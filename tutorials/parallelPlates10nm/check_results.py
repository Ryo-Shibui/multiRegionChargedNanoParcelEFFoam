#!/usr/bin/env python3
"""Verify native fields and tracking; export CSV comparisons and optional PNG plots."""
from pathlib import Path
import csv, json, math, re
import numpy as np
ROOT = Path(__file__).resolve().parent

def values(path, count):
    text = path.read_text()
    match = re.search(r"internalField\s+uniform\s+([^;]+);", text)
    if match:
        array = np.fromstring(match[1].replace("(", " ").replace(")", " "), sep=" ")
        return np.full(count, array[0]) if len(array) == 1 else np.tile(array, (count, 1))
    match = re.search(r"internalField\s+nonuniform\s+List<(\w+)>\s+\d+\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"Cannot read native field {path}")
    array = np.fromstring(match[2].replace("(", " ").replace(")", " "), sep=" ")
    return array.reshape(-1, 3) if match[1] == "vector" else array

def latest(case):
    paths = [p for p in case.iterdir() if p.is_dir() and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", p.name)]
    return max(paths, key=lambda p: float(p.name))

def patch_flux(path, name):
    text = re.search(r"\b" + name + r"\s*\{([^{}]*)\}", path.read_text(), re.S)[1]
    match = re.search(r"value\s+uniform\s+([^;]+);", text)
    if match:
        return 100 * float(match[1])
    match = re.search(r"value\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;", text, re.S)
    return np.fromstring(match[1], sep=" ").sum()

def read_rows(case, name):
    rows = []
    for path in (case / "postProcessing").rglob(name):
        with path.open() as stream:
            rows.extend(csv.DictReader(stream))
    return rows

def write_csv(path, names, rows):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(names)
        writer.writerows(rows)

def check():
    output = ROOT / "results"
    output.mkdir(exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        plt = None
        print("matplotlib unavailable: numerical checks and CSV outputs remain available")
    results = []
    for variant in ("noDielectric", "withDielectric"):
        pair = ROOT / variant
        steady, particles = pair / "steady", pair / "particles"
        meta = json.loads((particles / "expected.json").read_text())
        solid = meta["solid"]
        bottom = .001 if solid else 0
        gas_field = 1000 / (.02 + (.001 / 4 if solid else 0))
        # Parse the actual N2 transport dictionary used by this case.
        properties = (particles / "constant/fluid/foam.dat").read_text()
        nitrogen = properties.split("\nN2\n", 1)[1].split("\nO2\n", 1)[0]
        As = float(re.search(r"\bAs\s+([^;]+);", nitrogen)[1])
        Ts = float(re.search(r"\bTs\s+([^;]+);", nitrogen)[1])
        mu = As * 320**1.5 / (320 + Ts)
        diameter, mean_free_path = 1e-8, 65e-9
        slip = 1 + 2 * mean_free_path / diameter * (1.257 + .4 * math.exp(-1.1 * diameter / (2 * mean_free_path)))
        actual_gas_field = float(-values(particles / "0/fluid/E", 2000)[:, 1].mean())
        speed = 1.602176634e-19 * actual_gas_field * slip / (3 * math.pi * mu * diameter)
        deposition_theory = .01 / speed
        checks = {}
        def require(name, error, tolerance):
            checks[name] = float(error)
            if not np.isfinite(error) or error >= tolerance:
                raise AssertionError(f"{variant}: {name} = {error} >= {tolerance}")
        flow_output = latest(steady)
        c = values(steady / "0/fluid/C", 2000)
        u = values(flow_output / "fluid/U", 2000)
        eta = (c[:, 1] - bottom) / .02
        ua = -6 * .01 * eta * (1 - eta)
        require("flowRelativeL2Error", np.linalg.norm(u[:, 0] - ua) / np.linalg.norm(ua), .02)
        require("crossFlowMaximum", np.max(abs(u[:, 1:])), 2e-5)
        incoming = patch_flux(flow_output / "fluid/phi", "inlet")
        outgoing = patch_flux(flow_output / "fluid/phi", "outlet")
        require("massFluxRelativeImbalance", abs(incoming + outgoing) / abs(incoming), 1e-6)
        particle_output = latest(particles)
        profiles = []
        for region, count in [("fluid", 2000)] + ([("dielectric", 100)] if solid else []):
            cr = values(steady / "0" / region / "C", count)
            electric = values(flow_output / region / "E", count)
            potential = values(flow_output / region / "phiE", count)
            magnitude = gas_field if region == "fluid" else gas_field / 4
            require(region + "FieldRelativeError", np.max(abs(electric[:, 1] + magnitude)) / magnitude, .002)
            analytical = gas_field / 4 * cr[:, 1] if region == "dielectric" else gas_field * (cr[:, 1] - bottom) + (gas_field * .001 / 4 if solid else 0)
            require(region + "PotentialRelativeError", np.max(abs(potential - analytical)) / 1000, .002)
            fields = ("E", "phiE", "T", "p") + (("U", "p_rgh", "N2", "H2O", "rho", "h", "phi") if region == "fluid" else ("h",))
            for field in fields:
                initial = particles / "0" / region / field
                final = particle_output / region / field
                if initial.exists() and final.exists():
                    initial_values = values(initial, count)
                    difference = np.max(abs(values(final, count) - initial_values))
                    if field == "rho":
                        # heRhoThermo reconstructs rho from p,T at restart.
                        require(region + "/rhoRestartRelativeChange", difference / np.max(abs(initial_values)), 1e-8)
                    else:
                        require(region + "/" + field + "FrozenMaxChange", difference, 1e-10)
            # Average cells at the same height to make a compact physical profile.
            for y in sorted(set(cr[:, 1])):
                mask = cr[:, 1] == y
                profiles.append((region, y, potential[mask].mean(), analytical[mask].mean(), electric[mask, 1].mean(), -magnitude))
        write_csv(output / (variant + "_fieldProfiles.csv"), ["region", "y_m", "phiE_V", "phiE_theory_V", "Ey_V_per_m", "Ey_theory_V_per_m"], profiles)
        deposit = read_rows(particles, "deposition.csv")
        if len(deposit) != 1:
            raise AssertionError(f"Expected one deposited parcel: {deposit}")
        event = deposit[0]
        assert event["patch"] == ("fluid_to_dielectric" if solid else "lowerElectrode")
        require("depositionChargeError", abs(float(event["chargeNumber"]) - 1), 1e-12)
        observed = float(event["time"])
        require("depositionTimeRelativeError", abs(observed / deposition_theory - 1), .002)
        trajectory = read_rows(particles, "trajectory.csv")
        samples = []
        for row in trajectory:
            time = float(row["time"])
            if time > observed + 1e-12:
                continue
            x, y = float(row["x"]), float(row["y"])
            x_theory = .01013 - 1.5 * .01 * time + 2 * .01 * speed**2 / .02**2 * time**3
            y_theory = bottom + .01 - speed * time
            samples.append((time, x, y, x_theory, y_theory))
        data = np.array(samples)
        if not len(data):
            raise AssertionError("No native trajectory rows")
        require("verticalTrajectoryRelativeError", np.max(abs(data[:, 2] - data[:, 4])) / .01, .002)
        require("horizontalTrajectoryRelativeError", np.max(abs(data[:, 1] - data[:, 3])) / .01, .002)
        write_csv(output / (variant + "_trajectoryComparison.csv"), ["time_s", "x_m", "y_m", "x_theory_m", "y_theory_m"], samples)
        report = dict(variant=variant, temperature_K=320, mu_Pa_s=mu, slipCorrection=slip, gasField_V_per_m=gas_field, copiedGasField_V_per_m=actual_gas_field, solidField_V_per_m=gas_field/4 if solid else None, electricSpeed_m_per_s=speed, depositionTheory_s=deposition_theory, depositionNative_s=observed, inletMassFlux_kg_per_s=incoming, outletMassFlux_kg_per_s=outgoing, checks=checks)
        results.append(report)
        if plt:
            fig, axes = plt.subplots(1, 3, figsize=(13, 4), constrained_layout=True)
            for region in ("fluid", "dielectric"):
                subset = [r for r in profiles if r[0] == region]
                if subset:
                    a = np.array([r[1:] for r in subset])
                    axes[0].plot(a[:, 0]*1000, a[:, 1], label=region+" native")
                    axes[0].plot(a[:, 0]*1000, a[:, 2], "--", label=region+" theory")
            axes[0].set(xlabel="Height y (mm)", ylabel="Potential (V)")
            axes[0].legend(fontsize=7)
            axes[1].plot(data[:, 1]*1000, data[:, 2]*1000, label="Native parcel")
            axes[1].plot(data[:, 3]*1000, data[:, 4]*1000, "--", label="Terminal-drift theory")
            axes[1].set(xlabel="x (mm)", ylabel="y (mm)")
            axes[1].legend(fontsize=7)
            axes[2].plot(c[:, 1]*1000, u[:, 0], ".", markersize=2, label="Solved flow")
            axes[2].plot(c[:, 1]*1000, ua, ".", markersize=1, label="Poiseuille theory")
            axes[2].set(xlabel="Height y (mm)", ylabel="Ux (m/s)")
            axes[2].legend(fontsize=7)
            fig.suptitle(variant)
            fig.savefig(output / (variant + "_comparison.png"), dpi=180)
            plt.close(fig)
    if plt:
        from plotElectricField import plot_electric_field
        plot_electric_field(ROOT)
    (output / "acceptance.json").write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results, indent=2))
if __name__ == "__main__":
    check()
