#!/usr/bin/env python3
"""Native no-flow comparison of a single retained deposit and a flying probe."""
import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import make_force_cases as force_cases
import make_cases as cases


SOLVER = "multiRegionChargedNanoParcelEFFoam"
X0 = 0.01013
Z0 = 0.000513
DIAMETER = 48e-9
FIELD = 1.25e6
DT = 5e-11
END = 9e-8
PROBE_TIME = 1e-8
LAMBDA = force_cases.LAMBDA
CC = 1 + 2 * LAMBDA / DIAMETER * (
    1.257 + 0.4 * math.exp(-1.1 * DIAMETER / (2 * LAMBDA))
)
SPEED = 300 * force_cases.Q * FIELD * CC / (3 * math.pi * cases.MU * DIAMETER)


def injection(name, charge, time):
    return f"""
        {name}
        {{
            type manualInjection;
            positionsFile {name}Positions;
            SOI {time:.17g};
            parcelBasisType fixed;
            nParticle 1;
            massTotal 0;
            U0 (0 {-SPEED:.17g} 0);
            parcelCharge {charge};
            sizeDistribution
            {{
                type fixedValue;
                fixedValueDistribution {{ value {DIAMETER:.17g}; }}
            }}
        }}
"""


def write_control(case, end, latest=False):
    start = "latestTime" if latest else "startTime"
    cases.write(case / "system/controlDict", f"""
application {SOLVER};
startFrom {start}; startTime 0; stopAt endTime;
endTime {end:.17g}; deltaT {DT:.17g};
writeControl timeStep; writeInterval 10000; purgeWrite 0;
writeFormat ascii; writePrecision 17; writeCompression off;
timeFormat general; timePrecision 14; runTimeModifiable false;
adjustTimeStep no; maxCo .5; maxDi 10;
""")


def prepare(case, base, offset_nm, enabled, seed_charge):
    for part in ("0", "constant", "system"):
        shutil.copytree(base / part, case / part)
    electric = case / "0/fluid/E"
    electric.write_text(electric.read_text().replace("(50000 0 0)", f"(0 {-FIELD:g} 0)"))
    cloud = (cases.ROOT / "templates/reactingCloud1Properties").read_text()
    start = cloud.index("        injectionModels")
    stop = cloud.index("        dispersionModel", start)
    cloud = cloud[:start] + (
        "        injectionModels\n        {\n"
        + injection("seed", seed_charge, 0)
        + injection("probe", 300, PROBE_TIME)
        + "        }\n"
    ) + cloud[stop:]
    if enabled:
        cloud += "\ndepositedParticleInteractions { enabled true; maxRelativeTravel 0.1; }\n"

    constant = case / "constant/fluid"
    (constant / "reactingCloud1Properties").write_text(cloud)
    cases.write(constant / "seedPositions", f"1\n(({X0:.17g} 2e-7 {Z0:.17g}))\n", "vectorField")
    probe_x = X0 + offset_nm * 1e-9
    cases.write(constant / "probePositions", f"1\n(({probe_x:.17g} 1.2e-6 {Z0:.17g}))\n", "vectorField")


def run_solver(case, executable, log_name):
    path = case / log_name
    with path.open("w") as log:
        result = subprocess.run(
            [executable, "-case", str(case)], stdout=log, stderr=subprocess.STDOUT
        )
    if result.returncode:
        raise RuntimeError(f"{path}\n{path.read_text(errors='replace')[-5000:]}")


def read_rows(case, filename):
    rows = []
    for path in sorted((case / "postProcessing").rglob(filename)):
        with path.open(newline="") as stream:
            rows.extend(csv.DictReader(stream))
    return rows


def deposition_sites(case):
    times = [
        path for path in case.iterdir()
        if path.is_dir() and re.fullmatch(r"\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", path.name)
    ]
    checkpoint = max(times, key=lambda path: float(path.name))
    path = checkpoint / "fluid/lagrangian/reactingCloud1/depositionSite"
    body = path.read_text().split("}", 1)[1]
    match = re.search(r"\b(\d+)\s*\(([^()]*)\)", body)
    assert match, f"{path}: expected a native label list"
    sites = list(map(int, match.group(2).split()))
    assert len(sites) == int(match.group(1)), path
    return sorted(sites)


def run(case, base, executable, offset_nm, enabled, seed_charge=300, restart=False):
    prepare(case, base, offset_nm, enabled, seed_charge)
    if restart:
        write_control(case, 8e-9)
        run_solver(case, executable, "log.first")
    write_control(case, END, latest=restart)
    run_solver(case, executable, "log.solver")

    hits = sorted(read_rows(case, "deposition.csv"), key=lambda hit: float(hit["time"]))
    if len(hits) != 2:
        raise AssertionError(f"{case}: expected exactly two deposits, got {hits}")
    seed, probe = hits
    assert seed["patch"] == "lowerElectrode" and float(seed["time"]) < PROBE_TIME, hits
    assert float(probe["time"]) > PROBE_TIME, hits
    probe_id = (probe["origProc"], probe["origId"])
    trace = [
        row for row in read_rows(case, "trajectory.csv")
        if (row["origProc"], row["origId"]) == probe_id
    ]
    assert trace, f"{case}: probe trajectory missing"
    with (case / "probe_trajectory.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=trace[0].keys())
        writer.writeheader()
        writer.writerows(trace)
    return {
        "time_ns": (float(probe["time"]) - PROBE_TIME) * 1e9,
        "x_nm": (float(probe["x"]) - float(seed["x"])) * 1e9,
        "y_nm": float(probe["y"]) * 1e9,
        "patch": probe["patch"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "depositionInteraction")
    args = parser.parse_args()
    executable = shutil.which(SOLVER)
    if executable is None:
        parser.error(f"load OpenCFD and build {SOLVER} before running the native test")
    args.output.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="run-", dir=args.output)).resolve()
    base = cases.ROOT / "force/cases" / SOLVER / "positive"
    if not base.exists():
        force_cases.generate(SOLVER, "positive")

    results = []
    for offset, expected in ((250, 622), (500, 668), (1000, 1078)):
        off = run(output / f"off-{offset}", base, executable, offset, False)
        on = run(output / f"on-{offset}", base, executable, offset, True)
        assert off["patch"] == on["patch"] == "lowerElectrode", (off, on)
        assert abs(off["x_nm"] - offset) < 0.1, off
        assert abs(on["x_nm"] - expected) < 12 and on["time_ns"] > off["time_ns"], (off, on)
        results.append({
            "offset_nm": offset,
            "off": off,
            "on": on,
            "shift_nm": on["x_nm"] - off["x_nm"],
        })

    contact_case = output / "contact"
    contact = run(contact_case, base, executable, 0, True, seed_charge=0)
    assert contact["patch"] == "depositedParticle", contact
    assert abs(contact["y_nm"] - 72) < 0.1, contact
    original_sites = deposition_sites(contact_case)
    assert original_sites == [1, 2], (contact_case, original_sites)
    write_control(contact_case, END + 4 * DT, latest=True)
    run_solver(contact_case, executable, "log.contactRestart")
    restarted_sites = deposition_sites(contact_case)
    assert restarted_sites == original_sites, (contact_case, original_sites, restarted_sites)
    restart = run(output / "restart", base, executable, 500, True, restart=True)
    reference = results[1]["on"]
    assert restart["patch"] == reference["patch"], (restart, reference)
    assert abs(restart["x_nm"] - reference["x_nm"]) < 0.1, (restart, reference)

    report = {
        "comparisons": results,
        "neutral_seed_contact": contact,
        "retained_contact_after_restart": True,
        "restart": restart,
    }
    (output / "acceptance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    print(f"Native trajectories and acceptance report: {output}")


if __name__ == "__main__":
    main()
