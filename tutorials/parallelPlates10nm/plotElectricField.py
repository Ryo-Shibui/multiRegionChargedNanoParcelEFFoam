#!/usr/bin/env python3
"""Plot native electric-field strength on the two actual structured meshes."""
from pathlib import Path
import numpy as np
from check_results import latest, values, read_rows

def cell_edges(centres):
    """Recover faces of the tutorial's uniform Cartesian mesh from native C."""
    if len(centres) < 2:
        raise ValueError("At least two cells are required along the plotted axis")
    middle = (centres[:-1] + centres[1:]) / 2
    return np.r_[centres[0] - (centres[1] - centres[0]) / 2,
                 middle,
                 centres[-1] + (centres[-1] - centres[-2]) / 2]

def native_grid(case, time, region):
    """Return actual cell magnitudes; no analytical field is used for pixels."""
    count = 2000 if region == "fluid" else 100
    centres = values(case / "0" / region / "C", count)
    electric = values(time / region / "E", count)
    # Coordinates differing only at floating-point roundoff identify one column.
    x = np.round(centres[:, 0], 12)
    y = np.round(centres[:, 1], 12)
    ux, ix = np.unique(x, return_inverse=True)
    uy, iy = np.unique(y, return_inverse=True)
    if len(ux) * len(uy) != count:
        raise ValueError(f"{case}/{region}: expected one planar structured mesh")
    magnitude = np.linalg.norm(electric, axis=1) / 1000
    grid = np.full((len(uy), len(ux)), np.nan)
    grid[iy, ix] = magnitude
    if not np.isfinite(grid).all():
        raise ValueError(f"{case}/{region}: incomplete native cell data")
    return cell_edges(ux) * 1000, cell_edges(uy) * 1000, grid

def plot_electric_field(root=None):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.patheffects as effects
        from matplotlib.colors import Normalize
    except ImportError:
        print("matplotlib unavailable: skipping electric-field strength PNG")
        return None
    root = Path(root) if root is not None else Path(__file__).resolve().parent
    output = root / "results"
    output.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 6.2), constrained_layout=True)
    norm = Normalize(vmin=0, vmax=50)
    artist = None
    for axis, variant in zip(axes, ("noDielectric", "withDielectric")):
        pair = root / variant
        case, particles = pair / "steady", pair / "particles"
        time = latest(case)
        coated = variant == "withDielectric"
        regions = ["fluid", "dielectric"] if coated else ["fluid"]
        for region in regions:
            xe, ye, strength = native_grid(case, time, region)
            artist = axis.pcolormesh(xe, ye, strength, shading="flat",
                                     cmap="viridis", norm=norm, rasterized=True)
            if region == "fluid":
                axis.text(1, 17.9 if coated else 17,
                          f"Gas: mean |E| = {strength.mean():.3f} kV/m",
                          fontsize=10, color="#151515",
                          bbox=dict(facecolor="white", alpha=.85, edgecolor="none", pad=3))
            else:
                axis.text(10, .5, f"Dielectric: {strength.mean():.3f} kV/m",
                          ha="center", va="center", fontsize=8, color="white")
        top = 21 if coated else 20
        axis.plot([0, 20], [0, 0], color="black", linewidth=3)
        axis.plot([0, 20], [top, top], color="black", linewidth=3)
        if coated:
            axis.axhline(1, color="white", linewidth=.9)
        trajectory = read_rows(particles, "trajectory.csv")
        deposition = read_rows(particles, "deposition.csv")
        if trajectory and deposition:
            event = deposition[0]
            deposited_time = float(event["time"])
            rows = sorted((r for r in trajectory if float(r["time"]) <= deposited_time + 1e-12),
                          key=lambda r: float(r["time"]))
            path_x = np.array([float(r["x"]) for r in rows]) * 1000
            path_y = np.array([float(r["y"]) for r in rows]) * 1000
            line, = axis.plot(path_x, path_y, color="black", linewidth=1.8,
                              label="Native particle trajectory", zorder=5)
            line.set_path_effects([effects.Stroke(linewidth=3.4, foreground="white"), effects.Normal()])
            axis.plot(path_x[0], path_y[0], "o", markerfacecolor="white",
                      markeredgecolor="black", markersize=5, zorder=6)
            axis.plot(float(event["x"]) * 1000, float(event["y"]) * 1000,
                      "X", color="white", markeredgecolor="black", markersize=8, zorder=6)
            axis.annotate(f"Stick: {deposited_time:.6f} s",
                          xy=(float(event["x"]) * 1000, float(event["y"]) * 1000),
                          xytext=(12.1, 4.2), fontsize=9,
                          arrowprops=dict(arrowstyle="-", color="black"),
                          bbox=dict(facecolor="white", alpha=.85, edgecolor="none", pad=3))
        carrier = values(time / "fluid/U", 2000)
        mean_x = carrier[:, 0].mean()
        start_x, end_x = (17, 12) if mean_x < 0 else (12, 17)
        axis.annotate("", xy=(end_x, 15), xytext=(start_x, 15),
                      arrowprops=dict(arrowstyle="->", color="black", lw=1.8))
        axis.text(14.5, 15.5, f"Mean Ux = {mean_x:.4f} m/s", ha="center", fontsize=9)
        axis.text(.5, top + .25, "1000 V electrode", va="bottom", fontsize=9)
        axis.text(.5, -.35, "0 V electrode", va="top", fontsize=9)
        axis.set(xlim=(0, 20), ylim=(-1, 22.5), xlabel="x (mm)", ylabel="y (mm)")
        axis.set_aspect("equal")
        axis.set_title("Bare electrode" if not coated else "1 mm dielectric, epsilonR = 4", fontsize=12)
    bar = fig.colorbar(artist, ax=list(axes), shrink=.82, pad=.025)
    bar.set_label("Native electric-field strength |E| (kV/m)")
    fig.suptitle("Solved electric field and charged-particle deposition", fontsize=14)
    path = output / "electricFieldStrength.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    print(f"Wrote {path}")
    return path

if __name__ == "__main__":
    plot_electric_field()
