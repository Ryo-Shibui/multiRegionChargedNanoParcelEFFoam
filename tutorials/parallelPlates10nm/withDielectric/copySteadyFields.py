#!/usr/bin/env python3
"""Copy the solved carrier state and native mesh into physical time zero."""
from pathlib import Path
import json, re, shutil
pair = Path(__file__).resolve().parent
steady, particles = pair / "steady", pair / "particles"
times = [p for p in steady.iterdir() if p.is_dir() and re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", p.name)]
latest = max(times, key=lambda p: float(p.name))
if float(latest.name) < 2000:
    raise SystemExit("Run steady/Allrun through iteration 2000 first")
regions = ["fluid"] + (["dielectric"] if (steady / "constant/dielectric").exists() else [])
copied = []
for region in regions:
    shutil.copytree(steady / "constant" / region / "polyMesh", particles / "constant" / region / "polyMesh")
    for source in (latest / region).iterdir():
        if source.is_file():
            shutil.copy2(source, particles / "0" / region / source.name)
            copied.append(f"{region}/{source.name}")
    for name in ("C", "Cx", "Cy", "Cz"):
        source = steady / "0" / region / name
        if source.exists():
            shutil.copy2(source, particles / "0" / region / name)
(particles / "steady_copy.json").write_text(json.dumps({"sourceTime": latest.name, "fields": copied}, indent=2) + "\n")
print(f"Copied solved fields from iteration {latest.name} into particles/0")
