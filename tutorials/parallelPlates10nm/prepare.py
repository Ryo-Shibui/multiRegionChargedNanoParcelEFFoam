#!/usr/bin/env python3
"""Restore pristine inputs and remove only this tutorial's generated products."""
from pathlib import Path
import re, shutil, sys
case = Path(sys.argv[1]).resolve()
if not (case / "0.orig").is_dir():
    raise SystemExit("Expected 0.orig before preparing this tutorial")
for path in case.iterdir():
    if path.is_dir() and (re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", path.name) or path.name == "postProcessing"):
        shutil.rmtree(path)
for path in case.glob("constant/*/polyMesh"):
    shutil.rmtree(path)
for name in ("steady_copy.json", "acceptance.json"):
    (case / name).unlink(missing_ok=True)
shutil.copytree(case / "0.orig", case / "0")
