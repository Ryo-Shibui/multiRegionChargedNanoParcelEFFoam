#!/usr/bin/env python3
from pathlib import Path
import importlib.util,shutil,json,re
pair=Path(__file__).resolve().parent
steady=pair/'steady';particles=pair/'particles'
outputs=[p for p in steady.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)]
latest=max(outputs,key=lambda p:float(p.name))
if float(latest.name)<2000:raise SystemExit('Run steady/Allrun through iteration2000 before copying fields')
spec=importlib.util.spec_from_file_location('archive',pair.parent/'archiveRun.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
mod.archive(particles)
meta=json.loads((particles/'expected.json').read_text())
regions=['fluid']+(['dielectric'] if meta['solid'] else [])
copied=[]
for reg in regions:
    source=latest/reg
    if not source.is_dir():raise SystemExit(f'Missing native output region {source}')
    for path in source.iterdir():
        if path.is_file():
            shutil.copy2(path,particles/'0'/reg/path.name);copied.append(f'{reg}/{path.name}')
(particles/'steady_copy.json').write_text(json.dumps(dict(source=str(latest),fields=copied),indent=2)+'\n')
print('Copied',len(copied),'native fields from',latest,'to physical time0')
