#!/usr/bin/env python3
from pathlib import Path
import datetime,shutil,re,sys

def archive(case):
    case=Path(case).resolve()
    support=Path(__file__).resolve().parents[2]
    if not case.is_relative_to(support):raise RuntimeError(f'Unexpected run target: {case}')
    items=[p for p in case.iterdir() if (p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name) and float(p.name)>0) or p.name=='postProcessing' or p.name=='acceptance.json' or p.name.startswith('log.')]
    if not items:return
    dest=case/'previous_runs'/datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f')
    dest.mkdir(parents=True)
    for p in items:
        if not p.resolve().is_relative_to(case):raise RuntimeError('Unexpected archive input')
        shutil.move(str(p),str(dest/p.name))
    print('Archived previous outputs:',dest)
if __name__=='__main__':archive(sys.argv[1])
