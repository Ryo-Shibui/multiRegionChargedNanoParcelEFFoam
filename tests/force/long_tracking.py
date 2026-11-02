#!/usr/bin/env python3
"""CHT low-Co normal-progress tracking regression (>100 substeps in one timestep)."""
from pathlib import Path
import argparse,shutil,json,re,math,subprocess,datetime
import make_force_cases as force
import check_force as native
ROOT=Path(__file__).resolve().parent
CASE=ROOT/'longTracking/case'
SOLVER='chtReactingNanoParcelFoam'
MAXCO=.0001;END=.0005

def generate():
    if CASE.exists():return
    src=ROOT/'cases'/SOLVER/'positive'
    if not src.exists():force.generate(SOLVER,'positive')
    for part in ['0','constant','system']:shutil.copytree(src/part,CASE/part)
    m=json.loads((src/'expected.json').read_text());velocity=m['terminalVelocity']
    p=CASE/'constant/fluid/reactingCloud1Properties';text=p.read_text()
    text=re.sub(r'\bmaxCo\s+[^;]+;',f'maxCo {MAXCO};',text)
    text=re.sub(r'\bU0\s*\([^;]*\);',f'U0 ({velocity:.17g} 0 0);',text);p.write_text(text)
    p=CASE/'system/controlDict';text=p.read_text()
    text=re.sub(r'\bdeltaT\s+[^;]+;',f'deltaT {END};',text);text=re.sub(r'\bwriteInterval\s+[^;]+;','writeInterval 1;',text);p.write_text(text)
    m.update(variant='longTracking',dt=END,expectedDisplacement=velocity*END,initialVelocity=[velocity,0,0],maxCo=MAXCO,cellVolume=.001*.0002*.001)
    m['minimumNormalSubsteps']=max(math.ceil(1/MAXCO),math.ceil(abs(m['expectedDisplacement'])/(MAXCO*m['cellVolume']**(1/3))))
    assert m['minimumNormalSubsteps']>100
    (CASE/'expected.json').write_text(json.dumps(m,indent=2)+'\n')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['generate','run','check']);a=ap.parse_args();generate()
    if a.action=='generate':print(CASE);print((CASE/'expected.json').read_text());return
    if a.action=='run':
        old=[p for p in CASE.iterdir() if (p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?(?:e[-+]?\d+)?',p.name,re.I) and float(p.name)>0) or p.name.startswith('log.') or p.name=='acceptance.json']
        if old:
            dest=ROOT/'longTracking/prior_runs'/datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f');dest.mkdir(parents=True)
            for p in old:shutil.move(str(p),dest/p.name)
        with (CASE/'log.solver').open('w') as f:code=subprocess.run([shutil.which(SOLVER) or SOLVER,'-case',str(CASE)],stdout=f,stderr=subprocess.STDOUT).returncode
        if code:r=dict(solver=SOLVER,variant='longTracking',passed=False,exitCode=code,failures=['Solver process failed'],logTail=(CASE/'log.solver').read_text()[-3000:])
        else:
            try:r=native.check(CASE);r['exitCode']=code
            except Exception as err:r=dict(solver=SOLVER,variant='longTracking',passed=False,exitCode=code,failures=[str(err)])
    else:r=native.check(CASE)
    m=json.loads((CASE/'expected.json').read_text());r['minimumNormalSubsteps']=m['minimumNormalSubsteps'];r['maxCo']=MAXCO
    (CASE/'acceptance.json').write_text(json.dumps(r,indent=2)+'\n');dest=ROOT/'longTracking/evidence';dest.mkdir(exist_ok=True)
    for p in [CASE/'expected.json',CASE/'acceptance.json',CASE/'log.solver']:
        if p.exists():shutil.copy2(p,dest/p.name)
    print(json.dumps(r,indent=2));raise SystemExit(0 if r['passed'] else 1)
if __name__=='__main__':main()