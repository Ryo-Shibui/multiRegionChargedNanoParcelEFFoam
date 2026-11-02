#!/usr/bin/env python3
"""Two-rank native parcel ownership-transfer and charge-IO regression (MR solver)."""
from pathlib import Path
import argparse,shutil,re,json,subprocess,os,math,traceback
import numpy as np
import make_force_cases as force
import check_force as native
ROOT=Path(__file__).resolve().parent
SOLVER='multiRegionChargedNanoParcelEFFoam'
CASE=ROOT/'mpiCrossing/case'
START=np.array([.00999,.01,.000513]);FIRST=.00005;END=.0005

def generate():
    if CASE.exists():return
    source=ROOT/'cases'/SOLVER/'positive'
    if not source.exists():force.generate(SOLVER,'positive')
    for part in ['0','constant','system']:shutil.copytree(source/part,CASE/part)
    (CASE/'constant/fluid/particlePositions').write_text('FoamFile {version 2.0; format ascii; class vectorField; object particlePositions;}\n1\n((.00999 .01 .000513))\n')
    p=CASE/'system/controlDict';p.write_text(p.read_text().replace('writeInterval 250;','writeInterval 50;'))
    (CASE/'system/decomposeParDict').write_text('FoamFile {version 2.0; format ascii; class dictionary; object decomposeParDict;}\nnumberOfSubdomains 2; method simple; simpleCoeffs {n (2 1 1); delta .001;}\n')
    shutil.copy2(CASE/'system/decomposeParDict',CASE/'system/fluid/decomposeParDict')
    meta=json.loads((source/'expected.json').read_text());meta['initialPosition']=START.tolist();meta['variant']='mpiCrossing';meta['subdomains']=2;meta['firstCheckpoint']=FIRST
    (CASE/'expected.json').write_text(json.dumps(meta,indent=2)+'\n')

def record(time):
    rows=[]
    for rank in [0,1]:
        cloud=CASE/f'processor{rank}'/f'{time:g}'/'fluid/lagrangian/reactingCloud1'
        if not (cloud/'positions').exists():continue
        positions=native.positions(cloud/'positions')
        if not len(positions):continue
        charge=native.native(cloud/'chargeNumber');origId=native.native(cloud/'origId');origProc=native.native(cloud/'origProcId');dia=native.native(cloud/'d');vel=native.native(cloud/'U',True)
        for i,pos in enumerate(positions):rows.append(dict(rank=rank,position=pos.tolist(),chargeNumber=float(charge[i]),origId=int(origId[i]),origProcId=int(origProc[i]),diameter=float(dia[i]),velocity=vel[i].tolist()))
    return rows

def check():
    first=record(FIRST);final=record(END);failures=[];checks={}
    if len(first)!=1 or len(final)!=1:failures.append(f'Expected one global parcel before/after crossing, got {len(first)}/{len(final)}')
    else:
        a,b=first[0],final[0]
        if a['rank']!=0 or b['rank']!=1:failures.append('Expected native ownership transfer from processor0 to processor1')
        if not a['position'][0]<.01<b['position'][0]:failures.append('Parcel did not physically cross x=.01 partition')
        if a['chargeNumber']!=1 or b['chargeNumber']!=1:failures.append('Elementary charge did not survive MPI transfer/native IO')
        if (a['origId'],a['origProcId'])!=(b['origId'],b['origProcId']):failures.append('Native original identity changed during processor transfer')
        meta=json.loads((CASE/'expected.json').read_text());displacement=b['position'][0]-START[0]
        checks['driftDisplacementRelativeError']=abs(displacement-meta['expectedDisplacement'])/abs(meta['expectedDisplacement'])
        if checks['driftDisplacementRelativeError']>.01:failures.append('Parallel drift differs from Stokes-Cunningham theory by more than1%')
        checks['diameterRelativeError']=abs(b['diameter']-1e-8)/1e-8
        if checks['diameterRelativeError']>1e-12:failures.append('Particle diameter changed during transfer')
    return dict(solver=SOLVER,variant='mpiCrossing',passed=not failures,firstTime=FIRST,finalTime=END,first=first,final=final,checks=checks,failures=failures)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['generate','run','check']);a=ap.parse_args();generate()
    if a.action=='generate':print(CASE);return
    if a.action=='run':
        if (CASE/'processor0').exists():raise SystemExit('Existing MPI result retained: run check, or archive the case before a new run')
        with (CASE/'log.decomposePar').open('w') as f:code=subprocess.run(['decomposePar','-case',str(CASE),'-allRegions'],stdout=f,stderr=subprocess.STDOUT).returncode
        if code:raise SystemExit('decomposePar failed; see log.decomposePar')
        env=os.environ.copy()
        if os.geteuid()==0:env.update(OMPI_ALLOW_RUN_AS_ROOT='1',OMPI_ALLOW_RUN_AS_ROOT_CONFIRM='1')
        with (CASE/'log.parallel').open('w') as f:code=subprocess.run(['mpirun','-np','2',shutil.which(SOLVER) or SOLVER,'-parallel','-case',str(CASE)],stdout=f,stderr=subprocess.STDOUT,env=env).returncode
        if code:
            result=dict(solver=SOLVER,variant='mpiCrossing',passed=False,exitCode=code,failures=['Parallel solver failed'],logTail=(CASE/'log.parallel').read_text()[-4000:])
        else:result=check();result['exitCode']=code
    else:result=check()
    (CASE/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n');dest=ROOT/'mpiCrossing/evidence';dest.mkdir(exist_ok=True)
    for p in list(CASE.glob('log.*'))+[CASE/'expected.json',CASE/'acceptance.json']:shutil.copy2(p,dest/p.name)
    print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
if __name__=='__main__':main()