#!/usr/bin/env python3
"""Check native parcel positions, charge, force response and restart continuity."""
from pathlib import Path
import argparse,json,math,re,sys
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from check_cases import values

def payload(path):
    text=path.read_text();text=text[text.index('}')+1:]
    m=re.search(r'\b(\d+)\s*\((.*)\)',text,re.S)
    if not m:raise ValueError(f'Missing native list {path}')
    return int(m[1]),m[2]
def native(path,vector=False):
    count,text=payload(path);a=np.fromstring(text.replace('(',' ').replace(')',' '),sep=' ')
    if vector:a=a.reshape(-1,3)
    if len(a)!=count:raise ValueError(f'Native list length mismatch: {path}')
    return a
def positions(path):
    count,text=payload(path)
    xyz=np.array([list(map(float,v)) for v in re.findall(r'\(\s*([-+\d.eE]+)\s+([-+\d.eE]+)\s+([-+\d.eE]+)\s*\)',text)])
    if len(xyz)!=count:raise ValueError(f'Native position list mismatch: {path}')
    return xyz

def parcel(case,time=None):
    case=Path(case);m=json.loads((case/'expected.json').read_text())
    time=m['endTime'] if time is None else time
    cloud=case/f'{time:g}'/('' if m['singleRegion'] else 'fluid')/'lagrangian/reactingCloud1'
    return dict(position=positions(cloud/'positions'),velocity=native(cloud/'U',True),charge=native(cloud/'chargeNumber'),diameter=native(cloud/'d'))
def check(case):
    case=Path(case);m=json.loads((case/'expected.json').read_text());p=parcel(case)
    checks={};failures=[]
    def test(key,x,tol):
        checks[key]=float(x)
        if x>tol:failures.append(f'{key} {x:.12g} > {tol:.12g}')
    if len(p['position'])!=1:raise AssertionError(f'Expected exactly one native parcel, got {len(p["position"])}')
    displacement=p['position'][0]-np.array(m['initialPosition']);expect=m['expectedDisplacement']
    if abs(expect)>0:test('terminalDriftDisplacementRelativeError',abs(displacement[0]-expect)/abs(expect),.01)
    else:test('forceFreeDisplacementMax',np.max(abs(displacement)),1e-13)
    test('transverseDisplacementMax',np.max(abs(displacement[1:])),1e-13)
    if abs(m['terminalVelocity'])>0:test('terminalVelocityRelativeError',abs(p['velocity'][0,0]-m['terminalVelocity'])/abs(m['terminalVelocity']),1e-6)
    else:test('forceFreeVelocityMax',np.max(abs(p['velocity'])),1e-13)
    test('chargeNumberIOError',abs(p['charge'][0]-m['chargeNumber']),1e-14)
    test('diameterRelativeIOError',abs(p['diameter'][0]-m['diameter'])/m['diameter'],1e-14)
    region=Path('.') if m['singleRegion'] else Path('fluid')
    for field in ['U','T','p','p_rgh','N2','DEG','H2O','E','phiE']:
        f=case/f'{m["endTime"]:g}'/region/field;initial=case/'0'/region/field
        if f.exists() and initial.exists():test(field+'CarrierFrozenMaxChange',np.max(abs(values(f)-values(initial))),1e-10)
    if m['variant']=='restart':
        first=parcel(case,.00025);test('checkpointChargeNumberIOError',abs(first['charge'][0]-m['chargeNumber']),1e-14)
        continuous=case.parent/'positive'
        if (continuous/'.completed').exists():
            c=parcel(continuous)
            test('restartPositionMaxDifference',np.max(abs(p['position']-c['position'])),1e-12)
            test('restartVelocityMaxDifference',np.max(abs(p['velocity']-c['velocity'])),1e-12)
            test('restartChargeNumberDifference',np.max(abs(p['charge']-c['charge'])),1e-14)
    return dict(solver=m['solver'],variant=m['variant'],passed=not failures,position=p['position'][0].tolist(),velocity=p['velocity'][0].tolist(),displacement=displacement.tolist(),expectedDisplacement=expect,terminalVelocity=m['terminalVelocity'],checks=checks,failures=failures)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('case');a=ap.parse_args();r=check(a.case)
    (Path(a.case)/'acceptance.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2));raise SystemExit(0 if r['passed'] else 1)
if __name__=='__main__':main()