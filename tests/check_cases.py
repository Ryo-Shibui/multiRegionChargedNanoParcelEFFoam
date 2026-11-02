#!/usr/bin/env python3
"""Numerical acceptance checks for dilute-tracer and legacy diffusion modes."""
from pathlib import Path
import argparse,json,math,re
import numpy as np

def values(path,n=2000):
    text=path.read_text();u=re.search(r'internalField\s+uniform\s+([^;]+);',text)
    if u:
        a=np.fromstring(u[1].replace('(',' ').replace(')',' '),sep=' ')
        return np.full(n,a[0]) if len(a)==1 else np.tile(a,(n,1))
    u=re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    if not u:raise ValueError(f'Cannot parse internal field: {path}')
    a=np.fromstring(u[3].replace('(',' ').replace(')',' '),sep=' ')
    return a.reshape(-1,3) if u[1]=='vector' else a

def check(case):
    case=Path(case);meta=json.loads((case/'expected.json').read_text());n=meta['nx']*meta['ny']
    time=meta['endTime'];out=case/f'{time:g}';region=Path('.') if meta['singleRegion'] else Path('fluid')
    initial=case/'0'/region;final=out/region
    failures=[];checks={};results={}
    def test(key,value,tol):
        checks[key]=float(value)
        if value>tol:failures.append(f'{key}: {value:.12g} exceeds {tol:.12g}')
    if meta['variant']=='steadyField':
        for name in ['U','T','p','p_rgh','E','phiE','N2','DEG','H2O']:
            if name in ['E','phiE'] and not (final/name).exists():continue
            test(name+'FrozenMaxChange',np.max(abs(values(final/name,n)-values(initial/name,n))),1e-11)
    else:
        cosine=np.repeat(np.cos(math.pi*(np.arange(meta['ny'])+.5)/meta['ny']),meta['nx'])
        for name,p in meta['tracers'].items():
            y=values(final/name,n);initial_y=values(initial/name,n)
            continuum=p['mean']+p['continuumAmplitude']*cosine
            discrete=p['mean']+p['discreteEulerAmplitude']*cosine
            amplitude=float(np.dot(y-y.mean(),cosine)/np.dot(cosine,cosine))
            test(name+'MeanConservationError',abs(y.mean()-p['mean']),2e-12)
            test(name+'ContinuumMaxErrorOverInitialAmplitude',np.max(abs(y-continuum))/p['amplitude'],1e-4)
            test(name+'DiscreteEulerMaxErrorOverInitialAmplitude',np.max(abs(y-discrete))/p['amplitude'],1e-6)
            steps=round(time/meta['dt']);eigen=4*math.sin(math.pi/(2*meta['ny']))**2/(meta['length']/meta['ny'])**2
            inferred=((p['amplitude']/amplitude)**(1/steps)-1)/(meta['dt']*eigen)
            test(name+'InferredDRelativeError',abs(inferred-p['effectiveD'])/p['effectiveD'],1e-5)
            results[name]=dict(measuredAmplitude=amplitude,expectedContinuumAmplitude=p['continuumAmplitude'],expectedDiscreteEulerAmplitude=p['discreteEulerAmplitude'],inferredD=inferred,expectedD=p['effectiveD'])
        deg=values(final/'DEG',n);water=values(final/'H2O',n);background=values(final/'N2',n)
        test('massFractionSumMaxError',np.max(abs(deg+water+background-1)),2e-13)
        test('backgroundClosureMaxError',np.max(abs(background-(1-deg-water))),2e-13)
        test('backgroundMeanConservationError',abs(background.mean()-(1-meta['tracers']['DEG']['mean']-meta['tracers']['H2O']['mean'])),2e-12)
        if min(deg.min(),water.min(),background.min())<0:failures.append('Negative mass fraction')
        for name in ['U','T','p','p_rgh','E','phiE']:
            if name in ['E','phiE'] and not (final/name).exists():continue
            test(name+'FrozenMaxChange',np.max(abs(values(final/name,n)-values(initial/name,n))),1e-11)
    return dict(case=str(case),solver=meta['solver'],variant=meta['variant'],passed=not failures,checks=checks,tracers=results,failures=failures)

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('case');ap.add_argument('--output');a=ap.parse_args()
    result=check(a.case);path=Path(a.output) if a.output else Path(a.case)/'acceptance.json';path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['passed'] else 1)
if __name__=='__main__':main()