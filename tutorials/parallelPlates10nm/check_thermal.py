#!/usr/bin/env python3
"""Validate the generated thermal/species/switch/native-cloud acceptance runs."""
from pathlib import Path
import argparse,csv,json,math,re
import numpy as np
ROOT=Path(__file__).resolve().parent

def values(path,n):
    s=path.read_text()
    u=re.search(r'internalField\s+uniform\s+([^;]+);',s)
    if u:
        a=np.fromstring(u[1].replace('(',' ').replace(')',' '),sep=' ')
        return np.full(n,a[0]) if len(a)==1 else np.tile(a,(n,1))
    u=re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s+(\d+)\s*\((.*?)\)\s*;',s,re.S)
    if not u:raise ValueError(f'Cannot parse {path}')
    a=np.fromstring(u[3].replace('(',' ').replace(')',' '),sep=' ')
    return a.reshape(-1,3) if u[1]=='vector' else a

def latest(case):
    ds=[p for p in case.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)]
    return max(ds,key=lambda p:float(p.name))

def assert_error(result,key,error,tol):
    result[key]=float(error)
    assert error<tol,(key,error,tol)

def check(case):
    m=json.loads((case/'expected.json').read_text());out=latest(case)
    assert float(out.name)>0,f'{case}: no nonzero output'
    result={'case':case.name,'time':float(out.name),'checks':{}}
    r=result['checks'];n=2000
    cf=values(case/'constant/fluid/polyMesh/C',n)
    if m['mode']=='heat':
        q=m['heatFlux']
        tf=values(out/'fluid/T',n)
        ta=m['coldT']+q*(.021-cf[:,1])/m['gasKappa']
        assert_error(r,'fluidMaxTemperatureErrorK',np.max(abs(tf-ta)),.01)
        cs=values(case/'constant/dielectric/polyMesh/C',100)
        ts=values(out/'dielectric/T',100)
        sa=m['hotT']-q*cs[:,1]
        assert_error(r,'solidMaxTemperatureErrorK',np.max(abs(ts-sa)),.01)
        qf=m['gasKappa']*(tf[-20:].mean()-m['coldT'])/(.0001)
        qs=(m['hotT']-ts[:20].mean())/(.0001)
        assert_error(r,'fluidHeatFluxRelativeError',abs(qf-q)/q,.002)
        assert_error(r,'solidHeatFluxRelativeError',abs(qs-q)/q,.002)
        assert_error(r,'frozenVelocityMaxChange',np.max(abs(values(out/'fluid/U',n))),1e-14)
        assert_error(r,'speciesMaxChange',np.max(abs(values(out/'fluid/O2',n)-.21)),1e-14)
    elif m['mode']=='species':
        y=values(out/'fluid/O2',n)
        analytic=.21+m['expectedAmplitude']*np.cos(math.pi*cf[:,1]/.02)
        assert_error(r,'maxO2AnalyticError',np.max(abs(y-analytic)),1e-5)
        assert_error(r,'meanO2ConservationError',abs(y.mean()-.21),1e-12)
        total=sum(values(out/'fluid'/s,n) for s in ('N2','O2','H2O'))
        assert_error(r,'massFractionSumError',np.max(abs(total-1)),1e-13)
        for f in ['U','T','p','p_rgh','phiE','E']:
            assert_error(r,f+'FrozenMaxChange',np.max(abs(values(out/'fluid'/f,n)-values(case/'0/fluid'/f,n))),1e-10)
    elif m['mode']=='off':
        for reg,nn in [('fluid',2000)]+([('dielectric',100)] if m['solid'] else []):
            fields=['T','phiE','E']+(['U','p','p_rgh','N2','O2','H2O'] if reg=='fluid' else [])
            for f in fields:
                a=values(out/reg/f,nn);b=values(case/'0'/reg/f,nn)
                if case.name=='solidOffFluidEnergy' and reg=='fluid' and f=='T':
                    r['fluidTemperatureChangeK']=float(np.max(abs(a-b)))
                    assert r['fluidTemperatureChangeK']>1e-5
                else:assert_error(r,reg+'/'+f+'FrozenMaxChange',np.max(abs(a-b)),1e-10)
    else:
        electric=values(out/'fluid/E',n)
        assert_error(r,'electricFieldRelativeError',np.max(abs(electric[:,1]+m['expectedField']))/m['expectedField'],.001)
        eta=(cf[:,1]-(.001 if m['solid'] else 0))/.02
        ua=-6*.01*eta*(1-eta)
        uf=values(out/'fluid/U',n)
        assert_error(r,'prescribedParabolicVelocityError',np.max(abs(uf[:,0]-ua)),1e-12)
        if (case/'tracking.enabled').exists():
            files=list((case/'postProcessing').rglob('deposition.csv'))
            rows=[]
            for f in files:
                with f.open() as fh:rows+=list(csv.DictReader(fh))
            assert len(rows)==1,('Expected exactly one deposition event',files,rows)
            t=float(rows[0]['time']);r['depositionTime']=t
            assert_error(r,'depositionTimeAbsoluteError',abs(t-m['expectedDepositionTime']),1e-7)
            assert rows[0]['patch']==('fluid_to_dielectric' if m['solid'] else 'lowerElectrode'),rows
            assert abs(float(rows[0]['chargeNumber'])-1)<1e-12
    (case/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));return result

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('cases',nargs='*')
    a=ap.parse_args()
    names=a.cases or ['heatConduction','speciesDiffusion','electricBare','electricDielectric','allOff','solidOffFluidEnergy']
    results=[check(ROOT/n) for n in names]
    (ROOT/'acceptance_results.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':main()