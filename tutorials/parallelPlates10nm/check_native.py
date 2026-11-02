#!/usr/bin/env python3
from pathlib import Path
import csv,importlib.util,json,math,re,sys
import numpy as np
root=Path(__file__).resolve().parent
import check_thermal as v
case=Path(sys.argv[1]).resolve();m=json.loads((case/'expected.json').read_text());out=v.latest(case)
if float(out.name)<=0:raise SystemExit('No nonzero native output')
r={};n=2000;c=v.values(case/'constant/fluid/polyMesh/C',n);u=v.values(out/'fluid/U',n)
eta=(c[:,1]-(.001 if m['solid'] else 0))/.02;ua=-6*.01*eta*(1-eta)
v.assert_error(r,'velocityRelativeL2Error',np.linalg.norm(u[:,0]-ua)/np.linalg.norm(ua),.02)
v.assert_error(r,'maximumCrossFlow',np.max(abs(u[:,1:])),2e-5)
e=v.values(out/'fluid/E',n)
v.assert_error(r,'electricRelativeError',np.max(abs(e[:,1]+m['expectedField']))/m['expectedField'],.001)
if m['solid']:
    es=v.values(out/'dielectric/E',100)
    v.assert_error(r,'solidElectricRelativeError',np.max(abs(es[:,1]+m['expectedField']/4))/(m['expectedField']/4),.001)
def patchflux(path,patch):
    text=path.read_text();b=re.search(r'\b'+patch+r'\s*\{([^{}]*)\}',text,re.S)[1]
    q=re.search(r'\bvalue\s+uniform\s+([^;]+);',b)
    if q:return np.full(100,float(q[1]))
    q=re.search(r'\bvalue\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;',b,re.S)
    return np.fromstring(q[1],sep=' ')
fi=patchflux(out/'fluid/phi','inlet').sum();fo=patchflux(out/'fluid/phi','outlet').sum()
v.assert_error(r,'massFluxRelativeImbalance',abs(fi+fo)/abs(fi),1e-6)
if m['phase']=='particles':
    assert (case/'steady_copy.json').exists(),'Particle carrier fields must be copied from steady result'
    for reg,nn in [('fluid',2000)]+([('dielectric',100)] if m['solid'] else []):
        fields=['T','p','phiE','E']+(['U','p_rgh','N2','O2','H2O','rho','h','phi'] if reg=='fluid' else ['h'])
        for f in fields:
            if not (case/'0'/reg/f).exists():continue
            v.assert_error(r,f'{reg}/{f}FrozenMaxChange',np.max(abs(v.values(out/reg/f,nn)-v.values(case/'0'/reg/f,nn))),1e-10)
    rows=[]
    for f in (case/'postProcessing').rglob('deposition.csv'):
        with f.open() as fh:rows+=list(csv.DictReader(fh))
    assert len(rows)==1,('Expected one deposition',rows)
    time=float(rows[0]['time']);r['depositionTime']=time
    v.assert_error(r,'depositionTimeAbsoluteError',abs(time-m['expectedDepositionTime']),1e-7)
    assert rows[0]['patch']==('fluid_to_dielectric' if m['solid'] else 'lowerElectrode')
    assert abs(float(rows[0]['chargeNumber'])-1)<1e-12
result=dict(case=str(case),time=float(out.name),checks=r)
(case/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
