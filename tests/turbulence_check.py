#!/usr/bin/env python3
"""Frozen k-epsilon viscosity tests for molecular D plus turbulent mixing."""
from pathlib import Path
import sys,shutil,subprocess,json,math,re
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import check_cases
NAMES=['reactingChargedNanoParcelFoam','chtReactingNanoParcelFoam','multiRegionChargedNanoParcelEFFoam']
results=[]
for solver in (sys.argv[1:] or NAMES):
    case=ROOT/'turbulence'/solver
    if case.exists():
        archive=ROOT/'evidence'/'turbulence_archive'
        archive.mkdir(parents=True,exist_ok=True)
        import time
        shutil.move(str(case),archive/(solver+'_'+str(time.time_ns())))
    shutil.copytree(ROOT/'cases'/solver/'dilute',case,ignore=shutil.ignore_patterns('0.*','postProcessing','log.*','acceptance.json'))
    single=solver==NAMES[0]; region=Path('.') if single else Path('fluid')
    props=case/'constant'/region
    (props/'turbulenceProperties').write_text('''FoamFile { version 2.0; format ascii; class dictionary; object turbulenceProperties; }
simulationType RAS;
RAS { RASModel kEpsilon; turbulence on; printCoeffs off; }
''')
    for name,dims,val in [('alphat','1 -1 -1 0 0 0 0',1.204*9e-7/.85),('nut','0 2 -1 0 0 0 0',9e-7),('k','0 2 -2 0 0 0 0',.001),('epsilon','0 2 -3 0 0 0 0',.1)]:
        (case/'0'/region/name).write_text('FoamFile {version 2.0; format ascii; class volScalarField; object '+name+';}\ndimensions ['+dims+'];\ninternalField uniform '+str(val)+''';
boundaryField
{
inlet {type zeroGradient;}
outlet {type zeroGradient;}
lowerElectrode {type zeroGradient;}
upperElectrode {type zeroGradient;}
front {type symmetryPlane;}
back {type symmetryPlane;}
}
''')
    dp=props/'speciesDiffusionProperties'
    dp.write_text(re.sub(r'turbulentSchmidt\s+[^;]+;','turbulentSchmidt 0.3;',dp.read_text()))
    schemes=case/'system'/region/'fvSchemes'
    schemes.write_text(re.sub(r'laplacianSchemes\s*\{[^}]+\}', '''laplacianSchemes
{
 default none;
 laplacian(muEff,DEG) Gauss linear orthogonal;
 laplacian(muEff,H2O) Gauss linear orthogonal;
}''',schemes.read_text()))
    meta=json.loads((case/'expected.json').read_text());meta['variant']='turbulentMixing'
    eigen=4*math.sin(math.pi/(2*meta['ny']))**2/(meta['length']/meta['ny'])**2
    for p in meta['tracers'].values():
        d=p['D']+9e-7/.3;p['effectiveD']=d
        p['continuumAmplitude']=p['amplitude']*math.exp(-d*(math.pi/meta['length'])**2*meta['endTime'])
        p['discreteEulerAmplitude']=p['amplitude']/(1+d*meta['dt']*eigen)**round(meta['endTime']/meta['dt'])
    (case/'expected.json').write_text(json.dumps(meta,indent=2)+'\n')
    with (case/'log.solver').open('w') as log:
        run=subprocess.run([solver,'-case',str(case)],stdout=log,stderr=subprocess.STDOUT)
    try:r=check_cases.check(case)
    except Exception as e:r=dict(passed=False,failures=[str(e)])
    r.update(solver=solver,exitCode=run.returncode)
    r['passed']=r['passed'] and run.returncode==0
    results.append(r)
    print(solver,'turbulent mixing and explicit laplacian scheme', 'PASS' if r['passed'] else 'FAIL',r.get('failures'),flush=True)
dest=ROOT/'evidence'/'turbulence';dest.mkdir(parents=True,exist_ok=True)
(dest/'results.json').write_text(json.dumps(results,indent=2)+'\n')
raise SystemExit(0 if all(r['passed'] for r in results) else 1)
