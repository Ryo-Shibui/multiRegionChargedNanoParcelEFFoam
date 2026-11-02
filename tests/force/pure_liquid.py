#!/usr/bin/env python3
"""Generate/run/check pure native H2O liquidEvapFuchsKnudsen runtime cases."""
from pathlib import Path
import argparse,sys,re,shutil,json,subprocess,math,datetime,traceback
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
import make_force_cases as force
import check_force as native
sys.path.insert(0,str(ROOT.parent))
import make_cases as species
END=1e-6;DT=1e-8;DIAMETER=1e-6

def generate(solver):
    case=ROOT/'pureLiquid/cases'/solver
    if case.exists():return case
    source=ROOT/'cases'/solver/'noCoulomb'
    if not source.exists():force.generate(solver,'noCoulomb')
    for part in ['0','constant','system']:shutil.copytree(source/part,case/part)
    single=solver==force.SOLVERS[0];zero=case/'0' if single else case/'0/fluid';constant=case/'constant' if single else case/'constant/fluid'
    (zero/'DEG').unlink(missing_ok=True)
    for name,value in [('N2',.99),('H2O',.01)]:species.fields(zero,name,[value]*2000)
    p=zero/'E';p.write_text(re.sub(r'internalField\s+uniform\s+[^;]+;','internalField uniform (0 0 0);',p.read_text()))
    (constant/'foam.inp').write_text('species (N2 H2O);\nreactions {}\n')
    p=constant/'foam.dat';p.write_text(re.sub(r'\nDEG\s*\{\s*specie\s*\{[^{}]*\}\s*equationOfState\s*\{[^{}]*\}\s*thermodynamics\s*\{[^{}]*\}\s*transport\s*\{[^{}]*\}\s*\}','',p.read_text(),flags=re.S))
    p=constant/'thermophysicalProperties';p.write_text(re.sub(r'solids\s*\{\s*NaCl\s*\{[^{}]*\}\s*\}','solids {}',p.read_text()))
    p=constant/'reactingCloud1Properties';cloud=p.read_text()
    cloud=cloud.replace('rho0                1000;','rho0                997;').replace('Cp0                 850;','Cp0                 4180;')
    cloud=re.sub(r'parcelCharge\s+[^;]+;','parcelCharge 0;',cloud)
    cloud=cloud.replace('value 1e-8;','value 1e-6;').replace('phaseChangeModel none;','phaseChangeModel liquidEvapFuchsKnudsen;')
    a=cloud.index('        singleMixtureFractionCoeffs');b=cloud.index('    cloudFunctions',a)
    cloud=cloud[:a]+'''        singleMixtureFractionCoeffs
        {
            phases ( gas {} liquid { H2O 1; } solid {} );
            YGasTot0 0;
            YLiquidTot0 1;
            YSolidTot0 0;
        }
        liquidEvapFuchsKnudsenCoeffs
        {
            enthalpyTransfer latentHeat;
            solution (H2O);
            gamma 6.5e-8;
            alpham 1;
            activityCoefficient Hoff;
            ic 0;
        }
    }

'''+cloud[b:]
    p.write_text(cloud)
    species.write(case/'system/controlDict',f'''application {solver};
startFrom startTime; startTime 0; stopAt endTime; endTime {END}; deltaT {DT};
writeControl timeStep; writeInterval 100; purgeWrite 0;
writeFormat ascii; writePrecision 16; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
adjustTimeStep no; maxCo .5; maxDi 10;
''')
    expected=dict(solver=solver,singleRegion=single,model='liquidEvapFuchsKnudsen',liquid='H2O',solution=['H2O'],initialDiameter=DIAMETER,rhoP=997,T=298.15,gasH2O=.01,endTime=END,dt=DT,scope='Construction and benign pure-water evaporation runtime/finite-state regression; no DEG or predictive evaporation-accuracy claim')
    (case/'expected.json').write_text(json.dumps(expected,indent=2)+'\n')
    assert 'NaCl' not in (constant/'reactingCloud1Properties').read_text()
    assert 'NaCl' not in (constant/'thermophysicalProperties').read_text()
    assert 'DEG' not in (constant/'foam.dat').read_text()
    return case

def check(case):
    m=json.loads((case/'expected.json').read_text());region=Path('.') if m['singleRegion'] else Path('fluid')
    path=case/f'{END:g}'/region/'lagrangian/reactingCloud1'
    log=(case/'log.solver').read_text();failures=[];checks={}
    if 'liquidEvapFuchsKnudsen' not in log:failures.append('Custom phase change model was not selected')
    if not re.search(r'liquid:\s*H2O',log):failures.append('Pure water participating liquid was not reported')
    for name in ['d','rho','T','mass0','Yliquid','Ysolid','Ygas','YH2O(l)']:
        a=native.native(path/name)
        if len(a)!=1 or not np.isfinite(a).all():failures.append(f'{name}: expected one finite native parcel value')
        checks[name]=float(a[0]) if len(a) else None
    mass=checks['rho']*math.pi*checks['d']**3/6
    checks['finalMass']=mass;checks['massRatio']=mass/checks['mass0'];checks['diameterRatio']=checks['d']/DIAMETER
    if not 0<checks['d']<DIAMETER:failures.append('Droplet must remain positive and shrink in undersaturated gas')
    if not 0<mass<checks['mass0']:failures.append('Droplet must retain positive finite mass while evaporating')
    if not .95<checks['massRatio']<1:failures.append('Expected benign evaporation with less than5% mass change')
    if not 200<checks['T']<350:failures.append('Droplet temperature outside benign finite range')
    if abs(checks['Yliquid']-1)>1e-12 or abs(checks['YH2O(l)']-1)>1e-12:failures.append('Native parcel composition must remain pure liquid water')
    if abs(checks['Ysolid'])>1e-12 or abs(checks['Ygas'])>1e-12:failures.append('Unexpected solid or gas parcel phase')
    return dict(solver=m['solver'],variant='pureLiquidH2O',passed=not failures,checks=checks,failures=failures,scope=m['scope'])

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['generate','run','check']);ap.add_argument('--solver',choices=force.SOLVERS);a=ap.parse_args()
    results=[];evidence=ROOT/'pureLiquid/evidence';evidence.mkdir(parents=True,exist_ok=True)
    for solver in [a.solver] if a.solver else force.SOLVERS:
        case=generate(solver)
        if a.action=='generate':print(case);continue
        if a.action=='run':
            old=[p for p in case.iterdir() if (p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?(?:e[-+]?\d+)?',p.name,re.I) and float(p.name)>0) or p.name.startswith('log.') or p.name=='acceptance.json']
            if old:
                dest=evidence/'prior_runs'/solver/datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f');dest.mkdir(parents=True)
                for p in old:shutil.move(str(p),dest/p.name)
            with (case/'log.solver').open('w') as f:code=subprocess.run([shutil.which(solver) or solver,'-case',str(case)],stdout=f,stderr=subprocess.STDOUT).returncode
            if code:
                result=dict(solver=solver,variant='pureLiquidH2O',passed=False,exitCode=code,failures=['Solver process failed'],logTail=(case/'log.solver').read_text()[-4000:])
            else:
                try:result=check(case);result['exitCode']=code
                except Exception as err:result=dict(solver=solver,variant='pureLiquidH2O',passed=False,exitCode=code,failures=[str(err)],traceback=traceback.format_exc())
        else:result=check(case)
        (case/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
        dest=evidence/solver;dest.mkdir(exist_ok=True)
        for p in [case/'expected.json',case/'acceptance.json',case/'log.solver']:
            if p.exists():shutil.copy2(p,dest/p.name)
        results.append(result);print(solver,'pureLiquidH2O','PASS' if result['passed'] else 'FAIL',result['failures'],flush=True)
    if a.action!='generate':
        combined=[json.loads(p.read_text()) for p in evidence.glob('*/acceptance.json')]
        (evidence/'results.json').write_text(json.dumps(combined,indent=2)+'\n');raise SystemExit(0 if all(r['passed'] for r in results) else 1)
if __name__=='__main__':main()