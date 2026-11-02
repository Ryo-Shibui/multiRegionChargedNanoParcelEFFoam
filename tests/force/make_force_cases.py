#!/usr/bin/env python3
"""Generate native prescribed-field nanoparticle force and restart regressions."""
from pathlib import Path
import argparse,sys,shutil,re,json,math
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
import make_cases as species
SOLVERS=species.SOLVERS
VARIANTS=['positive','negative','neutral','noCoulomb','reverseField','restart','missingE']
D=1e-8;RHO=1000.;MU=1.81e-5;LAMBDA=6.5e-8;Q=1.602176634e-19;DT=1e-6;END=.0005;START=[.01013,.01,.000513]

def generate(solver,variant,force=False):
    case=ROOT/'cases'/solver/variant
    if case.exists():
        if not force:raise RuntimeError(f'{case} exists')
        assert case.resolve().is_relative_to((ROOT/'cases').resolve());shutil.rmtree(case)
    source=ROOT.parent/'cases'/solver/'legacyAbsent'
    if not source.exists():species.generate(solver,'legacyAbsent')
    for part in ['0','constant','system']:shutil.copytree(source/part,case/part)
    single=solver==SOLVERS[0];zero=case/'0' if single else case/'0/fluid';constant=case/'constant' if single else case/'constant/fluid';system=case/'system' if single else case/'system/fluid'
    # Carrier fields, chemistry and energy stay frozen; only standard cloud evolves.
    flags=(system/'fvSolution').read_text().replace('solveSpecies true;','solveSpecies false;').replace('solveParcels false;','solveParcels true;')
    (system/'fvSolution').write_text(flags)
    charge=-1 if variant=='negative' else 0 if variant=='neutral' else 1
    electric=-50000 if variant=='reverseField' else 50000
    p=zero/'E';text=p.read_text();text=re.sub(r'internalField\s+uniform\s+[^;]+;',f'internalField uniform ({electric} 0 0);',text);p.write_text(text)
    if variant=='missingE':p.unlink()
    cloud=(ROOT.parent/'templates/reactingCloud1Properties').read_text()
    cloud=re.sub(r'U0\s*\([^;]*\);','U0 (0 0 0);',cloud)
    cloud=re.sub(r'parcelCharge\s+[^;]+;',f'parcelCharge {charge};',cloud)
    cloud=cloud.replace('T0                  293.15;','T0                  298.15;')
    if variant=='noCoulomb':cloud=re.sub(r'CoulombSep\s*\{[^{}]*\}','',cloud)
    a=cloud.index('    cloudFunctions');cloud=cloud[:a]+'    cloudFunctions {}\n'
    (constant/'reactingCloud1Properties').write_text(cloud)
    species.write(constant/'particlePositions','1\n((.01013 .01 .000513))\n','vectorField')
    species.write(case/'system/controlDict',f'''application {solver};
startFrom startTime; startTime 0; stopAt endTime; endTime {END}; deltaT {DT};
writeControl timeStep; writeInterval 250; purgeWrite 0;
writeFormat ascii; writePrecision 16; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
adjustTimeStep no; maxCo .5; maxDi 10;
''')
    cc=1+2*LAMBDA/D*(1.257+.4*math.exp(-1.1*D/(2*LAMBDA)))
    beta=3*math.pi*MU*D/cc;mass=RHO*math.pi*D**3/6;tau=mass/beta
    terminal=0 if variant in ['neutral','noCoulomb'] else charge*Q*electric/beta
    expected=dict(solver=solver,variant=variant,singleRegion=single,initialPosition=START,diameter=D,rhoP=RHO,mu=MU,lambdaMeanFreePath=LAMBDA,Cunningham=cc,tau=tau,chargeNumber=charge,E=electric,terminalVelocity=terminal,endTime=END,dt=DT,expectedDisplacement=terminal*(END-tau*(1-math.exp(-END/tau))),expectedFailure=variant=='missingE')
    (case/'expected.json').write_text(json.dumps(expected,indent=2)+'\n');return case

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--solver',choices=SOLVERS);ap.add_argument('--variant',choices=VARIANTS);ap.add_argument('--force',action='store_true');a=ap.parse_args()
    for s in [a.solver] if a.solver else SOLVERS:
        for v in [a.variant] if a.variant else VARIANTS:
            if v=='missingE' and s==SOLVERS[-1]:continue
            print(generate(s,v,a.force))
if __name__=='__main__':main()