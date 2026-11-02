#!/usr/bin/env python3
"""Self-contained three-solver dilute-tracer acceptance case generator."""
from pathlib import Path
import argparse,shutil,json,math,re
ROOT=Path(__file__).resolve().parent
SOLVERS=['reactingChargedNanoParcelFoam','chtReactingNanoParcelFoam','multiRegionChargedNanoParcelEFFoam']
VARIANTS=['dilute','temperaturePressureLaw','legacyAbsent','legacyExplicit','missingD','missingSpecies','negativeD','steadyField']
RHO=1.204;MU=1.81e-5;L=.02;NX=20;NY=100;DT=.001;END=.2
TRACERS={'DEG':dict(mean=.0001,amplitude=.00002,D=7.24e-6,molWeight=106.12),'H2O':dict(mean=.002,amplitude=.0005,D=2e-5,molWeight=18.0153)}

def header(obj,cls='dictionary'):return f'FoamFile {{version 2.0; format ascii; class {cls}; object {obj};}}\n'
def write(path,text,cls='dictionary'):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(header(path.name,cls)+text)
def fields(path,name,values):
    b=''.join(f'{p} {{ type '+('symmetryPlane' if p in ['front','back'] else 'zeroGradient')+'; }\n' for p in ['inlet','outlet','lowerElectrode','upperElectrode','front','back'])
    write(path/name,'dimensions [0 0 0 0 0 0 0];\ninternalField nonuniform List<scalar> '+str(NX*NY)+' (\n'+'\n'.join(f'{v:.17g}' for v in values)+'\n);\nboundaryField {\n'+b+'}\n','volScalarField')
def generate(solver,variant,force=False):
    case=ROOT/'cases'/solver/variant
    if case.exists():
        if not force:raise RuntimeError(f'{case} exists; use --force to replace this generated case')
        assert case.resolve().is_relative_to((ROOT/'cases').resolve());shutil.rmtree(case)
    single=solver=='reactingChargedNanoParcelFoam'
    for part in ['0','constant','system']:
        target=case/part if single else case/part/'fluid'
        shutil.copytree(ROOT/'templates/baseFluid'/part,target)
    zero=case/'0' if single else case/'0/fluid'
    constant=case/'constant' if single else case/'constant/fluid'
    system=case/'system' if single else case/'system/fluid'
    (zero/'O2').unlink(missing_ok=True)
    vals={s:[p['mean']+p['amplitude']*math.cos(math.pi*(j+.5)/NY) for j in range(NY) for i in range(NX)] for s,p in TRACERS.items()}
    vals['N2']=[1-a-b for a,b in zip(vals['DEG'],vals['H2O'])]
    for name,v in vals.items():fields(zero,name,v)
    temperature=400 if variant=='temperaturePressureLaw' else 298.15
    pressure=202650 if variant=='temperaturePressureLaw' else 101325
    for name,value in [('T',temperature),('p',pressure),('p_rgh',pressure)]:
        p=zero/name;text=p.read_text()
        text=re.sub(r'internalField\s+uniform\s+[^;]+;',f'internalField uniform {value};',text)
        text=text.replace('value uniform 101325;',f'value uniform {pressure};')
        p.write_text(text)
    thermo=(constant/'thermophysicalProperties').read_text()
    if single:thermo=thermo.replace('<constant>/fluid/','<constant>/')
    (constant/'thermophysicalProperties').write_text(thermo)
    (constant/'foam.inp').write_text('species (N2 DEG H2O);\nreactions {}\n')
    write(constant/'foam.dat','\n'.join(f'''{name}
{{
    specie {{ molWeight {mw}; }}
    equationOfState {{ rho {RHO}; }}
    thermodynamics {{ Cp 1005; Hf 0; }}
    transport {{ mu {MU}; Pr .6996346153846154; }}
}}''' for name,mw in [('N2',28.0134),('DEG',106.12),('H2O',18.0153)]))
    diffusion='model dilute;\nturbulentSchmidt .9;\nspecies\n{\n'
    for name,p in TRACERS.items():
        if variant=='missingSpecies' and name=='H2O':continue
        if variant=='missingD' and name=='DEG':body=''
        else:body=f'D {-p["D"] if variant=="negativeD" and name=="DEG" else p["D"]};'
        if variant=='temperaturePressureLaw':body+=' Tref 298.15; pref 101325; temperatureExponent 1.75; pressureExponent 1;'
        diffusion+=f'    {name} {{ {body} }}\n'
    diffusion+='}\n'
    if variant=='legacyExplicit':diffusion='model legacy;\n'
    if variant!='legacyAbsent':write(constant/'speciesDiffusionProperties',diffusion)
    else:(constant/'speciesDiffusionProperties').unlink(missing_ok=True)
    # Standard single-region implementation always constructs its cloud. Disable it explicitly.
    if single:
        cloud=(ROOT/'templates/reactingCloud1Properties').read_text().replace('active          true;','active          false;')
        cloud=re.sub(r'particleForces\s*\{.*?\}\s*\}', 'particleForces {}',cloud, count=1,flags=re.S)
        a=cloud.index('        injectionModels');b=cloud.index('        dispersionModel',a)
        cloud=cloud[:a]+'        injectionModels {}\n\n'+cloud[b:]
        a=cloud.index('    cloudFunctions');cloud=cloud[:a]+'    cloudFunctions {}\n'
        (constant/'reactingCloud1Properties').write_text(cloud)
        write(constant/'dynamicMeshDict','dynamicFvMesh staticFvMesh;\n')
    else:
        write(case/'constant/regionProperties','regions (fluid (fluid) solid ());\n')
        write(case/'system/fvSolution','PIMPLE {nOuterCorrectors 1;}\ncoupledEnergyField false;\n')
    write(case/'constant/g','dimensions [0 1 -2 0 0 0 0]; value (0 0 0);','uniformDimensionedVectorField')
    flags=(system/'fvSolution').read_text().replace('solveSpecies true;',f'solveSpecies {"false" if variant=="steadyField" else "true"};')
    flags=flags.replace('|N2|O2|H2O|','|N2|DEG|H2O|')
    (system/'fvSolution').write_text(flags)
    schemes=(system/'fvSchemes').read_text()
    schemes=schemes.replace('default Euler;','default steadyState;' if variant=='steadyField' else 'default Euler;')
    (system/'fvSchemes').write_text(schemes)
    end=1 if variant=='steadyField' else .001 if variant in ['missingD','missingSpecies','negativeD'] else END
    dt=1 if variant=='steadyField' else DT
    steps=round(end/dt)
    write(case/'system/controlDict',f'''application {solver};
startFrom startTime; startTime 0; stopAt endTime; endTime {end}; deltaT {dt};
writeControl timeStep; writeInterval {steps}; purgeWrite 0;
writeFormat ascii; writePrecision 16; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
adjustTimeStep no; maxCo .5; maxDi 10;
''')
    expected=dict(solver=solver,variant=variant,singleRegion=single,rho=RHO,mu=MU,length=L,nx=NX,ny=NY,dt=dt,endTime=end,temperature=temperature,pressure=pressure,expectedFailure=variant in ['missingD','missingSpecies','negativeD'],tracers={})
    for name,p in TRACERS.items():
        actualD=MU/RHO if variant in ['legacyAbsent','legacyExplicit'] else p['D']*((temperature/298.15)**1.75*(101325/pressure) if variant=='temperaturePressureLaw' else 1)
        eigen=4*math.sin(math.pi/(2*NY))**2/(L/NY)**2
        expected['tracers'][name]=dict(**p,effectiveD=actualD,continuumAmplitude=p['amplitude']*math.exp(-actualD*(math.pi/L)**2*end),discreteEulerAmplitude=p['amplitude']/(1+dt*actualD*eigen)**steps)
    (case/'expected.json').write_text(json.dumps(expected,indent=2)+'\n')
    return case

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--solver',choices=SOLVERS);ap.add_argument('--variant',choices=VARIANTS);ap.add_argument('--force',action='store_true');a=ap.parse_args()
    for solver in [a.solver] if a.solver else SOLVERS:
        for variant in [a.variant] if a.variant else VARIANTS:print(generate(solver,variant,a.force))
if __name__=='__main__':main()