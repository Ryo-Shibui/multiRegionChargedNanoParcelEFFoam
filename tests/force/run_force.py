#!/usr/bin/env python3
"""Run upgraded native nanoparticle force cases after builds are ready."""
from pathlib import Path
import argparse,subprocess,shutil,json,re,datetime,traceback
import make_force_cases as gen
import check_force as check
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--solver',choices=gen.SOLVERS);ap.add_argument('--variants',nargs='+',choices=gen.VARIANTS);a=ap.parse_args()
    results=[];evidence=ROOT/'evidence';evidence.mkdir(exist_ok=True)
    for solver in [a.solver] if a.solver else gen.SOLVERS:
        for variant in a.variants or gen.VARIANTS:
            if variant=='missingE' and solver==gen.SOLVERS[-1]:continue
            case=ROOT/'cases'/solver/variant
            if not case.exists():gen.generate(solver,variant)
            old=[p for p in case.iterdir() if (p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name) and float(p.name)>0) or p.name=='postProcessing' or p.name.startswith('log.') or p.name in ['acceptance.json','.completed']]
            if old:
                dest=evidence/'prior_runs'/solver/variant/datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f');dest.mkdir(parents=True)
                for p in old:shutil.move(str(p),dest/p.name)
            ctrl=case/'system/controlDict';original=ctrl.read_text()
            original=re.sub(r'startFrom\s+[^;]+;','startFrom startTime;',original)
            original=re.sub(r'endTime\s+[^;]+;',f'endTime {gen.END};',original)
            ctrl.write_text(original)
            executable=shutil.which(solver) or solver
            def run(name):
                with (case/name).open('w') as f:return subprocess.run([executable,'-case',str(case)],stdout=f,stderr=subprocess.STDOUT).returncode
            if variant=='restart':
                ctrl.write_text(original.replace(f'endTime {gen.END};','endTime .00025;'))
                code=run('log.firstHalf')
                if code==0:
                    ctrl.write_text(original.replace('startFrom startTime;','startFrom latestTime;'))
                    code=run('log.restart')
            else:code=run('log.solver')
            logs='\n'.join(p.read_text() for p in case.glob('log.*'))
            if variant=='missingE':
                passed=code!=0 and ('cannot find file' in logs or 'not found' in logs) and ('/E' in logs or ' E ' in logs)
                r=dict(solver=solver,variant=variant,passed=passed,exitCode=code,failures=[] if passed else ['Missing prescribed E should fail with a clear field diagnostic'])
            elif code:r=dict(solver=solver,variant=variant,passed=False,exitCode=code,failures=['Solver process failed'],logTail=logs[-3500:])
            else:
                try:r=check.check(case);r['exitCode']=code
                except Exception as error:r=dict(solver=solver,variant=variant,passed=False,exitCode=code,failures=[str(error)],traceback=traceback.format_exc())
            if r['passed']:(case/'.completed').touch()
            (case/'acceptance.json').write_text(json.dumps(r,indent=2)+'\n')
            dest=evidence/solver/variant;dest.mkdir(parents=True,exist_ok=True)
            for p in list(case.glob('log.*'))+[case/'expected.json',case/'acceptance.json']:shutil.copy2(p,dest/p.name)
            results.append(r);print(solver,variant,'PASS' if r['passed'] else 'FAIL',r.get('failures',[]),flush=True)
    # Sign reversal is checked using independent force computations in separate native cases.
    for solver in [a.solver] if a.solver else gen.SOLVERS:
        subset={r['variant']:r for r in results if r['solver']==solver and r['passed'] and 'displacement' in r}
        if 'positive' in subset and 'negative' in subset:
            pos=subset['positive']['displacement'][0];neg=subset['negative']['displacement'][0]
            error=abs(pos+neg)/max(abs(pos),1e-30)
            r=dict(solver=solver,variant='signReversalComparison',passed=error<1e-8,relativeAntisymmetryError=error,failures=[] if error<1e-8 else ['Opposite charges did not reverse drift symmetrically'])
            results.append(r);print(solver,'signReversalComparison','PASS' if r['passed'] else 'FAIL',flush=True)
    for item in results:
        dest=evidence/item['solver']/item['variant'];dest.mkdir(parents=True,exist_ok=True)
        (dest/'acceptance.json').write_text(json.dumps(item,indent=2)+'\n')
    combined=[json.loads(p.read_text()) for p in evidence.glob('*/*/acceptance.json')]
    (evidence/'results.json').write_text(json.dumps(combined,indent=2)+'\n')
    raise SystemExit(0 if all(r['passed'] for r in results) else 1)
if __name__=='__main__':main()