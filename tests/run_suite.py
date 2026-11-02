#!/usr/bin/env python3
"""Run executable snapshots for RED evidence, or upgraded solvers for acceptance."""
from pathlib import Path
import argparse,subprocess,json,shutil,traceback,re,datetime
import make_cases,check_cases
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--stage',choices=['preupgrade','upgraded'],default='upgraded');ap.add_argument('--solver',choices=make_cases.SOLVERS);ap.add_argument('--variants',nargs='+',choices=make_cases.VARIANTS);a=ap.parse_args()
    results=[]
    variants=a.variants or (['dilute','legacyAbsent','missingD'] if a.stage=='preupgrade' else make_cases.VARIANTS)
    solvers=[a.solver] if a.solver else make_cases.SOLVERS
    evidence=ROOT/'evidence'/a.stage;evidence.mkdir(parents=True,exist_ok=True)
    for solver in solvers:
        executable=ROOT/'evidence/preupgrade_bin'/solver if a.stage=='preupgrade' else Path(shutil.which(solver) or solver)
        for variant in variants:
            case=ROOT/'cases'/solver/variant
            if not case.exists():make_cases.generate(solver,variant)
            # Preserve previous run artifacts before replacement; initial fields stay unchanged.
            old=[p for p in case.iterdir() if (p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name) and float(p.name)>0) or p.name=='postProcessing' or p.name.startswith('log.') or p.name=='acceptance.json']
            if old:
                archive=ROOT/'evidence/prior_runs'/solver/variant/datetime.datetime.now().strftime('%Y%m%d_%H%M%S_%f');archive.mkdir(parents=True)
                for p in old:shutil.move(str(p),archive/p.name)
            log=case/('log.'+a.stage)
            with log.open('w') as f:run=subprocess.run([str(executable),'-case',str(case)],stdout=f,stderr=subprocess.STDOUT)
            meta=json.loads((case/'expected.json').read_text());text=log.read_text()
            if meta['expectedFailure']:
                diagnostic=('FOAM FATAL' in text and ('DEG' in text or 'H2O' in text) and ('speciesDiffusionProperties' in text or 'diffus' in text.lower()))
                result=dict(solver=solver,variant=variant,exitCode=run.returncode,passed=run.returncode!=0 and diagnostic,failures=[] if run.returncode!=0 and diagnostic else ['Explicit dilute invalid/missing D must fail with a species/dictionary diagnostic'])
            elif run.returncode and not (case/f'{meta["endTime"]:g}').exists():
                result=dict(solver=solver,variant=variant,exitCode=run.returncode,passed=False,failures=['Solver exited before numerical acceptance'],logTail=text[-4000:])
            else:
                try:
                    result=check_cases.check(case);result['exitCode']=run.returncode
                    if run.returncode:
                        result['numericalChecksPassed']=result['passed'];result['passed']=False
                        result['failures'].append(f'Solver runtime/shutdown exitCode {run.returncode}')
                        result['logTail']=text[-2000:]
                except Exception as err:result=dict(solver=solver,variant=variant,passed=False,exitCode=run.returncode,failures=[str(err)],traceback=traceback.format_exc())
            dest=evidence/solver/variant;dest.mkdir(parents=True,exist_ok=True)
            shutil.copy2(log,dest/log.name);shutil.copy2(case/'expected.json',dest/'expected.json')
            if not meta['expectedFailure'] and (case/f'{meta["endTime"]:g}').exists():shutil.copytree(case/f'{meta["endTime"]:g}',dest/'output',dirs_exist_ok=True)
            (dest/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n');(case/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
            results.append(result)
            print(solver,variant,'PASS' if result['passed'] else 'FAIL',result.get('failures',[]),flush=True)
    (evidence/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    raise SystemExit(0 if all(r['passed'] for r in results) else 1)
if __name__=='__main__':main()