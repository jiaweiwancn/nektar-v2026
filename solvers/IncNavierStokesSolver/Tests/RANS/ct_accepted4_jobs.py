#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed CT4 independent-case controller; process completion never admits a case."""
import argparse,json,os,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted_feasibility import write_json
from ct_accepted_open import require
from rans_request_manifest import RequestManifest,sha,runtime_identity

CASES=[(m,6.,3) for m in ('A95','A93','E94E95')]+[(m,h,3) for h in (5.5,6.5) for m in ('A95','A93','E94E95')]+[('A95',6.,4)]
def label(m,h,p):return m+'-'+str(h)+'-P'+str(p)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--phase',choices=['initial','evolve'],required=True)
    for k in ('root','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--runtime',type=Path);p.add_argument('--initial-root',type=Path);p.add_argument('--evolution-root',type=Path)
    a=p.parse_args();r=a.root;rt=a.runtime or r/'build-v1/runtime';ir=a.initial_root or r/'initial-v1';er=a.evolution_root or r/'evolution-v1';s=Path(__file__).parent;refs=r/'reference-v1/results.json';contract=r/'source/docs/plans/cao-tamura-accepted-4-evidence/projection-contract.json'
    m=RequestManifest(a.output,dict(schema='CT4FixedIndependentCasesV1',phase=a.phase,cases=CASES,reference_sha256=sha(refs),contract_sha256=sha(contract),
        initial_root=str(ir),evolution_root=str(er),runtime=runtime_identity(rt/'IncNavierStokesSolver'),source_sha256={f:sha(s/f) for f in ('ct_accepted4_jobs.py','ct_accepted4_run.py')}))
    result=dict(diagnostic_execution='incomplete',cases={})
    if a.phase=='evolve':
        control=json.loads((r/'controls-v1/results.json').read_text());require(control['implementation']=='pass','common_implementation_controls_not_admitted')
    for model,height,degree in CASES:
        name=label(model,height,degree);row={};result['cases'][name]=row;initial=ir/name
        targets=[0.] if a.phase=='initial' else [.008,.02] if height==6 and degree==3 else [.008]
        parent=None
        for target in targets:
            if a.phase=='evolve':
                receipt=json.loads((initial/'results.json').read_text());require(receipt['case']['model']==model and receipt['case']['H_mm']==height,'foreign_model_height_initialization')
                if receipt['status']!='pass':row['not_admitted']='initialization gate failed';break
                if parent:
                    q=json.loads((parent/'results.json').read_text())
                    if q['qualification']['status']!='pass':row['continuation_not_admitted']='previous paired milestone failed';break
            output=initial if a.phase=='initial' else er/name/('T'+str(target))
            cmd=[sys.executable,'-B',str(s/'ct_accepted4_run.py'),'--mode','initial' if a.phase=='initial' else 'advance',
                '--runtime',str(rt),'--references',str(refs),'--contract',str(contract),'--model',model,'--height',str(height),'--degree',str(degree),'--output',str(output)]
            if a.phase=='evolve':cmd+=['--initial',str(initial),'--target',str(target)]+(['--parent',str(parent)] if parent else [])
            with (a.output/(name+'-'+str(target)+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
            file=output/'results.json';body=json.loads(file.read_text()) if file.exists() else None
            row[str(target)]=dict(exit_code=code,request=cmd,receipt_sha256=sha(file) if body else None,
                status=body.get('qualification',body).get('status') if body else 'diagnostic_execution_failed',directory=str(output))
            write_json(a.output/'progress.json',result)
            if code:break
            parent=output
    result['diagnostic_execution']='completed_independent_attempts';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT4',a.phase,'attempts complete')
    return 0 if all(v.get('status')=='pass' for row in result['cases'].values() for v in row.values() if isinstance(v,dict)) and not any('not_admitted' in k for row in result['cases'].values() for k in row) else 1
if __name__=='__main__':raise SystemExit(main())
