#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fixed COLAMD/AMD/Natural assessment using the actual guarded native factor."""
import argparse,os,re,resource
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('exports','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    data=sealed(a.exports);m=RequestManifest(a.output,dict(schema='CT13FixedOrderingAssessmentV1',source=sha(Path(__file__)),exports=sha(a.exports/'results.json'),contract=sha(a.contract),runtime=runtime_identity(a.runtime/'CT13SparseOrderingProbe'),orderings=['COLAMD','AMD','Natural']))
    results={};fill={k:[] for k in ('COLAMD','AMD','Natural')}
    for label,s in data['systems'].items():
        require(sha(Path(s['A']['file']))==s['A']['sha256'],'foreign canonical matrix')
        rhs=a.output/(label+'.rhs');rhs.write_text('\n'.join(format(x,'.17g') for x in s['rhs'])+'\n');results[label]={};control=None
        for ordering in fill:
            out=a.output/(label+'-'+ordering+'.json')
            e=run([a.runtime/'CT13SparseOrderingProbe',s['A']['file'],rhs,ordering,out],a.output,label+'-'+ordering+'.log')
            row=__import__('json').loads(out.read_text()) if out.is_file() else dict(status='native_failure',failure='no factor result; inspect native log')
            compact={k:v for k,v in row.items() if k not in ('solution','tiny_solution')};compact.update(execution=e,result_sha256=sha(out) if out.exists() else None,matrix_sha256=s['A']['sha256'])
            if e['exit_code']==0:
                require(row['factor_nnz']<=40000000 and row['true_residual']<=2e-12*row['rhs_maximum']+1e-25 and row['tiny_true_residual']<=2e-12*row['tiny_rhs_maximum']+1e-25,'pass label violates original factor/accuracy guards')
                compact['fill_margin']=40000000-row['factor_nnz'];fill[ordering].append(row['factor_nnz'])
                if ordering=='COLAMD':control=np.asarray(row['solution'])
                elif control is not None:
                    compact['maximum_solution_difference']=float(np.max(abs(np.asarray(row['solution'])-control),initial=0.))
                    require(compact['maximum_solution_difference']<=1e-8,'predictor solution disagreement')
            else:
                fill[ordering].append(None);match=re.search(r'factor_nnz=(\d+)',row.get('failure',''))
                if match:compact['observed_factor_nnz']=int(match[1])
            results[label][ordering]=compact;write_json(a.output/'progress.json',dict(results=results,fill=fill))
            print(label,ordering,compact['status'],compact.get('factor_nnz',compact.get('observed_factor_nnz')),flush=True)
    eligible=[k for k in ('AMD','Natural') if len(fill[k])==len(data['systems']) and all(v is not None and v<=40000000 for v in fill[k])]
    winner=min(eligible,key=lambda k:max(fill[k])) if eligible else None
    result=dict(schema='CT13FixedOrderingAssessmentV1',status='ordering_selected' if winner else 'blocked_capacity',selected=winner,results=results,fill=fill,selection='one fixed alternative with smallest maximum fill across all prospectively declared systems; original guards/residual/refinement unchanged',governing_component='not_enabled_by_read_only_assessment')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print(result['status'],winner,flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
