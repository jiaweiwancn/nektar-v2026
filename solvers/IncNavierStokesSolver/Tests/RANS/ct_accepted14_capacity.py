#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One prospectively fixed Scotch ordering on CT13's exact canonical systems."""
import argparse,json,os,resource
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from ct_accepted13_orderings import native_run
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import run,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct13','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--policy',choices=['ScotchND0V1','ScotchVectorND0V1','ScotchVectorNDKeepV1','VectorCOLAMDV1'],default='ScotchND0V1')
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    ex=sealed(a.ct13/'exports-v2');co=json.loads(a.contract.read_text());old=sealed(a.ct13/'orderings-v4')
    m=RequestManifest(a.output,dict(schema='CT14ScotchCapacityV1',source=sha(Path(__file__)),exports=sha(a.ct13/'exports-v2/results.json'),old_ordering=sha(a.ct13/'orderings-v4/results.json'),contract=sha(a.contract),runtime=runtime_identity(a.runtime/'CT13SparseOrderingProbe'),policy=a.policy,graph='canonical A+A^T union; optional exact native-DOF vector quotient before ND',seed=0,strategy=(dict(builder='Eigen::COLAMDOrdering<int>',knobs='unchanged defaults',directed_quotient=True) if a.policy=='VectorCOLAMDV1' else dict(builder='SCOTCH_stratGraphOrderBuild',flags=9,level=0,balance=.1)),permutation='old_to_new[old]=new; B[p[i],p[j]]=A[i,j]; x[i]=y[p[i]]'))
    algebra=run([a.runtime/'CT9SparseAlgebraChecks'],a.output,'algebra.log');require(algebra['exit_code']==0,'independent small sparse/scotch controls')
    results={}
    for label,s in ex['systems'].items():
        require(sha(Path(s['A']['file']))==s['A']['sha256'],'changed frozen canonical matrix')
        rhs=a.output/(label+'.rhs');rhs.write_text('\n'.join(format(x,'.17g') for x in s['rhs'])+'\n')
        out=a.output/(label+'.json');command=[a.runtime/'CT13SparseOrderingProbe',s['A']['file'],rhs,a.policy,out]
        if a.policy in ('ScotchVectorND0V1','ScotchVectorNDKeepV1','VectorCOLAMDV1'):
            require(sha(Path(s['packet']))==s['packet_sha256'],'changed canonical field map packet')
            frozen=json.loads(Path(s['packet']).read_text())['pre_factor_frozen'];n0=frozen['block_sizes'][0]
            require(frozen['block_sizes'][1]==n0,'velocity fields do not share native IDs')
            groups=a.output/(label+'.groups');groups.write_text('\n'.join(str(i%n0) for i in frozen['rows'][:frozen['velocity_size']])+'\n');command.append(groups)
        e=native_run(command,a.output,label+'.log',co['first_attempt']['per_factor_CPU_soft'],co['first_attempt']['per_factor_wall_seconds'])
        raw=json.loads(out.read_text()) if out.exists() else dict(status='unresolved_work_or_resource_limit')
        row={k:v for k,v in raw.items() if k not in ('solution','tiny_solution','old_to_new')};row.update(execution=e,matrix=s['A'],rhs_sha256=sha(rhs),result_sha256=sha(out) if out.exists() else None)
        if e['exit_code']==0:
            perm=np.asarray(raw['old_to_new']);require(sorted(perm.tolist())==list(range(s['A']['rows'])),'permutation inversion')
            row['permutation_inversion']='pass';row['fill_margin']=40000000-raw['factor_nnz']
            for prefix in ('','tiny_','cancelling_'):require(raw[prefix+'true_residual']<=2e-12*raw[prefix+'rhs_maximum']+1e-25,'original-system residual:'+prefix)
            if label.startswith('original'):
                prior=a.ct13/'orderings-v4'/(label+'-COLAMD.json');require(sha(prior)==old['results'][label]['COLAMD']['result_sha256'],'changed old original solution')
                before=json.loads(prior.read_text());x=np.asarray(raw['solution']);y=np.asarray(before['solution']);row['old_solution_scaled_difference']=float(np.max(abs(x-y),initial=0)/max(1.,np.max(abs(y),initial=0)))
                require(row['old_solution_scaled_difference']<=1e-8,'old admitted predictor equivalence')
            require(raw['factor_nnz']<=40000000 and e['peak_address_space_KiB']<=16*1024**2,'capacity primitive')
        results[label]=row;write_json(a.output/'progress.json',dict(results=results));print(label,row['status'],row.get('factor_nnz'),flush=True)
        if e['exit_code']!=0:break
    passed=len(results)==4 and all(r['execution']['exit_code']==0 for r in results.values())
    result=dict(schema='CT14ScotchCapacityV1',status='pass' if passed else 'unresolved',results=results,algebra=algebra,component_admitted=False)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if passed else 2
if __name__=='__main__':raise SystemExit(main())
