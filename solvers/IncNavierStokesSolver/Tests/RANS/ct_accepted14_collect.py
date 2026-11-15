#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT14 checkpoint evidence from sealed primitives and actual histories."""
import argparse,json,os,resource,time
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
from ct_accepted9_pressure_run import sealed
from ct_accepted11_target import target_pair
from ct_accepted_budget_gate import physical_budgets
from ct_accepted14_pair import numerical_receipts
from ct_accepted_open import require
from compare_cls_epawf_pimple import same_history
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def read(p):return json.loads(Path(p).read_text())
def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','capacity','component','runtime','output','self-control'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--static',type=Path);p.add_argument('--pair',type=Path)
    a=p.parse_args();started=time.monotonic();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    cap=sealed(a.capacity);component=sealed(a.component);runtime=runtime_identity(a.runtime/'IncNavierStokesSolver');probe=runtime_identity(a.runtime/'RANSPreservationProbe');require(cap['status']==component['status']=='pass','unqualified capacity/component')
    m=RequestManifest(a.output,dict(schema='CT14CheckpointEvidenceV1',source=sha(Path(__file__)),capacity=sha(a.capacity/'results.json'),component=sha(a.component/'results.json'),static=sha(a.static/'results.json') if a.static else None,pair=sha(a.pair/'results.json') if a.pair else None,runtime=runtime,probe=probe))
    require(set(cap['results'])=={'original-initial','original-evolved','candidate-alpha1000','candidate-alpha2000'},'capacity system coverage')
    fills={}
    for label,r in cap['results'].items():
        require(sha(Path(r['matrix']['file']))==r['matrix']['sha256'] and r['execution']['exit_code']==0 and 0<r['factor_nnz']<=40000000 and r['execution']['peak_address_space_KiB']<=16*1024**2,'actual capacity primitive')
        for prefix in ('','tiny_','cancelling_'):require(0<=r[prefix+'true_residual']<=2e-12*r[prefix+'rhs_maximum']+1e-25,'actual frozen true residual')
        fills[label]=r
    require(set(component['small'])=={'d0-p3','d0-p4','displaced-p3','displaced-p4'},'component fixture coverage')
    action=max(v for r in component['small'].values() for v in r['complete_independent_actions']['defects'].values());require(action<=2e-10,'component independent actions')
    for f,digest in component['histories']['reports'].items():require(sha(Path(f))==digest,'native history report changed')
    full,continued,one,rollback=[read(a.component/x/'case.evolution.json') for x in ('full','continued','one','rollback')];require(same_history(full,continued) and same_history(one,rollback),'full native restart/rollback history')
    update=max(x['maximum_scaled_field'] for x in component['original']['common_times']);traction=max(x['maximum_scaled_traction'] for x in component['original']['common_times']);require(max(update,traction)<=1e-8,'original full-domain update equivalence')
    out=dict(schema='CT14CheckpointEvidenceV1',component_identity='CT14-VectorCOLAMDV1',runtime=runtime,probe=probe,capacity=fills,component=dict(status='pass',maximum_action_defect=action,maximum_scaled_original_field_difference=update,maximum_scaled_original_traction_difference=traction,history=component['histories'],small=component['small'],original=component['original']),static=dict(status='UNATTEMPTED'),pair=dict(status='UNATTEMPTED'),closure='IN_PROGRESS')
    if a.static:
        s=sealed(a.static);sr=read(a.static/'request-manifest.json')['request'];require(s['status']=='pass' and sr['component']==sha(a.component/'results.json') and sr['capacity']==sha(a.capacity/'results.json') and sr['runtime']==runtime and sr['probe']==probe,'static exact-runtime prerequisite binding')
        require(set(s['controls'])=={'native','alpha2000','q48','q64'} and set(s['physical_scalar_controls'])=={'native','q48f16','q64f16','q64f32'} and set(s['canonical_exports'])=={'1000','2000'},'complete static coverage')
        for label,r in s['controls'].items():
            require(r['execution']['exit_code']==0 and r['predictor_factor_nnz']<=40000000 and r['pressure_C_nnz']==0,'static flow/capacity')
            if label in ('q48','q64'):require(r['pressure_change_scaled_maximum']<1e-8,'static governing quadrature')
        for r in s['physical_scalar_controls'].values():require(r['execution']['exit_code']==0 and r['execution']['peak_virtual_KiB']<=16*1024**2 and r['classification']['reconstruction_status']=='pass','static scalar capacity/reconstruction')
        out['static']=s
    if a.pair:
        require(a.static is not None,'physical pair needs its admitted static receipt');pair=sealed(a.pair);pr=read(a.pair/'request-manifest.json')['request'];require(pr['static']==sha(a.static/'results.json') and pr['component']==sha(a.component/'results.json') and pr['capacity']==sha(a.capacity/'results.json') and pr['runtime']==runtime and pr['probe']==probe,'pair exact-runtime prerequisite binding')
        for name,r in pair['members'].items():
            folder=Path(r['directory']);native=read(folder/'case.evolution.json');require(sha(folder/'case.evolution.json')==r['native_report_sha256'] and exact_vector_hash(native['state'])==r['state_sha256'] and exact_vector_hash(native['accepted_states'][0])==r['actual_initial_sha256'],'physical member actual state/report identity')
            require(numerical_receipts(folder/'sparse-receipts.jsonl')['status']=='pass','actual evolution linear receipts')
            require(physical_budgets(r['budgets'])['status']=='pass','actual consecutive-state budgets')
        decision=target_pair(pair);require(decision==pair['qualification'],'physical pass label differs from unchanged classifier')
        out['pair']=pair;out['closure']='PASS' if decision['status']=='pass' else 'UNRESOLVED'
    resources={};cpu=wall=0.;project=aggregate=0
    for folder in sorted(a.root.glob('resources-*')):
        if folder==a.self_control:continue
        f=folder/'execution.json';require(f.is_file(),'unclosed owned family:'+str(folder));e=read(f);r=e['resource'];require(r['project_status']==r['aggregate_status']=='pass','worker violation:'+str(folder))
        cpu+=e['process_cpu_seconds'];wall+=e['wall_seconds']-e['allocation_wait_seconds'];project=max(project,r['maximum_project']);aggregate=max(aggregate,r['maximum_aggregate']);resources[folder.name]=dict(execution=e,sha256=sha(f))
    # Whole elapsed package time is a conservative bound including Mac source,
    # transfer and inspection, even when no Linux family was active.
    whole_elapsed=time.time()-1794751358;output=sum(p.lstat().st_size for p in a.root.rglob('*') if p.is_file() or p.is_symlink());wall_upper=max(whole_elapsed,wall+time.monotonic()-started+240);cpu_upper=cpu+time.process_time()+120;output_upper=output+64*1024**2
    require(cpu_upper<=24000 and wall_upper<=28800 and output_upper<=16*1024**3,'finite CT14 package allocation')
    out['resources']=dict(closed_families=resources,closed_CPU_seconds=cpu,CPU_upper_with_collector_publication_reserve=cpu_upper,active_wall_upper=wall_upper,whole_elapsed_package_seconds=whole_elapsed,retained_new_bytes=output,new_bytes_upper_with_publication_reserve=output_upper,maximum_owned_family_project=project,maximum_aggregate=aggregate,reserve='120CPU seconds and64MiB; whole elapsed wall includes transport/source/inspection; add closed collector/publication receipts at final delivery')
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print(json.dumps(dict(stage=out['closure'],capacity='pass',component='pass',static=out['static']['status'],pair=out['pair'].get('qualification',out['pair']).get('status'),closed_CPU=cpu,wall_upper=wall_upper,new_bytes=output)),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
