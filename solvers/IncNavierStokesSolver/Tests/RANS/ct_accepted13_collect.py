#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate CT13 capacity, reconstruction, memory and physical admission."""
import argparse,copy,gc,json,math,os,subprocess,resource
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha

SYSTEMS={'original-initial','original-evolved','candidate-alpha1000','candidate-alpha2000'}
def finite(x):return type(x) in (int,float) and math.isfinite(x)
def factor_gate(row,matrix_sha):
    reasons=[]
    if row.get('matrix_sha256')!=matrix_sha:reasons.append('matrix_identity')
    if row.get('execution',{}).get('exit_code')!=0:reasons.append('native_completion')
    n=row.get('factor_nnz')
    if not finite(n) or not 0<n<=40000000:reasons.append('factor_fill')
    for prefix in ('','tiny_'):
        rhs=row.get(prefix+'rhs_maximum');r=row.get(prefix+'true_residual')
        if not finite(rhs) or rhs<0 or not finite(r) or not 0<=r<=2e-12*rhs+1e-25:reasons.append(prefix+'true_residual')
    return sorted(set(reasons))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','exports','orderings','memory','image','ct11','ct12','output','self-control'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));ex=sealed(a.exports);order=sealed(a.orderings);memory=sealed(a.memory);image=sealed(a.image)
    require(ex['schema']=='CT13PreFactorExportQualificationV2' and set(ex['controls'])==SYSTEMS,'missing/foreign no-factor controls')
    for name,r in ex['controls'].items():
        require(finite(r.get('whole_command_sparse_factorizations')) and r['whole_command_sparse_factorizations']==0,'startup factorization not excluded')
        require(r['execution']['exit_code']==0,'incomplete no-factor export')
        if name.startswith('original'):require(set(r.get('defects',{}))=={'system','mass','G','D','C','residual'},'missing no-factor reconstruction primitive')
        require(all(finite(v) and 0<=v<=2e-10 for v in r.get('defects',{}).values()),'invalid no-factor reconstruction primitive')
        if name.startswith('candidate'):require(r['slots']==62533 and r['exact_prototype_CSR_bytes'],'foreign selected mesh/coefficients')
        else:require(r['actual_history_differs_from_state'],'evolved history substituted with state')
    require(set(order['results'])==SYSTEMS,'ordering system coverage')
    policy={}
    for label in ('COLAMD','AMD','Natural'):
        policy[label]={name:factor_gate(order['results'][name].get(label,{}),ex['systems'][name]['A']['sha256']) for name in sorted(SYSTEMS)}
    eligible=[label for label in ('AMD','Natural') if all(not v for v in policy[label].values())]
    selected=min(eligible,key=lambda label:max(order['results'][name][label]['factor_nnz'] for name in SYSTEMS)) if eligible else None
    require(order['selected']==selected,'ordering pass label/selection differs from complete primitives')
    require(memory['observer_policy']=='NativeAuditReceiptV1' and set(memory['controls'])=={'48','64'},'missing/foreign memory policy/control')
    native_inventories={}
    for label,r in memory['controls'].items():
        require(r['execution']['exit_code']==0,'diagnostic capacity still blocked')
        if label=='48':require(r['mandatory_diagnostics_exact'] and r['all_observation_primitive_defect']<=2e-10 and max(r['all_row_defects'].values())<=2e-10,'supported observer equivalence failed')
        else:require(r['all_row_reconstruction']<=2e-10 and r['pressure_image']['maximum_scaled_value']<1e-8 and r['tighter_repeat']['maximum_scaled_value']<1e-8,'q64 action/pressure uncertainty failed')
        path=a.memory/('q'+label)/'loads.json';require(sha(path)==r['packet_sha256'],'foreign held observer state/packet')
        q=json.loads(path.read_text());audit=q['observation']['audit_alpha_zero']
        for key in ('accepted_awf_scalar_domain','accepted_awf_hard_admissible'):
            require(key in audit and len(audit[key])>0 and finite(audit[key][0]),'missing/nonfinite native observer domain flag')
        physical=all(audit[key][0]==1 for key in ('accepted_awf_scalar_domain','accepted_awf_hard_admissible'))
        require(q['physical_admissibility']==('pass' if physical else 'fail'),'observer domain label differs from actual audit')
        native_inventories[label]=dict(status='pass' if physical else 'fail',packet_sha256=r['packet_sha256'],state_sha256=r['state_sha256'],resolved_descriptor=q['observation']['resolved_descriptor'],numerical_operator=q['numerical_operator'],actual_quadrature=q['observation']['actual_quadrature'],audit=audit)
        del q;gc.collect()
    require(memory['controls']['48']['state_sha256']==memory['controls']['64']['state_sha256'],'foreign held state')
    require(set(memory['small_controls'])=={'d0','displaced','unknown_policy'},'missing small observer/operator controls')
    for k in ('d0','displaced'):require(max(memory['small_controls'][k]['defects'].values())<=2e-10,'small all-row/observation mismatch')
    require('unknown probe observer lifecycle' in memory['small_controls']['unknown_policy']['failure'],'wrong negative reason')
    require(set(image['cases'])=={'baseline','selected'},'missing initial image comparison')
    for r in image['cases'].values():
        require(r['row_reconstruction']<=2e-10 and abs(sum(r['pressure_volumes'])-1100)<2e-10,'initial load/mass/domain disagreement')
        for v in r['measurements'].values():require(finite(v['repeat_physical_max_upper']) and v['repeat_physical_max_upper']<1e-8,'initial pressure-image numerical uncertainty')
    base=image['cases']['baseline']['measurements']['inner_shear']['pressure']['whole']['weighted_RMS'];fine=image['cases']['selected']['measurements']['inner_shear']['pressure']['whole']['weighted_RMS']
    require(finite(base) and finite(fine) and base>1e-8 and fine>=0,'initial shear signal unresolved/nonfinite')
    ratio=fine/base;require(abs(ratio-image['selected_to_baseline_shear_RMS_ratio'])<1e-12,'image ratio not supported by primitives')
    uncertainty=max(1e-8,10*sum(r['measurements']['inner_shear']['repeat_physical_max_upper'] for r in image['cases'].values()))
    require(image['resolved_initial_image_reduction']==(base-fine>uncertainty),'image resolution label differs from primitives')
    good=copy.deepcopy(order['results']['original-initial']['COLAMD']);fingerprint=ex['systems']['original-initial']['A']['sha256'];negatives={}
    for name,key,value,word in [('pass_bad_fill','factor_nnz',40000001,'factor_fill'),('pass_bad_residual','true_residual',1.,'true_residual'),('pass_foreign_matrix','matrix_sha256','0'*64,'matrix_identity'),('pass_nan_rhs','rhs_maximum',float('nan'),'true_residual')]:
        r=copy.deepcopy(good);r.update(status='pass');r[key]=value;why=factor_gate(r,fingerprint);require(word in why,'nondiscriminating primitive negative:'+name);negatives[name]=dict(status='pass',recomputed_rejection=why)
    resources={};cpu=wall=0.;aggregate=0;roots=set()
    for d in sorted(a.root.glob('resources-*')):
        if d==a.self_control:continue
        require((d/'execution.json').is_file(),'unclosed family:'+d.name);e=json.loads((d/'execution.json').read_text());r=e['resource']
        require(r['project_status']==r['aggregate_status']=='pass','resource violation:'+d.name)
        hw=d/'worker-guard-highwater.jsonl';rows=[json.loads(x) for x in hw.read_text().splitlines()] if hw.exists() else []
        ev=d/'worker-guard-events.jsonl'
        resources[d.name]=dict(execution=e,sha256=sha(d/'execution.json'),highwater=rows,events=[json.loads(x) for x in ev.read_text().splitlines()] if ev.exists() else [])
        roots.update(e[k] for k in ('control_parent_pid','guard_pid','owned_family_pid'));cpu+=e['process_cpu_seconds'];wall+=max(0.,e['wall_seconds']-e['allocation_wait_seconds']);aggregate=max(aggregate,r['maximum_aggregate'])
    peak=0;witness=None
    for row in [r for v in resources.values() for r in v['highwater']]:
        owned=set(roots)
        while True:
            add={p['pid'] for p in row['processes'] if p['parent'] in owned}-owned
            if not add:break
            owned.update(add)
        threads=sum(p['threads'] for p in row['processes'] if p['pid'] in owned)
        if threads>peak:peak=threads;witness=dict(time_unix=row['time_unix'],threads=threads,processes=[p for p in row['processes'] if p['pid'] in owned])
    output=sum(f.stat().st_size for f in a.root.rglob('*') if f.is_file() and f.relative_to(a.root).parts[0] not in ('source','tmp'))
    require(peak<=16 and aggregate<=72 and cpu<24000 and wall<28800 and output<16*1024**3,'closed resource/package bound')
    old=sealed(a.ct11/'target-v1');require(old['qualification']['reasons']==['temporal:Pi'],'historical pressure classification changed')
    previous=json.loads((a.ct12/'publication-v1/prepublication-accounting.json').read_text())
    m=RequestManifest(a.output,dict(schema='CT13CompactDecisionV1',source=sha(Path(__file__)),inputs={k:dict(path=str(v),result_sha256=sha(v/'results.json'),seal_sha256=sha(v/'completed-artifacts.json')) for k,v in [('exports',a.exports),('orderings',a.orderings),('memory',a.memory),('initial_image',a.image)]},resource_receipts={k:v['sha256'] for k,v in resources.items()}))
    decision=dict(schema='CT13RequiredScopeV1',pre_factor_export=dict(status='pass',whole_command_global_factors=0,capacity_only=True),ordering=dict(status='qualified_candidate' if selected else 'no_qualified_policy',selected=selected,recomputed_reasons=policy),component_equivalence=dict(status='not_attempted',reason='no qualified ordering' if not selected else 'separate component qualification required'),diagnostic_memory=dict(status='pass',policy='NativeAuditReceiptV1',governing_change=False),initial_image=dict(status='diagnostic_only',selected_to_baseline_shear_RMS_ratio=ratio,reduction=fine<base,not_evolved_uncertainty=True),selected_static_admission=dict(status='blocked',first_gate='predictor_capacity_or_work',missing='component equivalence,complete static flow/scalar/alpha/quadrature and new-component history admission'),short_physical_pair=dict(status='not_attempted_or_admitted',candidate_horizon=0.,original_target_horizon=.002),historical_CT11=old['qualification'],campaign=dict(status='no_go_unchanged',governing_speedup=False,headroom=.25,existing_CPU_ceiling=24000,existing_active_wall_ceiling=28800,existing_output_ceiling_GiB=24),resources=dict(status='pass',closed_CPU_before_collection=cpu,active_wall_upper_before_collection=wall,output_bytes_before_collection=output,maximum_project=peak,maximum_aggregate=aggregate,witness=witness),terrain='unattempted/out_of_scope; mapped numbering,joined-map gradients/cofactors/transmission and historical curved gates remain',physical_applicability='not precursor or hill validation')
    bodies=dict(exports=ex,orderings=order,memory=memory,initial_image=image)
    decision['diagnostic_memory']['held_native_domain']={k:v['status'] for k,v in native_inventories.items()}
    decision['initial_image'].update(resolved_reduction=image['resolved_initial_image_reduction'],comparison_uncertainty_floor=uncertainty)
    # Heavy RHS/solutions/CSR/states stay in the sealed archive; retain their
    # hashes and every decisive numerical/classification primitive here.
    for s in bodies['exports']['systems'].values():s.pop('rhs',None)
    for name,v in bodies.items():write_json(a.output/(name+'.json'),v)
    write_json(a.output/'decision.json',decision);write_json(a.output/'resources.json',resources);write_json(a.output/'negative-controls.json',negatives)
    write_json(a.output/'native-observer-inventories.json',native_inventories)
    write_json(a.output/'source.json',dict(revision=subprocess.check_output(['git','-c','core.preloadIndex=false','rev-parse','HEAD'],cwd=a.root/'source',text=True).strip(),historical_root=str(a.ct12),previous_closed_accounting_sha256=sha(a.ct12/'publication-v1/prepublication-accounting.json'),physical_manifest_sha256=sha(a.root/'source/docs/plans/cao-tamura-accepted-8-evidence/physical-reference-manifest.json'),governing='unchanged pinned CT9 build-v8'))
    files=[f.name for f in a.output.iterdir() if f.is_file() and f.name not in ('request-manifest.json','completed-artifacts.json')];m.seal(a.output,files);print(json.dumps(decision,indent=2),flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
