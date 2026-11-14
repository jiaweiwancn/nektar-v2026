#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact CT10 scientific decisions; numerical/resource completion are separate."""
import argparse,copy,hashlib,json,os,subprocess
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair
from ct_accepted9_collect import compact,compact_pair,timing,LIBRARIES
from ct_accepted10_target_run import target_pair
from ct_accepted_feasibility import write_json,require
from rans_request_manifest import RequestManifest,sha,canonical


def diagnosis(q,vectors,scale):
    errors=[]
    for name in ('native_reconstruction','frozen_momentum_agreement'):
        x=q.get(name);errors+=[] if isinstance(x,(int,float)) and np.isfinite(x) and 0<=x<=2e-10 else [name]
    for name in ('signed_component_sum_error','signed_band_sum_error','tighter_repeat_difference'):
        x=q.get(name,{}).get('maximum_scaled_value');errors+=[] if isinstance(x,(int,float)) and np.isfinite(x) and 0<=x<1e-8 else [name]
    cond=q.get('pressure_conditions',{})
    if cond.get('maximum_C')!=0 or cond.get('unused_gauge_in_system') is not False or cond.get('absolute_pressure') is not True:errors.append('pressure_conditions')
    columns=np.array(list(vectors['components'].values()));total=np.array(vectors['total']['delta']);defect=float(np.max(abs(columns.sum(axis=0)-total)))/scale
    if not np.isfinite(defect) or defect>=1e-8:errors.append('recomputed_signed_coefficient_sum')
    for row in q.get('frozen_quadrature',[]):
        if row['pressure_image_change']['maximum_scaled_value']>=1e-8:errors.append('quadrature')
    if len(q.get('frozen_quadrature',[]))!=2:errors.append('missing_quadrature')
    ind=q.get('independent',{})
    if not ind.get('load_reconstruction') or any(v['maximum_defect']>2e-10 for v in ind['load_reconstruction'].values()):errors.append('independent_loads')
    if ind.get('partition_defect',1)>2e-12 or ind.get('constant_gradient_defect',1)>2e-10:errors.append('constant_modes')
    return dict(status='pass' if not errors else 'fail',reasons=errors,recomputed_signed_coefficient_sum_scaled=defect)


def resources(root):
    rows={}
    for file in sorted(root.glob('resources-*/execution.json')):
        q=n.read(file);q['high_water']={p.name:dict(sha256=sha(p),records=[json.loads(s) for s in p.read_text().splitlines()]) for p in file.parent.glob('*highwater.jsonl')};events=file.parent/'worker-guard-events.jsonl';q['suspension_events']=[json.loads(s) for s in events.read_text().splitlines()] if events.exists() else [];rows[file.parent.name]=q
    return rows


def member_profiles(q):
    return {dt:dict(actual_initial=r['actual_initial_traction_profile'],endpoint=r['traction_history'][-1],accepted_history=[dict(time=v['time'],mean_wall_shear=v['profile']['mean_wall_shear'],wall_error=v['profile']['normalized_wall_traction_error'],interface_error=v['profile']['normalized_interface_traction_RMS'],wall_standard_deviation=v['profile']['wall_spatial_standard_deviation'],state_sha256=v['profile']['state_sha256']) for v in r['traction_history']]) for dt,r in q['members'].items()}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct9','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--target',type=Path);a=p.parse_args();root=a.root
    labels=['controls-v1','static-v1','baseline-attribution-v3','candidate-attribution-v1','baseline-T008-v1','candidate-T008-v1','candidate-T020-v2','common-space-v2','target-preflight-v1','build-v2','build-v3']
    bodies={name:sealed(root/name) for name in labels};sources={name:dict(path=str(root/name),sha256=sha(root/name/'results.json'),request=sha(root/name/'request-manifest.json'),seal=sha(root/name/'completed-artifacts.json')) for name in labels}
    auditcap=sealed(root/'target-audit-capacity-v1');sources['required-audit-capacity']=dict(path=str(root/'target-audit-capacity-v1'),sha256=sha(root/'target-audit-capacity-v1/results.json'),request=sha(root/'target-audit-capacity-v1/request-manifest.json'),seal=sha(root/'target-audit-capacity-v1/completed-artifacts.json'))
    target=sealed(a.target) if a.target else None
    if a.target:sources['target-stepping']=dict(path=str(a.target),sha256=sha(a.target/'results.json'),request=sha(a.target/'request-manifest.json'),seal=sha(a.target/'completed-artifacts.json'))
    m=RequestManifest(a.output,dict(schema='CT10FinalEvidenceV1',sources=sources,source=sha(Path(__file__)),scope='primitive decisions; a successful collector does not admit a physical trajectory'))
    pre=bodies['static-v1'];scale=pre['cases']['baseline']['ustar']**2;dg={}
    for label in ('baseline-attribution-v3','candidate-attribution-v1'):
        vectors=n.read(root/label/'pressure-vectors.json');dg[label]=diagnosis(bodies[label],vectors,scale);write_json(a.output/(label+'.json'),compact(bodies[label]));write_json(a.output/(label+'-vectors.json'),vectors)
    require(bodies['baseline-attribution-v3']['selection']['knots']==pre['selection']['knots'],'postdiagnosis changed the prospectively frozen single contrast')
    decisions={}
    for label in ('baseline-T008-v1','candidate-T008-v1','candidate-T020-v2'):
        q=bodies[label];decisions[label]=pair(q);out=compact_pair(q);out['common_times']=q.get('common_times',[])
        if 'dense_fine_transfer' in q:out['dense_fine_transfer']=q['dense_fine_transfer']
        for dt,r in q['members'].items():out['members'][dt]['execution']=r.get('execution',{});out['members'][dt]['first_accuracy_crossing']=r.get('first_accuracy_crossing')
        write_json(a.output/(label+'.json'),compact(out));write_json(a.output/(label+'-tractions.json'),member_profiles(q))
    for label in ('controls-v1','common-space-v2','target-preflight-v1'):
        write_json(a.output/(label+'.json'),compact(bodies[label]))
    write_json(a.output/'target-audit-capacity.json',compact(auditcap))
    profiles={}
    for q in [bodies[k] for k in ('baseline-T008-v1','candidate-T008-v1','candidate-T020-v2')]+([target] if target else []):
        for row in q['members'].values():
            file=Path(row['directory'])/'timing.jsonl'
            if file.exists():profiles[str(file)]=timing(file)
    write_json(a.output/'cost.json',profiles)
    if target:
        out=compact_pair(target);out['qualification']=target_pair(target);out['first_stage_resource_plan']=target['first_stage_resource_plan'];out['members']={dt:dict(r,execution=target['members'][dt].get('execution',{})) for dt,r in out['members'].items()};write_json(a.output/'target-pair.json',compact(out));write_json(a.output/'target-tractions.json',member_profiles(target))
    res=resources(root);resource_ok=bool(res) and all(q['resource']['aggregate_status']==q['resource']['project_status']=='pass' for q in res.values());write_json(a.output/'resources.json',compact(res))
    attempts={}
    for name in ('build-v1','diagnosis-v1','common-space-v1','candidate-T020-v1'):
        folder=root/name;attempts[name]=dict(path=str(folder),files={str(f.relative_to(folder)):dict(sha256=sha(f),bytes=f.stat().st_size) for f in folder.glob('*') if f.is_file() and f.name in ('request-manifest.json','progress.json','results.json','1.log','independent.json')},scope='preserved diagnostic/controller preparation; completed native artifacts reused only by explicit checked lineage')
    write_json(a.output/'retained-attempts.json',attempts)
    pressure_reduction=bodies['candidate-attribution-v1']['total_images']['actual']['whole']['weighted_RMS']/bodies['baseline-attribution-v3']['total_images']['actual']['whole']['weighted_RMS']
    decision=dict(schema='CT10RequiredScopeV1',diagnostic={k:v for k,v in dg.items()},observational_controls=bodies['controls-v1']['status'],static_contrast=pre['status'],single_contrast=dict(knots=pre['selection']['knots'],native_slots=pre['initial']['candidate-0.001']['native_count']['stored_slots'],selection_signal=pre['selection']['signal_captured']),baseline_physical_uncertainty=decisions['baseline-T008-v1'],baseline_backend_transfer=bodies['baseline-T008-v1']['dense_fine_transfer']['status'],candidate_T008=decisions['candidate-T008-v1'],candidate_T020=decisions['candidate-T020-v2'],initial_momentum_pressure_image_ratio=pressure_reduction,target_static=bodies['target-preflight-v1']['static_admission'],target_pair=target_pair(target) if target else dict(status='unattempted'),actual_target_horizon=target['target'] if target else 0.,first_01_transit='not_started',resources='pass' if resource_ok else 'fail',source_artifacts=sources,terrain='not_admitted; mapped numbering,joined gradients/cofactors,transmission,curved temporal/spatial qualification and physical-interface/approach-flow admission remain',physical_applicability='central reference feasibility retained; no experimental validation; eight extraction failures and all historical qualifications unchanged')
    if target:decision['first_stage_resource_plan']=target['first_stage_resource_plan']
    decision['target_required_audit_capacity']=dict(status=auditcap['status'],first_failure=auditcap.get('first_failure'),scope=auditcap['physical_budget_qualification'])
    if not target:decision['target_blocker']='required physical observer capacity failed' if auditcap['status']!='pass' else 'target paired stepping still required'
    decision['investigation_completion']='complete' if all(v['status']=='pass' for v in dg.values()) and decision['observational_controls']==decision['static_contrast']==decision['baseline_backend_transfer']==decision['resources']=='pass' else 'incomplete'
    decision['target_continuation_admitted']=bool(target and target_pair(target)['status']=='pass' and target['first_stage_resource_plan']['resource_admitted'])
    require(not decision['target_continuation_admitted'],'target is admitted and affordable: execute its prospectively declared first stage before final collection')
    write_json(a.output/'decision.json',decision)
    # Exact source roles, immutable governing binary and physical-file identities.
    src=root/'source';git=['git','-c','core.preloadIndex=false','-c','index.threads=1','-c','core.filemode=false'];revision=subprocess.check_output(git+['rev-parse','HEAD'],cwd=src,text=True).strip();status=subprocess.check_output(git+['status','--porcelain'],cwd=src,text=True);require(not status,'dirty publication controller source')
    physical=['library/SolverUtils/Turbulence/RANSModels.cpp','library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSparse.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp','docs/plans/cao-tamura-accepted-8-evidence/physical-reference-manifest.json']
    hashes={}
    for name in physical:
        old=hashlib.sha256(subprocess.check_output(git+['show','af87414abc2b583492c3a9f3466b25cb0fe5285a:'+name],cwd=src)).hexdigest();require(old==sha(src/name),'governing/physical source unexpectedly changed:'+name);hashes[name]=old
    write_json(a.output/'source-scope.json',dict(start='af87414abc2b583492c3a9f3466b25cb0fe5285a',controller_revision=revision,diagnostic_compiled_revisions={k:n.read(root/k/'request-manifest.json')['request']['revision'] for k in ('build-v2','build-v3')},governing=dict(root=str(a.ct9/'build-v8/runtime'),solver_sha256=sha(a.ct9/'build-v8/runtime/IncNavierStokesSolver'),component='CSR1:EigenLU-LDLT-IR3-F40M:FGMRES80x1000'),unchanged_governing_and_physics=hashes,diagnostic_scope='read-only attribution and bounded face serialization; no governing runtime rebuilt',request_and_state_provenance=sources))
    write_json(a.output/'runtime.json',compact({k:bodies[k] for k in ('build-v2','build-v3')}));write_json(a.output/'loaded-libraries.json',LIBRARIES)
    files=[f.name for f in a.output.iterdir() if f.suffix=='.json' and f.name not in ('request-manifest.json','completed-artifacts.json')];m.seal(a.output,files);print('CT10 final investigation',decision['investigation_completion'],'target',decision['target_pair'],flush=True);return 0 if decision['investigation_completion']=='complete' else 1
if __name__=='__main__':raise SystemExit(main())
