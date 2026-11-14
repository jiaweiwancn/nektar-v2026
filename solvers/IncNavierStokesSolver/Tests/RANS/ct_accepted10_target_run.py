#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prospectively bounded actual-target pair after independent strip/static gates.

This is a paired stepping/cost check, not a flow-through claim. The old dense
guard, sparse component, physical functions, one-sweep method and tests remain.
"""
import argparse,json,math,os,resource,xml.etree.ElementTree as ET
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import ct_accepted8_native as n
import numpy as np
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair as strip_pair,FIELDS,finite
from ct_accepted9_controls import enable
from ct_accepted10_run import RUNTIME_SHA,comparison
from ct_accepted3_metrics import inventory
from ct_accepted_pilot import replay,constitutive_comparison,artifact_hash,budget_classification
from ct_accepted6_audit import captured_k_matrix
from ct_accepted_budget_gate import physical_budgets
from ct_accepted_open import run,probe,require
from ct_accepted5_traction import profiles
from ct_accepted_feasibility import write_json
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity


def target_pair(q):
    reasons=[]
    if q.get('schema')!='CT10TargetPhysicalPairV1' or q.get('target')!=.002 or q.get('required_dt')!=[.001,.0005]:reasons.append('foreign_target_request')
    if set(q.get('members',{}))!={'0.001','0.0005'}:reasons.append('missing_member')
    for name,r in q.get('members',{}).items():
        if r.get('dt')!=float(name) or not finite(r.get('time')) or abs(r['time']-.002)>1e-13 or r.get('failure'):reasons.append(name+':incomplete_or_rejected')
        expected=[round(j*float(name),12) for j in range(round(.002/float(name))+1)]
        metrics=r.get('accepted_metrics',[])
        if [round(v.get('time',-1),12) for v in metrics]!=expected:reasons.append(name+':accepted_coverage')
        if any(k not in r.get('signature','') for k in (':resolved=A95:OriginalHighRe-v1','SamePhysicalWallDisplacedAcceptedV1','DisplacedNeutralCompositeV1','AcceptedRichardsNorrisFlatV1','EigenLU-LDLT-IR3-F40M')):reasons.append(name+':foreign_physics_or_backend')
        for v in metrics:
            if set(v.get('drift',{}))!=FIELDS:reasons.append(name+':field_coverage')
            for f,a in v.get('drift',{}).items():
                x,s=a.get('weighted_RMS'),a.get('scale')
                if not finite(x) or x<0 or not finite(s) or s<=0 or x/s>(1e-4 if f=='w' else 1e-3):reasons.append(name+':drift:'+f)
            mass=v.get('mass',{});x,s=mass.get('outward_net'),mass.get('physical_scale')
            if not finite(x) or not finite(s) or s<=0 or abs(x)/s>1e-8:reasons.append(name+':mass')
            for k in ('normalized_wall_traction_error','normalized_interface_traction_RMS'):
                x=v.get('traction',{}).get(k)
                if not finite(x) or not 0<=x<=.001:reasons.append(name+':traction:'+k)
        if physical_budgets(r.get('budgets',{}))['status']!='pass':reasons.append(name+':budget_coverage_or_accuracy')
        binding=r.get('budgets',{}).get('budget_binding',{})
        if not binding or binding.get('state_sha256')!=r.get('state_sha256') or binding.get('configuration_sha256')!=r.get('configuration_sha256') or binding.get('native_report_sha256')!=r.get('native_report_sha256') or binding.get('dt')!=r.get('dt') or binding.get('time')!=r.get('time'):reasons.append(name+':budget_state_configuration_identity')
        cap=r.get('budgets',{}).get('captured_k_matrix',{});x,s=cap.get('maximum_defect'),cap.get('scale')
        if not finite(x) or x<0 or not finite(s) or s<=0 or x/s>2e-10:reasons.append(name+':captured_scalar_operator')
        b=r.get('initialization_receipt',{})
        if not b.get('open_explicit_implicit_agreement') or not b.get('explicit_boundary_sha256') or b['explicit_boundary_sha256']!=b.get('implicit_boundary_sha256'):reasons.append(name+':boundary_identity')
        if r.get('domain',{}).get('status')!='pass' or r.get('constitutive',{}).get('status')!='pass':reasons.append(name+':domain_or_constitutive')
        if r.get('corrections')!=0 or r.get('extra_solves')!=0:reasons.append(name+':correction')
    if q.get('same_actual_initial') is not True:reasons.append('unequal_actual_initial')
    if set(q.get('temporal',{}))!=FIELDS:reasons.append('missing_temporal_fields')
    for f,v in q.get('temporal',{}).items():
        x=v.get('difference',{}).get('weighted_RMS');s=v.get('difference',{}).get('scale');d=v.get('measured_drift')
        if not finite(x) or x<0 or not finite(s) or s<=0 or not finite(d) or d<0 or x/s>max(1e-7,.1*d):reasons.append('temporal:'+f)
    if set(q.get('traction_temporal',{}))!={'ground','interface','signed_mean'}:reasons.append('missing_traction_temporal')
    for f,v in q.get('traction_temporal',{}).items():
        x,d=v.get('spread'),v.get('measured_drift')
        if not finite(x) or x<0 or not finite(d) or d<0 or x>max(1e-7,.1*d):reasons.append('traction_temporal:'+f)
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)),scope='actual unchanged target,pairedT=.002 stepping only')


def bounded_budgets(rt,folder,ini,r,c):
    """Same existing physical integration; omit redundant exported volume rows."""
    attempt=next(x for x in reversed(r['pimple_attempts']) if x['scheduled_complete'])
    binding=dict(model=c['model'],policy=c['wall_policy'],dt=r['dt'],time=r['time'],state_sha256=exact_vector_hash(r['state']),configuration_sha256=sha(folder/'case.xml'),native_report_sha256=sha(folder/'case.evolution.json'),actual_k_stage_sha256=attempt['k_snapshot_sha256'])
    result=dict(budgets={},budget_binding=binding);path=Path(attempt['k_snapshot']);require(sha(path)==attempt['k_snapshot_sha256'],'changed actual target scalar stage');snap=n.read(path)
    for label,orders,face in n.CONTROLS:
        dest=folder/label;dest.mkdir();root=ET.parse(folder/'case.xml').getroot()
        for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
        write(root,dest/'case.xml');req=replay(ini,r['state'],r['time'],[],export_physical_scalar_budget=True,export_frozen_scalar_reconciliation=True,previous_state=r['accepted_states'][-2],stage_history=snap['history'],actual_flow_stage=snap['post_flow'],actual_pre_k_stage=snap['pre_k'],physical_face_order=face)
        try:
            q,e=probe(rt,dest,req,'budget');require(e['exit_code']==0,'native target physical budget failed:'+str(q.get('failure')));require(q['observation']['state']==r['state'] and q['observation']['identity']==ini['observation']['identity'],'target budget state/identity mismatch')
            b=budget_classification(q,c);b.update(execution=e,packet_sha256=sha(dest/'budget.json'),state_sha256=exact_vector_hash(r['state']),binding=binding,control_configuration_sha256=sha(dest/'case.xml'));result['budgets'][label]=b
            if label=='native':result['captured_k_matrix']=captured_k_matrix(q,snap)
        except Exception as err:result['budgets'][label]=dict(status='failed_required_control',reason=repr(err),directory=str(dest))
        write_json(folder/'budget-progress.json',result)
    result['gate']=physical_budgets(result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct9','runtime','probe','preflight','strip','contract','audit-capacity','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');n.reference=lru_cache(None)(n.reference)
    pre=sealed(a.preflight);strip=sealed(a.strip);require(strip_pair(strip)['status']=='pass' and strip['target']==.02,'complete strip continuation gate failed');require(pre['static_admission']['status']=='pass','target static gate failed');require(sha(a.runtime/'IncNavierStokesSolver')==RUNTIME_SHA,'foreign target governing runtime')
    c=pre['mesh'];scale=pre['mesh_scale_preflight'];require(scale['minimum_dx']>=.125 and abs(scale['minimum_dz']-.0029296875)<1e-15,'target mesh scales do not support the predeclared strip-informed timestep choice')
    cap=sealed(a.audit_capacity);cap_request=n.read(a.audit_capacity/'request-manifest.json')['request'];require(cap_request['case']==c and cap_request['actual_initial_sha256']==pre['actual_initial_sha256'],'foreign required-audit capacity fixture')
    require(set(cap.get('controls',{}))=={'native','q48f16','q64f16','q64f32'},'required audit capacity coverage incomplete')
    for label,q in cap['controls'].items():
        require(q.get('execution',{}).get('exit_code')==0 and q.get('status')=='pass' and set(q.get('equation_names',[]))=={'inner_k','outer_k','outer_epsilon','whole_k'},'required audit capacity failed:'+label)
        require(set(q.get('frozen_load_reconstruction',{}))=={'k','outer_epsilon'} and all(0<=r['maximum_defect']/r['normalization']<=2e-10 for r in q['frozen_load_reconstruction'].values()),'required audit reconstruction missing/failed:'+label)
    gov_identity=runtime_identity(a.runtime/'IncNavierStokesSolver');os.environ['LD_LIBRARY_PATH']=str(a.probe/'lib');probe_identity=runtime_identity(a.probe/'RANSPreservationProbe');os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    m=RequestManifest(a.output,dict(schema='CT10TargetSteppingV1',case=c,preflight=sha(a.preflight/'results.json'),audit_capacity=sha(a.audit_capacity/'results.json'),strip=sha(a.strip/'results.json'),contract=sha(a.contract),source=sha(Path(__file__)),runtime=gov_identity,probe=probe_identity,target=.002,dt=[.001,.0005],steps=[2,4],flow_through=40.621637857609905,scope='bounded actual target stepping and cost; no automatic .01transit admission'))
    out=dict(schema='CT10TargetPhysicalPairV1',required_dt=[.001,.0005],target=.002,case=c,members={});initials={};reports={};ini=n.read(a.preflight/'target/readiness.json');weights=n.inventory(c)[2];count=len(weights);first=None
    for dt in (.001,.0005):
        name=str(dt);zero=a.output/('zero-'+name);n.make(zero,a.runtime,c,dt,0);enable(zero);e=run([a.runtime/'IncNavierStokesSolver','case.xml'],zero,'solver.log');z=n.read(zero/'case.evolution.json');require(e['exit_code']==0 and not z['failure'] and z['time']==0. and exact_vector_hash(z['state'])==pre['actual_initial_sha256'],'fresh target zero-step state differs from admitted actual representation')
        require(z['state']==ini['observation']['state'],'initial target coefficient mismatch');require(first is None or first==z['state'],'timestep target initial coefficients differ');first=z['state'];initials[name]=ini
        folder=a.output/name;n.make(folder,a.runtime,c,dt,round(.002/dt),Path(z['checkpoint']),zero/'case.xml');os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);r=n.read(folder/'case.evolution.json');require(r['accepted_states'][0]==z['state'],'target did not consume native history')
        row=dict(dt=dt,time=r['time'],directory=str(folder),failure=r['failure'],signature=r['signature'],execution=e,actual_initial_sha256=exact_vector_hash(z['state']),state_sha256=exact_vector_hash(r['state']),configuration_sha256=sha(folder/'case.xml'),native_report_sha256=sha(folder/'case.evolution.json'),history_sha256=artifact_hash(r['checkpoint']),initialization_receipt=r['initialization_receipt'])
        row['accepted_metrics']=[dict(time=t,**n.physical_metrics(ini,s,c,weights,count)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['traction_history']=[dict(time=t,profile=profiles(ini,s,c)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['actual_initial_traction_profile']=profiles(ini,z['state'],c)
        fixed=inventory(ini['test_basis'],ini['observation']['grid'],np.asarray(r['state']),ini['observation']['offsets'],'A95',ini['observation']['grid_sides']);row['domain']=dict(status='pass' if not r['failure'] and fixed['domain_failure_count']==0 and fixed['minima']['covariance']['value']>=-1e-10*c['K0'] else 'fail',fixed=fixed,native_scope='mandatory production accepted-state and stage checks retained in native report',native_audited_extrema=audited_extrema(r))
        row['budgets']=bounded_budgets(a.runtime,folder,ini,r,c) if r['steps'] else dict(gate=dict(status='missing'))
        grid,sides,_,_=n.inventory(c);req=replay(ini,r['state'],r['time'],grid,grid_sides=sides,export_constitutive=True);q,pe=probe(a.runtime,folder,req,'constitutive');row['constitutive']=constitutive_comparison(q,c) if pe['exit_code']==0 else dict(status='fail',failure=q.get('failure'))
        row['corrections']=sum(x.get('bound_calls',0) for x in r['pimple_attempts']);row['extra_solves']=sum(x.get('additional_governing_equation_solves',0) for x in r['pimple_attempts']);row['actual_stage_checks']=[{k:x.get(k) for k in ('time','pressure_krylov_iterations','frozen_momentum_defect','corrected_continuity_gauge','bound_calls','additional_governing_equation_solves')} for x in r['pimple_attempts']];out['members'][name]=row;reports[name]=r;write_json(a.output/'progress.json',out);print('CT10 target',name,r['time'],r['failure'],flush=True)
    # Reuse the common arithmetic, with this separately versioned target gate.
    comparison(out,initials,reports,c,weights,count,classifier=target_pair)
    samples=[];overheads=[];timing=[]
    for r in out['members'].values():
        path=Path(r['directory'])/'timing.jsonl';by_step={}
        for line in path.read_text().splitlines():
            z=json.loads(line)
            if z['step']>0:by_step[z['step']]=by_step.get(z['step'],0.)+z['exclusive_process_cpu_seconds']
        require(len(by_step)==len(r['accepted_metrics'])-1,'accepted target timing coverage missing')
        samples.extend(by_step.values());overheads.append(max(0.,r['execution']['process_cpu_seconds']-sum(by_step.values())));timing.append(dict(path=str(path),sha256=sha(path),exclusive_step_CPU=by_step))
    per_step=max(samples);needed=round((.407-.002)/.001)+round((.407-.002)/.0005);segments=math.ceil(round((.407-.002)/.001)/32)+math.ceil(round((.407-.002)/.0005)/32);estimate=1.25*(per_step*needed+max(overheads)*segments)
    byte_rows=[dict(directory=r['directory'],bytes=sum(p.stat().st_size for p in Path(r['directory']).rglob('*') if p.is_file()),accepted_steps=len(r['accepted_metrics'])-1) for r in out['members'].values()];per_step_bytes=max(x['bytes']/max(1,x['accepted_steps']) for x in byte_rows);storage=1.25*per_step_bytes*needed
    out['first_stage_resource_plan']=dict(completed_horizon=.002,transit_fraction=.002/40.621637857609905,prospective_endpoint=.407,additional_pair_steps=needed,segments_at_most32steps=segments,maximum_measured_exclusive_accepted_step_CPU=per_step,maximum_measured_nonstep_native_CPU=max(overheads),CPU_with_25percent_headroom=estimate,CPU_limit=24000,measured_output=byte_rows,storage_bytes_with_25percent_headroom=storage,storage_limit_bytes=24*1024**3,resource_admitted=estimate<=24000 and storage<=24*1024**3,physical_admitted=out['qualification']['status']=='pass',timing=timing,scope='exclusive accepted-step CPU separated from startup/output; raw output measured separately,immutable histories retained; no flow-through claim')
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT10 target',out['qualification'],out['first_stage_resource_plan'],flush=True);return 0 if out['qualification']['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
