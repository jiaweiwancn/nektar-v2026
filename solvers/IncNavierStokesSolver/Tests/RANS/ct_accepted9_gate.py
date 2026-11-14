# SPDX-License-Identifier: MIT
"""Prospective CT9 pairs have real timestep identities; CT8 is immutable."""
import math
from ct_accepted_budget_gate import physical_budgets

FIELDS={'U','w','k','outer_epsilon','Pi'}
def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)

def backend_trajectory(result):
    """Recompute equivalence and common-time coverage, independently of labels."""
    reasons=[];members=result.get('members',{});comparisons=result.get('equivalence',{})
    if result.get('schema')!='CT9PhysicalPairV1' or result.get('required_dt')!=[.002,.001] or result.get('target')!=.008:reasons.append('foreign_equivalence_request')
    if set(members)!={'0.002','0.001'} or set(comparisons)!={'0.002','0.001'}:reasons.append('missing_backend_member')
    for name in ('0.002','0.001'):
        member=members.get(name,{});q=comparisons.get(name,{});rows=q.get('common_times',[])
        expected=[round(i*float(name),12) for i in range(round(.008/float(name))+1)]
        times=[round(r.get('time',-1),12) for r in rows if finite(r.get('time'))]
        if times!=expected:reasons.append(name+':missing_equal_time_coverage')
        if not q.get('dense_report') or q.get('sparse_report')!=member.get('native_report_sha256'):reasons.append(name+':report_identity')
        for row in rows:
            if set(row.get('fields',{}))!=FIELDS:reasons.append(name+':missing_equivalence_field')
            for f,v in row.get('fields',{}).items():
                x,s=v.get('weighted_RMS'),v.get('scale')
                if not finite(x) or x<0 or not finite(s) or s<=0 or x/s>1e-8:reasons.append(name+':equivalence:'+f)
            for f in ('ground','interface'):
                x=row.get('traction',{}).get('profiles',{}).get(f,{}).get('normalized_weighted_RMS')
                if not finite(x) or not 0<=x<=1e-8:reasons.append(name+':traction_equivalence:'+f)
    # The physical pressure failure is expected; every other physical primitive
    # must still pass. A foreign failure cannot masquerade as equivalence.
    if pair(result)['reasons']!=['temporal:Pi']:reasons.append('different_physical_classification')
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)))

def capacity(result):
    """Static capacity is evidence of actual actions, never physical admission."""
    reasons=[];rows=result.get('cases',{});expected={'intermediate-p3':8864,'target-p3':33267,'target-p4':60907}
    if result.get('schema')!='CT9FlatCapacityV1' or set(rows)!=set(expected):reasons.append('missing_or_foreign_capacity')
    for name,count in expected.items():
        row=rows.get(name,{});q=row.get('sparse_backend',{})
        if row.get('status')!='static_capacity_pass' or row.get('zero_time')!=0. or row.get('startup_physical_failure'):reasons.append(name+':incomplete_startup')
        if row.get('native_count',{}).get('stored_slots')!=count or q.get('stored_slots')!=count:reasons.append(name+':native_count')
        for phase in ('startup_execution','frozen_execution'):
            if row.get(phase,{}).get('exit_code')!=0:reasons.append(name+':'+phase)
        if row.get('reference_agreement',{}).get('status')!='pass':reasons.append(name+':reference_binding')
        for key in ('flow_rows','flow_nnz','velocity_rows','pressure_rows','predictor_factor_nnz','consistent_mass_factor_nnz','pressure_preconditioner_factor_nnz'):
            x=q.get(key)
            if not finite(x) or x<=0:reasons.append(name+':missing_allocation:'+key)
        if q.get('pressure_C_nnz')!=0 or not finite(q.get('pressure_iterations')) or not 0<=q['pressure_iterations']<=1000:reasons.append(name+':pressure_system')
        masses=q.get('masses',[])
        if len(masses)!=6 or any(not finite(m.get('rows')) or m['rows']<=0 or not finite(m.get('nnz')) or m['nnz']<=0 for m in masses):reasons.append(name+':mass_coverage')
        for error,rhs in [('predictor_residual','predictor_rhs_maximum'),('true_pressure_residual','pressure_rhs_maximum'),('consistent_mass_correction_true_residual','consistent_mass_correction_rhs_maximum')]:
            x,s=q.get(error),q.get(rhs)
            if not finite(x) or x<0 or not finite(s) or s<0 or x>2e-12*s+1e-25:reasons.append(name+':true_residual:'+error)
        scalars=q.get('scalars',[])
        if len(scalars)!=2 or {s.get('block') for s in scalars}!={2,3}:reasons.append(name+':scalar_coverage')
        for s in scalars:
            x,b=s.get('maximum_true_residual'),s.get('rhs_maximum')
            if not finite(x) or x<0 or not finite(b) or b<0 or x>2e-12*b+1e-25:reasons.append(name+':scalar_true_residual')
        for key in ('full_native_G_action_defect','true_continuity','splitting_identity'):
            x=q.get(key)
            if not finite(x) or not 0<=x<=1e-10:reasons.append(name+':native_action:'+key)
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)),scope='zero-step projection, complete frozen actions and resource capacity only; not target preservation')
def backend_controls(result):
    """Required primitive scope; a pass label cannot replace its evidence."""
    reasons=[];rows=result.get('checks',{});operators={'d0-p3','d0-p4','displaced-p3','displaced-p4'}
    if result.get('schema')!='CT9BackendControlsV2':reasons.append('foreign_control_schema')
    if set(rows)!=operators|{'restart_rollback_backend_identity','absolute_pressure_level'}:reasons.append('missing_required_control')
    required={'raw_residual','constrained_residual','external_constraints','traction_load','diffusive_load','advection_load','predictor','pressure_increment','consistent_mass_correction','scalar_increment_2','scalar_increment_3'}|{'mass_'+str(i) for i in range(6)}|{'block_'+str(i)+'_'+key for i in (0,2,3) for key in ('matrix','residual')}
    for name in operators:
        q=rows.get(name,{});d=q.get('defects',{})
        if q.get('status')!='pass' or not required<=set(d) or any(not finite(x) or not 0<=x<=2e-10 for x in d.values()):reasons.append(name+':operator_primitives')
        for key in ('negative_roughness_pressure_omission_signal','negative_sign_signal','negative_mass_replacement_signal','negative_stale_mask_signal'):
            x=q.get(key)
            if not finite(x) or x<=1e-9:reasons.append(name+':nondiscriminating_negative')
        artifacts=q.get('artifacts',[])
        if len(artifacts)!=2:reasons.append(name+':missing_native_member')
        for a in artifacts:
            r=a.get('independent',{});rank=r.get('pressure_minimum_singular_ratio');constant=r.get('constant_pressure_action');defects=r.get('defects',{})
            if r.get('status')!='pass' or not defects or any(not finite(x) or not 0<=x<=2e-10 for x in defects.values()) or not finite(rank) or rank<=1e-10 or not finite(constant) or constant<=1e-8 or r.get('additional_pressure_null_modes')!=0:reasons.append(name+':independent_rank_or_operator')
    h=rows.get('restart_rollback_backend_identity',{});x=h.get('dense_sparse_state_maximum')
    if h.get('status')!='pass' or not h.get('uninterrupted_state') or h.get('uninterrupted_state')!=h.get('continued_state') or h.get('rebuilt_dense_CT8_history_bitwise_equal') is not True or not finite(x) or not 0<=x<=1e-10:reasons.append('history_primitives')
    negative={q.get('name'):q for q in h.get('negative',[])}
    for name,word in [('backend_history','signature'),('stale_fit','signature'),('wrong_model','a95'),('unsupported_cache','caching')]:
        q=negative.get(name,{})
        if q.get('status')!='pass' or word not in q.get('reason','').lower():reasons.append('history_negative:'+name)
    offset=rows.get('absolute_pressure_level',{});x=offset.get('maximum_coefficient_defect')
    if offset.get('status')!='pass' or offset.get('shift')!=.125 or not finite(x) or not 0<=x<=2e-10:reasons.append('absolute_pressure_primitive')
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)))

def pair(result):
    reasons=[];expected=result.get('required_dt',[]);members=result.get('members',{});target=result.get('target')
    if result.get('schema')!='CT9PhysicalPairV1':reasons.append('foreign_schema')
    if expected not in ([.001,.0005],[.002,.001]):reasons.append('undeclared_pair')
    if not finite(target) or target not in (.008,.02):reasons.append('wrong_horizon');return dict(status='fail',reasons=reasons)
    if set(members)!={str(t) for t in expected}:reasons.append('missing_timestep_member')
    for name,row in members.items():
        if not finite(row.get('dt')) or str(row['dt'])!=name:reasons.append(name+':wrong_native_dt')
        if not finite(row.get('time')) or abs(row['time']-target)>1e-13 or row.get('failure'):reasons.append(name+':incomplete_or_rejected')
        signature=row.get('signature','')
        if any(k not in signature for k in (':resolved=A95:OriginalHighRe-v1','SamePhysicalWallDisplacedAcceptedV1','CT8PhysicalBelow25mmEqualIds012V1')):reasons.append(name+':foreign_identity')
        metrics=row.get('accepted_metrics',[])
        if not metrics or abs(metrics[-1].get('time',-1)-target)>1e-13:reasons.append(name+':missing_endpoint');continue
        if any(q.get('status')!='pass' for q in metrics):reasons.append(name+':accuracy_crossing')
        for q in metrics:
            fields=q.get('drift',{})
            if set(fields)!=FIELDS:reasons.append(name+':missing_field')
            for f,v in fields.items():
                s,x=v.get('scale'),v.get('weighted_RMS')
                if not finite(s) or s<=0 or not finite(x) or x<0 or x/s>(1e-4 if f=='w' else 1e-3):reasons.append(name+':field:'+f)
            mass=q.get('mass',{});x,s=mass.get('outward_net'),mass.get('physical_scale')
            if not finite(x) or not finite(s) or s<=0 or abs(x)/s>1e-8:reasons.append(name+':mass')
            for k in ('normalized_wall_traction_error','normalized_interface_traction_RMS'):
                x=q.get('traction',{}).get(k)
                if not finite(x) or x<0 or x>1e-3:reasons.append(name+':'+k)
        if physical_budgets(row.get('budgets',{}))['status']!='pass':reasons.append(name+':budget_primitives')
        cap=row.get('budgets',{}).get('captured_k_matrix',{});x,s=cap.get('maximum_defect'),cap.get('scale')
        if not finite(x) or not finite(s) or s<=0 or x<0 or x/s>2e-10:reasons.append(name+':captured_matrix')
        receipt=row.get('initialization_receipt',{})
        if not receipt.get('open_explicit_implicit_agreement') or not receipt.get('explicit_boundary_sha256') or receipt.get('explicit_boundary_sha256')!=receipt.get('implicit_boundary_sha256'):reasons.append(name+':boundary_identity')
        if row.get('domain',{}).get('status')!='pass' or row.get('constitutive',{}).get('status')!='pass':reasons.append(name+':domain')
        if row.get('corrections')!=0 or row.get('extra_solves')!=0:reasons.append(name+':correction')
    if set(result.get('temporal',{}))!=FIELDS:reasons.append('missing_temporal_field')
    for f,q in result.get('temporal',{}).items():
        d=q.get('measured_drift');v=q.get('difference',{});x,s=v.get('weighted_RMS'),v.get('scale')
        if not finite(d) or d<0 or not finite(x) or x<0 or not finite(s) or s<=0 or x/s>max(1e-7,.1*d):reasons.append('temporal:'+f)
    if set(result.get('traction_temporal',{}))!={'ground','interface','signed_mean'}:reasons.append('missing_traction_temporal')
    for f,q in result.get('traction_temporal',{}).items():
        x,d=q.get('spread'),q.get('measured_drift')
        if not finite(x) or x<0 or not finite(d) or d<0 or x>max(1e-7,.1*d):reasons.append('traction_temporal:'+f)
    if result.get('same_actual_initial') is not True:reasons.append('different_actual_initial')
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)),scope='prospective CT9 bounded pair; old CT8 gate and labels unchanged')
