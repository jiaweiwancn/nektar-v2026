# SPDX-License-Identifier: MIT
"""Prospective CT9 pairs have real timestep identities; CT8 is immutable."""
import math
from ct_accepted_budget_gate import physical_budgets

FIELDS={'U','w','k','outer_epsilon','Pi'}
def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
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
