# SPDX-License-Identifier: MIT
"""CT8 primitive admission; collection and process completion are not gates."""
import math
from ct_accepted_budget_gate import physical_budgets

FIELDS={'U','w','k','outer_epsilon','Pi'}
def finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)

def pair(result):
    reasons=[];target=result.get('target');members=result.get('members',{})
    if not finite(target) or target not in (.008,.02):reasons.append('wrong_horizon')
    if not finite(target):return dict(status='fail',reasons=reasons)
    if set(members)!={'0.002','0.001'}:reasons.append('missing_timestep_member')
    for name,row in members.items():
        if not finite(row.get('time')) or abs(row['time']-target)>1e-13 or row.get('failure'):reasons.append(name+':incomplete_or_rejected')
        if ':resolved=A95:OriginalHighRe-v1' not in row.get('signature','') or 'SamePhysicalWallDisplacedAcceptedV1' not in row.get('signature','') or 'CT8PhysicalBelow25mmEqualIds012V1' not in row.get('signature',''):reasons.append(name+':foreign_identity')
        metrics=row.get('accepted_metrics',[])
        if not metrics or abs(metrics[-1].get('time',-1)-target)>1e-13:reasons.append(name+':missing_endpoint_metrics');continue
        metric=metrics[-1];fields=metric.get('drift',{})
        if set(fields)!=FIELDS:reasons.append(name+':missing_field')
        for f,q in fields.items():
            scale=q.get('scale');rms=q.get('weighted_RMS')
            if not finite(scale) or scale<=0 or not finite(rms) or rms<0 or rms/scale>(1e-4 if f=='w' else 1e-3):reasons.append(name+':drift:'+f)
        mass=metric.get('mass',{});v,scale=mass.get('outward_net'),mass.get('physical_scale')
        if not finite(v) or not finite(scale) or scale<=0 or abs(v)/scale>1e-8:reasons.append(name+':mass')
        for key in ('normalized_wall_traction_error','normalized_interface_traction_RMS'):
            v=metric.get('traction',{}).get(key)
            if not finite(v) or v<0 or v>1e-3:reasons.append(name+':'+key)
        if physical_budgets(row.get('budgets',{}))['status']!='pass':reasons.append(name+':budget_primitives')
        if row.get('domain',{}).get('status')!='pass' or row.get('constitutive',{}).get('status')!='pass':reasons.append(name+':physical_domain_or_constitutive')
        if row.get('corrections')!=0 or row.get('extra_solves')!=0:reasons.append(name+':forbidden_correction')
    temporal=result.get('temporal',{})
    if set(temporal)!=FIELDS:reasons.append('missing_temporal_field')
    for f,q in temporal.items():
        drift=q.get('measured_drift');d=q.get('difference',{});rms,scale=d.get('weighted_RMS'),d.get('scale')
        if not finite(drift) or drift<0 or not finite(rms) or rms<0 or not finite(scale) or scale<=0 or rms/scale>max(1e-7,.1*drift):reasons.append('temporal:'+f)
    return dict(status='pass' if not reasons else 'fail',reasons=reasons,scope='bounded strip only; complete prerequisites and model/state/history binding are required independently')
