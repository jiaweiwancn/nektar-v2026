# SPDX-License-Identifier: MIT
"""CT11 prerequisite gates recomputed from complete finite primitives."""
from ct_accepted9_gate import finite
from ct_accepted_budget_gate import CONTROLS,EQUATIONS,RECONSTRUCTION
from ct_accepted11_equivalence import CASES,POLICY


def observer(q):
    reasons=[]
    if q.get('schema')!='CT11ObserverEquivalenceV1' or q.get('policy')!=POLICY:reasons.append('foreign_observer_request')
    if set(q.get('checks',{}))!=CASES:reasons.append('missing_fixture')
    for name,rows in q.get('checks',{}).items():
        if set(rows)!=CONTROLS:reasons.append(name+':missing_control')
        for label,r in rows.items():
            d=r.get('all_primitive_difference',{});v=d.get('maximum_scaled_difference');count=d.get('compared_numeric_entries')
            if not finite(v) or not 0<=v<=2e-10 or not finite(count) or count<=0:reasons.append(name+':operator_primitives')
            meta=r.get('volume_observer',{})
            if meta.get('policy')!=POLICY or meta.get('materialized_global_volume_rows') is not False or not meta.get('native_rows_provenance_fnv1a64') or any(not finite(meta.get(k)) or meta[k]<=0 for k in ('inner_points','outer_points')):reasons.append(name+':volume_identity')
            b=r.get('classification',{})
            if set(b.get('equations',{}))!=EQUATIONS or set(b.get('reconstruction',{}))!=RECONSTRUCTION:reasons.append(name+':coverage')
            for f,v in b.get('reconstruction',{}).items():
                x,s=v.get('maximum_defect'),v.get('normalization')
                if not finite(x) or x<0 or not finite(s) or s<=0 or x/s>2e-10:reasons.append(name+':reconstruction:'+f)
            if r.get('execution',{}).get('exit_code')!=0 or r.get('raw_state_unchanged') is not True:reasons.append(name+':execution_or_state')
            if not all(r.get(k) for k in ('state_sha256','previous_state_sha256','flow_stage_sha256','pre_k_stage_sha256','configuration_sha256','packet_sha256')):reasons.append(name+':missing_identity')
            if label=='native':
                p=r.get('native_provenance') or {};a=p.get('legacy',{});b=p.get('stream',{})
                if not a.get('native_rows_provenance_fnv1a64') or a.get('native_rows_provenance_fnv1a64')!=b.get('native_rows_provenance_fnv1a64') or a.get('physical_weights')!=b.get('physical_weights'):reasons.append(name+':native_row_provenance')
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)),scope='observer all-row/face/actual-stage equivalence,not physical admission')


def capacity(q):
    reasons=[]
    if q.get('schema')!='CT11AuditCapacityV1' or q.get('policy')!=POLICY or q.get('native_trajectory')!='not_run':reasons.append('foreign_capacity_request')
    if set(q.get('controls',{}))!=CONTROLS:reasons.append('missing_capacity_control')
    for name,r in q.get('controls',{}).items():
        e=r.get('execution',{})
        if e.get('exit_code')!=0 or e.get('control_wall_limit_reached') is not False:reasons.append(name+':execution')
        for k in ('peak_RSS_KiB','peak_virtual_KiB','process_cpu_seconds','wall_seconds'):
            if not finite(e.get(k)) or e[k]<=0:reasons.append(name+':resource_primitive')
        if finite(e.get('peak_virtual_KiB')) and e['peak_virtual_KiB']>16*1024**2:reasons.append(name+':memory_ceiling')
        if set(r.get('equation_names',[]))!=EQUATIONS or set(r.get('frozen_load_reconstruction',{}))!=RECONSTRUCTION:reasons.append(name+':coverage')
        for f,v in r.get('frozen_load_reconstruction',{}).items():
            x,s=v.get('maximum_defect'),v.get('normalization')
            if not finite(x) or x<0 or not finite(s) or s<=0 or x/s>2e-10:reasons.append(name+':reconstruction:'+f)
        meta=r.get('volume_observer',{})
        if meta.get('policy')!=POLICY or meta.get('materialized_global_volume_rows') is not False:reasons.append(name+':observer_identity')
        if not all(r.get(k) for k in ('request_sha256','configuration_sha256','packet_sha256')):reasons.append(name+':missing_identity')
    return dict(status='pass' if not reasons else 'fail',reasons=sorted(set(reasons)),scope='actual initialized target capacity/reconstruction only; no trajectory or finite-defect accuracy qualification')
