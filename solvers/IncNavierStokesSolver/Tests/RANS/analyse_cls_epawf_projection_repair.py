#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent identities, finite-QP checks and native restart/rollback evidence.

Read-only analysis of completed immutable packets. No governing solves and no
checkpoint migration. Run on Linux with one BLAS thread.
"""
import argparse
import collections
import json
from pathlib import Path
import platform
import numpy as np

from run_cls_epawf_picard import load,require,sha,state_hash,write_json
from continue_cls_epawf_pimple import header
from analyse_cls_epawf_constitutive_v2 import read_report
from compare_cls_epawf_pimple import same_history
from replay_cls_epawf_projection_limit import RAW_SHA,CP_SHA


def compact_bound(bound):
    return {k:v for k,v in bound.items() if k not in
            ('bounded_trace','active_constraints','current_inventory','inventory_counts')}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--uniform',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    require(platform.system()=='Linux','authoritative Linux analysis required')
    r=a.root.resolve();u=a.uniform.resolve();out=dict(schema=1,root=str(r),checks={},artifacts=[])
    def artifact(path,role):
        out['artifacts'].append(dict(path=str(path),sha256=sha(path),role=role))
    old=load(r/'replay-old-controls/manifest.json');new=load(r/'replay-total20000/manifest.json')
    result=load(r/'replay-total20000/result.json');q=load(r/'replay-total20000/archived-packet.json')
    d=load(r/'replay-total20000/projection-pool.json');n=len(d['free_rows'])
    mass=np.fromfile(d['mass_file'],dtype='<f8').reshape(n,n)
    require(sha(Path(d['mass_file']))==d['mass_sha256'],'projection metric hash differs')
    cuts=np.fromfile(d['cuts_file'],dtype='<f8').reshape(len(d['pool']),n)
    require(sha(Path(d['cuts_file']))==d['cuts_sha256'],'constraint pool hash differs')
    free=np.array(d['free_rows'],dtype=int);offset=int(result['offsets'][2]);index={int(row)-offset:i for i,row in enumerate(free)}
    E=np.zeros_like(cuts)
    for i,c in enumerate(d['pool']):
        for row,value in zip(c['physical_ids'],c['physical_phi']):
            if int(row) in index:E[i,index[int(row)]]+=value
    lower=np.linalg.cholesky(mass);white=np.linalg.solve(lower,E.T).T
    scales=np.linalg.norm(white,axis=1);white/=scales[:,None]
    white_error=float(np.max(np.abs(cuts-white)))
    require(white_error<1e-12,'exported physical mass/evaluation rows do not reconstruct native whitening')
    raw=np.array(result['raw']);corrected=np.array(result['corrected']);delta=(corrected-raw)[free]
    y=np.array(d['whitened_primal']);active=np.array(d['solved_active'],dtype=int);dual=np.array(d['solved_lambda'])
    rhs=np.array([c['rhs'] for c in d['pool']]);slack=cuts@y-rhs;action=cuts[active].T@dual
    objective=float(.5*delta@mass@delta);dual_objective=float(rhs[active]@dual-.5*action@action)
    stationarity=float(np.linalg.norm(y-action));complementarity=float(np.max(np.abs(dual*slack[active])))
    norm=float(np.sqrt(delta@mass@delta));b=result['k_bound'];scale=max(b['K0']**2,objective)
    require(np.max(rhs-cuts@y)<1e-12 and np.min(dual)>=-1e-12,'finite QP primal/dual infeasibility')
    require(stationarity<1e-10 and complementarity<1e-12 and abs(objective-dual_objective)<1e-12,'finite QP optimality failure')
    require(abs(norm-result['quadrature_correction_mass_norm'])<1e-12 and
            abs(b['integral_total']-result['quadrature_correction_integral'])<1e-12,'physical accounting differs')
    A=np.fromfile(r/'replay-total20000/transport_A.bin',dtype='<f8').reshape(n,n)
    raw_linear=np.array(q['transport_residual'])+A@np.array(q['raw_increment'])
    added=A@delta;storage=q['alpha']*(mass@delta);corrected_linear=raw_linear+added
    out['finite_qp']=dict(size=n,cuts=len(rhs),active=len(active),mass_sha256=d['mass_sha256'],cuts_sha256=d['cuts_sha256'],
        physical_rows_reconstruct_whitening_max=white_error,physical_objective=objective,dual_objective=dual_objective,
        primal_dual_gap=objective-dual_objective,scaled_gap=abs(objective-dual_objective)/scale,
        normalized_row_primal_violation=float(max(0.,np.max(-slack))),stationarity=stationarity,
        stationarity_over_K0=stationarity/b['K0'],complementarity=complementarity,
        scaled_complementarity=complementarity/scale,dual_violation=float(max(0.,-np.min(dual))),
        final_physical_inventory_violation=b['primal_floor_defect'],final_ground_violation=b['ground_floor_defect'],
        physical_quadrature_mass_norm=result['quadrature_correction_mass_norm'],physical_integral=b['integral_total'])
    out['frozen_accounting']=dict(raw_linear_max=float(np.max(np.abs(raw_linear))),
        raw_linear_readback_defect=float(np.max(np.abs(raw_linear-q['raw_linear_residual']))),
        corrected_linear_max=float(np.max(np.abs(corrected_linear))),added_A_delta_max=float(np.max(np.abs(added))),
        storage_alpha_M_delta_max=float(np.max(np.abs(storage))),
        A_delta_minus_alpha_M_delta_max=float(np.max(np.abs(added-storage))),
        transport_A_minus_alpha_M_max=float(np.max(np.abs(A-q['alpha']*mass))),
        definition='A is frozen transport, M is independently reconstructed PicardKMass; corrected residual is raw residual + A delta')
    for label,manifest in (('old',old),('larger_total',new)):
        out[label]=dict(status=manifest['status'],controls=manifest['controls'],corrected_sha256=manifest['corrected_sha256'],
            bound=compact_bound(manifest['bound']),runtime=manifest['probe'],execution=manifest['execution'])
    out['checks'].update(old_corrected_vector_bitwise_reproduced=old['corrected_sha256']==old['original_corrected_sha256'],
        complete_hard_audit_passed=result['complete_audit_failure']=='',
        exact_history_and_native_checkpoint=new['replay_context'],non_k_change=result['non_k_change'],
        essential_k_change=result['essential_k_change'],shared_cg_jump=result['shared_cg_jump'])
    require(out['checks']['old_corrected_vector_bitwise_reproduced'] and out['checks']['complete_hard_audit_passed'],'exact replay/audit check failed')
    events=[json.loads(x) for x in (r/'replay-total20000/trace.jsonl').read_text().splitlines()]
    separations=[v for v in events if v['kind']=='separation']
    out['progress']=dict(event_kinds=dict(collections.Counter(v['kind'] for v in events)),
        samples=[separations[i] for i in range(0,len(separations),100)]+[separations[-1]])
    baseline=load(u/'policy-checks/upstream6/case.evolution.json');earlier=load(r/'earlier-active6/case.evolution.json')
    out['checks']['earlier_active_bitwise_state_time_all_history_equal']=same_history(baseline,earlier)
    out['checks']['earlier_bound_result_hashes_equal']=[v['k_bound']['corrected_sha256'] for v in baseline['pimple_attempts']]==[v['k_bound']['corrected_sha256'] for v in earlier['pimple_attempts']]
    require(out['checks']['earlier_bound_result_hashes_equal'],'earlier correction changed under resource controls')
    cp=header(r/'native-inputs/checkpoint72',True);require(cp['sha256']==CP_SHA,'restart checkpoint changed')
    rollback=load(r/'rollback-active73/case.evolution.json')
    require(rollback['steps']==cp['steps'] and rollback['time']==cp['time'],'rejected step advanced physical clock')
    require(state_hash(rollback['state'])==state_hash(cp['coupled_state']),'rollback changed accepted coupled state')
    require([state_hash(x) for x in rollback['native_history']]==[state_hash(x) for x in cp['history']] and
            rollback['native_history_times']==cp['history_times'],'rollback changed native history')
    require('deliberate PIMPLE hard-audit rejection' in rollback['failure'],'wrong forced failure')
    require(state_hash(rollback['rejected_state'])==new['corrected_sha256'],'rejection saved wrong candidate')
    events=rollback['physical_range_events'];failure_event=events[-1]
    require(state_hash(failure_event['history'])==state_hash(q['history']),'rejection saved wrong callback history')
    out['checks']['activated_rejection_full_rollback_and_actual_candidate_history']=True
    segment=r/'repaired-campaign/trajectory/steps-000073-000073';native=read_report(segment)
    require(native['steps']==73 and native['failure']=='','native restart did not pass former failure')
    attempt=native['pimple_attempts'][0];packet=load(Path(attempt['k_snapshot']))
    require(state_hash(packet['raw'])==RAW_SHA and state_hash(packet['history'])==state_hash(q['history']),'native raw/history differs from archived attempt')
    require(state_hash(packet['pre_k'])==state_hash(q['pre_k']) and packet['transport_rows']==q['transport_rows'] and
            packet['matrix_sha256']==q['matrix_sha256'],'native pre-k/free rows/frozen transport A differs')
    require(state_hash(packet['corrected'])==new['corrected_sha256'],'native repair differs from diagnostic replay')
    require(native['initial_k_bound_check']['activations']==0,'startup bound changed accepted restart')
    require(state_hash(native['accepted_states'][0])==state_hash(cp['coupled_state']),'restart initial state changed')
    published=header(segment/'checkpoint.1',True)
    require(published['steps']==73 and published['time']==native['time'] and
            state_hash(published['coupled_state'])==state_hash(native['state']) and
            [state_hash(v) for v in published['history']]==[state_hash(v) for v in native['native_history']] and
            published['history_times']==native['native_history_times'],'native accepted publication differs from state/history')
    for key in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'):
        require(attempt[key]==1,'changed native governing schedule: '+key)
    require(attempt['additional_governing_equation_solves']==0 and native['implicit_callback_count']==1,'extra native callback/solve')
    require(abs(attempt['k_added_storage_residual']-out['frozen_accounting']['storage_alpha_M_delta_max'])<1e-12 and
            abs(attempt['k_corrected_frozen_linear_residual']-out['frozen_accounting']['corrected_linear_max'])<1e-12,'native corrected residual accounting differs')
    out['checks']['native_attempt73_raw_history_corrected_bitwise_matches_replay']=True
    out['checks']['native_restart_startup_no_op']=True
    out['checks']['native_pre_k_free_rows_transport_A_and_published_history_exact']=True
    samples=native['physical_samples'];w=np.array(samples['weights_2']);values=np.array(samples['values_2'])
    integral=float(w@values);field_norm=float(np.sqrt(w@(values*values)))
    out['step73_correction_scales']=dict(field_k_integral=integral,field_k_mass_norm=field_norm,
        correction_over_raw_step=b['mass_norm_correction']/b['mass_norm_raw_step'],
        correction_over_field_mass_norm=b['mass_norm_correction']/field_norm,
        correction_integral_over_field_integral=b['integral_total']/integral,
        signed_correction_integral_rate=b['integral_total']/.000625,
        inner_outer_near_positive_negative=b['integrals_inner_outer_near_positive_negative'])
    out['native_attempt73']=dict(time=native['time'],steps=native['steps'],attempt={k:v for k,v in attempt.items() if k!='k_bound'},
        startup=native['initial_k_bound_check'],signature=native['signature'],last_audit=native['state_audits'][-1])
    out['current_campaign_progress']=load(r/'repaired-campaign/trajectory/progress.json')
    corrector=load(r/'native-corrector-checks-h0/summary.json')
    require(corrector['actual_initial_sha256']==state_hash(baseline['accepted_states'][0]),'no-op fixture used final state instead of actual h0')
    out['native_corrector_checks']=dict(status=corrector['status'],actual_initial_sha256=corrector['actual_initial_sha256'])
    for path,role in ((r/'replay-total20000/archived-packet.json','immutable original failed packet'),
        (r/'native-inputs/checkpoint72/State.bin','native restart with all histories'),
        (r/'replay-total20000/projection-pool.json','complete accumulated finite pool and coherent primal/dual'),
        (Path(d['mass_file']),'actual consistent free k mass'),(Path(d['cuts_file']),'whitened finite constraints'),
        (Path(attempt['k_snapshot']),'native attempt73 raw/corrected/history/A packet'),
        (r/'qp-tests2.log','focused mass-QP and budget checks'),(r/'driver-tests2.log','ownership/recovery/prefix checks')):artifact(path,role)
    write_json(a.output,out);print(a.output)


if __name__=='__main__':main()
