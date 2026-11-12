#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact executed H2c evidence; heavy states/rows remain in the indexed run root."""
import argparse,json,math,hashlib,subprocess
from datetime import datetime,timezone,timedelta
from pathlib import Path
from run_cls_epawf_picard import load,sha,state_hash,require
from run_cao_tamura_budget_h2a import physical
from continue_cls_epawf_pimple import atomic

def small(x):
    if isinstance(x,dict):
        out={}
        for k,v in x.items():
            if k.endswith('evaluation_row') and isinstance(v,dict):out[k]=dict(geometry=v['geometry'],element=v['element'],content_sha256=hashlib.sha256(json.dumps(v,sort_keys=True).encode()).hexdigest())
            elif k not in ('physical_samples_initial','physical_samples_final','physical_samples'):out[k]=small(v)
        return out
    if isinstance(x,list):return [small(v) for v in x]
    return x

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--revision',default='working-tree');a=p.parse_args();r=a.root.resolve();a.output.mkdir(parents=True,exist_ok=False)
    def read(n):return load(r/n)
    def save(n,x):atomic(a.output/(n+'.json'),x)
    pre=read('preflight-final.json');save('calibration',pre);save('measurements',read('measurements.json'))
    screens=read('screen-v3/progress.json');adjusted=read('physical-height-primary-v1/progress.json');off=read('primary-uncached-v1/progress.json');off2=read('physical-height-uncached-v1/progress.json')
    require(screens['status']==adjusted['status']=='screen_completed','zero-step screens incomplete')
    for on,uncached in [(screens['results']['65-hill'],off['results']['65-hill']),(adjusted['results']['65-hill'],off2['results']['65-hill'])]:
        for key in ('state_sha256','audit','scalar_inventory','seed_constitutive','target','moving_foot_fixture'):
            require(on[key]==uncached[key],'cache changed '+key)
    for d in [screens,adjusted]:
        for k,v in d['results'].items():
            if k.endswith('-hill'):require(v['moving_foot_fixture']['passed'],'moving foot fixture failed')
    save('screens',small(dict(tangential=screens,physical_height=adjusted,cache_bitwise=True)))
    ana=read('analytical-witnesses-tangential-v2.json');ana2=read('analytical-witnesses-physical-height-v2.json')
    require(ana['status']==ana2['status']=='passed','analytical witness replay failed');save('witnesses',small(dict(tangential=ana,physical_height=ana2)))
    hist=read('historical-h2b-replay-v1/progress.json');old=hist['observations']['composite-hill']
    require(old['state_sha256']=='7adda6c42bd9d6e8fae1c5259055dbd60290888593c0fbc626d027b39b2203e4' and old['audit']['accepted_awf_hard_admissible']==[0],'historical H2b failure changed')
    native=read('native-seed-v1/hill-h0/case.evolution.json');raw=adjusted['results']['65-hill']
    require(native['steps']==native['time']==native['implicit_callback_count']==0 and native['failure'],'actual native initialization did not preserve rejection')
    require(state_hash(native['state'])==raw['state_sha256'],'native mapped seed differs from diagnostic seed')
    save('native',dict(capture_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),status='initialization_rejected',hill_steps=0,hill_attempted_time=0,
        state_sha256=raw['state_sha256'],report_sha256=sha(r/'native-seed-v1/hill-h0/case.evolution.json'),failure=native['failure'],
        reference_sha256=state_hash(native['fixed_reference']),defect_sha256=state_hash(native['fixed_defect']),
        history_sha256=[state_hash(v) for v in native['native_history']],history_times=native['native_history_times'],
        receipt=read('native-seed-v1/progress.json')['results']['hill-h0']['receipt'],
        conditional_8_16_pair='not executed: hard initial covariance gate failed',historical_H2b_preserved=hist,
        classification='chosen streamfunction initialization is outside covariance domain before projection; not a proof that every selected-model hill solution is inadmissible'))
    flat=read('flat-final-v2/progress.json');require(flat['status']=='completed_flat8','flat readback/identity validation incomplete')
    q=read('flat-v1/flat-h0/case.evolution.json');f=read('flat-v1/flat8/case.evolution.json')
    # Zero-step production reports have no refreshed physical-sample cache.
    # Both comparisons use the strict, unchanged native-checkpoint readbacks.
    qa=dict(q,physical_samples=flat['results']['flat-h0']['readback']['diagnostic']['physical_samples'])
    fa=dict(f,physical_samples=flat['results']['flat8']['readback']['diagnostic']['physical_samples'])
    diff=physical(qa,fa)
    H=.0065;K=pre['heights'][-1]['K0'];E=pre['u_tau']**3/(pre['kappa_CLS']*H);area=.8*.8
    normalization=dict(U=math.sqrt(area),W=math.sqrt(area),velocity=math.sqrt(area),k=K*math.sqrt(area),epsilon=E*math.sqrt(.8*(.8-H)),pi_inner=math.sqrt(.8*H),pi_outer=math.sqrt(.8*(.8-H)),pressure=math.sqrt(area),wall_stress=pre['u_tau']**2*math.sqrt(.8),wall_K=K*math.sqrt(.8),mean_traction=pre['u_tau']**2)
    norms={k:v/normalization[k] for k,v in diff.items()};passing={k:norms[k]<=1e-4 for k in ('velocity','k','epsilon','pressure')}
    readbacks={k:small(flat['results'][k]['readback']) for k in ('flat-h0','flat8')}
    maximum=max(max(x['diagnostic']['selected_spatial_block_norms'][:4]) for x in readbacks.values())
    attempts=[dict(time=t['time'],counts={k:t[k] for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')},solved_epsilon_retained_bitwise=t['solved_epsilon_retained_bitwise'],bound_active=t['k_bound'].get('activated',t['k_bound'].get('activations')),bound_signed_supply=t['k_bound'].get('integral_total')) for t in f['pimple_attempts']]
    require(all(t['k_bound']['activations']==0 and t['k_bound']['integral_total']==0 for t in f['pimple_attempts']),'unexpected physical correction supply')
    save('flat',dict(status=flat['status'],actual_time=f['time'],actual_steps=f['steps'],dt=f['dt'],tests=flat['tests'],
        initial_state_sha256=state_hash(q['state']),final_state_sha256=state_hash(f['state']),reference_sha256=state_hash(f['fixed_reference']),
        checkpoint=f['checkpoint'],checkpoint_sha256=sha(Path(f['checkpoint'])/'State.bin'),
        dimensional_drift=diff,normalization=normalization,normalization_definition='SI nominal U(Hhill)=1 held fixed; model K0 and epsilon(H_A); physical represented domains, absolute pressure',
        normalized_drift=norms,absolute_field_preservation_target=1e-4,field_target_pass=passing,
        selected_spatial_target=1e-3,maximum_first_four_spatial_blocks=maximum,spatial_target_pass=maximum<=1e-3,
        classification='hard-admissible bounded compatibility path; not a verified exact discrete equilibrium',
        initial_reference_errors=q['finite_inflow_reference']['flat_reference_errors'],final_reference_errors=f['finite_inflow_reference']['flat_reference_errors'],
        attempts=attempts,readbacks=readbacks,receipts={k:v['receipt'] for k,v in flat['results'].items()},
        boundary='reference-compatible both streamwise sides and top prescribed; separate from open hill free-slip/total-traction policy'))
    qp=read('new-inventory-projection-v1/result.json');check=qp['checks'];require(check['projection_only']['success'] and check['projection_only']['inventory_policy']=='AcceptedBudgetV2','new inventory projection failed')
    reg=read('regressions-final-v2/progress.json');require(reg['status']=='passed','old-mode regressions failed')
    final_runtime=read('final-runtime-verification.json');require(final_runtime['passed'],'final runtime restart differed')
    analytic=read('analytical-free-tests-v1/progress.json');require(analytic['status']=='completed','analytical free-test readback incomplete')
    save('verification',dict(final_runtime=final_runtime,analytical_free_tests=analytic,kernel_checks=(r/'kernel-checks.log').read_text(),old_mode_regressions=reg,
        moving_foot_fixture=adjusted['results']['65-hill']['moving_foot_fixture'],cache_on_off_state_audit_scalar_target_bitwise=True,
        actual_native_seed_matches_probe=True,old_H2b_seed_hash_retained=True,
        physical_mass_projection=check['projection_only'],projection_objective=check['projection_physical_mass_objective'],independent_GL32_objective=check['projection_independent_GL32_objective'],
        interface_actual_epsilon_rejection=check['interface_actual_epsilon_rejection'],zero_PDE_projection_fixture=True))
    timing={}
    for name in ('flat-h0','flat8','flat4','flat4-restart'):
        log=r/'flat-v1'/name/'usage.txt';timing[name]=dict(usage_wall_user_system_RSS_KiB=list(map(float,log.read_text().split())))
    timing['mapped_native_zero_steps']=read('native-seed-v1/progress.json')['results']['hill-h0']['receipt'];timing['mapped_physical_step_cost']='not measured: initialization rejected'
    save('cost',timing)
    source=r/'source';paths=sorted(set(list((source/'solvers/IncNavierStokesSolver/EquationSystems').glob('RANSWall*.cpp'))+list((source/'solvers/IncNavierStokesSolver/EquationSystems').glob('RANSWall*.h'))+list((source/'library/SolverUtils/Turbulence').glob('*'))+list((source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('*h2c*'))+[source/'solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp']))
    save('source-revisions',dict(base_revision='85ce59e95a2e44a5217216ae1842e248d8a24ec9',tested_revision=a.revision,
        source_hashes={str(v.relative_to(source)):sha(v) for v in paths if v.is_file()},final_solver=final_runtime['runtime'],final_probe=analytic['probe'],runtime=adjusted['runtime'],probe=adjusted['probe'],
        prototype_runtime=screens['runtime'],prototype_probe=screens['probe'],flat_runtime=flat['runtime'],
        data_snapshot_sha256=sha(source/'docs/plans/cao-tamura-validation-data/snapshot.json'),
        resource_contract='isolated build parallel 4; at most four single-rank native jobs simultaneously; all computational library threads 1; total below 72',
        actual_run_root=str(r),data_volume_UUID='5AB8A669B8A642FD',actual_mount='/media/jay/新加卷1'))
    files=[]
    for x in sorted(r.rglob('*')):
        if not x.is_file() or x.is_symlink() or not x.relative_to(r).parts or x.relative_to(r).parts[0] in ('source','build') or a.output in x.parents:continue
        files.append(dict(path=str(x.relative_to(r)),size=x.stat().st_size,sha256=sha(x)))
    save('artifact-index',dict(root=str(r),files=files,total_bytes=sum(v['size'] for v in files),scope='immutable native payloads, linked-library closures, sessions, receipts and logs; build/source are identified separately'))
if __name__=='__main__':main()
