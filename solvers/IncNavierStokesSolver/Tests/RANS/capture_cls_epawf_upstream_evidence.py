#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact evidence index; large native vectors/matrices stay in the run directory."""
import argparse
import collections
import datetime
from pathlib import Path
from run_cls_epawf_picard import load,write_json,sha,state_hash,require

def compact_status(status):
    """Keep process/checkpoint identities; index the full durable ledger on Linux."""
    result={k:v for k,v in status.items() if k!='cases'};result['cases']={}
    for name,record in status['cases'].items():
        copy=dict(record);progress=dict(copy['progress'])
        segments=progress.pop('segments',[])
        progress['completed_segment_count']=len(segments)
        progress['latest_completed_segment']=segments[-1] if segments else None
        copy['progress']=progress;result['cases'][name]=copy
    return result

def compact_observation(q):
    d=q['diagnostic']
    return {**{k:v for k,v in q.items() if k!='diagnostic'},'diagnostic':{
        k:v for k,v in d.items() if k not in ('profiles','native_scalar_domain')}}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--status',type=Path);a=p.parse_args()
    r=a.root.resolve();a.output.mkdir(parents=True,exist_ok=True);index={}
    def artifact(path,role):
        path=path.resolve();record=dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size,role=role)
        index[str(path)]=record;return record
    comparison=load(r/'comparison.json');budgets=load(r/'discrete-budget/summary.json')
    out=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),experiment=str(r),
        base='de3e0a1d6d9e24860a1957a1baa9505a97f47ec7',production_source=load(r/'screening/progress.json')['source'],
        policy='ConstitutiveV3: Q0 union intrinsically permitted Qpositive, current profiles, means, local/retained rows, continuous ground',
        positivity_scope='declared constitutive evaluation inventory and continuous cubic ground, not global volume positivity',
        temporal_contract='native IMEX1; one vector U, one pressure/gauge, one epsilon, one k and one logical bound call; no extra governing solves',
        cases={},paired=comparison['paired'],policy_comparisons=comparison['policy_comparisons'],source_balances={},runtime={},software={})
    for c,v in comparison['cases'].items():
        path=r/'screening'/c;native=load(path/'trajectory/case.evolution.json');result=load(path/'result.json')
        initial=load(path/'initial/case.evolution.json');attempts=native['pimple_attempts'];last=attempts[-1];b=last['k_bound']
        work=dict(v['work']);work.pop('increments',None)
        family=collections.Counter()
        for x in attempts:family.update(x['k_bound']['discovered_constraint_families'])
        health={k:max(abs(x['k_bound'].get(k,0.)) for x in attempts) for k in
            ['primal_floor_defect','ground_floor_defect','stationarity_mass_dual','complementarity','dual_feasibility_defect','primal_dual_gap_working_set']}
        health['maximum_continuity_gauge']=max(max(map(abs,x['corrected_continuity_gauge'])) for x in attempts)
        health['maximum_raw_k_linear_residual']=max(x['k_frozen_linear_residual'] for x in attempts)
        events=[{k:x[k] for k in ('time','alpha','reasons','stopped_before_history_advance') if k in x} for x in native['physical_range_events']]
        observations={}
        for n,q in v['observations'].items():
            d=q['diagnostic'];observations[n]=dict(time=q['time'],drift=q['drift_from_actual_h0'],error_to_own_reference=q['error_to_own_reference'],
                spatial_dimensional=d['dimensional_mass_dual_spatial_first4'],spatial_common_normalized=d['common_normalized_spatial_first4'],
                mean_wall_stress=d['mean_wall_stress'],bands=d['fixed_bands']['regions'])
            detail=r/'observations'/c/f'step-{int(n):06d}'/'height_comparison/result.json'
            if detail.exists():observations[n]['matching_trace_defect']=load(detail)['matching_trace_defect']
        budget=budgets[c];maxima=dict(closure=max(abs(t) for s in budget['steps'] for t in s['increments']['closure']),
            corrected_frozen_identity=max(s['corrected_identity_defect'] for s in budget['steps']))
        out['cases'][c]=dict(definition=v['definition'],status=v['status'],accepted_steps=v['steps'],accepted_time=v['time'],
            h0_sha256=result['actual_h0_sha256'],h0_native_history_hashes=result['initial_native_history_hashes'],
            initialization_work=initial['initialization_lift_calls_iterations'],initialization_seconds=result['initial_receipt']['seconds'],
            initial_error=v['initial_error_to_own_reference'],initial_reference=v['initialization']['reference'],
            initial_wall_stress=v['initialization']['mean_wall_stress'],initial_representation_errors=v['initialization']['representation_L2_velocity_gradient_production_referenceProduction'],
            degrees_of_freedom=v['initialization']['native_degrees_of_freedom'],observations=observations,
            trajectory_wall_seconds=result['trajectory_receipt']['seconds'],work=work,health=health,discovered_cut_families=dict(family),
            final_ground=b['corrected_ground'],final_inventory={k:last['k_domain'][k] for k in ['minimum_inner_k','minimum_outer_k','minimum_outer_epsilon','minimum_sampled_K','minimum_b','minimum_lambda']},
            qualification=dict(satisfied_at_final=native['state_audits'][-1]['r2_qualified'][0]==1,first_exit=events[0] if events else None,event_count=len(events)),
            frozen_budget=dict(bands=budget['bands'],totals=budget['totals'],endpoint_integral_change=budget['endpoint_integral_change'],
                initialized_projected_band_rates=budget['initialized_projected_band_rates'],max_defects=maxima))
        for f,role in [(path/'initial/case.xml','fresh physical session'),(path/'initial/case.evolution.json','actual native h0 and initialization histories'),
            (path/'trajectory/case.evolution.json','accepted trajectory and audits'),(Path(native['checkpoint'])/'State.bin','native accepted step-32 checkpoint')]:artifact(f,role)
        artifact(Path(initial['checkpoint'])/'State.bin','native initial history checkpoint')
        obs=r/'observations'/c
        out['source_balances'][c]={}
        for n in ('0','8','16','32'):
            f=obs/f'step-{int(n):06d}'/'source_balance/result.json'
            if not f.exists():continue
            d=load(f);points=[x for x in d['points'] if x.get('valid') and x.get('z')==0]
            keys=['x','variant','chi','spatial_gate','state_gate','admissibility_gate','margin_gate','reynolds_gate','sensor_d2',
                'local_k','sampled_K','production','reference_production','epsilon_A','D_pre','D','source','B','m','lambda',
                'saturation_active','linear_production_same_gradient','analytical_gradient_saturated_production','selected_native_spatial_rate']
            failures=[x for x in d['points'] if not x.get('valid')]
            out['source_balances'][c][n]=dict(ground=[{k:x.get(k) for k in keys} for x in points],
                invalid_extra_queries=failures,budget_identity_max=max(abs(q['budget_identity_defect']) for q in d['ray_budgets']),
                component_load_sum_defect=d['component_load_sum_defect'],native_mass_solve_defect=d['native_mass_solve_defect'])
            artifact(f,'read-only source/gate/counterfactual data; actual retained epsilon')
            branch=obs/f'step-{int(n):06d}'/'branch_coverage/result.json'
            if branch.exists():
                b=load(branch)['state'];out['source_balances'][c][n]['independent_vertical_minima']=[
                    dict(x=q['anchor'][0],**min(q['vertical_cubic_panels'],key=lambda v:v['minimum'])) for q in b['anchors']]
                artifact(branch,'independent continuous vertical cubic minima; finite-inventory scope diagnostic')
    for n in ['screening/runtime','screening/probe-runtime','detailed-probe-runtime','reviewed-validation-runtime','reviewed-probe-runtime','branch-entry-probe-runtime']:
        f=r/n/'identity.json';out['runtime'][n]=load(f)['frozen'];artifact(f,'frozen executable/library closure')
    for n in ['focused','active-followup','active-followup-detail','reviewed-followup']:
        f=r/n/'summary.json'
        if f.exists():
            d=load(f);out['software'][n]={k:v for k,v in d.items() if isinstance(v,bool)};artifact(f,'focused native verification')
            if 'strict_reader_rejections' in d:
                out['software']['strict_checkpoint_reader_rejections']={k:{x:q[x] for x in ('exit_code','input_checkpoint_unchanged')}
                    for k,q in d['strict_reader_rejections'].items()}
    for n in ['focused-ctest.log','ctest-reviewed.log']:
        f=r/n
        if f.exists():out['software'][n]=f.read_text();artifact(f,'focused Linux CTest results')
    f=r/'upstream-quadrature/summary.json';q=load(f);out['quadrature']=q;artifact(f,'frozen source integration order 16/32, all free modes')
    forensic=next(r/n for n in ('forensic-v3-entry','forensic-final','forensic103-native-history') if (r/n/'result.json').exists())
    f=load(forensic/'result.json');before=f['state']['anchors'][0];after=f['compare']['anchors'][0]
    def branch(a):return {k:a[k] for k in ['K','chi','spatial','state','d2','margin_gate','reynolds_gate','margin','reynolds','prospective_reason']}
    out['forensic']=dict(context=f['native_replay_context'],k_coefficient_difference=f['k_coefficient_difference'],
        before=branch(before),after=branch(after),minimum_vertical=min(before['vertical_cubic_panels'],key=lambda q:q['minimum']),
        baseline_minimum=before['profiles'][0]['minimum'],first_missing_negative=next(x for p in before['profiles'] if p['branch']=='prospective_active' for x in p['rows'] if x['k']<0))
    artifact(forensic/'result.json','strict read-only saved-state branch replay');artifact(forensic/'input.json','exact archived candidate and native history-formation request')
    out['forensic']['selected_inventory_before']=f['state']['selected_inventory']
    out['forensic']['selected_inventory_after']=f['compare']['selected_inventory']
    if 'prospective_inventory_only' in f:out['forensic']['v3_read_only_inventory']=f['prospective_inventory_only']
    for f in [r/'comparison.json',r/'discrete-budget/summary.json',r/'plots/plot-validation.json']:
        artifact(f,'physical norms, exact frozen budgets, or native plot verification')
    terminal=r/'terminal-U05-graded/manifest.json'
    if terminal.exists():
        t=load(terminal);b=t.pop('bound');t.pop('files',None)
        keys=('failure','success','activations','floor','arithmetic_reserve','raw_minimum','corrected_minimum',
            'raw_ground','corrected_ground','mass_norm_correction','mass_norm_raw_step','correction_to_step_ratio',
            'primal_floor_defect','ground_floor_defect','stationarity_mass_dual','complementarity','dual_feasibility_defect',
            'primal_dual_gap_working_set','current_inventory_passed','separation_passes','inventory_rebuilds',
            'active_set_iterations','kkt_solves','mass_factor_applications','maximum_work_one_add','maximum_active_size',
            'seconds_total','seconds_algebra','seconds_inventory','unfinished_working_set','bounded_trace',
            'last_completed_qr_diagonal_range','minimum_observed_qr_residual','finite_pool_normalized_violation',
            'completed_primal_dual_coherent','primal_objective','dual_objective','repeated_working_sets','dependent_entries',
            'unique_evaluation_rows','repeated_rows','repeated_identical_constraints','rediscovered_identical_constraints')
        t['bound']={k:b[k] for k in keys if k in b};out['continuation_terminal']=t
        for f in (terminal,terminal.with_name('attempt.json'),terminal.with_name('replay-input.json'),
                  terminal.parent/'checkpoint/State.bin',terminal.with_name('transport_A.bin')):
            artifact(f,'preserved terminal attempt; diagnostic replay prepared, not executed')
    monitor=r/'pair-monitor/summary.json'
    if monitor.exists():
        m=load(monitor);out['continuation_windows']={}
        for c,v in m['cases'].items():
            work=dict(v['work']);work.pop('increments',None)
            out['continuation_windows'][c]={**{k:v[k] for k in ('status','accepted_steps','accepted_time','completed_targets',
                'next_target','failure','initial_bulk_u','flow_through_time','H_over_ustar','fixed_time_windows')},'work':work,
                'observations':{n:compact_observation(q) for n,q in v['observations'].items() if int(float(n)) in (64,80,160) or int(float(n))%80==0}}
        artifact(monitor,'read-only fixed physical-window monitor; no governing solves')
    uniform=r.parent/'cls-epawf-height-sampler-20261107/pair-monitor/summary.json'
    if uniform.exists():
        m=load(uniform);out['existing_uniform_observations']={}
        for c,v in m['cases'].items():
            latest=sorted(v['observations'],key=float)[-2:]
            out['existing_uniform_observations'][c]={k:v[k] for k in ('status','accepted_steps','accepted_time','flow_through_time','H_over_ustar')}
            out['existing_uniform_observations'][c].update(observations={n:compact_observation(v['observations'][n]) for n in latest},
                fixed_time_windows=v['fixed_time_windows'][-3:],activations=v['work']['activations'],correction_integral=v['work']['correction_integral'])
        artifact(uniform,'read-only existing Uniform G/H context; frozen runtime and physical policy unchanged')
    for n in ('source-validation.json','probe-entry-build.log','build-reviewed.log','reviewed-followup.log'):
        f=r/n
        if f.exists():artifact(f,'tested source content or focused Linux build/verification log')
    f=r/'continuation-budgets.json'
    if f.exists():out['continuation_budgets']=load(f);artifact(f,'accepted-step-only cumulative and fixed-window correction rates')
    if a.status:
        out['latest_status']=compact_status(load(a.status));artifact(a.status,'dated full live process/ledger/checkpoint snapshot')
        write_json(a.output/'latest-status.json',out['latest_status'])
        for f in a.status.parent.glob('*/checkpoint/State.bin'):artifact(f,'immutable accepted frontier copy from this dated observation')
    out['historical_terminal_checkpoints']={}
    for n,f,expected in (
        ('72',r.parent/'cls-epawf-uniform-ep-flat-20261106/upstream-terminal/checkpoint-last-accepted/State.bin',
            '0ebcf803b486d63c42d6ed93db08efd4208d87b98e90dd918e139a77f8c26510'),
        ('103',r.parent/'cls-epawf-height-sampler-20261107/status-final/upstream-checkpoint/State.bin',
            '11d1c23373ed98a3f79c7ce6f242c23883398a3e369e1210e34a81a1c42445e9')):
        require(sha(f)==expected,'historical terminal checkpoint changed')
        out['historical_terminal_checkpoints'][n]=artifact(f,'unchanged historical terminal checkpoint; no new old-tail continuation')
    write_json(a.output/'native.json',out);write_json(a.output/'artifact-index.json',dict(root=str(r),artifacts=list(index.values()),
        per_step_packets='screening/<case>/trajectory/pimple-states/*.json: recoverable raw/corrected/pre-k/history plus A_k path/hash; A_k is not M_ff',
        mass='discrete-budget/<case>/mass-and-bands.json: native PicardKMass and identical free rows',
        continuations='selected-U05-{graded,one}/trajectory/progress.json, steps.jsonl, runtime/identity.json, steps-*/receipt.json'))
    print(a.output/'native.json')

if __name__=='__main__':main()
