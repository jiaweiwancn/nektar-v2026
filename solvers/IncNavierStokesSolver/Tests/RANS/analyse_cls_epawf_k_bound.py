#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read retained states, replay physical samples, and audit bounded experiments.

This performs no governing-equation solve. Native readbacks and complete
coefficient/matrix packets remain under the Linux experiment directory.
"""
import argparse
import math
from pathlib import Path
import shutil

from run_cls_epawf_picard import load, require, run, sha, state_hash, write_json
from run_cls_epawf_pimple import identity
from compare_cls_epawf_pimple import describe, differences, same_history
from test_cls_epawf_k_bound import replay


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiments', type=Path, required=True)
    p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.experiments = a.experiments.resolve(); a.probe = a.probe.resolve(strict=True)
    a.output = a.output.resolve(); a.output.mkdir(parents=True, exist_ok=False)
    reports = {x.parent.name: load(x) for x in a.experiments.glob('*/case.evolution.json')}
    summary = dict(schema=1, experiments=str(a.experiments), probe=identity(a.probe), runs={},
                   prefix_checks={}, readbacks={}, refinement={}, status='running')
    counts = ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves',
              'epsilon_transport_solves', 'k_transport_solves', 'additional_governing_equation_solves',
              'mass_inverse_applications', 'bound_calls', 'bound_activations', 'bound_active_set_iterations',
              'bound_mass_factor_applications', 'bound_kkt_solves')
    detail = ('time', 'alpha', 'governing_solves_complete', 'scheduled_complete', 'failure',
              'corrected_continuity_gauge', 'boundary_defect', 'imposed_trace_defect', 'refreshed_trace_defect',
              'k_frozen_linear_residual', 'k_corrected_frozen_linear_residual', 'k_added_frozen_residual',
              'k_added_storage_residual', 'refreshed_original_block_norms', 'history_sha256', 'k_snapshot',
              'k_snapshot_sha256', 'k_budget_last_accepted', 'k_budget_post_flow', 'k_budget_pre_k')
    for name, r in reports.items():
        m = load(a.experiments/name/'manifest.json')
        require(r['implicit_callback_count'] == len(r['pimple_attempts']), 'callback count differs from attempts')
        accepted = set(r['audit_times'][1:]); attempts = []; integral = 0.
        for row in r['pimple_attempts']:
            item = {k: row[k] for k in (*counts, *detail) if k in row}
            item['accepted'] = row['time'] in accepted
            if item['accepted']:
                require(all(row[k] == 1 for k in counts[:5]) and row['additional_governing_equation_solves'] == 0,
                        'accepted step changed the governing solve schedule')
            for phase in ('k_domain', 'k_ground'):
                if phase in row: item[phase] = row[phase]
            if row.get('k_snapshot'):
                packet = load(Path(row['k_snapshot']))
                item['raw_sha256'] = state_hash(packet['raw']); item['corrected_sha256'] = state_hash(packet['corrected'])
                first, last = int(r['offsets'][2]), int(r['offsets'][3])
                require(state_hash(packet['raw'][:first]+packet['raw'][last:]) ==
                        state_hash(packet['corrected'][:first]+packet['corrected'][last:]), 'production bound changed a non-k field')
                require(state_hash(packet['history']) == row['history_sha256'], 'snapshot history changed')
                require(sha(Path(row['k_snapshot'])) == row['k_snapshot_sha256'], 'snapshot changed')
                require(sha(Path(packet['matrix_file'])) == packet['matrix_sha256'], 'frozen matrix changed')
            if 'k_bound' in row:
                bound = dict(row['k_bound'])
                bound['active_constraint_count'] = len(bound['active_constraints'])
                bound['minimum_normalized_dual'] = min((v['normalized_dual_multiplier'] for v in bound['active_constraints']), default=0.)
                bound['active_locations'] = [dict(coordinate=v['coordinate'], geometry=v['geometry'], inventory=v['inventory'])
                                             for v in bound.pop('active_constraints')]
                item['bound'] = bound
                if item['accepted']: integral += bound['integral_total']
            item['cumulative_accepted_bound_integral'] = integral; attempts.append(item)
        summary['runs'][name] = dict(status=m['status'], steps=r['steps'], time=r['time'], dt=r['dt'],
            requested_steps=m['requested_steps'], seconds=m['execution']['seconds'], failure=r['failure'],
            state_sha256=state_hash(r['state']), history_sha256=[state_hash(v) for v in r['native_history']],
            history_times=r['native_history_times'], report_sha256=sha(a.experiments/name/'case.evolution.json'),
            manifest_sha256=sha(a.experiments/name/'manifest.json'), executable_sha256=m['executable']['sha256'],
            initial_check=r.get('initial_k_bound_check'), initialization=r['initialization_lift_calls_iterations'],
            signature=r['k_bound_signature'], implicit_callbacks=r['implicit_callback_count'],
            totals={k: sum(x.get(k, 0) for x in attempts) for k in counts},
            accepted_bound_integral=integral, attempts=attempts,
            final_audit={k: v for k, v in r['state_audits'][-1].items()
                         if k.startswith('r2_') or k in ('strong_continuity', 'robust_source', 'inner_stress')},
            qualification_events=[{k: v for k, v in e.items() if k not in ('state', 'history', 'audit')}
                                  for e in r['physical_range_events']])
    for i in range(3):
        off, bound = reports[f'off-{i}'], reports[f'bound-{i}']
        count = min(len(off['accepted_states']), len(bound['accepted_states']))
        equal = 0
        for x, y in zip(off['accepted_states'][:count], bound['accepted_states'][:count]):
            if state_hash(x) != state_hash(y): break
            equal += 1
        require(equal == len(off['accepted_states']), 'On changed an uncorrected accepted prefix')
        summary['prefix_checks'][str(i)] = dict(bitwise_equal_accepted_states=equal, off_time=off['time'])
        for x, y in zip(off['pimple_attempts'], bound['pimple_attempts']):
            if x.get('k_snapshot') and y.get('k_snapshot') and x['time'] <= off['time'] + off['dt']:
                px, py = load(Path(x['k_snapshot'])), load(Path(y['k_snapshot']))
                require(state_hash(px['history']) == state_hash(py['history']), 'incoming history changed before first activation')
                require(state_hash(px['raw']) == state_hash(py['raw']), 'raw transport changed before first activation')
        summary['prefix_checks'][str(i)]['raw_and_history_identical_through_first_rejection'] = True
    require(same_history(reports['inactive-on'], reports['inactive-off']), 'inactive On changed native transaction')
    summary['inactive_on_off_bitwise_state_time_history'] = True
    if 'restart' in reports:
        require(same_history(reports['bound-0'], reports['restart']), 'active restart differs')
        summary['active_restart_bitwise_state_time_all_native_history'] = True
    rejected, off = reports['hard-reject'], reports['off-0']
    require(same_history(rejected, off), 'active repair rejection changed committed transaction')
    event = next(e for e in reversed(rejected['physical_range_events'])
                 if 'PIMPLE_forced_hard_audit_rejection_test' in e['reasons'])
    packet = load(Path(rejected['pimple_attempts'][-1]['k_snapshot']))
    require(state_hash(event['state']) == state_hash(packet['corrected']) == state_hash(rejected['rejected_state']), 'wrong rejected candidate')
    require(state_hash(event['history']) == state_hash(packet['history']), 'wrong rejected incoming history')
    require(bool(event['audit']) and packet['k_bound']['activations'] == 1, 'audit missing or repair inactive')
    summary['hard_rejection'] = dict(committed_state_time_all_native_history_restored=True,
        actual_repaired_candidate_and_incoming_history_saved=True, audit_saved=True, reason=event['reasons'],
        rejected_sha256=state_hash(event['state']), incoming_history_sha256=state_hash(event['history']))

    def prefix(name, time):
        report = reports[name]
        idx = next(i for i, t in enumerate(report['audit_times']) if abs(t-time) < 1e-13)
        row = next(v for v in report['pimple_attempts'] if abs(v['time']-time) < 1e-13)
        packet = load(Path(row['k_snapshot'])); state = report['accepted_states'][idx]
        require(state_hash(state) == state_hash(packet['corrected']), 'snapshot is not accepted prefix')
        folder = a.output/(name+'-t'+str(time)); folder.mkdir()
        shutil.copy2(a.experiments/name/'case.xml', folder/'case.xml')
        write_json(folder/'input.json', dict(operation='samples', state=state, history=packet['history'],
            alpha=packet['alpha'], time=time, require_history_h0=False, inner_stress_model='SaturatedR2'))
        receipt = run([str(a.probe), 'case.xml', 'input.json', 'result.json'], folder, 'probe.log', 900)
        require(receipt['exit_code'] == 0, 'physical prefix replay failed')
        result = load(folder/'result.json'); result['initial'] = report['initial']
        short = dict(report, time=time, steps=idx, state_audits=report['state_audits'][:idx+1],
                     pimple_attempts=[v for v in report['pimple_attempts'] if v['time'] <= time])
        description = describe(short); description['ground_minimum'] = row['k_ground']['minimum']
        description['bound_integral'] = sum(v.get('k_bound', {}).get('integral_total', 0.) for v in short['pimple_attempts'])
        description['activations'] = sum(v.get('bound_activations', 0) for v in short['pimple_attempts'])
        summary['readbacks'][folder.name] = dict(receipt=receipt, state_sha256=state_hash(state), description=description,
                                               result_sha256=sha(folder/'result.json'))
        return result

    for time in (.0025, .005):
        states = [prefix(f'bound-{i}', time) for i in range(3)]
        delta = [differences(x, y) for x, y in zip(states, states[1:])]
        summary['refinement'][str(time)] = dict(physical_L2_differences=delta,
            successive_difference_orders={k: math.log2(delta[0][k]/delta[1][k]) for k in delta[0]
                if k != 'pressure_common_gauge_shift' and delta[0][k] > 0 and delta[1][k] > 0})
    summary['floor_trajectory_differences'] = {name: differences(reports[name], reports['bound-0'])
                                             for name in ('floor-low', 'floor-high')}
    raw_packet = load(Path(reports['off-0']['pimple_attempts'][-1]['k_snapshot']))
    summary['floor_replays'] = {}
    for floor in (1e-10, 1e-8, 1e-6):
        result = replay(a.probe, a.experiments/'off-0/case.xml', raw_packet, a.output/('floor-'+str(floor)), floor=floor)
        summary['floor_replays'][str(floor)] = result['k_bound']
    summary['status'] = 'passed'; write_json(a.output/'summary.json', summary)
    print(a.output/'summary.json')


if __name__ == '__main__':
    main()
