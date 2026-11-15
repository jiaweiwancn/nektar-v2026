#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict post-startup readback and physical norms for the open mapped hill.

Pressure is absolute: no separate regional or common-mean gauge is removed.
The readback uses saved native callback histories for physical storage. It
performs no physical solves, Lift, field transfer or fresh seed qualification.
"""
import argparse
from collections import defaultdict
from datetime import datetime
import math
from pathlib import Path
from zoneinfo import ZoneInfo

from cao_tamura_h2d_state import load
from compare_cls_epawf_pimple import same_history
from diagnose_cls_accepted_covariance import checkpoint, common_packet, run_probe
from run_cls_epawf_picard import require, sha, state_hash, write_json

FIELDS = ('U', 'W', 'k', 'outer_epsilon', 'inner_pressure', 'outer_pressure')


def norm(weights, values):
    require(len(weights) == len(values) and all(w > 0 for w in weights),
            'invalid physical integration weights')
    return math.sqrt(math.fsum(w*v*v for w, v in zip(weights, values)))


def difference(a, b, scales):
    result = {}
    for i, field in enumerate(FIELDS):
        w = a['weights_'+str(i)]
        require(state_hash(w) == state_hash(b['weights_'+str(i)]), 'quadrature changed')
        d = [x-y for x, y in zip(a['values_'+str(i)], b['values_'+str(i)])]
        raw = norm(w, d)
        result[field] = dict(physical_L2=raw, physical_RMS=raw/math.sqrt(math.fsum(w)),
                             fixed_scale_normalized_RMS=raw/math.sqrt(math.fsum(w))/scales[field])
    result['velocity_L2'] = math.hypot(result['U']['physical_L2'], result['W']['physical_L2'])
    result['absolute_pressure_L2'] = math.hypot(result['inner_pressure']['physical_L2'],
                                               result['outer_pressure']['physical_L2'])
    w = a['wall_weights']
    require(state_hash(w) == state_hash(b['wall_weights']), 'wall quadrature changed')
    for name in ('wall_stress', 'wall_K'):
        result[name+'_L2'] = norm(w, [x-y for x, y in zip(a[name], b[name])])
    return result


def timing(folder):
    rows = [__import__('json').loads(s) for s in (folder/'timing.jsonl').read_text().splitlines()]
    all_cost, continued = defaultdict(lambda: [0, 0., 0., 0., 0.]), defaultdict(lambda: [0, 0., 0., 0., 0.])
    for row in rows:
        for target in (all_cost, continued) if row['step'] >= 0 else (all_cost,):
            v = target[row['name']]
            v[0] += 1
            for j, key in enumerate(('wall_seconds', 'process_cpu_seconds',
                                     'exclusive_wall_seconds', 'exclusive_process_cpu_seconds'), 1):
                v[j] += row[key]
    return dict(columns=['count', 'inclusive_wall', 'inclusive_cpu', 'exclusive_wall', 'exclusive_cpu'],
                all=dict(all_cost), physical_step_context=dict(continued),
                definition='Exclusive times partition work; inclusive ancestors overlap and must not be added.')


def summarize(report, record, samples):
    audit = report['state_audits'][-1]
    rows = report['pimple_attempts']
    counter_names = ('outer_sweeps', 'momentum_predictor_solves', 'pressure_system_solves',
                     'epsilon_transport_solves', 'k_transport_solves', 'bound_calls',
                     'bound_activations', 'additional_governing_equation_solves',
                     'additional_correction_or_restoration_solves')
    counts = {key: sum(row[key] for row in rows) for key in counter_names}
    for row in rows:
        require(row['scheduled_complete'] and not row['failure'], 'unexpected physical rejection')
        require(all(row[k] == 1 for k in counter_names[:6]), 'changed one-sweep schedule')
        require(row['solved_epsilon_retained_bitwise'], 'solved epsilon overwritten')
    interface = audit['interface']
    selected = {key: audit[key] for key in (
        'accepted_awf_covariance', 'accepted_awf_inner_minimum_physical_point',
        'accepted_awf_outer_minimum_physical_point', 'accepted_awf_scalar_domain',
        'accepted_awf_sources', 'strong_continuity', 'coupled_norms',
        'native_inner_minimum_geometry_element_side_block',
        'native_outer_minimum_geometry_element_side_block')}
    stage_names = ('momentum_splitting_identity_defect', 'frozen_momentum_defect',
                   'corrected_continuity_gauge', 'boundary_defect', 'native_full_G_action_defect',
                   'imposed_trace_defect', 'refreshed_trace_defect',
                   'epsilon_same_state_source_storage_defect', 'k_frozen_linear_residual')
    final = {key: rows[-1].get(key) for key in stage_names} if rows else {}
    kb = [row['k_bound'] for row in rows]
    weight = samples['wall_weights']
    wall = {name: dict(minimum=min(samples[name]), maximum=max(samples[name]),
                       physical_mean=math.fsum(w*v for w, v in zip(weight, samples[name]))/math.fsum(weight))
            for name in ('wall_stress', 'wall_K')}
    return dict(time=report['time'], dt=report['dt'], global_accepted_steps=report['steps'],
                local_accepted_attempts=len(rows), counts=counts, final_audit=selected,
                named_interface=dict(normal_traction_jump_RMS=interface[2],
                                     full_traction_jump_RMS=interface[3],
                                     maximum_absolute_normal_jump=interface[4],
                                     k_flux_jump_RMS=interface[5]),
                wall=wall, final_stage=final,
                correction=dict(activations=sum(v['activations'] for v in kb),
                    signed_integral_sum=math.fsum(v['integral_total'] for v in kb),
                    region_near_positive_negative_integrals=[v['integrals_inner_outer_near_positive_negative'] for v in kb],
                    maximum_mass_correction=max([v['mass_norm_correction'] for v in kb] or [0.])),
                solver_driver_wall_seconds=record['seconds'], process_usage=record['usage'],
                implicit_callback_seconds=math.fsum(row['seconds'] for row in rows),
                checkpoint=report['checkpoint'], checkpoint_sha256=sha(Path(report['checkpoint'])/'State.bin'),
                state_sha256=state_hash(report['state']), history_sha256=[state_hash(v) for v in report['native_history']],
                history_times=report['native_history_times'],
                source_report=record['report'], source_report_sha256=sha(Path(record['report'])))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('root', 'probe', 'comparison', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(exist_ok=False)
    progress = load(a.root/'pilot-requalified/progress.json')
    require(progress['status'] == 'completed_bounded_pair_revalidated', 'startup protocol not qualified')
    records = progress['results']
    reports = {key: load(Path(records[key]['report'])) for key in
               ('coarse-h0', 'coarse8', 'fine-h0', 'fine16', 'coarse4', 'coarse4-restart')}
    require(same_history(reports['coarse8'], reports['coarse4-restart']), 'split histories differ')
    require(state_hash(reports['coarse-h0']['state']) == state_hash(reports['fine-h0']['state']), 'different seeds')
    # Native repeated addition need not give identical last bits for 8 versus
    # 16 steps. Preserve the actual labels; do not relabel/interpolate states.
    for label in ('coarse8', 'fine16'):
        require(abs(reports[label]['time']-1.6e-5) <= 32*math.ulp(1.6e-5), 'different physical endpoint')
    preflight = load(a.comparison)
    standard = next(v for v in preflight['models'] if v['name'] == 'StandardKEpsilon')
    U = preflight['u_tau']/preflight['kappa']*math.log(.04/preflight['z0_eff'])
    scales = dict(U=U, W=U, k=standard['K0'], outer_epsilon=standard['E_C_H'],
                  inner_pressure=U*U, outer_pressure=U*U)
    result = dict(version='H2kPostStartupAnalysisV1',
                  created_Asia_Shanghai=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),
                  probe=str(a.probe), probe_sha256=sha(a.probe), fixed_normalization=scales,
                  pressure='absolute mechanical pressure in both regions; no gauge alignment',
                  endpoint_time_labels={key: reports[key]['time'] for key in ('coarse8', 'fine16')},
                  endpoint_arithmetic_tolerance=32*math.ulp(1.6e-5),
                  norms='same native physical-volume quadrature with all represented coefficients; complete outer epsilon',
                  governing_solves=0, physical_commits=0, runs={}, readbacks={})
    result['limitations'] = 'Two levels give adjacent-step differences, not an observed order. Startup is not developed flow or experimental validation. Weak continuity does not imply pointwise zero divergence.'
    write_json(a.output/'analysis.json', result)
    native_samples = {}
    for label in ('coarse-h0', 'coarse8', 'fine16'):
        report = reports[label]
        cp = Path(report['checkpoint'])
        before = sha(cp/'State.bin')
        saved = checkpoint(cp)
        packet = dict(common_packet(cp), operation='h2c_checkpoint', inner_stress_model='Linear',
                      state=saved['coupled'], history=saved['coupled'], alpha=0., time=saved['time'],
                      native_report_file=records[label]['report'], native_report_sha256=sha(Path(records[label]['report'])))
        data, receipt = run_probe(a.probe, Path(records[label]['report']).parent/'case.xml', packet, a.output/(label+'-readback'))
        require(data['governing_solves'] == data['physical_commits'] == 0 and data['initialization_work'] == [0., 0.] and
                data['selected_state_unchanged'] and not data['native_history_advanced'], 'readback changed native state or advanced')
        require(before == sha(cp/'State.bin'), 'readback changed checkpoint bytes')
        require(data['state_sha256'] == state_hash(report['state']), 'readback selected another state')
        native_samples[label] = data['physical_samples']
        result['readbacks'][label] = dict(receipt=receipt, checkpoint_sha256=before,
            native_turbulence_budget=data['native_turbulence_budget'], developing_boundary=data['developing_boundary'],
            audit=data['audit'], native_divergence_GL16=data['native_divergence_GL16'],
            selected_spatial_block_norms=data['selected_spatial_block_norms'],
            native_replay_context=data['native_replay_context'],
            no_Lift_no_solve_no_advance_checkpoint_unchanged=True)
        write_json(a.output/'analysis.json', result)
    # A zero-step import does not populate the evolution reporter's samples.
    # Its strict checkpoint readback above supplies the same physical rules.
    for label, report in reports.items():
        samples = report['physical_samples'] or native_samples['coarse-h0']
        result['runs'][label] = summarize(report, records[label], samples)
        result['runs'][label]['timing'] = timing(Path(records[label]['report']).parent)
    h0 = native_samples['coarse-h0']
    result['differences'] = dict(
        coarse_minus_fine=difference(native_samples['coarse8'], native_samples['fine16'], scales),
        coarse_drift=difference(native_samples['coarse8'], h0, scales),
        fine_drift=difference(native_samples['fine16'], h0, scales))
    result['passed'] = True
    write_json(a.output/'analysis.json', result)


if __name__ == '__main__':
    main()
