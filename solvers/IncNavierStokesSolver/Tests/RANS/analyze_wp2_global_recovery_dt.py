#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact common-time recovery comparison from the original physical seed.

Only immutable progress snapshots are read. Neither accepted states nor
histories are advanced, projected, relabelled, or interpolated in time.
"""
import argparse, re, resource, time
from pathlib import Path
import numpy as np
from analyze_wp2_approach import native_states, paths
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, state_hash
from wp2_all_state import Rows, US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate as analytical_reference
from wp2_qualification import common_differences, quarter_budget_differences


def retained(root, progress_path):
    progress = load(progress_path)
    require(progress['model'] == 'standard', 'global control model differs')
    require(progress['mode'] in ('recover', 'recovery-trial'), 'not the original recovery trajectory')
    bindings = {str(progress_path): sha(progress_path)}
    reference = defect = context = None
    for path in paths(root, 'standard', 'recovery', progress):
        report = load(path)
        bindings[str(path)] = sha(path)
        if context is None:
            context = report
            reference, defect = (state_hash(report[k]) for k in ('fixed_reference', 'fixed_defect'))
        require(state_hash(report['fixed_reference']) == reference and
                state_hash(report['fixed_defect']) == defect, 'immutable data changed between reports')
    for entry in progress['results'].values():
        if entry.get('report') in bindings:
            require(bindings[entry['report']] == entry['report_sha256'], 'retained report changed')
    states, _ = native_states(root, 'standard', 'recovery', progress_path)
    require(states and min(states) == 0., 'global comparison must include the original pulse at t=0')
    return states, context, bindings, progress


def comparison_signature(report):
    # Ordinary restart still compares the entire signature. Only this
    # adjacent-dt readback excludes its explicitly validated timestep entry.
    pattern = r'\|TIMESTEP=([^|]+)'
    values = re.findall(pattern, report['signature'])
    require(len(values) == 1 and float(values[0]) == report['dt'],
            'signature timestep does not match the actual native report')
    return re.sub(pattern, '|TIMESTEP=<adjacent-dt comparison>', report['signature'])


def main():
    parser = argparse.ArgumentParser()
    for key in ('root', 'coarse', 'fine', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    coarse, context, bindings, pc = retained(args.root, args.coarse)
    fine, fine_context, fb, pf = retained(args.root, args.fine)
    bindings.update(fb)
    require(pc['dt'] == 2 * pf['dt'], 'not an adjacent timestep pair')
    for key in ('fixed_reference', 'fixed_defect'):
        require(state_hash(context[key]) == state_hash(fine_context[key]), 'pair changed ' + key)
    require(comparison_signature(context) == comparison_signature(fine_context),
            'physical signature differs beyond the declared timestep')
    seed = fine[0.][1]
    require(state_hash(seed) == state_hash(coarse[0.][1]), 'pair changed the original physical pulse')
    rows = Rows(args.root / 'geometry-rows-selected/result.json')
    K0, H = context['parameters'][5], context['parameters'][1]
    reference_values = None

    def evaluate(state):
        values = rows.tuple(state, context)
        dense = resample(dict(context,
            wp2_common_samples=dict(context['wp2_common_samples'], values=values['values']),
            physical_samples=dict(wall_K=values['endpoint_K'])), 28)
        sample = dict(columns=context['wp2_common_samples']['columns'], values=dense['values'],
                      weights=dense['weights'].tolist(), regions=dense['regions'].tolist())
        if reference_values is not None:
            sample['reference_values'] = reference_values
        return sample, dense

    seed_sample, seed_dense = evaluate(seed)
    reference_values = analytical_reference(context, seed_dense['physical'])['values']
    seed_sample['reference_values'] = reference_values
    history = []
    for key in sorted(set(coarse) & set(fine)):
        tc, sc = coarse[key]
        tf, sf = fine[key]
        require(abs(tc - tf) < 2e-13, 'actual common times differ')
        left, dl = evaluate(sc)
        right, dr = evaluate(sf)
        difference = common_differences(left, right, US, K0, H)
        cr = common_differences(left, seed_sample, US, K0, H)
        fr = common_differences(right, seed_sample, US, K0, H)
        response = {}
        for name in ('U', 'W', 'k', 'epsilon', 'pi', 'mechanical_pressure'):
            den = fr[name]['dimensional_L2']
            response[name] = None if den <= 1e-12 else difference[name]['dimensional_L2'] / den
        tensors, minima = {}, []
        for side in (0, 1):
            mask = dl['regions'] == side
            weight = dl['weights'][mask]
            for name in ('R', 'T'):
                delta = dl[name][mask] - dr[name][mask]
                tensors[str(side) + '_' + name] = float(np.sqrt(
                    np.sum(weight[:, None, None] * delta ** 2) / sum(weight)) / US ** 2)
        for dense in (dl, dr):
            R = dense['R']
            eig = np.minimum(R[:, 1, 1], .5 * (R[:, 0, 0] + R[:, 2, 2] -
                np.hypot(R[:, 0, 0] - R[:, 2, 2], 2 * R[:, 0, 2])))
            minima.append([float(min(eig[dense['regions'] == side])) for side in (0, 1)])
        history.append(dict(coarse_actual_time=tc, fine_actual_time=tf,
            coarse_state_sha256=state_hash(sc), fine_state_sha256=state_hash(sf),
            differences=difference, coarse_response_from_original_seed=cr,
            fine_response_from_original_seed=fr, difference_fraction_of_fine_response=response,
            quarter_budget=quarter_budget_differences(left, right, US, K0, H),
            full_tensor_RMS_difference_over_ustar2=tensors,
            full_tensor_quarter_budget_passed=max(tensors.values()) <= .0025,
            GL28_covariance_inner_outer_coarse_fine=minima))
    require(history, 'no exact common-time comparison')
    for path, expected in bindings.items():
        require(sha(Path(path)) == expected, 'analysis input changed')
    atomic(args.output, dict(version='WP2GlobalRecoveryAdjacentTimestepV1',
        coarse_status=pc['status'], fine_status=pf['status'], history=history,
        common_seed_sha256=state_hash(seed), fixed_reference_sha256=state_hash(context['fixed_reference']),
        fixed_defect_sha256=state_hash(context['fixed_defect']), bindings=bindings,
        response_denominator_floor_dimensional_L2=1e-12, pressure='absolute; no gauge realignment',
        quadrature='seen GL28 physical union partition; coefficient-derived fields and current endpoint K',
        source_sha256=sha(Path(__file__)), seconds=time.monotonic() - started,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Global adjacent-timestep differences from the original pulse; no temporal-order or future-recovery claim. '
              'Dependent inner epsilon excluded from transported epsilon norms. Zero governing solves.'))


if __name__ == '__main__':
    main()
