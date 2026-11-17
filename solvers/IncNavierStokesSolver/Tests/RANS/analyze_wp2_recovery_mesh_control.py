#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Same-dt original-pulse comparison on the qualified level4/level5 meshes."""
import argparse, resource, time
from pathlib import Path
import numpy as np
from analyze_wp2_global_recovery_dt import retained
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require, sha, state_hash
from wp2_all_state import Rows, US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate as reference
from wp2_qualification import quarter_budget_differences


def same_physical_problem(left, right):
    for key in ('parameters', 'model_specification', 'ustar'):
        require(left[key] == right[key], 'mesh control changes physical ' + key)
    require(left['developing_boundary']['definition'] == right['developing_boundary']['definition'],
            'mesh control changes physical boundary interpretation')


def compare_dense(left, right, reference_values, columns, K0):
    for key in ('physical', 'weights', 'regions'):
        require(np.array_equal(left[key], right[key]), 'mesh control physical integration grid differs: ' + key)
    base = dict(columns=columns, weights=left['weights'].tolist(),
        regions=left['regions'].tolist(), reference_values=reference_values)
    principal = quarter_budget_differences(dict(base, values=left['values']),
        dict(base, values=right['values']), US, K0, .0065)
    tensors = {}
    for region in (0, 1):
        mask = left['regions'] == region
        weight = left['weights'][mask]
        for name in ('R', 'T'):
            delta = left[name][mask] - right[name][mask]
            tensors[str(region) + '_' + name] = float(np.sqrt(
                np.sum(weight[:, None, None] * delta ** 2) / sum(weight)) / US ** 2)
    minima = []
    for dense in (left, right):
        R = dense['R']
        eig = np.minimum(R[:, 1, 1], .5 * (R[:, 0, 0] + R[:, 2, 2] -
            np.hypot(R[:, 0, 0] - R[:, 2, 2], 2 * R[:, 0, 2])))
        minima.append([float(min(eig[dense['regions'] == region])) for region in (0, 1)])
    return dict(principal=principal, full_tensor_RMS_over_ustar2=tensors,
        GL28_covariance_inner_outer_selected_confirmation=minima,
        passed=principal['passed'] and max(tensors.values()) <= .0025)


def main():
    parser = argparse.ArgumentParser()
    for key in ('root', 'selected', 'confirmation', 'output'):
        parser.add_argument('--' + key, type=Path, required=True)
    args = parser.parse_args()
    started = time.monotonic()
    coarse, lc, bindings, pc = retained(args.root, args.selected)
    fine, rc, other, pf = retained(args.root, args.confirmation)
    bindings.update(other)
    require(pc['dt'] == pf['dt'], 'mesh comparison changed timestep')
    require(pc['case_definition'] == dict(level=4, order=16, projection='RadialEnergyV1') and
            pf['case_definition'] == dict(level=5, order=16, projection='RadialEnergyV1'),
            'mesh comparison is not the qualified level4/level5 pair')
    require(pc['solver_sha256'] == pf['solver_sha256'], 'mesh comparison changed runtime')
    same_physical_problem(lc, rc)
    rows = [Rows(args.root / ('geometry-rows-' + name + '/result.json'))
            for name in ('selected', 'confirmation')]
    replay = [row.check_native(context) for row, context in zip(rows, (lc, rc))]
    history = []
    reference_values = None
    for key in sorted(set(coarse) & set(fine)):
        tl, sl = coarse[key]
        tr, sr = fine[key]
        require(abs(tl - tr) < 2e-13, 'mesh control has no exact common physical time')
        actual = []
        for row, context, state in zip(rows, (lc, rc), (sl, sr)):
            values = row.tuple(state, context)
            actual.append(resample(dict(context,
                wp2_common_samples=dict(context['wp2_common_samples'], values=values['values']),
                physical_samples=dict(wall_K=values['endpoint_K'])), 28))
        if reference_values is None:
            reference_values = reference(lc, actual[0]['physical'])['values']
            require(np.array_equal(reference_values, reference(rc, actual[0]['physical'])['values']),
                    'the two meshes do not share the same analytical reference')
        comparison = compare_dense(actual[0], actual[1], reference_values,
            lc['wp2_common_samples']['columns'], lc['parameters'][5])
        history.append(dict(selected_actual_time=tl, confirmation_actual_time=tr,
            selected_state_sha256=state_hash(sl), confirmation_state_sha256=state_hash(sr), **comparison))
    require(history and history[-1]['confirmation_actual_time'] >= pf['Treq'], 'mesh control misses required horizon')
    for path, expected in bindings.items():
        require(sha(Path(path)) == expected, 'mesh comparison input changed')
    atomic(args.output, dict(version='WP2RecoveryMatchedMeshControlV1', history=history,
        common_times=len(history), passed=all(row['passed'] for row in history),
        fixed_reference_sha256=[state_hash(c['fixed_reference']) for c in (lc, rc)],
        initial_seed_sha256=[state_hash(coarse[0.][1]), state_hash(fine[0.][1])],
        native_tuple_replay=replay, bindings=bindings,
        row_catalogues=[dict(path=str(row.path), sha256=sha(row.path)) for row in rows],
        source_sha256=sha(Path(__file__)), seconds=time.monotonic()-started,
        peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='Same model, physical data, continuum pulse and dt; distinct native projections/references/histories. '
              'Seen GL28 common physical panels, freshly recomputed selected tensors, absolute pressure, '
              'no field transfer or time interpolation. Two meshes do not establish spatial order. Zero governing solves.'))


if __name__ == '__main__':
    main()
