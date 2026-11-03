#!/usr/bin/env python3
"""Independent bounded all-mode audit of the selected curved M8 operator.

The independently reconstructed original M5/M6 weak operator is augmented by
the exact selected outer stress/production increments.  An explicit all-ray
plateau check is required before reusing the original inner source.  These
are unsolved generic states, not manufactured curved gamma solutions or
qualified terrain trajectories.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

import validate_cls_epawf_m5 as weak
import cls_epawf_m6_operator_oracle as transition
from validate_cls_epawf_m1 import indexed
from validate_cls_epawf_m2 import Profiles
from cls_epawf_m7_energy_oracle import Audit, _horizontal
from cls_epawf_m8_native_oracle import normal_rule, selected_stress
from cls_epawf_m8_protocol import sha
from validate_cls_epawf_m5_geometry import geometry


def plateau(raw, state, fields):
    """Every declared ray/source location must be in the open inactive set."""
    if not raw['robust_parameters'][6]:
        return dict(required=False, samples=0, minimum_k_over_K=None)
    d, offsets = raw['dimension'], raw['offsets']
    nu, H, rough, rK, rm, _, qstar, kappa = raw['parameters']
    minimum, count = float('inf'), 0
    for chart, _ in _horizontal(fields, d, H, int(raw['quadrature_orders'][0])):
        metric = geometry(raw['map_parameters'], chart)
        sampled = fields[d].evaluate(state, offsets[d], [chart[0], chart[1], rK/metric['scale']])[0]
        chi = transition.weight(raw, state, fields, chart) if 'transition' in raw else raw['chi']
        profile = Profiles(dict(geometry=[nu, metric['height'], rough, rK, rm],
            metric=[metric['kappa1'], metric['kappa2']], k=sampled, q=qstar, kappa=kappa, chi=chi))
        points, _ = normal_rule(profile, raw, fields, metric)
        for distance in points:
            point = [chart[0], chart[1], distance/metric['scale']]
            ratio = fields[d].evaluate(state, offsets[d], point)[0]/sampled
            minimum, count = min(minimum, ratio), count+1
    if not minimum > raw['robust_parameters'][2]:
        raise RuntimeError('Inner repaired source left its open inactive plateau; the original-source weak oracle cannot be reused')
    return dict(required=True, samples=count, minimum_k_over_K=minimum,
                theta=raw['robust_parameters'][2], source_repair_identically_inactive=True)


def residual(raw, state, fields, normal_order=24):
    source_plateau = plateau(raw, state, fields)
    view = dict(raw, lifted_state=np.asarray(state))
    options = dict(normal_order=normal_order, verify=False, return_loads=True,
                   tangential_order=int(raw['quadrature_orders'][0]),
                   outer_order=int(raw['quadrature_orders'][0]))
    if 'transition' in raw:
        options['weight_function'] = lambda chart, values, spaces: transition.weight(raw, values, spaces, chart)
    original = weak.residual_oracle(view, **options)
    result = np.asarray(original['residual'])
    change = np.zeros_like(result)
    d, offsets = raw['dimension'], raw['offsets']
    axes = (0, 2) if d == 2 else (0, 1, 2)
    stress_work = production_work = 0.
    for element in fields[d+3].elements:
        points, weights = fields[d+3].quadrature(element, int(raw['quadrature_orders'][0]))
        velocity = [weak.physical_rows(fields[0], state, offsets[block], points, element['geometry'])
                    for block in range(d)]
        ids, phi, derivative = velocity[0][:3]
        G = np.zeros((len(points), 3, 3))
        for block, axis in enumerate(axes):
            G[:, axis] = velocity[block][4]
        ki, kp, _, K, _ = weak.physical_rows(fields[d], state, offsets[d], points, element['geometry'])
        ei, ep, _, epsilon, _ = weak.physical_rows(fields[d+1], state, offsets[d+1], points, element['geometry'])
        delta_R = np.asarray([selected_stress(raw, g, k, e, point)[0]-indexed(g, k, e, 'physical')['stress'].real
                              for g, k, e, point in zip(G, K, epsilon, points)])
        delta_P = -np.einsum('qab,qab->q', delta_R, G)
        for block, axis in enumerate(axes):
            np.add.at(change, offsets[block]+ids,
                      np.einsum('qai,qa,q->i', derivative, -delta_R[:, axis], weights))
        np.add.at(change, offsets[d]+ki, kp.T@(-weights*delta_P))
        np.add.at(change, offsets[d+1]+ei, ep.T@(-weights*1.44*epsilon/K*delta_P))
        stress_work += float(weights@(-np.einsum('qab,qab->q', delta_R, G)))
        production_work += float(weights@delta_P)
    result += change
    return result, dict(source_plateau=source_plateau,
                        maximum_added_component=float(np.max(np.abs(change))),
                        modified_stress_weak_work=stress_work,
                        production_increment=production_work,
                        mean_plus_k_change=stress_work-production_work,
                        original_energy=original['energy_budget'])


def coupled_residual(raw, state, fields, trace_ids, normal_order=24):
    physical = transition._external_lift(raw, state, trace_ids)
    loads, detail = residual(raw, physical, fields, normal_order)
    offset = raw['offsets'][raw['dimension']+1]
    for row, value in raw['constraints']:
        row = int(row)
        if not (offset <= row < offset+fields[raw['dimension']+1].size and row-offset in trace_ids):
            loads[row] = state[row]-value
    baseline, correction = transition.trace_projection(raw, physical, fields)
    for index in trace_ids:
        loads[offset+index] = state[offset+index]-baseline[index]-correction[index]
    return loads, detail


def directions(raw, fields, audit):
    if 'transition' not in raw:
        return dict(independent_full_weak_actions=False, reason='Fixed-mode probe does not export coupled directional samples')
    if raw['dimension'] != 2:
        return dict(independent_full_weak_actions=False, reason='Bounded independent directional oracle is 2D; 3D complete weak values remain audited')
    state = np.asarray(raw['lifted_state'])
    maps, _ = transition._trace_projection(raw, state, fields)
    trace_ids = set(maps[0])
    result = []
    # Include velocity directions, where persistent gamma changes the
    # constitutive tangent, as well as the nonlocal turbulence/pressure rows.
    for case in raw['transition']['directions']:
        direction = np.asarray(case['vector'])
        audit.require(np.max(np.abs(direction)) > 0, 'Empty selected-operator direction')
        native = {float(action['step']): np.asarray(action['value']) for action in case['actions']}
        levels = []
        for step in (2e-4, 1e-4):
            plus, dp = coupled_residual(raw, state+step*direction, fields, trace_ids)
            minus, dm = coupled_residual(raw, state-step*direction, fields, trace_ids)
            action = (plus-minus)/(2*step)
            selected = min(native, key=lambda value: abs(value-step))
            audit.require(abs(selected-step) < 1e-14, 'Missing native selected-operator action level')
            audit.near(action, native[selected], 'm8_curved_complete_action_block_'+str(case['block'])+'_h_'+str(step), 1e-7)
            levels.append(dict(step=step, maximum_native_error=float(np.max(np.abs(action-native[selected]))),
                plus_minimum_k_over_K=dp['source_plateau']['minimum_k_over_K'],
                minus_minimum_k_over_K=dm['source_plateau']['minimum_k_over_K']))
        result.append(dict(block=case['block'], levels=levels,
                           includes_velocity_stress_increment=case['block'] < raw['dimension']))
    return dict(independent_full_weak_actions=True, blocks=result,
                scope='Central differences of independently reconstructed selected all-field weak loads and nonlinear matching rows')


def verify(raw, directional=True):
    audit = Audit()
    audit.require(raw.get('eta') == 0 and raw.get('model_variant') in ('Original', 'CubicScaled', 'Regularized'),
                  'Missing unanchored selected-model probe identity')
    audit.require(raw.get('robust_parameters') is not None and raw.get('ustar', 0) > 0,
                  'Missing robust parameter/reference metadata')
    fields, state = weak.spaces(raw), np.asarray(raw['lifted_state'])
    expected, detail = residual(raw, state, fields)
    audit.near(raw['raw_residual'], expected, 'm8_curved_selected_all_mode_weak_operator', 1e-10)
    audit.near(raw['energy_budget'], detail['original_energy'], 'm8_curved_selected_total_energy_exchange', 2e-10)
    audit.near(detail['mean_plus_k_change'], 0., 'm8_curved_same_stress_production_exchange', 1e-12)
    blocks = []
    actual = np.asarray(raw['raw_residual'])
    for block in range(raw['dimension']+4):
        start, end = raw['offsets'][block:block+2]
        difference = actual[start:end]-expected[start:end]
        blocks.append(dict(block=block, maximum_absolute_error=float(np.max(np.abs(difference), initial=0.))))
    action = directions(raw, fields, audit) if directional else dict(independent_full_weak_actions=False, reason='Disabled by caller')
    return dict(schema=1, success=True, checks=audit.checks, errors=audit.errors,
        model_variant=raw['model_variant'], robust_parameters=raw['robust_parameters'],
        dimension=raw['dimension'], map_parameters=raw['map_parameters'],
        complete_weak_blocks=blocks, operator_increment=detail, directional=action,
        state_scope='Unsolved generic bounded coefficient state with independently verified inactive source-repair plateau',
        variable_matching_height=any(value != 0 for value in raw['map_parameters'][3:6]),
        genuine_curved_gamma_mms=False,
        excluded=['exact curved gamma manufactured solution', 'admissible state-envelope coverage',
                  'positivity/covariance qualification', 'long-time or physical terrain validation'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--raw', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--skip-directional', action='store_true')
    a = p.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific numerical operator reconstruction must run on Linux')
    receipt = dict(started_utc=datetime.now(timezone.utc).isoformat(),
        input=dict(path=str(a.raw), sha256=sha(a.raw)), driver_sha256=sha(__file__),
        numerical_threads=THREADS)
    try:
        receipt.update(verify(json.loads(a.raw.read_text()), not a.skip_directional))
    except Exception as error:
        receipt.update(success=False, status='failed', error=str(error))
    if sha(a.raw) != receipt['input']['sha256']:
        receipt.update(success=False, status='input_changed')
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    a.output.parent.mkdir(parents=True, exist_ok=True)
    if a.output.exists():
        raise RuntimeError('Refusing to overwrite an existing selected-operator receipt')
    a.output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: receipt[key] for key in ('success', 'checks', 'errors') if key in receipt}))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
