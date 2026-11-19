#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replay a fixed-state pre-elimination A/B outlet test on Linux."""
import argparse
import json
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
import numpy as np


LIMIT = 1e-10
OUTLET_ORDER = 20
VOLUME_ORDERS = (6, 8)


def outlet_face(state, element):
    nodes, weights = po.leggauss(OUTLET_ORDER)
    _, right, bottom, top = element['bounds']
    panels = state.cuts(bottom, top)
    r = np.concatenate([(lo+hi)/2+(hi-lo)*nodes/2 for lo, hi in zip(panels[:-1], panels[1:])])
    q = np.concatenate([weights*(hi-lo)/2 for lo, hi in zip(panels[:-1], panels[1:])])
    return np.full_like(r, right), r, q


def static_volume(state, element, orders=VOLUME_ORDERS):
    xn, xw = po.leggauss(orders[0])
    rn, rw = po.leggauss(orders[1])
    left, right, bottom, top = element['bounds']
    panels = state.cuts(bottom, top)
    x = (left+right)/2+(right-left)*xn/2
    xs, rs, weights = [], [], []
    for lo, hi in zip(panels[:-1], panels[1:]):
        r = (lo+hi)/2+(hi-lo)*rn/2
        xx, rr = np.meshgrid(x, r)
        xs.extend(xx.ravel())
        rs.extend(rr.ravel())
        weights.extend(np.outer(rw, xw).ravel()*(hi-lo)*(right-left)/4)
    return np.asarray(xs), np.asarray(rs), np.asarray(weights)


def independent_actions(state, orders=VOLUME_ORDERS):
    """Rebuild the declared discrete actions with independent bases and kernels."""
    normal, tangent, discarded = np.zeros((3, state.n))
    for e, side in state.faces['u']:
        if side != 1:
            continue
        x, r, weight = outlet_face(state, e)
        b = po.basis(e, 'u', x, r)[0]
        prescribed, full = state.outlet_traction(r), state.ref.outlet(r)
        normal[e['u']] += b.T@(weight*prescribed[:, 0])
        tangent[e['u']] += b.T@(weight*prescribed[:, 1])
        discarded[e['u']] += b.T@(weight*full[:, 1])
    # Only neighboring elements with support on the constrained transverse
    # trace contribute. The flat sampled AWF force has no transverse term.
    momentum = np.zeros(2*state.n+state.np)
    touched = set(state.outlet_tangent_fixed)
    for e in state.elements:
        if not touched.intersection(e['u']):
            continue
        x, r, weight = static_volume(state, e, orders)
        a = state.sample(e, x, r)
        b, dx, dr = po.basis(e, 'u', x, r)
        advection = a['U']*a['G'][:, 2, 0]+a['W']*a['G'][:, 2, 2]
        momentum[state.n+np.asarray(e['u'])] += (
            dx.T@(weight*a['V'][:, 2, 0])+dr.T@(weight*(a['V'][:, 2, 2]-a['Pi']))+
            b.T@(weight*advection))
    momentum[state.n:2*state.n] -= tangent
    expected = dict(normal_applied=normal, tangent_applied=tangent,
                    tangent_reaction=state.outlet_tangent_reactions(momentum))
    return expected, discarded


def compare_actions(actual, expected, state):
    differences, vectors = {}, {}
    for key in ('normal_applied', 'tangent_applied'):
        want, got = expected[key], np.asarray(actual[key], dtype=float)
        po.require(got.shape == want.shape and np.all(np.isfinite(got)), 'invalid native '+key+' vector')
        vectors[key] = got
        differences[key] = float(np.max(abs(got-want))/state.ref.tau)
    rows = actual['tangent_reaction']
    po.require(all(type(row['dof']) is int for row in rows), 'noninteger native reaction DOF')
    got = {row['dof']: row for row in rows}
    want = {row['dof']: row for row in expected['tangent_reaction']}
    po.require(len(got) == len(rows) and set(got) == set(want), 'wrong native constrained reaction inventory')
    po.require(all(got[dof]['scope'] == row['scope'] and math.isfinite(got[dof]['value'])
                   for dof, row in want.items()), 'wrong native reaction attribution or nonfinite action')
    differences['tangent_reaction'] = max(abs(got[dof]['value']-row['value'])/state.ref.tau
                                           for dof, row in want.items())
    return differences, vectors


def pressure_level(state):
    shift, face = np.zeros((2, 2, state.n))
    for element in state.elements:
        state.clock()
        x, r, weight = static_volume(state, element)
        _, dx, dr = po.basis(element, 'u', x, r)
        ids = np.asarray(element['u'])
        shift[0, ids] -= dx.T@weight
        shift[1, ids] -= dr.T@weight
    for element, side in state.faces['u']:
        x, r, weight = state.face(element, side)
        b = po.basis(element, 'u', x, r)[0]
        for component, normal in enumerate(po.NORMALS[side]):
            face[component, element['u']] -= normal*(b.T@weight)
    free = [float(np.max(abs(shift[c, ids]))) for c, ids in enumerate(state.ufree)]
    defects = [float(np.max(abs(shift[c, ids]-face[c, ids]))) for c, ids in enumerate(state.ufree)]
    po.require(state.data['gauge'] is False and free[0] > 1e-8 and max(defects) <= LIMIT,
               'open-normal pressure level is not supported by the free tests')
    return dict(pressure_gauge=False, unit_Pi_free_action_MAX=free,
                volume_minus_face_free_action_MAX=defects)


def corners(state, vectors, expected):
    result = {}
    for adjacent, dof in state.outlet_corners.items():
        r = 0. if adjacent == 'ground' else state.ref.top
        owners = [e for e, side in state.faces['u'] if side == 1 and dof in e['u']]
        po.require(len(owners) == 1, 'missing unique outlet owner of shared vertex')
        owner = owners[0]
        x, r = np.array([owner['bounds'][1]]), np.array([r])
        b = po.basis(owner, 'u', x, r)[0][0]
        a = state.sample(owner, x, r)
        ids = np.asarray(owner['u'])
        free_tests = [float(max(abs(b[np.isin(ids, free)]), default=0.)) for free in state.ufree]
        po.require(abs(free_tests[0]-(adjacent == 'top')) <= LIMIT and free_tests[1] <= LIMIT,
                   'incorrect axial or transverse freedom at shared outlet corner')
        result[adjacent] = dict(dof=int(dof), owner=owner['geometry'], x=float(x[0]), r=float(r[0]),
            free_test_trace_MAX=free_tests,
            prescribed_normal_action=float(vectors['normal_applied'][dof]),
            prescribed_tangential_action=float(vectors['tangent_applied'][dof]),
            signed_tangent_reaction_action=next(row['value'] for row in expected['tangent_reaction']
                                                if row['dof'] == dof),
            raw_stress_traction=[float(a['V'][0, 0, 0]-a['Pi'][0]), float(a['V'][0, 2, 0])],
            Ux=float(a['G'][0, 0, 0]), Wr=float(a['G'][0, 2, 2]), D=float(a['div'][0]))
    if state.outlet_policy == po.ZERO_TRANSVERSE:
        for key in ('Ux', 'Wr', 'D'):
            po.require(abs(result['ground'][key]) <= LIMIT, 'B ground vertex violates its coupled traces')
    return result


def action_test(state):
    native = state.data['outlet_actions']
    expected, full_tangent = independent_actions(state)
    differences, vectors = compare_actions(native, expected, state)
    po.require(max(differences.values()) <= LIMIT, 'independent pre-elimination actions differ: '+str(differences))
    free_u = np.intersect1d(state.outlet_ids, state.ufree[0])
    free_w = np.intersect1d(state.outlet_ids, state.ufree[1])
    free_axial_load = float(np.max(abs(vectors['normal_applied'][free_u]))/state.ref.tau)
    po.require(free_axial_load > LIMIT, 'prescribed normal load has no exercised free outlet tests')
    result = dict(state_sha256=state.data['state_sha256'], policy=state.outlet_policy,
        policy_identity=state.outlet_identity, pressure_level=pressure_level(state),
        fixture_elements=len(state.elements),
        independent_source=po.identity(po.__file__), scaled_absolute_errors=differences,
        free_axial_normal_action_scaled_MAX=free_axial_load,
        outlet_trace_dofs=state.outlet_ids.tolist(), free_axial_outlet_dofs=free_u.tolist(),
        free_transverse_outlet_dofs=free_w.tolist(),
        corners=corners(state, vectors, expected),
        independent_actions={key: value.tolist() if isinstance(value, np.ndarray) else value
                             for key, value in expected.items()}, native_actions=native)
    if state.outlet_policy == po.FULL_TRACTION:
        free_tangent_load = float(np.max(abs(vectors['tangent_applied'][free_w]))/state.ref.tau)
        po.require(free_tangent_load > LIMIT, 'A transverse natural load has no exercised free outlet tests')
        result['free_transverse_tangent_action_scaled_MAX'] = free_tangent_load
    else:
        w = state.s['flow'][state.n:2*state.n]
        po.require(float(np.max(abs(w[state.outlet_ids]))) <= LIMIT and
                   float(np.max(abs(vectors['tangent_applied']))) == 0.,
                   'B transverse essential coefficients or omitted natural action are nonzero')
        mutant = dict(native, tangent_applied=vectors['tangent_applied']+full_tangent)
        mutated, _ = compare_actions(mutant, expected, state)
        free = float(np.max(abs(full_tangent[state.ufree[1]]))/state.ref.tau)
        work = float(abs(w@full_tangent)/state.ref.tau)
        corner_errors = {key: float(abs(full_tangent[dof])/state.ref.tau)
                         for key, dof in state.outlet_corners.items()}
        po.require(mutated['tangent_applied'] > LIMIT and min(corner_errors.values()) > LIMIT and
                   free <= LIMIT and work <= LIMIT, 'reintroduced traction negative control lacks coverage')
        result['reintroduced_full_tangent_negative_control'] = dict(status='DETECTED',
            pre_elimination_scaled_MAX=mutated['tangent_applied'], corner_scaled_errors=corner_errors,
            restricted_free_scaled_MAX=free, work_scaled_abs=work)
    return result


def load_state(path, reference, policy, deadline):
    data, metadata = po.load_native(path)
    po.require(metadata.get('normal_main') is True and metadata.get('complete_step_seconds') == [] and
               data['step'] == 0 and data['time'] == 0., 'outlet probe requires a registered zero-time state')
    po.require(data['outlet_policy'] == policy, 'wrong explicit outlet policy')
    po.require((data['quadratureX'], data['quadratureR']) == VOLUME_ORDERS,
               'different fixed-state fixture quadrature')
    return po.PhysicalOracle(data, reference, OUTLET_ORDER, deadline)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('a', 'b', 'table', 'out'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    parser.add_argument('--source', help='Expected compiled source identity')
    args = parser.parse_args()
    po.require(sys.platform.startswith('linux'), 'outlet action probe execution is Linux-only')
    po.require(math.isfinite(args.seconds) and args.seconds > 0, 'explicit finite wall cap required')
    po.require(not args.out.exists(), 'immutable outlet proof already exists')
    start = time.monotonic()
    report = dict(schema='WJH27OutletPreEliminationProbeV1', status='FAIL',
                  inputs={key: po.identity(getattr(args, key)) for key in ('a', 'b', 'table')},
                  source=po.identity(__file__), basis_source=po.identity(po.fixtures.__file__),
                  volume_gauss_orders=VOLUME_ORDERS, traction_panel_gauss_order=OUTLET_ORDER,
                  action_error_scale='tau0 times the fixed unit code length', tolerance=LIMIT)
    try:
        reference = po.load_reference(args.table)
        a = load_state(args.a, reference, po.FULL_TRACTION, start+args.seconds)
        b = load_state(args.b, reference, po.ZERO_TRANSVERSE, start+args.seconds)
        compiled = a.data['source_identity']
        po.require(len(compiled) == 64 and all(c in '0123456789abcdef' for c in compiled) and
                   compiled == b.data['source_identity'] and (args.source is None or compiled == args.source),
                   'probe cases require the same expected native source')
        po.require(a.geometry == b.geometry and a.elements == b.elements and a.outlet_identity != b.outlet_identity and
                   a.data['signature'] != b.data['signature'] and np.array_equal(a.ufree[0], b.ufree[0]) and
                   len(a.ufree[1]) > len(b.ufree[1]), 'fixture geometry, outlet identity or masks disagree')
        report['native_source_identity'] = compiled
        report['cases'] = {'A': action_test(a), 'B': action_test(b)}
        aa, bb = (np.asarray(report['cases'][key]['native_actions']['normal_applied']) for key in ('A', 'B'))
        report['between_policy_normal_action_scaled_MAX'] = float(np.max(abs(aa-bb))/reference.tau)
        po.require(report['between_policy_normal_action_scaled_MAX'] <= LIMIT, 'normal traction changed')
        report['status'] = 'PASS'
    except (KeyError, TypeError, ValueError, TimeoutError) as error:
        report['reason'] = str(error)
    report['elapsed_seconds'] = time.monotonic()-start
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'elapsed_seconds', 'reason') if key in report}))
    if report['status'] != 'PASS':
        raise SystemExit(2)


if __name__ == '__main__':
    main()
