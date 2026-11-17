#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only packet, ownership and independent basis controls for InitialFields.

Retains every valid/invalid fixture, native stdout/stderr and output in an
immutable sibling fixture directory. Does not launch a solver or alter history.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
from types import SimpleNamespace

import physical_oracle_test
import physical_oracle as po
import numpy as np


def state_hash(state):
    h = hashlib.sha256()
    for key in ('flow', 'k', 'epsilon'):
        h.update(np.asarray(state['coefficients'][key], dtype='<f8').tobytes())
    return h.hexdigest()


def check_samples(packet, result, path):
    state = packet['state']; g = state['geometry']; n = state['sizes'][0]
    po.require(result['schema'] == 'WJH27InitialFieldsProbeV1' and result['status'] == 'RECORDED', 'bad native probe schema/status')
    po.require(result['packet_sha256'] == po.identity(path)['sha256'] and result['state_sha256'] == state['state_sha256'], 'packet/coefficient identity mismatch')
    qualified = packet['schema'] == 'WJH27RegisteredCoupledRunV2'
    po.require(result['qualified_reference'] == qualified and result['role'] == ('QualifiedReference' if qualified else 'Candidate'), 'candidate was mislabeled or reference role lost')
    fields = (("u", "flow", 0), ("u", "flow", n), ("k", "k", 0), ("epsilon", "epsilon", 0), ("p", "flow", 2*n))
    maxima = [max(abs(v) for v in state['coefficients']['flow'][:n]), max(abs(v) for v in state['coefficients']['flow'][n:2*n]),
              max(abs(v) for v in state['coefficients']['k']), max(abs(v) for v in state['coefficients']['epsilon']),
              max(abs(v) for v in state['coefficients']['flow'][2*n:])]
    elements = {e['geometry']: e for e in state['elements']}
    maximum_error = 0.
    for row in result['samples']:
        x, r = row['x'], row['physical_z']-g['d0']
        def same(a, b): return abs(a-b) <= 64*np.finfo(float).eps*max(abs(a), abs(b), g['H'])
        for edge in (g['left'], g['right']):
            if same(x, edge): x = edge
        for edge in (0., g['H'], g['top']):
            if same(r, edge): r = edge
        def contains(e):
            xl, xr, rl, rr = e['bounds']
            return (x >= xl or same(x, xl)) and (x <= xr or same(x, xr)) and (r >= rl or same(r, rl)) and (r <= rr or same(r, rr))
        for f, (field, vector, offset) in enumerate(fields):
            for mode in ('owned', 'automatic'):
                value = row[mode][f]
                if mode == 'owned':
                    e = elements[row['owner']]
                    absent = not e[field]
                else:
                    absent = f == 3 and r < g['H']
                    owners = [e for e in elements.values() if e[field] and contains(e)]
                    e = min(owners, key=lambda a: a['geometry']) if owners else None
                if absent:
                    po.require(value is None, 'native probe fabricated inner epsilon'); continue
                po.require(value is not None and e is not None and contains(e), 'missing valid owned value')
                coefficients = np.asarray(state['coefficients'][vector])[offset+np.asarray(e[field])]
                expected = float(po.basis(e, field, np.array([x]), np.array([r]))[0]@coefficients)
                error = abs(value-expected)/maxima[f]
                maximum_error = max(maximum_error, error)
                po.require(error <= 2e-12, f'independent Modified_A field {f} value error {error}')
    names = {v['name'] for v in result['rejections']}
    po.require(names == {'inner_epsilon', 'outside_domain', 'invalid_field', 'absent_owner', 'nonfinite_coordinate',
        'explicit_inner_owner_at_H', 'wrong_explicit_owner', 'source_binding', 'reference_binding', 'model_binding',
        'manufactured_binding', 'variant_binding', 'geometry_binding'}, 'runtime rejection inventory differs')
    return dict(maximum_scaled_value_error=maximum_error, samples=len(result['samples']), role=result['role'], runtime_rejections=len(names))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    args = parser.parse_args()
    po.require(sys.platform.startswith('linux'), 'InitialFields numerical controls are Linux-only')
    po.require(np.isfinite(args.seconds) and args.seconds > 0 and not args.out.exists(), 'finite cap and immutable output required')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    directory = args.out.parent/(args.out.stem+'.fixtures')
    directory.mkdir(exist_ok=False)
    began = time.monotonic()
    result = dict(schema='WJH27InitialFieldsIndependentCheckV1', status='FAIL', trials=[],
        source=po.identity(__file__), coefficient_fixture_source=po.identity(physical_oracle_test.__file__),
        basis_source=po.identity(po.fixtures.__file__), fixture_directory=str(directory),
        scope='Packet role/identity/inventory/value transfer only; no physical solution, time history, admissibility or H27 qualification')
    try:
        ref = SimpleNamespace(Kin=.0088, H=.15, h=.108, top=.5, d0=.0875, A=.00035, tau=.002601)
        state, _ = physical_oracle_test.fixture(ref)
        state.update(step=0, time=0., dt=.731, basis='Modified_A', mms_variant='Historical',
                     source_identity=po.identity(physical_oracle_test.__file__)['sha256'])
        # Deliberately unrelated temporal metadata must never become a history.
        state['input_history'] = dict(step=999, dt=.937, time=936.063, current='not imported')
        candidate = dict(schema='WJH27SteadyCandidateV1', qualified_reference=False, state=state)
        registered = dict(schema='WJH27RegisteredCoupledRunV2', normal_main=True, factory='WJNativeCoupledV2',
                          steady_mode=True, steady_solver=dict(status='CONVERGED', converged=True), state=copy.deepcopy(state))
        packets = [('candidate', candidate, True), ('registered', registered, True)]
        reordered = copy.deepcopy(candidate); reordered['state']['elements'].reverse()
        packets.append(('reordered', reordered, True))
        def mutate(name, operation, base=candidate):
            value = copy.deepcopy(base); operation(value); packets.append((name, value, False))
        mutate('missing_basis', lambda p: p['state'].pop('basis'))
        mutate('wrong_basis', lambda p: p['state'].__setitem__('basis', 'Modified_B'))
        mutate('missing_identity', lambda p: p['state'].pop('source_identity'))
        mutate('wrong_coefficient_hash', lambda p: p['state'].__setitem__('state_sha256', 'f'*64))
        mutate('old_time', lambda p: p['state'].__setitem__('time', .1))
        mutate('old_step', lambda p: p['state'].__setitem__('step', 1))
        mutate('wrong_coefficient_size', lambda p: p['state']['coefficients']['k'].pop())
        mutate('nonfinite_coefficient', lambda p: p['state']['coefficients']['k'].__setitem__(0, float('nan')))
        mutate('invalid_sign', lambda p: p['state']['elements'][0]['us'].__setitem__(0, .5))
        mutate('repeated_local_id', lambda p: p['state']['elements'][0]['u'].__setitem__(0, p['state']['elements'][0]['u'][1]))
        mutate('out_of_range_id', lambda p: p['state']['elements'][0]['u'].__setitem__(0, p['state']['sizes'][0]))
        def unused(p):
            p['state']['sizes'][2] += 1; p['state']['coefficients']['k'].append(0.)
            p['state']['state_sha256'] = state_hash(p['state'])
        mutate('unused_global_id', unused)
        mutate('duplicate_geometry', lambda p: p['state']['elements'][1].__setitem__('geometry', p['state']['elements'][0]['geometry']))
        mutate('overlap', lambda p: p['state']['elements'][2]['bounds'].__setitem__(0, 0.))
        mutate('gap', lambda p: p['state']['elements'][1]['bounds'].__setitem__(2, ref.H+.01))
        mutate('curved_geometry', lambda p: p['state']['geometry'].__setitem__('curve', .01))
        mutate('missing_outer_epsilon', lambda p: (p['state']['elements'][1].__setitem__('epsilon', []), p['state']['elements'][1].__setitem__('epsilon_signs', [])))
        def inner(p):
            for key in ('epsilon', 'epsilon_signs'):
                p['state']['elements'][0][key] = list(p['state']['elements'][1][key])
        mutate('invented_inner_epsilon', inner)
        mutate('discontinuous_trace', lambda p: p['state']['elements'][0]['us'].__setitem__(1, -p['state']['elements'][0]['us'][1]))
        mutate('candidate_claims_convergence', lambda p: p.__setitem__('qualified_reference', True))
        mutate('registered_not_steady', lambda p: p.__setitem__('steady_mode', False), registered)
        mutate('registered_not_converged', lambda p: p['steady_solver'].__setitem__('converged', False), registered)
        mutate('contradictory_registered_status', lambda p: p['steady_solver'].__setitem__('status', 'NOT_CONVERGED'), registered)
        mutate('non_normal_factory', lambda p: p.__setitem__('normal_main', False), registered)
        valid_checks = {}
        for name, packet, expected in packets:
            remaining = args.seconds-(time.monotonic()-began)
            po.require(remaining > 0, 'InitialFields controls exhausted explicit wall allowance')
            path = directory/(name+'.json'); output = directory/(name+'.probe.json')
            path.write_text(json.dumps(packet, indent=2, allow_nan=True)+'\n')
            run = subprocess.run([str(args.probe.resolve()), str(path), str(output)], capture_output=True, text=True, timeout=remaining, check=False)
            (directory/(name+'.stdout')).write_text(run.stdout)
            (directory/(name+'.stderr')).write_text(run.stderr)
            trial = dict(name=name, expected_accepted=expected, returncode=run.returncode, input=po.identity(path), stderr=run.stderr)
            result['trials'].append(trial)
            po.require((run.returncode == 0) == expected, 'incorrect packet acceptance for '+name+': '+run.stderr)
            if expected:
                native = json.loads(output.read_text())
                valid_checks[name] = check_samples(packet, native, path)
                trial['output'] = po.identity(output)
            else:
                po.require(not output.exists(), 'rejected constructor left a positive probe artifact')
        result.update(status='PASS', valid_packets=valid_checks,
                      rejected_packet_count=sum(not item[2] for item in packets))
    except Exception as error:
        result.update(exception=type(error).__name__, reason=str(error))
    result['elapsed_seconds'] = time.monotonic()-began
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(result, allow_nan=False))
    if result['status'] != 'PASS': raise SystemExit(2)


if __name__ == '__main__':
    main()
