#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent same-inventory, same-time numerical stress-split comparison."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
            'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
from physical_oracle import basis


def identity(path):
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def check(paths):
    assert sys.platform.startswith('linux')
    runs = [json.loads(path.read_text()) for path in paths]
    states = [run['state'] for run in runs]
    for key in ('model', 'reference_identity', 'time', 'step', 'dt', 'sizes',
                'geometry', 'elements', 'free_velocity_ids', 'k_free_ids',
                'epsilon_free_ids'):
        assert states[0][key] == states[1][key], 'incompatible ' + key
    assert {s['stress_split'] for s in states} == {'IsotropicMuV1', 'PositiveShearV1'}
    assert all(not s['manufactured'] and s['step'] > 0 for s in states)
    n, np_, _, _ = states[0]['sizes']
    kin, cl, H = .008809015145515192, 2.3804029452632984, .15
    amplitude = kin**1.5 / cl
    maxima = dict(U=0., Pi=0., k=0., epsilon=0.)
    # |Modified_A mode| <= 1 on each owned reference cell. Absolute sums of
    # signed coefficient differences are rigorous field MAX upper bounds.
    bounds = dict(U=0., Pi=0., k=0., epsilon=0.)
    sample_nodes = np.unique(np.r_[np.linspace(-1, 1, 17),
                                  np.polynomial.legendre.leggauss(20)[0]])
    for e in states[0]['elements']:
        xl, xr, rl, rr = e['bounds']
        xx, radial = np.meshgrid((xl+xr)/2+(xr-xl)*sample_nodes/2,
                                 (rl+rr)/2+(rr-rl)*sample_nodes/2)
        x, r = xx.ravel(), radial.ravel()
        udelta = []
        for field, key, scale in (('u', 'U', 1.), ('v', 'U', 1.),
                                  ('p', 'Pi', .002601), ('k', 'k', kin),
                                  ('epsilon', 'epsilon', amplitude/rr)):
            if field == 'epsilon' and not e['epsilon']:
                continue
            native_key = 'u' if field == 'v' else field
            ids = np.asarray(e[native_key])
            if field in ('u', 'v', 'p'):
                offset = 0 if field == 'u' else n if field == 'v' else 2*n
                delta = (np.asarray(states[0]['coefficients']['flow']) -
                         np.asarray(states[1]['coefficients']['flow']))[offset+ids]
            else:
                delta = (np.asarray(states[0]['coefficients'][field]) -
                         np.asarray(states[1]['coefficients'][field]))[ids]
            signs = np.asarray(e[{'u': 'us', 'v': 'us', 'p': 'ps',
                                 'k': 'ks', 'epsilon': 'epsilon_signs'}[field]])
            bound = float(np.sum(abs(signs*delta))/scale)
            bounds[key] = max(bounds[key], bound)
            value = basis(e, native_key, x, r)[0] @ delta
            if field in ('u', 'v'):
                udelta.append(value)
            else:
                point_scale = amplitude/r if field == 'epsilon' else scale
                maxima[key] = max(maxima[key], float(np.max(abs(value)/point_scale)))
        maxima['U'] = max(maxima['U'], float(np.max(np.hypot(*udelta))))
    # Component-sum bounds imply the vector bound after a conservative sqrt(2).
    bounds['U'] *= np.sqrt(2)
    assert all(value <= 1e-8 for value in bounds.values()), bounds
    assert all(max(s['residuals']) <= s['nonlinear_tolerance'] for s in states)
    return dict(status='PASS_SAME_TIME_CURRENT_SPLITS_ONLY',
                full_coupled_gate_pass=False, sampled_maxima=maxima,
                rigorous_modified_basis_sum_bounds=bounds, tolerance=1e-8,
                scales=dict(U=1, Pi=.002601, k=kin, epsilon='Ain/r'),
                steps=states[0]['step'], physical_time=states[0]['time'],
                complete_seconds=[run['complete_advancing_seconds'] for run in runs],
                nonlinear_iterations=[s['nonlinear_iterations'] for s in states],
                residuals=[s['residuals'] for s in states],
                scope='Numerical split equivalence on this actual advancing state; no long-time accuracy or cost claim')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--old', type=Path, required=True)
    parser.add_argument('--new', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = dict(schema='WJH27CurrentSplitComparisonV1', status='FAIL',
                  inputs=[identity(args.old), identity(args.new)],
                  helper=identity(Path(__file__)))
    error = None
    try:
        report.update(check([args.old, args.new]))
    except Exception as caught:
        error = caught
        report.update(reason=str(caught), exception=type(caught).__name__)
    with args.out.open('x') as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if error is not None:
        raise error
    print(json.dumps(report, sort_keys=True))


if __name__ == '__main__':
    main()
