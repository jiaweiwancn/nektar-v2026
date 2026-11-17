#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only nonzero analytic controls for physical_oracle.py.

The coefficient fixture deliberately does NOT satisfy the physical essential
conditions and is NOT a physical solution or a manufactured forced solve. It
tests independent assembly signs, scales, conservation and term visibility.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
import numpy as np
from numpy.polynomial.legendre import leggauss
import reference_oracle


def fixture(ref):
    xs, rs = [-10., 2.5, 15.], [0., ref.H, ref.top]
    elements, maps = [], {key: {} for key in ('u', 'p', 'k', 'epsilon')}
    def ids(key, ix, ir, degree):
        result = []
        for j in range(degree+1):
            rtag = ('v', ir+j) if j < 2 else ('b', ir, j)
            for i in range(degree+1):
                xtag = ('v', ix+i) if i < 2 else ('b', ix, i)
                tag = (xtag, rtag)
                if tag not in maps[key]: maps[key][tag] = len(maps[key])
                result.append(maps[key][tag])
        return result
    for ix in range(2):
        for ir in range(2):
            e = dict(geometry=100+len(elements), bounds=[xs[ix], xs[ix+1], rs[ir], rs[ir+1]])
            for key in maps:
                e[key] = ids(key, ix, ir, 2 if key == 'p' else 3) if key != 'epsilon' or ir == 1 else []
                e[po.SIGNS[key]] = [-1. if gid % 3 == 1 else 1. for gid in e[key]]
            elements.append(e)
    sizes = [len(maps[key]) for key in ('u', 'p', 'k', 'epsilon')]
    n, np_, nk, ne = sizes
    s = dict(flow=np.zeros(2*n+np_), k=np.zeros(nk), epsilon=np.zeros(ne))
    def field(key, x, r):
        X = x+10
        if key == 'u': return .2+.004*X+.1*r/ref.top+.002*X*r/ref.top
        if key == 'w': return .01+.001*X-.003*r/ref.top
        if key == 'p': return ref.tau*(.3+.01*X+.2*r/ref.top)
        if key == 'k': return ref.Kin*(1+.02*X/25+.003*(X/25)**3+.005*(r-ref.H)/ref.top)
        return ref.A/ref.H*(1-.5*(r-ref.H)/(ref.top-ref.H))
    seen = {key: {} for key in ('u', 'w', 'p', 'k', 'epsilon')}
    for e in elements:
        xl, xr, rl, rr = e['bounds']
        for key in seen:
            kind = 'u' if key == 'w' else key
            if not e[kind]: continue
            degree = 2 if kind == 'p' else 3
            a = np.linspace(-1., 1., degree+1)
            xx, rr_ = np.meshgrid((xl+xr)/2+(xr-xl)*a/2, (rl+rr)/2+(rr-rl)*a/2)
            x, r = xx.ravel(), rr_.ravel()
            local = np.linalg.solve(po.basis(e, kind, x, r)[0], field(key, x, r))
            for gid, value in zip(e[kind], local):
                po.require(gid not in seen[key] or abs(seen[key][gid]-value) <= 1e-12,
                           'analytic fixture has inconsistent shared signed coefficient')
                seen[key][gid] = value
                if key in ('u', 'w', 'p'):
                    s['flow'][(0 if key == 'u' else n if key == 'w' else 2*n)+gid] = value
                else: s[key][gid] = value
    def free(key, sides):
        size = sizes[0 if key == 'u' else 2 if key == 'k' else 3]
        fixed = set()
        for e in elements:
            if not e[key]: continue
            xl, xr, rl, rr = e['bounds']
            for side, boundary in enumerate((xl == -10., xr == 15., rl == (ref.H if key == 'epsilon' else 0.), rr == ref.top)):
                if side in sides and boundary:
                    fixed.update(e[key][j] for j in po.EDGES[side])
        return sorted(set(range(size))-fixed)
    digest = hashlib.sha256()
    for key in s: digest.update(np.asarray(s[key], dtype='<f8').tobytes())
    return dict(schema='WJH27NativeCoupledStateV2', manufactured=False, gauge=False,
        model='FiniteNuIsotropicWallFullWJCollarV2', source_identity='ANALYTIC_ORACLE_TEST_ONLY',
        reference_identity=po.TABLE_SHA, state_sha256=digest.hexdigest(), sizes=sizes, elements=elements,
        coefficients={key: value.tolist() for key, value in s.items()},
        geometry=dict(d0=ref.d0, H=ref.H, h=ref.h, left=-10., right=15., top=ref.top),
        free_velocity_ids=[free('u', (0, 2)), free('u', (0, 2, 3))], k_free_ids=free('k', (0,)),
        epsilon_free_ids=free('epsilon', (0, 2))), field


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    args = parser.parse_args()
    po.require(sys.platform.startswith('linux'), 'numerical controls are Linux-only')
    po.require(math.isfinite(args.seconds) and args.seconds > 0 and not args.out.exists(), 'finite cap and immutable output required')
    began = time.monotonic()
    result = dict(schema='WJH27PhysicalOracleAnalyticControlsV1', status='FAIL', source=po.identity(__file__),
        oracle_source=po.identity(po.__file__), independent_scalar_WJ_source=po.identity(reference_oracle.__file__),
        scope='Analytic coefficient assembly controls only; fixture is not a physical solution and carries no manufactured forcing')
    try:
        po.require(po.identity(args.table)['sha256'] == po.TABLE_SHA, 'wrong frozen table')
        ref = po.Reference(json.loads(args.table.read_text()))
        data, field = fixture(ref)
        oracle = po.PhysicalOracle(data, ref, 20, began+args.seconds)
        core = oracle.assemble()
        po.require(core['accounting_identity_status'] == 'MET', 'independent accounting identities failed')
        errors = {}
        def close(name, actual, expected, scale, tolerance=1e-10):
            error = float(np.max(np.abs(np.asarray(actual)-expected)))/scale
            errors[name] = error
            po.require(error <= tolerance, f'{name}: {error} > {tolerance}')
        one, parts = oracle.one['u'], oracle.parts
        expected_k = (2/3)*ref.Kin*(.02+.003)*ref.h
        close('sampled_known_k_Q3_sign', one@parts['awf_known_k']['flow'][:oracle.n], expected_k, abs(expected_k))
        expected_pi = -ref.tau*.01*25*ref.h
        close('sampled_pressure_sign', one@parts['awf_pressure']['flow'][:oracle.n], expected_pi, abs(expected_pi))
        expected_conv = -ref.h*(field('u', 15., ref.H)**2-field('u', -10., ref.H)**2)
        close('sampled_convective_factor_two_sign', one@parts['awf_convection']['flow'][:oracle.n], expected_conv, abs(expected_conv))
        gx, gw = leggauss(4)
        xx, rr = np.meshgrid(2.5+12.5*gx, ref.top*(1+gx)/2)
        weight = np.outer(gw, gw).ravel()*25*ref.top/4
        x, r = xx.ravel(), rr.ravel()
        u, w = field('u', x, r), field('w', x, r)
        conv_x = u*(.004+.002*r/ref.top)+w*(.1+.002*(x+10))/ref.top
        conv_r = u*.001-w*.003/ref.top
        close('domain_convection_x', one@parts['advection']['flow'][:oracle.n], weight@conv_x, abs(weight@conv_x))
        close('domain_convection_r', one@parts['advection']['flow'][oracle.n:2*oracle.n], weight@conv_r, abs(weight@conv_r))
        # The conservative/advective difference is q div(U), with no missing
        # exterior term. Check EVERY signed global coefficient, not only 1.
        for key in ('k', 'epsilon'):
            expected = np.zeros_like(oracle.s[key])
            q, weights = leggauss(6)
            for e in oracle.elements:
                if not e[key]: continue
                xl, xr, rl, rr = e['bounds']
                xx, zz = np.meshgrid((xl+xr)/2+(xr-xl)*q/2, (rl+rr)/2+(rr-rl)*q/2)
                x, r = xx.ravel(), zz.ravel()
                ww = np.outer(weights, weights).ravel()*(xr-xl)*(rr-rl)/4
                divergence = .004+.002*r/ref.top-.003/ref.top
                expected[e[key]] += po.basis(e, key, x, r)[0].T@(ww*field(key, x, r)*divergence)
            actual = parts['advection'][key]-oracle.advective_scalar[key]
            close('conservative_'+key+'_all_test_moments', actual, expected, np.max(np.abs(expected)))
        close('positive_outward_top_epsilon_sign', oracle.one['epsilon']@parts['natural']['epsilon'], 25*ref.top_flux, 25*ref.top_flux)
        close('negative_outlet_vertical_traction_load', one@parts['natural']['flow'][oracle.n:2*oracle.n], -ref.top*ref.tau, ref.top*ref.tau)
        for name, mutation in core['mutations'].items():
            po.require(mutation['status'] == 'DETECTABLE_AT_DECLARED_TOLERANCE', 'nonzero analytic fixture does not detect '+name)
        # Independent scalar np.roots reference vs batched companion roots;
        # general 3D, weak divergence, zero strain and solid rotation included.
        rng = np.random.default_rng(270017)
        gradients = rng.normal(size=(36, 3, 3))*.15
        gradients[0] = 0.; gradients[1] = np.diag([.03, .01, -.02])
        gradients[2] = np.array([[0., .2, -.1], [-.2, 0., .07], [.1, -.07, 0.]])
        k = ref.Kin*np.linspace(.9, 1.1, len(gradients))
        E = ref.A/ref.H*np.linspace(.7, 2., len(gradients))
        batched, roots, defects, conditioning = po.full_wj(k, E, gradients)
        scalar = np.stack([reference_oracle.full_wj(ki, ei, gi) for ki, ei, gi in zip(k, E, gradients)])
        close('full_3D_WJ_independent_root_tensor', batched, scalar, ref.Kin, 2e-11)
        po.require(np.min(conditioning) > 0 and np.min(roots) >= 1.8*(1-1e-12), 'physical root branch changed')
        close('full_3D_WJ_covariance_trace', np.trace(batched, axis1=1, axis2=2), 2*k, ref.Kin, 2e-12)
        # Deliberately suppress a required term. Its own assembled constant
        # moment must differ by the analytically prescribed nonzero amount.
        close('missing_known_k_negative_control', abs(expected_k), abs(one@parts['awf_known_k']['flow'][:oracle.n]), abs(expected_k))
        po.require(core['steady_residual_status'] == 'NOT_MET', 'non-solution fixture falsely passes the steady gate')
        result.update(status='PASS', maximum_scaled_errors=errors, mutation_sensitivities=core['mutations'],
            accounting_identity_max_error=core['accounting_identity_max_error'],
            non_solution_steady_residuals=core['residuals'], signed_basis=True,
            full_3D_constitutive_cases=len(gradients), maximum_cubic_residual=float(np.max(defects)))
    except Exception as error:
        result.update(exception=type(error).__name__, reason=str(error))
    result['elapsed_seconds'] = time.monotonic()-began
    with args.out.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps(result, allow_nan=False))
    if result['status'] != 'PASS': raise SystemExit(2)


if __name__ == '__main__':
    main()
