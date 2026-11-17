#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent representable PolynomialScalarsV1 steady reference.

Constructs coefficients from continuous SymPy polynomials and native signed
inventories, never from the accepted coefficient values or a nonlinear Picard
iteration. Original independently derived continuous manufactured loads remain
unchanged. This verifies a known same-space steady solution and field agreement;
it is not a uniqueness, temporal, physical-reference, or full H27 proof.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import sys
import time

import coupled_oracle as co
import coupled_mms_check as mms_check
import numpy as np
import scipy.linalg as la
import sympy as sp
from numpy.polynomial.legendre import leggauss

TABLE_SHA = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'
FIELD_LIMIT = 1e-8
ORDERS = ((20, 24), (24, 32))


def identity(path):
    path = Path(path)
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def object_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     allow_nan=False).encode()).hexdigest()


def target_polynomials(m):
    """Derive t=0 fields independently; retain exact decimal input parameters."""
    x, r = sp.symbols('x r', real=True)
    rational = lambda value: sp.Rational(str(value))
    H, R, Kin, cl = map(rational, (m.H, m.R, m.Kin, m.cl))
    s = sp.sqrt(Kin)*(1+sp.Rational(8, 1000)*x)
    psi = sp.Rational(22, 1000)*r*r/2+sp.Rational(5, 100)*x*x*(1-x)*r*r*(R-r)
    expressions = dict(u=sp.diff(psi, r), w=-sp.diff(psi, x),
        Pi=sp.Rational(2, 10000)*(x*x+r*r/R**2-sp.Rational(2, 3)),
        k=s*s*(1+sp.Rational(5, 1000)*r*(r-H)/(H*R)), epsilon=s**3/(cl*H))
    expressions = {key: sp.expand(value) for key, value in expressions.items()}
    degrees = {key: [int(sp.Poly(value, x, r).degree(v)) for v in (x, r)]
               for key, value in expressions.items()}
    co.require(all(max(degree) <= (2 if key == 'Pi' else 3) for key, degree in degrees.items()),
               'continuous target is not representable in the declared tensor spaces')
    div = sp.simplify(sp.diff(expressions['u'], x)+sp.diff(expressions['w'], r))
    mean = sp.simplify(sp.integrate(expressions['Pi'], (x, 0, 1), (r, 0, R)))
    trace_k = sp.simplify(expressions['k'].subs(r, H)-s*s)
    trace_e = sp.simplify(expressions['epsilon'].subs(r, H)-s**3/(cl*H))
    co.require(div == mean == trace_k == trace_e == 0, 'symbolic divergence/gauge/Gamma identity failed')
    lower_factor = 1-sp.Rational(5, 1000)*H/(4*R)
    co.require(lower_factor > 0 and Kin > 0 and cl > 0 and H > 0,
               'symbolic target positivity branch is invalid')
    functions = {key: sp.lambdify((x, r), [value, sp.diff(value, x), sp.diff(value, r)],
                                 modules='numpy', cse=True) for key, value in expressions.items()}
    symbolic = dict(variables=['x', 'model_radius_r'], frozen_time=0.,
        expressions={key: sp.sstr(value) for key, value in expressions.items()},
        tensor_degrees=degrees, exact_divergence=str(div), exact_pressure_integral=str(mean),
        exact_K_H_minus_s_squared=str(trace_k), exact_epsilon_H_minus_s_cubed_over_clH=str(trace_e),
        positive_s_interval=[float(sp.sqrt(Kin)), float(sp.sqrt(Kin)*(1+sp.Rational(8, 1000)))],
        positive_k_factor_lower_bound=float(lower_factor),
        Gamma_branch='s>0 on 0<=x<=1, so (s^2)^(3/2)=s^3; no absolute-value or branch substitution')
    symbolic['sha256'] = object_hash(symbolic)
    return functions, symbolic


def jet(function, x, r):
    return np.vstack([np.broadcast_to(np.asarray(v, dtype=float), np.shape(x)) for v in function(x, r)])


def construct_target(data, m, functions, deadline):
    n, p, nk, ne = data['sizes']
    result = dict(flow=np.zeros(2*n+p), k=np.zeros(nk), epsilon=np.zeros(ne))
    seen = {key: {} for key in functions}
    errors = dict(shared_signed_coefficient=0., heldout_value=0., heldout_gradient=0., continuous_MMS_identity=0.)
    for e in sorted(data['elements'], key=lambda value: value['geometry']):
        co.require(time.monotonic() < deadline, 'polynomial reference construction exceeded elapsed cap')
        xl, xr, rl, rr = e['bounds']
        for name in functions:
            key = 'u' if name in ('u', 'w') else 'p' if name == 'Pi' else name
            if not e[key]: continue
            degree = 2 if key == 'p' else 3
            # A small polynomial interpolation problem, not a PDE or mass
            # solve. No native solution coefficients enter this construction.
            nodes = np.linspace(-1., 1., degree+1)
            xx, zz = np.meshgrid((xl+xr)/2+(xr-xl)*nodes/2, (rl+rr)/2+(rr-rl)*nodes/2)
            x, r = xx.ravel(), zz.ravel()
            b = co.basis(e, key, x, r)[0]
            local = la.solve(b, jet(functions[name], x, r)[0], assume_a='gen')
            scale = 1. if name in ('u', 'w') else m.tau0 if name == 'Pi' else m.Kin if name == 'k' else m.Ain/m.H
            for gid, value in zip(e[key], local):
                if gid in seen[name]:
                    error = abs(seen[name][gid]-value)/scale
                    errors['shared_signed_coefficient'] = max(errors['shared_signed_coefficient'], float(error))
                    co.require(error <= 2e-12, 'independent polynomial signed trace conflict')
                else:
                    seen[name][gid] = float(value)
                    if name in ('u', 'w', 'Pi'):
                        offset = 0 if name == 'u' else n if name == 'w' else 2*n
                        result['flow'][offset+gid] = value
                    else:
                        result[name][gid] = value
            # Heldout Gauss7 points plus every endpoint are disjoint from the
            # interpolation interior; check values and first derivatives.
            heldout = np.r_[-1., leggauss(7)[0], 1.]
            xx, zz = np.meshgrid((xl+xr)/2+(xr-xl)*heldout/2, (rl+rr)/2+(rr-rl)*heldout/2)
            x, r = xx.ravel(), zz.ravel()
            exact = jet(functions[name], x, r)
            reconstructed = np.vstack([operator@local for operator in co.basis(e, key, x, r)])
            scale = m.Ain/r if name == 'epsilon' else scale
            errors['heldout_value'] = max(errors['heldout_value'], float(np.max(abs(reconstructed[0]-exact[0])/scale)))
            errors['heldout_gradient'] = max(errors['heldout_gradient'], float(np.max(abs(reconstructed[1:]-exact[1:])/scale)))
            original = m.functions[name](x, r, 0.)
            original = np.vstack([np.broadcast_to(np.asarray(v), x.shape) for v in original[:3]])
            errors['continuous_MMS_identity'] = max(errors['continuous_MMS_identity'], float(np.max(abs(original-exact)/scale)))
    for name, size in (('u', n), ('w', n), ('Pi', p), ('k', nk), ('epsilon', ne)):
        co.require(set(seen[name]) == set(range(size)), 'incomplete independently constructed '+name+' inventory')
    co.require(errors['heldout_value'] <= 2e-12 and errors['heldout_gradient'] <= 2e-11 and
               errors['continuous_MMS_identity'] <= 2e-12, 'polynomial representability or continuous-field identity failed')
    co.require(all(np.all(np.isfinite(v)) for v in result.values()), 'nonfinite independently constructed target')
    return result, errors


def scaled_difference(o, native, reference):
    d = co.combine(native, 1., reference, -1.)
    n = o.nu
    energy = lambda value, matrix, volume: math.sqrt(max(0., float(value@matrix@value))/volume)
    return dict(U=energy(d['flow'][:2*n], o.Mv, o.volume),
        Pi=energy(d['flow'][2*n:], o.Q, o.volume)/o.m.tau0,
        k=energy(d['k'], o.Mk, o.volume)/o.m.Kin,
        epsilon=energy(d['epsilon'], o.Eerror, o.evolume))


def quadrature_delta(o, previous, current, old_mean, s):
    delta = co.combine(current, 1., previous, -1.)
    values = o.dual(delta); continuity = delta['flow'][2*o.nu:]
    values['continuity'] = math.sqrt(max(0., float(continuity@la.cho_solve(o.riesz_p, continuity)))/o.volume)
    values['gauge'] = abs(float((o.mean-old_mean)@s['flow'][2*o.nu:]))/(.01*o.volume)
    return values


def pointwise_gamma_error(o, s):
    result = 0.
    for e, side in o.faces['epsilon']:
        if side != 2: continue
        xl, xr, _, _ = e['bounds']
        nodes = np.r_[-1., leggauss(12)[0], 1.]
        x = (xl+xr)/2+(xr-xl)*nodes/2; r = np.full(len(x), o.m.H)
        K = o.trace(e, 'k', x, s)
        co.require(np.all(K > 0), 'actual H trace lost the positive Gamma branch')
        E = co.basis(e, 'epsilon', x, r)[0]@s['epsilon'][e['epsilon']]
        result = max(result, float(np.max(abs(E-K**1.5/(o.m.cl*o.m.H)))/(o.m.Ain/o.m.H)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'table', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seconds', type=float, default=180.)
    args = parser.parse_args()
    co.require(sys.platform.startswith('linux'), 'authoritative numerical reference construction is Linux-only')
    co.require(math.isfinite(args.seconds) and 0 < args.seconds <= 300, 'explicit bounded elapsed allowance must be in (0,300] seconds')
    co.require(not args.out.exists(), 'reference evidence output is immutable')
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    began = time.monotonic(); deadline = began+args.seconds
    report = dict(schema='WJH27IndependentPolynomialSteadyReferenceV1', status='NOT_MET',
        full_H27_qualification=False, physical_reference_gate_pass=False, temporal_order_qualified=False,
        source=identity(__file__), independent_oracle_source=identity(co.__file__),
        continuous_MMS_source=identity(mms_check.__file__), independent_WJ_source=identity(co.reference_oracle.__file__),
        independent_Jacobi_basis_source=identity(Path(co.__file__).resolve().parents[1]/'h26'/'fixtures.py'),
        native=identity(args.native), table=identity(args.table), wall_cap_seconds=args.seconds,
        memory_cap_bytes=cap, workers=1, field_RMS_limit=FIELD_LIMIT, residual_limits=co.LIMITS,
        scope='Exactly representable steady PolynomialScalarsV1 same-space reference at alpha=0, t=0; no nonlinear uniqueness or whole-H27 claim',
        prior_Picard_failure='Preserved as a distinct verifier-iteration failure. This construction neither runs nor repairs that iteration and does not relabel its result.',
        forcing_route='Original continuous SymPy manufactured fields integrated by parts by the independent oracle; analytic time derivatives removed at t=0. No projected A_h*target forcing.')
    try:
        packet = json.loads(args.native.read_text()); table = json.loads(args.table.read_text())
        co.require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True and
                   packet.get('factory') == 'WJNativeCoupledV2' and packet.get('steady_mode') is True,
                   'input is not an actual registered steady native run')
        co.require(packet['steady_solver']['converged'] is True and packet['steady_solver']['status'] == 'CONVERGED',
                   'native steady run does not report convergence')
        data = packet['state']
        co.require(data['schema'] == 'WJH27NativeCoupledStateV2' and data['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and
                   data.get('basis') == 'Modified_A' and data.get('mms_variant') == 'PolynomialScalarsV1' and
                   data['manufactured'] is True and data['gauge'] is True, 'wrong native state/model/basis/MMS variant')
        co.require(data['step'] == 0 and data['time'] == 0. and packet.get('physical_time_advanced', 0.) == 0.,
                   'steady reference requires unchanged step-zero/time-zero clock')
        co.require(report['table']['sha256'] == TABLE_SHA == data['reference_identity'] and
                   table['schema'] == 'WJH27V2ReferenceTablesV1' and table['status'] == 'PASS', 'qualified reference-table identity mismatch')
        co.require(len(data['source_identity']) == 64 and all(c in '0123456789abcdef' for c in data['source_identity']),
                   'native source identity must be a SHA256')
        native = co.state(data['coefficients'], data['sizes'])
        co.require(co.state_hash(native) == data['state_sha256'] == data['previous_sha256'], 'steady native state byte identity differs')
        cached = co.state(data['current_advection'], data['sizes'])
        cached_previous = co.state(data['previous_advection'], data['sizes'])
        co.require(co.state_hash(cached) == data['advection_sha256'] == data['previous_advection_sha256'] ==
                   co.state_hash(cached_previous), 'steady advection identities differ')
        co.require(np.all(cached['flow'][2*data['sizes'][0]:] == 0.), 'native explicit cache has a pressure load')
        m = mms_check.IndependentMMS(table, data.get('mms_body_amplitude', 0.), 'PolynomialScalarsV1')
        functions, symbolic = target_polynomials(m)
        reference, representation = construct_target(data, m, functions, deadline)
        inventory = {key: data[key] for key in ('sizes', 'geometry', 'elements', 'free_velocity_ids', 'k_free_ids', 'epsilon_free_ids', 'basis')}
        report.update(native_source_identity=data['source_identity'], native_state_sha256=data['state_sha256'],
            native_steady_solver=packet['steady_solver'], symbolic_target=symbolic, representation_errors=representation,
            reference_state_sha256=co.state_hash(reference), reference_coefficients={key: value.tolist() for key, value in reference.items()},
            signed_inventory_sha256=object_hash(inventory), native_reported_residuals=data['residuals'],
            mms_body_amplitude=m.body_amplitude, physical_time_advanced=0., alpha=0.,
            field_scales=dict(U='combined velocity RMS / 1', Pi='pressure RMS / tau0; native pressure is not shifted',
                              k='k RMS / Kin', epsilon='outer epsilon RMS weighted by r/Ain'))
        n, p, k, e = data['sizes']
        dense_bytes = 8*(8*n*n+3*p*p+4*k*k+6*e*e+8*n*p+4*n*k+2*p*k)
        co.require(dense_bytes < .75*cap, 'bounded dense independent MMS oracle exceeds declared memory reserve')
        report['estimated_dense_bytes'] = dense_bytes
        assemblies, residuals, means = [], [], []
        for nx, nr in ORDERS:
            o = co.Oracle(data, m, nx, nr, deadline)
            # ONE shared load per quadrature, constructed only from continuous
            # fields. The native state enters this API solely to add integration
            # breakpoints; it does not determine any manufactured source value.
            forcing = o.forcing(native, 0., steady=True)
            row = dict(quadrature=[nx, nr], continuous_forcing_sha256=co.state_hash(forcing))
            res = {}
            for role, s in (('native', native), ('reference', reference)):
                adv = o.advection(s, 0.)[0]
                res[role], norms = o.residual(s, 0., 0., co.combine(forcing, 1., adv, -1.))
                row[role] = dict(residuals=norms, residual_sha256=co.state_hash(res[role]),
                    current_advection_sha256=co.state_hash(adv), essential_Gamma_errors=o.trace_error(s, 0.),
                    dependent_H_pointwise_error=pointwise_gamma_error(o, s),
                    Gamma_note='Essential coefficient gate tests the current endpoint-constrained Gamma projection; pointwise error is reported separately.',
                    actual_pressure_mean=float(o.mean@s['flow'][2*o.nu:])/o.volume)
                if role == 'native': row['stored_current_advection_error'] = o.dual(co.combine(cached, 1., adv, -1.))
            row['reference_scaled_RMS_difference'] = scaled_difference(o, native, reference)
            assemblies.append(row); residuals.append(res); means.append(o.mean.copy())
            report['independent_assemblies'] = assemblies
            if len(assemblies) < len(ORDERS): del o
        differences = {role: quadrature_delta(o, residuals[0][role], residuals[1][role], means[0], s)
                       for role, s in (('native', native), ('reference', reference))}
        report['quadrature_difference_dual'] = differences
        residual_ok = all(value <= co.LIMITS[key] for row in assemblies for role in ('native', 'reference')
                          for key, value in row[role]['residuals'].items())
        boundary_ok = all(value <= 1e-10 for row in assemblies for role in ('native', 'reference')
                         for value in row[role]['essential_Gamma_errors'].values())
        quadrature_ok = all(value <= .2*co.LIMITS[key] for row in differences.values() for key, value in row.items())
        fields_ok = all(value <= FIELD_LIMIT for row in assemblies for value in row['reference_scaled_RMS_difference'].values())
        cache_ok = all(value <= 1e-10 for row in assemblies for value in row['stored_current_advection_error'].values())
        report.update(current_and_reference_residual_pass=residual_ok, actual_essential_Gamma_and_gauge_pass=boundary_ok and
                      all(row[role]['residuals']['gauge'] <= co.LIMITS['gauge'] for row in assemblies for role in ('native', 'reference')),
            independent_quadrature_pass=quadrature_ok, unchanged_reference_scaled_field_agreement_pass=fields_ok,
            stored_current_advection_pass=cache_ok, quadrature_reserve_fraction=.2,
            target_reference_binding_sha256=object_hash(dict(symbolic=symbolic['sha256'],
                coefficients=co.state_hash(reference), inventory=object_hash(inventory), table=TABLE_SHA,
                continuous_source=report['continuous_MMS_source']['sha256'])))
        passed = residual_ok and boundary_ok and quadrature_ok and fields_ok and cache_ok
        report['same_discrete_system_polynomial_reference_gate_pass'] = passed
        report['status'] = 'PASS_REPRESENTABLE_STEADY_REFERENCE_ONLY' if passed else 'NOT_MET'
        report['reason'] = ('Independent continuous polynomial target and accepted native state tested against the same original steady manufactured problem. '
                            'The previous optional Picard failure remains distinct and unchanged; no physical or whole-H27 qualification is inferred.')
    except Exception as error:
        report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-began
    report['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'reason', 'elapsed_seconds', 'exception',
                     'same_discrete_system_polynomial_reference_gate_pass') if key in report}))
    if report['status'] != 'PASS_REPRESENTABLE_STEADY_REFERENCE_ONLY': raise SystemExit(2)


if __name__ == '__main__':
    main()
