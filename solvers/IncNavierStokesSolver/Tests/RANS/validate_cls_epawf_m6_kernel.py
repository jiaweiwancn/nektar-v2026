#!/usr/bin/env python3
"""Independent transition/profile kernel oracle for the revised note's M6.

The reference evaluates explicit physical formulae with 80-digit mpmath and
uses mpmath's numerical differentiation. It does not import the C++ jets,
native gates, native quadrature, exported profile cutoffs, or solver helpers to
construct its answers. Native profile K derivatives are central differences
with a declared convergence sequence, not claimed exact analytic derivatives.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

sys.dont_write_bytecode = True
import mpmath as mp
import numpy as np
from scipy.integrate import quad

mp.mp.dps = 80


def real(value):
    return mp.mpf(str(value))


def decreasing(q, lo, hi):
    if q <= lo:
        return mp.mpf(1)
    if q >= hi:
        return mp.mpf(0)
    x = (q-lo)/(hi-lo)
    return 1-10*x**3+15*x**4-6*x**5


def transition(inputs, controls):
    ell, margin, reynolds = inputs[:3]
    s = inputs[3:]
    squares = [sum(x*x for x in s[:3]), sum(x*x for x in s[3:6]),
               s[6]**2, s[7]**2, s[8]**2,
               sum(x*x for x in s[9:12]), sum(x*x for x in s[12:15])]
    start, end, on, off, moff, mon, roff, ron = controls[:8]
    departure = sum(weight*square for weight, square in zip(controls[8:], squares))
    spatial = decreasing(ell, start, end)
    state = decreasing(departure, on**2, off**2)
    gm = 1-decreasing(margin, moff, mon)
    gr = 1-decreasing(reynolds, roff, ron)
    admissibility = gm*gr
    return [spatial*state*admissibility, spatial, state, admissibility, gm, gr, departure]


def physical_sensors(inputs):
    vectors = [inputs[3*j:3*j+3] for j in range(11)]
    n, uh, uh0, ua, ua0, ub, ub0, pressure, pressure0, convection, convection0 = vectors
    K, H, nu, Q, ka, ea, ea0 = inputs[33:]
    scale2 = K/Q+(nu/H)**2
    scale = mp.sqrt(scale2)

    def sub(a, b):
        return [x-y for x, y in zip(a, b)]

    def tangent(a):
        normal = sum(x*y for x, y in zip(a, n))
        return [x-normal*y for x, y in zip(a, n)]

    da, db = sub(ua, ua0), sub(ub, ub0)
    matching = [x/scale for x in tangent(sub(uh, uh0))]
    difference = [x/scale for x in tangent(sub(db, da))]
    normal = sum(x*y for x, y in zip(n, da))/scale
    return matching+difference+[normal, (ka-K)/K, (ea-ea0)/ea0]+[
        H*x/scale2 for x in sub(pressure, pressure0)]+[
        H*x/scale2 for x in sub(convection, convection0)]


def directional(function, point, direction):
    point, direction = list(map(real, point)), list(map(real, direction))

    def path(t):
        return function([value+t*delta for value, delta in zip(point, direction)])

    values = path(mp.mpf(0))
    first = [mp.diff(lambda t: path(t)[j], 0, 1) for j in range(len(values))]
    second = [mp.diff(lambda t: path(t)[j], 0, 2) for j in range(len(values))]
    return [[float(value) for value in result] for result in (values, first, second)]


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}
        self.context = ''

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError((self.context+': ' if self.context else '')+message)

    def near(self, actual, expected, name, tolerance=3e-11):
        actual, expected = np.asarray(actual, float), np.asarray(expected, float)
        self.require(actual.shape == expected.shape and np.all(np.isfinite(actual)) and
                     np.all(np.isfinite(expected)), 'invalid '+name)
        scale = np.maximum(1., np.abs(expected))
        error = float(np.max(np.abs(actual-expected)/scale, initial=0.))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: scaled error {error} > {tolerance}')


def validate_transitions(data, audit):
    cases = {case['name']: case for case in data['transition_cases']}
    audit.require(len(cases) == len(data['transition_cases']), 'duplicate transition cases')
    for case in cases.values():
        audit.context = case['name']
        inputs, controls = case['input'], case['controls']
        audit.require(len(inputs) == 18 and len(controls) == 15, 'transition input shape')
        audit.require(controls[0] < controls[1] and 0 < controls[2] < controls[3] and
                      0 < controls[4] < controls[5] and 0 < controls[6] < controls[7] and
                      min(controls[8:]) > 0, 'invalid exported transition controls')
        expected = directional(lambda x: transition(x[:18], x[18:]), inputs+controls,
                               case['direction']+case['control_direction'])
        audit.near(case['value'], expected[0], 'transition_value')
        audit.near(case['first'], expected[1], 'transition_first')
        audit.near(case['second'], expected[2], 'transition_second', 3e-10)
        audit.near(case['legacy_weight'], expected[0][0], 'legacy_compatibility')
        audit.require(0 <= case['value'][0] <= 1, 'chi outside [0,1]')
    audit.context = 'transition coverage'
    audit.near(cases['equal_reference']['value'][0], 1, 'projected_reference_plateau', 1e-15)
    for j in range(18):
        audit.require(f'input_direction_{j}' in cases, f'missing input derivative {j}')
    for j in range(15):
        audit.require(f'control_direction_{j}' in cases, f'missing control derivative {j}')
    for j in range(7):
        case = cases[f'independent_sensor_{j}']
        audit.near(case['value'][0], 0, 'independent_sensor_off', 1e-15)
        audit.require(case['value'][-1] > 0, 'sensor did not contribute independently')
    cancel = cases['pressure_convection_no_cancellation']
    audit.near(np.asarray(cancel['input'][12:15])+cancel['input'][15:18], [0, 0, 0], 'cancelling_forces')
    audit.require(cancel['value'][-1] > cancel['controls'][3]**2 and cancel['value'][0] == 0,
                  'pressure/convection sensors cancelled')
    audit.near(cases['inadmissible_inactive']['value'][0], 0, 'inadmissible_branch', 1e-15)
    downstream = [case for name, case in cases.items() if name.startswith('downstream_')]
    audit.require(len(downstream) == 3, 'missing one-way downstream controls')
    for case in downstream:
        audit.near(case['value'][0], 0, 'no_downstream_reactivation', 1e-15)
    endpoint_rates = []
    for gate in range(4):
        for endpoint in range(2):
            exact = cases[f'endpoint_{gate}_{endpoint}']
            audit.near([exact['first'][0], exact['second'][0]], [0, 0], 'C2_exact_endpoint', 1e-14)
            approach = [cases[f'approach_{gate}_{endpoint}_{level}'] for level in range(3)]
            first = [abs(case['first'][0]) for case in approach]
            second = [abs(case['second'][0]) for case in approach]
            audit.require(first[2] < first[1] < first[0] and second[2] < second[1] < second[0],
                          'native first/second derivatives do not approach zero')
            first_ratio, second_ratio = first[1]/first[2], second[1]/second[2]
            audit.require(3 < first_ratio < 5 and 1.5 < second_ratio < 2.5,
                          'native C2 endpoint convergence rates')
            endpoint_rates.append(dict(gate=gate, endpoint=endpoint,
                                       first_ratio=first_ratio, second_ratio=second_ratio))
    return dict(cases=len(cases), input_directions=18, control_directions=15,
                endpoint_convergence=endpoint_rates, no_pressure_convection_cancellation=True)


def validate_physical(data, audit):
    cases = {case['name']: case for case in data['physical_cases']}
    audit.require(len(cases) == 41, 'missing physical-input derivative coverage')
    for j in range(40):
        audit.require(f'physical_direction_{j}' in cases, f'missing physical derivative {j}')
    for case in cases.values():
        audit.context = case['name']
        audit.require(len(case['input']) == len(case['direction']) == 40, 'physical sample shape')
        expected = directional(physical_sensors, case['input'], case['direction'])
        audit.near(case['value'], expected[0], 'physical_sensor_value')
        audit.near(case['first'], expected[1], 'physical_sensor_first')
        audit.near(case['second'], expected[2], 'physical_sensor_second', 3e-10)
    audit.context = 'physical coverage'
    equal = list(map(real, data['equal_physical_input']))
    audit.near(data['equal_physical_sensors'], list(map(float, physical_sensors(equal))), 'equal_physical_formula')
    audit.near(data['equal_physical_sensors'], np.zeros(15), 'equal_reference_indicators', 1e-15)
    audit.require(abs(cases['physical_base']['value'][1]) > 0 and
                  abs(cases['physical_base']['value'][4]) > 0, 'second tangential sensor component omitted')
    return dict(cases=len(cases), scalar_input_directions=40,
                geometric_normal_derivatives='algebraic P=I-n*nT extension; geometric first variations are tangent')


def validate_coordinates(data, audit):
    names = set()
    for case in data['coordinates']:
        audit.context = case['name']
        names.add(case['name'])
        kind, dimension, H0, hx, hy, ha, hkx, hky, sx, sy, radius, amplitude, Lx, Ly = case['map_parameters']
        del dimension, H0, hx, hy, ha, hkx, hky, radius, Ly
        x, y, rho = case['chart']
        dx, dy, drho = case['direction']
        del rho, drho
        if case['affine']:
            audit.require(kind in (0, 1), 'affine coordinate on nonflat support')
            origin, covector = np.asarray(case['origin']), np.asarray(case['covector'])
            point = np.array([x, y, sx*x+sy*y])
            normal = np.array([-sx, -sy, 1.])/np.sqrt(1+sx*sx+sy*sy)
            audit.near(covector@normal, 0, 'affine_tangent_covector')
            coordinate = covector@(point-origin)
            gradient = np.array([covector[0]+sx*covector[2], covector[1]+sy*covector[2], 0.])
            second = 0.
        elif kind in (0, 1):
            coordinate = x*np.sqrt(1+sx*sx)
            gradient, second = np.array([np.sqrt(1+sx*sx), 0., 0.]), 0.
        else:
            audit.require(kind == 4 and Lx > 0, 'unsupported declared arc geometry')
            wave = np.pi/Lx
            integrand = lambda t: np.sqrt(1+(.5*amplitude*wave*np.sin(wave*t))**2)
            inside = min(abs(x), Lx)
            integral, error = quad(integrand, 0, inside, epsabs=2e-13, epsrel=2e-13)
            audit.require(error < 1e-10, 'independent arc quadrature unresolved')
            coordinate = np.sign(x)*(integral+abs(x)-inside)
            slope = -.5*amplitude*wave*np.sin(wave*x) if abs(x) < Lx else 0.
            curvature = -.5*amplitude*wave**2*np.cos(wave*x) if abs(x) < Lx else 0.
            gradient = np.array([np.sqrt(1+slope*slope), 0., 0.])
            second = slope*curvature/gradient[0]*dx*dx
        audit.near(case['coordinate'], coordinate, 'wall_coordinate')
        audit.near(case['chart_gradient'], gradient, 'wall_coordinate_gradient')
        audit.near(case['first'], gradient@np.array([dx, dy, 0.]), 'wall_coordinate_first')
        audit.near(case['second'], second, 'wall_coordinate_second')
        audit.near(case['ray_coordinate'], coordinate, 'coordinate_ray_constancy', 1e-14)
    audit.require({'flat', 'plane2d', 'plane3d', 'affine3d'}.issubset(names) and
                  len([name for name in names if name.startswith('cosine_')]) == 7,
                  'missing coordinate geometry/foot coverage')
    return dict(cases=len(names), physical_policy='prescribed geometry only; no flow direction; no downstream reactivation')


def validate_blends(data, audit):
    omission_controls = []
    for case in data['blends']:
        audit.context = 'blend'
        chi, dc, ddc = (real(case[key][0]) for key in ('transition_value', 'transition_first', 'transition_second'))
        p, dp, ddp = (list(map(real, case[key])) for key in ('profile', 'first_profile', 'second_profile'))

        def function(t):
            weight = chi+t*dc+t*t*ddc/2
            ea, ec, p0, production = [value+t*first+t*t*second/2 for value, first, second in zip(p, dp, ddp)]
            return [ea+weight*(ec-ea), weight*(ec-p0), production-ea+weight*(ea-p0)]

        expected = [[float(mp.diff(lambda t: function(t)[j], 0, order)) for j in range(3)] for order in range(3)]
        for name, value in zip(('value', 'first', 'second'), expected):
            audit.near(case[name], value, 'blend_'+name)
        missing = float((p[1]-p[0])*dc)
        source_missing = float((p[0]-p[2])*dc)
        audit.require(abs(missing) > 1e-8 and abs(source_missing) > 1e-8,
                      'chi-gradient omission controls are insensitive')
        omission_controls.append(dict(dissipation_missing_chi_derivative=missing,
                                      k_source_missing_chi_derivative=source_missing))
        audit.near(case['value'][0]*-1+case['value'][1]+float(p[3]), case['value'][2], 'source_assembled_once')
    safe = data['safe_inactive']
    audit.context = 'safe inactive blend'
    audit.require(safe['corrected_profile_evaluated'] is False, 'inactive corrected profile was evaluated')
    audit.near(safe['value'], [.5, 0, .1], 'safe_inactive_baseline')
    audit.near(safe['first']+safe['second'], np.zeros(6), 'safe_inactive_derivatives', 1e-15)
    return dict(cases=len(data['blends']), derivative_omission_controls=omission_controls,
                inactive_undefined_corrected_profile_bypassed=True)


def profile_reference(K, payload):
    nu, H, roughness, rK, rm = map(real, payload['geometry'])
    del rK
    Q, kappa, K0 = map(real, (payload['Q'], payload['kappa'], payload['K0']))
    root = mp.sqrt(K)
    rough_re = roughness*root/(70*nu)
    exponent = max(mp.mpf('.5')-mp.mpf('.4')*rough_re**mp.mpf('.7'),
                   1-mp.mpf('.79')*rough_re**mp.mpf('-.28'))
    zv = nu*mp.mpf('10.7')*(1-rough_re**exponent)/root
    slope = mp.mpf('.09')*mp.mpf('2.55')*root
    lower = max(mp.mpf(0), zv)
    if lower < H:
        d0, d1 = nu+slope*(lower-zv), nu+slope*(H-zv)
        production_budget = (K/Q)**2*(mp.log(d1/d0)+nu/d1-nu/d0)/slope
    else:
        production_budget = mp.mpf(0)
    numerator = (K/Q)**mp.mpf('1.5')/kappa
    minimum_budget = numerator*(1+mp.log(H/rm))
    margin = (production_budget-minimum_budget)/minimum_budget
    reynolds = H*root/nu
    cutoff = H*mp.exp(1-production_budget/numerator)
    inputs = [mp.mpf('-1.5'), margin, reynolds]+[mp.mpf(0)]*15
    inputs[3] = mp.mpf('.045')+mp.mpf('.02')*(K/K0-1)
    inputs[4] = mp.mpf('.012')
    weight = transition(inputs, list(map(real, payload['controls'])))
    dissipation_cutoff = mp.mpf('5.1')*nu/root
    values, profiles = [], []
    for r in map(real, payload['heights']):
        nut = slope*max(mp.mpf(0), r-zv)
        ea = K**mp.mpf('1.5')/(mp.mpf('2.55')*max(r, dissipation_cutoff))
        ec = numerator/max(r, cutoff)
        p0 = (K/Q)**2*nut/(nu+nut)**2
        production = mp.mpf('.27')+mp.mpf('.03')*r+mp.mpf('.02')*K
        profiles.append([ea, ec, p0, production])
        values.extend([ea+weight[0]*(ec-ea), weight[0]*(ec-p0),
                       production-ea+weight[0]*(ea-p0)])
    return dict(K=K, cutoff=cutoff, margin=margin, reynolds=reynolds,
                production_budget=production_budget, minimum_budget=minimum_budget,
                input=inputs, transition=weight, profiles=profiles, values=values, numerator=numerator)


def validate_profiles(payload, audit):
    K0 = real(payload['K0'])
    independent_base = profile_reference(K0, payload)
    controls = payload['controls']
    audit.near(controls[4:6], [float(independent_base['margin']*.5), float(independent_base['margin']*1.5)],
               'profile_margin_controls')
    audit.near(controls[6:8], [float(independent_base['reynolds']*.8), float(independent_base['reynolds']*1.2)],
               'profile_reynolds_controls')
    points = [payload['base']]+[point for level in payload['levels'] for point in (level['minus'], level['plus'])]
    for point in points:
        audit.context = 'profile K='+str(point['K'])
        reference = profile_reference(real(point['K']), payload)
        audit.require(0 < reference['cutoff'] < real(payload['geometry'][-1]), 'independent cutoff inadmissible')
        for key in ('cutoff', 'margin', 'reynolds', 'production_budget', 'minimum_budget'):
            audit.near(point[key], float(reference[key]), 'profile_'+key, 3e-10)
        audit.near(point['input'], list(map(float, reference['input'])), 'profile_transition_inputs', 3e-10)
        audit.near(point['transition'], list(map(float, reference['transition'])), 'profile_transition_values', 5e-10)
        audit.require(len(point['samples']) == 3, 'missing cutoff/interior/trace heights')
        for sample, expected in zip(point['samples'], reference['profiles']):
            audit.near(sample['profile'], list(map(float, expected)), 'actual_production_profile', 3e-10)
        actual_values = np.array([sample['value'] for sample in point['samples']]).reshape(-1)
        audit.near(actual_values, list(map(float, reference['values'])), 'actual_blended_profile', 5e-10)
        audit.near(point['redistribution_integral'], 0, 'ray_zero_net_redistribution', 3e-10)
        audit.near(point['flux_endpoints'], [0, 0], 'ray_redistribution_flux_endpoints', 3e-10)

    def function(K):
        reference = profile_reference(K, payload)
        return [reference['cutoff']]+reference['values']

    audit.context = 'profile current-K derivatives'

    analytic = [[float(mp.diff(lambda K: function(K)[j], K0, order))
                 for j in range(10)] for order in (1, 2)]
    native_base = np.array([payload['base']['cutoff']]+[
        value for sample in payload['base']['samples'] for value in sample['value']])
    levels = []
    for level in payload['levels']:
        h = level['h']
        plus, minus = [np.array([level[key]['cutoff']]+[
            value for sample in level[key]['samples'] for value in sample['value']]) for key in ('plus', 'minus')]
        derivatives = [(plus-minus)/(2*h), (plus-2*native_base+minus)/(h*h)]
        errors = [float(np.max(np.abs(actual-expected)/np.maximum(1, np.abs(expected))))
                  for actual, expected in zip(derivatives, map(np.asarray, analytic))]
        levels.append(dict(h=h, first_error=errors[0], second_error=errors[1],
                           native_first=derivatives[0].tolist(), native_second=derivatives[1].tolist()))
    audit.require(len(levels) == 3 and levels[-1]['h'] < levels[0]['h'], 'missing K derivative refinement')
    audit.require(levels[-1]['first_error'] < 3e-6 and levels[-1]['second_error'] < 3e-5,
                  'actual profile cutoff/source/trace K derivatives')
    for name in ('first_error', 'second_error'):
        if levels[0][name] > 1e-7:
            audit.require(levels[-1][name] < .45*levels[0][name], 'K derivative refinement did not converge')
    r = real(payload['heights'][0])
    ecK = mp.diff(lambda K: profile_reference(K, payload)['profiles'][0][1], K0)
    numerator_derivative = mp.diff(lambda K: profile_reference(K, payload)['numerator'], K0)
    cutoff_omission = ecK-numerator_derivative/independent_base['cutoff']
    audit.require(r < independent_base['cutoff'] and abs(cutoff_omission) > mp.mpf('1e-6'),
                  'below-cutoff K sensitivity omission control is insensitive')
    trace = independent_base['profiles'][-1]
    chiK = mp.diff(lambda K: profile_reference(K, payload)['transition'][0], K0)
    trace_omission = (trace[1]-trace[0])*chiK
    audit.require(abs(trace_omission) > mp.mpf('1e-6'), 'dependent epsilon trace omitted-chi control is insensitive')
    return dict(points=len(points), heights=payload['heights'], K_derivative_refinement=levels,
                analytic_K_first=analytic[0], analytic_K_second=analytic[1],
                omitted_cutoff_derivative=float(cutoff_omission), omitted_chi_trace_derivative=float(trace_omission),
                scope='local production profiles; native coupled trace/state Jacobian tested separately')


def validate_rejections(data, audit):
    audit.context = 'invalid-input coverage'
    expected = {f'unordered_{j}' for j in (0, 2, 4, 6)}
    expected.update(f'nonpositive_control_{j}' for j in range(2, 15))
    expected.update(f'nonpositive_physical_{j}' for j in range(33, 40))
    expected.update(('negative_reynolds', 'nonfinite_sensor', 'nonfinite_interval', 'nonunit_normal',
                     'transition_overlaps_hill', 'no_flat_gaussian', 'no_flat_cylinder',
                     'affine_normal_covector', 'affine_zero_covector', 'active_undefined_profile'))
    names = {case['name'] for case in data['rejections']}
    audit.require(names == expected and len(data['rejections']) == len(expected), 'invalid-input coverage')
    audit.require(data['rejected'] == len(expected) and all(case['message'] for case in data['rejections']),
                  'invalid inputs were not actually rejected')
    return dict(count=len(expected), names=sorted(expected))


def validate(data):
    audit = Audit()
    audit.require(data['schema'] == 1 and data['success'] is True and data['checks'] > 0,
                  'native transition probe failed')
    sections = dict(transition=validate_transitions(data, audit), physical=validate_physical(data, audit),
                    coordinates=validate_coordinates(data, audit), blends=validate_blends(data, audit),
                    profiles=validate_profiles(data['profiles'], audit), rejections=validate_rejections(data, audit))
    return audit, sections


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--probe', type=Path)
    inputs.add_argument('--report', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    output = args.output/'validation.json'
    receipt = dict(schema=1, success=False, verification_status='failed',
                   started_utc=datetime.now(timezone.utc).isoformat(), commands=[], failures=[])
    try:
        if args.probe:
            probe = args.probe.resolve()
            report = (args.output/'native-transition.json').resolve()
            log = (args.output/'native-transition.log').resolve()
            env = os.environ.copy()
            env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                       BLIS_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
            command = [str(probe), '--output', str(report)]
            started = time.monotonic()
            with log.open('w') as stream:
                result = subprocess.run(command, cwd=args.output, env=env, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=180, check=False)
            receipt['commands'].append(dict(argv=command, cwd=str(args.output.resolve()), exit_code=result.returncode,
                                             ranks=1, threads=1, elapsed_seconds=time.monotonic()-started,
                                             log=str(log), log_sha256=hashlib.sha256(log.read_bytes()).hexdigest()))
            if result.returncode:
                raise RuntimeError(f'native transition probe failed: exit {result.returncode}; see {log}')
        else:
            report = args.report.resolve()
        receipt['native_report'] = str(report)
        receipt['native_report_sha256'] = hashlib.sha256(report.read_bytes()).hexdigest()
        data = json.loads(report.read_text())
        audit, sections = validate(data)
        receipt.update(success=True, verification_status='passed', checks=audit.checks,
                       native_checks=data['checks'], native_rejections=data['rejected'],
                       errors=audit.errors, sections=sections)
    except Exception as error:
        receipt['failures'].append(str(error))
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    output.write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps(dict(success=receipt['success'], checks=receipt.get('checks', 0),
                          failures=receipt['failures'], report=str(output))))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
