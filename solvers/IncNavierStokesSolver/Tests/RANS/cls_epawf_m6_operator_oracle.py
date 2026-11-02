"""Independent M6 physical sensors, constrained trace and complete weak action.

Only native mesh/basis maps and coefficient vectors define the discrete space.
Independent Jacobi rows, analytic normal geometry, dimensional M2 profile
primitives and the indexed M5 weak oracle supply the physical answer. Analytic
chi is evaluated in physical space and is never projected as a polynomial.
"""
from functools import lru_cache
import math
import sys

sys.dont_write_bytecode = True
import numpy as np
from scipy.integrate import quad
from scipy.special import eval_legendre, roots_jacobi

import validate_cls_epawf_m5 as weak
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m3 import basis
from validate_cls_epawf_m5_geometry import geometry
from cls_epawf_m6_cases import arc_offset


def _gate(q, lo, hi):
    if q <= lo:
        return 1.
    if q >= hi:
        return 0.
    s = (q-lo)/(hi-lo)
    # Different polynomial evaluation from the native factored expression.
    return 1-10*s**3+15*s**4-6*s**5


@lru_cache(maxsize=64)
def _arc_excess(amplitude, length):
    return arc_offset(amplitude, length)


def coordinate(raw, chart):
    """Signed physical arc, independent of rho and instantaneous velocity."""
    m, x = raw['map_parameters'], float(chart[0])
    kind, sx = int(m[0]), float(m[8])
    if kind in (0, 1):
        return x*math.hypot(1., sx)
    if kind != 4:
        raise ValueError('M6 upstream reference requires a declared exactly flat fetch')
    amplitude, length = float(m[11]), float(m[12])
    excess = _arc_excess(amplitude, length)
    if x <= -length:
        return x-excess
    if x >= length:
        return x+excess
    wave = math.pi/length
    integral, error = quad(lambda q: math.hypot(1., .5*amplitude*wave*math.sin(wave*q)),
                           0., x, epsabs=1e-13, epsrel=1e-13)
    if error > 1e-10:
        raise RuntimeError('Independent physical wall arc was not resolved')
    return integral


def _budget(raw, K, metric):
    """J-weighted dimensional AWF budget; does not solve a cutoff problem."""
    nu, _, rough, rK, rm, _, qstar, kappa = map(np.longdouble, raw['parameters'])
    del rK
    K, H = np.longdouble(K), np.longdouble(metric['height'])
    k1, k2 = np.longdouble(metric['kappa1']), np.longdouble(metric['kappa2'])
    if not (K > 0 and nu > 0 and 0 < rm <= H):
        raise ValueError('Invalid independent AWF baseline domain')
    root = np.sqrt(K)
    if rough:
        x = rough*root/(70*nu)
        exponent = max(np.longdouble('.5')-np.longdouble('.4')*x**np.longdouble('.7'),
                       1-np.longdouble('.79')*x**np.longdouble('-.28'))
        intercept = np.longdouble('10.7')*(1-x**exponent)
    else:
        intercept = np.longdouble('10.7')
    zv, slope = nu*intercept/root, np.longdouble('.2295')*root
    ze, amplitude = np.longdouble('5.1')*nu/root, (K/qstar)**np.longdouble('1.5')/kappa
    s1, s2 = k1+k2, k1*k2

    def volume(r):
        return r+s1*r*r/2+s2*r*r*r/3

    def cap(numerator, cutoff):
        below, above = min(H, cutoff), max(H, cutoff)
        return numerator*(volume(below)/cutoff+np.log(above/cutoff)+
                          s1*(above-cutoff)+s2*(above*above-cutoff*cutoff)/2)

    lower, intercept_D = max(np.longdouble(0), zv), nu-slope*zv
    production = np.longdouble(0)
    if H > lower:
        w0, w1 = intercept_D+slope*lower, intercept_D+slope*H
        delta = slope*(H-lower)
        c0 = 1-s1*intercept_D/slope+s2*(intercept_D/slope)**2
        c1 = s1/slope-2*s2*intercept_D/slope**2
        c2 = s2/slope**2
        production = (K/qstar)**2/slope*((c0-nu*c1)*np.log1p(delta/w0)+
            (c1-nu*c2)*delta+c2*delta*(w1+w0)/2-nu*c0*delta/(w1*w0))
    minimum = cap(amplitude, rm)
    original_integral = cap(K**np.longdouble('1.5')/np.longdouble('2.55'), ze)
    return dict(K=float(K), H=float(H), nu=float(nu), rough=float(rough), rm=float(rm),
                qstar=float(qstar), kappa=float(kappa), zv=float(zv), slope=float(slope),
                ze=float(ze), amplitude=float(amplitude), production=float(production), minimum=float(minimum),
                margin=float((production-minimum)/minimum), reynolds=float(H*root/nu),
                original_integral=float(original_integral), kappas=[float(k1), float(k2)])


def _profile(raw, law, chi):
    nu, _, rough, rK, rm, _, qstar, kappa = raw['parameters']
    return Profiles(dict(geometry=[nu, law['H'], rough, rK, rm], metric=law['kappas'],
                         k=law['K'], q=qstar, kappa=kappa, chi=chi))


def _sample(fields, state, block, offset, point, side=0):
    return fields[block].evaluate(np.asarray(state), offset, point, side)


def forcing_parts(raw, state, fields, chart, zero_body=False):
    """Separate physical P(grad pi-2 grad k/3-b) and P div(ut tensor ut)."""
    d, off = raw['dimension'], raw['offsets']
    p = np.asarray(chart, float).copy()
    p[2] = raw['parameters'][1]
    metric = geometry(raw['map_parameters'], p)
    n, Br = metric['normal'], metric['Br']
    projection = np.eye(3)-np.outer(n, n)
    velocity, gradient = np.zeros(3), np.zeros((3, 3))
    for block, axis in enumerate((0, 2) if d == 2 else (0, 1, 2)):
        value, derivative = _sample(fields, state, 0, off[block], p, 1)
        velocity[axis], gradient[axis] = value, derivative
    _, pressure = _sample(fields, state, d+3, off[d+3], p, 1)
    _, grad_k = _sample(fields, state, d, off[d], p, 1)
    ut, normal_velocity = projection@velocity, n@velocity
    grad_ut = projection@gradient-np.outer(n, velocity@Br)-normal_velocity*Br
    body = np.zeros(3) if zero_body else weak.bodies([p], d, False)[0]
    pressure_body = projection@(pressure-2*grad_k/3-body)
    convection = projection@(grad_ut@ut+ut*np.trace(grad_ut))
    return np.array([pressure_body, convection])


def _diagnostic(raw, state, fields, chart, zero_body=False):
    d, off = raw['dimension'], raw['offsets']
    start, end, ra, rb, baseline_start, on, off_state, moff, mon, roff, ron = raw['transition']['controls']
    del baseline_start
    p = np.asarray(chart, float).copy()
    p[2] = raw['parameters'][1]
    metric = geometry(raw['map_parameters'], p)
    sample = np.array([p[0], p[1], raw['parameters'][3]/metric['scale']])
    K = _sample(fields, state, d, off[d], sample)[0]
    law = _budget(raw, K, metric)
    ell = coordinate(raw, p)
    sensors = np.zeros(15)
    if ell < end:
        reference = np.asarray(raw['transition']['reference'])
        n, H, nu, qstar = metric['normal'], metric['height'], raw['parameters'][0], raw['parameters'][6]
        projection, scale2 = np.eye(3)-np.outer(n, n), K/qstar+(nu/H)**2
        pa, pb = p.copy(), p.copy()
        pa[2], pb[2] = ra/metric['scale'], rb/metric['scale']
        U, U0, Ua, Ua0, Ub, Ub0 = [np.zeros(3) for _ in range(6)]
        for block, axis in enumerate((0, 2) if d == 2 else (0, 1, 2)):
            U[axis] = _sample(fields, state, 0, off[block], p, 1)[0]
            U0[axis] = _sample(fields, reference, 0, off[block], p, 1)[0]
            Ua[axis] = _sample(fields, state, 0, off[block], pa)[0]
            Ua0[axis] = _sample(fields, reference, 0, off[block], pa)[0]
            Ub[axis] = _sample(fields, state, 0, off[block], pb)[0]
            Ub0[axis] = _sample(fields, reference, 0, off[block], pb)[0]
        ka = _sample(fields, state, d, off[d], pa)[0]
        ea = _sample(fields, state, d+1, off[d+1], pa)[0]
        ea0 = _sample(fields, reference, d+1, off[d+1], pa)[0]
        if not (ea > 0 and ea0 > 0 and ka > 0):
            raise ValueError('Invalid independent outer sensor denominator/sample')
        force = forcing_parts(raw, state, fields, p, zero_body)
        reference_force = forcing_parts(raw, reference, fields, p, True)
        scale = math.sqrt(scale2)
        sensors[:3] = projection@(U-U0)/scale
        sensors[3:6] = projection@((Ub-Ub0)-(Ua-Ua0))/scale
        sensors[6:9] = [n@(Ua-Ua0)/scale, (ka-K)/K, (ea-ea0)/ea0]
        sensors[9:12] = H*(force[0]-reference_force[0])/scale2
        sensors[12:15] = H*(force[1]-reference_force[1])/scale2
    groups = [sensors[:3], sensors[3:6], sensors[6:7], sensors[7:8], sensors[8:9],
              sensors[9:12], sensors[12:15]]
    departure = sum(weight*float(values@values) for weight, values in zip(raw['transition']['weights'], groups))
    spatial, state_gate = _gate(ell, start, end), _gate(departure, on*on, off_state*off_state)
    margin, reynolds = 1-_gate(law['margin'], moff, mon), 1-_gate(law['reynolds'], roff, ron)
    chi = spatial*state_gate*margin*reynolds if raw.get('chi', 1) else 0.
    return dict(input=np.array([ell, law['margin'], law['reynolds']]), sensors=sensors,
                gates=np.array([chi, spatial, state_gate, margin*reynolds, margin, reynolds, departure]), law=law)


def weight(raw, state, fields, chart):
    """Public independent analytic chi callback for the full weak oracle."""
    # The inactive branch requires only baseline-domain validation and budgets;
    # it never constructs Profiles, hence never solves a forbidden cutoff.
    return float(_diagnostic(raw, state, fields, chart)['gates'][0])


@lru_cache(maxsize=32)
def _segment_projection(degree):
    """Independent constrained transform: GLL data, exact interior GL mass."""
    count = max(8, degree+3)
    middle, _ = roots_jacobi(count-2, 1, 1)
    nodes = np.r_[-1., middle, 1.]
    weights = 2/(count*(count-1)*eval_legendre(count-1, nodes)**2)
    B, _ = basis(degree, nodes)
    z, w = np.polynomial.legendre.leggauss(degree+2)
    P, _ = basis(degree, z)
    mass = P.T@(w[:, None]*P)
    transform = np.zeros((degree+1, count))
    transform[0, 0], transform[1, -1] = 1., 1.
    rhs = B[:, 2:].T*weights[None, :]-mass[2:, :2]@transform[:2]
    transform[2:] = np.linalg.solve(mass[2:, 2:], rhs)
    return nodes, transform


def _trace_projection(raw, state, fields, frozen_weight_state=None):
    d, off = raw['dimension'], raw['offsets']
    field, H = fields[d+1], raw['parameters'][1]
    baseline_start = raw['transition']['controls'][4]
    baseline, correction, baseline_ids = {}, {}, set()
    # Each input state has its own immutable cache. No id-based/global cache
    # can accidentally reuse a sensor after a Newton finite-difference column.
    samples = {}

    def evaluate(q):
        key = tuple(q)
        if key not in samples:
            value = _diagnostic(raw, state, fields, q)
            law, chi = value['law'], value['gates'][0]
            if frozen_weight_state is not None:
                chi = weight(raw, frozen_weight_state, fields, q)
            ea = law['K']**1.5/(2.55*max(law['H'], law['ze']))
            change = 0.
            if chi:
                profile = _profile(raw, law, chi)
                ec = law['amplitude']/max(law['H'], float(profile.rc))
                change = chi*(ec-ea)
            samples[key] = ea, change
        return samples[key]

    for element in field.elements:
        if abs(element['lower'][2]-H) > 1e-11:
            continue
        degree, modes = element['degree'], element['degree']+1
        nodes, transform = _segment_projection(degree)
        x = element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
        if d == 2:
            data = np.array([evaluate((float(xx), 0., H)) for xx in x]).T
            coefficients = [transform@data[part] for part in range(2)]
        else:
            y = element['lower'][1]+(nodes+1)*(element['upper'][1]-element['lower'][1])/2
            data = np.array([[evaluate((float(xx), float(yy), H)) for xx in x] for yy in y])
            coefficients = [(transform@data[:, :, part]@transform.T).reshape(-1) for part in range(2)]
        ids = np.asarray(element['ids'][:modes**(d-1)], int)
        signs = np.asarray(element['sign'][:modes**(d-1)])
        for table, values in zip((baseline, correction), coefficients):
            for index, value in zip(ids, signs*values):
                index, value = int(index), float(value)
                if index in table and abs(table[index]-value) > 2e-10*max(1., abs(value)):
                    raise RuntimeError('Independent trace transform has inconsistent shared edge/corner')
                table[index] = value
        if coordinate(raw, element['lower']) >= baseline_start-1e-11:
            baseline_ids.update(map(int, ids))
    for index in baseline_ids:
        correction[index] = 0.
    return [baseline, correction], baseline_ids


def trace_projection(raw, state, fields):
    """Return [baseline_coefficients, constrained_correction_coefficients]."""
    return _trace_projection(raw, state, fields)[0]


def _external_lift(raw, state, trace_ids):
    result, off = np.array(state, float), raw['offsets']
    epsilon = raw['dimension']+1
    for row, value in raw['constraints']:
        row = int(row)
        if not (off[epsilon] <= row < off[epsilon+1] and row-off[epsilon] in trace_ids):
            result[row] = value
    return result


def _residual(raw, state, fields, trace_ids, normal_order):
    """Complete coupled operator, including its nonlinear matching rows."""
    physical = _external_lift(raw, state, trace_ids)
    view = dict(raw, lifted_state=physical)
    loads = weak.residual_oracle(view, normal_order=normal_order, verify=False,
        weight_function=lambda chart, values, spaces: weight(raw, values, spaces, chart), return_loads=True)['residual']
    loads = np.asarray(loads)
    offset = raw['offsets'][raw['dimension']+1]
    for row, value in raw['constraints']:
        row = int(row)
        if not (offset <= row < offset+fields[raw['dimension']+1].size and row-offset in trace_ids):
            loads[row] = state[row]-value
    base, correction = trace_projection(raw, physical, fields)
    for index in trace_ids:
        loads[offset+index] = state[offset+index]-base[index]-correction[index]
    return loads


def _directional(raw, fields, near, require):
    if raw['dimension'] != 2:
        return dict(independent_full_weak_actions=False, reason='bounded independent action audit is 2D; full 3D weak rows audited separately')
    d, off, state = raw['dimension'], raw['offsets'], np.asarray(raw['lifted_state'])
    maps, _ = _trace_projection(raw, state, fields)
    trace_ids = set(maps[0])
    matching_rows = np.array(sorted(off[d+1]+index for index in trace_ids), int)
    results = []
    cases = {case['block']: case for case in raw['transition']['directions']}
    for block in (d, d+1, d+3):
        require(block in cases, 'missing independent coupled action block '+str(block))
        case, levels = cases[block], []
        direction = np.asarray(case['vector'])
        require(np.max(np.abs(direction), initial=0.) > 0, 'empty coupled action direction')
        native = {float(action['step']): np.asarray(action['value']) for action in case['actions']}
        for h in (2e-4, 1e-4):
            plus, minus = state+h*direction, state-h*direction
            a = _residual(raw, plus, fields, trace_ids, 24)
            b = _residual(raw, minus, fields, trace_ids, 24)
            action = (a-b)/(2*h)
            reference = min(native, key=lambda step: abs(step-h))
            require(abs(reference-h) < 1e-14, 'missing native action at independent step')
            near(action, native[reference], f'complete_transition_action_block_{block}_h_{h}', tolerance=1e-7)
            levels.append(dict(step=h, maximum=float(np.max(np.abs(action))),
                               native_error=float(np.max(np.abs(action-native[reference]))), action=action))
        native_fine = native[min(native)]
        near(levels[-1]['action'], native_fine, f'complete_transition_action_block_{block}_fine', tolerance=1e-7)
        drift = float(np.max(np.abs(levels[0]['action']-levels[1]['action'])))
        # Directly isolate the implicit trace rows; current independent epsilon
        # samples include these very native coefficients, rather than a lag.
        h = 1e-4
        plus, minus = _external_lift(raw, state+h*direction, trace_ids), _external_lift(raw, state-h*direction, trace_ids)
        tp, tm = trace_projection(raw, plus, fields), trace_projection(raw, minus, fields)
        frozen_plus = _trace_projection(raw, plus, fields, state)[0]
        frozen_minus = _trace_projection(raw, minus, fields, state)[0]
        derivative = np.array([(tp[0][i]+tp[1][i]-tm[0][i]-tm[1][i])/(2*h) for i in sorted(trace_ids)])
        frozen = np.array([(frozen_plus[0][i]+frozen_plus[1][i]-frozen_minus[0][i]-frozen_minus[1][i])/(2*h) for i in sorted(trace_ids)])
        expected_rows = direction[matching_rows]-derivative
        near(native_fine[matching_rows], expected_rows, 'implicit_trace_action_block_'+str(block), tolerance=1e-7)
        omission = float(np.max(np.abs(derivative-frozen), initial=0.))
        require(omission > 1e-12, 'coupled trace chi-derivative omission control is insensitive in block '+str(block))
        results.append(dict(block=block, independent_steps=[level['step'] for level in levels],
                            same_step_native_errors=[level['native_error'] for level in levels],
                            independent_step_drift=drift,
                            finest_native_step=float(min(native)),
                            independent_to_finest_error=float(np.max(np.abs(levels[-1]['action']-native_fine))),
                            trace_omitted_chi_derivative_maximum=omission, complete_rows=len(state),
                            nonlinear_trace_rows=len(trace_ids)))
    return dict(independent_full_weak_actions=True, blocks=results,
                interpretation='central finite differences of independently reconstructed complete operators; reported finite-step drift')


def verify(raw, near, require, directional=True):
    """Qualify an actual extended RANSCurvedWallProbe transition report."""
    require('transition' in raw and raw.get('eta') == 0, 'M6 needs an unanchored transition operator')
    fields, state = weak.spaces(raw), np.asarray(raw['lifted_state'])
    section, d = raw['transition'], raw['dimension']
    controls = section['controls']
    require(len(controls) == 11 and len(section['weights']) == 7 and min(section['weights']) > 0,
            'invalid physical sensor controls/positive weights')
    require(controls[0] < controls[1] < controls[4] and 0 < controls[5] < controls[6] and
            0 < controls[7] < controls[8] and 0 < controls[9] < controls[10], 'invalid gate/buffer thresholds')
    require(raw['parameters'][1] < controls[2] < controls[3] < raw['upper'][2], 'outer sample distances are not independent')
    if int(raw['map_parameters'][0]) == 4:
        foot = coordinate(raw, [-raw['map_parameters'][12], 0., 0.])
        require(controls[4] < foot, 'baseline projection buffer reaches the hill')
    active, inactive, reference_departure = 0, 0, 0.
    profile_records = []
    for index, sample in enumerate(section['samples']):
        value = _diagnostic(raw, state, fields, sample['chart'])
        law, chi = value['law'], value['gates'][0]
        near(sample['input'], value['input'], 'transition_physical_inputs')
        sensors = np.concatenate(sample['sensors'])
        near(sensors, value['sensors'], 'independent_seven_transition_sensors')
        near(sample['gates'], value['gates'], 'independent_transition_gates')
        near(sample['forcing_parts'], forcing_parts(raw, state, fields, sample['chart']), 'independent_separate_forcing')
        reference = _diagnostic(raw, np.asarray(section['reference']), fields, sample['chart'], True)
        reference_departure = max(reference_departure, reference['gates'][-1])
        near(sample['reference_d2'], reference['gates'][-1], 'actual_projected_reference_departure', tolerance=1e-20)
        near(sample['reference_state'], reference['gates'][2], 'actual_projected_reference_state')
        cutoff = 0.
        if chi:
            active += 1
            profile = _profile(raw, law, chi)
            cutoff = float(profile.rc)
            require(sample['profile'][2] > 0, 'active corrected cutoff was not solved')
        else:
            inactive += 1
            require(sample['profile'][2] == 0, 'inactive branch solved a corrected cutoff')
        ea0 = law['K']**1.5/(2.55*law['ze'])
        eaH = law['K']**1.5/(2.55*max(law['H'], law['ze']))
        E0 = ea0 if not chi else (1-chi)*ea0+chi*law['amplitude']/cutoff
        EH = eaH if not chi else (1-chi)*eaH+chi*law['amplitude']/max(law['H'], cutoff)
        near(np.asarray(sample['profile'])[[0, 1, 3, 4, 5]], [law['K'], cutoff, chi, E0, EH], 'actual_state_profile')
        expected_change = chi*(law['original_integral']-law['production'])
        near(sample['budgets'], [0., expected_change, expected_change, 0., 0.], 'transition_full_and_redistribution_budgets', tolerance=3e-10)
        profile_records.append(dict(sample=index, coordinate=value['input'][0], chi=chi, cutoff=cutoff,
                                    d2=value['gates'][-1], full_k_budget_change=expected_change))
    require(active > 0 and inactive > 0, 'operator did not exercise both active and one-way inactive regions')
    require(reference_departure < 1e-20, 'projected reference sensors did not vanish independently')
    low = np.asarray(section['inactive_low_K'])
    q = [(raw['lower'][0]-3.)/2, .37 if d == 3 else 0., raw['parameters'][3]]
    low_budget = _budget(raw, 1e-10, geometry(raw['map_parameters'], q))
    near(low, [1e-10, 0., 0., 0., low_budget['margin'], low_budget['reynolds']], 'inadmissible_correction_valid_baseline', tolerance=2e-10)
    maps, baseline_ids = _trace_projection(raw, state, fields)
    for part, expected in enumerate(maps):
        actual = {int(index): value for index, value in section['trace'][part]}
        require(set(actual) == set(expected), 'constrained trace ID/map mismatch')
        ids = sorted(expected)
        near([actual[index] for index in ids], [expected[index] for index in ids], 'independent_trace_projection_part_'+str(part))
    epsilon_offset = raw['offsets'][d+1]
    ids = sorted(maps[0])
    near(state[epsilon_offset+np.array(ids)], [maps[0][i]+maps[1][i] for i in ids], 'implicit_trace_lift_fixed_point')
    require(bool(baseline_ids), 'missing designated baseline support')
    near([maps[1][i] for i in sorted(baseline_ids)], np.zeros(len(baseline_ids)), 'designated_baseline_trace_support', tolerance=1e-15)
    previous_checks = weak.checks
    complete = weak.residual_oracle(raw, normal_order=24,
        weight_function=lambda chart, values, spaces: weight(raw, values, spaces, chart))
    finite_rule = weak.residual_oracle(raw, normal_order=16,
        weight_function=lambda chart, values, spaces: weight(raw, values, spaces, chart))
    identity = section['equal_input_baseline']
    near(identity['fixed_operator'], identity['disabled_transition_operator'], 'equal_input_adapted_baseline_operator')
    near([identity['maximum'], identity['constraint_maximum']], [0., 0.], 'equal_input_adapted_baseline_constraints')
    baseline_raw = dict(raw, lifted_state=identity['state'], chi=0.)
    baseline_answer = weak.residual_oracle(baseline_raw, normal_order=24, verify=False, return_loads=True)['residual']
    near(identity['fixed_operator'], baseline_answer, 'independent_equal_input_baseline_all_rows', tolerance=3e-10)
    directional_receipt = _directional(raw, fields, near, require) if directional else dict(independent_full_weak_actions=False, reason='disabled by caller')
    return dict(schema=1, dimension=d, active_samples=active, inactive_samples=inactive,
                reference_departure_maximum=float(reference_departure), physical_samples=profile_records,
                constrained_trace=dict(rows=len(ids), baseline_rows=len(baseline_ids),
                    transform='independent GLL sampled constrained Modified_A transform; exact interior GL mass; commuting tensor transforms',
                    analytic_chi_projected=False, dependent_trace_is_unknown=True),
                complete_weak=complete, independent_inner_quadrature_control=finite_rule,
                additional_weak_checks=weak.checks-previous_checks,
                equal_input_baseline=dict(rows=len(identity['fixed_operator']), native_operator_maximum=identity['maximum'],
                    native_constraint_maximum=identity['constraint_maximum'], independently_reconstructed=True),
                directional=directional_receipt)
