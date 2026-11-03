"""Independent source, stress and energy audits for M8 native reports.

Uses Jacobi coefficient reconstruction, dimensional AWF primitives and an
indexed CLS tensor.  It never calls the production robust kernel.  The source
normalization is reconstructed on the declared finite source Gauss rule;
this verifies that discrete operator, not an everywhere-positivity theorem.
"""
import math
import sys

sys.dont_write_bytecode = True
import numpy as np

from validate_cls_epawf_m1 import indexed
from cls_epawf_m7_energy_oracle import (
    Audit, _fields, _forcing, _horizontal, _mass, _profile, _sample,
    _velocity, _work,
)
from validate_cls_epawf_m5_geometry import geometry


def selected_stress(raw, gradient, k, epsilon, chart):
    original = indexed(gradient, k, epsilon, 'physical')
    stress = original['stress'].real.copy()
    if raw['model_variant'] == 'Original':
        return stress, float(original['nut'].real)
    beta, gamma, _, _, _, scale = raw['robust_parameters'][:6]
    distance = float(geometry(raw['map_parameters'], chart)['distance'])
    flow = np.asarray(raw['robust_reference_flow'])
    normal = np.asarray(raw['robust_reference_normal'])
    us, kappa = raw['ustar'], raw['parameters'][7]
    reference = us/(kappa*distance)*(np.outer(flow, normal)+np.outer(normal, flow))
    stress += (beta-1)*original['cubic'].real
    stress -= (gamma+scale*kappa*us*distance)*(gradient+gradient.T-reference)
    return stress, float(original['nut'].real)


def normal_rule(profile, raw, fields, metric):
    """Independent Legendre roots/weights on declared physical source panels."""
    height, scale = metric['height'], metric['scale']
    cutoff = float(profile.rc) if profile.chi > 0 else 0.
    breaks = {0., float(height)}
    for value in (profile.h, profile.zv, profile.ze, profile.rm, cutoff):
        if 0 < value < height:
            breaks.add(float(value))
    for element in fields[raw['dimension']+2].elements:
        for value in (element['lower'][2]*scale, element['upper'][2]*scale):
            if 0 < value < height:
                breaks.add(float(value))
    start = min(float(profile.rm), float(profile.nu/profile.a))
    if cutoff:
        start = min(start, cutoff)
    for _ in range(4096):
        if start >= height:
            break
        if start > 0:
            breaks.add(start)
        start *= 2
    z, w = np.polynomial.legendre.leggauss(int(raw['quadrature_orders'][1]))
    points, weights = [], []
    breaks = sorted(breaks)
    for lo, hi in zip(breaks, breaks[1:]):
        points.extend(lo+(z+1)*(hi-lo)/2)
        weights.extend(w*(hi-lo)/2)
    return np.asarray(points), np.asarray(weights)


def ray_source(raw, fields, state, chart):
    profile, metric = _profile(raw, fields, state, chart)
    r, dr = normal_rule(profile, raw, fields, metric)
    chart_points = np.column_stack([np.full_like(r, chart[0]), np.full_like(r, chart[1]), r/metric['scale']])
    local_k = np.array([fields[raw['dimension']].evaluate(state, raw['offsets'][raw['dimension']], point)[0]
                        for point in chart_points])
    values = np.asarray(profile.fields(r), float)
    pre_e, pre_d = values[2], values[2]-values[3]
    measure = dr*np.asarray(profile.jacobian(r), float)
    theta, delta, smooth = raw['robust_parameters'][2:5]
    capped = bool(raw['robust_parameters'][6])
    if capped:
        t = np.minimum(1., local_k/(theta*float(profile.k)))
        if np.any(t < 0):
            raise RuntimeError('Negative k is outside the repaired source domain')
        limiter = 10*t**3-15*t**4+6*t**5 if smooth else t
        budget = float(measure@pre_d)
        b = float(measure@(pre_d*limiter))/budget
        m = max(b, delta)
        lam = b/m
        e, destruction = lam*pre_e, pre_d*limiter/m
        rate = np.divide(destruction, local_k, out=np.zeros_like(local_k), where=local_k > 0)
        zero = local_k == 0
        rate[zero] = (0. if smooth else 1.)*pre_d[zero]/(m*theta*float(profile.k))
    else:
        limiter, e, destruction = np.ones_like(local_k), pre_e, pre_d
        lam, b, m = 1., 1., 1.
        rate = np.divide(destruction, local_k, out=np.zeros_like(local_k), where=local_k > 0)
    redistribution = e-destruction
    return dict(profile=profile, metric=metric, r=r, points=chart_points, measure=measure,
                local_k=local_k, pre_e=pre_e, pre_d=pre_d, e=e, destruction=destruction,
                redistribution=redistribution, limiter=limiter, b=b, m=m, lam=lam, rate=rate)


def direct(raw, fields, state, audit):
    d, offsets = raw['dimension'], raw['offsets']
    nu, H, rough = raw['parameters'][:3]
    stress_source = advection = rough_work = 0.
    source = np.zeros(8)
    source[6] = 1.
    wall = np.zeros(5)
    inner_volume = outer_volume = 0.
    inner_covariance_minimum = outer_covariance_minimum = math.inf
    modification = np.zeros(3)
    for chart, horizontal in _horizontal(fields, d, H, int(raw['quadrature_orders'][0])):
        ray = ray_source(raw, fields, state, chart)
        profile, metric = ray['profile'], ray['metric']
        force = _forcing(raw, fields, state, chart)
        weights = horizontal*metric['baseArea']*ray['measure']
        # Source panels include every inner element interface, so a pointwise
        # native element lookup is an independent coefficient reconstruction.
        U, G, pressure = [], [], []
        for point in ray['points']:
            u, g = _velocity(fields, state, offsets, d, [point])
            U.append(u[0]); G.append(g[0])
            pressure.append(_sample(fields, state, offsets, d+2, [point])[3][0])
        U, G, pressure = np.asarray(U), np.asarray(G), np.asarray(pressure)
        S = G+G.transpose(0, 2, 1)
        divergence = np.trace(G, axis1=1, axis2=2)
        nt = float(profile.a)*np.maximum(ray['r']-float(profile.zv), 0.)
        tensor = 2*ray['local_k'][:, None, None]/3*np.eye(3)-nt[:, None, None]*S
        production = -np.einsum('qab,qab->q', tensor, G)
        inner_covariance_minimum = min(inner_covariance_minimum, float(np.min(np.linalg.eigvalsh(tensor)[:, 0])))
        stress = (nu+nt)[:, None, None]*S-pressure[:, None, None]*np.eye(3)
        unreduced = np.einsum('qab,qab->q', G, stress)-(production-ray['destruction'])
        reduced = nu*np.sum(S*S, axis=(1, 2))/2+(2*ray['local_k']/3-pressure)*divergence+ray['destruction']
        audit.near(unreduced, reduced, 'repaired_inner_stress_production_cancellation', 3e-12)
        stress_source += float(weights@reduced)
        advection -= float(weights@np.einsum('qab,qa,qb->q', G, U, U))
        below = ray['r'] <= rough
        rough_work -= float(weights[below]@(U[below]@force))
        inner_volume += float(np.sum(weights))
        source[:4] += [float(weights@ray['pre_e']), float(weights@ray['e']),
                       float(weights@ray['destruction']), float(weights@ray['redistribution'])]
        if raw['robust_parameters'][6]:
            source[4] += float(np.sum(weights[ray['limiter'] < 1]))
            source[5] += float(np.sum(weights)) if ray['b'] < raw['robust_parameters'][3] else 0.
            source[6] = min(source[6], ray['lam'])
            source[7] = max(source[7], float(np.max(ray['rate'])))
        ea = float(profile.k)**1.5/(2.55*np.maximum(ray['r'], float(profile.ze)))
        wall += [float(weights@ea), float(weights@ray['e']), float(weights@ray['redistribution']),
                 float(weights@production), float(weights@(nu*np.sum(S*S, axis=(1, 2))/2))]
    for element in fields[d+3].elements:
        points, weights = fields[d+3].quadrature(element, int(raw['quadrature_orders'][0]))
        U, G = _velocity(fields, state, offsets, d, points, element['geometry'])
        K = _sample(fields, state, offsets, d, points, element['geometry'])[3]
        epsilon = _sample(fields, state, offsets, d+1, points, element['geometry'])[3]
        pressure = _sample(fields, state, offsets, d+3, points, element['geometry'])[3]
        tensor = np.array([selected_stress(raw, g, k, e, p)[0] for g, k, e, p in zip(G, K, epsilon, points)])
        original = np.array([indexed(g, k, e, 'physical')['stress'].real for g, k, e in zip(G, K, epsilon)])
        outer_covariance_minimum = min(outer_covariance_minimum, float(np.min(np.linalg.eigvalsh(tensor)[:, 0])))
        production = -np.einsum('qab,qab->q', tensor, G)
        stress = -tensor+(2*K/3-pressure)[:, None, None]*np.eye(3)
        reduced = (2*K/3-pressure)*np.trace(G, axis1=1, axis2=2)+epsilon
        unreduced = np.einsum('qab,qab->q', G, stress)-(production-epsilon)
        audit.near(unreduced, reduced, 'selected_outer_stress_production_cancellation', 3e-12)
        stress_source += float(weights@reduced)
        advection -= float(weights@np.einsum('qab,qa,qb->q', G, U, U))
        change = -np.einsum('qab,qab->q', tensor-original, G)
        modification[0] += float(weights@change)
        modification[1] += float(weights@(production+np.einsum('qab,qab->q', original, G)))
        outer_volume += float(np.sum(weights))
    modification[2] = modification[0]-modification[1]
    return dict(stress_work_minus_k_source=stress_source, resolved_advection=advection,
                roughness_force_work=rough_work, robust_source=source.tolist(),
                wall_sources=wall.tolist(), constitutive_modification=modification.tolist(),
                inner_covariance_minimum=inner_covariance_minimum,
                outer_covariance_minimum=outer_covariance_minimum,
                volumes=[inner_volume, outer_volume],
                quadrature='independently computed Legendre/Jacobi rows on declared finite native source panels')


def interface(raw, fields, state):
    d, offsets, H = raw['dimension'], raw['offsets'], raw['parameters'][1]
    nu = raw['parameters'][0]
    result = np.zeros(8)
    for point, horizontal in _horizontal(fields, d, H, 8):
        ray = ray_source(raw, fields, state, point)
        profile, metric = ray['profile'], ray['metric']
        area = metric['matching_area']; measure = float(np.linalg.norm(area))
        normal, weight = area/measure, horizontal*measure
        Ui, Gi = _velocity(fields, state, offsets, d, [point], side=-1)
        _, Go = _velocity(fields, state, offsets, d, [point], side=1)
        k = _sample(fields, state, offsets, d, [point])[3][0]
        epsilon = _sample(fields, state, offsets, d+1, [point], side=1)[3][0]
        pi = _sample(fields, state, offsets, d+2, [point], side=-1)[3][0]
        po = _sample(fields, state, offsets, d+3, [point], side=1)[3][0]
        gki = _sample(fields, state, offsets, d, [point], side=-1)[4][0]
        gko = _sample(fields, state, offsets, d, [point], side=1)[4][0]
        nt = float(profile.a)*max(metric['height']-float(profile.zv), 0.)
        tensor, outer_nut = selected_stress(raw, Go[0], k, epsilon, point)
        inner_stress = (nu+nt)*(Gi[0]+Gi[0].T)-pi*np.eye(3)
        outer_stress = -tensor+(2*k/3-po)*np.eye(3)
        jump = (inner_stress-outer_stress)@normal
        normal_jump = float(normal@jump)
        flux_jump = float(((nu+nt)*gki-outer_nut*gko)@normal)
        result[0] += weight; result[1] += weight*(pi-po)
        result[2] += weight*normal_jump**2; result[3] += weight*float(jump@jump)
        result[4] = max(result[4], abs(normal_jump)); result[5] += weight*flux_jump**2
        result[6] += weight*float(Ui[0]@normal); result[7] += weight*normal_jump
    for index in (2, 3, 5):
        result[index] = math.sqrt(max(0., result[index])/result[0])
    return result


def verify(raw, full_volume=True):
    audit = Audit()
    audit.require(raw.get('diagnostics_schema') == 1 and not raw.get('failure'), 'Missing successful M8 diagnostics')
    audit.require(raw.get('model_variant') in ('Original', 'CubicScaled', 'Regularized'), 'Missing explicit model identity')
    state = np.asarray(raw['state'], float)
    history = np.asarray(raw.get('last_stage_history') or raw['state'], float)
    alpha = float(raw.get('last_stage_alpha', 0.))
    fields, d, offsets = _fields(raw), raw['dimension'], raw['offsets']
    storage, stored = _mass(raw, fields, state, history, alpha)
    energy = np.asarray(raw['energy_budget'])
    audit.near(energy[[1, 8]], [storage, stored], 'm8_actual_physical_mass_storage_and_energy')
    loads = np.asarray(raw.get('final_stage_loads', []), float)
    constraints = raw.get('final_stage_constraints', [])
    if not len(loads):
        audit.near(state, raw['fixed_reference'], 'zero_step_reference_load_fallback', 1e-13)
        loads, constraints = np.asarray(raw['raw_reference_loads']), raw['reference_constraints']
    audit.require(loads.shape == state.shape, 'Missing complete unmasked M8 weak loads')
    unmasked = _work(state, loads, fields, offsets, d)
    reaction = _work(state, loads, fields, offsets, d, {int(row) for row, _ in constraints})
    anchor = raw['eta']*_work(state, np.asarray(raw['fixed_defect']), fields, offsets, d)
    audit.near(energy[0], unmasked, 'm8_actual_unmasked_weak_energy_work')
    audit.near(raw['anchor_energy_work'], anchor, 'm8_fixed_anchor_energy_work')
    audit.near(raw['state_audits'][-1]['energy'][9:12], [reaction, anchor, unmasked-reaction-anchor],
               'm8_essential_reaction_anchor_accounting')
    result = dict(schema=1, checks=0, errors={}, model_variant=raw['model_variant'],
                  physical_mass_storage=storage, stored_resolved_plus_k=stored,
                  unmasked_work=unmasked, reaction_work=reaction, anchor_work=anchor,
                  direct_volume_components_reconstructed=False,
                  scope='independent physical mass, complete captured weak loads and essential reactions')
    if full_volume:
        # Full direct account requires absence of prescribed body/scalar
        # functions.  Manufactured MMS records retain the limited scope.
        audit.require(not any('|'+prefix+'/' in raw.get('signature', '') for prefix in
            ('BodyForceInner', 'BodyForceOuter', 'ScalarSourceInner', 'ScalarSourceOuter')),
            'Full independent M8 energy oracle does not infer manufactured forcing')
        answer = direct(raw, fields, state, audit)
        native = raw['state_audits'][-1]
        audit.near(native['robust_source'], answer['robust_source'], 'independent_m8_repaired_ray_source_budget')
        audit.near(native['wall_sources'], answer['wall_sources'], 'independent_m8_wall_source_budget')
        audit.near(native['constitutive_modification'], answer['constitutive_modification'], 'independent_m8_constitutive_work')
        audit.near(native['source_volume_inner'][3], answer['inner_covariance_minimum'], 'independent_m8_inner_covariance')
        audit.near(native['source_volume_outer'][3], answer['outer_covariance_minimum'], 'independent_m8_outer_covariance')
        audit.near(energy[[2, 3, 4]], [answer['stress_work_minus_k_source'], answer['resolved_advection'],
                                      answer['roughness_force_work']], 'independent_m8_direct_energy_components')
        total = storage+answer['stress_work_minus_k_source']+answer['resolved_advection']+answer['roughness_force_work']
        audit.near(unmasked, total, 'independent_m8_complete_energy_identity')
        traction = interface(raw, fields, state)
        audit.near(native['interface'], traction, 'independent_m8_total_traction_k_flux')
        result.update(scope='independent selected indexed stress, physical source normalization, full direct energy and interface account',
                      direct_volume_components_reconstructed=True, direct=answer,
                      independently_reconstructed_total=total, independent_interface=traction.tolist())
    result.update(checks=audit.checks, errors=audit.errors, success=True)
    return result
