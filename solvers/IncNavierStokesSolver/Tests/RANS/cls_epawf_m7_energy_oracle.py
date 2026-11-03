"""Independent resolved-plus-k accounting for stored native M7 reports.

Flat Fixed CLS cases have no prescribed source/body force. Independent Jacobi
rows and physical mass matrices, dimensional AWF primitives and the indexed
CLS tensor construct the direct answer. Native weak loads are used only for
their explicitly labelled contraction/reaction audit, never to construct the
independent volume components. Upstream cases receive a limited mass/work
audit when manufactured forcing is present. Unforced neutral Upstream states
also receive independent physical sensor/chi and direct-volume accounting.
"""
import itertools
import math
import sys

sys.dont_write_bytecode = True
import numpy as np
from validate_cls_epawf_m1 import indexed
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m4 import spaces as flat_spaces
from validate_cls_epawf_m5 import spaces as mapped_spaces, physical_rows
from validate_cls_epawf_m5_geometry import geometry


class Audit:
    def __init__(self):
        self.checks, self.errors = 0, {}

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError(message)

    def near(self, actual, expected, name, tolerance=1e-10):
        actual, expected = np.asarray(actual), np.asarray(expected)
        self.require(actual.shape == expected.shape and np.all(np.isfinite(actual)) and
                     np.all(np.isfinite(expected)), 'Invalid M7 energy data: '+name)
        error = float(np.max(np.abs(actual-expected), initial=0.))/max(1., float(np.max(np.abs(expected), initial=0.)))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: {error} > {tolerance}')


def _fields(raw):
    return mapped_spaces(raw) if raw.get('mapped', False) else flat_spaces(raw)


def _sample(fields, state, offsets, block, points, geometry_id=None, side=0):
    return physical_rows(fields[block], state, offsets[block], np.atleast_2d(points), geometry_id, side)


def _velocity(fields, state, offsets, dimension, points, geometry_id=None, side=0):
    points = np.atleast_2d(points)
    U, G = np.zeros((len(points), 3)), np.zeros((len(points), 3, 3))
    for block, axis in enumerate((0, 2) if dimension == 2 else (0, 1, 2)):
        _, _, _, value, gradient = _sample(fields, state, offsets, block, points, geometry_id, side)
        U[:, axis], G[:, axis] = value, gradient
    return U, G


def _work(state, loads, fields, offsets, dimension, selected=None):
    result = 0.
    for block in range(dimension+1):
        ids = np.arange(fields[block].size)+offsets[block]
        test = fields[block].constant() if block == dimension else state[ids]
        if selected is not None:
            mask = np.array([int(index) in selected for index in ids])
            ids, test = ids[mask], test[mask]
        result += float(test@loads[ids])
    return result


def _mass(raw, fields, state, history, alpha):
    d, offsets = raw['dimension'], raw['offsets']
    storage, stored = 0., 0.
    for block in range(d+1):
        ids = np.arange(fields[block].size)+offsets[block]
        current, old, mass = state[ids], history[ids], fields[block].mass()
        test = fields[block].constant() if block == d else current
        storage += alpha*float(test@mass@(current-old))
        stored += (1. if block == d else .5)*float(test@mass@current)
    return storage, stored


def _columns(fields, dimension, H):
    values = {tuple(element['lower'][:2]+element['upper'][:2]) for element in fields[dimension+2].elements}
    if not values:
        raise ValueError('Missing native inner geometry')
    return sorted(values)


def _horizontal(fields, dimension, H, order):
    z, w = np.polynomial.legendre.leggauss(order)
    for x0, y0, x1, y1 in _columns(fields, dimension, H):
        for indices in itertools.product(range(order), repeat=dimension-1):
            point = np.array([x0+(z[indices[0]]+1)*(x1-x0)/2, 0., H])
            weight = w[indices[0]]*(x1-x0)/2
            if dimension == 3:
                point[1] = y0+(z[indices[1]]+1)*(y1-y0)/2
                weight *= w[indices[1]]*(y1-y0)/2
            yield point, weight


def _profile(raw, fields, state, chart):
    nu, H, rough, rK, rm, _, qstar, kappa = raw['parameters']
    metric = geometry(raw['map_parameters'], chart)
    sample = np.array([chart[0], chart[1], rK/metric['scale']])
    K = _sample(fields, state, raw['offsets'], raw['dimension'], [sample])[3][0]
    chi=raw['chi']
    if raw.get('upstream_transition'):
        # Independently reconstruct the seven physical sensor groups and
        # analytic gates. This path is used only with proven zero forcing.
        from cls_epawf_m6_operator_oracle import _diagnostic
        view=dict(raw,transition=dict(controls=raw['transition_controls'],weights=raw['transition_weights'],
                                     reference=raw['transition_reference']))
        chi=float(_diagnostic(view,state,fields,chart,zero_body=True)['gates'][0])
    return Profiles(dict(geometry=[nu, metric['height'], rough, rK, rm],
                         metric=[metric['kappa1'], metric['kappa2']], k=K,
                         q=qstar, kappa=kappa, chi=chi)), metric


def _forcing(raw, fields, state, point):
    d, offsets = raw['dimension'], raw['offsets']
    U, G = _velocity(fields, state, offsets, d, [point], side=1)
    U, G = U[0], G[0]
    metric = geometry(raw['map_parameters'], point)
    n, projection = metric['normal'], np.eye(3)-np.outer(metric['normal'], metric['normal'])
    ut, normal = projection@U, n@U
    grad_ut = projection@G-np.outer(n, U@metric['Br'])-normal*metric['Br']
    gp = _sample(fields, state, offsets, d+3, [point], side=1)[4][0]
    gk = _sample(fields, state, offsets, d, [point], side=1)[4][0]
    return projection@(grad_ut@ut+ut*np.trace(grad_ut)+gp-2*gk/3)


def _direct(raw, fields, state, audit):
    d, offsets = raw['dimension'], raw['offsets']
    nu, H, rough, rK, rm, _, qstar, kappa = raw['parameters']
    del rK, rm, qstar, kappa
    tangential_order = int(raw.get('quadrature_orders', [8, 16])[0])
    zn, wn = np.polynomial.legendre.leggauss(32)
    stress_source, advection, forcing = 0., 0., 0.
    inner_volume = outer_volume = 0.
    original_sink = corrected_sink = redistribution = 0.
    for point, horizontal in _horizontal(fields, d, H, tangential_order):
        profile, metric = _profile(raw, fields, state, point)
        height, scale = metric['height'], metric['scale']
        FA = _forcing(raw, fields, state, point)
        breaks = {0., float(height), rough, float(profile.zv), float(profile.ze), float(profile.rm), float(profile.rc)}
        start = min(float(profile.rc), nu/float(profile.a), float(profile.ze))
        while start < height:
            if start > 0:
                breaks.add(start)
            start *= 3
        for element in fields[d+2].elements:
            if element['lower'][0] <= point[0] <= element['upper'][0] and (
                    d == 2 or element['lower'][1] <= point[1] <= element['upper'][1]):
                breaks.update((element['lower'][2]*scale, element['upper'][2]*scale))
        breaks = sorted(value for value in breaks if 0 <= value <= height)
        for lo, hi in zip(breaks[:-1], breaks[1:]):
            r = lo+(zn+1)*(hi-lo)/2
            points = np.column_stack([np.full_like(r, point[0]), np.full_like(r, point[1]), r/scale])
            element = fields[d+2].element([point[0], point[1], (lo+hi)/(2*scale)])['geometry']
            weights = horizontal*wn*(hi-lo)/2*metric['baseArea']*np.asarray(profile.jacobian(r), float)
            U, G = _velocity(fields, state, offsets, d, points, element)
            K = _sample(fields, state, offsets, d, points, element)[3]
            pressure = _sample(fields, state, offsets, d+2, points, element)[3]
            S, divergence = G+G.transpose(0, 2, 1), np.trace(G, axis1=1, axis2=2)
            EA = float(profile.k)**1.5/(2.55*np.maximum(r, float(profile.ze)))
            p0 = np.asarray(profile.fields(r)[1], float)
            chi=float(profile.chi)
            direct = nu*np.sum(S*S, axis=(1, 2))/2+(2*K/3-pressure)*divergence+(1-chi)*EA+chi*p0
            stress_source += float(weights@direct)
            advection -= float(weights@np.einsum('qab,qa,qb->q', G, U, U))
            below = r <= rough
            forcing -= float(weights[below]@(U[below]@FA))
            inner_volume += float(np.sum(weights))
            original_sink += float(weights@EA)
            values = np.asarray(profile.fields(r), float)
            corrected_sink += float(weights@values[2])
            redistribution += float(weights@values[3])
    # Match the declared finite outer Gauss rule; CLS max(S,Omega) is not a
    # globally analytic function whose higher-order quadrature is identical.
    for element in fields[d+3].elements:
        points, weights = fields[d+3].quadrature(element, tangential_order)
        eid = element['geometry']
        U, G = _velocity(fields, state, offsets, d, points, eid)
        K = _sample(fields, state, offsets, d, points, eid)[3]
        epsilon = _sample(fields, state, offsets, d+1, points, eid)[3]
        pressure = _sample(fields, state, offsets, d+3, points, eid)[3]
        divergence = np.trace(G, axis1=1, axis2=2)
        closure = [indexed(g, k, e, 'physical') for g, k, e in zip(G, K, epsilon)]
        tensor = np.asarray([value['stress'].real for value in closure])
        production = -np.einsum('qab,qab->q', tensor, G)
        stress = -tensor+(2*K/3-pressure)[:, None, None]*np.eye(3)
        unreduced = np.einsum('qab,qab->q', G, stress)-(production-epsilon)
        direct = (2*K/3-pressure)*divergence+epsilon
        audit.near(unreduced, direct, 'indexed_outer_stress_production_cancellation', 3e-12)
        stress_source += float(weights@direct)
        advection -= float(weights@np.einsum('qab,qa,qb->q', G, U, U))
        outer_volume += float(np.sum(weights))
    return dict(stress_work_minus_k_source=stress_source, resolved_advection=advection,
                roughness_force_work=forcing, volumes=[inner_volume, outer_volume],
                original_dissipation_integral=original_sink, blended_dissipation_integral=corrected_sink,
                redistribution_integral=redistribution, quadrature=[tangential_order, 32, tangential_order])


def _interface(raw, fields, state):
    d, offsets, H = raw['dimension'], raw['offsets'], raw['parameters'][1]
    nu = raw['parameters'][0]
    result = np.zeros(8)
    for point, horizontal in _horizontal(fields, d, H, 8):
        profile, metric = _profile(raw, fields, state, point)
        area = metric['matching_area']
        measure = float(np.linalg.norm(area))
        normal, weight = area/measure, horizontal*measure
        Ui, Gi = _velocity(fields, state, offsets, d, [point], side=-1)
        Uo, Go = _velocity(fields, state, offsets, d, [point], side=1)
        k = _sample(fields, state, offsets, d, [point])[3][0]
        epsilon = _sample(fields, state, offsets, d+1, [point], side=1)[3][0]
        pi = _sample(fields, state, offsets, d+2, [point], side=-1)[3][0]
        po = _sample(fields, state, offsets, d+3, [point], side=1)[3][0]
        gki = _sample(fields, state, offsets, d, [point], side=-1)[4][0]
        gko = _sample(fields, state, offsets, d, [point], side=1)[4][0]
        nt = float(profile.a)*max(metric['height']-float(profile.zv), 0.)
        closure = indexed(Go[0], k, epsilon, 'physical')
        inner_stress = (nu+nt)*(Gi[0]+Gi[0].T)-pi*np.eye(3)
        outer_stress = -closure['stress'].real+(2*k/3-po)*np.eye(3)
        jump = (inner_stress-outer_stress)@normal
        normal_jump = float(normal@jump)
        k_flux_difference = float(((nu+nt)*gki-float(closure['nut'].real)*gko)@normal)
        result[0] += weight
        result[1] += weight*(pi-po)
        result[2] += weight*normal_jump*normal_jump
        result[3] += weight*float(jump@jump)
        result[4] = max(result[4], abs(normal_jump))
        result[5] += weight*k_flux_difference*k_flux_difference
        result[6] += weight*float(Ui[0]@normal)
        result[7] += weight*normal_jump
    for index in (2, 3, 5):
        result[index] = math.sqrt(max(0., result[index])/result[0])
    return result


def verify(raw):
    """Return a JSON-serializable independent final-stage energy receipt."""
    audit = Audit()
    audit.require(raw.get('diagnostics_schema') == 1 and not raw.get('failure'), 'M7 report is missing successful state diagnostics')
    state = np.asarray(raw['state'], float)
    history = np.asarray(raw.get('last_stage_history') or raw['state'], float)
    alpha = float(raw.get('last_stage_alpha', 0.))
    fields, d, offsets = _fields(raw), raw['dimension'], raw['offsets']
    audit.require(state.shape == history.shape and raw['eta'] in (0, 1), 'Invalid native history/anchor')
    storage, stored = _mass(raw, fields, state, history, alpha)
    energy = np.asarray(raw['energy_budget'])
    audit.near(energy[[1, 8]], [storage, stored], 'actual_physical_mass_storage_and_energy')
    loads = np.asarray(raw.get('final_stage_loads', []), float)
    constraints = raw.get('final_stage_constraints', [])
    work_origin = 'actual final unmasked stage loads'
    if not len(loads):
        audit.near(state, raw['fixed_reference'], 'zero_step_equal_reference_for_load_fallback', 1e-13)
        loads, constraints = np.asarray(raw['raw_reference_loads']), raw['reference_constraints']
        work_origin = 'zero-step state equals stored reference; original unmasked reference loads'
    audit.require(loads.shape == state.shape, 'Missing complete unmasked weak loads')
    unmasked = _work(state, loads, fields, offsets, d)
    selected = {int(row) for row, _ in constraints}
    reaction = _work(state, loads, fields, offsets, d, selected)
    anchor = raw['eta']*_work(state, np.asarray(raw['fixed_defect']), fields, offsets, d)
    audit.near(energy[0], unmasked, 'actual_unmasked_weak_energy_work')
    audit.near(raw['anchor_energy_work'], anchor, 'independent_fixed_anchor_work')
    last = raw['state_audits'][-1]['energy']
    audit.near(last[9:12], [reaction, anchor, unmasked-reaction-anchor], 'essential_reaction_and_anchor_accounting')
    result = dict(schema=1, checks=0, errors={}, scope='physical mass and actual unmasked-load/essential-reaction/anchor contracts',
                  dimension=d, mapped=raw.get('mapped', False), eta=raw['eta'], chi=raw['chi'],
                  stage_time=raw.get('last_stage_time'), alpha=alpha,
                  storage=storage, stored_resolved_plus_k=stored, unmasked_weak_work=unmasked,
                  essential_boundary_reaction_work=reaction, fixed_anchor_work=anchor,
                  free_row_balance=unmasked-reaction-anchor, unmasked_load_source=work_origin)
    unforced_transition=raw.get('upstream_transition') and raw.get('neutral_reference') and not any(
        '|'+prefix+'/' in raw.get('signature','') for prefix in (
            'BodyForceInner','BodyForceOuter','ScalarSourceInner','ScalarSourceOuter'))
    if raw.get('upstream_transition') and not unforced_transition:
        result.update(scope='limited Upstream audit: actual mapped mass/storage and unmasked work/reaction/anchor only',
                      direct_volume_components_reconstructed=False,
                      limitation='manufactured transition body/scalar forcing not reconstructed by this flat-neutral oracle')
    else:
        audit.require(int(raw['map_parameters'][0]) in (0, 1) and
                      raw['map_parameters'][3] == raw['map_parameters'][4] == raw['map_parameters'][5] == 0,
                      'Full energy oracle supports constant-H flat/plane neutral-data geometry')
        audit.require(abs(raw['parameters'][6]-3.1906750816242875) < 1e-10, 'Full energy oracle requires original CLS')
        direct = _direct(raw, fields, state, audit)
        if unforced_transition:
            from cls_epawf_m6_operator_oracle import _diagnostic
            view=dict(raw,transition=dict(controls=raw['transition_controls'],weights=raw['transition_weights'],
                                         reference=raw['transition_reference']))
            for row in raw['transition_samples']:
                independent=_diagnostic(view,state,fields,row['chart'],zero_body=True)
                audit.near(row['gates'],independent['gates'],'independent_unforced_transition_gates')
                audit.near(row['input'],independent['input'],'independent_unforced_transition_inputs')
            # The initial mixed-chi state exercises the active correction as
            # well as its baseline, before the physical transient changes it.
            initial=np.asarray(raw['accepted_states'][0],float)
            initial_direct=_direct(raw,fields,initial,audit)
            initial_energy=np.asarray(raw['state_audits'][0]['energy'])
            initial_storage,initial_stored=_mass(raw,fields,initial,initial,0.)
            audit.near(initial_energy[[1,8]],[initial_storage,initial_stored],'initial_unforced_physical_mass')
            terms=[initial_direct['stress_work_minus_k_source'],initial_direct['resolved_advection'],initial_direct['roughness_force_work']]
            audit.near(initial_energy[[2,3,4]],terms,'initial_mixed_chi_direct_energy_components')
            audit.near(initial_energy[0],sum(terms),'initial_mixed_chi_direct_equals_captured_weak_work')
            audit.near(raw['state_audits'][0]['wall_sources'][:3],
                [initial_direct['original_dissipation_integral'],initial_direct['blended_dissipation_integral'],initial_direct['redistribution_integral']],
                'initial_mixed_chi_wall_source_budgets')
            result['initial_unforced_direct']=initial_direct
        audit.near(energy[[2, 3, 4]], [direct['stress_work_minus_k_source'], direct['resolved_advection'],
                                      direct['roughness_force_work']], 'independent_direct_energy_components')
        total = storage+direct['stress_work_minus_k_source']+direct['resolved_advection']+direct['roughness_force_work']
        audit.near(unmasked, total, 'independent_direct_energy_equals_actual_weak_work')
        audit.near(energy[[6, 7]], [total, unmasked-total], 'independent_complete_energy_identity')
        audit.near(direct['redistribution_integral'], 0., 'independent_zero_net_redistribution')
        interface = _interface(raw, fields, state)
        audit.near(raw['state_audits'][-1]['interface'], interface, 'independent_total_matching_traction_and_k_flux')
        result.update(scope=('full unforced Upstream original-CLS direct energy, independent sensors/chi and matching-traction audit'
                             if unforced_transition else 'full Fixed original-CLS neutral-data direct energy and matching-traction audit'),
                      direct_volume_components_reconstructed=True, direct=direct,
                      independently_reconstructed_total=total, independent_interface=interface.tolist(),
                      interpretation='stage storage uses the actual BDF/IMEX history contract; not a new finite-difference definition of total-energy change')
    result.update(checks=audit.checks, errors=audit.errors, success=True)
    return result
