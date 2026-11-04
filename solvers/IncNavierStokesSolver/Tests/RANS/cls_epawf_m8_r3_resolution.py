#!/usr/bin/env python3
"""Read-only R3 residual-resolution and consistent-mass contact analysis.

Never invokes a native executable or advances a state.  Numerical invocation
belongs to Linux.  All models, equations, quadrature orders, saved coefficients
and the absolute2e-12 root criterion remain fixed.  Frozen-load accumulation
bounds and exact-polynomial mass bounds have deliberately separate scopes.
"""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

from cls_epawf_m8_r1_nonlinear import Audit
from cls_epawf_m8_r2_preflight import (Inventory, add, finite_ranges,
    trace_audit, weak_loads, native_domain_audit)
from cls_epawf_m8_r2_oracle import (selected_inner, selected_outer,
    indexed_outer, outer_evaluation_statistics)
from cls_epawf_m8_native_oracle import ray_source
from cls_epawf_m7_energy_oracle import _horizontal, _forcing
from cls_epawf_m6_operator_oracle import _segment_projection, _diagnostic
from validate_cls_epawf_m5 import physical_rows
from cls_epawf_m8_r3_mass import TensorMass, POLYNOMIALS, rational, decimal


TERMS = ('inner_momentum_stress', 'outer_momentum_stress', 'conservative_advection',
    'k_diffusion', 'epsilon_diffusion', 'inner_k_production', 'inner_k_destruction',
    'outer_k_production', 'outer_k_destruction', 'epsilon_production',
    'epsilon_destruction', 'epsilon_yap', 'prescribed_scalar_source',
    'roughness_force', 'prescribed_body', 'continuity', 'gauge')
TOLERANCE = 2e-10


class ResolutionAudit(Audit):
    """Same physical norm, with triangular solves of the cached Cholesky."""
    def norms(self, residual, constraints):
        residual = np.asarray(residual, float)
        values = []
        for block, field in enumerate(self.fields):
            free = np.array([i for i in range(field.size) if self.off[block]+i not in constraints], int)
            key = (block, tuple(free))
            if key not in self._factors:
                self._factors[key] = np.linalg.cholesky(self.mass[block][np.ix_(free, free)]) if len(free) else np.zeros((0, 0))
            load = residual[self.off[block]+free]
            inverse = solve_triangular(self._factors[key], load, lower=True, check_finite=False) if len(free) else np.zeros(0)
            rate = self.us/self.H if block >= self.d+2 else self.state_scale(block)*self.us/self.H
            values.append(float(np.linalg.norm(inverse))/(rate*math.sqrt(field.volume)))
        values.append(abs(float(residual[-1]))/self.us**2)
        bc = max((abs(float(residual[row]))/self.state_scale(next(b for b in range(len(self.fields))
            if self.off[b] <= row < self.off[b]+self.fields[b].size)) for row in constraints), default=0.)
        return np.asarray(values+[bc])


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def vector_sha(value):
    return hashlib.sha256(np.asarray(value, dtype='<f8').tobytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source),
        *arguments], text=True).strip()


def independent_terms(raw, fields, state):
    """Indexed full stress/P and independently generated Gauss/Jacobi rows."""
    d, offsets = raw['dimension'], raw['offsets']
    nu, H, rough = raw['parameters'][:3]
    result = {key: np.zeros(len(state)) for key in TERMS}
    for key in ('diagnostic_inner_Eplus', 'diagnostic_inner_minus_Cplus',
                'diagnostic_outer_cancelled_momentum',
                'diagnostic_outer_cancelled_k_production',
                'diagnostic_outer_cancelled_epsilon_production'):
        result[key] = np.zeros(len(state))
    meta = dict(source_rays=0, inner_volume=0., outer_volume=0.,
                minimum_b=1., minimum_lambda=1., maximum_source_identity_error=0.)

    def load(key, block, rows, weights, source=None, flux=None):
        ids, phi, derivative = rows[:3]
        zero = np.zeros(len(weights))
        add(result[key], offsets[block], ids, phi, derivative, weights,
            zero if source is None else source,
            np.zeros((len(weights), 3)) if flux is None else flux)

    def integrate(identity, points, weights, inner, ray=None, index=None, force=None):
        vr = [physical_rows(fields[0], state, offsets[b], points, identity) for b in range(d)]
        kr = physical_rows(fields[d], state, offsets[d], points, identity)
        pr = physical_rows(fields[d+2 if inner else d+3], state,
            offsets[d+2 if inner else d+3], points, identity)
        k, gk, pressure = kr[3], kr[4], pr[3]
        U, G = np.zeros((len(points), 3)), np.zeros((len(points), 3, 3))
        for block, axis in enumerate((0, 2)):
            U[:, axis], G[:, axis] = vr[block][3], vr[block][4]
        S = G+G.transpose(0, 2, 1)
        if inner:
            nt = float(ray['profile'].a)*np.maximum(points[:, 2]-float(ray['profile'].zv), 0.)
            selected = [selected_inner(raw, g, kval, viscosity) for g, kval, viscosity in zip(G, k, nt)]
            Fstress = np.asarray([row['turbulent'].real for row in selected])
            P = np.asarray([float(row['production'].real) for row in selected])
            stress = nu*S+Fstress-pressure[:, None, None]*np.eye(3)
            diffusion = nu+nt
            destruction, E, C = ray['destruction'][index], ray['e'][index], ray['redistribution'][index]
            meta['maximum_source_identity_error'] = max(meta['maximum_source_identity_error'],
                float(np.max(np.abs(E-C-destruction), initial=0.)))
            load('inner_k_production', d, kr, weights, -P)
            load('inner_k_destruction', d, kr, weights, destruction)
            load('diagnostic_inner_Eplus', d, kr, weights, E)
            load('diagnostic_inner_minus_Cplus', d, kr, weights, -C)
            body = (points[:, 2] <= rough)[:, None]*force
            meta['inner_volume'] += float(np.sum(weights))
        else:
            er = physical_rows(fields[d+1], state, offsets[d+1], points, identity)
            epsilon, ge = er[3], er[4]
            selected = [selected_outer(raw, g, kval, e, p)
                for g, kval, e, p in zip(G, k, epsilon, points)]
            R = np.asarray([row[0] for row in selected])
            nt = np.asarray([row[1] for row in selected])
            P = -np.einsum('qab,qab->q', R, G)
            stress = -R+(2*k/3-pressure)[:, None, None]*np.eye(3)
            diffusion = nt
            ratio = k**1.5/(2.5*epsilon*points[:, 2])
            yap = .83*epsilon**2/k*np.maximum(ratio-1, 0)*ratio**2
            load('outer_k_production', d, kr, weights, -P)
            load('outer_k_destruction', d, kr, weights, epsilon)
            load('epsilon_production', d+1, er, weights, -1.44*epsilon/k*P)
            load('epsilon_destruction', d+1, er, weights, 1.92*epsilon**2/k)
            load('epsilon_yap', d+1, er, weights, -yap)
            load('epsilon_diffusion', d+1, er, weights, flux=nt[:, None]/1.3*ge)
            load('conservative_advection', d+1, er, weights, flux=-U*epsilon[:, None])
            beta, uniform, _, _, _, scale = raw['robust_parameters'][:6]
            frame = np.outer(raw['robust_reference_flow'], raw['robust_reference_normal'])
            reference = (raw['ustar']/(raw['parameters'][7]*points[:, 2]))[:, None, None]*(frame+frame.T)
            basic = [indexed_outer(g, kval, e) for g, kval, e in zip(G, k, epsilon)]
            quadratic = np.asarray([row['quadratic'].real for row in basic])
            cubic = np.asarray([row['cubic'].real for row in basic])
            gamma = uniform+scale*raw['parameters'][7]*raw['ustar']*points[:, 2]
            B = nt[:, None, None]*S-quadratic-beta*cubic+gamma[:, None, None]*(S-reference)
            stable_stress = B-pressure[:, None, None]*np.eye(3)
            stable_P = np.einsum('qab,qab->q', B-2*k[:, None, None]/3*np.eye(3), G)
            for block, axis in enumerate((0, 2)):
                load('diagnostic_outer_cancelled_momentum', block, vr[block], weights, flux=stable_stress[:, axis])
            load('diagnostic_outer_cancelled_k_production', d, kr, weights, -stable_P)
            load('diagnostic_outer_cancelled_epsilon_production', d+1, er, weights, -1.44*epsilon/k*stable_P)
            body = np.zeros_like(U)
            meta['outer_volume'] += float(np.sum(weights))
        for block, axis in enumerate((0, 2)):
            load('inner_momentum_stress' if inner else 'outer_momentum_stress',
                block, vr[block], weights, flux=stress[:, axis])
            load('conservative_advection', block, vr[block], weights, flux=-U[:, axis, None]*U)
            load('roughness_force', block, vr[block], weights, -body[:, axis])
        load('k_diffusion', d, kr, weights, flux=diffusion[:, None]*gk)
        load('conservative_advection', d, kr, weights, flux=-U*k[:, None])
        load('continuity', d+2 if inner else d+3, pr, weights, np.trace(G, axis1=1, axis2=2))

    for chart, horizontal in _horizontal(fields, d, H, int(raw['quadrature_orders'][0])):
        ray = ray_source(raw, fields, state, chart)
        groups = defaultdict(list)
        for index, point in enumerate(ray['points']):
            groups[fields[0].element(point)['geometry']].append(index)
        force = _forcing(raw, fields, state, chart)
        for identity, index in groups.items():
            index = np.asarray(index, int)
            weights = horizontal*ray['metric']['baseArea']*ray['measure'][index]
            integrate(identity, ray['points'][index], weights, True, ray, index, force)
        meta['source_rays'] += 1
        meta['minimum_b'] = min(meta['minimum_b'], float(ray['b']))
        meta['minimum_lambda'] = min(meta['minimum_lambda'], float(ray['lam']))
    for element in fields[d+3].elements:
        points, weights = fields[d+3].quadrature(element, int(raw['quadrature_orders'][0]))
        integrate(element['geometry'], points, weights, False)
    for block in (d+2, d+3):
        mean = fields[block].mass()@fields[block].constant()/fields[0].volume
        result['gauge'][offsets[block]:offsets[block+1]] += state[-1]*mean
        result['gauge'][-1] += state[offsets[block]:offsets[block+1]]@mean
    return result, meta


def norm_record(audit, values, constraints):
    norms = audit.norms(values, constraints)
    return dict(blocks=norms.tolist(), maximum=float(max(norms)))


def physical_term_record(audit, values, constraints):
    # Strong-row physical loads are boundary reactions, not constraint
    # residuals.  Keep them explicitly instead of labelling their magnitude
    # as a constraint defect in the final generic Norms component.
    free_values = np.asarray(values, float).copy()
    strong = [abs(float(free_values[row])) for row in constraints]
    for row in constraints:
        free_values[row] = 0.
    return dict(free_physical_norms=norm_record(audit, free_values, constraints),
        maximum_strong_reaction_coefficient=max(strong, default=0.),
        scope='Free physical mass-dual blocks and gauge; constrained-row term loads retained as reactions, not a constraint residual')


def accumulation_bounds(raw):
    """Higham gamma_n bound for sums of frozen binary64 increments only."""
    u, offsets = F(1, 2**53), raw['offsets']
    result = {}
    for key in ('legacy_spatial',)+TERMS:
        counts = raw['addition_counts'][key]
        absolute = raw['absolute_accumulation'][key]
        bounds = []
        for count, magnitude in zip(counts, absolute):
            count = int(count)
            if count*u >= 1:
                raise ValueError('Unsupported gamma_n bound')
            bounds.append(count*u/(1-count*u)*rational(magnitude))
        blocks = []
        for block in range(len(offsets)-1):
            subset = bounds[offsets[block]:offsets[block+1]]
            blocks.append(float(max(subset, default=F(0))))
        result[key] = dict(maximum_coefficient_bound=float(max(bounds, default=F(0))),
            block_maximum_coefficient_bounds=blocks,
            maximum_number_of_additions=int(max(counts, default=0)))
    return dict(per_key=result,
        formula='|fl(sum_i x_i)-sum_i x_i| <= gamma_n sum_i |x_i|; gamma_n=n*u/(1-n*u),u=2^-53',
        absolute_sum='Native diagnostic supplies a conservative upper enclosure of the positive sum of frozen binary64 increments',
        scope='Summation error only. Does not bound basis/geometry/stress/profile/projection/product input errors, mass entries, or the entire residual. Not a roundoff-floor estimate')


def exact_row(field, point):
    element = field.element(point)
    degree = int(element['degree'])
    local = []
    for axis in (0, 2):
        xi = 2*(rational(point[axis])-rational(element['lower'][axis]))/(
            rational(element['upper'][axis])-rational(element['lower'][axis]))-1
        local.append([sum((coefficient*xi**power for power, coefficient in enumerate(poly)), F(0))
            for poly in POLYNOMIALS[:degree+1]])
    values = [F(0)]*field.size
    for j, (index, sign) in enumerate(zip(element['ids'], element['sign'])):
        values[int(index)] += int(sign)*local[0][j % (degree+1)]*local[1][j // (degree+1)]
    return values


def exact_norms(raw, audit, constraints, vectors):
    tables, factors = {}, []
    offsets, d = raw['offsets'], raw['dimension']
    for block, record in enumerate(raw['fields']):
        fixed = [row-offsets[block] for row in constraints if offsets[block] <= row < offsets[block+1]]
        factors.append(TensorMass(record, fixed))
    for name, values in vectors.items():
        records = []
        for block, mass in enumerate(factors):
            scale = audit.us/audit.H if block >= d+2 else audit.state_scale(block)*audit.us/audit.H
            records.append(mass.norm(values[offsets[block]:offsets[block+1]], scale))
        gauge = abs(rational(values[-1]))/rational(audit.us)**2
        bc = max((abs(rational(values[row]))/rational(audit.state_scale(next(b for b in range(len(factors))
            if offsets[b] <= row < offsets[b+1]))) for row in constraints), default=F(0))
        lower = max([Decimal(row['interval_decimal'][0]) for row in records]+[decimal(gauge), decimal(bc)])
        upper = max([Decimal(row['interval_decimal'][1]) for row in records]+[decimal(gauge), decimal(bc)])
        tables[name] = dict(blocks=records, gauge_decimal=str(decimal(gauge)),
            constraint_maximum_decimal=str(decimal(bc)), maximum_interval_decimal=[str(lower), str(upper)],
            below_original_2e_minus12_tolerance=bool(upper <= Decimal('2e-12')))
    return tables, factors


def contact_analysis(raw, fields, state, terms, constraints, audit, factors, inventory):
    if not raw['scalar_rates']:
        return dict(requested=False)
    d, offsets = raw['dimension'], raw['offsets']
    records, computed_rates = {}, {}
    for name, block in (('k', d), ('epsilon_with_induced_boundary_mass', d+1)):
        native = raw['scalar_rates'].get(name)
        if native is None:
            continue
        field, mass = fields[block], audit.mass[block]
        free, fixed = np.asarray(native['free_rows'], int), np.asarray(native['fixed_rows'], int)
        expected_fixed = sorted(row-offsets[block] for row in constraints if offsets[block] <= row < offsets[block+1])
        audit.require(fixed.tolist() == expected_fixed, 'Scalar contact fixed mass partition differs')
        audit.require(sorted(free.tolist()+fixed.tolist()) == list(range(field.size)), 'Scalar mass partition incomplete')
        boundary = np.asarray(native['boundary_rate'], float)
        factor = cho_factor(mass[np.ix_(free, free)], lower=True, check_finite=False)
        independent_spatial = sum((terms[key][offsets[block]:offsets[block+1]] for key in TERMS), np.zeros(field.size))
        rhs = -independent_spatial[free]-mass[np.ix_(free, fixed)]@boundary[fixed]
        rate = boundary.copy()
        rate[free] = cho_solve(factor, rhs, check_finite=False)
        computed_rates[name] = rate
        audit.near(native['rate'], rate, 'independent_consistent_mass_'+name+'_rate', TOLERANCE)
        balance = mass[np.ix_(free, np.arange(field.size))]@np.asarray(native['rate'])+np.asarray(native['spatial_load'])[free]
        audit.near(native['free_mass_residual'], balance, 'independent_native_mass_equation_'+name, TOLERANCE)
        component_sum = np.zeros(field.size)
        for key in TERMS:
            load = terms[key][offsets[block]:offsets[block+1]]
            component = np.zeros(field.size)
            component[free] = cho_solve(factor, -load[free], check_finite=False)
            audit.near(native['component_coefficient_rates'][key], component,
                'independent_mass_component_'+name+'_'+key, TOLERANCE)
            component_sum += component
        component_sum += np.asarray(native['component_coefficient_rates']['boundary_rate_mass'])
        audit.near(component_sum, native['rate'], 'component_mass_rate_sum_'+name, TOLERANCE)
        points = []
        for sample in native['physical_samples']:
            point = sample['coordinate']
            ids, phi, _ = field.rows(field.element(point, 1 if block == d+1 and point[2] == raw['parameters'][1] else 0), [point])
            value = float(phi[0]@state[offsets[block]+ids])
            physical_rate = float(phi[0]@rate[ids])
            audit.near(sample['value'], value, 'independent_scalar_contact_value_'+name, TOLERANCE)
            audit.near(sample['rate'], physical_rate, 'independent_scalar_contact_rate_'+name, TOLERANCE)
            row = dict(coordinate=point, value=value, native_value=sample['value'],
                independent_rate=physical_rate, native_rate=sample['rate'],
                component_rates=sample['component_rates'],
                largest_absolute_component=max(sample['component_rates'], key=lambda key: abs(sample['component_rates'][key])))
            if block == d and not np.any(boundary):
                weights = exact_row(field, point)
                rhs = [-float(x) for x in native['spatial_load']]
                row['fixed_native_load_exact_mass_rate'] = factors[block].functional(rhs, weights)
                row['exact_polynomial_value_decimal'] = str(decimal(sum((weights[i]*rational(state[offsets[block]+i]) for i in range(field.size)), F(0))))
            points.append(row)
        records[name] = dict(points=points, maximum_free_mass_balance=float(np.max(np.abs(balance), initial=0.)),
            physical_consistent_mass=True, lumped=False,
            boundary_mass_included=True)
    assumptions = raw['scalar_rates'].get('trace_rate_assumptions_verified', False)
    independent_assumptions = True
    trace_nodes = raw['scalar_rates'].get('trace_projection_node_checks', [])
    for node in trace_nodes:
        ray = inventory.ray(node['coordinate'], baseline=node['baseline'])
        audit.near([node['b'], node['lambda'], node['K']],
            [ray['b'], ray['lam'], float(ray['profile'].k)], 'independent_quiet_trace_normalization', TOLERANCE)
        independent_assumptions &= (ray['b'] > raw['robust_parameters'][3] and ray['lam'] == 1.)
        view = dict(raw, transition=dict(controls=raw['transition_controls'],
            weights=raw['transition_weights'], reference=raw['transition_reference']))
        gates = _diagnostic(view, state, fields, node['coordinate'], zero_body=True)['gates']
        audit.near([node['state_gate'], node['admissibility_gate'], node['d2']],
            gates[[2, 3, 6]], 'independent_quiet_trace_state_gates', TOLERANCE)
        independent_assumptions &= gates[2] == 1. and gates[3] == 1.
        independent_assumptions &= raw['parameters'][1] > float(ray['profile'].ze)
        independent_assumptions &= (ray['profile'].chi == 0 or raw['parameters'][1] > float(ray['profile'].rc))
    audit.require(not assumptions or independent_assumptions, 'Quiet trace assumptions failed independent normalization')
    derivative_check = None
    if assumptions and independent_assumptions:
        baseline, correction, designated = {}, {}, set()
        field, H = fields[d+1], raw['parameters'][1]
        for element in field.elements:
            if abs(element['lower'][2]-H) > 1e-11:
                continue
            nodes, transform = _segment_projection(element['degree'])
            values = []
            for node in nodes:
                point = np.array([element['lower'][0]+(node+1)*(element['upper'][0]-element['lower'][0])/2, 0., H])
                sample = point.copy(); sample[2] = raw['parameters'][3]
                dotK = fields[d].evaluate(computed_rates['k'], 0, sample)[0]
                derivatives = []
                for zero in (True, False):
                    ray = inventory.ray(point, baseline=zero)
                    E = float(np.asarray(ray['profile'].fields([H]), float)[2, 0])
                    derivatives.append(1.5*E/float(ray['profile'].k)*dotK)
                values.append([derivatives[0], derivatives[1]-derivatives[0]])
            projected = [transform@np.asarray(values).T[part] for part in range(2)]
            count = element['degree']+1
            for table, coefficients in zip((baseline, correction), projected):
                for index, sign, value in zip(element['ids'][:count], element['sign'][:count], coefficients):
                    index, value = int(index), float(sign*value)
                    if index in table:
                        audit.near(table[index], value, 'independent_shared_trace_rate_mode', TOLERANCE)
                    table[index] = value
            if element['lower'][0] >= raw['transition_controls'][4]-1e-11:
                designated.update(int(index) for index in element['ids'][:count])
        for index in designated:
            correction[index] = 0.
        expected = np.zeros(field.size)
        for index in baseline:
            expected[index] = baseline[index]+correction[index]
        native = np.asarray(raw['scalar_rates']['epsilon_with_induced_boundary_mass']['boundary_rate'])
        audit.near(native, expected, 'independent_complete_quiet_projected_trace_rate', TOLERANCE)
        derivative_check = dict(maximum_coefficient_error=float(np.max(np.abs(native-expected), initial=0.)),
            derivative='Separate baseline and current projection of1.5*Epre(H)/K*sampledKdot; designated correction modes zero; lambda1 and all state gates on their constant unit plateaus',
            g_E_zero_on_verified_branch=True, induced_boundary_mass_included=True)
    return dict(requested=True, scalar_mass_rates=records,
        full_scalar_DAE_rate_assessed=bool(assumptions and independent_assumptions),
        independent_projected_trace_rate=derivative_check,
        rate_scope=('Conditional fixed U and reduced Pi scalar subsystem. Both free scalar mass matrices are SPD; '
            'the verified quiet/unit/uncapped trace has g_E=0 and the prescribed sampled-K derivative is included through M_E,FB.'
            if assumptions else 'k-only instantaneous consistent-mass equation; current dependent-trace differential/subsystem solvability was not established'),
        pressure_scope='U and stored reduced Pi fixed; mechanical p/rho=Pi-2k/3 varies with k and remains in the forcing/sensors',
        original_fixed_contact=[1.6046496089512394, 0., 0.], constructed_zero_contact=[1., 0., 0.],
        no_positivity_or_implicit_root_claim=True)


def constructed_proof(raw, fields, state, constraints, trace, audit):
    if not raw['construct_contact']:
        return dict(constructed=False, saved_state_is_not_a_certified_zero_contact=True)
    block, offsets, K = raw['dimension'], raw['offsets'], raw['parameters'][5]
    field = fields[block]
    expected = {}
    for element in field.elements:
        audit.require(element['degree'] == 3, 'Constructed proof requires unchanged p3')
        one = [1, 1, 0, 0]
        if element['lower'][0] == -4. and element['upper'][0] == -3.:
            f = [0, 1, 1, 0]
        elif element['lower'][0] == 1. and element['upper'][0] == 2.:
            f = [1, 0, 1, 0]
        else:
            audit.require(-3. <= element['lower'][0] < element['upper'][0] <= 1., 'Unexpected contact plateau cell')
            f = one
        h = [1, 0, 1, 0] if element['lower'][2] == 0 else [0]*4
        if any(h):
            audit.require(element['upper'][2] == .01, 'Changed constructed radial support')
        for j, (index, sign) in enumerate(zip(element['ids'], element['sign'])):
            value = sign*K*(one[j % 4]*one[j // 4]-f[j % 4]*h[j // 4])
            if int(index) in expected:
                audit.require(expected[int(index)] == value, 'Constructed coefficient incompatible at a shared mode')
            expected[int(index)] = value
    audit.require(all(state[offsets[block]+i] == value for i, value in expected.items()),
        'Constructed state differs from exact0/+-K0 polynomial definition')
    audit.require(trace['relative_to_E0'] <= 1e-9, 'Constructed dependent trace is inconsistent')
    for row in constraints:
        if offsets[block] <= row < offsets[block+1]:
            audit.require(state[row] == constraints[row], 'Constructed state violates a static k condition')
    for block_id in (0, 1, 4, 5):
        audit.require(np.array_equal(state[offsets[block_id]:offsets[block_id+1]],
            np.asarray(raw['history'])[offsets[block_id]:offsets[block_id+1]]),
            'Constructed scalar diagnostic changed fixed velocity/reduced-pressure coefficients')
    return dict(constructed=True, coefficient_identity_exact=True,
        analytic_k_nonnegative=True, proof=(
            'On the left ramp f=3/4+xi/2-xi²/4 has derivative(1-xi)/2>=0 and endpoint values0,1; '
            'on the right f=3/4-xi/2-xi²/4 decreases from1to0; on the plateau f=1. '
            'g=((1+zeta)/2)² in the existing first radial cell and1above. Thus '
            'k=K0*((1-f)+f*g)>=0 throughout every polynomial element, with exact verified0/+-K0 coefficients.'),
        fixed_outer_K_and_sample_K=True,
        zero_contact_scope='Ground plateau including x=1. The original fixed row x=1.6046496089512394 is outside this plateau and must be reported separately.',
        trace_relative_error=trace['relative_to_E0'],
        trajectory_or_clipped_newton_update=False)


def analyze(raw, session, source, packet=None):
    audit = ResolutionAudit(raw, source)
    audit.require(raw['schema'] == 1 and raw['dimension'] == 2 and not raw['mapped'] and
        raw['eta'] == 0 and raw['upstream_transition'], 'Wrong R3 residual-resolution scope')
    audit.require(raw['robust_parameters'] == [.4, 0., .2, .1, 0., 1., 1.], 'Changed frozen outer/source law')
    audit.require(raw['model_variant'] == 'Regularized' and raw['inner_stress_model'] in ('Linear', 'SaturatedR2'),
        'Changed declared constitutive selectors')
    audit.require(raw['quadrature_orders'] == [8., 16.], 'Changed frozen source quadrature')
    audit.require(tuple(raw['spatial_sum_keys']) == TERMS, 'Native physical sum keys differ')
    state, history = np.asarray(raw['state']), np.asarray(raw['history'])
    audit.require(len(state) == 2165 and np.all(np.isfinite(state)), 'Wrong outer6 coefficient state')
    fields, offsets, d = audit.fields, raw['offsets'], raw['dimension']
    trace_ids = fields[d+1].boundary(lambda axis, side, e: axis == 2 and side == 0 and abs(e['lower'][2]-raw['parameters'][1]) < 1e-11)
    external = {int(row): float(value) for row, value in raw['reference_constraints']
                if int(row)-offsets[d+1] not in trace_ids or not offsets[d+1] <= int(row) < offsets[d+2]}
    effective = state.copy()
    for row, value in external.items():
        effective[row] = value
    audit.near(effective, raw['effective_state'], 'independent_external_lift', 1e-13)
    inventory = Inventory(raw, fields, effective, session)
    trace = trace_audit(raw, fields, effective, inventory)
    constraints = dict(external)
    constraints.update({offsets[d+1]+int(row): value for row, value in trace['projected_target_coefficients']})
    native_constraints = audit.constraints(raw['constraints'])
    audit.require(set(constraints) == set(native_constraints), 'Projected constraint rows differ')
    audit.near([constraints[i] for i in sorted(constraints)], [native_constraints[i] for i in sorted(constraints)],
        'independent_projected_constraint_targets', TOLERANCE)
    domain = native_domain_audit(raw, fields, effective, raw['scalar_domain'], trace, audit, inventory)
    terms, meta = independent_terms(raw, fields, effective)
    term_comparisons = {}
    for key, value in terms.items():
        native = np.asarray(raw['terms'][key])
        audit.near(native, value, 'independent_physical_term_'+key, TOLERANCE)
        term_comparisons[key] = dict(independent=physical_term_record(audit, value, constraints),
            native=physical_term_record(audit, native, constraints),
            difference=physical_term_record(audit, native-value, constraints),
            maximum_coefficient_error=float(np.max(np.abs(native-value), initial=0.)))
    spatial = sum((terms[key] for key in TERMS), np.zeros(len(state)))
    old_loads, weak_meta = weak_loads(raw, fields, effective)
    audit.near(spatial, old_loads, 'independent_term_sum_vs_complete_selected_operator', TOLERANCE)
    audit.near(raw['terms']['spatial_sum'], spatial, 'independent_complete_spatial_sum', TOLERANCE)
    audit.near(raw['native']['spatial_unmasked'], spatial, 'unchanged_native_spatial_residual', TOLERANCE)
    time_mass = np.zeros(len(state))
    for block in range(d+2):
        rows = np.arange(offsets[block], offsets[block+1])
        time_mass[rows] = raw['alpha']*(audit.mass[block]@(effective[rows]-history[rows]))
    audit.near(raw['terms']['time_mass'], time_mass, 'state_minus_history_before_physical_mass', TOLERANCE)
    total = spatial+time_mass
    for row, target in constraints.items():
        total[row] = state[row]-target
    audit.near(raw['native']['stage_constrained'], total, 'unchanged_native_coupled_stage_residual', TOLERANCE)
    native = np.asarray(raw['native']['stage_constrained'])
    constrained = np.asarray(raw['terms']['total_constrained'])
    cpp_wide = np.asarray(raw['extended_decimal']['total_constrained'], dtype=np.longdouble)
    cancelled = np.asarray(raw['extended_decimal']['total_constrained_cancelled'], dtype=np.longdouble)
    observations = dict(
        native=norm_record(audit, native, constraints),
        independent_terms=norm_record(audit, total, constraints),
        native_component_sum=norm_record(audit, constrained, constraints),
        native_minus_independent=norm_record(audit, native-total, constraints),
        native_minus_component_sum=norm_record(audit, native-constrained, constraints),
        native_minus_extended=norm_record(audit, native-np.asarray(cpp_wide, float), constraints),
        extended_cancellation_difference=norm_record(audit, np.asarray(cancelled-cpp_wide, float), constraints),
        native_time_minus_independent_time=physical_term_record(audit, np.asarray(raw['terms']['time_mass'])-time_mass, constraints))
    exact, factors = exact_norms(raw, audit, constraints, dict(native=native,
        extended=raw['extended_decimal']['total_constrained'],
        extended_cancelled=raw['extended_decimal']['total_constrained_cancelled']))
    exact_time = {}
    for block in range(d+2):
        rows = np.arange(offsets[block], offsets[block+1])
        # Native currently subtracts the binary64 state/history coefficients
        # before applying M. Exact rational matrix action freezes THAT delta.
        delta = effective[rows]-history[rows]
        mass_action = factors[block].action(delta)
        expected = np.array([float(rational(raw['alpha'])*value) for value in mass_action])
        exact_time[str(block)] = dict(maximum_native_stored_mass_action_difference=float(np.max(np.abs(
            np.asarray(raw['terms']['time_mass'])[rows]-expected), initial=0.)),
            frozen_delta_float64_LE_sha256=vector_sha(delta),
            scope='Observed native stored-mass/dot/cast difference from exact polynomial mass acting on the same already-rounded state-history delta; not an enclosing error bound')
    contact = contact_analysis(raw, fields, effective, terms, constraints, audit, factors, inventory)
    construction = constructed_proof(raw, fields, state, constraints, trace, audit)
    ranges = finite_ranges(raw, fields, effective, session, inventory) if raw['construct_contact'] else None
    packet_binding = None
    if packet is not None:
        payload = json.loads(Path(packet).read_text())
        audit.require(np.array_equal(history, payload['history']), 'Original h0 history packet changed')
        audit.require(raw['time'] == payload['time'] and raw['alpha'] == payload['alpha'], 'Target time/alpha changed')
        if not raw['construct_contact']:
            audit.require(np.array_equal(state, payload['state']), 'Saved diagnostic state changed')
        packet_binding = dict(path=str(Path(packet).resolve()), sha256=sha(packet))
        if not raw['construct_contact'] and payload.get('residual'):
            saved = np.asarray(payload['residual'])
            audit.near(native, saved, 'retained_and_current_same_state_equations', TOLERANCE)
            packet_binding['retained_residual_bitwise_equal'] = bool(np.array_equal(native, saved))
            packet_binding['current_minus_retained_residual'] = norm_record(audit, native-saved, constraints)
            packet_binding['retained_residual_physical_norms'] = norm_record(audit, saved, constraints)
    return dict(schema=1, phase='M8-R3 C/D fixed-state diagnostic', checks=audit.checks,
        success=True, failures=[], errors=audit.errors, source_packet=packet_binding,
        state_sha256=vector_sha(state), original_history_sha256=vector_sha(history),
        exact_target_tolerance=2e-12, inner_stress_model=raw['inner_stress_model'],
        time=raw['time'], alpha=raw['alpha'], term_comparisons=term_comparisons,
        observed_norms=observations, certified_exact_mass_norms=exact,
        frozen_accumulation_bounds=accumulation_bounds(raw),
        exact_polynomial_time_mass_comparisons=exact_time,
        independent_trace=trace, native_inventory_crosscheck=domain,
        physical_term_metadata=meta, independent_complete_weak_metadata=weak_meta,
        contact=contact, constructed_state_proof=construction,
        constructed_finite_ranges=ranges, outer_oracle=outer_evaluation_statistics(),
        scientific_limits=[
            'The2e-12 coupled root criterion is unchanged; no diagnostic result makes a historical failed candidate a root.',
            'Observed native/independent/extended arithmetic differences do not by themselves establish a roundoff floor or the cause of stagnation.',
            'The accumulation bound covers sums of frozen binary64 increments only; the exact-mass certificate fixes the declared residual coefficients.',
            'Native finite inventory and covariance/applicability observations do not prove a PDE or trajectory invariant region.',
            'A saved positive near-contact state or nonzero constraint residual is not a verified zero-contact positivity counterexample.',
            'The constructed zero contactx=1 differs from the retained original fixed contactx=1.6046496089512394; rate conclusions must retain that location scope.'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--packet', type=Path)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument('--producing-revision', help='Source epoch of the retained native executable; independent from analysis HEAD')
    parser.add_argument('--native-receipt', type=Path, help='Read-only parent runtime/source/binary fingerprint receipt')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Numerical resolution/contact analysis is authorized only on Linux')
    if args.output.exists():
        raise SystemExit('Refusing to overwrite an existing diagnostic result')
    started, wall = stamp(), time.monotonic()
    receipt = dict(schema=1, started_utc=started, source_revision=git(args.source, 'rev-parse', 'HEAD'),
        source_worktree_status=git(args.source, 'status', '--porcelain', '--untracked-files=no'),
        threads=THREADS, worker_count=1, native_report=dict(path=str(args.report.resolve()),
            sha256=sha(args.report), bytes=args.report.stat().st_size),
        session=dict(path=str(args.session.resolve()), sha256=sha(args.session)),
        producing_revision=args.producing_revision,
        native_runtime_receipt=(dict(path=str(args.native_receipt.resolve()), sha256=sha(args.native_receipt),
            bytes=args.native_receipt.stat().st_size) if args.native_receipt else None),
        analysis_sources={path.name: sha(path) for path in (Path(__file__),
            Path(__file__).with_name('cls_epawf_m8_r3_mass.py'),
            Path(__file__).with_name('cls_epawf_m8_r2_oracle.py'),
            Path(__file__).with_name('cls_epawf_m8_r2_preflight.py'))})
    try:
        if receipt['source_worktree_status']:
            raise ValueError('Analysis source has tracked working-tree changes')
        raw = json.loads(args.report.read_text())
        receipt.update(analyze(raw, args.session, args.source, args.packet))
    except Exception as error:
        receipt.update(success=False, failures=[dict(type=type(error).__name__,
            message=str(error), traceback=traceback.format_exc())])
    receipt.update(finished_utc=stamp(), elapsed_seconds=time.monotonic()-wall)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    save(args.output, receipt)
    print(json.dumps(dict(success=receipt['success'], checks=receipt.get('checks'),
        elapsed_seconds=receipt['elapsed_seconds'], output=str(args.output)), allow_nan=False))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
