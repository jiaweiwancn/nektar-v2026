#!/usr/bin/env python3
"""Independent M8-R1 t=0 representation audit; never launches a solver.

Reconstructs the complete unforced weak loads from native coefficient maps,
the indexed stress oracle and dimensional wall primitives. Source Gauss,
mass, matching/sensor, native GLL and dense reconstruction sets are finite
observations, not polynomial-everywhere or PDE invariant-region proofs.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)

# Import numerical modules only after fixing the numerical worker policy.
import numpy as np
from cls_epawf_m8_r1_cases import cost, RADIAL, physics, sha as case_sha
from cls_epawf_m8_state_coverage import (reconstruction_rules, point_sets,
    local_values, tensor_points, gll)
from cls_epawf_m8_native_oracle import (ray_source, selected_stress,
    direct as direct_energy, interface as interface_budget)
from cls_epawf_m7_energy_oracle import (Audit, _fields, _horizontal, _profile,
    _forcing, _sample, _velocity, _work)
from cls_epawf_m6_operator_oracle import _diagnostic, _segment_projection, coordinate
from validate_cls_epawf_m5 import physical_rows


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def vector_sha(value):
    return hashlib.sha256(np.asarray(value, dtype='<f8').tobytes()).hexdigest()


def selected_initial(raw):
    states, times = raw.get('accepted_states', []), raw.get('audit_times', [])
    if not states or len(states) != len(times) or times[0] != 0:
        raise ValueError('Missing actual first accepted t=0 coefficient state')
    state = np.asarray(states[0], float)
    if len(state) != raw['offsets'][-1]+1 or not np.all(np.isfinite(state)):
        raise ValueError('Invalid complete initial native state')
    if raw['dimension'] != 2 or raw['mapped'] or not raw['upstream_transition'] or raw['eta'] != 0:
        raise ValueError('Wrong physical preflight mode')
    if raw['model_variant'] != 'Regularized' or raw['robust_parameters'] != [.4, 0., .2, .1, 0., 1., 1.]:
        raise ValueError('Changed frozen M8-R1 model/source parameters')
    if any('|'+name+'/' in raw.get('signature', '') for name in
           ('BodyForceInner', 'BodyForceOuter', 'ScalarSourceInner', 'ScalarSourceOuter')):
        raise ValueError('Preflight requires the actual unforced control')
    return state


def transition_view(raw):
    return dict(raw, transition=dict(controls=raw['transition_controls'],
        weights=raw['transition_weights'], reference=raw['transition_reference']))


def add(load, offset, ids, phi, derivative, weight, source, flux):
    np.add.at(load, offset+ids, phi.T@(weight*source)+
              np.einsum('qai,qa,q->i', derivative, flux, weight))


def weak_loads(raw, fields, state, tangential=None, normal=None, outer=None):
    """Actual unforced alpha=0 operator; all pressure/gauge rows retained."""
    tangential = int(raw['quadrature_orders'][0]) if tangential is None else tangential
    normal = int(raw['quadrature_orders'][1]) if normal is None else normal
    outer = tangential if outer is None else outer
    view = deepcopy(raw)
    view['quadrature_orders'] = [tangential, normal]
    d, offsets = raw['dimension'], raw['offsets']
    nu, H, rough = raw['parameters'][:3]
    result = np.zeros(len(state))

    def integrate(identity, points, weights, inner, destruction=None, profile=None, force=None):
        velocity = [physical_rows(fields[0], state, offsets[b], points, identity) for b in range(d)]
        ids, phi, derivative = velocity[0][:3]
        U, G = np.zeros((len(points), 3)), np.zeros((len(points), 3, 3))
        for block, axis in enumerate((0, 2)):
            U[:, axis], G[:, axis] = velocity[block][3], velocity[block][4]
        ki, kp, kg, k, gk = physical_rows(fields[d], state, offsets[d], points, identity)
        block = d+2 if inner else d+3
        pi, pp, pg, pressure, _ = physical_rows(fields[block], state, offsets[block], points, identity)
        S = G+G.transpose(0, 2, 1)
        if inner:
            nt = float(profile.a)*np.maximum(points[:, 2]-float(profile.zv), 0.)
            R = 2*k[:, None, None]/3*np.eye(3)-nt[:, None, None]*S
            P = -np.einsum('qab,qab->q', R, G)
            stress = (nu+nt)[:, None, None]*S-pressure[:, None, None]*np.eye(3)
            dk, sk = nu+nt, P-destruction
            body = (points[:, 2] <= rough)[:, None]*force
        else:
            ei, ep, eg, epsilon, ge = physical_rows(fields[d+1], state, offsets[d+1], points, identity)
            chosen = [selected_stress(raw, g, kval, e, p) for g, kval, e, p in zip(G, k, epsilon, points)]
            R, nt = np.asarray([row[0] for row in chosen]), np.asarray([row[1] for row in chosen])
            P = -np.einsum('qab,qab->q', R, G)
            stress = -R+(2*k/3-pressure)[:, None, None]*np.eye(3)
            dk, sk = nt, P-epsilon
            ratio = k**1.5/(2.5*epsilon*points[:, 2])
            se = 1.44*epsilon/k*P-1.92*epsilon**2/k+.83*epsilon**2/k*np.maximum(ratio-1, 0)*ratio**2
            add(result, offsets[d+1], ei, ep, eg, weights, -se, nt[:, None]/1.3*ge-U*epsilon[:, None])
            body = np.zeros_like(U)
        for block, axis in enumerate((0, 2)):
            add(result, offsets[block], ids, phi, derivative, weights,
                -body[:, axis], stress[:, axis]-U[:, axis, None]*U)
        add(result, offsets[d], ki, kp, kg, weights, -sk, dk[:, None]*gk-U*k[:, None])
        add(result, offsets[d+2 if inner else d+3], pi, pp, pg, weights,
            np.trace(G, axis1=1, axis2=2), np.zeros_like(U))

    rays = 0
    for chart, horizontal in _horizontal(fields, d, H, tangential):
        ray = ray_source(view, fields, state, chart)
        force = _forcing(raw, fields, state, chart)
        groups = {}
        for index, point in enumerate(ray['points']):
            identity = fields[0].element(point)['geometry']
            groups.setdefault(identity, []).append(index)
        for identity, ids in groups.items():
            ids = np.asarray(ids)
            weights = horizontal*ray['metric']['baseArea']*ray['measure'][ids]
            integrate(identity, ray['points'][ids], weights, True,
                ray['destruction'][ids], ray['profile'], force)
        rays += 1
    for element in fields[d+3].elements:
        points, weights = fields[d+3].quadrature(element, outer)
        integrate(element['geometry'], points, weights, False)
    for block in (d+2, d+3):
        mean = fields[block].mass()@fields[block].constant()/fields[0].volume
        result[offsets[block]:offsets[block+1]] += state[-1]*mean
        result[-1] += state[offsets[block]:offsets[block+1]]@mean
    return result, dict(orders=[tangential, normal, outer], rays=rays,
                       scope='Complete independently reconstructed unforced alpha=0 weak loads')


def dual_norms(raw, fields, loads, state):
    d, offsets, us = raw['dimension'], raw['offsets'], raw['ustar']
    K, H, kappa = raw['parameters'][5], raw['parameters'][1], raw['parameters'][7]
    E0 = us**3/(kappa*H)
    constraints = {int(row): value for row, value in raw['reference_constraints']}
    result = []
    for block, field in enumerate(fields):
        free = np.array([index for index in range(field.size) if offsets[block]+index not in constraints], int)
        if not len(free):
            result.append(0.)
            continue
        value = loads[offsets[block]+free]
        mass = field.mass()[np.ix_(free, free)]
        inverse = np.linalg.solve(mass, value)
        square = float(value@inverse)
        if square < -1e-18:
            raise ValueError('Negative physical mass-dual square')
        scale = us if block < d else K if block == d else E0 if block == d+1 else 1.
        rate = us/H if block >= d+2 else scale*us/H
        result.append(math.sqrt(max(0., square))/(rate*math.sqrt(field.volume)))
    gauge = abs(loads[-1])/(us*us)
    bc = 0.
    for row, value in constraints.items():
        block = max(index for index, offset in enumerate(offsets) if offset <= row)
        scale = us if block < d else K if block == d else E0
        bc = max(bc, abs(state[row]-value)/scale)
    return result+[gauge, bc]


def eigen_example(R, points, gradients, k, sampled=None, nt=None):
    values, vectors = np.linalg.eigh((R+R.transpose(0, 2, 1))/2)
    index = int(np.argmin(values[:, 0]))
    vector = vectors[index, :, 0]
    if vector[np.argmax(np.abs(vector))] < 0:
        vector = -vector
    record = dict(chart=points[index].tolist(), covariance_minimum=float(values[index, 0]),
                  eigenvector=vector.tolist(), k=float(k[index]), gradient=gradients[index].tolist(),
                  strain=(gradients[index]+gradients[index].T).tolist())
    if sampled is not None:
        S = gradients[index]+gradients[index].T
        record.update(sampled_K=float(sampled[index]), k_over_sampled_K=float(k[index]/sampled[index]),
                      nu_A=float(nt[index]), strain_largest_eigenvalue=float(np.linalg.eigvalsh(S)[-1]),
                      r_R=float(3*nt[index]*np.linalg.eigvalsh(S)[-1]/(2*k[index])))
    return record


def inner_values(raw, fields, state, identity, points):
    d, offsets, rK = raw['dimension'], raw['offsets'], raw['parameters'][3]
    _, G = _velocity(fields, state, offsets, d, points, identity)
    k = _sample(fields, state, offsets, d, points, identity)[3]
    cache, sampled, viscosity = {}, [], []
    for point in points:
        key = float(point[0])
        if key not in cache:
            K = fields[d].evaluate(state, offsets[d], [point[0], 0., rK])[0]
            # Viscosity depends only on sampled K and the fixed wall law, not
            # chi or any transported-field floor. Flat geometry is explicit.
            root = math.sqrt(K)
            nu, rough = raw['parameters'][0], raw['parameters'][2]
            X = rough*root/(70*nu)
            exponent = max(.5-.4*X**.7, 1-.79*X**(-.28))
            zv = nu*10.7*(1-X**exponent)/root
            cache[key] = K, .2295*root, zv
        K, a, zv = cache[key]
        sampled.append(K)
        viscosity.append(a*max(point[2]-zv, 0.))
    sampled, nt = np.asarray(sampled), np.asarray(viscosity)
    R = 2*k[:, None, None]/3*np.eye(3)-nt[:, None, None]*(G+G.transpose(0, 2, 1))
    return G, k, sampled, nt, R


def inner_point_sets(raw, fields, state, rules):
    d, H, rK = raw['dimension'], raw['parameters'][1], raw['parameters'][3]
    for element in fields[d+2].elements:
        velocity_degree = fields[0].by_id[element['geometry']]['degree']
        sets = [('mass_volume_inner', np.polynomial.legendre.leggauss(8)[0]),
                ('reconstruction_native_inner', gll(rules[element['geometry']][0])),
                ('reconstruction_dense_inner', np.linspace(-1., 1., 2*velocity_degree+7))]
        for name, nodes in sets:
            yield name, element['geometry'], tensor_points(element, fields[0].axes, [nodes, nodes])
    for chart, _ in _horizontal(fields, d, H, int(raw['quadrature_orders'][0])):
        ray = ray_source(raw, fields, state, chart)
        groups = {}
        for point in ray['points']:
            identity = fields[0].element(point)['geometry']
            groups.setdefault(identity, []).append(point)
        for identity, points in groups.items():
            yield 'source_volume_inner', identity, np.asarray(points)
        for name, distance in (('physical_sample_K', rK), ('interface_inner', H)):
            point = np.asarray([[chart[0], 0., distance]])
            identity = fields[0].element(point[0], -1 if distance == H else 0)['geometry']
            yield name, identity, point


def finite_ranges(raw, fields, state, session):
    d, H, us, qstar, kappa = raw['dimension'], raw['parameters'][1], raw['ustar'], raw['parameters'][6], raw['parameters'][7]
    rules, issue = reconstruction_rules(raw, session)
    if issue:
        raise ValueError(issue)
    frame = np.column_stack([raw['robust_reference_flow'],
        np.cross(raw['robust_reference_normal'], raw['robust_reference_flow']), raw['robust_reference_normal']])
    center = np.zeros((3, 3)); center[0, 2] = 1.
    result = dict(outer={}, inner={}, outer_maximum_gradient_entry_deviation=0.,
                  outer_k_over_K0=[math.inf, -math.inf], outer_epsilon_over_epsilon0=[math.inf, -math.inf])
    for name, identity, points in point_sets(raw, fields, rules):
        G, k, epsilon, distance = local_values(raw, fields, state, identity, points)
        normalized = np.einsum('ai,qab,bj->qij', frame, G, frame)/(us/(kappa*distance))[:, None, None]
        deviation = np.max(np.abs(normalized-center), axis=(1, 2))
        kr, er = k/(qstar*us*us), epsilon/(us**3/(kappa*distance))
        strain, rotation = G+G.transpose(0, 2, 1), G-G.transpose(0, 2, 1)
        qs = k/epsilon*np.sqrt(np.sum(strain*strain, axis=(1, 2))/2)
        qo = k/epsilon*np.sqrt(np.sum(rotation*rotation, axis=(1, 2))/2)
        stress = np.asarray([selected_stress(raw, g, kval, e, p)[0] for g, kval, e, p in zip(G, k, epsilon, points)])
        eigen = eigen_example(stress, points, G, k)
        row = result['outer'].setdefault(name, dict(count=0, minimum_k=math.inf, minimum_epsilon=math.inf,
            maximum_gradient_entry_deviation=0., k_over_K0=[math.inf, -math.inf],
            epsilon_over_epsilon0=[math.inf, -math.inf], covariance_minimum=math.inf,
            strain_invariant=[math.inf, -math.inf], rotation_invariant=[math.inf, -math.inf],
            physical_max_invariant=[math.inf, -math.inf]))
        row['count'] += len(points)
        row['minimum_k'], row['minimum_epsilon'] = min(row['minimum_k'], float(min(k))), min(row['minimum_epsilon'], float(min(epsilon)))
        for key, values in (('k_over_K0', kr), ('epsilon_over_epsilon0', er),
                            ('strain_invariant', qs), ('rotation_invariant', qo), ('physical_max_invariant', np.maximum(qs, qo))):
            row[key] = [min(row[key][0], float(min(values))), max(row[key][1], float(max(values)))]
        row['maximum_gradient_entry_deviation'] = max(row['maximum_gradient_entry_deviation'], float(max(deviation)))
        if eigen['covariance_minimum'] < row['covariance_minimum']:
            row.update(covariance_minimum=eigen['covariance_minimum'], worst_covariance=dict(eigen, geometry=identity))
        result['outer_maximum_gradient_entry_deviation'] = max(result['outer_maximum_gradient_entry_deviation'], float(max(deviation)))
        for key, values in (('outer_k_over_K0', kr), ('outer_epsilon_over_epsilon0', er)):
            result[key] = [min(result[key][0], float(min(values))), max(result[key][1], float(max(values)))]
    for name, identity, points in inner_point_sets(raw, fields, state, rules):
        G, k, K, nt, stress = inner_values(raw, fields, state, identity, points)
        example = eigen_example(stress, points, G, k, K, nt)
        # Exact neutral AWF derivative, with the fixed incoming direction.
        reference = (us*us/(raw['parameters'][0]+nt))[:, None, None]*np.outer(
            raw['robust_reference_flow'], raw['robust_reference_normal'])
        relative = np.max(np.abs(G-reference), axis=(1, 2))/(us*us/(raw['parameters'][0]+nt))
        row = result['inner'].setdefault(name, dict(count=0, minimum_k=math.inf,
            minimum_k_over_sampled_K=math.inf, covariance_minimum=math.inf,
            maximum_relative_AWF_gradient_error=0., maximum_r_R=-math.inf))
        row['count'] += len(points)
        row['minimum_k'] = min(row['minimum_k'], float(min(k)))
        row['minimum_k_over_sampled_K'] = min(row['minimum_k_over_sampled_K'], float(min(k/K)))
        row['maximum_relative_AWF_gradient_error'] = max(row['maximum_relative_AWF_gradient_error'], float(max(relative)))
        largest = np.linalg.eigvalsh(G+G.transpose(0, 2, 1))[:, -1]
        row['maximum_r_R'] = max(row['maximum_r_R'], float(max(3*nt*largest/(2*k))))
        if example['covariance_minimum'] < row['covariance_minimum']:
            row.update(covariance_minimum=example['covariance_minimum'], worst_covariance=dict(example, geometry=identity))
    result['historical_radius001_kEPS99_101_covered'] = (
        result['outer_maximum_gradient_entry_deviation'] <= .001 and
        .99 <= result['outer_k_over_K0'][0] <= result['outer_k_over_K0'][1] <= 1.01 and
        .99 <= result['outer_epsilon_over_epsilon0'][0] <= result['outer_epsilon_over_epsilon0'][1] <= 1.01)
    result['scope'] = 'Finite source/mass/interface/transfer/sensor/native-GLL/dense outer sets and finite inner source/mass/K-sample/interface/native-GLL/dense sets; one-sided element gradients retained'
    return result


def trace_audit(raw, fields, state):
    view = transition_view(raw)
    d, offsets, H = raw['dimension'], raw['offsets'], raw['parameters'][1]
    baseline_start = raw['transition_controls'][4]
    baseline, correction, baseline_ids, samples = {}, {}, set(), []
    for element in fields[d+1].elements:
        if abs(element['lower'][2]-H) > 1e-12:
            continue
        nodes, transform = _segment_projection(element['degree'])
        x = element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
        data = []
        for value in x:
            point = np.array([value, 0., H])
            diagnostic = _diagnostic(view, state, fields, point, zero_body=True)
            ray = ray_source(raw, fields, state, point)
            if ray['lam'] != 1. or min(ray['local_k']/float(ray['profile'].k)) <= raw['robust_parameters'][2]:
                raise ValueError('Initial trace left the proven inactive source-cap plateau')
            profile, chi = ray['profile'], float(diagnostic['gates'][0])
            ea = float(profile.k)**1.5/(2.55*max(H, float(profile.ze)))
            target = float(profile.fields([H])[2, 0])
            epsilon = fields[d+1].evaluate(state, offsets[d+1], point, 1)[0]
            data.append([ea, target-ea])
            samples.append(dict(chart=point.tolist(), chi=chi, epsilon_over_epsilon0=epsilon/(raw['ustar']**3/(raw['parameters'][7]*H)),
                                lambda_value=ray['lam'], b=ray['b']))
        data = np.asarray(data).T
        maps = [transform@data[part] for part in range(2)]
        count = element['degree']+1
        for table, values in zip((baseline, correction), maps):
            for index, sign, value in zip(element['ids'][:count], element['sign'][:count], values):
                index, value = int(index), float(sign*value)
                if index in table and abs(table[index]-value) > 2e-10:
                    raise ValueError('Inconsistent shared initial epsilon trace projection')
                table[index] = value
        if coordinate(view, element['lower']) >= baseline_start-1e-11:
            baseline_ids.update(element['ids'][:count])
    for index in baseline_ids:
        correction[index] = 0.
    error = max(abs(state[offsets[d+1]+index]-baseline[index]-correction[index]) for index in baseline)
    zero = [row['epsilon_over_epsilon0'] for row in samples if row['chi'] == 0.]
    return dict(maximum_trace_coefficient_error=error, relative_to_E0=error/(raw['ustar']**3/(raw['parameters'][7]*H)),
                matching_rows=len(baseline), chi_zero_epsilon_ratio_range=[min(zero), max(zero)] if zero else None,
                samples=samples, source_plateau_verified=True,
                interpretation='Actual dependent epsilon trace with unforced sensor groups and constrained baseline-face correction; no pressure-continuity condition')


def verify(raw, session, ranges_only=False):
    state, audit = selected_initial(raw), Audit()
    fields, d, offsets = _fields(raw), raw['dimension'], raw['offsets']
    ranges = finite_ranges(raw, fields, state, session)
    trace = trace_audit(raw, fields, state)
    audit.require(trace['relative_to_E0'] <= 1e-9, 'Initial dependent epsilon trace is inconsistent')
    direct = direct_energy(raw, fields, state, audit)
    native = raw['state_audits'][0]
    audit.near(native['robust_source'], direct['robust_source'], 'R1_initial_independent_source_budget')
    audit.near(native['source_volume_inner'][3], direct['inner_covariance_minimum'], 'R1_initial_inner_source_covariance')
    audit.near(native['source_volume_outer'][3], direct['outer_covariance_minimum'], 'R1_initial_outer_source_covariance')
    traction = interface_budget(raw, fields, state)
    audit.near(native['interface'], traction, 'R1_initial_complete_traction_and_k_flux')
    minimum_inner = min([direct['inner_covariance_minimum']]+[row['covariance_minimum'] for row in ranges['inner'].values()])
    minimum_outer = min(row['covariance_minimum'] for row in ranges['outer'].values())
    positive = (min(row['minimum_k'] for row in ranges['inner'].values()) > 0 and
                min(row['minimum_k'] for row in ranges['outer'].values()) > 0 and
                min(row['minimum_epsilon'] for row in ranges['outer'].values()) > 0)
    result = dict(schema=1, state_selection='First stored accepted t=0 state; original report failure/status retained in provenance',
        coefficient_float64_LE_sha256=vector_sha(state), fields=[field.size for field in fields], unknowns=len(state),
        maximum_coefficient_difference_from_frozen_reference=float(max(abs(state-np.asarray(raw['fixed_reference'])))),
        reference_identity_note='Zero perturbation may still repeat the dependent-trace Lift after reference freezing. Actual accepted coefficients are used throughout; reference equality is measured, not assumed.',
        finite_ranges=ranges, independent_trace=trace, independent_source_energy=direct,
        independent_interface=traction.tolist(), native_initial_audit=native,
        initial_finite_admissibility=dict(positive_transported_k_and_outer_epsilon=positive,
            inner_covariance_minimum=minimum_inner, outer_covariance_minimum=minimum_outer,
            passed=positive and minimum_inner >= 0 and minimum_outer >= 0,
            scope='Strict finite covariance/positivity observations; no clipping or field floor, no continuum/polynomial-everywhere conclusion'),
        rigorous_joint_constitutive_qualification=False,
        scope='Bounded representation preflight, not a trajectory or physical terrain validation')
    if not ranges_only:
        loads, detail = weak_loads(raw, fields, state)
        native_loads = np.asarray(raw['raw_reference_loads'])
        audit.near(native_loads, loads, 'R1_independent_complete_unanchored_t0_weak_loads', 2e-10)
        norms = dual_norms(raw, fields, loads, state)
        audit.near(raw['raw_reference_norms'], norms, 'R1_independent_physical_block_dual_norms', 2e-9)
        result['raw_unanchored_residual'] = dict(native_block_norms=raw['raw_reference_norms'], independent_block_norms=norms,
            maximum_complete_weak_load_error=float(max(abs(loads-native_loads))),
            complete_weak_load_float64_LE_sha256=vector_sha(loads), gauge=norms[-2], essential_and_trace=norms[-1],
            block_names=['u', 'v', 'k', 'outer_epsilon', 'inner_pressure_divergence', 'outer_pressure_divergence', 'common_gauge', 'constraints'],
            quadrature=detail)
        changes = []
        native_orders = [int(raw['quadrature_orders'][0]), int(raw['quadrature_orders'][1]), int(raw['quadrature_orders'][0])]
        for orders in ((8, 16, 8), (8, 24, 8), (12, 24, 12)):
            if list(orders) == native_orders:
                continue
            other, metadata = weak_loads(raw, fields, state, *orders)
            changes.append(dict(orders=list(orders), maximum_weak_load_difference=float(max(abs(other-loads))),
                block_norms=dual_norms(raw, fields, other, state),
                free_load_difference_block_norms=dual_norms(raw, fields, other-loads, state),
                interpretation='Independent integration change at identical coefficient state; no new native solver/time-map execution'))
        result['same_state_quadrature_sensitivity'] = changes
    result.update(implementation_audit_passed=True, independent_checks=audit.checks, errors=audit.errors)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True, help='Clean archived source providing the original independent oracles')
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranges-only', action='store_true')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific preflight reconstruction belongs on the authoritative Linux server')
    if args.output.exists():
        raise FileExistsError('Preserve previous preflight receipts: '+str(args.output))
    raw_sha, xml_sha, plan_sha = sha(args.raw), sha(args.session), sha(args.plan)
    raw = json.loads(args.raw.read_text())
    receipt = dict(schema=1, started_utc=datetime.now(timezone.utc).isoformat(),
        raw=dict(path=str(args.raw), sha256=raw_sha, archived_failure=raw.get('failure'), archived_time=raw.get('time')),
        session=dict(path=str(args.session), sha256=xml_sha), driver_sha256=sha(__file__),
        plan=dict(path=str(args.plan), sha256=plan_sha),
        case_generator=dict(path=str(Path(sys.modules['cls_epawf_m8_r1_cases'].__file__).resolve()),
                            sha256=case_sha(sys.modules['cls_epawf_m8_r1_cases'].__file__)),
        numerical_threads=THREADS, affinity=sorted(os.sched_getaffinity(0)), solver_invoked=False,
        imported_oracles={})
    for name in ('cls_epawf_m8_state_coverage', 'cls_epawf_m8_native_oracle',
                 'cls_epawf_m7_energy_oracle', 'cls_epawf_m6_operator_oracle',
                 'validate_cls_epawf_m1', 'validate_cls_epawf_m2', 'validate_cls_epawf_m4',
                 'validate_cls_epawf_m5', 'validate_cls_epawf_m5_geometry'):
        path = Path(sys.modules[name].__file__).resolve()
        if not path.is_relative_to(args.source.resolve()):
            raise ValueError('Scientific oracle imported outside the frozen source: '+str(path))
        receipt['imported_oracles'][name] = dict(path=str(path), sha256=sha(path))
    try:
        source_head = subprocess.check_output(['git', '-C', str(args.source), 'rev-parse', 'HEAD'], text=True).strip()
        source_status = subprocess.check_output(['git', '-C', str(args.source), 'status', '--porcelain'], text=True).strip()
        if source_head != 'ab5c3c6fee19761e8d88692e8e6c331e5889234e' or source_status:
            raise ValueError('Independent oracle source is not the immutable reviewed ab5 baseline')
        receipt['oracle_source'] = dict(path=str(args.source), commit=source_head, clean=True)
        plan = json.loads(args.plan.read_text())
        if physics(ET.parse(args.session).getroot(), True) != plan['physics']:
            raise ValueError('Preflight changed an original physical/function/solver-policy input')
        radial = sorted({value for element in raw['fields'][0]['mesh'] for value in (element['lower'][2], element['upper'][2])})
        matching = [name for name, knots in RADIAL.items() if radial == knots]
        if len(matching) != 1:
            raise ValueError('Native report uses neither of the two frozen radial representations')
        expected = cost(RADIAL[matching[0]])
        if [field['size'] for field in raw['fields']] != expected['field_sizes'] or len(raw['state']) != expected['unknowns']:
            raise ValueError('Native unknown/field counts differ from the static cost contract')
        receipt['actual_representation'] = matching[0]
        receipt['cost'] = expected
        receipt['result'] = verify(raw, args.session, args.ranges_only)
        receipt['status'] = 'passed'
    except Exception as error:
        receipt.update(status='audit_failed', failure=dict(type=type(error).__name__, message=str(error)))
    receipt['input_hashes_unchanged'] = sha(args.raw) == raw_sha and sha(args.session) == xml_sha and sha(args.plan) == plan_sha
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=receipt['status'], output=str(args.output),
        finite_admissibility=receipt.get('result', {}).get('initial_finite_admissibility'),
        outer_ranges={key:value for key,value in receipt.get('result', {}).get('finite_ranges', {}).items()
                      if key.startswith('outer_')}, failure=receipt.get('failure'))), flush=True)
    return 0 if receipt['status'] == 'passed' and receipt['input_hashes_unchanged'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
