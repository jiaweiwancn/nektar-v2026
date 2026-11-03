#!/usr/bin/env python3
"""Independent R2 representation/selected-inner weak audit; no solver calls.

The unchanged R1 operator is reconstructed independently, then the indexed
selected-stress increment is assembled at every actual inner source point.
No anchor or artificial operator-reference correction is introduced.
Numerical invocation belongs to the authoritative Linux environment.
"""
import argparse
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

from cls_epawf_m8_r1_preflight import (weak_loads as linear_weak_loads,
    dual_norms, add, vector_sha, inner_values as linear_inner_values,
    inner_point_sets, eigen_example)
from cls_epawf_m8_r2_oracle import selected_inner, selected_outer as selected_stress, bind_outer_oracle
from cls_epawf_m8_native_oracle import (ray_source, normal_rule,
    direct as linear_direct)
from cls_epawf_m7_energy_oracle import (Audit, _fields, _horizontal, _profile,
    _forcing, _sample, _velocity)
from cls_epawf_m6_operator_oracle import _segment_projection, coordinate
from cls_epawf_m8_state_coverage import (reconstruction_rules, point_sets,
    local_values, tensor_points, gll)
from cls_epawf_m8_r2_cases import (RADIAL, cost, normalized_physics,
    REVIEWED_SOURCE, INNER_LAWS)
from validate_cls_epawf_m5 import physical_rows
from validate_cls_epawf_m2 import Profiles

# R2-only coefficient evaluation; original source/evidence are unchanged.
linear_weak_loads = bind_outer_oracle(linear_weak_loads)
linear_direct = bind_outer_oracle(linear_direct)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def real(value):
    array = np.asarray(value)
    if np.max(np.abs(array.imag), initial=0.) > 1e-12:
        raise ValueError('Unexpected complex value in a real-state audit')
    return array.real


def selected_values(raw, gradients, k, viscosity):
    values = [selected_inner(raw, g, kval, nt)
              for g, kval, nt in zip(gradients, k, viscosity)]
    return values


def weak_loads(raw, fields, state, tangential=None, normal=None, outer=None):
    """Complete unforced alpha=0 raw loads, including every pressure/gauge row.

    Importable by the independent active native probe. Arbitrary fixed-chi
    flat 2D states are permitted; the CLI adds the stricter R2 t0 controls.
    Scalar diffusivities, roughness force and capped preprofiles are unchanged.
    """
    if raw['dimension'] != 2 or raw.get('mapped'):
        raise ValueError('R2 preflight operator is explicitly flat and 2D')
    result, metadata = linear_weak_loads(raw, fields, state, tangential, normal, outer)
    if raw.get('inner_stress_model', 'Linear') == 'Linear':
        return result, dict(metadata, inner_stress_model='Linear', selected_increment=0.)
    tangential = int(raw['quadrature_orders'][0]) if tangential is None else tangential
    normal = int(raw['quadrature_orders'][1]) if normal is None else normal
    view = dict(raw, quadrature_orders=[tangential, normal])
    d, offsets, H = raw['dimension'], raw['offsets'], raw['parameters'][1]
    increment = np.zeros(len(state))
    active_volume, total_volume = 0., 0.
    for chart, horizontal in _horizontal(fields, d, H, tangential):
        ray = ray_source(view, fields, state, chart)
        nt = float(ray['profile'].a)*np.maximum(ray['r']-float(ray['profile'].zv), 0.)
        groups = {}
        for index, point in enumerate(ray['points']):
            groups.setdefault(fields[0].element(point)['geometry'], []).append(index)
        for identity, indices in groups.items():
            indices = np.asarray(indices)
            points = ray['points'][indices]
            weights = horizontal*ray['metric']['baseArea']*ray['measure'][indices]
            velocity = [physical_rows(fields[0], state, offsets[b], points, identity) for b in range(d)]
            ids, phi, derivative = velocity[0][:3]
            G = np.zeros((len(points), 3, 3))
            for block, axis in enumerate((0, 2)):
                G[:, axis] = velocity[block][4]
            ki, kp, kg, k, _ = physical_rows(fields[d], state, offsets[d], points, identity)
            values = selected_values(raw, G, k, nt[indices])
            F = np.asarray([real(row['turbulent']) for row in values])
            S = G+G.transpose(0, 2, 1)
            delta_F = F-nt[indices, None, None]*S
            delta_P = np.einsum('qab,qab->q', delta_F, G)
            for block, axis in enumerate((0, 2)):
                add(increment, offsets[block], ids, phi, derivative, weights,
                    np.zeros(len(points)), delta_F[:, axis])
            add(increment, offsets[d], ki, kp, kg, weights, -delta_P,
                np.zeros((len(points), 3)))
            total_volume += float(np.sum(weights))
            active_volume += float(np.sum(weights[[row['active'] for row in values]]))
    return result+increment, dict(metadata,
        inner_stress_model=raw.get('inner_stress_model', 'Linear'),
        selected_increment_float64_LE_sha256=vector_sha(increment),
        maximum_selected_weak_increment=float(np.max(np.abs(increment))),
        inner_source_volume=total_volume, saturation_active_source_volume=active_volume,
        source_diffusion_and_preprofiles='unchanged', anchor_added=False)


def selected_direct(raw, fields, state, audit):
    """Same-stress energy/source account with actual selected inner tensors."""
    result = linear_direct(raw, fields, state, audit)
    d, offsets, nu, H = raw['dimension'], raw['offsets'], raw['parameters'][0], raw['parameters'][1]
    selected_covariance, active_volume, positive_volume, zero_volume = math.inf, 0., 0., 0.
    ratios, radial_minimum, production_change = [], math.inf, 0.
    maximum_ratio, worst = -math.inf, None
    for chart, horizontal in _horizontal(fields, d, H, int(raw['quadrature_orders'][0])):
        ray = ray_source(raw, fields, state, chart)
        weights = horizontal*ray['metric']['baseArea']*ray['measure']
        G, pressure = [], []
        for point in ray['points']:
            G.append(_velocity(fields, state, offsets, d, [point])[1][0])
            pressure.append(_sample(fields, state, offsets, d+2, [point])[3][0])
        G, pressure, k = np.asarray(G), np.asarray(pressure), ray['local_k']
        nt = float(ray['profile'].a)*np.maximum(ray['r']-float(ray['profile'].zv), 0.)
        values = selected_values(raw, G, k, nt)
        R = np.asarray([real(row['stress']) for row in values])
        F = np.asarray([real(row['turbulent']) for row in values])
        P = np.asarray([float(real(row['production'])) for row in values])
        S, divergence = G+G.transpose(0, 2, 1), np.trace(G, axis1=1, axis2=2)
        stress = nu*S+F-pressure[:, None, None]*np.eye(3)
        direct = np.einsum('qab,qab->q', G, stress)-(P-ray['destruction'])
        reduced = nu*np.sum(S*S, axis=(1, 2))/2+(2*k/3-pressure)*divergence+ray['destruction']
        audit.near(direct, reduced, 'R2_selected_inner_full_strain_energy_cancellation', 3e-12)
        audit.near(P, np.einsum('qab,qab->q', F, G)-2*k/3*divergence,
            'R2_selected_inner_weak_divergence_production', 3e-12)
        linear_P = np.einsum('qab,qab->q', nt[:, None, None]*S, G)-2*k/3*divergence
        production_change += float(weights@(P-linear_P))
        selected_covariance = min(selected_covariance, float(np.min(np.linalg.eigvalsh(R)[:, 0])))
        for index, row in enumerate(values):
            radial_minimum = min(radial_minimum, float(real(row['radial'])))
            if row['active']:
                active_volume += float(weights[index])
            if k[index] > 0:
                positive_volume += float(weights[index])
                ratio = float(real(row['a']/row['c']))
                ratios.append(ratio)
                if ratio > maximum_ratio:
                    maximum_ratio = ratio
                    worst = dict(chart=ray['points'][index].tolist(), k=float(k[index]),
                        sampled_K=float(ray['profile'].k), nu_A=float(nt[index]),
                        G=G[index].tolist(), S=S[index].tolist(), a_over_c=ratio,
                        t=float(real(row['t'])), radial_tangent=float(real(row['radial'])))
            else:
                zero_volume += float(weights[index])
    result['wall_sources'][3] += production_change
    result['inner_covariance_minimum'] = selected_covariance
    result.update(inner_stress_model=raw.get('inner_stress_model', 'Linear'),
        inner_stress=[active_volume, positive_volume, min(ratios) if ratios else 0.,
            max(ratios) if ratios else 0., zero_volume, radial_minimum],
        inner_constitutive_modification=[production_change, production_change, 0.],
        inner_stress_max_ratio_point=worst,
        energy_scope='Actual selected full stress; production includes -(2k/3)*div(U); unchanged scalar diffusion and capped preprofiles')
    return result


def selected_interface(raw, fields, state):
    """Full total-traction jump with nu_w=nu, nu_o=0, plus scalar k flux."""
    d, offsets, H, nu = raw['dimension'], raw['offsets'], raw['parameters'][1], raw['parameters'][0]
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
        selected = selected_inner(raw, Gi[0], k, nt)
        tensor, outer_nut = selected_stress(raw, Go[0], k, epsilon, point)
        inner_stress = nu*(Gi[0]+Gi[0].T)+real(selected['turbulent'])-pi*np.eye(3)
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


class Inventory:
    """Independent finite inventory with freshly state-dependent full rays.

    Repeated source/trace anchors preserve their call-set identity. Normal
    points are reconstructed from the actual candidate state's profile and
    cutoff. A frozen seed ray is included for prescribed-load quadrature.
    """
    def __init__(self, raw, fields, state, session):
        self.raw, self.fields, self.state = raw, fields, np.asarray(state, float)
        self.d, self.H, self.offsets = raw['dimension'], raw['parameters'][1], raw['offsets']
        self.rules, issue = reconstruction_rules(raw, session)
        if issue:
            raise ValueError(issue)
        self.ray_cache = {}
        self.ray_summary = {}

    def ray(self, chart, baseline=False, seed=False):
        key = (float(chart[0]), float(chart[1]), bool(baseline), bool(seed))
        if key not in self.ray_cache:
            view, state = self.raw, self.state
            if baseline or seed:
                view = dict(view, upstream_transition=False, chi=0.)
            if seed:
                state = np.zeros_like(self.state)
                state[self.offsets[self.d]:self.offsets[self.d+1]] = (
                    self.raw['parameters'][5]*self.fields[self.d].constant())
            self.ray_cache[key] = ray_source(view, self.fields, state,
                np.array([chart[0], chart[1], self.H]))
        return self.ray_cache[key]

    def grouped(self, name, points, inner=True):
        groups = {}
        for point in points:
            side = (-1 if inner else 1) if point[2] == self.H else 0
            identity = self.fields[0].element(point, side)['geometry']
            groups.setdefault(identity, []).append(point)
        for identity, values in groups.items():
            yield name, identity, np.asarray(values), inner

    def trace_charts(self):
        points = set()
        for element in self.fields[self.d+1].elements:
            if abs(element['lower'][2]-self.H) > 1e-12:
                continue
            nodes, _ = _segment_projection(element['degree'])
            for value in nodes:
                points.add((element['lower'][0]+(value+1)*(element['upper'][0]-element['lower'][0])/2, 0., self.H))
        for point in sorted(points):
            yield np.asarray(point)

    def ray_sets(self, name, chart, baseline=False, seed=False):
        ray = self.ray(chart, baseline, seed)
        profile = ray['profile']
        row = self.ray_summary.setdefault(name, dict(rays=0, points=0,
            minimum_sampled_K=math.inf, lambda_range=[math.inf, -math.inf],
            b_range=[math.inf, -math.inf], chi_range=[math.inf, -math.inf],
            cutoff_range=[math.inf, -math.inf], minimum_dependent_epsilon_H=math.inf))
        row['rays'] += 1; row['points'] += len(ray['points'])
        row['minimum_sampled_K'] = min(row['minimum_sampled_K'], float(profile.k))
        dependent = float(ray['lam']*np.asarray(profile.fields([self.H]), float)[2, 0])
        row['minimum_dependent_epsilon_H'] = min(row['minimum_dependent_epsilon_H'], dependent)
        for field, value in (('lambda_range', ray['lam']), ('b_range', ray['b']),
                             ('chi_range', float(profile.chi)),
                             ('cutoff_range', float(profile.rc) if profile.chi > 0 else 0.)):
            row[field] = [min(row[field][0], value), max(row[field][1], value)]
        yield from self.grouped(name, ray['points'])
        # These explicit state-dependent endpoints supplement the actual
        # Gauss ray. They never replace any normalization point.
        boundaries = {0., self.H}
        for value in (profile.h, profile.zv, profile.ze, profile.rm,
                      profile.rc if profile.chi > 0 else 0.):
            if 0 <= value <= self.H:
                boundaries.add(float(value))
        endpoints = np.array([[chart[0], chart[1], value] for value in sorted(boundaries)])
        yield from self.grouped(name+'_profile_breakpoints', endpoints)

    def fixed_sets(self):
        for name, identity, points in point_sets(self.raw, self.fields, self.rules):
            yield name, identity, points, False
        for element in self.fields[self.d+2].elements:
            degree = self.fields[0].by_id[element['geometry']]['degree']
            for name, nodes in (
                ('mass_volume_inner', np.polynomial.legendre.leggauss(8)[0]),
                ('reconstruction_native_inner', gll(self.rules[element['geometry']][0])),
                ('reconstruction_dense_inner', np.linspace(-1., 1., 2*degree+7))):
                yield name, element['geometry'], tensor_points(element, self.fields[0].axes, [nodes, nodes]), True
        for chart, _ in _horizontal(self.fields, self.d, self.H, 8):
            yield from self.grouped('interface_inner', [chart])

    def sets(self):
        static = list(self.fixed_sets())
        anchors = {float(point[0]) for _, _, points, _ in static for point in points}
        yield from static
        for chart, _ in _horizontal(self.fields, self.d, self.H, int(self.raw['quadrature_orders'][0])):
            anchors.add(float(chart[0]))
            ray = self.ray(chart)
            yield from self.grouped('source_volume_inner', ray['points'])
            yield from self.ray_sets('normalization_volume_current', chart)
            # Prescribed-load quadrature was frozen at construction using
            # the neutral K0 and baseline profile before transition readiness.
            seed = self.ray(chart, seed=True)
            yield from self.grouped('prescribed_load_frozen_seed', seed['points'])
        for chart in self.trace_charts():
            anchors.add(float(chart[0]))
            for baseline in (True, False):
                yield from self.ray_sets('normalization_trace_'+('baseline' if baseline else 'current'),
                    chart, baseline=baseline)
        for x in sorted(anchors):
            yield from self.grouped('complete_physical_matching_transfer', [[x, 0., self.H]], inner=False)
            yield from self.grouped('complete_physical_sample_K', [[x, 0., self.raw['parameters'][3]]])
            for index, distance in enumerate(self.raw['transition_controls'][2:4]):
                yield from self.grouped('complete_physical_sensor_transfer_'+str(index),
                    [[x, 0., distance]], inner=False)
        # The strict constitutive audit also evaluates dependent epsilon at
        # native/dense inner locations. Cover those complete normalization
        # rays as an explicit diagnostic superset, not old-state anchors.
        anchors = set()
        for element in self.fields[self.d+2].elements:
            degree = self.fields[0].by_id[element['geometry']]['degree']
            for nodes in (gll(self.rules[element['geometry']][0]), np.linspace(-1., 1., 2*degree+7)):
                for value in nodes:
                    anchors.add(element['lower'][0]+(value+1)*(element['upper'][0]-element['lower'][0])/2)
        for x in sorted(anchors):
            yield from self.ray_sets('normalization_constitutive_reconstruction_current', [x, 0., self.H])


def finite_ranges(raw, fields, state, session, inventory=None):
    inventory = Inventory(raw, fields, state, session) if inventory is None else inventory
    us, qstar, kappa = raw['ustar'], raw['parameters'][6], raw['parameters'][7]
    K0 = qstar*us*us
    frame = np.column_stack([raw['robust_reference_flow'],
        np.cross(raw['robust_reference_normal'], raw['robust_reference_flow']), raw['robust_reference_normal']])
    center = np.zeros((3, 3)); center[0, 2] = 1.
    result = dict(sets={}, covariance_inner_minimum=math.inf, covariance_outer_minimum=math.inf,
        inner_minimum_k=math.inf, outer_minimum_k=math.inf, outer_minimum_epsilon=math.inf,
        minimum_sampled_K=math.inf, outer_gradient_entry_departure=0.,
        outer_k_over_K0=[math.inf, -math.inf], outer_epsilon_over_epsilon0=[math.inf, -math.inf],
        sampled_K_over_K0=[math.inf, -math.inf], inner_a_over_c=[math.inf, -math.inf])
    for name, identity, points, inner in inventory.sets():
        row = result['sets'].setdefault(name, dict(count=0, inner=inner,
            minimum_k=math.inf, covariance_minimum=math.inf))
        row['count'] += len(points)
        if inner:
            G, k, K, nt, _ = linear_inner_values(raw, fields, state, identity, points)
            if not np.all(np.isfinite(k)) or np.any(k < 0) or np.any(K <= 0):
                raise ValueError('R2 independent inventory outside inner scalar domain: '+name)
            selected = selected_values(raw, G, k, nt)
            R = np.asarray([real(value['stress']) for value in selected])
            row['minimum_sampled_K'] = min(row.get('minimum_sampled_K', math.inf), float(np.min(K)))
            row['minimum_k_over_K'] = min(row.get('minimum_k_over_K', math.inf), float(np.min(k/K)))
            result['minimum_sampled_K'] = min(result['minimum_sampled_K'], float(np.min(K)))
            kr = K/K0
            result['sampled_K_over_K0'] = [min(result['sampled_K_over_K0'][0], float(np.min(kr))),
                max(result['sampled_K_over_K0'][1], float(np.max(kr)))]
            result['inner_minimum_k'] = min(result['inner_minimum_k'], float(np.min(k)))
            ratios = np.asarray([float(real(value['a']/value['c'])) for value in selected if value['c'] != 0])
            row['saturation_active_count'] = row.get('saturation_active_count', 0)+sum(value['active'] for value in selected)
            row['zero_k_count'] = row.get('zero_k_count', 0)+int(np.count_nonzero(k == 0))
            if len(ratios):
                previous = row.get('a_over_c', [math.inf, -math.inf])
                row['a_over_c'] = [min(previous[0], float(np.min(ratios))), max(previous[1], float(np.max(ratios)))]
                result['inner_a_over_c'] = [min(result['inner_a_over_c'][0], float(np.min(ratios))),
                    max(result['inner_a_over_c'][1], float(np.max(ratios)))]
            eigen = eigen_example(R, points, G, k)
            index = int(np.argmin(np.linalg.eigvalsh(R)[:, 0]))
            largest = float(np.linalg.eigvalsh(G[index]+G[index].T)[-1])
            eigen.update(sampled_K=float(K[index]), k_over_sampled_K=float(k[index]/K[index]),
                nu_A=float(nt[index]), strain_largest_eigenvalue=largest,
                r_R=float(3*nt[index]*largest/(2*k[index])) if k[index] > 0 else None)
        else:
            G, k, epsilon, distance = local_values(raw, fields, state, identity, points)
            if not np.all(np.isfinite(k+epsilon)) or np.any(k <= 0) or np.any(epsilon <= 0):
                raise ValueError('R2 independent inventory outside outer scalar domain: '+name)
            normalized = np.einsum('ai,qab,bj->qij', frame, G, frame)/(us/(kappa*distance))[:, None, None]
            departure = np.max(np.abs(normalized-center), axis=(1, 2))
            kr, er = k/K0, epsilon/(us**3/(kappa*distance))
            strain, rotation = G+G.transpose(0, 2, 1), G-G.transpose(0, 2, 1)
            qs = k/epsilon*np.sqrt(np.sum(strain*strain, axis=(1, 2))/2)
            qo = k/epsilon*np.sqrt(np.sum(rotation*rotation, axis=(1, 2))/2)
            row['minimum_epsilon'] = min(row.get('minimum_epsilon', math.inf), float(np.min(epsilon)))
            for key, values in (('k_over_K0', kr), ('epsilon_over_epsilon0', er),
                                ('strain_invariant', qs), ('rotation_invariant', qo)):
                previous = row.get(key, [math.inf, -math.inf])
                row[key] = [min(previous[0], float(np.min(values))), max(previous[1], float(np.max(values)))]
            row['maximum_gradient_entry_departure'] = max(row.get('maximum_gradient_entry_departure', 0.), float(np.max(departure)))
            index = int(np.argmax(departure))
            if float(departure[index]) >= result['outer_gradient_entry_departure']:
                result['outer_gradient_entry_departure'] = float(departure[index])
                result['outer_gradient_worst'] = dict(set=name, geometry=identity,
                    chart=points[index].tolist(), gradient=G[index].tolist(),
                    normalized_gradient=normalized[index].tolist(), k_over_K0=float(kr[index]),
                    epsilon_over_epsilon0=float(er[index]))
            for key, values in (('outer_k_over_K0', kr), ('outer_epsilon_over_epsilon0', er)):
                result[key] = [min(result[key][0], float(np.min(values))), max(result[key][1], float(np.max(values)))]
            result['outer_minimum_k'] = min(result['outer_minimum_k'], float(np.min(k)))
            result['outer_minimum_epsilon'] = min(result['outer_minimum_epsilon'], float(np.min(epsilon)))
            R = np.asarray([selected_stress(raw, g, kval, e, p)[0] for g, kval, e, p in zip(G, k, epsilon, points)])
            eigen = eigen_example(R, points, G, k)
            for key, margin in (
                ('gradient_box_margin', .001-float(np.max(departure))),
                ('k_box_margin', min(float(np.min(kr))-.99, 1.01-float(np.max(kr)))),
                ('epsilon_box_margin', min(float(np.min(er))-.98, 1.02-float(np.max(er))))):
                row[key] = min(row.get(key, math.inf), margin)
        row['minimum_k'] = min(row['minimum_k'], float(np.min(k)))
        if eigen['covariance_minimum'] < row['covariance_minimum']:
            row.update(covariance_minimum=eigen['covariance_minimum'], worst_covariance=dict(eigen, geometry=identity))
        key = 'covariance_inner_minimum' if inner else 'covariance_outer_minimum'
        result[key] = min(result[key], eigen['covariance_minimum'])
    result['normalization_rays'] = inventory.ray_summary
    result['minimum_lambda'] = min(row['lambda_range'][0] for row in inventory.ray_summary.values())
    result['fixed_outer_box_margins'] = dict(
        gradient=.001-result['outer_gradient_entry_departure'],
        k=min(result['outer_k_over_K0'][0]-.99, 1.01-result['outer_k_over_K0'][1]),
        epsilon=min(result['outer_epsilon_over_epsilon0'][0]-.98, 1.02-result['outer_epsilon_over_epsilon0'][1]))
    result['finite_scalar_admissible'] = (result['inner_minimum_k'] >= 0 and
        min(result['outer_minimum_k'], result['outer_minimum_epsilon'], result['minimum_sampled_K']) > 0)
    result['finite_covariance_nonnegative'] = min(result['covariance_inner_minimum'], result['covariance_outer_minimum']) >= 0
    result['fixed_outer_box_covered'] = min(result['fixed_outer_box_margins'].values()) >= 0
    result['lambda_guard_satisfied'] = result['minimum_lambda'] >= .5
    result['finite_joint_preflight_passed'] = (result['finite_scalar_admissible'] and
        result['finite_covariance_nonnegative'] and result['fixed_outer_box_covered'] and result['lambda_guard_satisfied'])
    result['scope'] = 'Explicit finite source/mass/body/sampler/interface/transfer/GLL/dense sets, full current volume rays, full baseline/current trace rays and diagnostic reconstruction rays; one-sided element gradients retained. Not a polynomial-everywhere or forward-invariance proof.'
    return result


def trace_audit(raw, fields, state, inventory):
    """Independently project baseline/current actual capped EPS-H values."""
    d, offsets, H = raw['dimension'], raw['offsets'], raw['parameters'][1]
    baseline, correction, baseline_ids, samples = {}, {}, set(), []
    for element in fields[d+1].elements:
        if abs(element['lower'][2]-H) > 1e-12:
            continue
        nodes, transform = _segment_projection(element['degree'])
        x = element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
        data = []
        for value in x:
            point = np.array([value, 0., H])
            base, current = inventory.ray(point, baseline=True), inventory.ray(point)
            base_E = float(base['lam']*np.asarray(base['profile'].fields([H]), float)[2, 0])
            current_E = float(current['lam']*np.asarray(current['profile'].fields([H]), float)[2, 0])
            actual = fields[d+1].evaluate(state, offsets[d+1], point, 1)[0]
            E0 = raw['ustar']**3/(raw['parameters'][7]*H)
            data.append([base_E, current_E-base_E])
            samples.append(dict(chart=point.tolist(), chi=float(current['profile'].chi),
                baseline_lambda=float(base['lam']), current_lambda=float(current['lam']),
                baseline_target=base_E, current_target=current_E, actual_epsilon=actual,
                actual_epsilon_over_epsilon0=actual/E0,
                current_target_over_epsilon0=current_E/E0))
        maps = [transform@np.asarray(data).T[part] for part in range(2)]
        count = element['degree']+1
        for table, values in zip((baseline, correction), maps):
            for index, sign, value in zip(element['ids'][:count], element['sign'][:count], values):
                index, value = int(index), float(sign*value)
                if index in table and abs(table[index]-value) > 2e-10:
                    raise ValueError('Inconsistent independently shared R2 trace projection')
                table[index] = value
        if coordinate(dict(raw, transition=dict(controls=raw['transition_controls'])), element['lower']) >= raw['transition_controls'][4]-1e-11:
            baseline_ids.update(element['ids'][:count])
    for index in baseline_ids:
        correction[index] = 0.
    targets = {index: baseline[index]+correction[index] for index in baseline}
    differences = [abs(state[offsets[d+1]+index]-value) for index, value in targets.items()]
    zero = [row['actual_epsilon_over_epsilon0'] for row in samples if row['chi'] == 0.]
    return dict(maximum_trace_coefficient_error=max(differences, default=0.),
        relative_to_E0=max(differences, default=0.)/(raw['ustar']**3/(raw['parameters'][7]*H)),
        matching_rows=len(targets), projected_target_coefficients=[[index, value] for index, value in sorted(targets.items())],
        chi_zero_epsilon_ratio_range=[min(zero), max(zero)] if zero else None,
        actual_epsilon_ratio_range=[min(row['actual_epsilon_over_epsilon0'] for row in samples),
            max(row['actual_epsilon_over_epsilon0'] for row in samples)],
        baseline_and_current_lambda_range=[min(min(row['baseline_lambda'], row['current_lambda']) for row in samples),
            max(max(row['baseline_lambda'], row['current_lambda']) for row in samples)],
        samples=samples,
        scope='Independent constrained native projection of lambda*Epre(H), including chi=0; no source-cap plateau assumption')


def native_domain_audit(raw, fields, state, domain, trace, audit, inventory):
    """Verify every exported actual scalar observation from independent rows.

    The native inventory supplies coordinates and declared field/side identity;
    it supplies no value used to construct the independent answer. Independently
    rebuilt full rays are checked separately, including their count and lambda.
    """
    if domain is None:
        return dict(available=False, complete_crosscheck=False,
            reason='Native ScalarDomain observations were not supplied; independent finite inventory retained')
    observations = domain.get('observations', [])
    audit.require(bool(observations), 'Missing complete native ScalarDomain observations')
    effective = np.asarray(domain['effective_state'], float)
    target = effective.copy()
    offsets, d = raw['offsets'], raw['dimension']
    for index, value in trace['projected_target_coefficients']:
        target[offsets[d+1]+index] = value
    native_target = np.asarray(domain['trace_target_state'], float)
    audit.near(effective, state, 'R2_initial_external_lift_identity', 2e-10)
    trace_ids = np.asarray([offsets[d+1]+int(row[0]) for row in trace['projected_target_coefficients']], int)
    audit.near(native_target[trace_ids], target[trace_ids], 'R2_native_independent_projected_trace_state', 2e-10)
    groups, counts, margins, projected_ratios = {}, {}, [], []
    for row in observations:
        counts[row['inventory']] = counts.get(row['inventory'], 0)+1
        groups.setdefault((int(row['block']), int(row['geometry']), bool(row['dependent_trace'])), []).append(row)
    maximum_value_error, projected_covariance = 0., math.inf
    for (block, identity, dependent), rows in groups.items():
        points = np.asarray([row['coordinate'] for row in rows], float)
        coefficients = target if dependent else effective
        values = physical_rows(fields[block], coefficients, offsets[block], points, identity)[3]
        expected_scale = np.full(len(rows), raw['parameters'][5]) if block == d else (
            raw['ustar']**3/(raw['parameters'][7]*points[:, 2]))
        native_values, native_scales = np.asarray([row['value'] for row in rows]), np.asarray([row['scale'] for row in rows])
        audit.near(native_values, values, 'R2_native_scalar_observation_values', 2e-10)
        audit.near(native_scales, expected_scale, 'R2_native_scalar_observation_scales', 2e-12)
        maximum_value_error = max(maximum_value_error, float(np.max(np.abs(values-native_values))))
        margins.extend((values/expected_scale).tolist())
        if dependent:
            projected_ratios.extend((values/expected_scale).tolist())
            G = _velocity(fields, effective, offsets, d, points, identity)[1]
            k = physical_rows(fields[d], effective, offsets[d], points, identity)[3]
            tensors = np.asarray([selected_stress(raw, g, kval, epsilon, point)[0]
                for g, kval, epsilon, point in zip(G, k, values, points)])
            projected_covariance = min(projected_covariance, float(np.min(np.linalg.eigvalsh(tensors)[:, 0])))
        audit.require(all(value > 0 if row['strict'] else value >= 0 for value, row in zip(values, rows)),
            'An exported native scalar observation is inadmissible')
    audit.require(counts == domain['inventory_counts'], 'Native ScalarDomain count/list mismatch')
    audit.require(domain['admissible'] and domain['phase'] == 'complete', 'Native initial ScalarDomain did not complete')
    audit.near(domain['minimum_dimensionless_margin'], min(margins), 'R2_native_scalar_minimum_margin', 2e-10)
    ray_comparison = {}
    for name in ('volume_normalization_ray', 'trace_baseline_normalization_ray', 'trace_current_normalization_ray'):
        actual = [row for row in observations if row['inventory'] == name]
        expected_points, expected_weights = [], []
        anchors = (_horizontal(fields, d, raw['parameters'][1], int(raw['quadrature_orders'][0]))
            if name == 'volume_normalization_ray' else ((chart, 0.) for chart in inventory.trace_charts()))
        for chart, horizontal in anchors:
            ray = inventory.ray(chart, baseline=name == 'trace_baseline_normalization_ray')
            expected_points.extend(ray['points'])
            expected_weights.extend(horizontal*ray['metric']['baseArea']*ray['measure']
                if name == 'volume_normalization_ray' else np.zeros(len(ray['points'])))
        actual_points = np.asarray([row['coordinate'] for row in actual])
        actual_weights = np.asarray([row['weight'] for row in actual])
        expected_points, expected_weights = np.asarray(expected_points), np.asarray(expected_weights)
        ai = np.lexsort((actual_points[:, 2], actual_points[:, 1], actual_points[:, 0]))
        ei = np.lexsort((expected_points[:, 2], expected_points[:, 1], expected_points[:, 0]))
        audit.near(actual_points[ai], expected_points[ei], 'R2_independent_dynamic_ray_coordinates_'+name, 2e-11)
        audit.near(actual_weights[ai], expected_weights[ei], 'R2_independent_dynamic_ray_weights_'+name, 2e-11)
        ray_comparison[name] = dict(count=len(actual),
            maximum_coordinate_error=float(np.max(np.abs(actual_points[ai]-expected_points[ei]))),
            maximum_weight_error=float(np.max(np.abs(actual_weights[ai]-expected_weights[ei]))))
    return dict(available=True, complete_crosscheck=True, counts=counts,
        observations=len(observations), maximum_value_error=maximum_value_error,
        minimum_dimensionless_margin=min(margins), minimum_lambda=domain['minimum_lambda'],
        projected_trace_over_epsilon0=[min(projected_ratios), max(projected_ratios)],
        projected_trace_fixed_box_margin=min(min(projected_ratios)-.98, 1.02-max(projected_ratios)),
        projected_trace_covariance_minimum=projected_covariance,
        effective_state_float64_LE_sha256=vector_sha(effective),
        trace_target_state_float64_LE_sha256=vector_sha(native_target),
        independent_dynamic_ray_coordinates_and_weights=ray_comparison,
        scope='Every exported actual native scalar value independently reconstructed; dynamic ray construction is independently compared separately')


def original_derivative_witness(raw, fields, state):
    point = np.array([-3.2041499092834287, 0., 8/3])
    us, kappa = raw['ustar'], raw['parameters'][7]
    frame = np.column_stack([raw['robust_reference_flow'],
        np.cross(raw['robust_reference_normal'], raw['robust_reference_flow']), raw['robust_reference_normal']])
    center = np.zeros((3, 3)); center[0, 2] = 1.
    result = []
    for element in fields[0].elements:
        if element['lower'][2] < raw['parameters'][1] or not (
                np.all(point >= np.asarray(element['lower'])-1e-12) and np.all(point <= np.asarray(element['upper'])+1e-12)):
            continue
        G, k, epsilon, _ = local_values(raw, fields, state, element['geometry'], [point])
        normalized = frame.T@G[0]@frame/(us/(kappa*point[2]))
        result.append(dict(geometry=element['geometry'], chart=point.tolist(),
            normal_cell=[element['lower'][2], element['upper'][2]],
            gradient=G[0].tolist(), normalized_gradient=normalized.tolist(),
            gradient_entry_departure=float(np.max(np.abs(normalized-center))),
            k=float(k[0]), epsilon=float(epsilon[0])))
    if len(result) != 2:
        raise ValueError('Frozen r=8/3 interface witness did not retain both element sides')
    return result


def verify(raw, session, native_domain=None, ranges_only=False):
    audit = Audit()
    audit.require(raw['dimension'] == 2 and not raw['mapped'] and raw['upstream_transition'],
        'R2 preflight requires the frozen flat 2D upstream control')
    audit.require(raw['eta'] == 0 and raw['dt'] == .0025 and raw['steps'] == 0 and raw['time'] == 0,
        'R2 preflight is zero-step with target startup dt=.0025 and eta=0')
    audit.require(raw['robust_parameters'] == [.4, 0., .2, .1, 0., 1., 1.], 'Changed frozen R2 model/source coefficients')
    info = {node.get('PROPERTY'): node.get('VALUE') for node in ET.parse(session).getroot().findall('CONDITIONS/SOLVERINFO/I')}
    audit.require(raw.get('inner_stress_model') == info['RANSWallInnerStress'], 'Missing/mismatched explicit inner-stress identity')
    audit.require(raw.get('nonlinear_policy') == info['RANSWallNonlinearPolicy'], 'Missing/mismatched nonlinear-policy identity')
    states, times = raw.get('accepted_states', []), raw.get('audit_times', [])
    audit.require(len(states) == len(times) == 1 and times == [0], 'Preflight has missing or evolved native history')
    state, fields = np.asarray(states[0], float), _fields(raw)
    audit.require(len(state) == raw['offsets'][-1]+1 and np.all(np.isfinite(state)), 'Invalid full initial coefficient vector')
    inventory = Inventory(raw, fields, state, session)
    ranges = finite_ranges(raw, fields, state, session, inventory)
    trace = trace_audit(raw, fields, state, inventory)
    audit.require(trace['relative_to_E0'] <= 1e-9, 'Initial dependent epsilon trace is inconsistent')
    direct = selected_direct(raw, fields, state, audit)
    native = raw['state_audits'][0]
    audit.near(native['robust_source'], direct['robust_source'], 'R2_independent_capped_source_budgets')
    audit.near(native['wall_sources'], direct['wall_sources'], 'R2_independent_actual_inner_production_and_molecular_work')
    audit.near(native['source_volume_inner'][3], direct['inner_covariance_minimum'], 'R2_independent_selected_inner_source_covariance')
    audit.near(native['source_volume_outer'][3], direct['outer_covariance_minimum'], 'R2_independent_outer_source_covariance')
    if raw['inner_stress_model'] == 'SaturatedR2':
        audit.near(native['inner_stress'], direct['inner_stress'], 'R2_independent_saturation_ranges_and_physical_volumes')
        audit.near(native['inner_constitutive_modification'], direct['inner_constitutive_modification'],
            'R2_selected_inner_mean_work_and_production_exchange')
    traction = selected_interface(raw, fields, state)
    audit.near(native['interface'], traction, 'R2_independent_all_component_total_traction_and_scalar_flux')
    domain_result = native_domain_audit(raw, fields, state, native_domain, trace, audit, inventory)
    if domain_result['complete_crosscheck']:
        aliases = {'volume_normalization_ray': 'normalization_volume_current',
            'trace_baseline_normalization_ray': 'normalization_trace_baseline',
            'trace_current_normalization_ray': 'normalization_trace_current'}
        for native_name, independent_name in aliases.items():
            audit.require(domain_result['counts'][native_name] == ranges['normalization_rays'][independent_name]['points'],
                'Actual/independent full normalization ray count differs: '+native_name)
        expected_lambda = min(ranges['normalization_rays'][name]['lambda_range'][0] for name in aliases.values())
        audit.near(domain_result['minimum_lambda'], expected_lambda, 'R2_independent_actual_volume_and_trace_minimum_lambda', 2e-10)
    result = dict(schema=2, state_selection='Actual sole stored t=0 accepted state; no solver call by this audit',
        inner_stress_model=raw['inner_stress_model'], nonlinear_policy=raw['nonlinear_policy'],
        coefficient_float64_LE_sha256=vector_sha(state), field_sizes=[field.size for field in fields], unknowns=len(state),
        maximum_coefficient_difference_from_fixed_reference=float(np.max(np.abs(state-np.asarray(raw['fixed_reference'])))),
        finite_ranges=ranges, independent_trace=trace, independent_source_energy=direct,
        independent_total_traction=traction.tolist(), native_scalar_domain=domain_result,
        old_r8_over3_one_sided_derivative_witness=original_derivative_witness(raw, fields, state),
        finite_observed_joint_preflight=(ranges['finite_joint_preflight_passed'] and
            domain_result.get('projected_trace_fixed_box_margin', -math.inf) >= 0 and
            domain_result.get('projected_trace_covariance_minimum', -math.inf) >= 0),
        stricter_implementation_sampled_K_box_margin=min(ranges['sampled_K_over_K0'][0]-.99,
            1.01-ranges['sampled_K_over_K0'][1]),
        native_inventory_independently_crosschecked=domain_result['complete_crosscheck'],
        rigorous_discrete_or_trajectory_qualification=False,
        scope='Bounded finite representation/native implementation audit; no source/model/certificate/field replacement')
    if not ranges_only:
        loads, detail = weak_loads(raw, fields, state)
        norms = dual_norms(raw, fields, loads, state)
        reference = np.asarray(raw['fixed_reference'], float)
        reference_loads = loads if np.array_equal(reference, state) else weak_loads(raw, fields, reference)[0]
        audit.near(raw['raw_reference_loads'], reference_loads, 'R2_independent_complete_unanchored_fixed_reference_weak_loads', 2e-10)
        audit.near(raw['raw_reference_norms'], dual_norms(raw, fields, reference_loads, reference),
            'R2_independent_native_fixed_reference_mass_dual_block_norms', 2e-9)
        result['raw_unanchored_residual'] = dict(actual_initial_block_norms=norms,
            native_fixed_reference_block_norms=raw['raw_reference_norms'],
            maximum_native_fixed_reference_load_error=float(np.max(np.abs(reference_loads-np.asarray(raw['raw_reference_loads'])))),
            actual_initial_weak_load_float64_LE_sha256=vector_sha(loads), quadrature=detail,
            block_names=['u', 'v', 'k', 'outer_epsilon', 'inner_pressure_divergence', 'outer_pressure_divergence', 'gauge', 'constraints'],
            interpretation='Unanchored alpha=0 loads; actual initial state and frozen reference distinguished explicitly')
        result['same_state_quadrature_sensitivity'] = []
        native_orders = [int(raw['quadrature_orders'][0]), int(raw['quadrature_orders'][1]), int(raw['quadrature_orders'][0])]
        for orders in ((8, 16, 8), (8, 24, 8), (12, 24, 12)):
            if list(orders) == native_orders:
                continue
            other, _ = weak_loads(raw, fields, state, *orders)
            result['same_state_quadrature_sensitivity'].append(dict(orders=list(orders),
                maximum_weak_load_difference=float(np.max(np.abs(other-loads))),
                block_norms=dual_norms(raw, fields, other, state),
                difference_block_norms=dual_norms(raw, fields, other-loads, state),
                identical_coefficient_state=True, native_solve=False))
    result.update(implementation_audit_passed=True, independent_checks=audit.checks, errors=audit.errors)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--expected-source', required=True)
    parser.add_argument('--native-domain', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ranges-only', action='store_true')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Scientific numerical reconstruction belongs on authoritative Linux')
    if args.output.exists():
        raise FileExistsError('Preserve previous R2 audit: '+str(args.output))
    inputs = {name: dict(path=str(value.resolve()), sha256=sha(value)) for name, value in
        (('raw', args.raw), ('session', args.session), ('plan', args.plan))}
    if args.native_domain:
        inputs['native_domain'] = dict(path=str(args.native_domain.resolve()), sha256=sha(args.native_domain))
    receipt = dict(schema=2, milestone='M8-R2', role='E independent representation preflight',
        started_utc=datetime.now(timezone.utc).isoformat(), inputs=inputs,
        driver=dict(path=str(Path(__file__).resolve()), sha256=sha(__file__)),
        numerical_threads={key: os.environ[key] for key in THREADS},
        affinity=sorted(os.sched_getaffinity(0)), solver_invoked=False, imported_oracles={})
    try:
        source = args.source.resolve()
        head = subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
        status = subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), 'status', '--porcelain'], text=True).strip()
        if head != args.expected_source or status:
            raise ValueError('Oracle source differs from the explicitly frozen clean revision')
        receipt['source'] = dict(path=str(source), commit=head, clean=True)
        for name, module in sorted(sys.modules.items()):
            if not name.startswith(('cls_epawf_', 'validate_cls_epawf_')) or not getattr(module, '__file__', None):
                continue
            path = Path(module.__file__).resolve()
            if not path.is_relative_to(source):
                raise ValueError('Imported scientific helper outside the frozen source: '+str(path))
            receipt['imported_oracles'][name] = dict(path=str(path), sha256=sha(path))
        raw, plan = json.loads(args.raw.read_text()), json.loads(args.plan.read_text())
        if normalized_physics(ET.parse(args.session).getroot()) != plan['protected_physics']:
            raise ValueError('Protected physical inputs/functions changed')
        rows = [row for row in plan['cases'] if row['xml_sha256'] == inputs['session']['sha256']]
        if len(rows) != 1:
            raise ValueError('Exact session is not unique in the frozen case plan')
        expected = cost(rows[0]['representation'])
        if [field['size'] for field in raw['fields']] != expected['field_sizes'] or len(raw['state']) != expected['unknowns']:
            raise ValueError('Actual native unknown/field count differs from the static representation contract')
        radial = sorted({value for element in raw['fields'][0]['mesh'] for value in (element['lower'][2], element['upper'][2])})
        if radial != RADIAL[rows[0]['representation']] or list(raw['quadrature_orders']) != rows[0]['source_orders']:
            raise ValueError('Actual radial mesh or source quadrature differs from frozen plan')
        native_domain = json.loads(args.native_domain.read_text()) if args.native_domain else raw.get('initial_scalar_domain')
        receipt['cost'] = expected
        receipt['result'] = verify(raw, args.session, native_domain, args.ranges_only)
        receipt['status'] = 'passed'
    except Exception as error:
        receipt.update(status='audit_failed', failure=dict(type=type(error).__name__, message=str(error)))
    receipt['input_hashes_unchanged'] = all(sha(value['path']) == value['sha256'] for value in inputs.values())
    receipt['finished_utc'] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=receipt['status'], output=str(args.output),
        observed_joint_preflight=receipt.get('result', {}).get('finite_observed_joint_preflight'),
        failure=receipt.get('failure'))), flush=True)
    return 0 if receipt['status'] == 'passed' and receipt['input_hashes_unchanged'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
