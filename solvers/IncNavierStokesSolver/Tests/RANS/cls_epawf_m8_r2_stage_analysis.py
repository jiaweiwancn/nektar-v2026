#!/usr/bin/env python3
"""Read-only independent M8-R2 matched-stage and Newton-path analysis.

This program never calls a native executable. It reconstructs the unchanged
unforced equations, scalar evaluation domain, projected epsilon trace, physical
mass norms and matrix actions from retained native coefficients. A collected
nonlinear failure, an inaccurate action and a physical-range failure remain
different results. Numerical invocation belongs to the Linux environment.
"""
import argparse
from collections import Counter, defaultdict, OrderedDict
from datetime import datetime, timezone
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
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

from cls_epawf_m8_r1_nonlinear import Audit
from cls_epawf_m8_r2_preflight import (Inventory, finite_ranges, selected_direct,
    selected_interface, trace_audit, weak_loads)
from cls_epawf_m7_energy_oracle import _horizontal, _mass, _profile, _velocity, _work
from cls_epawf_m8_state_coverage import gll, local_values, tensor_points
from cls_epawf_m8_r2_oracle import selected_inner, selected_outer as selected_stress
from cls_epawf_m8_native_oracle import normal_rule
from validate_cls_epawf_m8_r2_pressure import PressureOracle

ACTION_TOLERANCE = 1e-7
LEVELS = [3200., 1600., 800., 400.]
LADDER = [2e-5, 2e-6, 5e-7]


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        allow_nan=False).encode()).hexdigest()


def vector_hash(vector):
    return hashlib.sha256(np.asarray(vector, dtype='<f8').tobytes()).hexdigest()


def rle(values):
    result = []
    for value in values:
        if result and result[-1][0] == value:
            result[-1][1] += 1
        else:
            result.append([value, 1])
    return result


def branch_difference(base, endpoint):
    result = {}
    for name in ('outer_max', 'inner_saturation', 'source_limiter'):
        a, b = base.get(name), endpoint.get(name)
        if a is None or b is None:
            result[name] = dict(available=False)
            continue
        result[name] = dict(available=True, identical=a['pattern_sha256'] == b['pattern_sha256'],
            same_length=a['count'] == b['count'], base_minimum_switch_distance=a['minimum_switch_distance'],
            endpoint_minimum_switch_distance=b['minimum_switch_distance'],
            interpretation='Fixed outer points' if name == 'outer_max' else
                'Literal actual source-ray ordering; moving normal coordinates need not represent identical physical points')
    return result


def save(path, payload):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def git(source, *arguments):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), *arguments], text=True).strip()


class Analyzer:
    def __init__(self, raw, session, source, matrix_root=None):
        self.raw, self.session, self.source = raw, Path(session), Path(source)
        self.audit = Audit(raw, self.source)
        self.audit.gauge_variable_scale = raw['ustar']/raw['parameters'][1]
        self.fields, self.offsets = self.audit.fields, raw['offsets']
        self.constants = [field.constant() for field in self.fields]
        self.row_cache, self.row_cache_bytes = OrderedDict(), 0
        self.d, self.H, self.K0, self.E0 = self.audit.d, self.audit.H, self.audit.K0, self.audit.E0
        self.size = len(raw['state'])
        self.trace_ids = self.fields[self.d+1].boundary(
            lambda axis, side, element: axis == 2 and side == 0 and abs(element['lower'][2]-self.H) < 1e-11)
        self.trace_rows = {self.offsets[self.d+1]+int(index) for index in self.trace_ids}
        initial_constraints = self.audit.constraints(raw['reference_constraints'])
        self.external = {row: value for row, value in initial_constraints.items() if row not in self.trace_rows}
        self.static = self.static_locations()
        self.bind_actual_static_locations()
        self.matrix_root = Path(matrix_root).resolve() if matrix_root else None
        self.state_cache, self.matrix_files, self.gates = {}, {}, []
        self._continuity, self._continuity_lifted, self._continuity_factors = None, None, None
        self.continuity_gates = []
        self.failures = []

    def attempt(self, label, function):
        try:
            return function()
        except Exception as error:
            record = dict(label=label, type=type(error).__name__, message=str(error), traceback=traceback.format_exc())
            self.failures.append(record)
            return dict(analysis_failed=True, **record)

    def external_lift(self, state):
        result = np.asarray(state, float).copy()
        self.audit.require(result.shape == (self.size,) and np.all(np.isfinite(result)), 'Invalid full coupled state')
        for row, value in self.external.items():
            result[row] = value
        return result

    def values(self, state, block, identity, points):
        """Centered long-double values, with independent Jacobi basis rows."""
        field = self.fields[block]
        points = np.atleast_2d(points)
        key = (int(block), int(identity), points.tobytes())
        if key not in self.row_cache:
            ids, phi, _ = field.rows(field.by_id[int(identity)], points)
            self.row_cache[key] = (ids, phi)
            self.row_cache_bytes += len(key[2])+ids.nbytes+phi.nbytes
            while self.row_cache_bytes > 128*1024*1024 and len(self.row_cache) > 1:
                old_key, (old_ids, old_phi) = self.row_cache.popitem(last=False)
                self.row_cache_bytes -= len(old_key[2])+old_ids.nbytes+old_phi.nbytes
        else:
            self.row_cache.move_to_end(key)
        ids, phi = self.row_cache[key]
        coefficients = np.asarray(state)[self.offsets[block]+ids].astype(np.longdouble)
        constant = self.constants[block][ids].astype(np.longdouble)
        pivot = next((i for i, value in enumerate(constant) if value != 0), None)
        self.audit.require(pivot is not None, 'Native tensor restriction lacks its constant')
        centre = coefficients[pivot]/constant[pivot]
        delta = coefficients-centre*constant
        value = np.asarray(centre+phi.astype(np.longdouble)@delta, float)
        bound = 64*np.finfo(float).eps*(abs(float(centre))+np.abs(phi)@np.abs(np.asarray(delta, float)))
        return value, np.asarray(bound, float)

    def grouped(self, name, points, inner, side=0, weight=None):
        groups = defaultdict(list)
        for index, point in enumerate(np.atleast_2d(points)):
            identity = self.fields[self.d].element(point, side)['geometry']
            groups[identity].append(index)
        points = np.atleast_2d(points)
        weights = np.zeros(len(points)) if weight is None else np.asarray(weight)
        for identity, indices in groups.items():
            indices = np.asarray(indices, int)
            yield dict(name=name, geometry=int(identity), points=points[indices], block=self.d,
                strict=not inner, side=side, dependent=False, weights=weights[indices])
            if not inner:
                yield dict(name=name, geometry=int(identity), points=points[indices], block=self.d+1,
                    strict=True, side=side, dependent=False, weights=weights[indices])

    def static_locations(self):
        """The common residual-domain set, not the wider constitutive audit."""
        locations, samplers, trace = [], set(), set()
        def add(name, points, inner, side=0):
            points = np.atleast_2d(points)
            locations.extend(self.grouped(name, points, inner, side))
            samplers.update((float(point[0]), float(point[1])) for point in points)
        order = int(self.raw['quadrature_orders'][0])
        seed = Inventory(self.raw, self.fields, self.external_lift(self.raw['initial']), self.session)
        for element in self.fields[self.d].elements:
            inner = element['upper'][2] <= self.H
            points, _ = self.fields[self.d].quadrature(element, max(8, element['degree']+3))
            add('mass_volume', points, inner)
            if not inner:
                points, _ = self.fields[self.d].quadrature(element, order)
                add('source_volume_outer', points, False)
            for kind, nodes in (('reconstruction_native', [gll(n) for n in seed.rules[element['geometry']]]),
                    ('reconstruction_dense', [np.linspace(-1., 1., 2*element['degree']+(7 if inner else 3))]*2)):
                points = tensor_points(element, (0, 2), nodes)
                # Preserve each element's one-sided rows at internal seams.
                for point in points:
                    side = (-1 if inner else 1) if abs(point[2]-self.H) < 1e-11 else 0
                    locations.append(dict(name=kind, geometry=int(element['geometry']), points=np.asarray([point]),
                        block=self.d, strict=not inner, side=side, dependent=False, weights=np.zeros(1)))
                    if not inner:
                        epsilon_geometry = self.fields[self.d+1].element(point,
                            1 if abs(point[2]-self.H) < 1e-11 else 0)['geometry']
                        locations.append(dict(name=kind, geometry=int(epsilon_geometry), points=np.asarray([point]),
                            block=self.d+1, strict=True, side=side, dependent=False, weights=np.zeros(1)))
                    samplers.add((float(point[0]), float(point[1])))
        for chart, _ in _horizontal(self.fields, self.d, self.H, order):
            add('prescribed_source_inner', seed.ray(chart, seed=True)['points'], True)
            samplers.add((float(chart[0]), float(chart[1])))
        for chart, _ in _horizontal(self.fields, self.d, self.H, 8):
            add('interface_inner', [chart], True, -1)
            add('interface_outer', [chart], False, 1)
        for element in self.fields[self.d+1].elements:
            if abs(element['lower'][2]-self.H) > 1e-11:
                continue
            for node in gll(max(8, element['degree']+3)):
                chart = np.array([element['lower'][0]+(node+1)*(element['upper'][0]-element['lower'][0])/2, 0., self.H])
                add('matching_trace_nodes', [chart], False, 1)
                trace.add(tuple(chart))
        # Roundoff-identical horizontal coordinates should not become extra
        # mathematical anchors merely because two independent maps differ by
        # a last bit. Exported counts and coordinate multisets remain checked.
        unique = {}
        for x, y in sorted(samplers):
            unique.setdefault((round(x, 13), round(y, 13)), (x, y))
        self.sampler_anchors = list(unique.values())
        self.trace_anchors = [np.asarray(point) for point in sorted(trace)]
        for x, y in self.sampler_anchors:
            add('physical_matching_transfer', [[x, y, self.H]], False, 1)
            point = np.array([x, y, self.raw['parameters'][3]])
            identity = self.fields[self.d].element(point, -1 if abs(point[2]-self.H) < 1e-11 else 0)['geometry']
            locations.append(dict(name='physical_sampled_K', geometry=int(identity), points=np.asarray([point]),
                block=self.d, strict=True, side=-1 if abs(point[2]-self.H) < 1e-11 else 0,
                dependent=False, weights=np.zeros(1)))
            for radius in self.raw['transition_controls'][2:4]:
                add('physical_sensor_transfer', [[x, y, radius]], False)
        return locations

    @staticmethod
    def point_key(row):
        return (row['inventory'], int(row['geometry']), int(row['block']), int(row['side']),
            bool(row['strict']), bool(row['dependent_trace']), *(round(float(x), 12) for x in row['coordinate']))

    def bind_actual_static_locations(self):
        """Check geometric coverage, then preserve every actual duplicate call.

        The compiled tensor rules and NumPy rules can differ in the last bit.
        A C++ set of binary64 sampler coordinates can consequently contain
        duplicates at one mathematical location. Coordinates carry no field
        answer; all values are reconstructed independently at every native
        call, including these literal duplicates.
        """
        observations = self.raw['initial_scalar_domain'].get('observations', [])
        self.audit.require(bool(observations), 'Missing full initial common-inventory observations')
        static_names = {group['name'] for group in self.static}
        actual = [row for row in observations if row['inventory'] in static_names]
        expected = []
        for group in self.static:
            for point in group['points']:
                expected.append(dict(inventory=group['name'], coordinate=point.tolist(), geometry=group['geometry'],
                    block=group['block'], side=group['side'], strict=group['strict'], dependent_trace=False))
        actual_keys, expected_keys = Counter(map(self.point_key, actual)), Counter(map(self.point_key, expected))
        self.audit.require(set(actual_keys) == set(expected_keys), 'Native/common independent static geometric coverage differs')
        self.static_coverage = dict(native_calls=len(actual), independent_rule_calls=len(expected),
            unique_descriptor_locations=len(actual_keys),
            native_literal_duplicate_calls=len(actual)-len(actual_keys),
            independent_rule_duplicate_calls=len(expected)-len(expected_keys),
            differing_literal_multiplicities=sum(actual_keys[key] != expected_keys[key] for key in actual_keys),
            coordinate_comparison_rounding_digits=12,
            interpretation='Rounding identifies geometric locations only; actual call coordinates and values are retained without rounding')
        groups = defaultdict(list)
        for row in actual:
            groups[(row['inventory'], int(row['geometry']), int(row['block']), int(row['side']), bool(row['strict']))].append(row)
        self.static = [dict(name=name, geometry=identity, points=np.asarray([row['coordinate'] for row in rows]),
            block=block, strict=strict, side=side, dependent=False, weights=np.zeros(len(rows)))
            for (name, identity, block, side, strict), rows in groups.items()]
        source = sorted({(float(row['coordinate'][0]), float(row['coordinate'][1])) for row in observations
            if row['inventory'] == 'volume_normalization_ray'})
        horizontal = list(_horizontal(self.fields, self.d, self.H, int(self.raw['quadrature_orders'][0])))
        self.audit.require(len(source) == len(horizontal), 'Native/independent source ground-ray count differs')
        self.source_anchors = []
        for x, y in source:
            matches = [(chart, weight) for chart, weight in horizontal
                if abs(chart[0]-x) < 1e-12 and abs(chart[1]-y) < 1e-12]
            self.audit.require(len(matches) == 1, 'Native source anchor is not an independent Gauss anchor')
            self.source_anchors.append((np.array([x, y, self.H]), float(matches[0][1])))
        self.trace_anchors = [np.asarray(point) for point in sorted({tuple(row['coordinate']) for row in actual
            if row['inventory'] == 'matching_trace_nodes'})]

    def projected(self, state):
        state = self.external_lift(state)
        inventory = Inventory(self.raw, self.fields, state, self.session)
        trace = trace_audit(self.raw, self.fields, state, inventory)
        constraints = dict(self.external)
        constraints.update({self.offsets[self.d+1]+int(index):float(value)
            for index, value in trace['projected_target_coefficients']})
        self.audit.require(set(constraints)-set(self.external) == self.trace_rows, 'Independent matching-row identity changed')
        target = state.copy()
        for row, value in constraints.items():
            target[row] = value
        return state, target, constraints, trace, inventory

    def domain(self, state, native=None, retain=False):
        """Reconstruct actual moving source/trace rays before calling a residual."""
        effective = self.external_lift(state)
        retain = retain or bool(native and native.get('observations'))
        records, counts, margins, minimum = [], Counter(), [], {}
        def observe(group, values):
            block, dependent = group['block'], group['dependent']
            points = group['points']
            numbers, arithmetic = self.values(values, block, group['geometry'], points)
            for index, (point, value, bound) in enumerate(zip(points, numbers, arithmetic)):
                scale = self.K0 if block == self.d else self.audit.us**3/(self.audit.kappa*point[2])
                strict = group['strict']; valid = bool(value > 0 if strict else value >= 0)
                name = group['name']; counts[name] += 1; margins.append(float(value/scale))
                row = dict(inventory=name, coordinate=point.tolist(), geometry=group['geometry'], block=block,
                    side=group['side'], value=float(value), scale=float(scale), weight=float(group['weights'][index]),
                    strict=bool(strict), dependent_trace=bool(dependent), admissible=valid,
                    arithmetic_roundoff_scale=float(bound), sign_resolved=bool(abs(value) > bound))
                if retain or not valid:
                    records.append(row)
                key = 'minimum_sampled_K' if name == 'physical_sampled_K' else (
                    'minimum_inner_k' if block == self.d and point[2] <= self.H and group['side'] <= 0 else
                    'minimum_outer_k' if block == self.d else
                    'minimum_projected_trace_epsilon' if dependent else 'minimum_outer_epsilon')
                minimum[key] = min(minimum.get(key, math.inf), float(value))
                if not valid:
                    return False
            return True
        for group in self.static:
            if not observe(group, effective):
                return self.domain_result(False, 'fixed_fields_and_samplers', records, counts, margins, minimum, native)
        inventory = Inventory(self.raw, self.fields, effective, self.session)
        rays, branch = [], []
        for chart, horizontal in self.source_anchors:
            rays.append(('volume_normalization_ray', chart, False, horizontal))
        for chart in self.trace_anchors:
            rays.extend([('trace_baseline_normalization_ray', chart, True, 0.),
                         ('trace_current_normalization_ray', chart, False, 0.)])
        normalized = []
        try:
            for name, chart, baseline, horizontal in rays:
                view = dict(self.raw, upstream_transition=False, chi=0.) if baseline else self.raw
                profile, metric = _profile(view, self.fields, effective, chart)
                radial, radial_weight = normal_rule(profile, self.raw, self.fields, metric)
                points = np.column_stack([np.full_like(radial, chart[0]), np.full_like(radial, chart[1]), radial])
                weights = horizontal*radial_weight*np.asarray(profile.jacobian(radial), float) if horizontal else np.zeros(len(points))
                for group in self.grouped(name, points, True, weight=weights):
                    if not observe(group, effective):
                        return self.domain_result(False, 'moving_ray_k_before_normalization', records, counts, margins, minimum, native)
                # Only after all transported k on this actual moving ray has
                # passed may the nonlocal capped normalization be evaluated.
                ray = inventory.ray(chart, baseline=baseline)
                branch.append('B' if float(ray['profile'].chi) == 0 else 'E' if float(ray['profile'].chi) == 1 else 'M')
                normalized.append(ray)
            for ray in normalized:
                branch.append('C' if ray['b'] < self.raw['robust_parameters'][3] else 'P')
            effective, target, constraints, trace, inventory = self.projected(state)
            for group in self.static:
                if group['block'] != self.d+1:
                    continue
                actual = dict(group, name='actual_lifted_outer_epsilon')
                if not observe(actual, effective):
                    return self.domain_result(False, 'normalization_and_dependent_trace', records, counts, margins, minimum, native)
                at = np.abs(group['points'][:, 2]-self.H) < 1e-11
                if np.any(at):
                    dependent = dict(group, name='projected_dependent_trace', points=group['points'][at],
                        weights=group['weights'][at], side=1, dependent=True)
                    if not observe(dependent, target):
                        return self.domain_result(False, 'normalization_and_dependent_trace', records, counts, margins, minimum, native)
        except Exception as error:
            result = self.domain_result(False, 'moving_profile_or_projected_trace', records, counts, margins, minimum, native)
            result['independent_exception'] = dict(type=type(error).__name__, message=str(error))
            return result
        result = self.domain_result(True, 'complete', records, counts, margins, minimum, native)
        result.update(minimum_b=float(min(ray['b'] for ray in normalized)),
            minimum_lambda=float(min(ray['lam'] for ray in normalized)),
            cutoff_range=[float(min(float(ray['profile'].rc) if ray['profile'].chi > 0 else 0. for ray in normalized)),
                float(max(float(ray['profile'].rc) if ray['profile'].chi > 0 else 0. for ray in normalized))],
            independent_trace=trace, profile_branch_pattern=''.join(branch),
            independent_branch_context=self.branch_context(effective, normalized[:len(self.source_anchors)]),
            effective_state_float64_LE_sha256=vector_hash(effective), trace_target_state_float64_LE_sha256=vector_hash(target))
        if native is not None:
            self.audit.require(bool(native['admissible']), 'Native rejected an independently admissible scalar state')
            for key in ('minimum_b', 'minimum_lambda', 'cutoff_range'):
                self.audit.near(result[key], native[key], 'independent_domain_'+key, 2e-10)
            # Branch ordering can differ only through independent anchor
            # ordering. Compare branch populations; literal native patterns
            # are retained separately for each actual FD endpoint.
            self.audit.require(Counter(result['profile_branch_pattern']) == Counter(native['profile_branch_pattern']),
                'Independent/native source/trace profile branch populations differ')
            if native.get('effective_state'):
                self.audit.near(effective, native['effective_state'], 'independent_external_lift', 2e-10)
                self.audit.near(target, native['trace_target_state'], 'independent_projected_trace_state', 2e-10)
            if native.get('observations'):
                self.observations(native['observations'], effective, target)
                a = sorted(native['observations'], key=self.point_key)
                b = sorted(records, key=self.point_key)
                self.audit.require(len(a) == len(b), 'Native/independent complete common inventory lengths differ')
                self.audit.require([self.point_key(row) for row in a] == [self.point_key(row) for row in b],
                    'Native/independent moving-ray inventory geometry differs')
                self.audit.near([row['coordinate'] for row in a], [row['coordinate'] for row in b],
                    'independent_complete_dynamic_ray_coordinates', 2e-11)
                self.audit.near([row['weight'] for row in a], [row['weight'] for row in b],
                    'independent_complete_dynamic_ray_weights', 2e-11)
        return result

    def domain_result(self, admitted, phase, records, counts, margins, minimum, native):
        result = dict(admissible=bool(admitted), phase=phase, inventory_counts=dict(counts),
            minimum_dimensionless_margin=min(margins) if margins else None, **minimum,
            invalid_points=[row for row in records if not row['admissible']],
            scope='Common actual residual-evaluation inventory; additional diagnostic reconstruction rays are separate')
        if native is not None:
            result['native_phase'] = native['phase']; result['native_failure'] = native['failure']
            result['native_limiting_inventory'] = native.get('limiting_inventory')
            result['native_limiting_point'] = native.get('limiting_point')
            if admitted:
                self.audit.require(counts == Counter(native['inventory_counts']), 'Common scalar inventory counts differ')
                self.audit.near(result['minimum_dimensionless_margin'], native['minimum_dimensionless_margin'],
                    'independent_scalar_domain_minimum_margin', 2e-10)
                for key in minimum:
                    self.audit.near(minimum[key], native[key], 'independent_scalar_domain_'+key, 2e-10)
            elif not native['admissible']:
                result['same_rejection_class'] = True
            else:
                self.audit.require(False, 'Native admitted an independently invalid scalar state')
        if records and admitted:
            result['minimum_points'] = sorted(records, key=lambda row: row['value']/row['scale'])[:4]
        return result

    def branch_context(self, state, source_rays):
        """Independent switch data, not constitutive stability conclusions."""
        patterns = dict(outer_max=[], inner_saturation=[], source_limiter=[])
        distance = {name:math.inf for name in patterns}
        uncertain = 0
        for group in self.static:
            if group['name'] != 'source_volume_outer' or group['block'] != self.d:
                continue
            G, k, epsilon, _ = local_values(self.raw, self.fields, state, group['geometry'], group['points'])
            S, O = G+G.transpose(0, 2, 1), G-G.transpose(0, 2, 1)
            qs = k/epsilon*np.sqrt(np.sum(S*S, axis=(1, 2))/2)
            qo = k/epsilon*np.sqrt(np.sum(O*O, axis=(1, 2))/2)
            gap = qs-qo
            patterns['outer_max'].extend('S' if value >= 0 else 'O' for value in gap)
            distance['outer_max'] = min(distance['outer_max'], float(np.min(np.abs(gap))))
            uncertain += int(np.count_nonzero(np.abs(gap) <= 64*np.finfo(float).eps*np.maximum(1., np.maximum(qs, qo))))
        for ray in source_rays:
            profile = ray['profile']; k = np.asarray(ray['local_k'])
            _, G = _velocity(self.fields, state, self.offsets, self.d, ray['points'])
            S = G+G.transpose(0, 2, 1)
            viscosity = float(profile.a)*np.maximum(ray['r']-float(profile.zv), 0.)
            a = viscosity*np.sqrt(np.sum(S*S, axis=(1, 2)))
            c = 2*k/3
            gap = (a-.8*c)/self.K0
            patterns['inner_saturation'].extend('Z' if kval == 0 else 'A' if value > 0 else 'I'
                for kval, value in zip(k, gap))
            distance['inner_saturation'] = min(distance['inner_saturation'], float(np.min(np.abs(gap))))
            ratio = k/(self.raw['robust_parameters'][2]*float(profile.k))
            patterns['source_limiter'].extend('Z' if value == 0 else 'L' if value < 1 else 'P' for value in ratio)
            distance['source_limiter'] = min(distance['source_limiter'], float(np.min(np.minimum(np.abs(ratio), np.abs(ratio-1)))))
        result = {}
        for name, pattern in patterns.items():
            text = ''.join(pattern)
            result[name] = dict(count=len(pattern), population=dict(Counter(pattern)), lossless_RLE=rle(pattern),
                pattern_sha256=hashlib.sha256(text.encode()).hexdigest(), minimum_switch_distance=distance[name])
        result['outer_max']['binary64_tie_uncertainty_count'] = uncertain
        result['interpretation'] = 'Actual finite evaluation signatures and switch distances; a matching endpoint signature is not a proof that no switch occurs between endpoints'
        return result

    def observations(self, observations, effective, target):
        groups = defaultdict(list)
        for row in observations:
            groups[(int(row['block']), int(row['geometry']), bool(row['dependent_trace']))].append(row)
        for (block, identity, dependent), rows in groups.items():
            values, _ = self.values(target if dependent else effective, block, identity,
                [row['coordinate'] for row in rows])
            self.audit.near(values, [row['value'] for row in rows], 'independent_exported_scalar_observations', 2e-10)
            scales = [self.K0 if block == self.d else self.audit.us**3/(self.audit.kappa*row['coordinate'][2]) for row in rows]
            self.audit.near(scales, [row['scale'] for row in rows], 'independent_exported_scalar_scales', 2e-10)

    def constraints(self, state, entries=None):
        effective, target, constraints, trace, inventory = self.projected(state)
        if entries is not None:
            native = self.audit.constraints(entries)
            self.audit.require(set(native) == set(constraints), 'Current projected constraint-row set differs')
            ids = sorted(constraints)
            self.audit.near([constraints[i] for i in ids], [native[i] for i in ids],
                'independent_current_projected_constraints', 2e-10)
        return constraints

    def time_mass(self, state, history, alpha, constraints=None):
        result = np.zeros(self.size)
        for block in range(self.d+2):
            ids = self.offsets[block]+np.arange(self.fields[block].size)
            result[ids] = alpha*self.audit.mass[block]@(np.asarray(state)[ids]-np.asarray(history)[ids])
        if constraints is not None:
            for row in constraints:
                result[row] = 0.
        return result

    def affine_continuity_matrix(self):
        """Exact polynomial weak-divergence and existing mean/gauge blocks.

        This independent matrix is a comparison, never a replacement for the
        native finite-difference Jacobian. Flat, canonical tensor geometry and
        fixed external velocity values make this block affine in the raw
        coefficients. No pressure sensor, source, roughness force or projected
        epsilon target contributes to these rows.
        """
        if self._continuity is not None:
            return self._continuity, self._continuity_lifted
        pressure_fields = self.fields[self.d+2:self.d+4]
        count = sum(field.size for field in pressure_fields)
        matrix = np.zeros((count+1, self.size), np.longdouble)
        cursor, degrees, native_adequate = 0, [], True
        for block, pressure in enumerate(pressure_fields, self.d+2):
            for element in pressure.elements:
                velocity = self.fields[0].by_id[element['geometry']]
                degree = velocity['degree']+element['degree']
                order = max(8, (degree+2)//2)
                native_adequate = native_adequate and all(2*int(n)-1 >= degree for n in self.raw['quadrature_orders'])
                degrees.append(degree)
                points, weights = self.fields[0].quadrature(velocity, order)
                pi, phi, _ = pressure.rows(element, points)
                vi, _, gradient = self.fields[0].rows(velocity, points)
                for component, axis in enumerate((0, 2)):
                    local = phi.astype(np.longdouble).T@(
                        weights.astype(np.longdouble)[:, None]*gradient[:, axis].astype(np.longdouble))
                    matrix[np.ix_(cursor+pi, self.offsets[component]+vi)] += local
            mean = self.audit.mass[block]@self.constants[block]/self.fields[0].volume
            local = np.arange(pressure.size)
            matrix[cursor+local, -1] = mean
            matrix[-1, self.offsets[block]+local] = mean
            cursor += pressure.size
        matrix = np.asarray(matrix, float)
        lifted = matrix.copy()
        for row in self.external:
            if row < self.offsets[self.d]:
                lifted[:, row] = 0.
        self._continuity, self._continuity_lifted = matrix, lifted
        self._continuity_factors = [np.linalg.cholesky(self.audit.mass[block]) for block in (self.d+2, self.d+3)]
        self.continuity_metadata = dict(shape=list(matrix.shape),
            matrix_float64_LE_sha256=vector_hash(matrix),
            external_lifted_matrix_float64_LE_sha256=vector_hash(lifted),
            polynomial_integrand_degree_maximum=max(degrees), independent_minimum_Gauss_order=8,
            native_declared_rules_polynomially_adequate=bool(native_adequate),
            gauge_means='Physical regional pressure mass times the exact constant, divided by total velocity volume',
            scope='Separately assembled weak continuity/gauge rows; native nonpressure Jacobian columns remain finite differences')
        self.audit.require(native_adequate, 'Native rule cannot integrate the affine weak-divergence polynomials exactly')
        return matrix, lifted

    def continuity_rows(self):
        return np.r_[np.arange(self.offsets[self.d+2], self.offsets[self.d+4]), self.size-1].astype(int)

    def continuity_error_norms(self, error):
        """Pressure mass-dual and gauge norms, with batched column support."""
        self.affine_continuity_matrix()
        error = np.asarray(error, float)
        vector = error.ndim == 1
        if vector:
            error = error[:, None]
        values, cursor = [], 0
        for index, block in enumerate((self.d+2, self.d+3)):
            count = self.fields[block].size
            inverse = np.linalg.solve(self._continuity_factors[index], error[cursor:cursor+count])
            values.append(np.linalg.norm(inverse, axis=0)/((self.audit.us/self.H)*math.sqrt(self.fields[block].volume)))
            cursor += count
        values.append(np.abs(error[-1])/self.audit.us**2)
        result = np.asarray(values)
        return result[:, 0] if vector else result

    def affine_continuity_values(self, state, residual):
        matrix, _ = self.affine_continuity_matrix()
        expected = matrix@self.external_lift(state)
        actual = np.asarray(residual)[self.continuity_rows()]
        self.audit.near(actual, expected, 'independent_exact_affine_weak_continuity_and_gauge_values', 2e-10)
        norms = self.continuity_error_norms(actual-expected)
        self.continuity_gates.append(dict(kind='actual_residual_value', maximum_absolute_physical_error=float(max(norms)),
            gate_met=bool(max(norms) <= ACTION_TOLERANCE)))
        return dict(exact_values=expected.tolist(), absolute_physical_errors=norms.tolist(),
            maximum_absolute_physical_error=float(max(norms)),
            fixed_1e_minus7_gate_met=bool(max(norms) <= ACTION_TOLERANCE),
            mean_gauge_included=True, external_lifting_included=True)

    def affine_continuity_columns(self, matrix, scalar_only):
        _, expected = self.affine_continuity_matrix()
        if scalar_only:
            frozen = np.r_[np.arange(self.offsets[self.d]), np.arange(self.offsets[self.d+2], self.size)]
            actual = np.asarray(matrix[frozen, :])
            identity = np.zeros_like(actual)
            identity[np.arange(len(frozen)), frozen] = 1.
            self.audit.near(actual, identity, 'declared_frozen_velocity_pressure_gauge_identity_rows', 0.)
            return dict(checked_as_coupled=False, declared_frozen_identities_checked=True,
                frozen_rows=len(frozen), maximum_identity_error=float(np.max(np.abs(actual-identity))),
                reason='Frozen-scalar diagnostic identity rows are verified as declared; no coupled continuity derivative claim')
        native = np.asarray(matrix[self.continuity_rows(), :])
        dimensioned = (native-expected)*self.audit.scales()[None, :]
        errors = self.continuity_error_norms(dimensioned)
        maximum = np.max(errors, axis=0)
        self.continuity_gates.append(dict(kind='current_J_columns', maximum_absolute_physical_error=float(max(maximum)),
            gate_met=bool(max(maximum) <= ACTION_TOLERANCE)))
        worst = np.argsort(maximum)[-12:][::-1]
        velocity = np.arange(self.offsets[self.d])
        return dict(checked_as_coupled=True, native_rows_are_finite_differences_for_nonpressure_columns=True,
            dimensionless_column_scale='us, K0, E0, us^2, us/H by variable block',
            maximum_absolute_physical_error=float(np.max(maximum)),
            maximum_velocity_column_error=float(np.max(maximum[velocity])),
            fixed_1e_minus7_gate_met=bool(np.max(maximum) <= ACTION_TOLERANCE),
            all_column_absolute_physical_errors=maximum.tolist(),
            worst_columns=[dict(column=int(index), physical_block_errors=errors[:, index].tolist(),
                input_scale=float(self.audit.scales()[index]),
                externally_fixed_velocity=bool(index in self.external and index < self.offsets[self.d])) for index in worst],
            interpretation='The independent B/mean block is exact. Any native FD continuity-row discrepancy is retained, not replaced or called analytic.')

    def affine_continuity_action(self, direction, action, scalar_only=False):
        _, matrix = self.affine_continuity_matrix()
        expected = np.asarray(direction)[self.continuity_rows()] if scalar_only else matrix@direction
        actual = np.asarray(action)[self.continuity_rows()]
        errors = self.continuity_error_norms(actual-expected)
        self.continuity_gates.append(dict(kind='current_J_direction', maximum_absolute_physical_error=float(max(errors)),
            gate_met=bool(max(errors) <= ACTION_TOLERANCE)))
        return dict(exact_action=expected.tolist(), absolute_physical_errors=errors.tolist(),
            maximum_absolute_physical_error=float(max(errors)),
            fixed_1e_minus7_gate_met=bool(max(errors) <= ACTION_TOLERANCE),
            coupled_continuity_derivative_claim=not bool(scalar_only),
            scope=('Declared frozen pressure/gauge identity rows only' if scalar_only else
                'Only weak-divergence and existing pressure/gauge mean rows; separate from pressure-gradient/roughness and nonlinear sensor/trace actions'))

    def weak(self, state, history, alpha, native_residual, entries=None):
        effective = self.external_lift(state)
        constraints = self.constraints(state, entries)
        loads, details = weak_loads(self.raw, self.fields, effective)
        unmasked = loads+self.time_mass(effective, history, alpha)
        full = unmasked.copy()
        for row, value in constraints.items():
            full[row] = state[row]-value
        self.audit.near(native_residual, full, 'independent_true_unforced_coupled_residual', 2e-10)
        norms = self.audit.norms(full, constraints)
        difference = self.audit.norm_record(full-np.asarray(native_residual), constraints)
        return dict(physical_norms=norms.tolist(), maximum_physical_norm=float(max(norms)),
            maximum_component_error=float(np.max(np.abs(full-np.asarray(native_residual)))),
            independent_difference_norm=difference, quadrature=details,
            raw_unmasked_load_float64_LE_sha256=vector_hash(unmasked),
            exact_affine_weak_continuity=self.affine_continuity_values(state, native_residual),
            masked_true_residual_float64_LE_sha256=vector_hash(full)), unmasked, constraints

    def continuity(self, state):
        result = np.zeros(5)
        field = self.fields[0]
        nodes, weights = np.polynomial.legendre.leggauss(8)
        for element in field.elements:
            points, measure = field.quadrature(element, max(8, element['degree']+3))
            _, G = _velocity(self.fields, state, self.offsets, self.d, points, element['geometry'])
            divergence = np.trace(G, axis1=1, axis2=2)
            region = 0 if element['upper'][2] <= self.H else 1
            result[region] += measure@divergence**2; result[region+2] += measure@divergence
            for axis in (0, 2):
                other = 2 if axis == 0 else 0
                for side in (0, 1):
                    endpoint = element['upper' if side else 'lower'][axis]
                    if abs(endpoint-(field.upper if side else field.lower)[axis]) > 1e-11:
                        continue
                    points = np.zeros((8, 3)); points[:, axis] = endpoint
                    width = element['upper'][other]-element['lower'][other]
                    points[:, other] = element['lower'][other]+(nodes+1)*width/2
                    U, _ = _velocity(self.fields, state, self.offsets, self.d, points, element['geometry'])
                    result[4] += (1 if side else -1)*float(weights@U[:, axis])*width/2
        result[:2] = np.sqrt(result[:2])
        return result

    def restricted_ranges(self, state, native_domain, native_audit=None):
        """Independent physical qualification on the actual common inventory.

        The dependent target and retained actual epsilon-H state both enter
        this inventory, with their own native call identities. The broader E
        reconstruction audit is reported separately and never substituted.
        """
        observations = native_domain.get('observations', [])
        self.audit.require(bool(observations), 'Restricted physical audit lacks its common-inventory coordinates')
        effective, target, _, _, _ = self.projected(state)
        groups = defaultdict(list)
        for row in observations:
            if int(row['block']) == self.d or row['dependent_trace']:
                groups[(int(row['geometry']), int(row['side']), bool(row['dependent_trace']))].append(row)
        outer = np.array([0., math.inf, -math.inf, math.inf, -math.inf, math.inf, -math.inf])
        covariance = [math.inf, math.inf]
        worst = [None, None]
        profiles = {}
        flow, normal = np.asarray(self.raw['robust_reference_flow']), np.asarray(self.raw['robust_reference_normal'])
        for (identity, side, dependent), rows in groups.items():
            points = np.asarray([row['coordinate'] for row in rows])
            _, G = _velocity(self.fields, effective, self.offsets, self.d, points, identity)
            k, _ = self.values(effective, self.d, identity, points)
            for index, point in enumerate(points):
                inner = bool(point[2] < self.H-1e-11 or abs(point[2]-self.H) < 1e-11 and side < 0)
                if inner:
                    key = (float(point[0]), float(point[1]))
                    if key not in profiles:
                        baseline = dict(self.raw, upstream_transition=False, chi=0.)
                        profiles[key] = _profile(baseline, self.fields, effective, point)[0]
                    profile = profiles[key]
                    viscosity = float(profile.a)*max(float(point[2])-float(profile.zv), 0.)
                    law = selected_inner(self.raw, G[index], float(k[index]), viscosity)
                    stress = np.asarray(law['stress']).real
                    value = float(np.linalg.eigvalsh(stress)[0])
                    if value < covariance[0]:
                        covariance[0] = value
                        worst[0] = dict(chart=point.tolist(), geometry=identity, side=side,
                            k=float(k[index]), sampled_K=float(profile.k), nu_A=viscosity,
                            a=float(np.asarray(law['a']).real), c=float(np.asarray(law['c']).real),
                            G=G[index].tolist(), S=(G[index]+G[index].T).tolist(),
                            covariance_minimum=value, inventory=rows[index]['inventory'])
                    continue
                epsilon_state = target if dependent else effective
                epsilon = self.fields[self.d+1].evaluate(epsilon_state, self.offsets[self.d+1], point,
                    1 if abs(point[2]-self.H) < 1e-11 else 0)[0]
                self.audit.require(k[index] > 0 and epsilon > 0, 'Common physical audit exited the scalar domain')
                stress, _ = selected_stress(self.raw, G[index], float(k[index]), epsilon, point)
                value = float(np.linalg.eigvalsh(stress)[0])
                if value < covariance[1]:
                    covariance[1] = value
                    worst[1] = dict(chart=point.tolist(), geometry=identity, side=side, dependent_trace=dependent,
                        k=float(k[index]), epsilon=float(epsilon), G=G[index].tolist(),
                        covariance_minimum=value, inventory=rows[index]['inventory'])
                gradient_reference = self.audit.us/(self.audit.kappa*point[2])
                epsilon_reference = self.audit.us**3/(self.audit.kappa*point[2])
                outer[0] = max(outer[0], float(np.max(np.abs(G[index]/gradient_reference-np.outer(flow, normal)))))
                kr, er = float(k[index]/self.K0), float(epsilon/epsilon_reference)
                outer[1:5] = [min(outer[1], kr), max(outer[2], kr), min(outer[3], er), max(outer[4], er)]
        sampled = [row for row in observations if row['inventory'] == 'physical_sampled_K']
        samples = []
        for row in sampled:
            value, _ = self.values(effective, self.d, int(row['geometry']), [row['coordinate']])
            samples.append(float(value[0]/self.K0))
        outer[5], outer[6] = min(samples), max(samples)
        box = bool(outer[0] <= .001 and outer[1] >= .99 and outer[2] <= 1.01 and
            outer[3] >= .98 and outer[4] <= 1.02 and outer[5] >= .99 and outer[6] <= 1.01)
        covariance_guard = bool(min(covariance) >= -1e-10*self.K0)
        independent_domain = self.domain(state, native_domain)
        qualified = bool(box and covariance_guard and independent_domain['minimum_lambda'] >= .5)
        if native_audit:
            self.audit.near(native_audit['r2_outer_range'], outer, 'independent_common_inventory_fixed_outer_range', 3e-10)
            self.audit.near(native_audit['r2_covariance'], covariance, 'independent_common_inventory_selected_covariance', 3e-10)
            self.audit.require(bool(native_audit['r2_qualified'][0]) == qualified,
                'Native restricted qualification disagrees with independent common-inventory gates')
        return dict(outer_range=outer.tolist(), covariance_minima=covariance, worst_points=worst,
            fixed_outer_box_covered=box, strict_covariance_nonnegative=bool(min(covariance) >= 0),
            declared_native_covariance_roundoff_guard=covariance_guard,
            covariance_guard_absolute_tolerance=1e-10*self.K0,
            minimum_lambda=independent_domain['minimum_lambda'], qualified_on_common_inventory=qualified,
            rigorous_or_polynomial_everywhere_qualification=False,
            scope='Actual common scalar inventory including both retained epsilon-H and independently projected target; one-sided geometry gradients preserved')

    def physical(self, state, history, alpha, native=None, unmasked=None, constraints=None, native_domain=None):
        effective = self.external_lift(state)
        inventory = Inventory(self.raw, self.fields, effective, self.session)
        direct = selected_direct(self.raw, self.fields, effective, self.audit)
        interface = selected_interface(self.raw, self.fields, effective)
        continuity = self.continuity(effective)
        storage, stored = _mass(self.raw, self.fields, effective, history, alpha)
        subtotal = storage+direct['stress_work_minus_k_source']+direct['resolved_advection']+direct['roughness_force_work']
        energy = dict(storage=float(storage), stored_energy=float(stored),
            stress_work_minus_k_source=direct['stress_work_minus_k_source'],
            resolved_advection=direct['resolved_advection'], roughness_force_work=direct['roughness_force_work'],
            prescribed_body_work=0., direct_sum=float(subtotal), anchor_work=0.)
        if unmasked is not None:
            work = _work(effective, unmasked, self.fields, self.offsets, self.d)
            reaction = _work(effective, unmasked, self.fields, self.offsets, self.d, set(constraints))
            energy.update(actual_unmasked_weak_work=float(work), essential_reaction=float(reaction),
                free_weak_work=float(work-reaction), direct_remainder=float(work-subtotal))
        if native:
            for key in ('robust_source', 'wall_sources', 'constitutive_modification'):
                if key in native:
                    self.audit.near(native[key], direct[key], 'independent_stage_'+key, 3e-10)
            if 'inner_stress' in native:
                self.audit.near(native['inner_stress'], direct['inner_stress'], 'independent_stage_selected_inner_stress', 3e-10)
            if 'inner_constitutive_modification' in native:
                self.audit.near(native['inner_constitutive_modification'], direct['inner_constitutive_modification'],
                    'independent_stage_inner_production_transfer', 3e-10)
            if 'interface' in native:
                self.audit.near(native['interface'], interface, 'independent_stage_full_traction_and_k_flux', 3e-10)
            if 'strong_continuity' in native:
                self.audit.near(native['strong_continuity'], continuity, 'independent_stage_strong_continuity', 3e-10)
            if 'energy' in native and unmasked is not None:
                expected = [energy[key] for key in ('actual_unmasked_weak_work', 'storage', 'stress_work_minus_k_source',
                    'resolved_advection', 'roughness_force_work', 'prescribed_body_work', 'direct_sum',
                    'direct_remainder', 'stored_energy', 'essential_reaction', 'anchor_work', 'free_weak_work')]
                self.audit.near(native['energy'], expected, 'independent_stage_complete_energy_and_reaction', 3e-10)
            for key, independent in (('source_volume_inner', 'inner_covariance_minimum'),
                    ('source_volume_outer', 'outer_covariance_minimum')):
                if key in native:
                    self.audit.near(native[key][3], direct[independent], 'independent_stage_'+key+'_covariance', 3e-10)
        try:
            ranges = finite_ranges(self.raw, self.fields, effective, self.session, inventory)
        except Exception as error:
            ranges = dict(complete=False, type=type(error).__name__, message=str(error),
                finite_joint_preflight_passed=False,
                scope='Independent additional physical reconstruction inventory; failure is not relabelled as a common scalar-domain rejection')
        return dict(source_energy=direct, energy=energy, full_total_traction_and_k_flux=interface.tolist(),
            tangential_traction_RMS=float(math.sqrt(max(0., interface[3]**2-interface[2]**2))),
            strong_continuity=continuity.tolist(), physical_ranges=ranges,
            pressure_jump_is_not_a_pressure_continuity_condition=True,
            constraints_and_projected_trace_are_separate_from_physical_covariance=True,
            common_inventory_physical_qualification=self.restricted_ranges(state, native_domain, native)
                if native_domain and native_domain.get('observations') else None)

    def state(self, state, history, alpha, native_residual, entries, domain, full=False, native_audit=None):
        key = (vector_hash(state), vector_hash(history), float(alpha))
        if key in self.state_cache and (not full or 'physical' in self.state_cache[key]):
            result = self.state_cache[key]
            if domain is not None:
                self.domain(state, domain, retain=full)
            if entries is not None:
                self.constraints(state, entries)
            if native_residual is not None and len(native_residual):
                self.audit.near(native_residual, result['residual'], 'same_state_same_alpha_true_residual')
            return result
        result = dict(state_float64_LE_sha256=key[0], history_float64_LE_sha256=key[1], alpha=float(alpha),
            state=np.asarray(state).tolist(), scalar_domain=self.domain(state, domain, retain=full))
        if result['scalar_domain']['admissible'] and native_residual is not None and len(native_residual):
            residual, unmasked, constraints = self.weak(state, history, alpha, native_residual, entries)
            result.update(weak_residual=residual, residual=np.asarray(native_residual).tolist(),
                constraints=[[row, value] for row, value in sorted(constraints.items())])
            if full:
                result['physical'] = self.physical(state, history, alpha, native_audit, unmasked, constraints)
        self.state_cache[key] = result
        return result

    def matrix(self, record):
        path = Path(record['jacobian_binary_file'])
        if self.matrix_root:
            path = self.matrix_root/path.name
        path = path.resolve()
        self.audit.require(path.is_file(), 'Missing retained current Jacobian: '+str(path))
        self.audit.require(record['jacobian_binary_shape'] == [self.size, self.size], 'Unexpected current Jacobian shape')
        self.audit.require(path.stat().st_size == 8*self.size*self.size, 'Truncated/different-shape current Jacobian')
        self.matrix_files[str(path)] = dict(bytes=path.stat().st_size, sha256=sha(path),
            shape=[self.size, self.size], dtype='float64', byteorder=sys.byteorder,
            retention='Dense matrix stays in the authoritative Linux run directory')
        matrix = np.memmap(path, dtype=np.float64, mode='r', shape=(self.size, self.size), order='C')
        self.audit.require(bool(np.all(np.isfinite(matrix))), 'Current Jacobian contains nonfinite entries')
        return matrix

    def action(self, item, matrix, state, history, alpha, constraints, scalar_only):
        direction = np.asarray(item['direction'], float)
        action = matrix@direction
        self.audit.near(action, item['jacobian_action'], 'independent_dense_current_J_direction')
        exact_mass = self.time_mass(direction, np.zeros(self.size), alpha, constraints)
        self.audit.near(exact_mass, item['exact_mass_action'], 'independent_exact_alpha_mass_action')
        mask = np.ones(self.size, bool)
        if scalar_only:
            mask[:self.offsets[self.d]] = False; mask[self.offsets[self.d+2]:] = False
        base = np.asarray(item['spatial_base'])
        spatial_expected = np.asarray(item['_full_residual'])-self.time_mass(
            self.external_lift(state), history, alpha, constraints)
        self.audit.near(base, spatial_expected, 'independent_spatial_plus_exact_mass_residual_partition', 2e-10)
        action_norm = self.audit.norm_record(action, constraints)
        scale = max(1., float(action_norm['maximum']))
        levels, strict_errors = [], []
        self.audit.near([level['step'] for level in item['levels']], LADDER, 'fixed_dimensionless_current_action_ladder', 0.)
        for level in item['levels']:
            h, sides = float(level['step']), {}
            for name, sign in (('plus', 1), ('minus', -1)):
                side, endpoint = level[name], np.asarray(state)+sign*h*direction
                scalar = self.domain(endpoint, side['domain'])
                result = dict(evaluated=bool(side['evaluated']), scalar_domain=scalar,
                    native_profile_branch_pattern=side['domain'].get('profile_branch_pattern', ''),
                    base_to_endpoint_branch_pattern_changed=bool(side['domain'].get('profile_branch_pattern', '') !=
                        item.get('_base_branches', '')),
                    exception_type=side.get('exception_type'), exception_message=side.get('exception_message'))
                result['independent_branch_comparison'] = branch_difference(item.get('_base_independent_branches', {}),
                    scalar.get('independent_branch_context', {}))
                if side['evaluated']:
                    self.audit.require(scalar['admissible'], 'Directional residual evaluated outside independent scalar domain')
                    c = self.constraints(endpoint, side.get('constraints'))
                    residual = np.asarray(side['spatial_residual'])
                    for row, value in c.items():
                        self.audit.near(residual[row], endpoint[row]-value, 'independent_directional_endpoint_constraint_rows', 2e-10)
                    one_sided = sign*(residual-base)/h+exact_mass
                    difference = one_sided-action; difference[~mask] = 0.
                    result['one_sided_J_error'] = self.audit.norm_record(difference, constraints)
                    result['relative_one_sided_J_error'] = float(result['one_sided_J_error']['maximum']/scale)
                sides[name] = result
            receipt = dict(step=h, **sides)
            if all(level[name]['evaluated'] for name in ('plus', 'minus')):
                central = (np.asarray(level['plus']['spatial_residual'])-
                    np.asarray(level['minus']['spatial_residual']))/(2*h)+exact_mass
                difference = central-action; difference[~mask] = 0.
                norms = self.audit.norm_record(difference, constraints)
                relative = float(norms['maximum']/scale)
                met = bool(norms['maximum'] <= ACTION_TOLERANCE)
                receipt.update(central_J_error=norms, relative_central_J_error=relative,
                    fixed_1e_minus7_gate_met=met,
                    relative_1e_minus7_gate_met=bool(relative <= ACTION_TOLERANCE),
                    complete_constraint_row_error=float(norms['maximum_constraint']))
                strict_errors.append((h, float(norms['maximum']), relative, met))
            levels.append(receipt)
        finest = min(strict_errors, default=None)
        gate = dict(kind=item['kind'], available=finest is not None, threshold=ACTION_TOLERANCE,
            finest_feasible_central_step=finest[0] if finest else None,
            finest_absolute_physical_error=finest[1] if finest else None,
            finest_relative_error=finest[2] if finest else None,
            gate_met=bool(finest[3]) if finest else False,
            threshold_definition='Absolute dimensionless physical mass-dual/block/constraint norm; relative errors are reported separately',
            qualification='Current represented action; branch/admissibility evidence retained. No smoothness claim across a switch.')
        self.gates.append(gate)
        analytic_pressure = None
        if item['kind'] == 'pressure_and_gauge' and not scalar_only:
            view = dict(self.raw, state=np.asarray(state).tolist(), constraints=[[i, v] for i, v in sorted(constraints.items())])
            oracle = PressureOracle(view, self.external_lift(state))
            exact_affine = oracle.affine(direction)
            nonlinear, records = oracle.source_and_trace(direction)
            analytic = exact_affine+nonlinear
            norm = self.audit.norm_record(action-analytic, constraints)
            analytic_pressure = dict(independent_exact_complete_action_error=norm,
                fixed_1e_minus7_gate_met=bool(norm['maximum'] <= ACTION_TOLERANCE),
                nonlinear_source_and_projected_trace_norm=self.audit.norm_record(nonlinear, constraints),
                ray_records=records, interpretation='Independent exact affine plus analytic sensor/chi/profile/cap/trace chain; not a full-residual finite difference')
        return dict(kind=item['kind'], dimensionless_direction_maximum=float(np.max(np.abs(direction/self.audit.scales()))),
            native_dense_product_checked=True, exact_physical_mass_checked=True,
            action_norm=action_norm, constraint_row_actions=action[sorted(constraints)].tolist(),
            matching_epsilon_row_actions=action[sorted(self.trace_rows)].tolist(), levels=levels, gate=gate,
            independent_affine_weak_continuity=self.affine_continuity_action(direction, action, scalar_only),
            frozen_scalar_masked_rows_not_claimed_as_coupled=bool(scalar_only), independent_analytic_pressure=analytic_pressure)

    def independent_outer_k_action(self, matrix, state, history, alpha, constraints, scalar_only, base_domain):
        """Additional outer-k evidence absent from the native inner bubble.

        Free matching-height coefficients belong to both native regions; their
        inner CG support is retained rather than forcibly discontinuous. Every
        endpoint is independently screened and assembled, with no native call.
        """
        field, offset = self.fields[self.d], self.offsets[self.d]
        outer_ids = {int(index) for element in field.elements if element['lower'][2] >= self.H-1e-11
            for index in element['ids']}
        selected = sorted(index for index in outer_ids if offset+index not in self.external)
        self.audit.require(bool(selected), 'No free represented outer-k coefficient exists')
        direction = np.zeros(self.size)
        for index in selected:
            direction[offset+index] = self.K0*math.sin(.917*(offset+index+1)+.271)
        direction *= self.K0/float(np.max(np.abs(direction)))
        matching = field.boundary(lambda axis, side, element:
            axis == 2 and side == 0 and abs(element['lower'][2]-self.H) < 1e-11)
        shared = sorted(set(selected).intersection(matching))
        self.audit.require(bool(shared), 'Additional outer-k direction omitted the matching trace')
        action = matrix@direction
        mass = self.time_mass(direction, np.zeros(self.size), alpha, constraints)
        mask = np.ones(self.size, bool)
        if scalar_only:
            mask[:self.offsets[self.d]] = False; mask[self.offsets[self.d+2]:] = False
        def endpoint(point):
            effective, _, current, trace, _ = self.projected(point)
            loads, details = weak_loads(self.raw, self.fields, effective)
            result = loads.copy()
            for row, value in current.items():
                result[row] = point[row]-value
            self.audit.require(bool(np.all(np.isfinite(result))), 'Nonfinite independent outer-k endpoint residual')
            return result, current, trace, details
        base, _, _, _ = endpoint(np.asarray(state))
        action_norm = self.audit.norm_record(np.where(mask, action, 0.), constraints)
        scale = max(1., float(action_norm['maximum']))
        levels, errors = [], []
        for h in LADDER:
            sides, residuals = {}, {}
            for name, sign in (('plus', 1), ('minus', -1)):
                point = np.asarray(state)+sign*h*direction
                domain = self.domain(point)
                receipt = dict(evaluated=False, endpoint_source='independent weak oracle; not a native recorded endpoint',
                    state_float64_LE_sha256=vector_hash(point), scalar_domain=domain,
                    branches=branch_difference(base_domain.get('independent_branch_context', {}),
                        domain.get('independent_branch_context', {})))
                if domain['admissible']:
                    try:
                        residual, current, trace, details = endpoint(point)
                        one_sided = sign*(residual-base)/h+mass-action
                        one_sided[~mask] = 0.
                        receipt.update(evaluated=True, spatial_residual=residual.tolist(),
                            spatial_residual_float64_LE_sha256=vector_hash(residual),
                            constraints=[[int(row), float(value)] for row, value in sorted(current.items())],
                            independent_trace=trace, quadrature=details,
                            one_sided_J_error=self.audit.norm_record(one_sided, constraints))
                        residuals[name] = residual
                    except Exception as error:
                        receipt['exception'] = dict(type=type(error).__name__, message=str(error))
                sides[name] = receipt
            row = dict(step=h, **sides)
            if all(sides[name]['evaluated'] for name in ('plus', 'minus')):
                central = (residuals['plus']-residuals['minus'])/(2*h)+mass
                difference = central-action; difference[~mask] = 0.
                norm = self.audit.norm_record(difference, constraints)
                absolute = float(norm['maximum'])
                amplitude = max(float(self.audit.norm_record(residuals[name], constraints)['maximum']) for name in residuals)
                row.update(central_J_error=norm, relative_central_J_error=absolute/scale,
                    fixed_1e_minus7_gate_met=bool(absolute <= ACTION_TOLERANCE),
                    relative_1e_minus7_gate_met=bool(absolute/scale <= ACTION_TOLERANCE),
                    arithmetic_subtraction_scale_proxy=float(64*np.finfo(float).eps*max(1., amplitude)/h),
                    uncertainty='The subtraction scale is a proxy, not an error bound; branch changes and all absolute failures remain recorded')
                errors.append((h, absolute, bool(absolute <= ACTION_TOLERANCE)))
            levels.append(row)
        transfers = dict(matching_H=[], physical_sensor_A=[], physical_sensor_B=[])
        maximum_tangential_gradient = 0.
        for chart, _ in self.source_anchors:
            for name, radial in (('matching_H', self.H), ('physical_sensor_A', self.raw['transition_controls'][2]),
                    ('physical_sensor_B', self.raw['transition_controls'][3])):
                point = np.array([chart[0], chart[1], radial])
                value, gradient = field.evaluate(direction, offset, point, 1 if radial == self.H else 0)
                transfers[name].append(float(value/self.K0))
                if name == 'matching_H':
                    maximum_tangential_gradient = max(maximum_tangential_gradient, float(abs(gradient[0])))
        finest = min(errors, default=None)
        gate = dict(kind='outer_k_independent_oracle', available=finest is not None, threshold=ACTION_TOLERANCE,
            finest_feasible_central_step=finest[0] if finest else None,
            finest_absolute_physical_error=finest[1] if finest else None,
            gate_met=bool(finest[2]) if finest else False,
            scope='Independent weak-oracle differences of a fixed outer-support direction; no native probe or numerical search')
        self.gates.append(gate)
        return dict(kind=gate['kind'], direction=direction.tolist(), selected_outer_support_DOFs=selected,
            selected_shared_matching_DOFs=shared,
            matching_trace_has_unavoidable_inner_CG_support=True,
            independent_endpoint_source=True, native_endpoint_calls=0,
            exact_mass_action=mass.tolist(), action_norm=action_norm,
            actual_current_J_action=action.tolist(), levels=levels, gate=gate,
            outer_k_transfer_dimensionless_values=transfers,
            maximum_matching_tangential_k_gradient=maximum_tangential_gradient,
            maximum_mechanical_forcing_minus_two_thirds_grad_k_direction=2*maximum_tangential_gradient/3,
            scalar_block_only_comparison=bool(scalar_only),
            native_five_families_not_claimed_to_cover_this_direction=True,
            independent_affine_weak_continuity=self.affine_continuity_action(direction, action, scalar_only))

    def fraction_bound(self, state, step, matrix, domain):
        actual = self.external_lift(np.asarray(state)+step)-self.external_lift(state)
        predicted = actual.copy(); jacstep = matrix@step
        for row in self.trace_rows:
            predicted[row] = step[row]-jacstep[row]
        maximum, limiting = 1., None
        self.audit.require(bool(domain.get('observations')), 'Missing actual current scalar observations for fraction bound')
        for point in domain['observations']:
            values, _ = self.values(predicted if point['dependent_trace'] else actual,
                int(point['block']), int(point['geometry']), [point['coordinate']])
            delta = float(values[0])
            if delta < 0 and point['value']/(-delta) < maximum:
                maximum = max(0., float(point['value']/(-delta)))
                limiting = dict(point, independent_directional_change=delta)
        return maximum, limiting

    def column_controls(self, columns, state, constraints, scalar_only):
        self.audit.require([column['column'] for column in columns] == list(range(self.size)),
            'Current Jacobian column provenance is incomplete')
        scales = self.audit.scales()
        stencils, rejected, one_sided = Counter(), Counter(), []
        for column in columns:
            index, stencil = int(column['column']), column['stencil']
            block = next((b for b, field in enumerate(self.fields)
                if self.offsets[b] <= index < self.offsets[b]+field.size), len(self.fields))
            stencils[stencil] += 1
            if scalar_only and block not in (self.d, self.d+1):
                self.audit.require(stencil == 'frozen_scalar_diagnostic_not_a_coupled_column',
                    'Frozen diagnostic column was mislabelled as a coupled derivative')
                continue
            if block >= self.d+2:
                self.audit.require(stencil == 'exact_affine_and_full_nonlinear_pressure_action',
                    'A pressure/gauge column replaced or omitted its full analytic action')
                continue
            if index in constraints and index not in self.trace_rows:
                self.audit.require(stencil == 'lifted_constraint_identity', 'External fixed column is not the lift identity')
                continue
            self.audit.require(stencil in ('scaled_central', 'admissible_one_sided_first_order',
                'admissible_one_sided_second_order'), 'Unknown scalar-feasible finite-difference stencil')
            initial = 2e-6*max(float(scales[index]), abs(float(state[index])))
            self.audit.require(any(abs(column['step']-initial*shrink) <= 2e-13*initial
                for shrink in (1., .25, .0625, .015625, .00390625)),
                'Column derivative step differs from its fixed dimensionless shrink ladder')
            probes = column['probes']
            self.audit.require(1 <= len(probes) <= 5, 'Missing or excessive column feasibility probes')
            self.audit.near([probe['step']/initial for probe in probes],
                [1., .25, .0625, .015625, .00390625][:len(probes)], 'fixed_column_shrink_ladder', 2e-13)
            for probe in probes:
                for name in ('plus', 'minus'):
                    domain = probe[name+'_domain']
                    self.audit.require(not probe[name+'_evaluated'] or domain['admissible'],
                        'Column residual evaluated after its scalar-domain rejection')
                    if not probe[name+'_evaluated']:
                        rejected[domain['phase']+': '+domain['failure']] += 1
            if stencil == 'scaled_central':
                self.audit.require(probes[-1]['plus_evaluated'] and probes[-1]['minus_evaluated'],
                    'Central column has no two feasible residual evaluations')
                self.audit.near(column['step'], probes[-1]['step'], 'selected_central_column_step', 2e-13*initial)
            else:
                self.audit.require(column['side'] in (-1, 1), 'One-sided column has no declared admissible side')
                one_sided.append(dict(column=index, stencil=stencil, side=column['side'], step=column['step'],
                    probes=probes))
        return dict(stencils=dict(stencils), rejected_probe_classifications=dict(rejected),
            one_sided_columns=one_sided,
            scope='All column topology, dimensions, steps and declared pre-residual checks inspected; current five-family matrix actions are independently checked separately')

    def iteration(self, record, history, alpha, scalar_only):
        state, residual, step = (np.asarray(record[key], float) for key in ('state', 'residual', 'step'))
        current = self.state(state, history, alpha, residual, record.get('constraints'), record['domain'], full=True)
        constraints = self.constraints(state, record.get('constraints'))
        norms = self.audit.norms(residual, constraints)
        self.audit.near(norms, record['scaled_block_norms'], 'independent_iteration_physical_mass_norms')
        matrix = self.matrix(record); linear = matrix@step+residual
        if scalar_only:
            linear[:self.offsets[self.d]] -= residual[:self.offsets[self.d]]
            linear[self.offsets[self.d+2]:] -= residual[self.offsets[self.d+2]:]
        linear_norms = self.audit.norms(linear, constraints)
        self.audit.near(linear_norms, record['linear_scaled_block_norms'], 'independent_current_linear_residual_norms', 2e-10)
        dimensionless = float(np.max(np.abs(step/self.audit.scales())))
        self.audit.near(dimensionless, record['dimensionless_step'], 'independent_current_dimensionless_Newton_step')
        error = float(max(norms[self.d], norms[self.d+1], norms[-1]) if scalar_only else max(norms))
        self.audit.require(max(linear_norms) < 1e-8*max(1., error), 'Current dense linear solve exceeds its fixed residual guard')
        maximum, limiting = self.fraction_bound(state, step, matrix, record['domain'])
        self.audit.near(maximum, record['fixed_location_maximum_fraction'], 'independent_fraction_to_boundary_prediction', 2e-10)
        self.audit.require(bool(record['boundary_limited']) == bool(maximum < 1), 'Fraction-limit classification disagrees with actual bounds')
        if maximum < 1:
            possible = [point for point in record['domain']['observations'] if point['inventory'] == record['limiting_inventory']
                and np.max(np.abs(np.asarray(point['coordinate'])-record['limiting_point'])) < 2e-11]
            self.audit.require(bool(possible), 'Native fraction-limiting point is absent from the actual scalar inventory')
            self.audit.require(any(abs(point['value']-record['limiting_value']) < 2e-10 for point in possible),
                'Native fraction-limiting scalar value differs from its actual observation')
        fraction = .995*maximum if maximum < 1 else 1.
        trials, accepted = [], None
        for index, trial in enumerate(record['trials']):
            h = float(trial['fraction'])
            self.audit.near(h, fraction*2.**(-index), 'actual_fraction_to_boundary_backtracking_ladder', 2e-10)
            endpoint = state+h*step
            if trial.get('state'):
                self.audit.near(endpoint, trial['state'], 'actual_current_line_search_trial_state', 2e-10)
            scalar = self.domain(endpoint, trial.get('domain'))
            classification = trial['classification']
            entry = dict(fraction=h, classification=classification, scalar_domain=scalar,
                native_profile_branches_changed=trial.get('profile_branches_changed'),
                exception_type=trial.get('type'), exception_message=trial.get('message'))
            if 'residual' in trial:
                c = self.constraints(endpoint, trial.get('constraints'))
                n = self.audit.norms(trial['residual'], c)
                self.audit.near(n, trial['scaled_block_norms'], 'independent_feasible_trial_physical_norms', 2e-10)
                next_error = float(max(n[self.d], n[self.d+1], n[-1]) if scalar_only else max(n))
                descent = bool(next_error < (1.-1e-4*h)*error)
                admissible_descent = descent and scalar['admissible']
                self.audit.require((classification == 'feasible_accepted_descent') == admissible_descent,
                    'Actual true-norm Armijo/domain decision disagrees with line-search classification')
                for row, value in c.items():
                    self.audit.near(trial['residual'][row], endpoint[row]-value, 'actual_trial_projected_constraint_rows', 2e-10)
                entry.update(physical_norms=n.tolist(), tested_residual_maximum=next_error,
                    Armijo_target=(1.-1e-4*h)*error, strict_descent=descent)
                if admissible_descent:
                    self.audit.require(accepted is None, 'More than one accepted trial in a Newton iteration')
                    accepted = endpoint
                    entry['independent_accepted_state'] = self.state(endpoint, history, alpha,
                        trial['residual'], trial.get('constraints'), trial.get('domain'), full=True)
            elif classification == 'scalar_domain_rejection_before_residual':
                self.audit.require(not scalar['admissible'], 'Native scalar rejection has no independent invalid domain')
            elif classification == 'evaluable_scalar_state_residual_exception':
                self.audit.require(scalar['admissible'], 'Residual exception was mislabelled as a scalar-valid candidate')
            else:
                self.audit.require(False, 'Unknown/missing line-search trial evidence')
            trials.append(entry)
        self.audit.require(bool(accepted is not None) == bool(record['accepted']), 'Accepted trial/iteration flags disagree')
        self.audit.near(record['accepted_fraction'], trials[-1]['fraction'] if accepted is not None else 0.,
            'actual_accepted_fraction', 2e-10)
        columns = self.column_controls(record['columns'], state, constraints, scalar_only)
        actions = []
        for item in record['current_directional_checks']:
            item = dict(item, _base_branches=record['domain'].get('profile_branch_pattern', ''), _full_residual=residual,
                _base_independent_branches=current['scalar_domain'].get('independent_branch_context', {}))
            actions.append(self.attempt('current action '+item['kind'], lambda item=item:
                self.action(item, matrix, state, history, alpha, constraints, scalar_only)))
        actions.append(self.attempt('independent outer-k current action', lambda:
            self.independent_outer_k_action(matrix, state, history, alpha, constraints, scalar_only,
                current['scalar_domain'])))
        return dict(iteration=record['iteration'], state=current, accepted=bool(record['accepted']),
            dimensionless_step=dimensionless, linear_physical_norms=linear_norms.tolist(),
            linear_residual_float64_LE_sha256=vector_hash(linear),
            fraction_to_boundary=dict(maximum=maximum, limiting_point=limiting,
                fixed_location_prediction_is_not_an_admissibility_proof=True),
            current_column_controls=columns, column_metadata_sha256=canonical(record['columns']),
            independent_affine_weak_continuity_J_columns=self.affine_continuity_columns(matrix, scalar_only),
            actions=actions, trials=trials,
            next_accepted_state_float64_LE_sha256=vector_hash(accepted) if accepted is not None else None)

    def solve(self, record):
        history, state = np.asarray(record['history'], float), np.asarray(record['state'], float)
        alpha = float(record['alpha']); scalar_only = bool(record.get('frozen_velocity_pressure_diagnostic'))
        self.audit.near(history, self.raw['initial'], 'fixed_actual_initial_stage_history', 0.)
        self.audit.near(record['time'], self.raw['stage_test_target_time'], 'fixed_physical_stage_time', 0.)
        self.audit.require(record['history_advanced'] is False and record['residual_unchanged'] is True,
            'Auxiliary solve advanced history or changed the residual')
        iterations = []
        for index, iteration in enumerate(record['iterations']):
            self.audit.require(iteration['iteration'] == index, 'Newton iteration numbering is incomplete')
            if index and iterations[-1].get('next_accepted_state_float64_LE_sha256'):
                self.audit.require(iterations[-1]['next_accepted_state_float64_LE_sha256'] == vector_hash(iteration['state']),
                    'Current Newton base is not the preceding accepted internal state')
            iterations.append(self.attempt('alpha '+str(alpha)+' iteration '+str(index), lambda iteration=iteration:
                self.iteration(iteration, history, alpha, scalar_only)))
        final = self.attempt('alpha '+str(alpha)+' final state', lambda: self.state(state, history, alpha,
            record.get('residual'), record.get('constraints'), record['domain'], full=True))
        if record.get('residual'):
            constraints = self.constraints(state, record.get('constraints'))
            norms = self.audit.norms(record['residual'], constraints)
            self.audit.near(norms, record['scaled_block_norms'], 'independent_final_true_mass_dual_norms', 2e-10)
            threshold = float(record['tolerance'])
            solved_error = float(max(norms[self.d], norms[self.d+1], norms[-1]) if scalar_only else max(norms))
            self.audit.require(not record['converged'] or solved_error < threshold,
                'Native converged flag does not meet the actual solved residual tolerance')
        counts = Counter(trial['classification'] for iteration in record['iterations'] for trial in iteration['trials'])
        return dict(alpha=alpha, time=record['time'], terminal=record['terminal'], converged=bool(record['converged']),
            frozen_velocity_pressure_diagnostic=scalar_only, native_history_advanced=False,
            history_float64_LE_sha256=vector_hash(history), line_search_classifications=dict(counts),
            iterations=iterations, final=final,
            failure_interpretation='A bounded unresolved search is not a proof that the unchanged stage has no root')

    def protocol(self, solves):
        raw = self.raw
        self.audit.require(raw['dimension'] == 2 and not raw['mapped'] and raw['upstream_transition'],
            'R2 stage analyzer supports the frozen serial flat 2D Upstream family')
        self.audit.require(raw['eta'] == 0 and raw['dt'] == .0025 and raw['steps'] == 0 and raw['time'] == 0.,
            'Restricted stage problem changed dt/anchor or advanced physical time')
        self.audit.require(raw['nonlinear_policy'] == 'ScalarFeasibleR2', 'Unexpected stage nonlinear policy')
        self.audit.require(self.size <= 3000 and all(record['maximum_iterations'] <= 12 for record in solves),
            'Restricted dense/Newton capacity exceeded its declared bounds')
        self.audit.require(raw['robust_parameters'] == [.4, 0., .2, .1, 0., 1., 1.], 'Changed frozen model/source parameters')
        root = ET.parse(self.session).getroot()
        from cls_epawf_m8_r2_cases import parameters
        params = parameters(root)
        self.audit.require(all(record['tolerance'] == float(params.get('RANSWALLNONLINEARTOLERANCE', '2e-12'))
            and record['maximum_iterations'] == int(params.get('RANSWALLMAXITERATIONS', '20')) for record in solves),
            'Restricted native tolerance/iteration cap differs from the actual XML/default')
        functions = {node.get('NAME', '') for node in root.findall('CONDITIONS/FUNCTION')}
        self.audit.require(not functions.intersection({'BodyForceInner', 'BodyForceOuter', 'ScalarSourceInner', 'ScalarSourceOuter'}),
            'Restricted stage contains manufactured/body/scalar forcing')
        self.audit.require(raw['stage_test_native_history_advanced'] is False, 'Auxiliary native history advanced')
        self.audit.near(raw['state'], raw['initial'], 'restricted_stage_kept_actual_initial_state', 0.)
        self.audit.require(raw['audit_times'] == [0.] and len(raw['accepted_states']) == 1, 'Auxiliary state became an accepted physical state')
        self.audit.near(raw['accepted_states'][0], raw['initial'], 'sole_native_accepted_initial_state', 0.)
        initial = raw['stage_test_initial_native_history']
        self.audit.require(initial['steps'] == 0 and initial['time'] == 0., 'Initial startup native snapshot is not step zero')
        self.audit.near(initial['state'], raw['state'][:len(initial['state'])], 'initial_native_differential_state', 0.)
        self.audit.near(initial['history'], raw['native_history'], 'all_native_history_slots_unchanged', 0.)
        self.audit.near(initial['times'], raw['native_history_times'], 'all_native_history_times_unchanged', 0.)
        levels = raw['stage_test_levels']; expected = [400.] if raw['stage_test'] == 'FrozenScalarR2' else LEVELS
        self.audit.near([level['alpha'] for level in levels], expected, 'preregistered_auxiliary_alpha_order', 0.)
        attempted = [level for level in levels if level['attempted']]
        self.audit.require(len(attempted) == len(solves), 'JSONL solves do not match attempted native stage levels')
        self.audit.near([row['alpha'] for row in solves], [level['alpha'] for level in attempted], 'trace_stage_alpha_correspondence', 0.)
        for record, level in zip(solves, attempted):
            self.audit.near(record['state'], level['state'], 'trace_final_state_matches_native_stage_report', 0.)
            self.audit.near(record['history'], level['history'], 'trace_history_matches_native_stage_report', 0.)
            self.audit.require(bool(record['converged']) == bool(level['converged']), 'Trace/stage convergence disagrees')
            self.audit.require(level['native_history_advanced'] is False, 'Stage level advanced native history')
            if level.get('residual') and record.get('residual'):
                self.audit.near(record['residual'], level['residual'], 'trace_stage_true_residual_identity', 0.)
        self.audit.require(attempted[-1]['alpha'] == 400., 'The actual alpha=400 target was not attempted')
        return dict(physical_target_time=raw['stage_test_target_time'], target_alpha=400.,
            native_time=raw['time'], native_steps=raw['steps'], native_history_advanced=False,
            all_native_solution_and_explicit_history_slots_unchanged=True,
            alpha_levels=[level['alpha'] for level in levels],
            history_float64_LE_sha256=vector_hash(raw['initial']),
            mode=raw['stage_test'], target_converged=bool(raw['stage_test_target_converged']),
            target_physical_range=bool(raw['stage_test_target_physical_range']),
            frozen_scalar_is_not_a_coupled_stage=raw['stage_test'] == 'FrozenScalarR2')

    def levels(self, solves):
        result = []
        for level in self.raw['stage_test_levels']:
            if not level['attempted']:
                result.append(dict(alpha=level['alpha'], attempted=False, reason=level['reason']))
                continue
            state, history, alpha = np.asarray(level['state']), np.asarray(level['history']), level['alpha']
            record = self.attempt('stage level '+str(alpha), lambda: self.state(state, history, alpha,
                level.get('residual'), level.get('constraints'), level['scalar_domain'], full=True, native_audit=level.get('audit')))
            # Repeated final state uses a cache, so check its newly exported
            # postqualification audit independently even when it is cached.
            if level.get('audit') and level.get('residual'):
                weak, unmasked, constraints = self.weak(state, history, alpha, level['residual'], level.get('constraints'))
                record['post_solve_physical'] = self.attempt('native postqualification audit '+str(alpha), lambda:
                    self.physical(state, history, alpha, level['audit'], unmasked, constraints, level['scalar_domain']))
            result.append(dict(alpha=alpha, attempted=True, native_converged=bool(level['converged']),
                native_true_coupled_residual_converged=bool(level['true_coupled_residual_converged']),
                native_physical_range_passed=bool(level['physical_range_passed']),
                residual_verification_failure=level.get('residual_verification_failure', ''),
                qualification_audit_failure=level.get('qualification_audit_failure', ''),
                nonlinear_failure=level.get('failure', ''), independent=record))
        return result


def read_jsonl(path):
    records = []
    with Path(path).open() as stream:
        for number, line in enumerate(stream, 1):
            if line.strip():
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError as error:
                    raise ValueError('Incomplete native solve trace line '+str(number)) from error
    if not records:
        raise ValueError('Missing native scalar-feasible solve trace')
    return records


def xml_structure(node):
    if node is None:
        return None
    return [node.tag, sorted(node.attrib.items()), ' '.join((node.text or '').split()),
        [xml_structure(child) for child in node]]


def verify_auxiliary(preflight, record, preflight_session, stage_session, source, matrix_root=None):
    """Analyze one completed native JSONL solve before the final case report.

    The actual stage XML must reproduce the source-bound preflight geometry,
    expansion definitions and protected physics, and its actual history must
    equal the preflight h0 byte for byte. Native field descriptors are reused
    from that preflight, not invented for an unfinished evolution report.
    """
    from cls_epawf_m8_r2_cases import normalized_physics, parameters, solver_info
    before, stage = (ET.parse(path).getroot() for path in (preflight_session, stage_session))
    if normalized_physics(before) != normalized_physics(stage):
        raise ValueError('Interim stage XML changed protected preflight physics')
    for name in ('GEOMETRY', 'EXPANSIONS'):
        if xml_structure(before.find(name)) != xml_structure(stage.find(name)):
            raise ValueError('Interim stage changed source-bound native '+name)
    info, params = solver_info(stage), parameters(stage)
    mode = info.get('RANSWALLSTAGETEST')
    if mode not in ('MatchedR2', 'FrozenScalarR2'):
        raise ValueError('Interim analysis requires the actual restricted-stage XML')
    if bool(record.get('frozen_velocity_pressure_diagnostic')) != (mode == 'FrozenScalarR2'):
        raise ValueError('Interim native record and actual stage XML disagree on the frozen-scalar mode')
    if mode == 'FrozenScalarR2' and record['alpha'] != 400.:
        raise ValueError('Frozen-scalar diagnosis is restricted to the original alpha=400 stage')
    if info.get('RANSWALLNONLINEARPOLICY') != 'ScalarFeasibleR2' or preflight['nonlinear_policy'] != 'ScalarFeasibleR2':
        raise ValueError('Interim stage/preflight did not select the required scalar-feasible policy')
    if float(params['TIMESTEP']) != .0025 or int(params['NUMSTEPS']) != 0:
        raise ValueError('Interim stage changed dt or requested physical advancement')
    if info['RANSWALLINNERSTRESS'] != preflight['inner_stress_model']:
        raise ValueError('Interim preflight and stage selected different inner laws')
    if list(map(int, (params['RANSWALLTANGENTIALORDER'], params['RANSWALLNORMALORDER']))) != preflight['quadrature_orders']:
        raise ValueError('Interim preflight/stage quadrature differs')
    if vector_hash(record['history']) != vector_hash(preflight['initial']):
        raise ValueError('Actual stage h0 differs from the source-bound preflight initial coefficients')
    if not np.array_equal(np.asarray(record['history']), np.asarray(preflight['initial'])):
        raise ValueError('Actual stage h0 is not byte-equivalent to the represented preflight state')
    if record['alpha'] not in LEVELS or record['time'] != .0025 or record['maximum_iterations'] > 12:
        raise ValueError('Interim record changed the preregistered stage time/alpha/capacity')
    if record['tolerance'] != float(params.get('RANSWALLNONLINEARTOLERANCE', '2e-12')):
        raise ValueError('Interim native tolerance differs from the actual stage XML/default')
    if record['maximum_iterations'] != int(params.get('RANSWALLMAXITERATIONS', '20')):
        raise ValueError('Interim native iteration cap differs from the actual stage XML/default')
    if record['inner_stress_model'] != preflight['inner_stress_model'] or record['policy'] != 'ScalarFeasibleR2:v1':
        raise ValueError('Interim record has a different inner law or nonlinear policy')
    view = dict(preflight, stage_test_target_time=record['time'])
    analyzer = Analyzer(view, stage_session, source, matrix_root)
    initial = analyzer.attempt('source-bound initial common inventory', lambda:
        analyzer.domain(view['initial'], view['initial_scalar_domain'], retain=True))
    result = analyzer.attempt('completed native auxiliary solve', lambda: analyzer.solve(record))
    return dict(schema=1, milestone='M8-R2 interim single-record analysis', success=not bool(analyzer.failures),
        analysis_success=not bool(analyzer.failures), checks=analyzer.audit.checks,
        errors=analyzer.audit.errors, failures=analyzer.failures, solved_record=result,
        selected_alpha=record['alpha'], physical_stage_time=record['time'],
        source_bound_h0_float64_LE_sha256=vector_hash(preflight['initial']),
        native_field_descriptors_sha256=canonical(preflight['fields']),
        protected_physics_sha256=canonical(normalized_physics(stage)),
        native_field_map_binding='Byte-identical preflight descriptors reused after identical source-bound geometry/expansions and actual h0 checks; final report comparison remains required',
        independent_common_static_coverage=analyzer.static_coverage, initial_scalar_domain=initial,
        dense_matrix_files=analyzer.matrix_files, current_action_gates=analyzer.gates,
        exact_affine_weak_continuity_assembly=getattr(analyzer, 'continuity_metadata', None),
        exact_affine_weak_continuity_gates=analyzer.continuity_gates,
        native_calls=0, only_one_completed_native_record_analyzed=True,
        full_native_case_protocol_assessed=False, final_native_history_slots_assessed=False,
        target_stage_status_assessed=False, physical_trajectory_qualification=False, M8_closed=False,
        qualification='Interim auxiliary-record implementation/finite-set diagnosis only. No final target result, no full-case native-history claim and no root-absence inference.')


def verify(raw, solves, session, source, matrix_root=None):
    """Importable analysis of retained data; no solver/subprocess is invoked."""
    analyzer = Analyzer(raw, session, source, matrix_root)
    protocol = analyzer.protocol(solves)
    initial = analyzer.attempt('initial actual common scalar domain', lambda:
        analyzer.domain(raw['initial'], raw['initial_scalar_domain'], retain=True))
    results = [analyzer.attempt('solve '+str(index), lambda record=record: analyzer.solve(record))
        for index, record in enumerate(solves)]
    levels = analyzer.levels(solves)
    gates = analyzer.gates
    action_passed = bool(gates) and all(gate['available'] and gate['gate_met'] for gate in gates)
    return dict(schema=1, milestone='M8-R2 independent stage/path analysis',
        success=not bool(analyzer.failures), analysis_success=not bool(analyzer.failures),
        checks=analyzer.audit.checks, errors=analyzer.audit.errors, failures=analyzer.failures,
        protocol=protocol, initial_scalar_domain=initial, solves=results, stage_levels=levels,
        independent_common_static_coverage=analyzer.static_coverage,
        exact_affine_weak_continuity_assembly=getattr(analyzer, 'continuity_metadata', None),
        exact_affine_weak_continuity_gates=analyzer.continuity_gates,
        exact_affine_weak_continuity_gate_passed=bool(analyzer.continuity_gates) and
            all(gate['gate_met'] for gate in analyzer.continuity_gates),
        dense_matrix_files=analyzer.matrix_files, current_action_gates=gates,
        fixed_action_tolerance=ACTION_TOLERANCE, current_action_gate_passed=action_passed,
        target_converged=bool(raw['stage_test_target_converged']),
        target_physical_range_passed=bool(raw['stage_test_target_physical_range']),
        rigorous_discrete_stability=False, M8_closed=False,
        qualification='Read-only implementation/finite-set and bounded nonlinear diagnosis. Root convergence, physical range, action accuracy and rigorous stability are separate gates.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--trace', type=Path, required=True)
    parser.add_argument('--session', type=Path, required=True)
    parser.add_argument('--runtime-pin', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source', type=Path, default=Path(__file__).resolve().parents[4])
    parser.add_argument('--expected-commit')
    parser.add_argument('--matrix-root', type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Scientific retained-stage analysis must run on Linux')
    args.output.mkdir(parents=True, exist_ok=False)
    source = args.source.resolve()
    inputs = [args.report.resolve(), args.trace.resolve(), args.session.resolve(), args.runtime_pin.resolve()]
    receipt = dict(schema=1, started_utc=stamp(), success=False, native_calls=0,
        driver_sha256=sha(__file__), source=str(source), source_commit=git(source, 'rev-parse', 'HEAD'),
        source_clean_start=not bool(git(source, 'status', '--porcelain')), numerical_threads=THREADS,
        byteorder=sys.byteorder, python=sys.version, cpu_affinity=sorted(os.sched_getaffinity(0)),
        inputs={str(path):dict(bytes=path.stat().st_size, sha256=sha(path)) for path in inputs}, failures=[])
    receipt['imported_source_files_start'] = {str(Path(module.__file__).resolve()):sha(module.__file__)
        for module in list(sys.modules.values()) if getattr(module, '__file__', None)
        and str(Path(module.__file__).resolve()).startswith(str(source))}
    save(args.output/'analysis-receipt.json', receipt)
    started = time.monotonic()
    try:
        if not receipt['source_clean_start'] or (args.expected_commit and receipt['source_commit'] != args.expected_commit):
            raise ValueError('Analyzer source is not the requested clean pinned revision')
        pin = json.loads(args.runtime_pin.read_text())
        if not pin['success'] or not pin['clean']:
            raise ValueError('Producing runtime fingerprint failed or source was dirty')
        receipt['producing_source_commit'] = pin['source_commit']
        producing_source = Path(pin['source_root'])
        if git(producing_source, 'rev-parse', 'HEAD') != pin['source_commit'] or git(producing_source, 'status', '--porcelain'):
            raise ValueError('Retained producing source no longer matches its clean runtime pin')
        solver_keys = [key for key in pin['runtime_artifacts'] if Path(key).name == 'IncNavierStokesSolver']
        if len(solver_keys) != 1:
            raise ValueError('Producing runtime has no unique native evolution image')
        solver_key = solver_keys[0]
        solver = Path(pin['runtime_resolved_paths'][solver_key])
        if sha(solver) != pin['runtime_artifacts'][solver_key]:
            raise ValueError('Retained producing solver changed after the runtime pin')
        receipt['producing_solver'] = dict(path=str(solver), sha256=sha(solver), bytes=solver.stat().st_size)
        receipt['producing_runtime_sections'] = {key:canonical(pin[key]) for key in
            ('runtime_artifacts', 'runtime_bindings', 'runtime_environment', 'python_packages', 'tool_images', 'source_inputs')}
        raw = json.loads(args.report.read_text()); solves = read_jsonl(args.trace)
        result = verify(raw, solves, args.session, source, args.matrix_root)
        for path, item in result['dense_matrix_files'].items():
            receipt['inputs'][path] = item
        receipt['analysis_success'] = result['analysis_success']; receipt['success'] = result['analysis_success']
        result['input_receipt'] = str(args.output/'analysis-receipt.json')
        save(args.output/'validation.json', result)
    except Exception as error:
        receipt['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
    for path, expected in receipt['inputs'].items():
        if sha(path) != expected['sha256'] or Path(path).stat().st_size != expected['bytes']:
            receipt['success'] = False
            receipt['failures'].append(dict(type='ProvenanceError', message='Retained analysis input changed: '+path))
    receipt['source_commit_end'] = git(source, 'rev-parse', 'HEAD')
    receipt['source_clean_end'] = not bool(git(source, 'status', '--porcelain'))
    if receipt['source_commit_end'] != receipt['source_commit'] or not receipt['source_clean_end']:
        receipt['success'] = False
        receipt['failures'].append(dict(type='ProvenanceError', message='Analyzer source changed during read-only analysis'))
    receipt['imported_source_files'] = {str(Path(module.__file__).resolve()):sha(module.__file__)
        for module in list(sys.modules.values()) if getattr(module, '__file__', None)
        and str(Path(module.__file__).resolve()).startswith(str(source))}
    for path, expected in receipt['imported_source_files_start'].items():
        if receipt['imported_source_files'].get(path) != expected:
            receipt['success'] = False
            receipt['failures'].append(dict(type='ProvenanceError', message='Imported analysis source changed: '+path))
    receipt['finished_utc'] = stamp(); receipt['seconds'] = time.monotonic()-started
    save(args.output/'analysis-receipt.json', receipt)
    print(json.dumps(dict(success=receipt['success'], output=str(args.output/'validation.json'))))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
