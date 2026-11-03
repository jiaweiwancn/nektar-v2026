#!/usr/bin/env python3
"""Read-only M8-R1 accepted-prefix reconstruction; never launches a solver.

Reconstructs the archived native finite evaluation sets from independent
Jacobi rows, dimensional AWF profiles and indexed CLS contractions.  Range
membership is distinct from a counterexample, and a Newton stop without a
saved rejected iterate is not a physical inadmissibility witness.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import traceback
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
    'MKL_NUM_THREADS', 'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np

BASELINE = 'ab5c3c6fee19761e8d88692e8e6c331e5889234e'


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                    allow_nan=False).encode()).hexdigest()


def save(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n')
    temporary.replace(path)


def git(source, *args):
    return subprocess.check_output(['git', '-c', 'gc.auto=0', '-C', str(source), *args], text=True).strip()


class Reconstruction:
    def __init__(self, raw, source, xml):
        sys.path.insert(0, str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
        from cls_epawf_m7_energy_oracle import _fields, _profile, _forcing, _horizontal
        from cls_epawf_m8_native_oracle import normal_rule, selected_stress
        from cls_epawf_m6_operator_oracle import _segment_projection
        from validate_cls_epawf_m1 import indexed
        self.raw, self.fields = raw, _fields(raw)
        self.profile, self.forcing, self.horizontal = _profile, _forcing, _horizontal
        self.normal_rule, self.stress = normal_rule, selected_stress
        self.segment, self.indexed = _segment_projection, indexed
        self.d, self.off = raw['dimension'], raw['offsets']
        if self.d != 2 or raw.get('mapped'):
            raise ValueError('This archived diagnostic is explicitly flat two-dimensional')
        if raw['eta'] != 0 or not raw['upstream_transition']:
            raise ValueError('Wrong archived unanchored transition family')
        functions = [node.attrib['NAME'] for node in ET.fromstring(xml).findall('./CONDITIONS/FUNCTION')]
        if any(name.startswith(('BodyForce', 'ScalarSource')) for name in functions):
            raise ValueError('This independent weak operator requires the archived unforced case')
        self.nu, self.H, self.h, self.rK, _, self.K0, self.Q, self.kappa = raw['parameters']
        self.us = raw['ustar']
        self.mass = [field.mass() for field in self.fields]
        self.constants = [field.constant() for field in self.fields]
        self.reference_constraints = {int(row): float(value) for row, value in raw['reference_constraints']}
        self.points_native = []
        expansion = ET.fromstring(xml).find('./EXPANSIONS/E')
        point_counts = [int(value) for value in expansion.attrib.get('NUMPOINTS', '8,8').split(',')]
        if len(point_counts) == 1:
            point_counts *= 2
        self.native_counts = point_counts
        self.errors, self.checks = {}, 0

    def near(self, actual, expected, name, tolerance=2e-9):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        error = float(np.max(np.abs(a-b), initial=0.))/max(1., float(np.max(np.abs(b), initial=0.)))
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.checks += 1
        if a.shape != b.shape or not np.all(np.isfinite(a)) or error > tolerance:
            raise ValueError(f'{name}: scaled error {error} exceeds {tolerance}; actual={a.tolist()}, expected={b.tolist()}')

    def evaluate(self, block, points, side=0, eid=None):
        points = np.atleast_2d(points)
        field = self.fields[block]
        value, gradient = np.empty(len(points)), np.empty((len(points), 3))
        groups = {}
        for i, point in enumerate(points):
            element = field.by_id[eid] if eid is not None else field.element(point, side)
            groups.setdefault(element['geometry'], []).append(i)
        for geometry_id, indices in groups.items():
            ids, phi, derivative = field.rows(field.by_id[geometry_id], points[indices])
            coefficients = self.state[self.off[block]+ids]
            value[indices] = phi@coefficients
            gradient[indices] = np.einsum('qaj,j->qa', derivative, coefficients)
        return value, gradient

    def velocity(self, points, side=0, eid=None):
        points = np.atleast_2d(points)
        u, g = np.zeros((len(points), 3)), np.zeros((len(points), 3, 3))
        for block, axis in enumerate((0, 2)):
            u[:, axis], g[:, axis] = self.evaluate(block, points, side, eid)
        return u, g

    def ray(self, point, baseline=False, seed=False):
        key = (float(point[0]), float(point[1]), baseline, seed)
        if key in self.rays:
            return self.rays[key]
        chart = np.asarray([point[0], point[1], self.H])
        raw = dict(self.raw, upstream_transition=False, chi=0.) if baseline else self.raw
        state = self.state
        if seed:
            # Upstream constructor seed quadrature precedes transition
            # readiness; its default profile uses the baseline weight.
            raw = dict(self.raw, upstream_transition=False, chi=0.)
            state = np.zeros_like(self.state)
            state[self.off[self.d]:self.off[self.d+1]] = self.K0*self.constants[self.d]
        profile, metric = self.profile(raw, self.fields, state, chart)
        r, dr = self.normal_rule(profile, raw, self.fields, metric)
        points = np.column_stack([np.full_like(r, chart[0]), np.zeros_like(r), r])
        k, _ = self.evaluate(self.d, points)
        values = np.asarray(profile.fields(r), float)
        pre_e, pre_d = values[2], values[2]-values[3]
        measure = dr*np.asarray(profile.jacobian(r), float)
        theta, delta, smooth = raw['robust_parameters'][2:5]
        t = np.minimum(1., k/(theta*float(profile.k)))
        if np.any(t < 0):
            raise ValueError('An accepted coefficient state has negative source-point k')
        limiter = 10*t**3-15*t**4+6*t**5 if smooth else t
        budget = float(measure@pre_d)
        b = float(measure@(pre_d*limiter))/budget
        m = max(b, delta)
        lam = b/m
        answer = dict(profile=profile, metric=metric, points=points, r=r, measure=measure,
                      k=k, pre_e=pre_e, pre_d=pre_d, limiter=limiter, b=b, m=m, lam=lam,
                      e=lam*pre_e, destruction=pre_d*limiter/m)
        self.rays[key] = answer
        return answer

    def add(self, block, points, weights, source, flux=None):
        field = self.fields[block]
        groups = {}
        for i, point in enumerate(points):
            groups.setdefault(field.element(point)['geometry'], []).append(i)
        for eid, positions in groups.items():
            ids, phi, derivative = field.rows(field.by_id[eid], points[positions])
            w = weights[positions]
            load = phi.T@(w*np.asarray(source)[positions])
            if flux is not None:
                load += np.einsum('qaj,qa,q->j', derivative, np.asarray(flux)[positions], w)
            self.loads[self.off[block]+ids] += load

    def outer_ranges(self, points, k, epsilon, G):
        G0 = self.us/(self.kappa*points[:, 2])
        center = np.zeros((3, 3)); center[0, 2] = 1.
        entry = np.max(np.abs(G/G0[:, None, None]-center), axis=(1, 2))
        kr = k/self.K0
        er = epsilon/(self.us**3/(self.kappa*points[:, 2]))
        tr = np.trace(G, axis1=1, axis2=2)/G0
        qs = k/epsilon*np.sqrt(np.sum((G+G.transpose(0, 2, 1))**2, axis=(1, 2))/2)
        qo = k/epsilon*np.sqrt(np.sum((G-G.transpose(0, 2, 1))**2, axis=(1, 2))/2)
        return dict(gradient_entry_departure=entry, k_over_K0=kr, epsilon_over_epsilon0=er,
                    normalized_divergence=tr, q=np.maximum(qs, qo), q_strain=qs, q_rotation=qo)

    def sample(self, name, points, inner, constitutive=False, ray=None, side=0, eid=None):
        points = np.atleast_2d(points)
        k, _ = self.evaluate(self.d, points, side, eid)
        epsilon = None
        if not inner:
            epsilon = self.evaluate(self.d+1, points, side)[0]
        elif constitutive:
            p = ray['profile']
            epsilon = ray['lam']*np.asarray(p.fields(points[:, 2]), float)[2]
        record = self.sets.setdefault(name, dict(count=0, minimum_k=math.inf, minimum_epsilon=None,
            minimum_covariance=None, minimum_yx_tangent=None, invalid_count=0, negative_covariance_count=0))
        record['count'] += len(points)
        record['minimum_k'] = min(record['minimum_k'], float(np.min(k)))
        if epsilon is not None:
            minimum = float(np.min(epsilon))
            record['minimum_epsilon'] = minimum if record['minimum_epsilon'] is None else min(record['minimum_epsilon'], minimum)
        invalid = ~np.isfinite(k) | (k < 0 if inner else k <= 0)
        if epsilon is not None:
            invalid |= ~np.isfinite(epsilon) | (epsilon <= 0)
        record['invalid_count'] += int(np.count_nonzero(invalid))
        if not constitutive:
            return
        _, G = self.velocity(points, side, eid)
        S = G+G.transpose(0, 2, 1)
        if inner:
            profile = ray['profile']
            nt = float(profile.a)*np.maximum(points[:, 2]-float(profile.zv), 0.)
            tensor = 2*k[:, None, None]/3*np.eye(3)-nt[:, None, None]*S
            rR = 3*nt*np.linalg.eigvalsh(S)[:, -1]/(2*k)
        else:
            tensor = np.asarray([self.stress(self.raw, g, kk, ep, p)[0]
                                 for g, kk, ep, p in zip(G, k, epsilon, points)])
            tangents = []
            for g, kk, ep, p in zip(G, k, epsilon, points):
                step = 1e-5*max(1., abs(g[0, 2])); plus, minus = g.copy(), g.copy()
                plus[1, 0] += step; minus[1, 0] -= step
                tangents.append(-(self.stress(self.raw, plus, kk, ep, p)[0][1, 0]-
                                  self.stress(self.raw, minus, kk, ep, p)[0][1, 0])/(2*step))
            minimum = float(min(tangents))
            record['minimum_yx_tangent'] = minimum if record['minimum_yx_tangent'] is None else min(record['minimum_yx_tangent'], minimum)
            ranges = self.outer_ranges(points, k, epsilon, G)
            for radius in (.001, .01):
                # These are the previously certified local boxes. The .01
                # box needs GammaScale=3 for the joint lower bound; the actual
                # archived GammaScale=1 does not meet that momentum certificate.
                inside = ((ranges['gradient_entry_departure'] <= radius) &
                          (ranges['k_over_K0'] >= .99) & (ranges['k_over_K0'] <= 1.01) &
                          (ranges['epsilon_over_epsilon0'] >= .99) & (ranges['epsilon_over_epsilon0'] <= 1.01))
                row = self.range_records.setdefault(str(radius), dict(count=0, outside_count=0, first_outside=None,
                    maximum_gradient_departure=0., k_ratio=[math.inf, -math.inf], epsilon_ratio=[math.inf, -math.inf],
                    maximum_normalized_divergence=0., maximum_q=0.))
                row['count'] += len(points); row['outside_count'] += int(np.count_nonzero(~inside))
                row['maximum_gradient_departure'] = max(row['maximum_gradient_departure'], float(np.max(ranges['gradient_entry_departure'])))
                row['maximum_normalized_divergence'] = max(row['maximum_normalized_divergence'], float(np.max(np.abs(ranges['normalized_divergence']))))
                row['maximum_q'] = max(row['maximum_q'], float(np.max(ranges['q'])))
                for key, value in (('k_ratio', ranges['k_over_K0']), ('epsilon_ratio', ranges['epsilon_over_epsilon0'])):
                    row[key] = [min(row[key][0], float(np.min(value))), max(row[key][1], float(np.max(value)))]
                if row['first_outside'] is None and np.any(~inside):
                    i = int(np.flatnonzero(~inside)[0])
                    row['first_outside'] = dict(set=name, coordinate=points[i].tolist(),
                        **{key: float(value[i]) for key, value in ranges.items()})
        eigenvalues, vectors = np.linalg.eigh(tensor)
        values = eigenvalues[:, 0]
        record['negative_covariance_count'] += int(np.count_nonzero(values < 0))
        i = int(np.argmin(values))
        if record['minimum_covariance'] is None or values[i] < record['minimum_covariance']:
            vector = vectors[i, :, 0]
            if vector[np.argmax(np.abs(vector))] < 0:
                vector = -vector
            point = points[i]
            detail = dict(coordinate=point.tolist(), eigenvalue=float(values[i]), eigenvector=vector.tolist(),
                eigenvalues=eigenvalues[i].tolist(), covariance=tensor[i].tolist(), k=float(k[i]),
                epsilon=float(epsilon[i]), G=G[i].tolist(), S=S[i].tolist())
            if inner:
                reference_state = self.state
                self.state = np.asarray(self.raw['initial'])
                _, initial_G = self.velocity([point], side)
                self.state = reference_state
                seed_profile = self.ray(point, seed=True)['profile']
                analytic_nt = float(seed_profile.a)*max(point[2]-float(seed_profile.zv), 0.)
                detail.update(sampled_K=float(profile.k), k_over_K=float(k[i]/float(profile.k)),
                    nu_A=float(nt[i]), r_R=float(rR[i]), maximum_strain_eigenvalue=float(np.linalg.eigvalsh(S[i])[-1]),
                    initial_native_G=initial_G[0].tolist(),
                    initial_native_strain_norm=float(np.linalg.norm(initial_G[0]+initial_G[0].T)),
                    strain_norm=float(np.linalg.norm(S[i])),
                    analytic_neutral_shear=float(self.us**2/(self.nu+analytic_nt)),
                    depletion_relative_K0=float(k[i]/self.K0), b=ray['b'], lambda_value=ray['lam'], chi=float(profile.chi))
            else:
                detail.update(**{key: float(value[i]) for key, value in ranges.items()})
            record['minimum_covariance'] = float(values[i]); record['worst_point'] = detail

    def trace(self):
        field, expected, baseline_ids, data = self.fields[self.d+1], {}, set(), []
        baseline_start = self.raw['transition_controls'][4]
        parts = [{}, {}]
        for element in field.elements:
            if abs(element['lower'][2]-self.H) > 1e-12:
                continue
            nodes, transform = self.segment(element['degree'])
            xs = element['lower'][0]+(nodes+1)*(element['upper'][0]-element['lower'][0])/2
            values = [[], []]
            for x in xs:
                point = np.array([x, 0., self.H]); current, baseline = self.ray(point), self.ray(point, True)
                e = current['lam']*float(np.asarray(current['profile'].fields(np.array([self.H])), float)[2, 0])
                eb = baseline['lam']*float(np.asarray(baseline['profile'].fields(np.array([self.H])), float)[2, 0])
                values[0].append(eb); values[1].append(e-eb)
                data.append(dict(x=float(x), K=float(current['profile'].k), chi=float(current['profile'].chi),
                    b=current['b'], lambda_value=current['lam'], baseline_lambda=baseline['lam'], analytic_trace=e))
            count = element['degree']+1
            for table, value in zip(parts, values):
                for index, coefficient in zip(element['ids'][:count], np.asarray(element['sign'][:count])*(transform@value)):
                    if index in table:
                        self.near(coefficient, table[index], 'shared_trace_endpoints')
                    table[int(index)] = float(coefficient)
            if element['lower'][0] >= baseline_start-1e-11:
                baseline_ids.update(map(int, element['ids'][:count]))
        for index in baseline_ids:
            parts[1][index] = 0.
        for index in parts[0]:
            expected[self.off[self.d+1]+index] = parts[0][index]+parts[1][index]
        return expected, data

    def state_record(self, state, history, alpha, time, native):
        self.state, self.rays = np.asarray(state, float), {}
        self.sets, self.range_records, self.loads = {}, {}, np.zeros(len(state))
        robust, wall, modification = np.zeros(8), np.zeros(5), np.zeros(3)
        robust[6] = 1.
        direct_energy, advection, rough_work = 0., 0., 0.
        ray_records = []
        for chart, horizontal in self.horizontal(self.fields, self.d, self.H, self.raw['quadrature_orders'][0]):
            ray = self.ray(chart); points, p = ray['points'], ray['profile']
            weights = horizontal*ray['measure']
            k, gk = self.evaluate(self.d, points); U, G = self.velocity(points)
            pi = self.evaluate(self.d+2, points)[0]
            S, divergence = G+G.transpose(0, 2, 1), np.trace(G, axis1=1, axis2=2)
            nt = float(p.a)*np.maximum(ray['r']-float(p.zv), 0.)
            R = 2*k[:, None, None]/3*np.eye(3)-nt[:, None, None]*S
            production = -np.einsum('qab,qab->q', R, G)
            stress = (self.nu+nt)[:, None, None]*S-pi[:, None, None]*np.eye(3)
            force = self.forcing(self.raw, self.fields, self.state, chart)
            for block, axis in enumerate((0, 2)):
                self.add(block, points, weights, -(ray['r'] <= self.h).astype(float)*force[axis], stress[:, axis]-U[:, axis, None]*U)
            self.add(self.d, points, weights, -(production-ray['destruction']), (self.nu+nt)[:, None]*gk-U*k[:, None])
            self.add(self.d+2, points, weights, divergence)
            reduced = self.nu*np.sum(S*S, axis=(1, 2))/2+(2*k/3-pi)*divergence+ray['destruction']
            direct_energy += float(weights@reduced)
            advection -= float(weights@np.einsum('qab,qa,qb->q', G, U, U))
            rough_work -= float(weights[(ray['r'] <= self.h)]@(U[(ray['r'] <= self.h)]@force))
            robust[:4] += [float(weights@ray['pre_e']), float(weights@ray['e']),
                float(weights@ray['destruction']), float(weights@(ray['e']-ray['destruction']))]
            robust[4] += float(np.sum(weights[ray['limiter'] < 1]))
            robust[5] += float(np.sum(weights)) if ray['b'] < self.raw['robust_parameters'][3] else 0.
            robust[6] = min(robust[6], ray['lam'])
            robust[7] = max(robust[7], float(np.max(ray['destruction']/k)))
            ea = float(p.k)**1.5/(2.55*np.maximum(ray['r'], float(p.ze)))
            wall += [float(weights@ea), float(weights@ray['e']), float(weights@(ray['e']-ray['destruction'])),
                     float(weights@production), float(weights@(self.nu*np.sum(S*S, axis=(1, 2))/2))]
            ray_records.append(dict(x=float(chart[0]), K=float(p.k), chi=float(p.chi), b=ray['b'],
                m=ray['m'], lambda_value=ray['lam'], cutoff=float(p.rc) if p.chi > 0 else 0.,
                source_points=len(points), minimum_k_over_K=float(np.min(k)/float(p.k))))
            self.sample('source_volume_inner', points, True, True, ray)
            self.sample('physical_sample_K', [chart[0], 0., self.rK], True, True, ray)
            for height in self.raw['transition_controls'][2:4]:
                self.sample('physical_sample_outer', [chart[0], 0., height], False, True)
            seed = self.ray(chart, seed=True)
            self.sample('prescribed_source_inner', seed['points'], True)
        for element in self.fields[self.d+3].elements:
            points, weights = self.fields[self.d+3].quadrature(element, self.raw['quadrature_orders'][0])
            k, gk = self.evaluate(self.d, points); eps, ge = self.evaluate(self.d+1, points)
            po = self.evaluate(self.d+3, points)[0]; U, G = self.velocity(points)
            tensor, nuts, original = [], [], []
            for g, kk, e, q in zip(G, k, eps, points):
                R, nt = self.stress(self.raw, g, kk, e, q)
                tensor.append(R); nuts.append(nt); original.append(self.indexed(g, kk, e, 'physical')['stress'].real)
            tensor, nuts, original = np.asarray(tensor), np.asarray(nuts), np.asarray(original)
            production = -np.einsum('qab,qab->q', tensor, G)
            stress = -tensor+(2*k/3-po)[:, None, None]*np.eye(3)
            for block, axis in enumerate((0, 2)):
                self.add(block, points, weights, np.zeros(len(points)), stress[:, axis]-U[:, axis, None]*U)
            self.add(self.d, points, weights, -(production-eps), nuts[:, None]*gk-U*k[:, None])
            ratio = k**1.5/(2.5*eps*points[:, 2])
            yap = eps*eps/k*.83*np.maximum(ratio-1, 0)*ratio**2
            self.add(self.d+1, points, weights, -(1.44*eps/k*production-1.92*eps*eps/k+yap),
                     (nuts/1.3)[:, None]*ge-U*eps[:, None])
            divergence = np.trace(G, axis1=1, axis2=2)
            self.add(self.d+3, points, weights, divergence)
            direct_energy += float(weights@((2*k/3-po)*divergence+eps))
            advection -= float(weights@np.einsum('qab,qa,qb->q', G, U, U))
            change = -np.einsum('qab,qab->q', tensor-original, G)
            modification[:2] += float(weights@change)
            self.sample('source_volume_outer', points, False, True)
        strong = np.zeros(5)
        for element in self.fields[0].elements:
            points, weights = self.fields[0].quadrature(element, max(8, element['degree']+3))
            _, G = self.velocity(points, eid=element['geometry'])
            div = np.trace(G, axis1=1, axis2=2); region = int(element['lower'][2] >= self.H)
            strong[region] += float(weights@(div*div)); strong[region+2] += float(weights@div)
            self.sample('mass_volume', points, not region, eid=element['geometry'])
            for dense in (False, True):
                counts = [2*element['degree']+3]*2 if dense else self.native_counts
                axes = []
                for axis, count in zip((0, 2), counts):
                    z = np.linspace(-1., 1., count) if dense else np.r_[-1., np.polynomial.legendre.Legendre.basis(count-1).deriv().roots(), 1.]
                    axes.append(element['lower'][axis]+(z+1)*(element['upper'][axis]-element['lower'][axis])/2)
                points = np.array([[x, 0., z] for z in axes[1] for x in axes[0]])
                self.sample('reconstruction_dense' if dense else 'reconstruction_native', points,
                            not region, eid=element['geometry'])
        strong[:2] = np.sqrt(strong[:2])
        for axis in (0, 2):
            tangent = 2 if axis == 0 else 0
            z, w = np.polynomial.legendre.leggauss(8)
            for side in (0, 1):
                bound = self.fields[0].upper[axis] if side else self.fields[0].lower[axis]
                for element in self.fields[0].elements:
                    if abs((element['upper'] if side else element['lower'])[axis]-bound) > 1e-11:
                        continue
                    points = np.zeros((len(z), 3)); points[:, axis] = bound
                    width = element['upper'][tangent]-element['lower'][tangent]
                    points[:, tangent] = element['lower'][tangent]+(z+1)*width/2
                    U, _ = self.velocity(points, eid=element['geometry'])
                    strong[4] += float((w*width/2)@U[:, axis])*(1 if side else -1)
        face = np.zeros(8); tangential_square = 0.; interface_records = []
        for point, horizontal in self.horizontal(self.fields, self.d, self.H, 8):
            ray = self.ray(point); p = ray['profile']
            self.sample('interface_inner', point, True, True, ray, side=-1)
            self.sample('interface_outer', point, False, True, side=1)
            U, Gi = self.velocity([point], -1); _, Go = self.velocity([point], 1)
            k, gki = self.evaluate(self.d, [point], -1); _, gko = self.evaluate(self.d, [point], 1)
            eps = self.evaluate(self.d+1, [point], 1)[0][0]
            pi = self.evaluate(self.d+2, [point], -1)[0][0]; po = self.evaluate(self.d+3, [point], 1)[0][0]
            nt = float(p.a)*max(self.H-float(p.zv), 0.)
            R, nto = self.stress(self.raw, Go[0], k[0], eps, point)
            ti = (self.nu+nt)*(Gi[0]+Gi[0].T)-pi*np.eye(3)
            to = -R+(2*k[0]/3-po)*np.eye(3)
            jump, normal = (ti-to)[:, 2], (ti-to)[2, 2]
            kflux = (self.nu+nt)*gki[0, 2]-nto*gko[0, 2]
            face += horizontal*np.array([1., pi-po, normal**2, jump@jump, 0., kflux**2, U[0, 2], normal])
            face[4] = max(face[4], abs(normal)); tangential_square += horizontal*(jump[0]**2+jump[1]**2)
            interface_records.append(dict(x=float(point[0]), pi=float(pi), po=float(po),
                mechanical_pressure_jump=float(pi-po), inner_traction=ti[:, 2].tolist(),
                outer_traction=to[:, 2].tolist(), full_traction_jump=jump.tolist(),
                normal_jump=float(normal), tangential_jump=jump[:2].tolist(), k_flux_jump=float(kflux)))
        for i in (2, 3, 5):
            face[i] = math.sqrt(max(0., face[i])/face[0])
        storage = stored = 0.
        for block in range(self.d+2):
            ids = np.arange(self.fields[block].size)+self.off[block]
            delta = self.state[ids]-np.asarray(history)[ids]
            self.loads[ids] += alpha*(self.mass[block]@delta)
            if block <= self.d:
                test = self.constants[block] if block == self.d else self.state[ids]
                storage += alpha*float(test@self.mass[block]@delta)
                stored += (1. if block == self.d else .5)*float(test@self.mass[block]@self.state[ids])
        pressure_mean = 0.
        for block in (self.d+2, self.d+3):
            ids = np.arange(self.fields[block].size)+self.off[block]
            mean = self.mass[block]@self.constants[block]/self.fields[0].volume
            self.loads[ids] += self.state[-1]*mean
            pressure_mean += float(mean@self.state[ids])
        self.loads[-1] = pressure_mean
        constraints, trace_records = self.trace()
        all_constraints = dict(self.reference_constraints); all_constraints.update(constraints)
        residual = self.loads.copy()
        for row, value in all_constraints.items():
            residual[row] = self.state[row]-value
        norms = []
        for block, field in enumerate(self.fields):
            free = np.array([i for i in range(field.size) if self.off[block]+i not in all_constraints], int)
            load = residual[self.off[block]+free]
            square = float(load@np.linalg.solve(self.mass[block][np.ix_(free, free)], load)) if len(free) else 0.
            scale = self.us if block < self.d else self.K0 if block == self.d else self.us**3/(self.kappa*self.H) if block == self.d+1 else 1.
            rate = self.us/self.H if block >= self.d+2 else scale*self.us/self.H
            norms.append(math.sqrt(max(0., square))/(rate*math.sqrt(field.volume)))
        norms.append(abs(pressure_mean)/(self.us*self.us))
        bc = 0.
        for row in all_constraints:
            block = next(i for i in range(len(self.off)-1) if self.off[i] <= row < self.off[i+1])
            scale = self.us if block < self.d else self.K0 if block == self.d else self.us**3/(self.kappa*self.H)
            bc = max(bc, abs(residual[row])/scale)
        norms.append(bc)
        unmasked = reaction = 0.
        balance = np.zeros(2*(self.d+1))
        for block in range(self.d+1):
            ids = np.arange(self.fields[block].size)+self.off[block]
            test = self.constants[block] if block == self.d else self.state[ids]
            constrained = np.array([row in all_constraints for row in ids])
            unmasked += float(test@self.loads[ids]); reaction += float(test[constrained]@self.loads[ids[constrained]])
            balance[block] = float(self.constants[block]@self.loads[ids])
            balance[self.d+1+block] = float(self.constants[block][constrained]@self.loads[ids[constrained]])
        energy = [unmasked, storage, direct_energy, advection, rough_work, 0.,
                  storage+direct_energy+advection+rough_work, unmasked-storage-direct_energy-advection-rough_work,
                  stored, reaction, 0., unmasked-reaction]
        for key, value in self.sets.items():
            expected = [value['count'], value['minimum_k'], value['minimum_epsilon'] or 0.,
                        value['minimum_covariance'] or 0., value['minimum_yx_tangent'] or 0., value['invalid_count']]
            self.near(expected, native[key], 'finite_set_'+key)
        self.near(robust, native['robust_source'], 'independent_robust_source')
        self.near(wall, native['wall_sources'], 'independent_wall_sources')
        self.near(modification, native['constitutive_modification'], 'independent_constitutive_modification')
        self.near(strong, native['strong_continuity'], 'independent_strong_continuity')
        self.near(face, native['interface'], 'independent_full_traction_kflux')
        self.near(norms, native['coupled_norms'], 'independent_complete_weak_scaled_norms')
        self.near(balance, native['momentum_k_load_and_boundary_reaction'], 'independent_boundary_reactions')
        self.near(energy, native['energy'], 'independent_complete_energy_accounting')
        return dict(time=time, alpha=alpha, state_sha256=canonical(state), history_sha256=canonical(history),
            sets=self.sets, constitutive_range_membership=self.range_records,
            robust_source_native=native['robust_source'], robust_source_independent=robust.tolist(),
            source_rays=ray_records, matching_trace=trace_records,
            pressure=dict(gauge_multiplier=float(self.state[-1]), common_pressure_mean=pressure_mean,
                scaled_pressure_rows=norms[self.d+2:self.d+4], scaled_gauge_residual=norms[-2],
                interpretation='One common pressure gauge; mechanical pressures need not be continuous'),
            continuity=dict(strong=strong.tolist(), weak_pressure_rows=norms[self.d+2:self.d+4]),
            interface=dict(native_layout=face.tolist(), tangential_rms=math.sqrt(tangential_square/face[0]),
                points=interface_records, interpretation='Complete traction and diffusive k flux, not mechanical pressure matching'),
            constraints=dict(count=len(all_constraints), scaled_maximum=bc, trace_rows=len(constraints)),
            energy=dict(independent=energy, native=native['energy'], boundary_reaction=reaction, anchor=0.),
            complete_scaled_weak_norms=norms, coefficient_load_sha256=canonical(self.loads.tolist()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--archive', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise SystemExit('Numerical reconstruction must run on the designated Linux host')
    args.output.mkdir(parents=True, exist_ok=False)
    report = dict(schema=1, milestone='M8-R1', phase='A archived accepted prefix', started_utc=stamp(),
        success=False, trajectory_launched=False, numerical_threads=THREADS, source_commit=git(args.source, 'rev-parse', 'HEAD'),
        source_clean_start=not bool(git(args.source, 'status', '--porcelain')), source=str(args.source),
        archive=str(args.archive), driver_sha256=sha(__file__), states=[], events=[], failures=[])
    output = args.output/'archive-ledger.json'
    try:
        if report['source_commit'] != BASELINE or not report['source_clean_start']:
            raise ValueError('Archive reconstruction source is not the clean reviewed baseline')
        inputs = {}
        for name in ('case.evolution.json', 'case.xml', 'run.log'):
            path = args.archive/name
            data = path.read_bytes()  # complete reads precede all reconstruction
            inputs[name] = dict(path=str(path), bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
        raw = json.loads((args.archive/'case.evolution.json').read_text())
        xml, log = (args.archive/'case.xml').read_text(), (args.archive/'run.log').read_text()
        flow = args.archive.parent.parent
        runtime_path, validation_path = flow/'runtime-before.json', args.archive.parent/'validation.json'
        runtime = json.loads(runtime_path.read_text()); validation = json.loads(validation_path.read_text())
        report['inputs'] = inputs
        report['producing_runtime'] = dict(source_commit=runtime['source_commit'],
            program=validation['program'], runtime_pin_path=str(runtime_path), runtime_pin_sha256=sha(runtime_path),
            validation_path=str(validation_path), validation_sha256=sha(validation_path),
            actual_command=next(c for c in validation['commands'] if c['name'] == args.archive.name),
            interpretation='Historical producing binary, distinct from the later capped Functions startup correction')
        report['history'] = {key: dict(entries=len(raw[key]), canonical_sha256=canonical(raw[key])) for key in
            ('initial', 'state', 'accepted_states', 'audit_times', 'native_history', 'native_history_times',
             'last_stage_history', 'fixed_reference', 'fixed_defect', 'raw_reference_loads', 'transition_reference')}
        report['archived_metadata'] = {key: raw[key] for key in ('time', 'steps', 'dt', 'order', 'ranks',
            'parameters', 'robust_parameters', 'transition_controls', 'quadrature_orders', 'offsets', 'failure')}
        report['archive_limits'] = dict(rejected_newton_vector_exported=bool(raw.get('rejected_state')),
            final_stage_loads_exported=bool(raw.get('final_stage_loads')), checkpoint_exported=bool(raw.get('checkpoint')),
            inference='No rejected trial vector, line-search trace or Jacobian age is archived; Newton cause is unresolved')
        proof_path = args.source/'docs/plans/cls-epawf-revised-evidence/m8/proof-local-joint-11aa.json'
        proof = json.loads(proof_path.read_text())
        report['certified_ranges'] = dict(path=str(proof_path), sha256=sha(proof_path),
            finite_boxes=[{key: box[key] for key in ('beta', 'gradient_entry_radius', 'k_over_k0', 'epsilon_over_epsilon0',
                'normalization', 'branch_scope', 'qualified_property')} | dict(gamma_candidates=[
                    {key: c[key] for key in ('gamma_over_nut0', 'joint_constitutive_gate', 'certified_positive',
                                            'certified_principal_lower_over_nut0')} for c in box['fixed_gamma_candidates']])
                for box in proof['finite_boxes']],
            interpretation='Outer constitutive membership only; the radius .01 joint box requires GammaScale3, not this run\'s 1. '
                'Tracefree principal-form certificates do not assert strong discrete incompressibility or inner AWF realizability. '
                'The disproven broad box is not a certified admissible envelope.')
        if len(raw['accepted_states']) != len(raw['audit_times']) or len(raw['state_audits']) != len(raw['accepted_states']):
            raise ValueError('Accepted-state audit histories have incompatible lengths')
        engine = Reconstruction(raw, args.source, xml)
        engine.near(raw['accepted_states'][0], raw['initial'], 'initial_archive_state', 0.)
        engine.near(raw['accepted_states'][-1], raw['state'], 'last_accepted_archive_state', 0.)
        engine.near(raw['time'], raw['audit_times'][-1], 'last_accepted_time', 0.)
        differential = raw['offsets'][raw['dimension']+2]
        for i, (state, time, native) in enumerate(zip(raw['accepted_states'], raw['audit_times'], raw['state_audits'])):
            if i == 0:
                history, alpha = list(state), 0.
            else:
                history = np.asarray(raw['accepted_states'][i-1], float).copy()
                alpha = 1/raw['dt'] if i == 1 else 1.5/raw['dt']
                if i > 1:
                    history[:differential] = (4*np.asarray(raw['accepted_states'][i-1])[:differential]-
                        np.asarray(raw['accepted_states'][i-2])[:differential])/3
            record = engine.state_record(state, history, alpha, time, native)
            record['accepted_step'] = i
            report['states'].append(record)
            for region in ('inner', 'outer'):
                candidates = [(key, value) for key, value in record['sets'].items()
                    if value['minimum_covariance'] is not None and key in (
                        ('source_volume_inner', 'physical_sample_K', 'interface_inner') if region == 'inner' else
                        ('source_volume_outer', 'physical_sample_outer', 'interface_outer'))]
                key, minimum = min(candidates, key=lambda item: item[1]['minimum_covariance'])
                if minimum['minimum_covariance'] < 0 and not any(event['type'] == 'first_negative_'+region+'_covariance' for event in report['events']):
                    report['events'].append(dict(type='first_negative_'+region+'_covariance', time=time,
                        accepted_step=i, set=key, detail=minimum['worst_point'], scope='actual archived constitutive finite sets'))
            for radius, membership in record['constitutive_range_membership'].items():
                event_type = 'first_outer_box_exit_radius_'+radius
                if membership['outside_count'] and not any(event['type'] == event_type for event in report['events']):
                    report['events'].append(dict(type=event_type, time=time, accepted_step=i,
                        detail=membership['first_outside'], interpretation='Certificate membership lost; not a failure proof'))
            for index, event_type in ((4, 'first_limiter_active_volume'), (5, 'first_cap_attenuation_volume')):
                if record['robust_source_independent'][index] > 1e-12 and not any(event['type'] == event_type for event in report['events']):
                    report['events'].append(dict(type=event_type, time=time, accepted_step=i,
                        volume=record['robust_source_independent'][index], minimum_lambda=record['robust_source_independent'][6]))
            report['checks'], report['errors'] = engine.checks, engine.errors
            save(output, report)
            print(f'accepted step {i:02d} time={time:.6g} independently reconstructed', flush=True)
        engine.near(raw['native_history'][0], raw['state'][:differential], 'native_history_current')
        engine.near(raw['native_history'][1], raw['accepted_states'][-2][:differential], 'native_history_previous')
        for slot in raw['native_history'][2:]:
            engine.near(slot, np.zeros(differential), 'native_zero_prescribed_body_history', 0.)
        engine.near(np.asarray(raw['last_stage_history'])[:differential], np.asarray(history)[:differential], 'independent_IMEX2_history')
        engine.near(raw['explicit_audits'], [[t, 0.] for t in raw['audit_times']], 'zero_explicit_body_at_every_epoch')
        report['events'].extend([dict(type='last_accepted_stage', time=raw['time'], accepted_step=raw['steps']),
            dict(type='first_failed_stage', time=raw['time']+raw['dt'], alpha=1.5/raw['dt'],
                 message=raw['failure'], classification='unqualified numerical stop; physical rejection cause not exported')])
        report['interpretation'] = dict(initial_inner_failure_is_preexisting=True,
            cap_attenuation_observed=any(s['robust_source_independent'][5] > 1e-12 for s in report['states']),
            outer_negative_covariance_observed=any(e['type'] == 'first_negative_outer_covariance' for e in report['events']),
            Newton_cause='Unresolved without failed iterate, line-search or controlled fresh-J comparison',
            physical_model_qualified=False, trajectory_status='Diagnostic accepted prefix only; already outside inner realizability at initialization')
        for name, receipt in inputs.items():
            engine.near(float(sha(args.archive/name) == receipt['sha256']), 1., 'archive_bytes_unchanged', 0.)
        report['source_clean_end'] = not bool(git(args.source, 'status', '--porcelain'))
        report['source_commit_end'] = git(args.source, 'rev-parse', 'HEAD')
        if not report['source_clean_end'] or report['source_commit_end'] != BASELINE:
            raise ValueError('Pinned reconstruction source changed during archive analysis')
        report['imported_source_files'] = {str(Path(module.__file__).resolve()): sha(module.__file__) for module in
            list(sys.modules.values()) if getattr(module, '__file__', None) and str(Path(module.__file__).resolve()).startswith(str(args.source))}
        report['checks'], report['errors'], report['success'] = engine.checks, engine.errors, True
        report['finished_utc'] = stamp(); save(output, report)
        lines = ['| Event | Time | Observation |', '|---|---:|---|']
        for event in sorted(report['events'], key=lambda value: value['time']):
            detail = event.get('detail', {})
            observation = ('eigenvalue='+str(detail['eigenvalue'])+' at '+str(detail['coordinate'])) if 'eigenvalue' in detail else event.get('interpretation', event.get('message', str(event.get('volume', ''))))
            lines.append(f"| {event['type']} | {event['time']:.6g} | {observation} |")
        (args.output/'archive-events.md').write_text('\n'.join(lines)+'\n')
        print(json.dumps(dict(success=True, checks=engine.checks, events=report['events'], output=str(output)), allow_nan=False))
        return 0
    except Exception as error:
        report['failures'].append(dict(type=type(error).__name__, message=str(error), traceback=traceback.format_exc()))
        report['finished_utc'] = stamp(); save(output, report)
        print(traceback.format_exc(), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
