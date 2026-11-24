#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent steady physical Q3/Q2 residual and budget oracle.

No native matrices, sampled bases, constitutive kernels, manufactured forcing,
time-mass terms, or stored native residuals enter the calculation. Global Riesz
operators are sparse. Element-local gradients retain their original owner.
An output is evidence, never an implicit full H27 qualification.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'h26'))
import numpy as np
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from scipy.optimize import brentq
import fixtures
from fixtures import modified

SCHEMA = 'WJH27IndependentPhysicalOracleV1'
TABLE_SHA = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'
EXACT_TABLE_SHA = '3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1'
EXACT_SCHEMA = 'WJH27ExactInnerFiniteNuTablesV1'
FULL_TRACTION = 'FullTractionV1'
ZERO_TRANSVERSE = 'ZeroTransverseV1'
SIGNS = dict(u='us', p='ps', k='ks', epsilon='epsilon_signs')
EDGES = ((0, 4, 8, 12), (1, 5, 9, 13), (0, 1, 2, 3), (4, 5, 6, 7))
NORMALS = ((-1., 0.), (1., 0.), (0., -1.), (0., 1.))
SIDES = ('inlet', 'outlet', 'bottom', 'top')
STATIONS = sorted(set([-10, -7.5, -5, -2.5, -1.25, 0, 1.25, 2.5, 4, 5.5,
                       6.5, 7.5, 10, 12.5, 15, 3.75, 5, 6.25] +
                      [3.75 + .25 * i for i in range(18)]))
RESIDUAL_LIMITS = dict(momentum=1e-9, continuity=1e-9, k=1e-8, epsilon=1e-8)
ENVELOPES = dict(U=(.01, .02), W=(.01, .02), k=(.05, .10), epsilon=(.10, .20))
NUMERICAL = dict(U=1e-4, W=1e-4, k=1e-3, epsilon=1e-2, R=1e-3,
                 Pi=1e-3, G=1e-2, traction=1e-3)


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def identity(path):
    path = Path(path)
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def outlet_identity(boundary_identity, policy):
    return hashlib.sha256(f'{boundary_identity}|{policy}'.encode()).hexdigest()


def coordinate_aliases(values, tolerance=1e-12):
    """Group geometrically coincident edge coordinates for topology only."""
    aliases = {}
    anchor = -float('inf')
    for coordinate in sorted(set(values)):
        if coordinate-anchor > tolerance:
            anchor = coordinate
        aliases[coordinate] = anchor
    return aliases


def basis(e, key, x, r):
    degree = 2 if key == 'p' else 3
    xl, xr, rl, rr = e['bounds']
    bx, dx = modified(degree, 2 * (x-xl)/(xr-xl)-1)
    br, dr = modified(degree, 2 * (r-rl)/(rr-rl)-1)
    signs = np.asarray(e[SIGNS[key]])
    n = (degree+1)**2
    return ((br[:, :, None]*bx[:, None, :]).reshape(len(x), n)*signs,
            (br[:, :, None]*dx[:, None, :]).reshape(len(x), n)*signs*2/(xr-xl),
            (dr[:, :, None]*bx[:, None, :]).reshape(len(x), n)*signs*2/(rr-rl))


def full_wj(k, epsilon, gradient):
    """Independent batched companion-eigenvalue root and full 3D WJ tensor.

    Uses no iterative native root or native derivative. Retains trace(G)/3,
    including S_yy for a weakly divergence-free planar finite element field.
    """
    k, epsilon = np.asarray(k), np.asarray(epsilon)
    require(np.all(k > 0) and np.all(epsilon > 0), 'nonpositive WJ input')
    I = np.eye(3)
    S = (k/epsilon)[:, None, None] * ((gradient+gradient.transpose(0, 2, 1))/2 -
          np.trace(gradient, axis1=1, axis2=2)[:, None, None]*I/3)
    O = (k/epsilon)[:, None, None] * (gradient-gradient.transpose(0, 2, 1))/2
    S2, O2 = S@S, O@O
    s = np.trace(S2, axis1=1, axis2=2)
    w = np.trace(O2, axis1=1, axis2=2)
    i4 = np.trace(S@O2, axis1=1, axis2=2)
    companion = np.zeros((len(k), 3, 3))
    companion[:, 1, 0] = companion[:, 2, 1] = 1.
    companion[:, 0, 2] = -3.6*w
    companion[:, 1, 2] = 2*w+2.7*s
    companion[:, 2, 2] = 1.8
    roots = np.linalg.eigvals(companion)
    selected = (np.abs(roots.imag) <= 1e-10*np.maximum(1., np.abs(roots.real))) & (roots.real >= 1.8*(1-1e-12))
    require(np.all(selected.sum(axis=1) == 1), 'independent WJ root absent or ambiguous')
    N = roots.real[np.arange(len(k)), np.argmax(selected, axis=1)]
    defect = (N-1.8)*(N*N-2*w)-2.7*s*N
    denominator = N**3+1.8*N*N+np.abs((2*w+2.7*s)*N)+np.abs(3.6*w)
    require(np.all(np.abs(defect) <= 5e-12*denominator), 'independent cubic root residual failed')
    conditioning = 3*N*N-3.6*N-2*w-2.7*s
    require(np.all(conditioning > 0), 'independent physical cubic root is not simple')
    q = (5/6)*(N*N-2*w)*(2*N*N-w)
    b1, b3 = -N*(2*N*N-7*w)/q, -12*i4/(N*q)
    b4, b6, b9 = -2*(N*N-2*w)/q, -6*N/q, 6/q
    a = (b1[:, None, None]*S + b3[:, None, None]*(O2-w[:, None, None]*I/3) +
         b4[:, None, None]*(S@O-O@S) + b6[:, None, None]*(S@O2+O2@S-2*i4[:, None, None]*I/3) +
         b9[:, None, None]*(O@S@O2-O2@S@O))
    R = k[:, None, None]*(2*I/3+a)
    require(np.all(np.isfinite(R)), 'nonfinite independently evaluated WJ tensor')
    return R, N, np.abs(defect)/denominator, conditioning


def load_reference(path):
    table = json.loads(path.read_text())
    reference = Reference(table)
    require(identity(path)['sha256'] == reference.boundary_identity, 'wrong immutable qualified reference table')
    return reference


class Reference:
    def __init__(self, table):
        self.exact_inner = table['schema'] == EXACT_SCHEMA
        require(table['schema'] in ('WJH27V2ReferenceTablesV1', EXACT_SCHEMA) and table['status'] == 'PASS' and
                table['model'] == 'IsotropicLowerWallFullWJCollarV2', 'unqualified v2 reference table')
        if self.exact_inner:
            require(table['policy'] == 'ExactInnerFiniteNuV1' and table['transport'] == 'ScalarCdV2',
                    'wrong exact-inner policy or transport')
            for component, digest in table['bindings'].items():
                encoded = json.dumps(table[component], sort_keys=True, separators=(',', ':'), allow_nan=False).encode()
                require(hashlib.sha256(encoded).hexdigest() == digest, 'wrong exact-inner binding: '+component)
        p, c, inlet = table['parameters'], table['constants'], table['inlet_binding']
        self.table = table
        self.identity = EXACT_TABLE_SHA if self.exact_inner else TABLE_SHA
        self.boundary_identity = self.identity
        self.inner_k_transport = 'ScalarCdV2'
        self.Kin, self.H, self.h, self.nu = p['K'], p['H'], p['h'], p['nu']
        self.cd, self.sk, self.se = p['cd'], p['sigma_k'], table['guard']['sigma_e']
        self.cl, self.qe, self.Ne, self.kappa = c['c_length'], c['q_e'], c['N_e'], c['kappa']
        self.d0, self.top, self.utau = inlet['d0'], inlet['top_model_radius'], inlet['u_tau']
        self.A, self.tau = self.Kin**1.5/self.cl, self.Kin/self.qe
        self.z0 = .005676270321216798
        if self.exact_inner:
            require(inlet['z0e'] == self.z0, 'outer calibration changed')
        self.lo, self.hi = table['cutoff_table']['domain']
        self.coefficients = table['cutoff_table']['coefficients']
        derivative = np.polynomial.chebyshev.chebder(self.coefficients)*2/(self.hi-self.lo)
        reserve = 32*(len(derivative)+1)*np.finfo(float).eps*float(np.sum(abs(derivative)))
        self.log_cutoff_slope_bounds = (float(derivative[0]-np.sum(abs(derivative[1:]))-reserve),
                                       float(derivative[0]+np.sum(abs(derivative[1:]))+reserve))
        require(self.log_cutoff_slope_bounds[1] < 0, 'independent cutoff inverse lacks a decreasing certificate')
        self.cutoff_radius_range = tuple(self.cutoff(np.asarray([self.hi,self.lo])))
        self.a = self.kappa*math.sqrt(self.Kin/self.qe)
        self.zv = float(self.intercept(np.array([self.Kin]))[0])
        require(self.zv < 0, 'historical inlet requires the fixed negative roughness intercept')
        self.D0 = self.nu-self.a*self.zv
        self.old_c = self.H*math.exp(1-self.old_integral_p(self.H))
        self.inlet_cutoff = table['inlet']['cutoff'] if self.exact_inner else self.old_c
        self.reference_shear = self.tau
        self.top_tangent = self.tau
        if self.exact_inner:
            self.reference_shear += self.nu*self.utau/(self.kappa*self.H)
            self.top_tangent += self.nu*self.utau/(self.kappa*self.top)
            self.inlet_panels = table['inlet']['panels']
            self.inlet_panel_bounds = np.array([p['domain'][1] for p in self.inlet_panels])
        self.top_flux = self.nu*self.A/self.top**2+self.cd*self.Kin**2/(self.se*self.top)
        self.inlet_quadrature = (leggauss(64), leggauss(96))
        require(0 < self.old_c < self.H, 'invalid historical inlet cutoff')

    def cutoff(self, K):
        K = np.asarray(K, dtype=float)
        require(np.all(np.isfinite(K)) and np.all((K >= self.lo) & (K <= self.hi)),
                f'current interface K range [{float(np.min(K))},{float(np.max(K))}] left qualified [{self.lo},{self.hi}]; no clamp')
        c = self.H*np.exp(chebval((2*K-self.lo-self.hi)/(self.hi-self.lo), self.coefficients))
        require(np.all((c > 0) & (c < .4*self.H)), 'conditional current cutoff branch failed')
        return c

    def intercept(self, K):
        rough = self.h*np.sqrt(K)/(70*self.nu)
        exponent = np.maximum(.5-.4*rough**.7, 1-.79*rough**(-.28))
        return -10.7*self.nu*np.expm1(exponent*np.log(rough))/np.sqrt(K)

    def aux(self, K, r):
        zv = self.intercept(K)
        require(np.all(zv < 0), 'dynamic roughness branch changed')
        return self.kappa*np.sqrt(K/self.qe)*(r-zv)

    def k_diffusivity(self, k, epsilon, auxiliary, outer):
        return self.nu+self.cd*k*k/(self.sk*epsilon)

    def old_integral_p(self, r):
        return np.log1p(self.a*r/self.D0)+self.nu/(self.D0+self.a*r)-self.nu/self.D0

    def old_Q(self, r):
        r = np.asarray(r, dtype=float)
        sink = np.empty_like(r)
        lower = r <= self.old_c
        sink[lower] = r[lower]/self.old_c
        sink[~lower] = 1+np.log(r[~lower]/self.old_c)
        return self.A*(sink-self.old_integral_p(r))

    def inlet(self, r):
        r = np.asarray(r, dtype=float)
        require(np.all((r >= -1e-13) & (r <= self.top+1e-12)), 'inlet physical radius outside domain')
        if self.exact_inner:
            return self.exact_inlet(r)
        U, k, G = np.empty_like(r), np.full_like(r, self.Kin), np.empty_like(r)
        inner = r <= self.H
        U[inner] = self.tau/self.a*np.log1p(self.a*r[inner]/self.D0)
        G[inner] = self.tau/(self.D0+self.a*r[inner])
        U[~inner] = self.utau/self.kappa*np.log(r[~inner]/self.z0)
        G[~inner] = self.utau/(self.kappa*r[~inner])
        # Direct split physical quadrature of the historical profile, never a
        # counterfactual-column substitution or native InletK call.
        ri = r[inner]
        values = []
        for a, w in self.inlet_quadrature:
            integral = np.zeros_like(ri)
            for lower, upper in ((ri, np.maximum(ri, self.old_c)),
                                 (np.maximum(ri, self.old_c), np.full_like(ri, self.H))):
                z = (lower[:, None]+upper[:, None])/2+(upper-lower)[:, None]*a/2
                integral += (upper-lower)/2*((self.old_Q(z)/(self.D0+self.a*z))@w)
            values.append(integral)
        require(np.max(np.abs(values[1]-values[0]), initial=0.) <= 2e-12*self.Kin,
                'historical inlet independent split quadrature unresolved')
        k[inner] -= values[1]
        return dict(U=U, W=np.zeros_like(r), k=k, G=G)

    def exact_inlet(self, r):
        inner = r <= self.H
        ri = r[inner]
        indices = np.searchsorted(self.inlet_panel_bounds, ri)
        U, G = np.empty_like(r), np.empty_like(r)
        k, kr, Pi = np.full_like(r, self.Kin), np.zeros_like(r), np.zeros_like(r)
        U[~inner] = self.utau/self.kappa*np.log(r[~inner]/self.z0)
        G[~inner] = self.utau/(self.kappa*r[~inner])
        inner_ids = np.flatnonzero(inner)
        for index in np.unique(indices):
            ids = inner_ids[indices == index]
            panel = self.inlet_panels[index]
            lo, hi = panel['domain']
            x = (2*r[ids]-lo-hi)/(hi-lo)
            for target, key in ((U, 'U'), (k, 'k'), (Pi, 'Pi')):
                target[ids] = chebval(x, panel[key])
            for target, key in ((G, 'U'), (kr, 'k')):
                target[ids] = chebval(x, np.polynomial.chebyshev.chebder(panel[key]))*2/(hi-lo)
        return dict(U=U, W=np.zeros_like(r), k=k, G=G, kr=kr, Pi=Pi)

    def outlet(self, r):
        r = np.asarray(r, dtype=float)
        if self.exact_inner:
            tz = np.full_like(r, self.reference_shear)
            outer = r > self.H
            tz[outer] = self.tau+self.nu*self.utau/(self.kappa*r[outer])
            return np.column_stack((-self.Kin/self.Ne-2*self.exact_inlet(r)['Pi'], tz))
        tx = np.full_like(r, self.Kin/self.Ne)
        for j, segment in enumerate(self.table['outlet_traction']['segments']):
            mask = ((r > .4*self.H) & (r <= .75*self.H) if j == 0 else
                    (r > .75*self.H) & (r <= self.H) if j == 1 else r > self.H)
            arg = np.log(r[mask]/self.H) if j == 2 else r[mask]
            lo, hi = segment['domain']
            require(np.all((arg >= lo-1e-13) & (arg <= hi+1e-13)), 'outlet reference domain error')
            tx[mask] = chebval((2*arg-lo-hi)/(hi-lo), segment['tx_coefficients'])
        return np.column_stack((tx, np.full_like(r, self.tau)))


class SparseMass:
    def __init__(self, size):
        self.size, self.rows, self.cols, self.values = size, [], [], []

    def add(self, ids, matrix):
        self.rows.append(np.repeat(ids, len(ids)))
        self.cols.append(np.tile(ids, len(ids)))
        self.values.append(matrix.ravel())

    def finish(self):
        result = coo_matrix((np.concatenate(self.values),
                             (np.concatenate(self.rows), np.concatenate(self.cols))),
                            shape=(self.size, self.size)).tocsc()
        result.sum_duplicates()
        self.rows = self.cols = self.values = []
        return result


class Norms:
    def __init__(self):
        self.data = {}

    def add(self, region, name, difference, weights):
        d = np.asarray(difference)
        if d.ndim > 1:
            # Component-wise RMS/MAX, then the most demanding component.
            d = d.reshape(len(weights), -1)
        else:
            d = d[:, None]
        key = region + '/' + name
        a = self.data.setdefault(key, dict(measure=0., square=np.zeros(d.shape[1]), maximum=np.zeros(d.shape[1])))
        a['measure'] += float(np.sum(weights))
        a['square'] += weights @ (d*d)
        a['maximum'] = np.maximum(a['maximum'], np.max(np.abs(d), axis=0))

    def report(self, limits):
        out = {}
        for key, a in self.data.items():
            name = key.split('/')[-1]
            rms = np.sqrt(a['square']/a['measure']) if a['measure'] > 0 else np.zeros_like(a['square'])
            limit = limits.get(name)
            if isinstance(limit, (int, float)):
                limit = (limit, limit)
            out[key] = dict(RMS=float(np.max(rms)), MAX=float(np.max(a['maximum'])),
                            component_RMS=rms.tolist(), component_MAX=a['maximum'].tolist(), measure=a['measure'])
            if limit is not None:
                out[key].update(limits=list(limit), status='MET' if max(rms) <= limit[0] and
                                max(a['maximum']) <= limit[1] and a['measure'] > 0 else 'NOT_MET')
        return out


class PhysicalOracle:
    def __init__(self, data, reference, order, deadline, *, outlet_x=15.):
        self.data, self.ref, self.order, self.deadline = data, reference, order, deadline
        self.n, self.np, self.nk, self.ne = map(int, data['sizes'])
        self.elements, self.geometry = data['elements'], data['geometry']
        require(self.geometry.get('mapping', 'Flat') == 'Flat',
                'use the normal-coordinate audit for a mapped WJ state')
        require(data['schema'] == 'WJH27NativeCoupledStateV2' and data['manufactured'] is False and
                data['gauge'] is False and data['model'] == 'FiniteNuIsotropicWallFullWJCollarV2',
                'physical oracle requires actual open-outlet v2 state, no MMS or pressure gauge')
        require(data['reference_identity'] == reference.identity and
                data.get('boundary_identity', TABLE_SHA) == reference.boundary_identity and
                data.get('inner_k_transport', 'ScalarCdV2') == reference.inner_k_transport,
                'native reference, boundary or inner transport identity mismatch')
        self.outlet_policy = data.get('outlet_policy', FULL_TRACTION)
        require(self.outlet_policy in (FULL_TRACTION, ZERO_TRANSVERSE), 'unknown physical outlet policy')
        require(self.outlet_policy == FULL_TRACTION or reference.exact_inner,
                'zero-transverse outlet requires the qualified exact-inner reference')
        require(outlet_x == 15. or (self.outlet_policy == ZERO_TRANSVERSE
                                    and math.isfinite(outlet_x) and outlet_x > 15.),
                'a downstream diagnostic requires the explicitly bound B domain')
        self.left, self.right = -10., outlet_x
        self.length = self.right-self.left
        self.outlet_identity = outlet_identity(reference.boundary_identity, self.outlet_policy)
        if 'outlet_policy' in data or 'outlet_identity' in data:
            require(data.get('outlet_policy') == self.outlet_policy and
                    data.get('outlet_identity') == self.outlet_identity,
                    'native outlet policy identity mismatch')
        self.s = {key: np.asarray(data['coefficients'][key], dtype=float) for key in ('flow', 'k', 'epsilon')}
        for key, length in zip(self.s, (2*self.n+self.np, self.nk, self.ne)):
            require(self.s[key].shape == (length,) and np.all(np.isfinite(self.s[key])), 'invalid '+key+' coefficients')
        digest = hashlib.sha256()
        for key in self.s:
            digest.update(np.asarray(self.s[key], dtype='<f8').tobytes())
        require(digest.hexdigest() == data['state_sha256'], 'native coefficient state hash mismatch')
        expected = dict(d0=reference.d0, H=reference.H, h=reference.h,
                        left=self.left, right=self.right, top=reference.top)
        for key, value in expected.items():
            require(abs(self.geometry[key]-value) <= 1e-12, 'retained physical geometry changed: '+key)
        self.volume, self.evolume = self.length*reference.top, self.length*(reference.top-reference.H)
        self.ufree = [np.asarray(v, dtype=int) for v in data['free_velocity_ids']]
        self.free = dict(k=np.asarray(data['k_free_ids'], dtype=int), epsilon=np.asarray(data['epsilon_free_ids'], dtype=int))
        self.owner, self.faces, self.internal = {}, dict(u=[], k=[], epsilon=[]), []
        self.validate_topology()
        self.gauss = leggauss(order)
        self.masses = None

    def clock(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('physical oracle exhausted explicit cooperative wall cap')

    def validate_topology(self):
        p = self.ref
        require(self.elements and len({e['geometry'] for e in self.elements}) == len(self.elements), 'empty/duplicate geometry inventory')
        visited = {key: set() for key in ('u', 'p', 'k', 'epsilon')}
        columns, edges = {}, {}
        for e in self.elements:
            xl, xr, rl, rr = e['bounds']
            require(xr > xl and rr > rl and xl >= self.left-1e-12 and xr <= self.right+1e-12 and
                    rl >= -1e-12 and rr <= p.top+1e-12 and not (rl < p.H-1e-12 and rr > p.H+1e-12),
                    'invalid physical rectangle or missing H mesh edge')
            # Snap trusted geometry endpoints only, never a field/state value.
            for j, value in enumerate(e['bounds']):
                for edge in ((self.left, self.right) if j < 2 else (0., p.H, p.top)):
                    if abs(value-edge) <= 64*np.finfo(float).eps*max(p.H, abs(value), abs(edge)):
                        e['bounds'][j] = edge
        x_alias = coordinate_aliases(e['bounds'][i] for e in self.elements for i in (0,1))
        r_alias = coordinate_aliases(e['bounds'][i] for e in self.elements for i in (2,3))
        for e in self.elements:
            xl, xr, rl, rr = e['bounds']
            outer = rl >= p.H
            require(bool(e['epsilon']) == outer, 'transported epsilon ownership mismatch')
            ax,bx,ar,br = x_alias[xl],x_alias[xr],r_alias[rl],r_alias[rr]
            require(ax<bx and ar<br, 'topology tolerance collapses an element')
            columns.setdefault((ax,bx), []).append(e)
            for key, size in (('u', self.n), ('p', self.np), ('k', self.nk), ('epsilon', self.ne)):
                ids, signs = e[key], e[SIGNS[key]]
                if key == 'epsilon' and not outer:
                    require(not signs, 'inner epsilon signs exist'); continue
                count = 9 if key == 'p' else 16
                require(len(ids) == len(signs) == count and len(set(ids)) == count and
                        all(0 <= i < size for i in ids) and all(s in (-1., 1.) for s in signs), 'invalid '+key+' signed map')
                visited[key].update(ids)
            for side, coordinates in enumerate(((ax, ax, ar, br), (bx, bx, ar, br),
                                                 (ax, bx, ar, ar), (ax, bx, br, br))):
                edges.setdefault(coordinates, []).append((e, side))
            for key in self.faces:
                if key == 'epsilon' and not outer:
                    continue
                bottom = p.H if key == 'epsilon' else 0.
                for side, present in enumerate((xl == self.left, xr == self.right, rl == bottom, rr == p.top)):
                    if present:
                        self.faces[key].append((e, side))
        ordered = sorted(columns)
        require(abs(ordered[0][0]-self.left) < 1e-12 and abs(ordered[-1][1]-self.right) < 1e-12 and
                all(abs(a[1]-b[0]) < 1e-12 for a, b in zip(ordered[:-1], ordered[1:])), 'x columns overlap or leave gaps')
        for column in columns.values():
            column.sort(key=lambda e: e['bounds'][2])
            require(column[0]['bounds'][2] == 0 and column[-1]['bounds'][3] == p.top and
                    all(abs(a['bounds'][3]-b['bounds'][2]) < 1e-12 for a, b in zip(column[:-1], column[1:])),
                    'radial column overlaps or leaves gaps')
            owners = [e for e in column if e['bounds'][2] == p.H]
            require(len(owners) == 1, 'missing or ambiguous outer H owner')
            for e in column:
                self.owner[e['geometry']] = owners[0]
        for key, size in (('u', self.n), ('p', self.np), ('k', self.nk), ('epsilon', self.ne)):
            require(visited[key] == set(range(size)), 'incomplete '+key+' global coefficient inventory')
        for owners in edges.values():
            require(len(owners) in (1, 2), 'nonmanifold physical edge')
            if len(owners) == 2:
                self.internal.append(owners)
            else:
                require(owners[0] in self.faces['u'], 'unpaired interior geometry edge')
        transverse_sides = (0, 1, 2, 3) if self.outlet_policy == ZERO_TRANSVERSE else (0, 2, 3)
        self.fixed = {}
        for key, size, free, sides in (('u', self.n, self.ufree[0], (0, 2)),
                ('w', self.n, self.ufree[1], transverse_sides), ('k', self.nk, self.free['k'], (0,)),
                ('epsilon', self.ne, self.free['epsilon'], (0, 2))):
            field = 'u' if key == 'w' else key
            fixed = {e[field][j] for e, side in self.faces[field] if side in sides for j in EDGES[side]}
            require(len(set(free)) == len(free) and set(free) == set(range(size))-fixed, 'wrong physical '+key+' essential mask')
            self.fixed[key] = np.asarray(sorted(fixed), dtype=int)
        outlets = [e for e, side in self.faces['u'] if side == 1]
        bottom = [e for e in outlets if e['bounds'][2] == 0.]
        top = [e for e in outlets if e['bounds'][3] == p.top]
        require(len(bottom) == len(top) == 1, 'wrong shared outlet corner ownership')
        self.outlet_corners = dict(ground=bottom[0]['u'][EDGES[1][0]], top=top[0]['u'][EDGES[1][1]])
        trace_ids = {e['u'][j] for e in outlets for j in EDGES[1]}
        require(trace_ids & set(self.fixed['u']) == {self.outlet_corners['ground']},
                'wrong axial outlet freedom')
        fixed_transverse = trace_ids if self.outlet_policy == ZERO_TRANSVERSE else set(self.outlet_corners.values())
        require(trace_ids & set(self.fixed['w']) == fixed_transverse, 'wrong transverse outlet freedom')
        self.outlet_ids = np.asarray(sorted(trace_ids), dtype=int)
        self.outlet_tangent_fixed = np.asarray(sorted(fixed_transverse), dtype=int)

    def outlet_traction(self, r):
        traction = self.ref.outlet(r)
        if self.outlet_policy == ZERO_TRANSVERSE:
            traction[:, 1] = 0.
        return traction

    def outlet_tangent_reactions(self, full_flow_residual):
        """Signed coefficient actions of the constrained transverse outlet tests."""
        corners = set(self.outlet_corners.values())
        return [dict(dof=int(dof), value=float(full_flow_residual[self.n+dof]),
                     scope='shared_corner' if dof in corners else 'outlet')
                for dof in self.outlet_tangent_fixed]

    def trace(self, e, key, x, derivative=0):
        owner = self.owner[e['geometry']]
        b = basis(owner, key, x, np.full(len(x), self.ref.H))[derivative]
        values = (self.s['flow'][2*self.n+np.asarray(owner['p'])] if key == 'p' else
                  self.s['flow'][owner['u']] if key == 'u' else self.s[key][owner[key]])
        return b@values

    def cuts(self, lo, hi, extras=()):
        eps = 64*np.finfo(float).eps*max(self.ref.H, abs(lo), abs(hi))
        values = [lo, hi, .4*self.ref.H, .75*self.ref.H, self.ref.h, self.ref.H, *extras]
        return sorted(set([lo, hi]+[float(v) for v in values if lo+eps < v < hi-eps]))

    def cutoff_x_splits(self, e, bounds=None):
        """Independent direct-c(K) roots; no native inverse levels or panels."""
        xl,xr,rl,rr = bounds if bounds is not None else e['bounds']
        if rr < self.ref.cutoff_radius_range[0] or rl > self.ref.cutoff_radius_range[1]: return []
        fit_nodes = np.asarray([-1.,-.5,.5,1.]);points = (xl+xr)/2+(xr-xl)*fit_nodes/2
        values = self.trace(e,'k',points)
        polynomial = np.polynomial.Chebyshev.fit(points,values,3,domain=[xl,xr])
        held = (xl+xr)/2+(xr-xl)*np.asarray([-.81,-.27,.19,.71])/2
        require(float(np.max(abs(polynomial(held)-self.trace(e,'k',held)))) <= 1e-12*self.ref.Kin,
                'independent interface trace is not the declared cubic')
        critical = []
        for value in polynomial.deriv().roots():
            if abs(np.imag(value)) <= 64*np.finfo(float).eps*max(abs(xl),abs(xr),xr-xl):
                point = float(np.real(value))
                if xl < point < xr: critical.append(point)
        panels = sorted(set([xl,xr,*critical]));cut = self.ref.cutoff(self.trace(e,'k',np.asarray(panels)))
        if not any(float(np.min(cut)) <= face <= float(np.max(cut)) for face in (rl,rr)): return []
        roots = list(critical)
        for face in (rl,rr):
            def residual(point):return float(self.ref.cutoff(self.trace(e,'k',np.asarray([point])))[0]-face)
            for left,right in zip(panels,panels[1:]):
                a,b = residual(left),residual(right)
                if (a < 0 < b) or (b < 0 < a):
                    roots.append(brentq(residual,left,right,xtol=8*np.finfo(float).eps*max(abs(left),abs(right),right-left),
                                        rtol=4*np.finfo(float).eps))
        return sorted(set(roots))

    def rule(self, e, bounds=None, extra_cutoff=None, extra_x_splits=()):
        a, w = self.gauss
        radial, radial_weights = getattr(self, 'radial_gauss', self.gauss)
        xl, xr, rl, rr = bounds if bounds is not None else e['bounds']
        xs, rs, weights = [], [], []
        panels = sorted(set([xl,xr,*self.cutoff_x_splits(e,bounds),*extra_x_splits]))
        require(panels[0] == xl and panels[-1] == xr and all(right > left for left,right in zip(panels,panels[1:])),
                'independent cutoff x partitions left the original owner')
        for left,right in zip(panels,panels[1:]):
            xnodes = (left+right)/2+(right-left)*a/2
            cutoff = self.ref.cutoff(self.trace(e, 'k', xnodes))
            extra = extra_cutoff(xnodes) if extra_cutoff else cutoff
            for x, wx, c, c2 in zip(xnodes, w, cutoff, extra):
                cuts = self.cuts(rl, rr, (c, c2))
                for lo, hi in zip(cuts[:-1], cuts[1:]):
                    xs.extend(np.full(len(radial), x)); rs.extend((lo+hi)/2+(hi-lo)*radial/2)
                    weights.extend(wx*radial_weights*(right-left)*(hi-lo)/4)
        return np.asarray(xs), np.asarray(rs), np.asarray(weights)

    def face(self, e, side, endpoints=False):
        a, w = self.gauss
        xl, xr, rl, rr = e['bounds']
        if side < 2:
            x = xl if side == 0 else xr
            cutoff = self.ref.cutoff(self.trace(e, 'k', np.array([x])))[0]
            cuts = self.cuts(rl, rr, (cutoff, self.ref.inlet_cutoff))
            rs = np.concatenate([(lo+hi)/2+(hi-lo)*a/2 for lo, hi in zip(cuts[:-1], cuts[1:])])
            weights = np.concatenate([w*(hi-lo)/2 for lo, hi in zip(cuts[:-1], cuts[1:])])
            if endpoints:
                rs = np.r_[rs, cuts]; weights = np.r_[weights, np.zeros(len(cuts))]
            return np.full(len(rs), x), rs, weights
        x, weights = (xl+xr)/2+(xr-xl)*a/2, w*(xr-xl)/2
        if endpoints:
            x, weights = np.r_[x, xl, xr], np.r_[weights, 0., 0.]
        return x, np.full(len(x), rl if side == 2 else rr), weights

    def sample(self, e, x, r):
        self.clock()
        ref = self.ref
        b, dx, dr = basis(e, 'u', x, r); ids = np.asarray(e['u'])
        u, w = b@self.s['flow'][ids], b@self.s['flow'][self.n+ids]
        G = np.zeros((len(x), 3, 3))
        G[:, 0, 0], G[:, 0, 2] = dx@self.s['flow'][ids], dr@self.s['flow'][ids]
        G[:, 2, 0], G[:, 2, 2] = dx@self.s['flow'][self.n+ids], dr@self.s['flow'][self.n+ids]
        kb, kx, kr = basis(e, 'k', x, r)
        k, gradk = kb@self.s['k'][e['k']], np.column_stack((kx@self.s['k'][e['k']], kr@self.s['k'][e['k']]))
        K = self.trace(e, 'k', x); cutoff = ref.cutoff(K)
        outer = bool(e['epsilon'])
        if outer:
            eb, ex, er = basis(e, 'epsilon', x, r)
            E = eb@self.s['epsilon'][e['epsilon']]
            grade = np.column_stack((ex@self.s['epsilon'][e['epsilon']], er@self.s['epsilon'][e['epsilon']]))
        else:
            E = K**1.5/ref.cl/np.maximum(r, cutoff)
            grade = np.zeros((len(x), 2))
        bad = np.flatnonzero((k <= 0) | (E <= 0) | ~np.isfinite(k) | ~np.isfinite(E))
        if len(bad):
            q = bad[0]
            raise ValueError(f'sampled scalar positivity/nonfinite failure: owner={e["geometry"]}, x={x[q]}, r={r[q]}, k={k[q]}, E={E[q]}')
        s = (r/ref.H-.4)/.35
        chi = np.where(r <= .4*ref.H, 0., np.where(r >= .75*ref.H, 1., s**3*(10-15*s+6*s*s)))
        strain = G+G.transpose(0, 2, 1)
        R = 2*k[:, None, None]*np.eye(3)/3
        aux = np.zeros(len(x)) if outer else ref.aux(K, r)
        R[:, 0, 2] = R[:, 2, 0] = -aux*strain[:, 0, 2]
        R[:, 1, 2] = R[:, 2, 1] = -aux*strain[:, 1, 2]
        root, root_error = np.full(len(x), 1.8), np.zeros(len(x))
        active = chi > 0
        if np.any(active):
            full, root[active], root_error[active], _ = full_wj(k[active], E[active], G[active])
            R[active] = (1-chi[active, None, None])*R[active]+chi[active, None, None]*full
        V = ref.nu*strain-R+2*k[:, None, None]*np.eye(3)/3
        Pi = basis(e, 'p', x, r)[0]@self.s['flow'][2*self.n+np.asarray(e['p'])]
        P = -np.einsum('qij,qij->q', R, G)
        Dk = ref.k_diffusivity(k,E,aux,outer)
        De = ref.nu+ref.cd*k*k/(ref.se*E) if outer else np.zeros(len(x))
        div = G[:, 0, 0]+G[:, 2, 2]
        return dict(U=u, W=w, G=G, k=k, E=E, K=K, cutoff=cutoff, chi=chi, R=R, V=V, Pi=Pi,
                    P=P, Dk=Dk, De=De, gradk=gradk, grade=grade, div=div, root=root,
                    root_error=root_error, covariance=np.linalg.eigvalsh(R)[:, 0]/k, outer=outer)

    def zeros(self):
        return dict(flow=np.zeros(2*self.n+self.np), k=np.zeros(self.nk), epsilon=np.zeros(self.ne))

    def build_riesz(self):
        self.clock()
        masses = {key: SparseMass(size) for key, size in (('u', self.n), ('p', self.np), ('k', self.nk), ('epsilon', self.ne))}
        self.one = {key: np.zeros(mass.size) for key, mass in masses.items()}
        seen = {key: {} for key in masses}
        for e in self.elements:
            self.clock(); x, r, weight = self.rule(e)
            for key, mass in masses.items():
                if not e[key]: continue
                b = basis(e, key, x, r)[0]
                scale = self.ref.Kin if key == 'k' else self.ref.A/r if key == 'epsilon' else 1.
                ids = np.asarray(e[key]); mass.add(ids, b.T@((weight*scale**2)[:, None]*b))
                width = 3 if key == 'p' else 4
                for j, gid in enumerate(ids):
                    value = 1/e[SIGNS[key]][j] if j % width < 2 and j//width < 2 else 0.
                    require(gid not in seen[key] or seen[key][gid] == value, 'inconsistent constant signed-map representation')
                    seen[key][gid] = value; self.one[key][gid] = value
        self.masses = {key: mass.finish() for key, mass in masses.items()}
        self.factors = {key: splu(self.masses[key][ids][:, ids].tocsc()) for key, ids in self.free.items()}
        self.factors['u'] = [splu(self.masses['u'][ids][:, ids].tocsc()) for ids in self.ufree]
        self.factors['p'] = splu(self.masses['p'])

    def dual(self, residual):
        require(self.masses is not None, 'Riesz operators have not been built')
        energy = 0.
        for component, ids in enumerate(self.ufree):
            f = residual['flow'][component*self.n+ids]
            energy += float(f@self.factors['u'][component].solve(f))
        require(energy >= -1e-25, 'negative momentum sparse Riesz energy')
        result = dict(momentum=math.sqrt(max(0., energy)/self.volume))
        f = residual['flow'][2*self.n:]
        result['continuity'] = math.sqrt(max(0., float(f@self.factors['p'].solve(f)))/self.volume)
        for key in self.free:
            f = residual[key][self.free[key]]
            energy = float(f@self.factors[key].solve(f))
            require(energy >= -1e-25, 'negative scalar sparse Riesz energy')
            result[key] = math.sqrt(max(0., energy)/(self.volume if key == 'k' else self.evolume))
        return result

    def assemble(self):
        if self.masses is None:
            self.build_riesz()
        names = ('diffusion', 'pressure', 'advection', 'awf_convection', 'awf_pressure',
                 'awf_known_k', 'source', 'natural', 'continuity')
        parts = {name: self.zeros() for name in names}
        advective_scalar = self.zeros()
        quantities = ('k_inventory', 'epsilon_inventory', 'P', 'E', 'k_source', 'epsilon_source',
                      'kinetic_advection', 'stress_work', 'molecular_work', 'pressure_divergence',
                      'isotropic_k_divergence', 'kinetic_energy_divergence', 'roughness_work', 'divergence', 'strong_divergence_square')
        bulk = {region: {key: 0. for key in quantities} for region in ('whole', 'inner', 'outer')}
        norms, shear_diagnostic = Norms(), Norms()
        minima = dict(k=float('inf'), epsilon=float('inf'), covariance_over_k=float('inf'))
        max_root_error = 0.
        for e in self.elements:
            self.clock(); x, r, weight = self.rule(e)
            a = self.sample(e, x, r)
            b, dx, dr = basis(e, 'u', x, r); ui = np.asarray(e['u'])
            pb = basis(e, 'p', x, r)[0]; pi = np.asarray(e['p'])
            div, G = a['div'], a['G']
            hmask = r < self.ref.h
            Uh, Uxh = self.trace(e, 'u', x), self.trace(e, 'u', x, 1)
            force_conv = np.where(hmask, 2*Uh*Uxh, 0.)
            force_pi = np.where(hmask, self.trace(e, 'p', x, 1), 0.)
            force_k = np.where(hmask, -2*self.trace(e, 'k', x, 1)/3, 0.)
            force = force_conv+force_pi+force_k
            for component, index in ((0, 0), (1, 2)):
                ids = component*self.n+ui
                convection = a['U']*G[:, index, 0]+a['W']*G[:, index, 2]
                parts['advection']['flow'][ids] += b.T@(weight*convection)
                parts['diffusion']['flow'][ids] += dx.T@(weight*a['V'][:, index, 0])+dr.T@(weight*a['V'][:, index, 2])
                parts['pressure']['flow'][ids] -= (dx if component == 0 else dr).T@(weight*a['Pi'])
            for key, f in (('awf_convection', force_conv), ('awf_pressure', force_pi), ('awf_known_k', force_k)):
                parts[key]['flow'][ui] -= b.T@(weight*f)
            parts['continuity']['flow'][2*self.n+pi] += pb.T@(weight*div)
            ksource = a['P']-a['E']
            esource = (1.44*a['P']-1.82*a['E'])*a['E']/a['k'] if a['outer'] else np.zeros(len(x))
            for key, values, gradient, diffusivity, source in (
                    ('k', a['k'], a['gradk'], a['Dk'], ksource),
                    ('epsilon', a['E'], a['grade'], a['De'], esource)):
                if not e[key]: continue
                qb, qx, qr = basis(e, key, x, r); ids = np.asarray(e[key])
                parts['advection'][key][ids] -= qx.T@(weight*a['U']*values)+qr.T@(weight*a['W']*values)
                parts['diffusion'][key][ids] += qx.T@(weight*diffusivity*gradient[:, 0])+qr.T@(weight*diffusivity*gradient[:, 1])
                parts['source'][key][ids] -= qb.T@(weight*source)
                advective_scalar[key][ids] += qb.T@(weight*(a['U']*gradient[:, 0]+a['W']*gradient[:, 1]))
            conv_work = a['U']*(a['U']*G[:, 0, 0]+a['W']*G[:, 0, 2])+a['W']*(a['U']*G[:, 2, 0]+a['W']*G[:, 2, 2])
            values = dict(k_inventory=a['k'], epsilon_inventory=a['E'] if a['outer'] else np.zeros(len(x)),
                P=a['P'], E=a['E'], k_source=ksource, epsilon_source=esource,
                kinetic_advection=conv_work, stress_work=np.einsum('qij,qij->q', a['V'], G),
                molecular_work=self.ref.nu*np.einsum('qij,qij->q', G+G.transpose(0, 2, 1), G),
                pressure_divergence=a['Pi']*div, isotropic_k_divergence=2*a['k']*div/3,
                kinetic_energy_divergence=.5*(a['U']**2+a['W']**2)*div,
                roughness_work=a['U']*force, divergence=div, strong_divergence_square=div*div)
            region = 'outer' if a['outer'] else 'inner'
            for target in ('whole', region):
                for key, value in values.items():
                    bulk[target][key] += float(weight@value)
            self.inlet_norms(norms, 'domain/'+region, a, r, weight)
            self.inlet_norms(norms, 'domain/whole', a, r, weight)
            minima['k'] = min(minima['k'], float(np.min(a['k'])))
            minima['epsilon'] = min(minima['epsilon'], float(np.min(a['E'])))
            minima['covariance_over_k'] = min(minima['covariance_over_k'], float(np.min(a['covariance'])))
            max_root_error = max(max_root_error, float(np.max(a['root_error'])))

        face_rows, boundary = [], dict(mass=0., kinetic_advection=0., traction_actual_work=0.,
            traction_prescribed_work=0., traction_essential_actual_work=0.,
            k_advection=0., k_diffusion_actual=0., k_diffusion_essential=0., k_diffusion_prescribed=0.,
            epsilon_advection=0., epsilon_diffusion_actual=0., epsilon_diffusion_essential=0., epsilon_diffusion_prescribed=0.)
        for key in ('u', 'k', 'epsilon'):
            for e, side in self.faces[key]:
                self.clock(); x, r, weight = self.face(e, side, endpoints=True)
                a = self.sample(e, x, r); nx, nr = NORMALS[side]
                un = nx*a['U']+nr*a['W']
                row = dict(field=key, geometry=e['geometry'], side=SIDES[side], bounds=e['bounds'])
                if key == 'u':
                    ui = np.asarray(e['u']); b = basis(e, 'u', x, r)[0]
                    traction = np.column_stack((nx*(a['V'][:, 0, 0]-a['Pi'])+nr*a['V'][:, 0, 2],
                                                nx*a['V'][:, 2, 0]+nr*(a['V'][:, 2, 2]-a['Pi'])))
                    prescribed = np.zeros((len(x), 2))
                    if side == 1:
                        prescribed = self.outlet_traction(r)
                    elif side == 3:
                        prescribed[:, 0] = self.ref.top_tangent
                    for component in (0, 1):
                        parts['natural']['flow'][component*self.n+ui] -= b.T@(weight*prescribed[:, component])
                    natural_work = a['U']*prescribed[:, 0]+a['W']*prescribed[:, 1]
                    actual_work = a['U']*traction[:, 0]+a['W']*traction[:, 1]
                    essential_work = np.zeros(len(x))
                    if side in (0, 2):
                        essential_work = actual_work
                    elif side == 3 or (side == 1 and self.outlet_policy == ZERO_TRANSVERSE):
                        essential_work = a['W']*traction[:, 1]
                    metrics = dict(mass=float(weight@un), kinetic_advection=float(weight@(.5*(a['U']**2+a['W']**2)*un)),
                        traction_actual_work=float(weight@actual_work), traction_prescribed_work=float(weight@natural_work),
                        traction_essential_actual_work=float(weight@essential_work))
                    for name, value in metrics.items(): boundary[name] += value
                    row.update(metrics, minimum_normal_velocity=float(np.min(un)),
                        actual_traction_integral=(weight@traction).tolist(), prescribed_traction_integral=(weight@prescribed).tolist())
                    if side in (1, 3):
                        components = 2 if side == 1 and self.outlet_policy == FULL_TRACTION else 1
                        compared = (traction[:, :components]-prescribed[:, :components])/self.ref.tau
                        row['natural_traction_scaled_RMS'] = float(np.sqrt(np.max(weight@(compared*compared)/weight.sum())))
                        row['natural_traction_scaled_MAX'] = float(np.max(np.abs(compared)))
                    if side == 1:
                        row.update(prescribed_normal_work=float(weight@(a['U']*prescribed[:, 0])),
                                   prescribed_tangential_work=float(weight@(a['W']*prescribed[:, 1])))
                    if side == 2:
                        norms.add('ground', 'shear', (a['V'][:, 0, 2]-self.ref.tau)/self.ref.tau, weight)
                        if self.ref.exact_inner:
                            shear_diagnostic.add('ground', 'shear_vs_reference',
                                      (a['V'][:, 0, 2]-self.ref.reference_shear)/self.ref.tau, weight)
                else:
                    b = basis(e, key, x, r)[0]; ids = np.asarray(e[key])
                    q, gradient, D = (a['k'], a['gradk'], a['Dk']) if key == 'k' else (a['E'], a['grade'], a['De'])
                    adv = un*q
                    actual = -D*(nx*gradient[:, 0]+nr*gradient[:, 1])
                    essential = side == 0 or (key == 'epsilon' and side == 2)
                    prescribed = np.full(len(x), self.ref.top_flux if key == 'epsilon' and side == 3 else 0.)
                    parts['advection'][key][ids] += b.T@(weight*adv)
                    parts['natural'][key][ids] += b.T@(weight*prescribed)
                    metrics = {key+'_advection': float(weight@adv), key+'_diffusion_actual': float(weight@actual),
                               key+'_diffusion_essential': float(weight@actual) if essential else 0.,
                               key+'_diffusion_prescribed': float(weight@prescribed)}
                    for name, value in metrics.items(): boundary[name] += value
                    row.update(metrics, essential=essential,
                        owner_gradient_diffusive_flux_MIN=float(np.min(actual)), owner_gradient_diffusive_flux_MAX=float(np.max(actual)))
                face_rows.append(row)
        total = {key: sum((part[key] for part in parts.values()), np.zeros_like(self.s[key])) for key in self.s}
        # Retained for the separately executable analytic mutation fixture;
        # these are independent arrays, never copied from native operators.
        self.parts, self.total, self.advective_scalar = parts, total, advective_scalar
        residual = self.dual(total)
        reactions = {}
        for key in ('k', 'epsilon'):
            fixed = self.fixed[key]
            reactions[key] = -float(self.one[key][fixed]@total[key][fixed])
        fixed_work = 0.
        for component, key in ((0, 'u'), (1, 'w')):
            ids = component*self.n+self.fixed[key]
            fixed_work += float(self.s['flow'][ids]@total['flow'][ids])
        bulk_whole = bulk['whole']
        budget = {}
        for key, source in (('k', 'k_source'), ('epsilon', 'epsilon_source')):
            scale = self.length*self.ref.A if key == 'k' else self.length*(
                self.ref.nu*self.ref.A/self.ref.H**2+self.ref.cd*self.ref.Kin**2/(self.ref.se*self.ref.H))
            adv, prescribed, physical, essential = (boundary[key+suffix] for suffix in
                ('_advection', '_diffusion_prescribed', '_diffusion_actual', '_diffusion_essential'))
            weak = adv+prescribed+reactions[key]-bulk_whole[source]
            free_action = float(self.one[key][self.free[key]]@total[key][self.free[key]])
            budget[key] = dict(storage=0., steady_defect_only=True, advective_outflow=adv,
                prescribed_outward_diffusion=prescribed, essential_outward_diffusion_reaction=reactions[key],
                essential_reaction_scope='inlet' if key == 'k' else 'combined inlet and H; corner reactions are not uniquely separable',
                owner_gradient_outward_diffusion_all=physical, owner_gradient_outward_diffusion_essential=essential,
                source=bulk_whole[source], reaction_balance=weak, free_test_action=free_action,
                assembly_identity_error_over_scale=abs(weak-free_action)/scale,
                mixed_physical_balance=adv+prescribed+essential-bulk_whole[source],
                all_owner_gradient_balance=adv+physical-bulk_whole[source], fixed_scale=scale,
                status='OBSERVED_REQUIRES_RESOLUTION_CHECK', actual_production_minus_destruction_not_constrained=True)
        energy_volume = (bulk_whole['kinetic_advection']+bulk_whole['stress_work']-
                         bulk_whole['pressure_divergence']-bulk_whole['roughness_work'])
        energy_free = sum(float(self.s['flow'][component*self.n+ids]@total['flow'][component*self.n+ids])
                          for component, ids in enumerate(self.ufree))
        energy_balance = energy_volume-boundary['traction_prescribed_work']-fixed_work
        energy_scale = self.length*self.ref.tau
        stress_identity = bulk_whole['stress_work']-bulk_whole['molecular_work']-bulk_whole['P']-bulk_whole['isotropic_k_divergence']
        budget['mechanical'] = dict(storage=0., steady_defect_only=True,
            kinetic_advection=bulk_whole['kinetic_advection'], kinetic_advective_boundary=boundary['kinetic_advection'],
            molecular_work=bulk_whole['molecular_work'], raw_turbulent_production_work=bulk_whole['P'],
            isotropic_k_divergence_work=bulk_whole['isotropic_k_divergence'], modified_pressure_divergence_work=-bulk_whole['pressure_divergence'],
            full_stress_work=bulk_whole['stress_work'], sampled_roughness_work=bulk_whole['roughness_work'],
            prescribed_boundary_traction_work=boundary['traction_prescribed_work'],
            essential_traction_reaction_work=fixed_work, essential_owner_gradient_traction_work=boundary['traction_essential_actual_work'],
            all_owner_gradient_traction_work=boundary['traction_actual_work'],
            reaction_balance=energy_balance, free_velocity_test_action=energy_free,
            assembly_identity_error_over_scale=abs(energy_balance-energy_free)/energy_scale,
            full_stress_raw_production_identity_error_over_scale=abs(stress_identity)/energy_scale,
            mixed_physical_balance=energy_volume-boundary['traction_prescribed_work']-boundary['traction_essential_actual_work'],
            all_owner_gradient_balance=energy_volume-boundary['traction_actual_work'], fixed_scale=energy_scale,
            status='OBSERVED_REQUIRES_RESOLUTION_CHECK')
        budget['mechanical']['advective_integration_identity_error_over_scale'] = abs(
            bulk_whole['kinetic_advection']-boundary['kinetic_advection']+
            bulk_whole['kinetic_energy_divergence'])/energy_scale
        budget['combined_mechanical_plus_k'] = dict(
            reaction_balance=energy_balance+budget['k']['reaction_balance'],
            molecular_work=bulk_whole['molecular_work'], destruction=bulk_whole['E'],
            raw_production_cancels_algebraically=True,
            pressure_and_kinetic_divergence_correction=-bulk_whole['pressure_divergence']+
                bulk_whole['isotropic_k_divergence']-bulk_whole['kinetic_energy_divergence'],
            boundary_advective_energy=boundary['kinetic_advection']+boundary['k_advection'],
            essential_scalar_outward_reaction=reactions['k'], sampled_roughness_work=bulk_whole['roughness_work'],
            prescribed_traction_work=boundary['traction_prescribed_work'], essential_traction_reaction_work=fixed_work,
            scope='Current actual stress production; no equilibrium or instantaneous P=E constraint')
        budget['mass'] = dict(boundary_outflow=boundary['mass'], integrated_divergence=bulk_whole['divergence'],
            assembly_identity_error=abs(boundary['mass']-bulk_whole['divergence'])/self.volume,
            weak_continuity=residual['continuity'], strong_divergence_RMS=math.sqrt(bulk_whole['strong_divergence_square']/self.volume))
        mutations = {}
        for name, differences in (
                ('missing_domain_momentum_convection', dict(flow=-parts['advection']['flow'])),
                ('missing_sampled_AWF_convection', dict(flow=-parts['awf_convection']['flow'])),
                ('missing_sampled_AWF_pressure', dict(flow=-parts['awf_pressure']['flow'])),
                ('missing_Q3_known_k_force', dict(flow=-parts['awf_known_k']['flow'])),
                ('wrong_known_k_sign', dict(flow=-2*parts['awf_known_k']['flow'])),
                ('nonconservative_scalar_advection', {key: advective_scalar[key]-parts['advection'][key] for key in ('k', 'epsilon')}),
                ('flipped_natural_boundary_sign', {key: -2*parts['natural'][key] for key in self.s})):
            delta = self.zeros(); delta.update(differences)
            changed = {key: total[key]+delta[key] for key in total}
            sensitivity = self.dual(delta)
            mutations[name] = dict(delta_dual=sensitivity, mutated_residual=self.dual(changed),
                status='DETECTABLE_AT_DECLARED_TOLERANCE' if any(sensitivity[key] > RESIDUAL_LIMITS[key] for key in sensitivity)
                else 'INACTIVE_IN_THIS_PHYSICAL_STATE; independent nonzero fixture still required')
        envelope = norms.report(dict(**ENVELOPES, shear=(.05, .10)))
        identities = [budget[key]['assembly_identity_error_over_scale'] for key in ('k', 'epsilon', 'mechanical')]
        identities += [budget['mechanical']['full_stress_raw_production_identity_error_over_scale'], budget['mass']['assembly_identity_error']]
        identities.append(budget['mechanical']['advective_integration_identity_error_over_scale'])
        return dict(residuals=residual, residual_limits=RESIDUAL_LIMITS,
                    steady_residual_status='MET' if all(residual[key] <= limit for key, limit in RESIDUAL_LIMITS.items()) else 'NOT_MET',
                    bulk=bulk, boundary=boundary, boundary_faces=face_rows, budgets=budget,
                    outlet=self.outlet_budget(face_rows),
                    accounting_identity_status='MET' if max(identities) <= 1e-10 else 'NOT_MET',
                    accounting_identity_max_error=max(identities), physical_envelope=envelope,
                    ground_shear_diagnostic=shear_diagnostic.report({}),
                    sampled_minima=minima, maximum_independent_cubic_relative_residual=max_root_error,
                    mutations=mutations, sparse_Riesz_nonzeros={key: matrix.nnz for key, matrix in self.masses.items()},
                    principal_admissibility=dict(status='NOT_MET', reason='Independent continuum principal certification is a separate required gate; sampled covariance is not a substitute'))

    def outlet_budget(self, face_rows):
        velocity = [row for row in face_rows if row['side'] == 'outlet' and row['field'] == 'u']
        work = {key: sum(row[key] for row in velocity) for key in (
            'traction_actual_work', 'traction_essential_actual_work',
            'prescribed_normal_work', 'prescribed_tangential_work')}
        rows = self.outlet_tangent_reactions(self.total['flow'])
        actions = {}
        for scope in ('outlet', 'shared_corner'):
            selected = [row for row in rows if row['scope'] == scope]
            actions[scope] = dict(
                unit_displacement_action=sum(self.one['u'][row['dof']]*row['value'] for row in selected),
                coefficient_work=sum(self.s['flow'][self.n+row['dof']]*row['value'] for row in selected))
        scalar = {}
        for field in ('k', 'epsilon'):
            faces = [row for row in face_rows if row['side'] == 'outlet' and row['field'] == field]
            scalar[field] = dict(
                advective_outflow=sum(row[field+'_advection'] for row in faces),
                prescribed_outward_diffusive_flux=sum(row[field+'_diffusion_prescribed'] for row in faces),
                raw_outward_diffusive_flux=sum(row[field+'_diffusion_actual'] for row in faces))
        return dict(policy=self.outlet_policy, identity=self.outlet_identity, normal=[1., 0.],
            component_order=['x', 'r'], mass=sum(row['mass'] for row in velocity),
            kinetic_advection=sum(row['kinetic_advection'] for row in velocity),
            sampled_minimum_normal_velocity=min(row['minimum_normal_velocity'] for row in velocity),
            prescribed_traction_integral=np.sum([row['prescribed_traction_integral'] for row in velocity], axis=0).tolist(),
            raw_stress_traction_integral=np.sum([row['actual_traction_integral'] for row in velocity], axis=0).tolist(),
            traction_work=work, scalar_fluxes=scalar,
            tangent_reaction=dict(sign='Positive outward traction action on the positive signed transverse test',
                definition='Full constrained steady momentum row after subtraction of applied loads',
                shared_corner_scope='One assembled row per shared corner, including its ground or top support',
                shared_corners=self.outlet_corners, rows=rows, actions=actions))

    def inlet_norms(self, norms, region, a, r, weight):
        inlet = self.ref.inlet(r)
        for name, scale in (('U', 1.), ('W', 1.), ('k', self.ref.Kin)):
            norms.add(region, name, (a[name]-inlet[name])/scale, weight)
        if a['outer']:
            norms.add(region, 'epsilon', a['E']/(self.ref.A/r)-1, weight)

    def boundary_checks(self):
        result, scalar_targets = {}, {}
        for key in ('k', 'epsilon'):
            target, assigned = np.zeros_like(self.s[key]), set()
            fixed_sides = (0,) if key == 'k' else (0, 2)
            scale = self.ref.Kin if key == 'k' else self.ref.A/self.ref.H
            for e, side in self.faces[key]:
                if side not in fixed_sides: continue
                self.clock()
                xl, xr, rl, rr = e['bounds']; a, w = leggauss(16)
                if side < 2:
                    cuts = self.cuts(rl, rr)
                    r = np.concatenate([(lo+hi)/2+(hi-lo)*a/2 for lo, hi in zip(cuts[:-1], cuts[1:])])
                    weights = np.concatenate([w*(hi-lo)/2 for lo, hi in zip(cuts[:-1], cuts[1:])])
                    x = np.full(len(r), xl if side == 0 else xr)
                    ex, er = np.array([xl if side == 0 else xr]*2), np.array([rl, rr])
                else:
                    x, r, weights = (xl+xr)/2+(xr-xl)*a/2, np.full(16, rl), w*(xr-xl)/2
                    ex, er = np.array([xl, xr]), np.array([rl, rl])
                def desired(xx, radius):
                    if key == 'k': return self.ref.inlet(radius)['k']
                    return (self.trace(e, 'k', xx)**1.5/self.ref.cl/self.ref.H if side == 2 else self.ref.A/radius)
                modes = EDGES[side]; b = basis(e, key, x, r)[0]
                endpoints = desired(ex, er)/np.asarray(e[SIGNS[key]])[list(modes[:2])]
                interior = b[:, list(modes[2:])]
                remainder = desired(x, r)-b[:, list(modes[:2])]@endpoints
                coefficients = np.r_[endpoints, np.linalg.solve(interior.T@(weights[:, None]*interior), interior.T@(weights*remainder))]
                for mode, value in zip(modes, coefficients):
                    gid = e[key][mode]
                    require(gid not in assigned or abs(target[gid]-value) <= 2e-11*scale,
                            'independent essential scalar corner incompatibility')
                    target[gid] = value; assigned.add(gid)
            require(assigned == set(self.fixed[key]), 'incomplete independently reconstructed scalar essential lift')
            error = float(np.max(np.abs((self.s[key]-target)[self.fixed[key]]))/scale)
            scalar_targets[key] = target
            result[key+'_essential_coefficient_error'] = dict(value=error, limit=1e-10, status='MET' if error <= 1e-10 else 'NOT_MET',
                projection='Endpoint constrained physical Gauss16 per fixed material panel; basis independent of native operators')
        normals = Norms(); trace_max = 0.
        for e, side in self.faces['u']:
            x, r, weight = self.face(e, side, endpoints=True)
            a = self.sample(e, x, r)
            if side == 0:
                inlet = self.ref.inlet(r)
                normals.add('inlet', 'U_representation', a['U']-inlet['U'], weight)
                normals.add('inlet', 'W_essential', a['W'], weight)
            if side == 1 and self.outlet_policy == ZERO_TRANSVERSE:
                normals.add('outlet', 'W_essential', a['W'], weight)
                normals.add('outlet', 'Wr_essential', a['G'][:, 2, 2], weight)
            if side == 2:
                normals.add('ground', 'U_essential', a['U'], weight)
                normals.add('ground', 'W_essential', a['W'], weight)
            if side == 3:
                normals.add('top', 'W_essential', a['W'], weight)
        for e, side in self.faces['epsilon']:
            if side != 2: continue
            x, r, weight = self.face(e, side, endpoints=True)
            a = self.sample(e, x, r)
            exact = a['K']**1.5/self.ref.cl/self.ref.H
            trace_max = max(trace_max, float(np.max(np.abs(a['E']-exact)/(self.ref.A/self.ref.H))))
        result['pointwise_boundary_representation'] = normals.report(
            dict(U_essential=1e-10, W_essential=1e-10, Wr_essential=1e-10))
        result['homogeneous_velocity_boundary_check'] = dict(status='MET' if all(value.get('status', 'MET') == 'MET'
            for value in result['pointwise_boundary_representation'].values()) else 'NOT_MET')
        result['dependent_H_pointwise_representation_error'] = trace_max
        result['dependent_H_note'] = 'Algebraic condition is the current endpoint-constrained projection; nonzero pointwise projection error is a spatial error, not an extra boundary equation.'
        if 'essential_flow_lift' in self.data:
            target = np.asarray(self.data['essential_flow_lift'])
            require(target.shape == self.s['flow'].shape and np.all(np.isfinite(target)), 'invalid exported essential velocity lift')
            error = max(float(np.max(np.abs(self.s['flow'][component*self.n+ids]-target[component*self.n+ids])))
                        for component, ids in ((0, self.fixed['u']), (1, self.fixed['w'])))
            result['velocity_essential_coefficient_error'] = dict(value=error, limit=1e-10,
                status='MET' if error <= 1e-10 else 'NOT_MET', scope='Actual coefficients versus exported fresh native lift; historical pointwise comparison above is independent')
        else:
            result['velocity_essential_coefficient_error'] = dict(status='NOT_MET', reason='Native boundary projection point inventory/fresh lift is not exported; pointwise representation evidence is separate')
        result['status'] = 'MET' if all(value['status'] == 'MET' for value in result.values() if isinstance(value, dict) and 'status' in value) else 'NOT_MET'
        return result

    def interfaces(self):
        norms, h_rows = Norms(), []
        for owners in self.internal:
            self.clock()
            (e, side), (other, other_side) = owners
            x, r, weight = self.face(e, side, endpoints=True)
            a, b = self.sample(e, x, r), self.sample(other, x, r)
            nx, nr = NORMALS[side]
            una, unb = nx*a['U']+nr*a['W'], nx*b['U']+nr*b['W']
            tk = (-a['Dk']*(nx*a['gradk'][:, 0]+nr*a['gradk'][:, 1])+
                   b['Dk']*(nx*b['gradk'][:, 0]+nr*b['gradk'][:, 1]))
            traction = np.column_stack((nx*(a['V'][:, 0, 0]-a['Pi']-b['V'][:, 0, 0]+b['Pi'])+nr*(a['V'][:, 0, 2]-b['V'][:, 0, 2]),
                nx*(a['V'][:, 2, 0]-b['V'][:, 2, 0])+nr*(a['V'][:, 2, 2]-a['Pi']-b['V'][:, 2, 2]+b['Pi'])))
            group = 'H' if a['outer'] != b['outer'] else 'interior'
            norms.add(group, 'mass_shared', una-unb, weight)
            norms.add(group, 'k_advection_shared', (una*a['k']-unb*b['k'])/self.ref.Kin, weight)
            norms.add(group, 'k_owner_diffusion_jump', tk/self.ref.A, weight)
            norms.add(group, 'total_traction_owner_jump', traction/self.ref.tau, weight)
            if a['outer'] and b['outer']:
                ej = (-a['De']*(nx*a['grade'][:, 0]+nr*a['grade'][:, 1])+
                       b['De']*(nx*b['grade'][:, 0]+nr*b['grade'][:, 1]))
                norms.add(group, 'epsilon_owner_diffusion_jump', ej/(self.ref.A/self.ref.H), weight)
            if group == 'H':
                h_rows.append(dict(owners=[e['geometry'], other['geometry']], side=SIDES[side],
                    k_diffusion_jump_integral=float(weight@tk), total_traction_jump_integral=(weight@traction).tolist(),
                    dependent_epsilon_trace_jump_MAX=float(np.max(np.abs(a['E']-b['E'])))/(self.ref.A/self.ref.H)))
        return dict(norms=norms.report({}), H_faces=h_rows,
            assembled_common_flux='One global conforming test space; no artificial shared gradient is substituted for either physical owner.',
            physical_gradient_status='OBSERVED_REQUIRES_RESOLUTION_CHECK')

    @staticmethod
    def observable(a):
        result = dict(U=a['U'], W=a['W'], k=a['k'], Pi=a['Pi'], mechanical_pressure=a['Pi']-2*a['k']/3,
                    R=a['R'][:, (0, 0, 0, 1, 1, 2), (0, 1, 2, 1, 2, 2)],
                    G=a['G'][:, (0, 0, 2, 2), (0, 2, 0, 2)], traction=a['V'][:, 0, 2],
                    K=a['K'], cutoff=a['cutoff'], chi=a['chi'], P=a['P'], production_over_destruction=a['P']/a['E'])
        result['epsilon' if a['outer'] else 'Ew'] = a['E']
        return result

    def profiles(self, *, detailed=False, stations=STATIONS):
        path = Path(__file__).resolve().parents[1]/'cao-tamura-fig3/rough-fig3.csv'
        fig = sorted(set(float(row['y_mm'])/40-self.ref.d0 for row in csv.DictReader(path.open())))
        extra = [0., self.ref.H, self.ref.h, .4*self.ref.H, .75*self.ref.H, .0875, self.ref.top, *fig]
        rows, norms, ground = [], Norms(), []
        for station in stations:
            self.clock()
            owners = [e for e in self.elements if e['bounds'][0]-1e-12 <= station <= e['bounds'][1]+1e-12]
            widths = len({(e['bounds'][0], e['bounds'][1]) for e in owners})
            require(widths in (1, 2), 'station has absent or ambiguous x-side owners')
            for e in owners:
                xl, xr, rl, rr = e['bounds']
                cutoff = self.ref.cutoff(self.trace(e, 'k', np.array([station])))[0]
                cuts = self.cuts(rl, rr, (cutoff, self.ref.inlet_cutoff, *extra))
                gauss, weight = self.gauss
                for lo, hi in zip(cuts[:-1], cuts[1:]):
                    r = (lo+hi)/2+(hi-lo)*gauss/2
                    w = weight*(hi-lo)/(2*widths)
                    a = self.sample(e, np.full(len(r), station), r)
                    region = 'outer' if a['outer'] else 'inner'
                    self.inlet_norms(norms, f'x={station:g}/{region}', a, r, w)
                    self.inlet_norms(norms, f'x={station:g}/whole', a, r, w)
                # Finite MAX inventory includes both exact owner endpoints and
                # inward offsets at EVERY material/cutoff/native face break.
                observed = sorted(set(cuts+[(a+b)/2 for a, b in zip(cuts[:-1], cuts[1:])] +
                                      [a+1e-8*(b-a) for a, b in zip(cuts[:-1], cuts[1:])] +
                                      [b-1e-8*(b-a) for a, b in zip(cuts[:-1], cuts[1:])]))
                r = np.asarray(observed); x = np.full(len(r), station)
                a = self.sample(e, x, r); fields = self.observable(a); inlet = self.ref.inlet(r)
                if detailed:
                    fields.update(strong_divergence=a['div'],
                                  raw_traction_e_x=np.column_stack((a['V'][:,0,0]-a['Pi'],a['V'][:,2,0])))
                region = 'outer' if a['outer'] else 'inner'
                self.inlet_norms(norms, f'x={station:g}/{region}', a, r, np.zeros(len(r)))
                self.inlet_norms(norms, f'x={station:g}/whole', a, r, np.zeros(len(r)))
                for j, radius in enumerate(r):
                    rows.append(dict(x=station, r=float(radius), z=float(radius+self.ref.d0), x_mm=40*station,
                        z_mm=40*(float(radius)+self.ref.d0), owner=e['geometry'], bounds=e['bounds'], outer=a['outer'],
                        values={key: np.asarray(value[j]).tolist() for key, value in fields.items()},
                        retained_inlet={key: float(value[j]) for key, value in inlet.items() if key != 'G'}))
                if rl == 0:
                    ground.append(dict(x=station, owner=e['geometry'], upward_shear=float(a['V'][0, 0, 2]),
                                       scaled_deviation=float(a['V'][0, 0, 2]/self.ref.tau-1),
                                       scaled_deviation_vs_reference=float((a['V'][0, 0, 2]-self.ref.reference_shear)/self.ref.tau)))
        result = norms.report(ENVELOPES)
        return dict(stations=stations, norms=result, rows=rows, ground_stations=ground,
            status='MET' if all(v.get('status') == 'MET' for v in result.values()) and all(abs(v['scaled_deviation']) <= .10 for v in ground) else 'NOT_MET',
            comparator=('Exact-inner and retained outer inlet. Ground shear gate uses tau0.'
                        if self.ref.exact_inner else 'Retained imposed historical inlet.'),
            sampling='Physical Gauss integration per native/material/cutoff/Fig3 panel; MAX also includes exact endpoints, panel midpoints and relative 1e-8 inward offsets, both original owners.',
            physical_7mm=dict(z_code=.175, r_code=.0875, region='inner', conversion_mm_per_code=40), fig3_source=identity(path))

    def compare(self, finer):
        norms, count = Norms(), 0
        xbreaks = sorted(set(e['bounds'][i] for collection in (self.elements, finer.elements)
                            for e in collection for i in (0, 1)))
        for e in self.elements:
            xl, xr, rl, rr = e['bounds']
            for fine in finer.elements:
                fl, fr, fb, ft = fine['bounds']
                left, right, bottom, top = max(xl, fl), min(xr, fr), max(rl, fb), min(rr, ft)
                if right-left <= 1e-12 or top-bottom <= 1e-12: continue
                self.clock(); count += 1
                bounds = (left, right, bottom, top)
                x, r, weight = self.rule(e, bounds, lambda x: finer.ref.cutoff(finer.trace(fine, 'k', x)),
                                         finer.cutoff_x_splits(fine,bounds))
                # Include both side limits for all common subcell endpoints.
                # Each pair has a common nonempty interior, defining the side;
                # no nearest-cell search averages opposite gradients.
                xx, radial, ww = list(x), list(r), list(weight)
                for xend in (left, right):
                    cuts = self.cuts(bottom, top, (self.ref.cutoff(self.trace(e, 'k', np.array([xend])))[0],
                        finer.ref.cutoff(finer.trace(fine, 'k', np.array([xend])))[0]))
                    xx.extend([xend]*len(cuts)); radial.extend(cuts); ww.extend([0.]*len(cuts))
                x, r, weight = np.asarray(xx), np.asarray(radial), np.asarray(ww)
                a, b = self.sample(e, x, r), finer.sample(fine, x, r)
                require(a['outer'] == b['outer'], 'comparison crossed distinct epsilon owner domains')
                aa, bb = self.observable(a), self.observable(b)
                for key in NUMERICAL:
                    if key == 'epsilon' and not a['outer']: continue
                    scale = self.ref.Kin if key == 'k' else self.ref.tau if key in ('R', 'Pi', 'traction') else self.ref.A/r if key == 'epsilon' else 1.
                    delta = (aa[key]-bb[key])/scale
                    for region in ('whole', 'outer' if a['outer'] else 'inner'):
                        norms.add(region, key, delta, weight)
                for station in STATIONS:
                    if station < left-1e-12 or station > right+1e-12: continue
                    multiplicity = 2 if -10 < station < 15 and any(abs(station-edge) < 1e-12 for edge in xbreaks) else 1
                    cuts = self.cuts(bottom, top, (self.ref.cutoff(self.trace(e, 'k', np.array([station])))[0],
                        finer.ref.cutoff(finer.trace(fine, 'k', np.array([station])))[0]))
                    nodes, weights = self.gauss
                    r = np.concatenate([(lo+hi)/2+(hi-lo)*nodes/2 for lo, hi in zip(cuts[:-1], cuts[1:])])
                    weight = np.concatenate([weights*(hi-lo)/(2*multiplicity) for lo, hi in zip(cuts[:-1], cuts[1:])])
                    r, weight = np.r_[r, cuts], np.r_[weight, np.zeros(len(cuts))]
                    x = np.full(len(r), station)
                    a, b = self.sample(e, x, r), finer.sample(fine, x, r)
                    aa, bb = self.observable(a), self.observable(b)
                    for key in NUMERICAL:
                        if key == 'epsilon' and not a['outer']: continue
                        scale = self.ref.Kin if key == 'k' else self.ref.tau if key in ('R', 'Pi', 'traction') else self.ref.A/r if key == 'epsilon' else 1.
                        delta = (aa[key]-bb[key])/scale
                        for region in ('whole', 'outer' if a['outer'] else 'inner'):
                            norms.add(f'x={station:g}/'+region, key, delta, weight)
        require(count > 0, 'no common refinement cells')
        result = norms.report(NUMERICAL)
        return dict(status='MET' if all(v['status'] == 'MET' for v in result.values()) else 'NOT_MET',
                    norms=result, common_owner_subcells=count, common_owner_sampling=True,
                    gradient_scale='Uref/L0=1 in code coordinates; L0=40mm',
                    scope='Single common-physical h/p/q comparison; temporal order, whole principal and finite perturbations remain separate gates')


def load_native(path):
    # Match native coefficient/checkpoint bytes when a double is emitted -0.
    data = json.loads(Path(path).read_text(),
                      parse_int=lambda token: -0.0 if token == '-0' else int(token))
    if 'state' in data:
        require(data.get('normal_main') is True and data.get('factory') == 'WJNativeCoupledV2',
                'input is not a registered normal-main coupled run')
        return data['state'], {key: value for key, value in data.items() if key != 'state'}
    return data, dict(scope='Standalone exported native scientific state')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    parser.add_argument('--compare-finer', type=Path)
    parser.add_argument('--order', type=int, default=20)
    args = parser.parse_args()
    require(sys.platform.startswith('linux'), 'numerical oracle execution is Linux-only')
    require(math.isfinite(args.seconds) and args.seconds > 0, 'explicit finite wall cap required')
    require(20 <= args.order <= 48, 'declared independent Gauss order must be in [20,48]')
    require(not args.out.exists(), 'immutable oracle output already exists')
    # Resource policy is additional to the parent runner, never permission to
    # exceed its allocation. Preserve any stricter inherited virtual-memory cap.
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = 8*1024**3
    if soft != resource.RLIM_INFINITY: cap = min(cap, soft)
    if hard != resource.RLIM_INFINITY: cap = min(cap, hard)
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    began = time.monotonic()
    report = dict(schema=SCHEMA, status='NOT_MET', full_H27_qualification=False,
                  source=identity(__file__), basis_source=identity(fixtures.__file__),
                  native=identity(args.native), table=identity(args.table),
                  quadrature=[args.order, args.order], wall_cap_seconds=args.seconds,
                  memory_cap_bytes=cap, workers=1,
                  unresolved_gates=['Whole-state admissibility certification', 'Temporal and finite perturbation qualification',
                                    'Stationary uniqueness/existence is not inferred from one residual evaluation'])
    try:
        reference = load_reference(args.table)
        data, metadata = load_native(args.native)
        report['native_run_metadata'] = metadata
        report['native_source_identity'] = data['source_identity']
        oracle = PhysicalOracle(data, reference, args.order, began+args.seconds)
        report['state_sha256'] = data['state_sha256']
        report['core'] = oracle.assemble()
        report['boundary_checks'] = oracle.boundary_checks()
        report['owner_interfaces'] = oracle.interfaces()
        report['profiles'] = oracle.profiles()
        if args.compare_finer:
            finer_data, finer_metadata = load_native(args.compare_finer)
            finer = PhysicalOracle(finer_data, reference, args.order, began+args.seconds)
            finer_core = finer.assemble()
            comparison = oracle.compare(finer)
            physical_budget_changes = {}
            for key in ('k', 'epsilon', 'mechanical'):
                coarse_budget = report['core']['budgets'][key]
                fine_budget = finer_core['budgets'][key]
                physical_budget_changes[key] = dict(
                    coarse_mixed_physical_balance=coarse_budget['mixed_physical_balance'],
                    finer_mixed_physical_balance=fine_budget['mixed_physical_balance'],
                    finer_all_owner_gradient_balance=fine_budget['all_owner_gradient_balance'],
                    balance_difference_over_fixed_scale=(fine_budget['mixed_physical_balance']-
                        coarse_budget['mixed_physical_balance'])/coarse_budget['fixed_scale'],
                    absolute_mixed_balance_reduced=abs(fine_budget['mixed_physical_balance']) < abs(coarse_budget['mixed_physical_balance']),
                    status='OBSERVED; no field-comparison result substitutes for flux-budget convergence')
            report['refinement'] = dict(input=identity(args.compare_finer), native_run_metadata=finer_metadata,
                **comparison, finer_current_residuals=finer_core['residuals'],
                finer_steady_residual_status=finer_core['steady_residual_status'],
                finer_accounting_identity_status=finer_core['accounting_identity_status'],
                physical_budget_changes=physical_budget_changes)
            if finer_core['steady_residual_status'] != 'MET' or finer_core['accounting_identity_status'] != 'MET':
                report['refinement']['status'] = 'NOT_MET'
        else:
            report['refinement'] = dict(status='NOT_MET', reason='No independent finer h/p/q state was supplied')
        core = report['core']
        report['core_steady_equations_status'] = 'MET' if core['steady_residual_status'] == 'MET' and core['accounting_identity_status'] == 'MET' else 'NOT_MET'
        envelope = core['physical_envelope']
        report['physical_development_status'] = 'MET' if all(value['status'] == 'MET' for value in envelope.values()) and report['profiles']['status'] == 'MET' else 'NOT_MET'
        report['steady_residual_boundary_and_field_refinement_status'] = 'MET' if all(value == 'MET' for value in (
            report['core_steady_equations_status'], report['boundary_checks']['status'], report['refinement']['status'])) else 'NOT_MET'
        report['same_inlet_resolved_steady_evidence_status'] = 'NOT_MET'
        report['unresolved_gates'].append('Owner-gradient physical budget convergence must be assessed with the recorded resolution evidence; weak reaction identities alone cannot close it')
        report['reason'] = 'Independent evidence recorded. Overall NOT_MET retains the separate unresolved qualification gates; no positive default.'
    except Exception as error:
        report['exception'] = type(error).__name__
        report['reason'] = str(error)
    report['elapsed_seconds'] = time.monotonic()-began
    report['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: value for key, value in report.items() if key in
        ('status', 'reason', 'exception', 'elapsed_seconds', 'maximum_resident_bytes',
         'core_steady_equations_status', 'physical_development_status', 'same_inlet_resolved_steady_evidence_status')}))
    if 'exception' in report:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
