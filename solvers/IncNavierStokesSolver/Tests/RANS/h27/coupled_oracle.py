#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux-only current coupled MMS step or true steady audit.

Native coefficient maps are data; Jacobi basis, physical Gauss integration,
constitutive evaluation and history reconstruction are independent. Continuous
MMS loads are integrated by parts from ANALYTIC fields with all exterior terms;
they are never formed as a native discrete matrix applied to target coefficients.
An optional dense block-direct nonlinear reference starts from the prior state.
Top-level steady_mode=true instead uses alpha=0, frozen t=0 continuous forcing,
and current advection; its optional direct reference starts from the original
fresh seed and refreshes advection before each block. Neither mode qualifies
temporal order, physical development or the full H27 gate.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'h26'))
import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.chebyshev import chebval
from fixtures import modified
import coupled_mms_check
import reference_oracle

EDGES = ((0, 4, 8, 12), (1, 5, 9, 13), (0, 1, 2, 3), (4, 5, 6, 7))
NORMALS = ((-1., 0.), (1., 0.), (0., -1.), (0., 1.))
LIMITS = dict(momentum=1e-9, continuity=1e-9, gauge=1e-9, k=1e-8, epsilon=1e-8)
SIGN = dict(u='us', p='ps', k='ks', epsilon='epsilon_signs')


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def state(packet, sizes):
    result = {key: np.asarray(packet[key], dtype=float) for key in ('flow', 'k', 'epsilon')}
    for key, size in zip(result, (2*sizes[0]+sizes[1], sizes[2], sizes[3])):
        require(result[key].shape == (size,) and np.all(np.isfinite(result[key])),
                'invalid '+key+' native coefficients/history')
    return result


def state_hash(s):
    digest = hashlib.sha256()
    for key in ('flow', 'k', 'epsilon'):
        digest.update(np.asarray(s[key], dtype='<f8').tobytes())
    return digest.hexdigest()


def combine(a, wa, b, wb):
    return {key: wa*a[key]+wb*b[key] for key in a}


def basis(e, key, x, r):
    degree = 2 if key == 'p' else 3
    xl, xr, rl, rr = e['bounds']; signs = np.asarray(e[SIGN[key]])
    bx, dx = modified(degree, 2*(x-xl)/(xr-xl)-1)
    bz, dz = modified(degree, 2*(r-rl)/(rr-rl)-1)
    n = (degree+1)**2
    return ((bz[:, :, None]*bx[:, None, :]).reshape(len(x), n)*signs,
            (bz[:, :, None]*dx[:, None, :]).reshape(len(x), n)*signs*2/(xr-xl),
            (dz[:, :, None]*bx[:, None, :]).reshape(len(x), n)*signs*2/(rr-rl))


class Oracle:
    def __init__(self, data, mms, nx, nr, deadline):
        self.data, self.m, self.nx, self.nr, self.deadline = data, mms, nx, nr, deadline
        require(data.get('mms_variant', 'Historical') == mms.variant,
                'native state and independent MMS scalar variants differ')
        self.nu, self.np, self.nk, self.ne = map(int, data['sizes'])
        self.sizes = data['sizes']; self.elements = data['elements']; self.g = data['geometry']
        self.ufree = [np.asarray(v, dtype=int) for v in data['free_velocity_ids']]
        self.kfree = np.asarray(data['k_free_ids'], dtype=int)
        self.efree = np.asarray(data['epsilon_free_ids'], dtype=int)
        self.vfree = np.r_[self.ufree[0], self.nu+self.ufree[1]]
        self.volume = (self.g['right']-self.g['left'])*self.g['top']
        self.evolume = (self.g['right']-self.g['left'])*(self.g['top']-mms.H)
        require(data['manufactured'] is True and data['gauge'] is True,
                'bounded oracle requires all-essential-velocity MMS with pressure gauge')
        for name, expected in dict(d0=mms.d0, H=mms.H, h=mms.h, left=0., right=1., top=.5).items():
            require(abs(self.g[name]-expected) <= 1e-12, 'MMS geometry mismatch: '+name)
        self.owner, self.faces = {}, {'u': [], 'k': [], 'epsilon': []}
        visited = {key: set() for key in ('u', 'p', 'k', 'epsilon')}
        require(len({e['geometry'] for e in self.elements}) == len(self.elements), 'duplicate geometry')
        for e in self.elements:
            xl, xr, rl, rr = e['bounds']
            require(xr > xl and rr > rl and not (rl < mms.H-1e-12 and rr > mms.H+1e-12),
                    'invalid rectangular/native H element inventory')
            require(xl >= -1e-12 and xr <= 1+1e-12 and rl >= -1e-12 and rr <= .5+1e-12,
                    'element lies outside declared MMS rectangle')
            outer = rl >= mms.H-1e-12
            require(bool(e['epsilon']) == outer, 'outer epsilon ownership differs')
            for key, size in (('u', self.nu), ('p', self.np), ('k', self.nk), ('epsilon', self.ne)):
                ids, signs = e[key], e[SIGN[key]]
                expected = 9 if key == 'p' else 16
                if key == 'epsilon' and not outer:
                    require(not signs, 'inner geometry has epsilon signs'); continue
                require(len(ids) == len(signs) == expected and len(set(ids)) == expected and
                        all(0 <= v < size for v in ids) and all(v in (-1., 1.) for v in signs),
                        'invalid independent signed '+key+' map')
                visited[key].update(ids)
            owners = [a for a in self.elements if abs(a['bounds'][2]-mms.H) < 1e-12 and
                      abs(a['bounds'][0]-xl) < 1e-12 and abs(a['bounds'][1]-xr) < 1e-12]
            require(len(owners) == 1, 'missing/ambiguous outer H owner')
            self.owner[e['geometry']] = owners[0]
            for key in self.faces:
                if key == 'epsilon' and not outer:
                    continue
                bottom = mms.H if key == 'epsilon' else 0.
                for side, exterior in enumerate((abs(xl) < 1e-12, abs(xr-1) < 1e-12,
                                                  abs(rl-bottom) < 1e-12, abs(rr-.5) < 1e-12)):
                    if exterior:
                        self.faces[key].append((e, side))
        for key, size in (('u', self.nu), ('p', self.np), ('k', self.nk), ('epsilon', self.ne)):
            require(visited[key] == set(range(size)), 'incomplete '+key+' global map')
        area = lambda e: (e['bounds'][1]-e['bounds'][0])*(e['bounds'][3]-e['bounds'][2])
        require(abs(sum(area(e) for e in self.elements)-self.volume) < 1e-11*self.volume and
                abs(sum(area(e) for e in self.elements if e['epsilon'])-self.evolume) < 1e-11*self.evolume,
                'full/outer geometry inventories do not cover their declared volumes')
        self.essential = {}
        for key, free, size, sides in (('u', self.ufree[0], self.nu, (0, 1, 2, 3)),
                                      ('k', self.kfree, self.nk, (0,)),
                                      ('epsilon', self.efree, self.ne, (0, 2))):
            fixed = {e[key][j] for e, side in self.faces[key] if side in sides for j in EDGES[side]}
            require(len(set(free)) == len(free) and set(free) == set(range(size))-fixed,
                    'wrong independently derived '+key+' essential mask')
            self.essential[key] = np.asarray(sorted(fixed), dtype=int)
        require(set(self.ufree[1]) == set(self.ufree[0]), 'MMS velocity component masks differ')
        self.static()

    def clock(self):
        if time.monotonic() >= self.deadline:
            raise TimeoutError('independent coupled oracle reached its admitted elapsed cap')

    def cutoff(self, K):
        K = np.asarray(K)
        require(np.all((K >= self.m.lo) & (K <= self.m.hi)), 'current K outside qualified interval')
        return self.m.H*np.exp(chebval((2*K-self.m.lo-self.m.hi)/(self.m.hi-self.m.lo), self.m.coefficients))

    def aux(self, K, r):
        m = self.m; root = np.sqrt(K); rough = m.h*root/(70*m.nu)
        exponent = np.maximum(.5-.4*rough**.7, 1-.79*rough**(-.28))
        zv = -10.7*m.nu*np.expm1(exponent*np.log(rough))/root
        require(np.all(zv < 0), 'reference roughness branch changed')
        return m.kappa*np.sqrt(K/m.qe)*(r-zv)

    def trace(self, e, key, x, s, derivative=False):
        owner = self.owner[e['geometry']]
        b = basis(owner, key, np.asarray(x), np.full(len(x), self.m.H))[int(derivative)]
        values = s['flow'][2*self.nu+np.asarray(owner['p'])] if key == 'p' else s[key][owner[key]]
        return b@values

    def rule(self, e, s=None, when=None):
        ax, wx = leggauss(self.nx); ar, wr = leggauss(self.nr)
        xl, xr, rl, rr = e['bounds']; xs, rs, weights = [], [], []
        for a, wa in zip(ax, wx):
            x = (xl+xr)/2+(xr-xl)*a/2
            cuts = [rl, rr, .4*self.m.H, .75*self.m.H, self.m.h, self.m.H]
            if s is not None:
                cuts.append(float(self.cutoff(self.trace(e, 'k', np.array([x]), s)[0])))
            if when is not None:
                cuts.append(float(self.cutoff(self.m.functions['K'](x, self.m.H, when)[0])))
            cuts = sorted(set([rl, rr]+[v for v in cuts if rl < v < rr]))
            for lo, hi in zip(cuts[:-1], cuts[1:]):
                xs.extend(np.full(self.nr, x)); rs.extend((lo+hi)/2+(hi-lo)*ar/2)
                weights.extend(wa*wr*(xr-xl)*(hi-lo)/4)
        return np.asarray(xs), np.asarray(rs), np.asarray(weights)

    def face(self, e, side):
        a, w = leggauss(max(self.nx, self.nr)); xl, xr, rl, rr = e['bounds']
        if side < 2:
            return np.full(len(a), xl if side == 0 else xr), (rl+rr)/2+(rr-rl)*a/2, w*(rr-rl)/2
        return (xl+xr)/2+(xr-xl)*a/2, np.full(len(a), rl if side == 2 else rr), w*(xr-xl)/2

    def exact(self, x, r, when, outer):
        return {name: np.vstack([np.broadcast_to(np.asarray(v), x.shape) for v in f(x, r, when)])
                for name, f in self.m.functions.items() if name != 'epsilon' or outer}

    def constitutive(self, k, E, G, auxiliary, r, need_mu=False):
        m = self.m; count = len(r); R = np.zeros((count, 3, 3)); mu = np.full(count, m.nu)
        coordinate = np.clip((r/m.H-.4)/.35, 0., 1.)
        chi = coordinate**3*(10-15*coordinate+6*coordinate**2)
        M = G+G.transpose(0, 2, 1)
        require(np.min(k) > 0 and np.min(E) > 0, 'nonpositive current scalar/destruction')
        for q in range(count):
            if q % 256 == 0: self.clock()
            iso = 2*k[q]*np.eye(3)/3
            iso[0, 2] = iso[2, 0] = -auxiliary[q]*M[q, 0, 2]
            iso[1, 2] = iso[2, 1] = -auxiliary[q]*M[q, 1, 2]
            R[q] = (1-chi[q])*iso
            if chi[q] > 0:
                R[q] += chi[q]*reference_oracle.full_wj(k[q], E[q], G[q])
                if need_mu:
                    st = (k[q]/E[q])*(M[q]/2-np.trace(G[q])*np.eye(3)/3)
                    ot = (k[q]/E[q])*(G[q]-G[q].T)/2; iw = float(np.trace(ot@ot))
                    root = reference_oracle.physical_root(float(np.trace(st@st)), iw)
                    mu[q] += chi[q]*3*root/(5*(root*root-2*iw))*k[q]*k[q]/E[q]
            mu[q] += (1-chi[q])*auxiliary[q]
        V = m.nu*M-R+(2*k/3)[:, None, None]*np.eye(3)
        return dict(V=V, production=-np.einsum('qij,qij->q', R, G), mu=mu, M=M)

    def sampled(self, e, x, r, s, frozen=None, need_mu=False):
        bu, ux, ur = basis(e, 'u', x, r); ui = np.asarray(e['u'])
        u, v = bu@s['flow'][ui], bu@s['flow'][self.nu+ui]
        G = np.zeros((len(x), 3, 3))
        G[:, 0, 0], G[:, 0, 2] = ux@s['flow'][ui], ur@s['flow'][ui]
        G[:, 2, 0], G[:, 2, 2] = ux@s['flow'][self.nu+ui], ur@s['flow'][self.nu+ui]
        scalar = frozen if frozen is not None else s
        bk, kx, kr = basis(e, 'k', x, r); k = bk@scalar['k'][e['k']]
        outer = bool(e['epsilon'])
        K = self.trace(e, 'k', x, scalar)
        E = basis(e, 'epsilon', x, r)[0]@scalar['epsilon'][e['epsilon']] if outer else K**1.5/self.m.cl/np.maximum(r, self.cutoff(K))
        auxiliary = np.zeros(len(x)) if outer else self.aux(K, r)
        law = self.constitutive(k, E, G, auxiliary, r, need_mu)
        law.update(u=u, v=v, G=G, k=k, E=E, K=K, Dk=self.m.nu+self.m.cd*k*k/(self.m.sk*E),
                   De=self.m.nu+self.m.cd*k*k/(self.m.se*E) if outer else np.zeros(len(x)))
        return law

    def exact_law(self, x, r, when, outer):
        j = self.exact(x, r, when, outer); k, K = j['k'][0], j['K'][0]
        E = j['epsilon'][0] if outer else K**1.5/self.m.cl/np.maximum(r, self.cutoff(K))
        G = np.zeros((len(x), 3, 3)); G[:, 0, 0], G[:, 0, 2] = j['u'][1:3]
        G[:, 2, 0], G[:, 2, 2] = j['w'][1:3]
        law = self.constitutive(k, E, G, np.zeros(len(x)) if outer else self.aux(K, r), r)
        law.update(j=j, E=E, Dk=self.m.nu+self.m.cd*k*k/(self.m.sk*E),
                   De=self.m.nu+self.m.cd*k*k/(self.m.se*E) if outer else np.zeros(len(x)))
        return law

    def zero(self):
        return dict(flow=np.zeros(2*self.nu+self.np), k=np.zeros(self.nk), epsilon=np.zeros(self.ne))

    def static(self):
        n, p = self.nu, self.np
        self.M = np.zeros((n, n)); self.Q = np.zeros((p, p)); self.Mk = np.zeros((self.nk, self.nk))
        self.Me = np.zeros((self.ne, self.ne)); self.Re = np.zeros_like(self.Me); self.Eerror = np.zeros_like(self.Me)
        self.D = np.zeros((p, 2*n)); self.C = np.zeros((2*n, p)); self.Ak = np.zeros((2*n, self.nk))
        self.cross = np.zeros((p, self.nk)); self.mean = np.zeros(p)
        for e in self.elements:
            self.clock(); x, r, w = self.rule(e)
            b, dx, dr = basis(e, 'u', x, r); P, _, _ = basis(e, 'p', x, r)
            k, _, _ = basis(e, 'k', x, r); ui, pi, ki = map(np.asarray, (e['u'], e['p'], e['k']))
            self.M[np.ix_(ui, ui)] += b.T@(w[:, None]*b)
            self.Q[np.ix_(pi, pi)] += P.T@(w[:, None]*P)
            self.Mk[np.ix_(ki, ki)] += k.T@(w[:, None]*k)
            self.cross[np.ix_(pi, ki)] += P.T@(w[:, None]*k)
            self.mean[pi] += P.T@w
            self.D[np.ix_(pi, ui)] += P.T@(w[:, None]*dx)
            self.D[np.ix_(pi, n+ui)] += P.T@(w[:, None]*dr)
            if e['epsilon']:
                ei = np.asarray(e['epsilon']); eb = basis(e, 'epsilon', x, r)[0]
                self.Me[np.ix_(ei, ei)] += eb.T@(w[:, None]*eb)
                self.Re[np.ix_(ei, ei)] += eb.T@((w*(self.m.Ain/r)**2)[:, None]*eb)
                self.Eerror[np.ix_(ei, ei)] += eb.T@((w*(r/self.m.Ain)**2)[:, None]*eb)
            mask = r < self.m.h
            if np.any(mask):
                owner = self.owner[e['geometry']]; hh = np.full(mask.sum(), self.m.H)
                px = basis(owner, 'p', x[mask], hh)[1]
                kx = basis(owner, 'k', x[mask], hh)[1]
                self.C[np.ix_(ui, owner['p'])] += b[mask].T@(w[mask, None]*px)
                self.Ak[np.ix_(ui, owner['k'])] += b[mask].T@(w[mask, None]*kx)
        self.G = -self.D.T; self.Mv = la.block_diag(self.M, self.M)
        self.riesz_u = [la.cho_factor(self.M[np.ix_(ids, ids)], lower=True) for ids in self.ufree]
        self.riesz_p = la.cho_factor(self.Q, lower=True)
        self.riesz_k = la.cho_factor((self.m.Kin**2*self.Mk)[np.ix_(self.kfree, self.kfree)], lower=True)
        self.riesz_e = la.cho_factor(self.Re[np.ix_(self.efree, self.efree)], lower=True)
        self.one = {}
        for key, size in (('u', n), ('p', p), ('k', self.nk), ('epsilon', self.ne)):
            value = np.zeros(size); seen = {}; width = 3 if key == 'p' else 4
            for e in self.elements:
                for j, gid in enumerate(e[key]):
                    expected = 1/e[SIGN[key]][j] if j % width < 2 and j//width < 2 else 0.
                    require(gid not in seen or abs(seen[gid]-expected) < 1e-14, 'constant signed-map conflict')
                    seen[gid] = expected; value[gid] = expected
            self.one[key] = value
        require(la.norm((self.G-self.C)[self.vfree]@self.one['p']) < 1e-11,
                'actual all-essential pressure action lacks the declared gauge')

    def dual(self, load):
        flow = load['flow']; energy = 0.
        for component, ids in enumerate(self.ufree):
            values = flow[component*self.nu+ids]
            energy += float(values@la.cho_solve(self.riesz_u[component], values))
        result = dict(momentum=math.sqrt(max(0., energy)/self.volume))
        for key, ids, factor, area in (('k', self.kfree, self.riesz_k, self.volume),
                                       ('epsilon', self.efree, self.riesz_e, self.evolume)):
            v = load[key][ids]; energy = float(v@la.cho_solve(factor, v))
            result[key] = math.sqrt(max(0., energy)/area)
        return result

    def advection(self, s, when):
        out = self.zero(); advective = self.zero(); awf = np.zeros(2*self.nu)
        for e in self.elements:
            self.clock(); x, r, w = self.rule(e, s, when)
            b, dx, dr = basis(e, 'u', x, r); ids = np.asarray(e['u'])
            u, v = b@s['flow'][ids], b@s['flow'][self.nu+ids]
            if np.any(r < self.m.h):
                owner = self.owner[e['geometry']]; hb, hx, _ = basis(owner, 'u', x, np.full(len(x), self.m.H))
                hu = s['flow'][owner['u']]; conv = np.where(r < self.m.h, 2*(hb@hu)*(hx@hu), 0.)
                awf[ids] += b.T@(w*conv)
            for component in (0, 1):
                coeff = s['flow'][component*self.nu+ids]
                out['flow'][component*self.nu+ids] += b.T@(w*(u*(dx@coeff)+v*(dr@coeff)))
            for key in ('k', 'epsilon'):
                if not e[key]: continue
                qids = np.asarray(e[key]); qb, qx, qr = basis(e, key, x, r); coeff = s[key][qids]; q = qb@coeff
                out[key][qids] -= qx.T@(w*u*q)+qr.T@(w*v*q)
                advective[key][qids] += qb.T@(w*(u*(qx@coeff)+v*(qr@coeff)))
        out['flow'][:2*self.nu] -= awf
        for key in ('k', 'epsilon'):
            for e, side in self.faces[key]:
                x, r, w = self.face(e, side); b = basis(e, 'u', x, r)[0]
                qb = basis(e, key, x, r)[0]; nx, nr = NORMALS[side]
                un = nx*(b@s['flow'][e['u']])+nr*(b@s['flow'][self.nu+np.asarray(e['u'])])
                out[key][e[key]] += qb.T@(w*un*(qb@s[key][e[key]]))
        return out, advective, awf

    def forcing(self, accepted, when, steady=False):
        require(not steady or when == 0., 'steady MMS forcing is frozen at exactly t=0')
        out = self.zero()
        bulk_body = self.zero(); sampled_body = self.zero()
        removed_time = self.zero() if steady else None
        # This is integration by parts of the independently analytic continuous
        # source. It uses no projected target/native matrix and includes exterior
        # terms. Scalar natural terms cancel the matching exact diffusive term
        # on every native callback face; epsilon H retains its Dirichlet reaction.
        for e in self.elements:
            self.clock(); x, r, w = self.rule(e, accepted, when); outer = bool(e['epsilon'])
            law = self.exact_law(x, r, when, outer); j = law['j']
            b, dx, dr = basis(e, 'u', x, r); ids = np.asarray(e['u'])
            atH = self.exact(x, np.full(len(x), self.m.H), when, True)
            awf = np.where(r < self.m.h, 2*atH['u'][0]*atH['u'][1]+atH['Pi'][1]-2*atH['k'][1]/3, 0.)
            if self.m.body_amplitude != 0:
                body = self.m.body(x, r, when)
                body_H = self.m.body(x, np.full(len(x), self.m.H), when)
                sampled_x = np.where(r < self.m.h, body_H[0], 0.)
                awf -= sampled_x
                # Native FullRHS must add b and subtract I_h*b_x(H)*e_x
                # separately from manufactured f_M. These actions are retained
                # individually for missing-term controls; do not simplify them
                # against the body terms in the continuous verification source.
                sampled_body['flow'][ids] -= b.T@(w*sampled_x)
            for component, key, axis in ((0, 'u', 0), (1, 'w', 2)):
                f = j[key][3]+j['u'][0]*j[key][1]+j['w'][0]*j[key][2]+j['Pi'][1 if component == 0 else 2]
                if steady:
                    f -= j[key][3]
                    removed_time['flow'][component*self.nu+ids] += b.T@(w*j[key][3])
                if component == 0: f -= awf
                if self.m.body_amplitude != 0:
                    f -= body[axis]
                    bulk_body['flow'][component*self.nu+ids] += b.T@(w*body[axis])
                out['flow'][component*self.nu+ids] += (b.T@(w*f)+dx.T@(w*law['V'][:, axis, 0])+
                                                       dr.T@(w*law['V'][:, axis, 2]))
            for key in ('k', 'epsilon'):
                if not e[key]: continue
                qb, qx, qr = basis(e, key, x, r); q = j[key]; idsq = np.asarray(e[key])
                D = law['Dk'] if key == 'k' else law['De']
                source = law['production']-law['E'] if key == 'k' else (1.44*law['production']-1.82*q[0])*q[0]/j['k'][0]
                f = q[3]+j['u'][0]*q[1]+j['w'][0]*q[2]+q[0]*(j['u'][1]+j['w'][2])-source
                if steady:
                    f -= q[3]
                    removed_time[key][idsq] += qb.T@(w*q[3])
                out[key][idsq] += qb.T@(w*f)+qx.T@(w*D*q[1])+qr.T@(w*D*q[2])
        for e, side in self.faces['u']:
            x, r, w = self.face(e, side); law = self.exact_law(x, r, when, bool(e['epsilon']))
            b = basis(e, 'u', x, r)[0]; ids = np.asarray(e['u']); nx, nr = NORMALS[side]
            for component, axis in ((0, 0), (1, 2)):
                out['flow'][component*self.nu+ids] -= b.T@(w*(nx*law['V'][:, axis, 0]+nr*law['V'][:, axis, 2]))
        for e, side in self.faces['epsilon']:
            if side != 2: continue
            x, r, w = self.face(e, side); law = self.exact_law(x, r, when, True)
            b = basis(e, 'epsilon', x, r)[0]
            # -int_H phi D grad(e).n with n=-ez; natural H flux is NOT imposed.
            out['epsilon'][e['epsilon']] += b.T@(w*law['De']*law['j']['epsilon'][2])
        self.prescribed_body_actions = dict(bulk=bulk_body, sampled_H=sampled_body)
        if steady:
            self.removed_analytic_time_load = removed_time
        if self.m.body_amplitude != 0:
            out['flow'] += bulk_body['flow']+sampled_body['flow']
        return out

    def current(self, s, when, frozen=None, matrices=False):
        out = self.zero(); Kk = np.zeros((self.nk, self.nk)); Ke = np.zeros((self.ne, self.ne))
        reaction = np.zeros_like(Ke); sources = self.zero(); Ku = np.zeros((2*self.nu, 2*self.nu)) if matrices else None
        for e in self.elements:
            self.clock(); x, r, w = self.rule(e, frozen if frozen is not None else s, when)
            a = self.sampled(e, x, r, s, frozen, matrices); b, dx, dr = basis(e, 'u', x, r); ids = np.asarray(e['u'])
            for component, axis in ((0, 0), (1, 2)):
                out['flow'][component*self.nu+ids] += dx.T@(w*a['V'][:, axis, 0])+dr.T@(w*a['V'][:, axis, 2])
            if matrices:
                gradients = (dx, dr); wm = w*a['mu']; lap = dx.T@(wm[:, None]*dx)+dr.T@(wm[:, None]*dr)
                for i in (0, 1):
                    for j in (0, 1):
                        block = gradients[j].T@(wm[:, None]*gradients[i])
                        if i == j: block += lap
                        Ku[np.ix_(i*self.nu+ids, j*self.nu+ids)] += block
            for key in ('k', 'epsilon'):
                if not e[key]: continue
                qb, qx, qr = basis(e, key, x, r); idsq = np.asarray(e[key]); D = a['Dk'] if key == 'k' else a['De']
                local = qx.T@((w*D)[:, None]*qx)+qr.T@((w*D)[:, None]*qr)
                if key == 'k':
                    Kk[np.ix_(idsq, idsq)] += local
                    sources[key][idsq] += qb.T@(w*(a['production']-a['E']))
                else:
                    Ke[np.ix_(idsq, idsq)] += local
                    reaction[np.ix_(idsq, idsq)] += qb.T@((w*1.82*a['E']/a['k'])[:, None]*qb)
                    sources[key][idsq] += qb.T@(w*1.44*a['production']*a['E']/a['k'])
        out['k'] = Kk@s['k']-sources['k']
        out['epsilon'] = (Ke+reaction)@s['epsilon']-sources['epsilon']
        return dict(action=out, Ku=Ku, Kk=Kk, Ke=Ke+reaction, sources=sources)

    def lift(self, s, when, include_epsilon=True):
        result = self.zero(); seen = {key: {} for key in ('u', 'w', 'k', 'epsilon')}
        for key, scalar_key, sides in (('u', 'u', (0, 1, 2, 3)), ('w', 'u', (0, 1, 2, 3)),
                                       ('k', 'k', (0,)), ('epsilon', 'epsilon', (0, 2))):
            if key == 'epsilon' and not include_epsilon: continue
            def value(e, side, x, r):
                if key == 'epsilon' and side == 2:
                    return self.trace(e, 'k', x, s)**1.5/(self.m.cl*self.m.H)
                f = self.m.functions[key](x, r, when)[0]
                return np.broadcast_to(np.asarray(f), x.shape)
            for e, side in self.faces[scalar_key]:
                if side not in sides: continue
                x, r, w = self.face(e, side); b = basis(e, scalar_key, x, r)[0]
                xl, xr, rl, rr = e['bounds']; modes = EDGES[side]
                xe = np.array([xl, xl] if side == 0 else [xr, xr] if side == 1 else [xl, xr])
                re = np.array([rl, rr] if side < 2 else [rl, rl] if side == 2 else [rr, rr])
                coefficients = value(e, side, xe, re)/np.asarray(e[SIGN[scalar_key]])[list(modes[:2])]
                interior = b[:, modes[2:]]
                remaining = value(e, side, x, r)-b[:, modes[0]]*coefficients[0]-b[:, modes[1]]*coefficients[1]
                coefficients = np.r_[coefficients, la.solve(interior.T@(w[:, None]*interior), interior.T@(w*remaining), assume_a='pos')]
                for local, val in zip(modes, coefficients):
                    gid = e[scalar_key][local]; old = seen[key].get(gid)
                    scale = 1. if key in ('u', 'w') else self.m.Kin if key == 'k' else self.m.Ain/self.m.H
                    require(old is None or abs(old-val) <= 1e-11*scale,
                            'independent essential corner/edge is incompatible')
                    if old is None: seen[key][gid] = float(val)
            for gid, val in seen[key].items():
                if key in ('u', 'w'): result['flow'][(key == 'w')*self.nu+gid] = val
                else: result[key][gid] = val
        return result

    def impose(self, s, when):
        result = {key: value.copy() for key, value in s.items()}
        # Update k's inlet first: Gamma must see that NEW corner value before
        # it is compared with the outer epsilon inlet at the new time.
        lift = self.lift(result, when, include_epsilon=False)
        for component in (0, 1):
            ids = component*self.nu+self.essential['u']; result['flow'][ids] = lift['flow'][ids]
        result['k'][self.essential['k']] = lift['k'][self.essential['k']]
        lift = self.lift(result, when)
        result['epsilon'][self.essential['epsilon']] = lift['epsilon'][self.essential['epsilon']]
        return result

    def trace_error(self, s, when):
        lift = self.lift(s, when); errors = {}
        for key in ('k', 'epsilon'):
            ids = self.essential[key]; scale = self.m.Kin if key == 'k' else self.m.Ain/self.m.H
            errors[key] = float(np.max(abs(s[key][ids]-lift[key][ids]), initial=0.)/scale)
        ids = np.r_[self.essential['u'], self.nu+self.essential['u']]
        errors['velocity'] = float(np.max(abs(s['flow'][ids]-lift['flow'][ids]), initial=0.))
        return errors

    def mass(self, s):
        result = self.zero(); result['flow'][:2*self.nu] = self.Mv@s['flow'][:2*self.nu]
        result['k'] = self.Mk@s['k']; result['epsilon'] = self.Me@s['epsilon']; return result

    def residual(self, s, when, alpha, base, assembled=None):
        terms = assembled if assembled is not None else self.current(s, when)
        out = ({key: value.copy() for key, value in terms['action'].items()} if alpha == 0.
               else combine(self.mass(s), alpha, terms['action'], 1.))
        for key in out: out[key] -= base[key]
        pressure = s['flow'][2*self.nu:]
        out['flow'][:2*self.nu] += (self.G-self.C)@pressure+(2./3)*self.Ak@s['k']
        out['flow'][2*self.nu:] = self.D@s['flow'][:2*self.nu]
        norm = self.dual(out); cont = out['flow'][2*self.nu:]
        norm['continuity'] = math.sqrt(max(0., float(cont@la.cho_solve(self.riesz_p, cont)))/self.volume)
        norm['gauge'] = abs(float(self.mean@pressure))/(.01*self.volume)
        return out, norm

    def direct(self, seed, accepted, when, alpha, base, seconds, steady=False):
        require(not steady or (alpha == 0. and when == 0.), 'steady direct reference forbids time-mass or advancing clock')
        stop = min(self.deadline, time.monotonic()+seconds); rows = []; s = self.impose(seed, when)
        self.direct_trace = rows
        for iteration in range(40):
            if time.monotonic() >= stop:
                return dict(status='NOT_MET_DIRECT_TIME_CAP', iterations=rows)
            a = self.current(s, when, matrices=True); nv = len(self.vfree); n = self.nu; p = self.np
            active_base = combine(base, 1., self.advection(s, when)[0], -1.) if steady else base
            velocity_matrix = a['Ku'] if steady else alpha*self.Mv+a['Ku']; gp = self.G-self.C
            system = np.zeros((nv+p+1, nv+p+1)); system[:nv, :nv] = velocity_matrix[np.ix_(self.vfree, self.vfree)]
            system[:nv, nv:nv+p] = gp[self.vfree]; system[nv:nv+p, :nv] = self.D[:, self.vfree]
            system[nv:nv+p, -1] = self.mean; system[-1, nv:nv+p] = self.mean
            lift = self.lift(s, when)['flow'][:2*n]
            remainder_rhs = a['Ku']@s['flow'][:2*n]-a['action']['flow'][:2*n]
            rhs = np.r_[(active_base['flow'][:2*n]+remainder_rhs-(2./3)*self.Ak@s['k']-velocity_matrix@lift)[self.vfree],
                        -self.D@lift, 0.]
            solution = la.solve(system, rhs, assume_a='gen')
            s['flow'][:2*n] = lift; s['flow'][self.vfree] = solution[:nv]; s['flow'][2*n:] = solution[nv:nv+p]
            a = self.current(s, when); A = a['Kk'] if steady else alpha*self.Mk+a['Kk']; lift = self.lift(s, when)['k']
            active_base = combine(base, 1., self.advection(s, when)[0], -1.) if steady else base
            s['k'] = lift; s['k'][self.kfree] = la.solve(A[np.ix_(self.kfree, self.kfree)],
                (active_base['k']+a['sources']['k']-A@lift)[self.kfree], assume_a='pos')
            s = self.impose(s, when); a = self.current(s, when); A = a['Ke'] if steady else alpha*self.Me+a['Ke']; lift = self.lift(s, when)['epsilon']
            active_base = combine(base, 1., self.advection(s, when)[0], -1.) if steady else base
            s['epsilon'] = lift; s['epsilon'][self.efree] = la.solve(A[np.ix_(self.efree, self.efree)],
                (active_base['epsilon']+a['sources']['epsilon']-A@lift)[self.efree], assume_a='pos')
            active_base = combine(base, 1., self.advection(s, when)[0], -1.) if steady else base
            _, norm = self.residual(s, when, alpha, active_base); rows.append(norm)
            if max(norm.values()) <= 1e-10:
                d = {key: s[key]-accepted[key] for key in s}
                difference = dict(U=math.sqrt(float(d['flow'][:2*n]@self.Mv@d['flow'][:2*n])/self.volume),
                    Pi=math.sqrt(float(d['flow'][2*n:]@self.Q@d['flow'][2*n:])/self.volume)/self.m.tau0,
                    k=math.sqrt(float(d['k']@self.Mk@d['k'])/self.volume)/self.m.Kin,
                    epsilon=math.sqrt(float(d['epsilon']@self.Eerror@d['epsilon'])/self.evolume))
                return dict(status='PASS_INDEPENDENT_DENSE_NONLINEAR_REFERENCE' if max(difference.values()) <= 1e-8 else 'NOT_MET_DIRECT_FIELD_AGREEMENT',
                    iterations=rows, reference_scaled_RMS_difference=difference, final_state_sha256=state_hash(s),
                    coefficients={key: value.tolist() for key, value in s.items()},
                    method=('Independent dense block-direct steady nonlinear iteration from original fresh seed; current advection refreshed before every block and final residual, alpha=0, t=0.' if steady else
                            'Independent dense block-direct current nonlinear iteration from pre-step state; same spaces with independently resolved quadrature.'))
        return dict(status='NOT_MET_DIRECT_ITERATION_CAP', iterations=rows)


def audit(data, table, args, report, deadline):
    require('mms_variant' in data, 'current coupled oracle requires explicit state.mms_variant')
    m = coupled_mms_check.IndependentMMS(table, data.get('mms_body_amplitude', 0.), data['mms_variant'])
    report['mms_variant'] = m.variant
    report['mms_body_amplitude'] = m.body_amplitude
    report['body_amplitude_provenance'] = ('state.mms_body_amplitude' if 'mms_body_amplitude' in data
                                          else 'Historical absent-field default is exactly zero.')
    sizes = data['sizes']; s = state(data['coefficients'], sizes)
    h = data['input_history']; current = state(h['current'], sizes); previous = state(h['previous'], sizes)
    advection = state(h['advection'], sizes); old_advection = state(h['previous_advection'], sizes)
    require(data['step'] == h['step']+1 and data['dt'] == h['dt'] > 0 and
            abs(data['time']-data['step']*data['dt']) < 1e-13 and abs(h['time']-h['step']*h['dt']) < 1e-13,
            'accepted/input step time or unchanged dt differs')
    require(state_hash(s) == data['state_sha256'] and state_hash(current) == data['previous_sha256'] and
            state_hash(advection) == data['previous_advection_sha256'], 'accepted state/history byte identities differ')
    startup = h['step'] == 0; dt = h['dt']; alpha = (1. if startup else 1.5)/dt
    history = combine(current, 1/dt if startup else 2/dt, previous, 0. if startup else -.5/dt)
    report['history_scheme'] = 'BE_STARTUP_EX2_NOT_EXERCISED' if startup else 'BDF2_EX2'
    report['ex2_exercised'] = not startup
    comparisons = []; oracles = []; residuals = []; bases = []; input_adv = []
    for nx, nr in ((args.nx, args.nr), (args.refine_nx, args.refine_nr)):
        o = Oracle(data, m, nx, nr, deadline); oracles.append(o)
        new_adv, conv_new, rough_new = o.advection(current, h['time'])
        prev_adv, conv_old, rough_old = o.advection(previous, max(0, h['step']-1)*dt)
        errors = {label: o.dual(combine(actual, 1., expected, -1.)) for label, actual, expected in
                  (('current', advection, new_adv), ('previous', old_advection, prev_adv))}
        require(all(np.max(abs(v['flow'][2*o.nu:]), initial=0.) == 0 for v in (advection, old_advection)),
                'stored explicit history has a pressure load')
        explicit = combine(new_adv, 1. if startup else 2., prev_adv, 0. if startup else -1.)
        base = combine(o.forcing(s, data['time']), 1., o.mass(history), 1.)
        base = combine(base, 1., explicit, -1.); bases.append(base)
        assembled = o.current(s, data['time']); residual, norms = o.residual(s, data['time'], alpha, base, assembled)
        residuals.append(residual); input_adv.append((new_adv, prev_adv, conv_new, conv_old, rough_new, rough_old))
        comparisons.append(dict(quadrature=[nx, nr], current_residual=norms, stored_history_errors=errors,
            current_trace=o.trace_error(s, data['time']), input_trace=o.trace_error(current, h['time']),
            previous_trace=o.trace_error(previous, max(0, h['step']-1)*dt)))
    report['independent_assemblies'] = comparisons
    o = oracles[-1]; base = bases[-1]; norms = comparisons[-1]['current_residual']
    difference = combine(residuals[1], 1., residuals[0], -1.); qerror = o.dual(difference)
    cont = difference['flow'][2*o.nu:]
    qerror['continuity'] = math.sqrt(max(0., float(cont@la.cho_solve(o.riesz_p, cont)))/o.volume)
    qerror['gauge'] = abs(float((oracles[1].mean-oracles[0].mean)@s['flow'][2*o.nu:]))/(.01*o.volume)
    report['quadrature_difference_dual'] = qerror
    history_ok = all(value <= 1e-10 for row in comparisons[-1]['stored_history_errors'].values() for value in row.values())
    traces_ok = max(value for name in ('current_trace', 'input_trace', 'previous_trace')
                    for value in comparisons[-1][name].values()) <= 1e-10
    residual_ok = all(norms[key] <= threshold for key, threshold in LIMITS.items())
    quadrature_ok = all(qerror[key] <= .2*threshold for key, threshold in LIMITS.items())
    report.update(current_residual_pass=residual_ok, stored_history_pass=history_ok,
                  essential_trace_pass=traces_ok, independent_quadrature_pass=quadrature_ok)
    new_adv, prev_adv, conv_new, conv_old, rough_new, rough_old = input_adv[-1]
    controls = {}
    def control(name, change):
        values = o.dual(change)
        controls[name] = dict(change_dual=values,
            status='DETECTED' if any(values[k] > 10*LIMITS[k] for k in values) else 'NOT_EXERCISED_AT_THIS_STATE')
    a = o.current(s, data['time']); frozen = o.current(s, data['time'], frozen=current)
    control('frozen_prior_scalar_and_source', combine(frozen['action'], 1., a['action'], -1.))
    v = o.zero(); projection = la.solve(o.Q, o.cross@s['k'], assume_a='pos')
    v['flow'][:2*o.nu] = (2./3)*(o.C@projection-o.Ak@s['k']); control('replace_actual_Q3_k_by_Q2_projection', v)
    v = o.zero(); v['flow'][:2*o.nu] = rough_new if startup else 2*rough_new-rough_old; control('missing_sampled_AWF_convection', v)
    v = o.zero(); v['flow'][:2*o.nu] = o.C@s['flow'][2*o.nu:]; control('missing_sampled_AWF_pressure', v)
    v = o.zero(); v['flow'][:2*o.nu] = -(2./3)*o.Ak@s['k']; control('missing_sampled_AWF_Q3_k', v)
    # Removing a signed RHS action changes LHS-RHS by that same signed action.
    # Bulk b and -I_h*b_x(H) are deliberately tested independently, although
    # their sum cancels the corresponding change in the manufactured source.
    control('missing_bulk_physical_body', o.prescribed_body_actions['bulk'])
    control('missing_sampled_AWF_body', o.prescribed_body_actions['sampled_H'])
    v = o.zero()
    for key in ('k', 'epsilon'):
        v[key] = conv_new[key]-new_adv[key] if startup else 2*(conv_new[key]-new_adv[key])-(conv_old[key]-prev_adv[key])
    control('replace_conservative_scalar_transport_by_advective_form', v)
    control('stale_EX1_instead_of_EX2', o.zero() if startup else combine(prev_adv, 1., new_adv, -1.))
    report['negative_controls'] = controls
    body_coverage_ok = m.body_amplitude == 0 or all(controls[name]['status'] == 'DETECTED' for name in
        ('missing_bulk_physical_body', 'missing_sampled_AWF_body'))
    report['nonzero_prescribed_body_mutation_coverage_pass'] = m.body_amplitude != 0 and body_coverage_ok
    report['prescribed_body_actions'] = dict(
        bulk_rhs=o.prescribed_body_actions['bulk']['flow'][:2*o.nu].tolist(),
        sampled_H_rhs=o.prescribed_body_actions['sampled_H']['flow'][:2*o.nu].tolist(),
        signs='FullRHS adds integral phi*b and subtracts integral_(r<h) phi*b_x(x,H,t)*e_x, distinct from f_M.',
        source='Independent analytic prescribed body, not current state or native cancellation data.')
    native_body = {}
    for field, action in (('mms_bulk_body_rhs', 'bulk'), ('mms_sampled_body_rhs', 'sampled_H')):
        if field not in data:
            continue
        values = np.asarray(data[field], dtype=float)
        require(values.shape == (2*o.nu,) and np.all(np.isfinite(values)),
                'native separate body action has wrong canonical shape or nonfinite entries')
        delta = o.zero(); delta['flow'][:2*o.nu] = values-o.prescribed_body_actions[action]['flow'][:2*o.nu]
        error = o.dual(delta)['momentum']
        native_body[field] = dict(momentum_dual_error=error, passed=error <= 1e-10)
    actions_present = len(native_body) == 2
    actions_pass = actions_present and all(value['passed'] for value in native_body.values())
    report['native_prescribed_body_action_comparison'] = dict(
        status=('PASS' if actions_pass else 'NOT_MET' if actions_present else 'NOT_PROVIDED'),
        actions=native_body, momentum_dual_tolerance=1e-10,
        limitation='Net manufactured forcing cancels body changes; the final state alone cannot establish that the native separate body paths ran.')
    # Historical zero-body data need no new ABI. Missing separate actions leave
    # body-path closure explicitly open; provided wrong actions fail the step.
    body_coverage_ok &= all(value['passed'] for value in native_body.values())
    report['nonzero_prescribed_body_coverage_pass'] = m.body_amplitude != 0 and body_coverage_ok and actions_pass
    if 'current_advection' in data:
        actual = state(data['current_advection'], sizes); expected = o.advection(s, data['time'])[0]
        require(state_hash(actual) == data['advection_sha256'], 'accepted advection hash differs')
        report['accepted_advection_error'] = o.dual(combine(actual, 1., expected, -1.))
        history_ok &= max(report['accepted_advection_error'].values()) <= 1e-10
        report['stored_history_pass'] = history_ok
    if 'previous_advection' in data:
        require(state_hash(state(data['previous_advection'], sizes)) == data['previous_advection_sha256'],
                'reported previous advection array/hash differs')
    report['essential_reaction_diagnostics'] = {
        key: dict(constant_test_full_residual=float(o.one[key]@residuals[-1][key]),
                  essential_coefficient_reaction_norm=float(la.norm(residuals[-1][key][o.essential[key]])),
                  note='Essential reactions are reported, not constrained to zero.') for key in ('k', 'epsilon')}
    report['direct_nonlinear'] = dict(status='NOT_PERFORMED', reason='Optional separate --direct admission required; no direct-reference gate is claimed.')
    if args.direct:
        original_deadline = o.deadline
        o.deadline = min(original_deadline, time.monotonic()+args.direct_seconds)
        try:
            report['direct_nonlinear'] = o.direct(current, s, data['time'], alpha, base, args.direct_seconds)
        except Exception as error:
            report['direct_nonlinear'] = dict(status='NOT_MET_DIRECT_REFERENCE',
                exception=type(error).__name__, reason=str(error), iterations=getattr(o, 'direct_trace', []))
        finally:
            o.deadline = original_deadline
    passed = residual_ok and history_ok and traces_ok and quadrature_ok and body_coverage_ok
    report['current_coupled_step_gate_pass'] = passed
    report['status'] = 'PASS_CURRENT_COUPLED_STEP_ONLY' if passed else 'NOT_MET_CURRENT_COUPLED_STEP'
    report['direct_reference_gate_pass'] = report['direct_nonlinear']['status'] == 'PASS_INDEPENDENT_DENSE_NONLINEAR_REFERENCE'
    report['norm_limits'] = LIMITS
    report['quadrature_reserve_fraction'] = .2
    report['native_reported_residuals'] = data['residuals']
    report['norms'] = 'Free-test mass dual / sqrt(domain volume): U force and continuity scale1; k uses integral Kin² phi phi; epsilon uses integral (Ain/r)² phi phi on outer volume; gauge |integral Pi|/(.01 volume). No alpha normalization.'


def audit_steady(data, table, args, report, deadline):
    """Current true steady MMS equations, distinct from the time-step audit."""
    require('mms_variant' in data, 'steady coupled oracle requires explicit state.mms_variant')
    m = coupled_mms_check.IndependentMMS(table, data.get('mms_body_amplitude', 0.), data['mms_variant'])
    sizes = data['sizes']; s = state(data['coefficients'], sizes)
    h = data['input_history']; seed = state(h['current'], sizes)
    require(data['time'] == h['time'] == 0. and data['step'] == h['step'] == 0,
            'steady evidence must retain the fresh zero-time/zero-step clock')
    require(state_hash(s) == data['state_sha256'] == data['previous_sha256'],
            'steady current/previous state hashes must both bind the solved state')
    require('current_advection' in data and 'previous_advection' in data,
            'steady evidence needs actual solved-state current/previous advection arrays')
    cached = state(data['current_advection'], sizes)
    previous_cached = state(data['previous_advection'], sizes)
    require(state_hash(cached) == data['advection_sha256'] == data['previous_advection_sha256'] and
            state_hash(previous_cached) == data['previous_advection_sha256'],
            'steady actual current advection identity differs')
    report.update(mms_variant=m.variant, mms_body_amplitude=m.body_amplitude,
        steady_mode=True, physical_time_advanced=0., alpha=0., ex2_exercised=False,
        temporal_order_qualified=False, physical_reference_gate_pass=False,
        history_scheme='STEADY_CURRENT_ADVECTION_NO_TIME_HISTORY',
        steady_seed_state_sha256=state_hash(seed),
        seed_binding='input_history.current is the original fresh seed; previous_sha256 binds the solved state, not this seed.',
        forcing_route='Continuous analytic t=0 MMS forcing integrated by parts, with analytic u_t,w_t,k_t,epsilon_t removed. No time-mass, EX2 load or projected A_h*target source.',
        scope='Current true steady MMS residual and optional independent dense steady reference only. No advancing-time, temporal-order, physical steady-reference or full H27 claim.')
    comparisons = []; residuals = []; oracles = []; forcings = []; advections = []
    for nx, nr in ((args.nx, args.nr), (args.refine_nx, args.refine_nr)):
        o = Oracle(data, m, nx, nr, deadline); oracles.append(o)
        forcing = o.forcing(s, 0., steady=True); forcings.append(forcing)
        current_adv, advective, rough = o.advection(s, 0.); advections.append((current_adv, advective, rough))
        require(np.max(abs(cached['flow'][2*o.nu:]), initial=0.) == 0 and
                np.max(abs(previous_cached['flow'][2*o.nu:]), initial=0.) == 0,
                'steady advection cache contains a pressure load')
        base = combine(forcing, 1., current_adv, -1.)
        residual, norms = o.residual(s, 0., 0., base); residuals.append(residual)
        comparisons.append(dict(quadrature=[nx, nr], current_residual=norms,
            stored_current_advection_error=o.dual(combine(cached, 1., current_adv, -1.)),
            current_trace=o.trace_error(s, 0.), seed_trace=o.trace_error(seed, 0.)))
    report['independent_assemblies'] = comparisons
    o = oracles[-1]; forcing = forcings[-1]; norms = comparisons[-1]['current_residual']
    difference = combine(residuals[1], 1., residuals[0], -1.); qerror = o.dual(difference)
    cont = difference['flow'][2*o.nu:]
    qerror['continuity'] = math.sqrt(max(0., float(cont@la.cho_solve(o.riesz_p, cont)))/o.volume)
    qerror['gauge'] = abs(float((oracles[1].mean-oracles[0].mean)@s['flow'][2*o.nu:]))/(.01*o.volume)
    residual_ok = all(norms[key] <= limit for key, limit in LIMITS.items())
    quadrature_ok = all(qerror[key] <= .2*limit for key, limit in LIMITS.items())
    advection_ok = max(comparisons[-1]['stored_current_advection_error'].values()) <= 1e-10
    traces_ok = max(value for label in ('current_trace', 'seed_trace')
                    for value in comparisons[-1][label].values()) <= 1e-10
    report.update(quadrature_difference_dual=qerror, current_residual_pass=residual_ok,
        independent_quadrature_pass=quadrature_ok, stored_current_advection_pass=advection_ok,
        essential_trace_pass=traces_ok,
        continuity_check='Physical D*U and pressure mean are checked directly; no gauge multiplier can absorb a continuity defect.')
    current_adv, advective, rough = advections[-1]
    controls = {}
    def control(name, change):
        values = o.dual(change)
        controls[name] = dict(change_dual=values,
            status='DETECTED' if any(values[key] > 10*LIMITS[key] for key in values) else 'NOT_EXERCISED_AT_THIS_STATE')
    a = o.current(s, 0.); frozen = o.current(s, 0., frozen=seed)
    control('frozen_seed_scalar_and_source', combine(frozen['action'], 1., a['action'], -1.))
    seed_adv = o.advection(seed, 0.)[0]
    control('stale_seed_instead_of_current_advection', combine(seed_adv, 1., current_adv, -1.))
    control('missing_current_advection', combine(current_adv, -1., o.zero(), 0.))
    control('retaining_unsteady_MMS_time_derivatives', combine(o.removed_analytic_time_load, -1., o.zero(), 0.))
    v = o.zero(); projection = la.solve(o.Q, o.cross@s['k'], assume_a='pos')
    v['flow'][:2*o.nu] = (2./3)*(o.C@projection-o.Ak@s['k']); control('replace_actual_Q3_k_by_Q2_projection', v)
    v = o.zero(); v['flow'][:2*o.nu] = rough; control('missing_sampled_AWF_convection', v)
    v = o.zero(); v['flow'][:2*o.nu] = o.C@s['flow'][2*o.nu:]; control('missing_sampled_AWF_pressure', v)
    v = o.zero(); v['flow'][:2*o.nu] = -(2./3)*o.Ak@s['k']; control('missing_sampled_AWF_Q3_k', v)
    v = o.zero()
    for key in ('k', 'epsilon'): v[key] = advective[key]-current_adv[key]
    control('replace_conservative_scalar_transport_by_advective_form', v)
    control('missing_bulk_physical_body', o.prescribed_body_actions['bulk'])
    control('missing_sampled_AWF_body', o.prescribed_body_actions['sampled_H'])
    report['negative_controls'] = controls
    body_ok = m.body_amplitude == 0 or all(controls[key]['status'] == 'DETECTED' for key in
        ('missing_bulk_physical_body', 'missing_sampled_AWF_body'))
    report['nonzero_prescribed_body_mutation_coverage_pass'] = m.body_amplitude != 0 and body_ok
    native_body = {}
    for field, action in (('mms_bulk_body_rhs', 'bulk'), ('mms_sampled_body_rhs', 'sampled_H')):
        if field not in data: continue
        values = np.asarray(data[field], dtype=float)
        require(values.shape == (2*o.nu,) and np.all(np.isfinite(values)), 'invalid native steady body action')
        delta = o.zero(); delta['flow'][:2*o.nu] = values-o.prescribed_body_actions[action]['flow'][:2*o.nu]
        error = o.dual(delta)['momentum']
        native_body[field] = dict(momentum_dual_error=error, passed=error <= 1e-10)
    actions_present = len(native_body) == 2
    actions_pass = actions_present and all(value['passed'] for value in native_body.values())
    body_ok &= all(value['passed'] for value in native_body.values())
    report['native_prescribed_body_action_comparison'] = dict(actions=native_body,
        status='PASS' if actions_pass else 'NOT_MET' if actions_present else 'NOT_PROVIDED',
        momentum_dual_tolerance=1e-10,
        limitation='The net manufactured source cancels body additions; a solved state alone cannot verify separate body paths.')
    report['nonzero_prescribed_body_coverage_pass'] = m.body_amplitude != 0 and body_ok and actions_pass
    report['essential_reaction_diagnostics'] = {
        key: dict(constant_test_full_residual=float(o.one[key]@residuals[-1][key]),
                  essential_coefficient_reaction_norm=float(la.norm(residuals[-1][key][o.essential[key]])),
                  note='Essential reactions are reported, not constrained to zero.') for key in ('k', 'epsilon')}
    report['direct_nonlinear'] = dict(status='NOT_PERFORMED',
        reason='Optional --direct is required for independent dense steady reference evidence.')
    if args.direct:
        original_deadline = o.deadline
        o.deadline = min(original_deadline, time.monotonic()+args.direct_seconds)
        try:
            report['direct_nonlinear'] = o.direct(seed, s, 0., 0., forcing, args.direct_seconds, steady=True)
        except Exception as error:
            report['direct_nonlinear'] = dict(status='NOT_MET_DIRECT_REFERENCE',
                exception=type(error).__name__, reason=str(error), iterations=getattr(o, 'direct_trace', []))
        finally:
            o.deadline = original_deadline
    passed = residual_ok and quadrature_ok and advection_ok and traces_ok and body_ok
    report.update(current_coupled_step_gate_pass=False, current_coupled_steady_gate_pass=passed,
        status='PASS_CURRENT_COUPLED_STEADY_ONLY' if passed else 'NOT_MET_CURRENT_COUPLED_STEADY',
        direct_reference_gate_pass=report['direct_nonlinear']['status'] == 'PASS_INDEPENDENT_DENSE_NONLINEAR_REFERENCE',
        norm_limits=LIMITS, quadrature_reserve_fraction=.2, native_reported_residuals=data['residuals'],
        norms='Same free-test reference-weighted mass duals and volume factors as the unsteady audit; alpha=0 and no mass normalization.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--nx', type=int, default=20); parser.add_argument('--nr', type=int, default=24)
    parser.add_argument('--refine-nx', type=int, default=24); parser.add_argument('--refine-nr', type=int, default=32)
    parser.add_argument('--seconds', type=float, default=240.)
    parser.add_argument('--direct', action='store_true'); parser.add_argument('--direct-seconds', type=float, default=150.)
    args = parser.parse_args()
    if args.out.exists(): parser.error('immutable oracle output already exists')
    started = time.monotonic()
    report = dict(schema='WJH27IndependentCoupledStepOracleV1', status='FAIL_COUPLED_STEP_ORACLE',
        full_coupled_gate_pass=False, current_coupled_step_gate_pass=False,
        current_coupled_steady_gate_pass=False, direct_reference_gate_pass=False,
        platform=sys.platform, native=str(args.native), table=str(args.table),
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        continuous_MMS_source_sha256=hashlib.sha256(Path(coupled_mms_check.__file__).read_bytes()).hexdigest(),
        independent_WJ_source_sha256=hashlib.sha256(Path(reference_oracle.__file__).read_bytes()).hexdigest(),
        independent_basis_sha256=hashlib.sha256((Path(__file__).resolve().parents[1]/'h26'/'fixtures.py').read_bytes()).hexdigest(),
        scope='Accepted current coupled MMS step only; startup is identified explicitly. No temporal-order, physical-development, principal-certification or Task3 claim.',
        forcing_route='Continuous analytic fields integrated by parts with complete exterior momentum stress terms and exact scalar natural-flux cancellation; no native/projected A_h*target forcing.')
    failure = None
    try:
        require(sys.platform.startswith('linux'), 'authoritative numerical audit requires Linux')
        require(args.nx >= 20 and args.nr >= 24 and args.refine_nx > args.nx and args.refine_nr > args.nr,
                'independent quadrature/refinement is insufficient')
        require(0 < args.seconds <= 300 and 0 < args.direct_seconds <= args.seconds,
                'oracle elapsed budget exceeds bounded admission')
        packet = json.loads(args.native.read_text()); table = json.loads(args.table.read_text())
        require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True,
                'input is not a registered normal-main coupled run')
        steady_mode = packet.get('steady_mode', False)
        require(type(steady_mode) is bool, 'top-level steady_mode must be an explicit boolean')
        report['steady_mode'] = steady_mode
        if steady_mode:
            require(packet.get('physical_time_advanced', 0.) == 0., 'steady artifact reports physical-time advancement')
            report['schema'] = 'WJH27IndependentCoupledSteadyOracleV1'
            report['status'] = 'FAIL_COUPLED_STEADY_ORACLE'
            report['scope'] = 'Current true steady MMS audit only; no EX2, advancing-time, temporal-order, physical-reference or full H27 claim.'
        data = packet['state']; require(data['schema'] == 'WJH27NativeCoupledStateV2', 'wrong state ABI')
        require(table['status'] == 'PASS' and data['reference_identity'] == hashlib.sha256(args.table.read_bytes()).hexdigest(),
                'qualified reference identity differs')
        require(data['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and len(data['source_identity']) == 64,
                'unbound or wrong model/source identity')
        configuration = ({key: data[key] for key in ('quadratureX', 'quadratureR',
                          'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations')}
                         if all(key in data for key in ('quadratureX', 'quadratureR',
                            'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations'))
                         else dict(quadratureX=8, quadratureR=12, linear_tolerance=1e-12,
                            nonlinear_tolerance=1e-10, maximum_iterations=60,
                            provenance='Coordinator-declared source-v3 smoke defaults; not present in this report.'))
        report.update(native_sha256=hashlib.sha256(args.native.read_bytes()).hexdigest(),
                      table_sha256=hashlib.sha256(args.table.read_bytes()).hexdigest(),
                      native_source_identity=data['source_identity'], native_signature=data['signature'],
                      native_configuration=configuration)
        configuration['mms_body_amplitude'] = data.get('mms_body_amplitude', 0.)
        configuration['mms_variant'] = data.get('mms_variant')
        if steady_mode: audit_steady(data, table, args, report, started+args.seconds)
        else: audit(data, table, args, report, started+args.seconds)
    except Exception as error:
        failure = error; report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream: stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if failure is not None: raise failure
    print(json.dumps(dict(status=report['status'], current_coupled_step_gate_pass=report['current_coupled_step_gate_pass'],
                          current_coupled_steady_gate_pass=report['current_coupled_steady_gate_pass'],
                          direct_reference_gate_pass=report['direct_reference_gate_pass'], elapsed_seconds=report['elapsed_seconds']), sort_keys=True))
    current_pass = (report['current_coupled_steady_gate_pass'] if report.get('steady_mode', False)
                    else report['current_coupled_step_gate_pass'])
    if not current_pass or (args.direct and not report['direct_reference_gate_pass']):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
