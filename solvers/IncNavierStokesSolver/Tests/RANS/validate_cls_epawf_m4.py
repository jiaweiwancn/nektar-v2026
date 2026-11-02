#!/usr/bin/env python3
"""Linux qualification of the flat composite reference and native IMEX map.

Independent Jacobi/physical quadrature checks use native maps as structural
metadata. The constitutive oracle is the dimensional indexed M1 implementation;
no exported native stress, production, quadrature or mass matrix is an oracle.
Raw sessions, logs, states and complete checkpoints remain in the run directory.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
sys.dont_write_bytecode = True
THREADS = {k: '1' for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import scipy
from validate_cls_epawf_m1 import indexed
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m3 import basis, wall_integrals, manufactured_values, residual_oracle
from cls_epawf_m4_cases import neutral_session, evolving_session, tensor_mesh, write, set_parameter, set_info

checks, errors = 0, {}

def require(value, message):
    global checks
    checks += 1
    if not value:
        raise RuntimeError(message)

def near(actual, expected, name, scale=1., tolerance=1e-10):
    a, b = np.asarray(actual), np.asarray(expected)
    require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'Invalid '+name)
    error = float(np.max(np.abs(a-b), initial=0.))/max(scale, float(np.max(np.abs(b), initial=0.)))
    errors[name] = max(error, errors.get(name, 0.))
    require(error <= tolerance, f'{name}: {error} > {tolerance}')
    return error

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class Space:
    def __init__(self, record, dimension):
        self.elements, self.size, self.dimension = record['mesh'], record['size'], dimension
        self.axes = (0, 2) if dimension == 2 else (0, 1, 2)
        self.by_id = {e['geometry']: e for e in self.elements}
        self.lower = np.min([e['lower'] for e in self.elements], axis=0)
        self.upper = np.max([e['upper'] for e in self.elements], axis=0)
        self._mass = None

    def element(self, point, side=0):
        possible = []
        p = np.asarray(point)
        for e in self.elements:
            if np.all(p >= np.array(e['lower'])-1e-11) and np.all(p <= np.array(e['upper'])+1e-11):
                if side > 0 and p[2] >= e['upper'][2]-1e-12 or side < 0 and p[2] <= e['lower'][2]+1e-12:
                    continue
                possible.append(e)
        require(bool(possible), 'Independent sampler missed a physical element')
        return min(possible, key=lambda e: e['geometry'])

    def rows(self, e, points):
        points = np.atleast_2d(points)
        b, db = [], []
        for axis in self.axes:
            length = e['upper'][axis]-e['lower'][axis]
            p, d = basis(e['degree'], 2*(points[:, axis]-e['lower'][axis])/length-1)
            b.append(p); db.append(2*d/length)
        def tensor(factors):
            product = factors[-1]
            for factor in factors[-2::-1]:
                product = (product[:, :, None]*factor[:, None, :]).reshape(len(points), -1)
            return product*np.asarray(e['sign'])
        values = tensor(b)
        gradients = np.zeros((len(points), 3, values.shape[1]))
        for i, axis in enumerate(self.axes):
            factors = list(b); factors[i] = db[i]
            gradients[:, axis] = tensor(factors)
        return np.asarray(e['ids'], dtype=int), values, gradients

    def evaluate(self, state, offset, point, side=0):
        e = self.element(point, side)
        ids, p, g = self.rows(e, [point]); c = np.asarray(state)[offset+ids]
        return float((p@c)[0]), g[0]@c

    def quadrature(self, e, order=None):
        z, w = np.polynomial.legendre.leggauss(order or e['degree']+2)
        points, weights = [], []
        for index in itertools.product(range(len(z)), repeat=self.dimension):
            p, W = [0., 0., 0.], 1.
            for axis, i in zip(self.axes, index):
                width = e['upper'][axis]-e['lower'][axis]
                p[axis] = e['lower'][axis]+(z[i]+1)*width/2
                W *= w[i]*width/2
            points.append(p); weights.append(W)
        return np.asarray(points), np.asarray(weights)

    @property
    def volume(self):
        return sum(math.prod(e['upper'][a]-e['lower'][a] for a in self.axes) for e in self.elements)

    def mass(self):
        if self._mass is None:
            self._mass = np.zeros((self.size, self.size))
            for e in self.elements:
                x, w = self.quadrature(e)
                ids, p, _ = self.rows(e, x)
                self._mass[np.ix_(ids, ids)] += p.T@(w[:, None]*p)
        return self._mass

    def constant(self):
        c = np.zeros(self.size)
        for e in self.elements:
            for local, modes in enumerate(itertools.product(range(e['degree']+1), repeat=self.dimension)):
                # Native tensor index has physical x as its fastest coordinate.
                if all(j < 2 for j in modes):
                    c[e['ids'][local]] = e['sign'][local]
        return c

    def boundary(self, predicate):
        ids = set()
        for e in self.elements:
            n = e['degree']+1
            for axis_id, axis in enumerate(self.axes):
                for side in (0, 1):
                    if predicate(axis, side, e):
                        for local in range(n**self.dimension):
                            if (local//(n**axis_id)) % n == side:
                                ids.add(e['ids'][local])
        return ids


def spaces(raw):
    return [Space(record, raw['dimension']) for record in raw['fields']]

def constrained(raw, field, block):
    d, H = raw['dimension'], raw['parameters'][1]
    if block >= d+2:
        return set()
    def face(axis, side, e):
        if axis != 2:
            return abs((e['upper'] if side else e['lower'])[axis]-(field.upper if side else field.lower)[axis]) < 1e-11
        z = (e['upper'] if side else e['lower'])[2]
        return (side == 1 and abs(z-field.upper[2]) < 1e-11 or
                side == 0 and (block < d and abs(z) < 1e-11 or block == d+1 and abs(z-H) < 1e-11))
    return field.boundary(face)

def state_difference(a, b, name, tolerance=1e-8):
    require(a['dimension'] == b['dimension'] and a['offsets'] == b['offsets'], 'Incompatible comparison spaces')
    A, B = np.asarray(a['state']), np.asarray(b['state'])
    values = []
    d, us = a['dimension'], a.get('ustar', .4)
    K, kappa, H = a['parameters'][5], a['parameters'][7], a['parameters'][1]
    for block, f in enumerate(spaces(a)):
        offset = a['offsets'][block]
        delta = A[offset:offset+f.size]-B[offset:offset+f.size]
        scale = us if block < d else K if block == d else us**3/(kappa*H) if block == d+1 else us*us
        error = math.sqrt(max(0., float(delta@f.mass()@delta)))/(scale*math.sqrt(f.volume))
        values.append(error)
        require(error <= tolerance, f'{name} block {block}: {error}')
    near(A[-1], B[-1], name+'_gauge', tolerance=tolerance)
    errors[name] = max(values)
    return values


def check_native_spatial(raw):
    f = spaces(raw)
    for field, record in zip(f, raw['native_evaluations']):
        for sample in record['samples']:
            value, gradient = field.evaluate(np.asarray(record['coefficients']), 0, sample[:3])
            near([value, *gradient], sample[3:7], 'independent_tensor_rows')
            near(sample[3:7], sample[7:11], 'native_LocalRegions_evaluation')
    B = np.asarray(raw['pressure_divergence'])
    s = np.linalg.svd(B, compute_uv=False)
    rank = int(np.count_nonzero(s > 1e-10*s[0]))
    require(rank == B.shape[0]-1, f'Unqualified pressure pairing: {rank}, {B.shape}')
    # The common constant pressure is the only null vector, across both regions.
    pressure_constant = np.r_[f[-2].constant(), f[-1].constant()]
    near(pressure_constant@B, np.zeros(B.shape[1]), 'one_common_pressure_nullspace')
    for case in raw['generic_cases']:
        for control in case['constraint_columns']:
            unit=np.zeros(len(case['state']));unit[control['column']]=1.
            near(control['derivative'],unit,'lifted_constraint_Jacobian_columns',tolerance=1e-9)
        body=case['body_partition'];derivative=np.asarray(body['explicit_derivative'])
        near(derivative,body['different_K_derivative'],'prescribed_body_is_state_independent')
        require(abs(body['sampled_K_pair'][1]-body['sampled_K_pair'][0])>1e-3*raw['parameters'][5],'Body control did not move the physical K sampler')
        difference=np.asarray(body['full_minus_implicit']);expected=np.zeros_like(difference)
        for b in range(raw['dimension']):
            start=raw['offsets'][b];expected[start:start+f[b].size]=-f[b].mass()@derivative[start:start+f[b].size]
        for row,_ in case['constraints']:expected[int(row)]=0.
        near(difference,expected,'exact_discrete_IMEX_body_partition',tolerance=1e-10)
    if raw['dimension'] == 2:
        layout = []
        for b in (0, 0, 0, 2, 3, 4, 5):
            r = raw['fields'][b]
            layout.append(dict(global_coefficients=r['size'], mesh=[dict(geometry=e['geometry'], degree=e['degree'], ids=e['ids'], sign=e['sign'],
                bounds=[e['lower'][0], e['upper'][0], e['lower'][2], e['upper'][2]]) for e in r['mesh']]))
        structural = dict(parameters=raw['parameters'], layout=layout, offsets=raw['offsets'], stage_unknowns=len(raw['generic_cases'][0]['state']))
        for case in raw['generic_cases']:
            converted = dict(case, forcing_samples=[[v[0], v[2], v[3]] for v in case['forcing']], potential_cancellation=[])
            residual_oracle(structural, converted)
    else:
        for case in raw['generic_cases']:
            residual_3d(raw, case)
    return dict(dimension=raw['dimension'], pressure_rows=B.shape[0], pressure_rank=rank, singular_values=s.tolist(), generic_chi=[c['chi'] for c in raw['generic_cases']])


def residual_3d(raw, case):
    f, offsets = spaces(raw), raw['offsets']
    v, kf, ef, pw, po = f[0], f[3], f[4], f[5], f[6]
    state, history = np.asarray(case['lifted_state']), np.asarray(case['history'])
    nu, H, rough, rK, rm, K0, qstar, kappa = raw['parameters']
    chi, alpha = case['chi'], case['alpha']
    result = np.zeros(len(state))
    def values(field, geometry, points, block):
        ids, p, g = field.rows(field.by_id[geometry], points)
        c = state[offsets[block]+ids]
        return ids, p, g, p@c, np.einsum('qai,i->qa', g, c)
    def add(block, ids, p, g, w, source, flux):
        np.add.at(result, offsets[block]+ids, p.T@(w*source)+np.einsum('qai,qa,q->i', g, flux, w))
    def integrate(geometry, points, w, inner, profile=None, force=None):
        velocity = [values(v, geometry, points, b) for b in range(3)]
        ids, p, g = velocity[0][:3]
        U, G = np.column_stack([r[3] for r in velocity]), np.stack([r[4] for r in velocity], axis=1)
        ki, kp, kg, K, gK = values(kf, geometry, points, 3)
        ps, pb = (pw, 5) if inner else (po, 6)
        pi, pp, pg, pressure, _ = values(ps, geometry, points, pb)
        x, y, z = points.T
        strain = G+G.transpose(0, 2, 1)
        if inner:
            nt = float(profile.a)*np.maximum(z-float(profile.zv), 0.)
            reynolds = 2*K[:, None, None]/3*np.eye(3)-nt[:, None, None]*strain
            production = -np.sum(reynolds*G, axis=(1, 2))
            stress = (nu+nt)[:, None, None]*strain-pressure[:, None, None]*np.eye(3)
            reference = np.asarray(profile.fields(z)[1], dtype=float)
            original = float(profile.k)**1.5/(2.55*np.maximum(z, float(profile.ze)))
            sk = production-(1-chi)*original-chi*reference+.03+.01*x+.04*y+.02*z
            dk = nu+nt
            body = np.column_stack([.01+.02*x, .015-.01*y+.02*z, -.02+.03*z])
            body[:, :2] += (z <= rough)[:, None]*np.asarray(force)
        else:
            ei, ep, eg, E, gE = values(ef, geometry, points, 4)
            closure = [indexed(a, k, e, 'physical') for a, k, e in zip(G, K, E)]
            reynolds = np.asarray([c['stress'].real for c in closure])
            nt = np.asarray([c['nut'].real for c in closure])
            production = -np.sum(reynolds*G, axis=(1, 2))
            stress = -reynolds+(2*K/3-pressure)[:, None, None]*np.eye(3)
            dk, de = nt, nt/1.3
            sk = production-E-.02+.01*x+.02*y+.03*z
            ratio = K**1.5/(2.5*E*z)
            se = 1.44*E/K*production-1.92*E*E/K+.83*E*E/K*np.maximum(ratio-1, 0)*ratio**2+.01*(1+x+y+z)
            body = np.column_stack([.04*x+.01*z, .02*y-.01*z, -.03+.02*z])
        for b in range(3):
            store = alpha*(U[:, b]-p@history[offsets[b]+ids])
            add(b, ids, p, g, w, store-body[:, b], stress[:, b]-U[:, b, None]*U)
        store = alpha*(K-kp@history[offsets[3]+ki])
        add(3, ki, kp, kg, w, store-sk, dk[:, None]*gK-U*K[:, None])
        add(pb, pi, pp, pg, w, np.trace(G, axis1=1, axis2=2), np.zeros_like(U))
        if not inner:
            store = alpha*(E-ep@history[offsets[4]+ei])
            add(4, ei, ep, eg, w, store-se, de[:, None]*gE-U*E[:, None])
    zt, wt = np.polynomial.legendre.leggauss(10)
    zn, wn = np.polynomial.legendre.leggauss(20)
    columns = sorted({tuple(e['lower'][:2]+e['upper'][:2]) for e in pw.elements})
    for x0, y0, x1, y1 in columns:
        for ix, iy in itertools.product(range(len(zt)), repeat=2):
            x, y = x0+(zt[ix]+1)*(x1-x0)/2, y0+(zt[iy]+1)*(y1-y0)/2
            wh = wt[ix]*wt[iy]*(x1-x0)*(y1-y0)/4
            K, _ = kf.evaluate(state, offsets[3], [x, y, rK], -1 if rK == H else 0)
            profile = Profiles(dict(geometry=[nu,H,rough,rK,rm], metric=[0,0], k=K, q=qstar, kappa=kappa, chi=chi))
            velocities, gradients = [], []
            for b in range(3):
                value, gradient = v.evaluate(state, offsets[b], [x,y,H], 1)
                velocities.append(value); gradients.append(gradient)
            _, gp = po.evaluate(state, offsets[6], [x,y,H], 1)
            _, gk = kf.evaluate(state, offsets[3], [x,y,H], 1)
            force = [sum(velocities[b]*gradients[a][b]+velocities[a]*gradients[b][b] for b in range(2))+gp[a]-2*gk[a]/3-
                     ([.04*x+.01*H,.02*y-.01*H][a]) for a in range(2)]
            breaks = {0.,H,rough,float(profile.rc),float(profile.ze)}
            if profile.zv > 0:breaks.add(float(profile.zv))
            start = min(float(profile.rc), nu/float(profile.a), float(profile.ze))
            while start < H:breaks.add(start); start *= 3
            for e in pw.elements:breaks.update((e['lower'][2],e['upper'][2]))
            breaks = sorted(r for r in breaks if 0 <= r <= H)
            for a,b in zip(breaks[:-1],breaks[1:]):
                points = np.column_stack([np.full_like(zn,x),np.full_like(zn,y),a+(zn+1)*(b-a)/2])
                geometry = pw.element([x,y,(a+b)/2])['geometry']
                integrate(geometry,points,wh*wn*(b-a)/2,True,profile,force)
    for e in po.elements:
        points,w = po.quadrature(e,12)
        integrate(e['geometry'],points,w,False)
    volume = v.volume
    for block, field in ((5,pw),(6,po)):
        load = field.mass()@field.constant()/volume
        result[offsets[block]:offsets[block]+field.size] += state[offsets[7]]*load
        result[offsets[7]] += state[offsets[block]:offsets[block]+field.size]@load
    for row,value in case['constraints']:result[int(row)] = case['state'][int(row)]-value
    for block in range(7):
        near(np.asarray(case['residual'])[offsets[block]:offsets[block+1]],result[offsets[block]:offsets[block+1]],'independent_3D_residual_'+str(block),tolerance=2e-10)
    near(case['residual'][-1],result[-1],'independent_3D_gauge')

def radial_projection(raw, block, profile):
    """Independent constrained one-dimensional physical L2 minimization."""
    field = spaces(raw)[block]
    H, us, kappa = raw['parameters'][1], raw['ustar'], raw['parameters'][7]
    Z = field.upper[2]
    intervals = sorted({(e['lower'][2],e['upper'][2],e['degree']) for e in field.elements})
    mapping, local = {}, []
    for a,b,p in intervals:
        ids=[]
        for j in range(p+1):
            key=(a,a,0) if j==0 else (b,b,0) if j==1 else (a,b,j)
            if key not in mapping:mapping[key]=len(mapping)
            ids.append(mapping[key])
        local.append(ids)
    n=len(mapping);M=np.zeros((n,n));rhs=np.zeros(n)
    z,w=np.polynomial.legendre.leggauss(64)
    velocity=block < raw['dimension']
    def target(r):
        r=np.asarray(r)
        if velocity:
            j0,_=wall_integrals(profile,np.minimum(r,H))
            return np.asarray(float(profile.t0)*j0+np.where(r>H,us/kappa*np.log(np.maximum(r,H)/H),0.),dtype=float)
        return us**3/(kappa*r)
    for (a,b,p),ids in zip(intervals,local):
        breaks={a,b}
        start=1e-5
        while start < b:
            if start>a:breaks.add(start)
            start*=3
        for v in (float(profile.zv),H):
            if a<v<b:breaks.add(v)
        cuts=sorted(breaks)
        for lower,upper in zip(cuts[:-1],cuts[1:]):
            r=lower+(z+1)*(upper-lower)/2;W=w*(upper-lower)/2
            phi,_=basis(p,2*(r-a)/(b-a)-1)
            M[np.ix_(ids,ids)] += phi.T@(W[:,None]*phi)
            rhs[ids] += phi.T@(W*target(r))
    D_H=float(profile.nu+profile.a*max(H-profile.zv,0))
    constraints=([(0,0,0),(H,0,float(profile.UH)),(H,-1,us*us/D_H),(H,1,us/(kappa*H)),
                  (Z,0,float(target([Z])[0])),(Z,-1,us/(kappa*Z))] if velocity else
                 [(H,0,us**3/(kappa*H)),(Z,0,us**3/(kappa*Z)),(Z,-1,-us**3/(kappa*Z*Z))])
    C,values=[],[]
    for height,side,value in constraints:
        i=next(i for i,(a,b,p) in enumerate(intervals) if a-1e-12<=height<=b+1e-12 and
               not(side>0 and height>=b-1e-12 or side<0 and height<=a+1e-12))
        a,b,p=intervals[i];phi,dp=basis(p,[2*(height-a)/(b-a)-1]);row=np.zeros(n)
        row[local[i]]=phi[0] if side==0 else dp[0]*2/(b-a)
        norm=np.linalg.norm(row);C.append(row/norm);values.append(value/norm)
    C=np.asarray(C);values=np.asarray(values)
    # Scale the objective by the physical interval volume; this changes neither
    # its minimizer nor the imposed physical trace/derivative constraints.
    volume=Z-field.lower[2]
    system=np.block([[M/volume,C.T],[C,np.zeros((len(C),len(C)))]])
    solution=np.linalg.solve(system,np.r_[rhs/volume,values])[:n]
    x,y=(field.lower[:2]+field.upper[:2])/2
    direction=raw['direction'][block] if velocity else 1.
    samples=[]
    for (a,b,p),ids in zip(intervals,local):
        for r in a+(np.linspace(-.91,.91,5)+1)*(b-a)/2:
            phi,_=basis(p,[2*(r-a)/(b-a)-1])
            expected=float((phi@solution[ids])[0])*direction
            actual,_=field.evaluate(np.asarray(raw['fixed_reference']),raw['offsets'][block],[x,y,r])
            samples.append([r,actual,expected])
    near(np.asarray(samples)[:,1],np.asarray(samples)[:,2],'constrained_physical_L2_projection',scale=us if velocity else us**3/(kappa*H),tolerance=1e-9)
    return dict(radial_unknowns=n,constraints=len(C),max_sample_difference=float(np.max(np.abs(np.asarray(samples)[:,1]-np.asarray(samples)[:,2]))))


def check_reference(raw, projection=False, dual=False):
    d=raw['dimension'];f=spaces(raw);o=raw['offsets'];s=np.asarray(raw['fixed_reference']);defect=np.asarray(raw['fixed_defect'])
    nu,H,rough,rK,rm,K,q,kappa=raw['parameters'];us=raw['ustar'];Z=f[0].upper[2]
    profile=Profiles(dict(geometry=[nu,H,rough,rK,rm],metric=[0,0],k=K,q=q,kappa=kappa,chi=1))
    near(s[o[d]:o[d]+f[d].size],K*f[d].constant(),'exact_constant_k')
    near(s[o[d-1]:o[d-1]+f[d-1].size],np.zeros(f[d-1].size),'exact_zero_normal_velocity')
    near(defect[o[d+2]:],np.zeros(len(defect)-o[d+2]),'anchor_excludes_all_algebraic_rows')
    exported_constraints={int(row):value for row,value in raw['reference_constraints']}
    raw_load=np.asarray(raw['raw_reference_loads'])
    removed=[]
    integrals=[];work=0.;independent_norms=[]
    for b,field in enumerate(f):
        c=field.constant();load=defect[o[b]:o[b]+field.size]
        boundary=constrained(raw,field,b)
        require({o[b]+i for i in boundary}=={i for i in exported_constraints if o[b]<=i<o[b]+field.size},'Native/independent essential row sets differ')
        if boundary:near(load[sorted(boundary)],np.zeros(len(boundary)),'anchor_excludes_essential_and_dependent_rows')
        if b<d+2:
            mask=np.asarray([i for i in range(field.size) if i not in boundary],dtype=int)
            near(load[mask],raw_load[o[b]+mask],'only_fixed_free_differential_loads_are_stored')
            boundary_load=np.asarray([raw_load[o[b]+i] for i in sorted(boundary)])
            removed.append(dict(block=b,row_count=len(boundary),maximum_raw_load=float(np.max(np.abs(boundary_load),initial=0.)),
                                constant_test_load=float(sum(field.constant()[i]*raw_load[o[b]+i] for i in boundary)),
                                raw_load_sha256=hashlib.sha256(boundary_load.tobytes()).hexdigest()))
        if b<d+2:integrals.append(float(c@load))
        if b<d:work += float(np.asarray(raw['state'])[o[b]:o[b]+field.size]@load)
        if b==d:work += float(c@load)
        if dual:
            free=np.asarray([i for i in range(field.size) if i not in boundary],dtype=int)
            M=field.mass()[np.ix_(free,free)];L=load[free]
            scale=us if b<d else K if b==d else us**3/(kappa*H) if b==d+1 else 1.
            rate=us/H if b>=d+2 else scale*us/H
            independent_norms.append(math.sqrt(max(0.,float(L@np.linalg.solve(M,L))))/(rate*math.sqrt(field.volume)) if len(free) else 0.)
    near(integrals,raw['defect_integrals'],'free_defect_integrals',tolerance=1e-9)
    near(raw['eta']*work,raw['anchor_energy_work'],'anchor_energy_work',tolerance=1e-9)
    if dual:
        near(independent_norms[:d+2],raw['raw_reference_norms'][:d+2],'physical_mass_dual_defect_norms',tolerance=1e-9)
    require(max(raw['raw_reference_norms'][d+2:])<=1e-9,'Reference algebraic constraints failed before anchoring')
    if raw['eta']==1:require(max(raw['corrected_reference_norms'])<=1e-9,'Anchored spatial reference is not a fixed point')
    else:near(raw['corrected_reference_norms'],raw['raw_reference_norms'],'eta0_leaves_operator_unchanged')
    gauge=0.;max_div=0.;Uerror=Eerror=0.
    for b in (d+2,d+3):
        gauge += float((f[b].mass()@f[b].constant())@s[o[b]:o[b]+f[b].size])/f[0].volume
    near(gauge,0.,'physical_global_pressure_gauge',scale=us*us,tolerance=1e-9)
    axes=(0,2) if d==2 else (0,1,2)
    for b,field in enumerate(f):
        for e in field.elements:
            n=e['degree']+1;local=s[o[b]+np.asarray(e['ids'])]*np.asarray(e['sign'])
            expected=np.zeros_like(local)
            for j in range(n**d):
                radial=j//n**(d-1)
                horizontal=j%n**(d-1)
                if all((horizontal//n**i)%n<2 for i in range(d-1)):
                    expected[j]=local[radial*n**(d-1)]
            near(local,expected,'native_exact_horizontal_invariance',tolerance=1e-9)
    for e in f[0].elements:
        p=[(e['lower'][a]+e['upper'][a])/2 for a in range(3)]
        G=np.zeros((3,3))
        for b,a in enumerate(axes):
            _,G[a]=f[b].evaluate(s,o[b],p)
        max_div=max(max_div,abs(np.trace(G)))
    require(max_div/(us/H)<=1e-9,'Reference discrete divergence failed')
    x=(f[0].lower[0]+f[0].upper[0])/2;y=(f[0].lower[1]+f[0].upper[1])/2
    for xi in (.17,.51,.83):
        p=[f[0].lower[0]+xi*(f[0].upper[0]-f[0].lower[0]),y,H]
        Gi,Go=np.zeros((3,3)),np.zeros((3,3))
        for b,a in enumerate(axes):
            ui,Gi[a]=f[b].evaluate(s,o[b],p,-1);uo,Go[a]=f[b].evaluate(s,o[b],p,1)
            expected=raw['direction'][b]*float(profile.UH) if b<d-1 else 0.
            near([ui,uo],[expected,expected],'shared_matching_velocity_trace',scale=us,tolerance=1e-9)
        ki,gki=f[d].evaluate(s,o[d],p,-1);ko,gko=f[d].evaluate(s,o[d],p,1)
        E,_=f[d+1].evaluate(s,o[d+1],p,1)
        pi,_=f[d+2].evaluate(s,o[d+2],p,-1);po,_=f[d+3].evaluate(s,o[d+3],p,1)
        closure=indexed(Go,ko,E,'physical');R=closure['stress'].real
        stress_i=float(nu+profile.a*max(H-profile.zv,0))*(Gi+Gi.T)-pi*np.eye(3)
        stress_o=-R+(2*ko/3-po)*np.eye(3)
        near(stress_i[:,2],stress_o[:,2],'full_matching_traction',scale=us*us,tolerance=1e-9)
        near(pi-po,R[2,2]-2*K/3,'required_relative_pressure_offset',scale=us*us,tolerance=1e-9)
        near(pi-po,-11*K/75,'CLS_reference_pressure_jump',scale=us*us,tolerance=1e-9)
        near([ki,ko,E],[K,K,us**3/(kappa*H)],'matching_k_and_dependent_epsilon',scale=K,tolerance=1e-9)
        near(gki,np.zeros(3),'zero_inner_k_flux',tolerance=1e-9);near(gko,np.zeros(3),'zero_outer_k_flux',tolerance=1e-9)
        for b in range(d):
            ground,_=f[b].evaluate(s,o[b],[p[0],p[1],0.])
            near(ground,0.,'strong_ground_no_slip',scale=us,tolerance=1e-9)
        top=[p[0],p[1],Z];Et,gEt=f[d+1].evaluate(s,o[d+1],top,-1)
        nt=us*kappa*Z
        near([Et,nt*gEt[2]/1.3],[us**3/(kappa*Z),-us**4/(1.3*Z)],'compatible_top_epsilon_trace_and_flux',scale=us**3/(kappa*H),tolerance=1e-9)
    # Physical profile errors use a separate high-order rule, including graded
    # panels near the molecular AWF scale. They measure projection, not drift.
    for e in f[0].elements:
        cuts={e['lower'][2],e['upper'][2]};start=1e-5
        area_fraction=math.prod((e['upper'][a]-e['lower'][a])/(f[0].upper[a]-f[0].lower[a]) for a in axes if a!=2)
        while start<e['upper'][2]:
            if start>e['lower'][2]:cuts.add(start)
            start*=3
        z,w=np.polynomial.legendre.leggauss(32)
        for a,b in zip(sorted(cuts)[:-1],sorted(cuts)[1:]):
            r=a+(z+1)*(b-a)/2
            points=np.column_stack([np.full_like(r,x),np.full_like(r,y),r])
            _,phi,_=f[0].rows(e,points)
            actual=phi@s[o[0]+np.asarray(e['ids'])]
            j0,_=wall_integrals(profile,np.minimum(r,H))
            exact=raw['direction'][0]*(float(profile.t0)*j0+np.where(r>H,us/kappa*np.log(np.maximum(r,H)/H),0.))
            Uerror+=area_fraction*float(np.dot(w*(b-a)/2,(actual-exact)**2))
            if a>=H-1e-12:
                ee=f[d+1].by_id[e['geometry']];ids,ep,_=f[d+1].rows(ee,points)
                Eerror+=area_fraction*float(np.dot(w*(b-a)/2,(ep@s[o[d+1]+ids]-us**3/(kappa*r))**2))
    result=dict(raw_norms=raw['raw_reference_norms'],corrected_norms=raw['corrected_reference_norms'],defect_integrals=integrals,anchor_energy_work=raw['anchor_energy_work'],
                velocity_profile_RMS=math.sqrt(Uerror/Z),epsilon_profile_RMS=math.sqrt(Eerror/(Z-H)),maximum_divergence=max_div,unknowns=len(s),removed_row_audit=removed)
    if projection:
        result['independent_velocity_projection']=radial_projection(raw,0,profile)
        result['independent_epsilon_projection']=radial_projection(raw,d+1,profile)
    return result


def manufactured_error(raw):
    result=[];s=np.asarray(raw['state']);o=raw['offsets'];scales=[.4,.4,raw['parameters'][5],.4**3/(raw['parameters'][7]*2),.16,.16]
    for block,field in enumerate(spaces(raw)):
        square=0.
        for e in field.elements:
            points,w=field.quadrature(e,12);ids,p,_=field.rows(e,points)
            exact=np.array([manufactured_values(x,r,raw['time'],raw['parameters'])[block] for x,_,r in points])
            delta=p@s[o[block]+ids]-exact;square+=float(np.dot(w,delta*delta))
        result.append(math.sqrt(square/field.volume)/scales[block])
    return result


def history_difference(a,b,name):
    near(a['native_history_times'],b['native_history_times'],name+'_times',tolerance=1e-8)
    near(a['native_history'],b['native_history'],name+'_all_native_slots',tolerance=1e-8)
    near(a['fixed_reference'],b['fixed_reference'],name+'_stored_reference',tolerance=1e-10)
    near(a['fixed_defect'],b['fixed_defect'],name+'_stored_free_defect',tolerance=1e-10)
    return state_difference(a,b,name,tolerance=1e-8)

class Runner:
    def __init__(self, build, output):
        self.build, self.output = build.resolve(), output.resolve()
        self.output.mkdir(parents=True, exist_ok=False)
        self.commands=[]
        self.environment=dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
        self.environment['LD_LIBRARY_PATH']=str(self.build/'geometry-lib')+':'+self.environment.get('LD_LIBRARY_PATH','')
        self.template=Path(__file__).resolve().parent/'EquilibriumWallFlat.xml'

    def command(self, name, binary, arguments, ranks=1, timeout=1800, expect_failure=False):
        directory=self.output/name
        directory.mkdir(exist_ok=True)
        cmd=([str(self.build/binary)]+list(map(str,arguments)))
        if ranks>1:cmd=['mpirun','--report-bindings','--bind-to','core','--map-by','core','-np',str(ranks),*cmd]
        start=time.monotonic()
        try:
            with (directory/'run.log').open('w') as log:
                completed=subprocess.run(cmd,cwd=directory,env=self.environment,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
            code=completed.returncode
        except subprocess.TimeoutExpired:code=124
        self.commands.append(dict(name=name,command=cmd,ranks=ranks,threads=THREADS,exit_code=code,seconds=time.monotonic()-start,
                                  log=str(directory/'run.log'),log_sha256=digest(directory/'run.log'),expected_failure=expect_failure))
        (self.output/'commands.json').write_text(json.dumps(self.commands,indent=2)+'\n')
        require((code!=0 and code!=124) if expect_failure else code==0, f'{name} exit {code}: '+(directory/'run.log').read_text()[-1800:])
        print(f'{name}: '+('rejected as required' if expect_failure else 'passed')+f' ({time.monotonic()-start:.2f}s)',flush=True)
        return directory

    def evolution(self, name, generator=neutral_session, ranks=1, expect_failure=False, **kwargs):
        directory=self.output/name;directory.mkdir()
        xml=generator(self.template,directory/'case.xml',**kwargs)
        self.command(name,'IncNavierStokesSolver',[xml],ranks,expect_failure=expect_failure)
        if expect_failure:return None
        raw=json.loads((directory/'case.evolution.json').read_text())
        require(raw['native_time_map'] and raw['ranks']==ranks,'Missing native map/collective execution')
        require(len(raw['native_history'])==2*raw['order'],'Incomplete exported native history')
        for stage in raw['stage_audits']:
            require(stage[4]<=1e-9 and max(stage[6:])<=1e-9,'Unconverged coupled native IMEX stage')
        near(raw['native_history'][0],raw['state'][:raw['offsets'][raw['dimension']+2]],'actual_native_output_state')
        return raw

    def spatial(self,name,dimension):
        directory=self.output/name;directory.mkdir()
        root=ET.parse(self.template).getroot()
        if dimension==3:tensor_mesh(root,[0.,.5,2.,3.,4.],3,3,1.,1.)
        set_parameter(root,'RANSWallTangentialOrder',12)
        xml=write(root,directory/'case.xml');report=directory/'probe.json'
        self.command(name,'RANSEquilibriumEvolutionProbe',[xml,report])
        return json.loads(report.read_text())


def graded_knots(profile, factor):
    ell=float(profile.nu/profile.a-profile.zv)
    values={0.,float(profile.h),float(profile.H)};r=ell
    while r<float(profile.H):values.add(r);r=(r+ell)*factor-ell
    return sorted(values)


def qualification(run, section):
    report={}
    if section in ('all','spatial'):
        report['native_spatial']=[check_native_spatial(run.spatial('spatial-'+str(d),d)) for d in (2,3)]
    if section in ('all','reference'):
        # Mandatory unanchored evaluation precedes every anchored run.
        base=run.evolution('unanchored-coarse',eta=0,xknots=[0,150,300])
        report['unanchored_coarse']=check_reference(base,projection=True,dual=True)
        nu,H,h,rK,rm,K,q,kappa=base['parameters']
        profile=Profiles(dict(geometry=[nu,H,h,rK,rm],metric=[0,0],k=K,q=q,kappa=kappa,chi=1))
        outer=[H]
        while outer[-1]*1.6<30:outer.append(outer[-1]*1.6)
        outer.append(30.)
        p_results=[]
        for p in (3,5,7):
            raw=run.evolution('unanchored-p'+str(p),degree=p,inner_knots=graded_knots(profile,2.5),outer_knots=outer,eta=0,tangential_order=max(8,p+2))
            p_results.append(dict(degree=p,**check_reference(raw)))
        h_results=[]
        for i,factor in enumerate((4.,2.5,1.6)):
            raw=run.evolution('unanchored-h'+str(i),degree=4,inner_knots=graded_knots(profile,factor),outer_knots=outer,eta=0,tangential_order=8)
            h_results.append(dict(grading_factor=factor,**check_reference(raw)))
        for name,values in (('p',p_results),('h',h_results)):
            require(values[-1]['velocity_profile_RMS']<values[0]['velocity_profile_RMS'],name+' projection lacks refinement')
            require(max(values[-1]['raw_norms'][:4])<max(values[0]['raw_norms'][:4]),name+' unanchored residual lacks refinement')
        require(p_results[-1]['epsilon_profile_RMS']<p_results[0]['epsilon_profile_RMS']/10,'Outer dissipation projection lacks p refinement')
        report['unanchored_p_refinement']=p_results;report['unanchored_h_refinement']=h_results
        anchored=run.evolution('anchored-spatial',eta=1,xknots=[0,150,300])
        report['anchored_spatial']=check_reference(anchored,dual=True)
        near(base['fixed_reference'],anchored['fixed_reference'],'eta_independent_projected_target')
        near(base['fixed_defect'],anchored['fixed_defect'],'fixed_global_anchor_loads')
        require(max(abs(v) for v in anchored['defect_integrals'])>1e-10,'Budget fixture has no retained constant test')
        require(abs(anchored['anchor_energy_work'])>1e-10,'Budget fixture has no retained energy work')
    if section in ('all','unanchored-time'):
        baseline=[]
        for dt,steps in ((.01,1),(.01,3),(.005,6)):
            raw=run.evolution(f'unanchored-dt{dt}-n{steps}',eta=0,order=1,dt=dt,steps=steps,xknots=[0,150,300])
            near(raw['anchor_energy_work'],0.,'unanchored_zero_anchor_work')
            near(raw['raw_reference_norms'],raw['corrected_reference_norms'],'unanchored_original_operator')
            require(max(row[2] for row in raw['step_audits'])>1e-8,'Unanchored target was silently held fixed')
            baseline.append(dict(order=1,dt=dt,steps=steps,time=raw['time'],final_scaled_field_drift=raw['step_audits'][-1][4:],
                                 max_stage_residual=max(row[4] for row in raw['stage_audits']),anchor_energy_work=raw['anchor_energy_work']))
        report['unanchored_native_time_maps']=baseline
    if section in ('all','stationary'):
        cases=[];baseline=None
        for order in (1,2):
            for dt in (.1,.05,.025):
                raw=run.evolution(f'stationary-o{order}-dt{dt}',eta=1,steps=round(1/dt),order=order,dt=dt,xknots=[0,150,300])
                check_reference(raw)
                require(max(row[2] for row in raw['step_audits'])<=1e-9,'Repeated native timestep drifts from anchored reference')
                cases.append(dict(order=order,dt=dt,steps=raw['steps'],time=raw['time'],max_drift=max(row[2] for row in raw['step_audits']),max_stage_residual=max(row[4] for row in raw['stage_audits'])))
                if order==2 and dt==.05:baseline=raw
        long=run.evolution('stationary-one-turnover',eta=1,steps=100,dt=.05,condensed=True,xknots=[0,150,300])
        require(max(row[2] for row in long['step_audits'])<=1e-9,'One-turnover native timestep drifts')
        report['stationary_time_maps']=cases;report['one_turnover']=dict(steps=long['steps'],time=long['time'],max_drift=max(row[2] for row in long['step_audits']))
        def neutral_body(template,destination,**kwargs):
            xml=neutral_session(template,destination,**kwargs);root=ET.parse(xml).getroot()
            for name,expressions in (('BodyForceInner',('0.01*cos(0.03*x)+0.004*sin(31*y)','0.002*cos(7*y)')),
                                     ('BodyForceOuter',('0.005*sin(0.02*x)+0.003*cos(17*y)','0.001*sin(3*y)'))):
                node=ET.SubElement(root.find('CONDITIONS'),'FUNCTION',{'NAME':name})
                for variable,expression in zip(('u','v'),expressions):ET.SubElement(node,'E',{'VAR':variable,'VALUE':expression})
            return write(root,xml)
        body_cases=[]
        for order in (1,2):
            body=run.evolution('stationary-prescribed-body-o'+str(order),generator=neutral_body,eta=1,steps=6,order=order,dt=.05,xknots=[0,150,300])
            require(max(row[2] for row in body['step_audits'])<=1e-9,'Explicit prescribed force breaks the anchored stationary map')
            require(max(row[1] for row in body['explicit_audits'])>1e-8,'Stationary force case has no actual explicit derivative history')
            body_cases.append(dict(order=order,max_drift=max(row[2] for row in body['step_audits']),
                                   explicit_derivative_maximum=max(row[1] for row in body['explicit_audits']),**check_reference(body,dual=True)))
        report['stationary_nonzero_explicit_body']=body_cases
        variants=[]
        for d,direction in ((3,(1.,0.)),(3,(.6,.8))):
            label='3D-straight' if direction[1]==0 else '3D-rotated'
            serial=run.evolution(label,dimension=d,direction=direction,eta=1,steps=8,dt=.05,xknots=[0,150,300])
            check_reference(serial,projection=direction[1]!=0)
            require(max(row[2] for row in serial['step_audits'])<=1e-9,'3D stationary timestep drifts')
            entries=[]
            for ranks in (2,4):
                parallel=run.evolution(label+'-mpi'+str(ranks),dimension=d,direction=direction,eta=1,steps=8,dt=.05,xknots=[0,150,300],ranks=ranks)
                entries.append(dict(ranks=ranks,state_difference=history_difference(parallel,serial,label+'-mpi'+str(ranks))))
                require(max(row[2] for row in parallel['step_audits'])<=1e-9,'MPI stationary timestep drifts')
            variants.append(dict(dimension=d,direction=direction,max_drift=max(row[2] for row in serial['step_audits']),MPI=entries))
        report['full_3D_and_rotation']=variants
        # A stationary anchored state also restores the stored nonzero defect.
        saved=run.evolution('stationary-save',eta=1,steps=3,dt=.05,xknots=[0,150,300],checkpoint=run.output/'stationary-history')
        resumed=run.evolution('stationary-restart-mpi4',eta=1,steps=17,dt=.05,xknots=[0,150,300],restart=saved['checkpoint'],ranks=4)
        report['stationary_restart']=history_difference(resumed,baseline,'anchored_restart')
    if section in ('all','evolving'):
        rows=[];baselines={}
        for order in (1,2):
            errors_by_dt=[]
            for dt in (.1,.05,.025):
                steps=round(.4/dt)
                raw=run.evolution(f'evolving-o{order}-dt{dt}',generator=evolving_session,order=order,dt=dt,steps=steps)
                field_errors=manufactured_error(raw);errors_by_dt.append(field_errors)
                baselines[(order,dt)]=raw
                require(max(v[1] for v in raw['explicit_audits'])>1e-8,'Evolving case does not exercise native explicit history')
            rates=(np.log(np.asarray(errors_by_dt[:-1])/np.asarray(errors_by_dt[1:]))/math.log(2)).tolist()
            # Check every differential field. Pressure is an algebraic multiplier;
            # its measured rates are reported separately, with a single gauge.
            minimum=.8 if order==1 else 1.7
            require(min(np.asarray(rates)[:,:4].ravel())>=minimum,f'Native IMEX{order} lacks its differential-field order: {rates}')
            rows.append(dict(order=order,dt=[.1,.05,.025],scaled_physical_field_errors=errors_by_dt,observed_orders=rates))
        report['manufactured_time_refinement']=rows
        base=baselines[(2,.1)]
        condensed=run.evolution('evolving-condensed',generator=evolving_session,steps=4,dt=.1,condensed=True)
        report['evolving_condensation']=history_difference(condensed,base,'full_condensed_time_maps')
        mpi=[]
        for ranks in (2,4):
            raw=run.evolution('evolving-mpi'+str(ranks),generator=evolving_session,steps=4,dt=.1,ranks=ranks)
            mpi.append(dict(ranks=ranks,state_difference=history_difference(raw,base,'evolving_mpi'+str(ranks)),iterations=[r[3] for r in raw['stage_audits']]))
        report['evolving_MPI']=mpi
        restarts=[]
        for order,split,save_ranks,ranks in ((1,1,1,1),(2,0,1,1),(2,1,1,4),(2,2,4,1),(2,3,1,4)):
            saved=run.evolution(f'evolving-o{order}-save{split}',generator=evolving_session,order=order,steps=split,dt=.1,ranks=save_ranks,
                                checkpoint=run.output/f'evolving-history-o{order}-{split}')
            resumed=run.evolution(f'evolving-o{order}-restart{split}',generator=evolving_session,order=order,steps=4-split,dt=.1,restart=saved['checkpoint'],ranks=ranks)
            restarts.append(dict(order=order,split=split,save_ranks=save_ranks,restart_ranks=ranks,
                                 state_difference=history_difference(resumed,baselines[(order,.1)],f'evolving_o{order}_restart{split}')))
        report['full_history_restarts']=restarts
        # Source/time/eta compatibility are exercised at the coupled system in
        # addition to corruption/shape/order controls in the native time probe.
        path=saved['checkpoint']
        run.evolution('restart-invalid-dt',generator=evolving_session,order=2,steps=1,dt=.05,restart=path,expect_failure=True)
        def wrong_source(template,destination,**kwargs):
            xml=evolving_session(template,destination,**kwargs);root=ET.parse(xml).getroot()
            node=root.find("CONDITIONS/FUNCTION[@NAME='BodyForceOuter']/E[@VAR='u']")
            node.set('VALUE','('+node.get('VALUE')+')+0.0001');return write(root,xml)
        run.evolution('restart-invalid-source',generator=wrong_source,steps=1,dt=.1,restart=path,expect_failure=True)
        def illegal_anchor(template,destination,**kwargs):
            xml=evolving_session(template,destination,**kwargs);root=ET.parse(xml).getroot();set_parameter(root,'RANSWallEta',1);return write(root,xml)
        run.evolution('manufactured-invalid-anchor',generator=illegal_anchor,steps=1,dt=.1,expect_failure=True)
        report['coupled_rejected_controls']=['dt','source','eta1_outside_flat_neutral_verification']
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',required=True,type=Path);p.add_argument('--output',type=Path)
    p.add_argument('--section',choices=('all','spatial','reference','unanchored-time','stationary','evolving'),default='all')
    args=p.parse_args()
    require(platform.system()=='Linux','Native numerical qualification requires Linux')
    source=Path(__file__).resolve().parents[4]
    revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    source_status=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
    require(not source_status,'Qualification requires a clean source worktree')
    output=args.output or args.build/('cls-epawf-m4-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f'))
    run=Runner(args.build,output)
    started=time.monotonic()
    result=dict(schema=1,milestone='M4',source_commit=revision,source_status=source_status,started_utc=datetime.now(timezone.utc).isoformat(),
                platform=platform.platform(),python=sys.version,numpy=np.__version__,scipy=scipy.__version__,thread_environment=THREADS,
                scope='native flat 2D/full3D CG, replicated global coefficients, cell-owned MPI assembly, root dense solve',
                section=args.section,milestone_closure=False,verification_status='running')
    scope=[source/'docs/plans/cls-epawf-revised-contract-v1.json',source/'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex']
    scope+=list((source/'solvers/IncNavierStokesSolver/EquationSystems').glob('RANS*'))
    scope+=list((source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('*m[1234]*.py'))
    scope+=list((source/'library/SolverUtils/Turbulence').glob('*.*'))
    scope+=list((source/'library/UnitTests/SolverUtils').glob('RANSEquilibrium*Probe.cpp'))
    result['source_sha256']={str(path.relative_to(source)):digest(path) for path in scope if path.is_file()}
    result['binary_sha256']={name:digest(args.build/name) for name in ('IncNavierStokesSolver','RANSEquilibriumEvolutionProbe','RANSEquilibriumTimeProbe')}
    try:
        result.update(qualification(run,args.section))
        # Native scalar recurrence/startup/history tests are independent of the
        # spatial operator and must accompany a complete M4 qualification.
        if args.section=='all':
            native=run.output/'native-time';native.mkdir()
            run.command('native-time','RANSEquilibriumTimeProbe',['--output',native/'time.json','--checkpoint-root',native/'history'])
            result['native_time_probe']=json.loads((native/'time.json').read_text())
        result['verification_status']='passed';result['milestone_closure']=args.section=='all'
    except Exception as error:
        result['verification_status']='failed';result['failure']=str(error)
        raise
    finally:
        result.update(checks=checks,maximum_scaled_errors=errors,commands=run.commands,seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat())
        result['source_end_commit']=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
        result['source_end_status']=subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
        (run.output/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    require(result['source_end_commit']==revision and not result['source_end_status'],'Source changed during qualification')
    print(f'CLS EP-AWF M4: {checks} checks; report {run.output/"validation.json"}',flush=True)

if __name__=='__main__':main()
