#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent continuous P1/P2 pressure projection and tensor-polynomial action.

Native value rows supply the represented polynomial, never its expected
derivative. The rough term uses the containing OUTER element at (actual x,H).
"""
import os
from rans_qualification import THREADS
os.environ.update(THREADS)
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.polynomial import polyvander2d
from ct_accepted_open import require, row_value
from ct_accepted2_diagnosis import analytic, WITNESS


def lagrange(t, degree):
    nodes = np.linspace(0., 1., degree+1)
    return np.array([np.prod([(t-b)/(a-b) for b in nodes if b != a], axis=0)
                     for a in nodes])


def lagrange_derivative(t, degree):
    nodes = np.linspace(0., 1., degree+1)
    return np.array([sum(np.prod([(t-b)/(a-b) for j,b in enumerate(nodes)
                                  if j != i and j != m], axis=0)/(a-nodes[m])
                         for m in range(len(nodes)) if m != i)
                     for i,a in enumerate(nodes)])


def projection(root, c, degree, o, witness_row):
    z = np.array([x for x in c['knots'] if x >= c['H']])
    n = (len(z)-1)*degree+1
    mass, load = np.zeros((n,n)), np.zeros(n)
    g, w = leggauss(8)
    for j,(a,b) in enumerate(zip(z[:-1], z[1:])):
        t = (g+1)/2
        ph = lagrange(t, degree)
        ids = np.arange(j*degree, j*degree+degree+1)
        p = analytic(root, 'Outer', 'p', a+(b-a)*t)
        mass[np.ix_(ids,ids)] += (b-a)/2*(ph*w)@ph.T
        load[ids] += (b-a)/2*(ph@(w*p))
    coeff = np.linalg.solve(mass, load)
    j = np.searchsorted(z, WITNESS[2], side='right')-1
    t = (WITNESS[2]-z[j])/(z[j+1]-z[j])
    grad = lagrange_derivative(t, degree)@coeff[j*degree:j*degree+degree+1]/(z[j+1]-z[j])
    native = row_value(witness_row['pressure'], np.array(o['state']), o['offsets'][5])[1][2]
    require(abs(grad-native) < 2e-11, 'degree_correct_pressure_projection_mismatch')
    return dict(status='pass', degree=degree, independent_gradient_at_witness=float(grad),
                native_gradient_at_witness=float(native), maximum_defect=abs(float(grad-native)),
                coefficients=coeff.tolist(), nodes=z.tolist(), projection_gauss_order=8,
                scope='continuous nodal P1/P2 mass projection of unchanged analytical Pi; no P1 shortcut for P2')


def projected_normal_load(root, c, projected, order):
    """Independent radial integral, analytical Qzz=-Pi vs represented pressure.

    This is the initial supplied continuum/pressure-projection diagnostic, not
    the actual evolved stress load or an artificial governing residual.
    """
    z=np.array(projected['nodes']);coeff=np.array(projected['coefficients']);degree=projected['degree']
    g,w=leggauss(order);total=0.
    for i,(a,b) in enumerate(zip(z[:-1],z[1:])):
        t=(g+1)/2;y=a+(b-a)*t
        native=lagrange(t,degree).T@coeff[i*degree:i*degree+degree+1]
        derivative=-2*(y-6.25)/1.5**2*np.exp(-((y-6.25)/1.5)**2)
        total+=(b-a)/2*(w@((analytic(root,'Outer','p',y)-native)*derivative))
    return float(total*(c['x'][-1]-c['x'][0])**3/6)


class TensorPressure:
    def __init__(self, o, state, degree):
        require(degree in (1,2), 'pressure_oracle_only_P1_P2')
        self.degree, self.elements = degree, {e['geometry']:e for e in o['elements']}
        self.coefficients = {}
        off = o['offsets']
        groups = {}
        for q in o['pressure_mass']:
            groups.setdefault(q['pressure']['geometry'], []).append(q)
        for gid, rows in groups.items():
            e = self.elements[gid]
            xy = np.array([[2*(q['coordinate'][a]-e['lower'][a])/(e['upper'][a]-e['lower'][a])-1
                            for a in (0,2)] for q in rows])
            V = polyvander2d(xy[:,0],xy[:,1],[degree,degree]).reshape(len(rows),-1)
            y = np.array([row_value(q['pressure'],state,off[4+q['region']])[0] for q in rows])
            v,_,rank,_ = np.linalg.lstsq(V,y,rcond=None)
            require(rank == (degree+1)**2, 'pressure_tensor_fit_rank')
            require(np.max(abs(V@v-y)) < 2e-10*max(1.,np.max(abs(y))), 'pressure_tensor_not_represented')
            self.coefficients[gid] = v.reshape(degree+1,degree+1)

    def evaluate(self, gid, x, z):
        e = self.elements[gid]
        # Known containing geometry, physical coordinates unchanged. Endpoint
        # arithmetic allowed only at the declared H trace, not neighbor search.
        require(e['lower'][0] <= x <= e['upper'][0] and e['lower'][2]-1e-14 <= z <= e['upper'][2]+1e-14,
                'pressure_polynomial_wrong_owner')
        tx = 2*(x-e['lower'][0])/(e['upper'][0]-e['lower'][0])-1
        tz = 2*(z-e['lower'][2])/(e['upper'][2]-e['lower'][2])-1
        p,dx,dz = 0.,0.,0.
        for i in range(self.degree+1):
            for j in range(self.degree+1):
                c = self.coefficients[gid][i,j]
                p += c*tx**i*tz**j
                if i: dx += c*i*tx**(i-1)*tz**j*2/(e['upper'][0]-e['lower'][0])
                if j: dz += c*j*tx**i*tz**(j-1)*2/(e['upper'][2]-e['lower'][2])
        return p,dx,dz


def pressure_action(o, direction, hr, degree):
    state = np.asarray(direction)
    off = o['offsets']
    poly = TensorPressure(o,state,degree)
    result, rough = np.zeros(len(state)), np.zeros(len(state))
    H = min(q['coordinate'][2] for q in o['faces'] if q['interface'])
    gids = {q['pressure']['geometry'] for q in o['faces'] if q['interface'] and q['region']==1}
    provenance = []
    for q in o['volume']:
        v, pr = q['velocity_k'], q['pressure']
        ids = np.array(v['ids'],int)
        p = row_value(pr,state,off[4+q['region']])[0]
        for b,a in enumerate((0,2)):
            result[ids+off[b]] -= q['weight']*p*np.array(v['derivative'][a])
        if q['region']==0 and q['coordinate'][2] <= hr:
            x = q['coordinate'][0]
            owners = [gid for gid in gids if poly.elements[gid]['lower'][0] < x < poly.elements[gid]['upper'][0]]
            require(len(owners)==1, 'rough_volume_outer_trace_owner_ambiguous')
            gid = owners[0]
            dp = poly.evaluate(gid,x,H)[1]
            rough[ids+off[0]] -= q['weight']*dp*np.array(v['phi'])
            provenance.append(dict(coordinate=q['coordinate'], generating_geometry=q['generating_geometry'],
                                   trace_coordinate=[x,0.,H], trace_geometry=gid, trace_side=1, dpdx=float(dp)))
    return result+rough, rough, provenance
