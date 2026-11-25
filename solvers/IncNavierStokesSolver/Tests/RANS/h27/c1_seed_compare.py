#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare source and target chart velocities and independently test seed continuity."""

import argparse
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu

from c1_geometry import Hill
from physical_oracle import SparseMass, basis, load_native, require
from physical_transient_oracle import state_hash
from resources import identity


def smooth(value, low, high):
    t = np.clip((value-low)/(high-low), 0., 1.)
    return t**3*(10+t*(-15+6*t))


def fade(value, low, high):
    return smooth(high-value, 0., high-low)


def gradient(element, a, r, metric, flow, n):
    _, da, dr = basis(element, 'u', a, r)
    ids = np.asarray(element['u'])
    u = np.column_stack((flow[ids], flow[n+ids]))
    return ((da@u)[:, :, None]*metric['da'][:, None, :] +
            (dr@u)[:, :, None]*metric['dr'][:, None, :])


def rotate_gradient(gradient, source, target):
    old = np.stack((source['t'], source['n']), axis=2)
    new = np.stack((target['t'], target['n']), axis=2)
    rotation = new@old.swapaxes(1, 2)
    return rotation@gradient@rotation.swapaxes(1, 2)


def coefficients(state, prepared):
    fields = state['coefficients']
    expected = state['state_sha256']
    if prepared and state['step'] > 0:
        require(state['step'] == 1, 'prepared coefficients require the first BE frame')
        fields = state['input_history']['current']
        expected = state['previous_sha256']
    arrays = {key: np.asarray(values,dtype=float) for key, values in fields.items()}
    require(state_hash(arrays) == expected, 'seed coefficient identity differs')
    return arrays


def compare(source, target, order, prepared, point):
    n, np_, _, _ = target['sizes']
    require(source['sizes'] == target['sizes'], 'seed comparison requires the same signed grid')
    s, t = coefficients(source, False), coefficients(target, prepared)
    sh, th = Hill(source['geometry']), Hill(target['geometry'])
    q, w = leggauss(order)
    qa, qr = (x.ravel() for x in np.meshgrid(q,q,indexing='ij'))
    product = np.outer(w,w).ravel()/4
    prev = {e['geometry']: e for e in source['elements']}
    divergence = np.zeros(np_)
    mass = SparseMass(np_)
    area, strong = 0., 0.
    moments = np.zeros((2, 3, 2))
    stations = []
    for e in target['elements']:
        src = prev[e['geometry']]
        require(all(e[key] == src[key] for key in ('bounds','u','us','p','ps')),
                'seed comparison found a changed signed element')
        xl, xr, rl, rr = e['bounds']
        a, r = xl+(xr-xl)*(qa+1)/2, rl+(rr-rl)*(qr+1)/2
        m, sm = th.metrics(e,a,r), sh.metrics(src,a,r)
        weight = product*(xr-xl)*(rr-rl)*m['det']
        G = gradient(e,a,r,m,t['flow'],n)
        before = gradient(src,a,r,sm,s['flow'],n)
        d = G[:,0,0]+G[:,1,1]
        p, _, _ = basis(e, 'p', a,r)
        np.add.at(divergence,e['p'],p.T@(weight*d))
        mass.add(e['p'],p.T@(weight[:,None]*p))
        area += weight.sum()
        strong += weight@(d*d)
        if rl >= th.H-1e-13:
            patch = 0
            window = (fade(abs(a-th.L),.25,.80)*smooth(r,4*th.H,5.5*th.H)*
                      fade(r,9*th.H,11*th.H))
        else:
            patch = 1
            before = rotate_gradient(before,sm,m)
            window = (fade(abs(a),.2*th.L,.6*th.L)*smooth(r,.15*th.H,.4*th.H)*
                      fade(r,.8*th.H,th.H))
        error = (G-before)**2
        measure = window*weight
        moments[patch,0] += (measure@error.sum(axis=(1,2)), measure.sum())
        moments[patch,1] += (measure@error[:,:,0].sum(axis=1), measure.sum())
        moments[patch,2] += (measure@error[:,:,1].sum(axis=1), measure.sum())
        if xl <= point[0] <= xr and rl <= point[1] <= rr:
            a, r = (np.asarray([value]) for value in point)
            m, sm = th.metrics(e,a,r), sh.metrics(src,a,r)
            before = gradient(src,a,r,sm,s['flow'],n)
            if patch:
                before = rotate_gradient(before,sm,m)
            stations.append(dict(element=e['geometry'], physical=m['position'][0].tolist(),
                source_gradient=before[0].tolist(), target_gradient=gradient(e,a,r,m,t['flow'],n)[0].tolist()))
    riesz = splu(mass.finish()).solve(divergence)
    lift = np.asarray(target['essential_flow_lift'])
    fixed = [np.setdiff1d(np.arange(n), free)+i*n for i,free in enumerate(target['free_velocity_ids'])]
    return dict(weak_continuity=float(np.sqrt(max(0.,divergence@riesz)/area)),
        strong_divergence_rms=float(np.sqrt(strong/area)),
        gradient_rms={name: dict(zip(('full','x_derivative','z_derivative'),
                     np.sqrt(row[:,0]/row[:,1]).tolist()))
                      for name,row in zip(('foot','collar'),moments)},
        essential_velocity_max_error=float(np.max(np.abs((t['flow']-lift)[np.concatenate(fixed)]))),
        pressure_max_change=float(np.max(np.abs(t['flow'][2*n:]-s['flow'][2*n:]))),
        k_max_change=float(np.max(np.abs(t['k']-s['k']))),
        chart_point=point, station=stations)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source','target','out'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--order', type=int, default=24)
    parser.add_argument('--prepared', action='store_true')
    parser.add_argument('--point', type=float, nargs=2, default=(2.4884753964041995,1.1168173532711385))
    args = parser.parse_args()
    require(args.order >= 12, 'seed comparison requires at least 12 Gauss points per direction')
    source, _ = load_native(args.source)
    target, _ = load_native(args.target)
    require(source['geometry']['mapping']==target['geometry']['mapping']=='NormalCosineHillV1',
            'seed comparison requires physical normal charts')
    report = dict(source=identity(args.source), target=identity(args.target), audit_source=identity(__file__),
                  order=args.order, prepared=args.prepared, frame_step=target['step'],
                  source_height=source['geometry']['hill_height'], target_height=target['geometry']['hill_height'],
                  **compare(source,target,args.order,args.prepared,args.point))
    with args.out.open('x') as stream:
        json.dump(report,stream,indent=2,allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
