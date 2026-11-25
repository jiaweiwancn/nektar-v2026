#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently sample the planar principal response on a normal hill."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.polynomial import polyroots

from c1_oracle import Audit, load_c1_state, regularization_viscosity
import physical_oracle as po
from resources import identity


DIRECTIONS = np.zeros((3, 3, 3))
DIRECTIONS[0, 0, 0] = 1.
DIRECTIONS[0, 2, 2] = -1.
DIRECTIONS[1, 0, 2] = 1.
DIRECTIONS[2, 2, 0] = 1.


def angular_minimum(action, shift):
    # For t = tan(theta), the coefficients of the perpendicular
    # polarization tensor are (-t, -t*t, 1)/(1+t*t).
    q = action + action.swapaxes(1, 2)
    p = np.column_stack((action[:, 2, 2], -q[:, 0, 2],
                         action[:, 0, 0]-q[:, 1, 2], q[:, 0, 1], action[:, 1, 1]))
    phi = np.arange(32)*2*np.pi/32
    c, s = np.cos(phi), np.sin(phi)
    v = np.column_stack((-s/2, (c-1)/2, (c+1)/2))
    grid = np.einsum('ij,njk,ik->ni', v, action, v).min(axis=1) + shift
    # Taylor's theorem at an angular stationary point bounds how far the
    # grid could overestimate the minimum. Use the bound to screen points.
    v0, vc, vs = np.array([[0., -.5, .5], [0., .5, .5], [-.5, 0., 0.]])
    b = v0@q@vc
    c = v0@q@vs
    d = (vc@action@vc-vs@action@vs)/2
    e = (vc@q@vs)/2
    lower = grid - .5*(np.pi/32)**2*(abs(b)+abs(c)+4*(abs(d)+abs(e)))
    best = float(grid.min())
    winner = int(grid.argmin())
    for index in np.flatnonzero(lower <= best+5e-15):
        # Each chart is evaluated only for |t| <= 1. The reversed
        # polynomial includes angles where cos(theta) vanishes.
        for coeff in (p[index], p[index, ::-1]):
            a0, a1, a2, a3, a4 = coeff
            roots = polyroots([a1, 2*a2-4*a0, 3*a3-3*a1, 4*a4-2*a2, -a3])
            t = np.r_[-1., 0., 1., roots.real[(abs(roots.imag)<1e-9)&(abs(roots.real)<=1.)]]
            value = np.polynomial.polynomial.polyval(t, coeff)/(1+t*t)**2 + shift[index]
            if value.min() < best:
                best, winner = float(value.min()), int(index)
    return best, winner


def response(ref, G, k, E, radius, metric, outer, added_viscosity=None):
    normal = np.column_stack((metric['n'][:, 0], np.zeros(len(k)), metric['n'][:, 1]))
    tangent = np.column_stack((metric['t'][:, 0], np.zeros(len(k)), metric['t'][:, 1]))
    size = len(k)
    full = np.zeros((size, 3, 3))
    full[:, 0, 0], full[:, 0, 2] = G[:, 0, 0], G[:, 0, 1]
    full[:, 2, 0], full[:, 2, 2] = G[:, 1, 0], G[:, 1, 1]
    collar = (radius/ref.H-.4)/.35
    chi = np.where(collar <= 0., 0., np.where(collar >= 1., 1.,
                   collar**3*(10-15*collar+6*collar**2)))
    aux = np.zeros(size) if outer else ref.aux(metric['K'], radius)
    mu = ref.nu+(1-chi)*aux
    added = np.zeros(size) if added_viscosity is None else added_viscosity
    # The lower wall's exact angular minimum is attained at 45 degrees
    # to its frame. The collar uses independent full-WJ differences.
    lower = np.flatnonzero(chi == 0.)
    best = (float('inf'), 0)
    difference = 0.
    if len(lower):
        values = (ref.nu+added[lower])/mu[lower]
        winner = int(values.argmin())
        best = float(values[winner]), int(lower[winner])
    active = np.flatnonzero(chi > 0.)
    if len(active):
        g, ka, ea = full[active], k[active], E[active]
        weight = chi[active]
        _, root, _, _ = po.full_wj(ka, ea, g)
        omega = (.5*ka/ea)[:, None, None]*(g-g.swapaxes(1, 2))
        inv_w = np.einsum('nij,nji->n', omega, omega)
        mu[active] += weight*.6*ka*ka/ea*root/(root*root-2*inv_w)
        t, n = tangent[active], normal[active]
        forms = []
        for relative in (1e-3, 5e-4):
            delta = (relative*ea/ka)[:, None, None]
            rows = []
            for direction in DIRECTIONS:
                derivative = sum(w*po.full_wj(ka, ea, g+j*delta*direction)[0]
                                 for j, w in ((-2, 1), (-1, -8), (1, 8), (2, -1)))/(12*delta)
                M = direction+direction.T
                tangent_strain = np.einsum('ni,ij,nj->n', t, M, n)
                wall = (t[:, :, None]*n[:, None, :]+n[:, :, None]*t[:, None, :])
                scale = ((1-weight)*aux[active]*tangent_strain)[:, None, None]
                stress = ref.nu*M-weight[:, None, None]*derivative+scale*wall
                rows.append(np.einsum('aij,nij->na', DIRECTIONS, stress)/mu[active, None])
            forms.append(np.stack(rows, axis=2))
        difference = float(np.max(np.sum(abs(forms[1]-forms[0]), axis=(1, 2))))
        value, winner = angular_minimum(forms[1], added[active]/mu[active])
        if value < best[0]:
            best = value, int(active[winner])
    value, index = best
    return value, index, difference, float(aux[index]), float(chi[index])


def audit_state(data, ref, order, radial_order, deadline, regularization=None, one_sided_faces=False):
    audit = Audit(data, ref, order, deadline, radial_order)
    results = {}
    for element in audit.mesh:
        audit.clock()
        a, r, _ = audit.rule(element)
        outer = bool(element['epsilon'])
        groups = [('outer' if outer else 'inner', a, r)]
        if min(abs(element['bounds'][j]-ref.H) for j in (2, 3)) < 1e-12:
            lo, hi = element['bounds'][:2]
            face_a = lo+(hi-lo)*(audit.gauss[0]+1)/2
            groups.append(('matching_outer' if outer else 'matching_inner',
                           face_a, np.full_like(face_a, ref.H)))
        if one_sided_faces:
            _, _, rl, rr = element["bounds"]
            name = "faces_outer" if outer else "faces_inner"
            for side in range(4):
                along, radial, _ = audit.face(element, side, endpoints=True)
                groups.append((name, along, radial))
            if not outer:
                cutoff = ref.cutoff(audit.trace(element, "k", along))
                inside = (cutoff >= rl) & (cutoff <= rr)
                if inside.any():
                    groups.append(("cutoff", along[inside], cutoff[inside]))
        for group, a, r in groups:
            metric = audit.hill.metrics(element, a, r, outer=outer)
            _, _, G, _, _, _, scalars, _, _, covariance = audit.sample(element, a, r, metric)
            k, E = scalars['k'], scalars['epsilon']
            metric['K'] = audit.trace(element, 'k', a)
            added = regularization_viscosity(audit.geometry, regularization, metric["position"], r)
            principal, index, difference, aux, chi = response(ref, G, k, E, r, metric, outer, added)
            result = results.setdefault(group, dict(samples=0, minimum_k=float('inf'),
                minimum_epsilon=float('inf'), covariance_over_k=float('inf'),
                principal_over_mu=float('inf'), max_directional_difference=0.))
            result['samples'] += len(a)
            result['minimum_k'] = min(result['minimum_k'], float(k.min()))
            result['minimum_epsilon'] = min(result['minimum_epsilon'], float(E.min()))
            result['covariance_over_k'] = min(result['covariance_over_k'], float(covariance))
            result['max_directional_difference'] = max(result['max_directional_difference'], difference)
            if principal < result['principal_over_mu']:
                result['principal_over_mu'] = principal
                result['witness'] = dict(geometry=element['geometry'], a=float(a[index]),
                    rho=float(r[index]), position=metric['position'][index].tolist(),
                    element_bounds=element["bounds"],
                    gradient=G[index].tolist(), k=float(k[index]), epsilon=float(E[index]),
                    auxiliary=aux, chi=chi)
    for result in results.values():
        result['sample_pass'] = (result['minimum_k'] > 0. and result['minimum_epsilon'] > 0.
            and result['covariance_over_k'] > 1e-10
            and result['principal_over_mu']-result['max_directional_difference'] > 1e-8)
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('native', 'table', 'out'):
        parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--order', type=int, default=20)
    parser.add_argument('--radial-order', type=int)
    parser.add_argument('--seconds', type=float, default=180.)
    parser.add_argument('--diagnostic', action='store_true')
    parser.add_argument('--regularized', action='store_true')
    parser.add_argument('--one-sided-faces', action='store_true',
                        help='also scan both polynomial owners of faces, corners, and the cutoff')
    args = parser.parse_args()
    began = time.monotonic()
    ref = po.load_reference(args.table)
    data, metadata = load_c1_state(args.native, args.diagnostic)
    if args.regularized and (not args.diagnostic or "regularization" not in metadata):
        parser.error("--regularized requires a diagnostic steady candidate")
    regularization = metadata["regularization"] if args.regularized else None
    groups = audit_state(data, ref, args.order, args.radial_order, began+args.seconds,
                         regularization, args.one_sided_faces)
    report = dict(native=identity(args.native), table=identity(args.table), audit_source=identity(__file__),
                  time=data['time'], step=data['step'], height=data['geometry']['hill_height'],
                  order=args.order, radial_order=args.radial_order or args.order,
                  one_sided_faces=args.one_sided_faces, regions=groups,
                  elapsed=time.monotonic()-began)
    if args.diagnostic:
        report.update(candidate_role=data['diagnostic_role'], regularization=regularization)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps({key: value['principal_over_mu'] for key, value in groups.items()}))


if __name__ == '__main__':
    main()
