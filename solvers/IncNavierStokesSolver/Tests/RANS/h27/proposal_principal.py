#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""PROPOSAL ONLY: independent isotropic-normal-wall/full-WJ-collar screen.

Run on Linux with python3 -B under an admitted analysis allowance. This script
does not modify the active PDE or qualify Task3. Positive finite direction/point
coverage is not a continuum certificate. Output is immutable exclusive-create.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from wj_h5_neutral import NeutralReference, KAPPA

IDENTITY = np.eye(3)
WALL_NORMAL = np.array([0., 0., 1.])
PROJECTOR = IDENTITY-np.outer(WALL_NORMAL, WALL_NORMAL)
CP = 1.8


def full_wj(k, e, gradient):
    tau = k/e
    s = tau*((gradient+gradient.T)/2-np.trace(gradient)*IDENTITY/3)
    o = tau*(gradient-gradient.T)/2
    s2, o2 = s@s, o@o
    invs, invw, inv4 = np.trace(s2), np.trace(o2), np.trace(s@o2)
    roots = np.roots([1., -CP, -2*invw-2.7*invs, 2*CP*invw])
    physical = [float(v.real) for v in roots
                if abs(v.imag) <= 1e-10*max(1., abs(v.real)) and v.real >= CP*(1-1e-12)]
    if len(physical) != 1:
        raise ValueError('proposal-only WJ physical root absent or ambiguous')
    n = physical[0]
    denominator = (5/6)*(n*n-2*invw)*(2*n*n-invw)
    b1 = -n*(2*n*n-7*invw)/denominator
    b3 = -12*inv4/(n*denominator)
    b4 = -2*(n*n-2*invw)/denominator
    b6 = -6*n/denominator
    b9 = 6/denominator
    a = (b1*s+b3*(o2-invw*IDENTITY/3)+b4*(s@o-o@s)
         +b6*(s@o2+o2@s-2*inv4*IDENTITY/3)
         +b9*(o@s@o2-o2@s@o))
    return k*(2*IDENTITY/3+a)


def chi(y):
    if y <= .4:
        return 0.
    if y >= .75:
        return 1.
    s = (y-.4)/.35
    return s**3*(10-15*s+6*s*s)


def proposal_tensor(k, e, gradient, auxiliary, blend):
    b = -auxiliary*PROJECTOR@(gradient+gradient.T)@WALL_NORMAL
    iso = 2*k*IDENTITY/3+np.outer(b, WALL_NORMAL)+np.outer(WALL_NORMAL, b)
    return iso if blend == 0 else (1-blend)*iso+blend*full_wj(k, e, gradient)


def normal_inventory(count):
    normals = []
    for i in range(count):
        nz = 1-2*(i+.5)/count
        angle = i*math.pi*(3-math.sqrt(5))
        normals.append([math.sqrt(1-nz*nz)*math.cos(angle),
                        math.sqrt(1-nz*nz)*math.sin(angle), nz])
    return np.array(normals)


def polarizations(normal):
    axis = IDENTITY[int(np.argmin(abs(normal)))]
    a = np.cross(normal, axis)
    a /= np.linalg.norm(a)
    return np.array([a, np.cross(normal, a)])


def derivative_tensor(c, k, e, gradient, auxiliary, blend, relative):
    result = np.zeros((3, 3, 3, 3))
    step = relative*e/k
    for a in range(3):
        for b in range(3):
            h = np.zeros((3, 3)); h[a, b] = 1.
            dr = (proposal_tensor(k, e, gradient+step*h, auxiliary, blend)
                  -proposal_tensor(k, e, gradient-step*h, auxiliary, blend))/(2*step)
            result[:, :, a, b] = c.nu*(h+h.T)-dr
    return result


def screen_point(c, r, cutoff, owner, sphere, planar):
    blend = chi(r/c.H)
    auxiliary = c.nt(r)
    if owner == 'outer':
        k, e, g = c.K, c.A/r, c.u_tau/(KAPPA*r)
    else:
        k, e, g = c.k_inner(r), c.A/max(r, cutoff), c.u_tau**2/c.D(r)
    gradient = np.zeros((3, 3)); gradient[0, 2] = g
    stress = proposal_tensor(k, e, gradient, auxiliary, blend)
    wj = full_wj(k, e, gradient)
    scale = c.nu+(1-blend)*auxiliary+blend*(-wj[0, 2]/g)
    b = -auxiliary*PROJECTOR@(gradient+gradient.T)@WALL_NORMAL
    expected_parallel = ((1-blend)*auxiliary*g-blend*wj[0, 2])*g
    steps = []
    for relative in (3e-4, 1e-4):
        action = derivative_tensor(c, k, e, gradient, auxiliary, blend, relative)
        records = {}
        identity_error = 0.
        for name, normals in [('planar', planar), ('full3d', sphere)]:
            minimum = (float('inf'), None)
            in_plane = (float('inf'), None)
            for normal in normals:
                basis = polarizations(normal)
                matrix = np.zeros((2, 2))
                for j, a in enumerate(basis):
                    h = np.outer(a, normal)
                    av = np.einsum('ijkl,kl->ij', action, h)
                    matrix[:, j] = basis@(av@normal)
                eigenvalues, eigenvectors = np.linalg.eigh((matrix+matrix.T)/2)
                quality = float(eigenvalues[0]/scale)
                if quality < minimum[0]:
                    minimum = (quality, dict(normal=normal.tolist(),
                        polarization=(eigenvectors[:, 0]@basis).tolist(),
                        normalized_principal=quality))
                if name == 'planar':
                    a = np.array([-normal[2], 0., normal[0]])
                    h = np.outer(a, normal)
                    value = float(np.sum(h*np.einsum('ijkl,kl->ij', action, h))/scale)
                    if value < in_plane[0]:
                        in_plane = (value, dict(normal=normal.tolist(), polarization=a.tolist(),
                                               normalized_principal=value))
                if blend == 0:
                    sh = [PROJECTOR@(np.outer(a, normal)+np.outer(normal, a))@WALL_NORMAL for a in basis]
                    exact = c.nu*np.eye(2)+auxiliary*np.array([[v@w for w in sh] for v in sh])
                    identity_error = max(identity_error, float(np.max(abs(matrix-exact))/scale))
            records[name+'_two_mode'] = minimum[1]
            if name == 'planar':
                records['planar_in_plane'] = in_plane[1]
        steps.append(dict(relative=relative, minima=records,
                          chi_zero_exact_coercivity_matrix_error=identity_error if blend == 0 else None))
    difference = max(abs(steps[0]['minima'][key]['normalized_principal']
                        -steps[1]['minima'][key]['normalized_principal'])
                     for key in steps[0]['minima'])
    return dict(owner=owner, r=float(r), physical_z=float(r+c.d0), chi=blend,
                k=k, epsilon=e, G=g, q=k*g/e, auxiliary=auxiliary, scale=scale,
                covariance_eigenvalues_over_k=(np.linalg.eigvalsh(stress)/k).tolist(),
                iso_covariance_margin_over_k=float((2*k/3-np.linalg.norm(b))/k),
                trace_error_over_k=float(abs(np.trace(stress)-2*k)/k),
                parallel_production_error_over_e=float(abs(-np.sum(stress*gradient)-expected_parallel)/e),
                principal_steps=steps, minimum_step_difference=difference)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--cutoff', type=float, action='append', default=[])
    parser.add_argument('--reference-json', type=Path)
    parser.add_argument('--radius-count', type=int, default=25)
    parser.add_argument('--sphere-directions', type=int, default=256)
    parser.add_argument('--planar-directions', type=int, default=128)
    args = parser.parse_args()
    if min(args.radius_count, args.sphere_directions, args.planar_directions) < 8:
        parser.error('at least eight points/directions per inventory')
    repo = Path(__file__).resolve().parents[5]
    calibration = repo/'docs/plans/wj-h5-linux/calibration.json'
    historical = next(v for v in json.loads(calibration.read_text())['candidates'] if v['H_w_mm'] == 6.)
    c = NeutralReference(historical['H_w'], historical['z0_e'])
    for name, value in [('K', c.K), ('A', c.A), ('h_s', c.h), ('nu', c.nu), ('z0_e', c.z0_e)]:
        if float(value).hex() != float(historical[name]).hex():
            raise ValueError('retained inlet calibration identity mismatch '+name)
    cutoffs = list(args.cutoff)
    if args.reference_json:
        for case in json.loads(args.reference_json.read_text())['cases']:
            if abs(case.get('parameters', {}).get('K', -1)/c.K-1) < 1e-12:
                cutoffs.extend(v['audit']['c'] for v in case.get('roots', []) if 'audit' in v)
    if not cutoffs or any(not 0 < v <= c.H for v in cutoffs):
        parser.error('supply actual positive cutoff <=H or a reference-oracle JSON with matching K')
    sphere = normal_inventory(args.sphere_directions)
    planar = np.array([[math.cos(t), 0., math.sin(t)]
                       for t in np.linspace(0., math.pi, args.planar_directions, endpoint=False)])
    cases = []
    for cutoff in sorted(set(cutoffs)):
        radii = sorted(set(np.linspace(0, c.H, args.radius_count).tolist()
                          +[.0015, .003, .006, .4*c.H, .5*c.H, .6*c.H, .75*c.H, cutoff, c.h, c.w['rc']]))
        rows = [screen_point(c, r, cutoff, 'inner', sphere, planar) for r in radii]
        rows.extend(screen_point(c, float(r), cutoff, 'outer', sphere, planar)
                    for r in np.geomspace(c.H, c.top-c.d0, 17))
        minima = {key: min((dict(r=row['r'], owner=row['owner'], chi=row['chi'], **step['minima'][key])
                           for row in rows for step in row['principal_steps']),
                          key=lambda v: v['normalized_principal'])
                  for key in rows[0]['principal_steps'][0]['minima']}
        covariance = min(v['covariance_eigenvalues_over_k'][0] for v in rows)
        stable_steps = max(v['minimum_step_difference'] for v in rows) <= 1e-6
        positive = all(v['normalized_principal'] > 0 for v in minima.values()) and covariance >= 0
        cases.append(dict(c=cutoff, c_over_H=cutoff/c.H, sampled_points=len(rows), minima=minima,
                          minimum_covariance_over_k=covariance, fd_step_agreement_1e_minus_6=stable_steps,
                          status='POSITIVE_FINITE_COVERAGE_ONLY' if positive and stable_steps else 'NOT_MET_OR_UNRESOLVED',
                          states=rows))
    output = dict(schema='WJH27ProposalOnlyPrincipalScreenV1', proposal_only=True,
                  active_model_modified=False, task3_ready=False,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  calibration_sha256=hashlib.sha256(calibration.read_bytes()).hexdigest(),
                  tensor='(1-chi)*(2kI/3+b_A tensor n+n tensor b_A)+chi*R_W',
                  chi='unchanged geometry C2: zero<=.4H, full WJ>=.75H',
                  sphere_directions=args.sphere_directions, planar_directions=args.planar_directions,
                  caveat='Finite frozen-scalar constitutive coverage only, no continuum or coupled-PDE qualification.',
                  cases=cases)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(json.dumps(output, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(proposal_only=True, output=str(args.output),
                          cases=[dict(c=v['c'], status=v['status'], minima=v['minima']) for v in cases]), sort_keys=True))


if __name__ == '__main__':
    main()
