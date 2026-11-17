#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent tiny native Q3 scalar algebra verifier.

Uses the H26 independent Jacobi basis definition, never native sampled basis or
operators. Physical Gauss20x24 rules split the wall-distance coordinate at the
moving cutoff, collar boundaries, roughness height, interface and cell bounds.
This is a scalar foundation gate only, not coupled-model qualification.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
FIXTURE_SOURCE = Path(__file__).resolve().parents[1]/'h26'/'fixtures.py'
sys.path.insert(0, str(FIXTURE_SOURCE.parent))
import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
from fixtures import modified

EDGE_MODES = ((0, 4, 8, 12), (1, 5, 9, 13), (0, 1, 2, 3), (4, 5, 6, 7))
NORMALS = ((-1., 0.), (1., 0.), (0., -1.), (0., 1.))


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def vector(data, name, size):
    result = np.asarray(data[name], dtype=float)
    require(result.shape == (size,) and np.all(np.isfinite(result)),
            name+' vector dimension/data differ')
    return result


def matrix(packet, size):
    require(packet['rows'] == size and packet['cols'] == size, 'matrix shape differs')
    result = np.zeros((size, size))
    for i, j, value in packet['entries']:
        require(0 <= i < size and 0 <= j < size, 'matrix index outside own scalar map')
        result[i, j] += value
    require(np.all(np.isfinite(result)), 'nonfinite native matrix')
    return result


def relative(actual, reference):
    scale = la.norm(reference)
    require(scale > 0, 'relative comparison needs a nonzero reference')
    return float(la.norm(actual-reference)/scale)


def comparison(checks, name, actual, reference, tolerance=1e-10):
    error = relative(actual, reference)
    checks[name] = dict(relative_error=error, tolerance=tolerance)
    require(error <= tolerance, name+' independent comparison failed')
    return error


def basis(element, x, r, signs):
    xl, xr, rl, rr = element['bounds']
    bx, dx = modified(3, 2*(x-xl)/(xr-xl)-1)
    bz, dz = modified(3, 2*(r-rl)/(rr-rl)-1)
    b = (bz[:, :, None]*bx[:, None, :]).reshape(len(x), 16)*signs
    gx = (bz[:, :, None]*dx[:, None, :]).reshape(len(x), 16)*signs*2/(xr-xl)
    gz = (dz[:, :, None]*bx[:, None, :]).reshape(len(x), 16)*signs*2/(rr-rl)
    return b, gx, gz


def analytic_q(x, r):
    return .01*(1+.2*x+.1*r+.3*x**3*(1+.1*r**2))


def analytic_gradient(x, r):
    return .01*(.2+.9*x*x*(1+.1*r*r)), .01*(.1+.06*x**3*r)


def volume_rule(element, height, roughness, nx, nz):
    ax, wx = leggauss(nx); ar, wr = leggauss(nz)
    xl, xr, rl, rr = element['bounds']
    xx, zz, ww = [], [], []
    for xi, wi in zip(ax, wx):
        x = (xl+xr)/2+(xr-xl)*xi/2
        cutoff = .03*(1+.03*x)
        cuts = sorted(set([rl, rr]+[v for v in
                      (cutoff, .4*height, .75*height, roughness, height) if rl < v < rr]))
        for lo, hi in zip(cuts[:-1], cuts[1:]):
            xx.extend(np.full(nz, x)); zz.extend((lo+hi)/2+(hi-lo)*ar/2)
            ww.extend(wi*wr*(xr-xl)*(hi-lo)/4)
    return np.asarray(xx), np.asarray(zz), np.asarray(ww)


def face_rule(element, side, order):
    abscissa, weight = leggauss(order)
    xl, xr, rl, rr = element['bounds']
    if side < 2:
        return (np.full(order, xl if side == 0 else xr),
                (rl+rr)/2+(rr-rl)*abscissa/2, weight*(rr-rl)/2)
    return ((xl+xr)/2+(xr-xl)*abscissa/2,
            np.full(order, rl if side == 2 else rr), weight*(xr-xl)/2)


def inventory(data, case):
    ns, nu = data['scalar_size'], data['velocity_size']
    require(ns > 0 and nu > 0, 'empty native scalar/velocity map')
    flow = {e['geometry']: e for e in data['flow_elements']}
    cells = data['scalar_elements']
    require(len(flow) == len(data['flow_elements']) and len(cells) > 0,
            'duplicate/empty geometry inventory')
    require(len({e['geometry'] for e in cells}) == len(cells), 'duplicate scalar geometry')
    seen_scalar, seen_velocity = set(), set()
    for elements, key, sign_key, size, seen in (
            (data['flow_elements'], 'u', 'us', nu, seen_velocity),
            (cells, 's', 'ss', ns, seen_scalar)):
        for e in elements:
            ids = e[key]; signs = np.asarray(e[sign_key])
            require(len(ids) == 16 and len(set(ids)) == 16 and
                    all(isinstance(i, int) and 0 <= i < size for i in ids),
                    'invalid own native Q3 coefficient map')
            require(signs.shape == (16,) and np.all(np.isin(signs, [-1., 1.])),
                    'invalid own native Q3 orientation signs')
            xl, xr, rl, rr = e['bounds']
            require(np.all(np.isfinite(e['bounds'])) and xr > xl and rr > rl,
                    'nonpositive/nonfinite rectangular geometry')
            seen.update(ids)
    require(seen_scalar == set(range(ns)) and seen_velocity == set(range(nu)),
            'native coefficient map has unvisited global IDs')
    for e in cells:
        require(e['geometry'] in flow and
                np.max(abs(np.asarray(e['bounds'])-flow[e['geometry']]['bounds'])) < 1e-12,
                'scalar geometry differs from its flow owner')
    expected_geometries = {g for g, e in flow.items()
                           if case == 'full' or e['bounds'][2] >= data['H']-1e-12}
    require({e['geometry'] for e in cells} == expected_geometries,
            'scalar domain does not match the selected full/outer geometry set')
    xmin = min(e['bounds'][0] for e in cells); xmax = max(e['bounds'][1] for e in cells)
    rmin = min(e['bounds'][2] for e in cells); rmax = max(e['bounds'][3] for e in cells)
    require(abs(xmin) < 1e-12, 'fixture requires x_left=0 for compatible gamma corner')
    require(abs(rmin-(0. if case == 'full' else data['H'])) < 1e-12,
            'wrong ground/interface wall distance')
    area = sum((e['bounds'][1]-e['bounds'][0])*(e['bounds'][3]-e['bounds'][2]) for e in cells)
    require(abs(area-(xmax-xmin)*(rmax-rmin)) < 1e-11*area,
            'scalar subset is not a complete rectangle')
    faces, essential = [], set()
    essential_sides = (0,) if case == 'full' else (0, 2)
    for e in cells:
        xl, xr, rl, rr = e['bounds']
        exterior = (abs(xl-xmin) < 1e-12, abs(xr-xmax) < 1e-12,
                    abs(rl-rmin) < 1e-12, abs(rr-rmax) < 1e-12)
        for side, boundary in enumerate(exterior):
            if boundary:
                faces.append((e, side))
                if side in essential_sides:
                    essential.update(e['s'][j] for j in EDGE_MODES[side])
    free = np.asarray(data['free_ids'], dtype=int)
    require(len(set(free)) == len(free) and set(free) == set(range(ns))-essential,
            'explicit essential/free mask differs from independent scalar edge map')
    require(data['essential_count'] == len(essential) > 0, 'essential scalar count differs')
    if case == 'outer':
        require(data['native_dirichlet_count'] == 0,
                'outer fixture must have no native XML Dirichlet mask')
        require(all(e['bounds'][2] >= data['H']-1e-12 for e in cells),
                'outer epsilon includes an inner geometry')
    else:
        require(data['native_dirichlet_count'] == len(essential),
                'full fixture native left boundary count differs')
    return flow, cells, faces, np.asarray(sorted(essential), dtype=int), free, area


def assemble(data, case, nx, nz, flow, cells, faces, target, velocity):
    ns, nu = data['scalar_size'], data['velocity_size']
    m = np.zeros((ns, ns)); k = np.zeros_like(m); riesz = np.zeros_like(m); transport = np.zeros_like(m)
    adv = np.zeros(ns); advective = np.zeros(ns); strong = np.zeros(ns)
    boundary_adv = np.zeros(ns); source = np.zeros(ns); flux = np.zeros(ns)
    one = np.zeros(ns); one_seen = {}; pointwise = dict(q=0., qx=0., qr=0., u=0., w=0., divergence=0.)
    volume = work = source_work = boundary_flow = natural_work = 0.
    for e in cells:
        ids = np.asarray(e['s']); signs = np.asarray(e['ss'])
        for j, gid in enumerate(ids):
            value = 1/signs[j] if j % 4 < 2 and j // 4 < 2 else 0.
            require(gid not in one_seen or abs(one_seen[gid]-value) < 1e-14,
                    'constant scalar shared-map inconsistency')
            one_seen[gid] = value; one[gid] = value
        x, r, weight = volume_rule(e, data['H'], data['h'], nx, nz)
        b, gx, gr = basis(e, x, r, signs)
        owner = flow[e['geometry']]; ui = np.asarray(owner['u'])
        ub, ux, ur = basis(owner, x, r, np.asarray(owner['us']))
        q = b@target[ids]; qx = gx@target[ids]; qr = gr@target[ids]
        u, w = ub@velocity[ui], ub@velocity[nu+ui]
        divergence = ux@velocity[ui]+ur@velocity[nu+ui]
        exact_q = analytic_q(x, r); exact_qx, exact_qr = analytic_gradient(x, r)
        exact_u, exact_w = .2+.1*x+.03*r, .04+.02*r
        for name, actual, exact, scale in (
                ('q', q, exact_q, .01), ('qx', qx, exact_qx, .01),
                ('qr', qr, exact_qr, .01), ('u', u, exact_u, 1.),
                ('w', w, exact_w, 1.), ('divergence', divergence, .12, 1.)):
            pointwise[name] = max(pointwise[name], float(np.max(abs(actual-exact))/scale))
        cutoff = .03*(1+.03*x)
        diffusion = .003*(1+x+2*r)*np.maximum(r, cutoff)
        reaction = .03*(1+.1*x)
        forcing = 1+x*x+r**3+(1+.02*x)/np.maximum(r, cutoff)
        require(np.min(diffusion) > 0, 'independent diffusion is not positive')
        scale = np.full(len(x), data['Kin']) if case == 'full' else data['Ain']/r
        m[np.ix_(ids, ids)] += b.T@(weight[:, None]*b)
        riesz[np.ix_(ids, ids)] += b.T@((weight*scale*scale)[:, None]*b)
        k[np.ix_(ids, ids)] += (gx.T@((weight*diffusion)[:, None]*gx) +
                               gr.T@((weight*diffusion)[:, None]*gr) +
                               b.T@((weight*reaction)[:, None]*b))
        adv[ids] -= gx.T@(weight*exact_u*exact_q)+gr.T@(weight*exact_w*exact_q)
        transport[np.ix_(ids, ids)] -= gx.T@((weight*exact_u)[:, None]*b)+gr.T@((weight*exact_w)[:, None]*b)
        advective[ids] += b.T@(weight*(exact_u*exact_qx+exact_w*exact_qr))
        strong[ids] += b.T@(weight*(exact_u*exact_qx+exact_w*exact_qr+.12*exact_q))
        source[ids] += b.T@(weight*forcing)
        volume += float(weight.sum())
        work += float(weight@(diffusion*(exact_qx**2+exact_qr**2)+reaction*exact_q**2))
        source_work += float(weight@(exact_q*forcing))
    for e, side in faces:
        x, r, weight = face_rule(e, side, max(nx, nz))
        ids = np.asarray(e['s']); b, _, _ = basis(e, x, r, np.asarray(e['ss']))
        nx_face, nr_face = NORMALS[side]
        normal_velocity = nx_face*(.2+.1*x+.03*r)+nr_face*(.04+.02*r)
        q = analytic_q(x, r); outward_j = .002*(1+x+.5*r)
        boundary_adv[ids] += b.T@(weight*normal_velocity*q)
        transport[np.ix_(ids, ids)] += b.T@((weight*normal_velocity)[:, None]*b)
        flux[ids] -= b.T@(weight*outward_j)
        boundary_flow += float(weight@(normal_velocity*q))
        natural_work -= float(weight@(outward_j*q))
    adv += boundary_adv
    return dict(mass=m, K=k, riesz=riesz, transport=transport, advection=adv, advective_mutation=advective,
                conservative_strong=strong, boundary_advection=boundary_adv,
                source_load=source, natural_flux_rhs=flux, one=one, volume=volume,
                work=work, source_work=source_work, exterior_advective_flux=boundary_flow,
                natural_work=natural_work, pointwise_errors=pointwise)


def gamma_checks(data, case, faces, essential, free, order):
    ns = data['scalar_size']; gamma = vector(data, 'dependent_gamma_lift', ns)
    if case == 'full':
        require(la.norm(gamma) == 0., 'full fixture has an unexpected dependent gamma lift')
        require(data['corner_conflict_rejected'] is False,
                'compatible full single-side lift was rejected')
        return dict(applicable=False, compatible_single_side_not_rejected=True)
    require(data['corner_conflict_rejected'] is True,
            'incompatible inlet/interface corner was not rejected')
    require(np.max(abs(gamma[free]), initial=0.) == 0., 'gamma lift includes free coefficients')
    expected = np.zeros(ns); assigned = {}; bottom_error = endpoint_error = 0.
    trace_owners = []
    def value(x, r, side):
        return data['Ain']*(1+.02*x)**3/data['H'] if side == 2 else data['Ain']/r
    def assign(gid, coefficient):
        if gid in assigned:
            require(abs(assigned[gid]-coefficient) <= 1e-12*data['Ain']/data['H'],
                    'independent gamma corner is incompatible')
        else:
            assigned[gid] = coefficient; expected[gid] = coefficient
    for e, side in faces:
        if side not in (0, 2):
            continue
        ids = np.asarray(e['s']); signs = np.asarray(e['ss']); modes = EDGE_MODES[side]
        xl, xr, rl, rr = e['bounds']
        endpoints_x = np.array([xl, xl]) if side == 0 else np.array([xl, xr])
        endpoints_r = np.array([rl, rr]) if side == 0 else np.array([rl, rl])
        endpoint_basis, _, _ = basis(e, endpoints_x, endpoints_r, signs)
        endpoint_values = value(endpoints_x, endpoints_r, side)
        endpoint_error = max(endpoint_error, float(np.max(abs(
            endpoint_basis@gamma[ids]-endpoint_values))/(data['Ain']/data['H'])))
        coefficients = endpoint_values/signs[np.asarray(modes[:2])]
        for j in range(2):
            assign(ids[modes[j]], coefficients[j])
        x, r, weight = face_rule(e, side, order)
        b, _, _ = basis(e, x, r, signs)
        rhs_values = value(x, r, side)-b[:, modes[0]]*coefficients[0]-b[:, modes[1]]*coefficients[1]
        interior = b[:, modes[2:]]
        interior_coefficients = la.solve(interior.T@(weight[:, None]*interior),
                                         interior.T@(weight*rhs_values), assume_a='pos')
        for j in range(2):
            assign(ids[modes[j+2]], interior_coefficients[j])
        if side == 2:
            require(abs(rl-data['H']) < 1e-12 and rr > data['H'],
                    'dependent epsilon trace lacks outer H+ ownership')
            # An independent non-quadrature grid also checks the complete cubic.
            xx = np.linspace(xl, xr, 31); zz = np.full_like(xx, data['H'])
            bb, _, _ = basis(e, xx, zz, signs)
            bottom_error = max(bottom_error, float(np.max(abs(
                bb@gamma[ids]-value(xx, zz, side)))/(data['Ain']/data['H'])))
            trace_owners.append(dict(geometry=e['geometry'], bounds=e['bounds']))
    require(set(assigned) == set(essential), 'independent gamma mask coverage differs')
    coefficient_error = relative(gamma, expected)
    require(coefficient_error <= 1e-10 and endpoint_error <= 1e-10 and bottom_error <= 1e-10,
            'dependent gamma edge projection/endpoint/cubic H+ trace mismatch')
    return dict(applicable=True, corner_conflict_rejected=True,
                independent_endpoint_constrained_projection_error=coefficient_error,
                reference_scaled_endpoint_max_error=endpoint_error,
                reference_scaled_cubic_H_trace_max_error=bottom_error,
                H_plus_trace_owners=trace_owners)


def analyze(path, case, nx, nz, report):
    data = json.loads(path.read_text())
    require(data['case'] == case and data['foundation_only'] is True and
            data['full_coupled_gate_pass'] is False, 'wrong case/qualification scope')
    require(data['alpha'] == 2. and abs(data['d0']-.0875) < 1e-15 and
            abs(data['H']-.15) < 1e-15 and 0 < data['h'] < data['H'],
            'unexpected scalar fixture units/geometry/time coefficient')
    require(abs(data['Kin']-.008809015145515192) < 1e-16 and
            abs(data['Ain']-data['Kin']**1.5/2.3804029452632984) < 1e-16,
            'unexpected physical residual reference scales')
    flow, cells, faces, essential, free, area = inventory(data, case)
    ns, nu = data['scalar_size'], data['velocity_size']; nf = len(free)
    target = vector(data, 'target', ns); lift = vector(data, 'lift', ns)
    velocity = np.asarray(data['velocity'], dtype=float)
    require(velocity.ndim == 1 and len(velocity) >= 2*nu and np.all(np.isfinite(velocity)),
            'invalid canonical full velocity vector')
    checks = report['checks'] = {}
    ref = assemble(data, case, nx, nz, flow, cells, faces, target, velocity)
    require(abs(ref['volume']-area) < 1e-11*area, 'independent rule domain coverage differs')
    report['polynomial_reference_scaled_pointwise_errors'] = ref['pointwise_errors']
    require(max(ref['pointwise_errors'].values()) <= 1e-10,
            'native q/velocity/divergence differs from independent analytic polynomial')
    comparison(checks, 'one', vector(data, 'one', ns), ref['one'])
    for name in ('mass', 'K'):
        comparison(checks, name, matrix(data[name], ns), ref[name])
    if 'transport' in data:
        comparison(checks, 'transport_matrix', matrix(data['transport'], ns), ref['transport'])
        comparison(checks, 'transport_matrix_current_action', matrix(data['transport'], ns)@target, ref['advection'])
    for name in ('advection', 'source_load', 'natural_flux_rhs'):
        comparison(checks, name, vector(data, name, ns), ref[name])
    comparison(checks, 'conservative_weak_vs_strong', ref['advection'], ref['conservative_strong'])
    require(np.max(abs(lift[free]), initial=0.) == 0. and
            np.max(abs(lift[essential]-target[essential]), initial=0.) <= 1e-12,
            'essential lift does not use the independently derived scalar mask')
    full = data['alpha']*ref['mass']+ref['K']; system = full[np.ix_(free, free)]
    comparison(checks, 'system', matrix(data['system'], nf), system)
    load = full@target
    rhs = (load-full@lift)[free]
    comparison(checks, 'rhs', vector(data, 'rhs', nf), rhs)
    weighted_free = ref['riesz'][np.ix_(free, free)]
    factor = la.cho_factor(weighted_free, lower=True)
    def dual(residual):
        square = float(residual@la.cho_solve(factor, residual))
        require(np.isfinite(square) and square >= 0, 'invalid independent scalar dual norm')
        return float(np.sqrt(square/area))
    solves = {}
    for changed in (False, True):
        prefix = 'changed_' if changed else ''
        current_full = data['alpha']*ref['mass']+(1.2 if changed else 1.)*ref['K']
        current = current_full[np.ix_(free, free)]
        current_rhs = (load-current_full@lift)[free]
        comparison(checks, prefix+'rhs_current', vector(data, prefix+'rhs', nf), current_rhs)
        solution = vector(data, prefix+'solution', ns)
        require(np.max(abs(solution[essential]-lift[essential]), initial=0.) <= 1e-12,
                prefix+'solution changes essential scalar values')
        direct = la.solve(current, current_rhs, assume_a='pos')
        agreement = comparison(checks, prefix+'direct_solution', solution[free], direct, 1e-8)
        residual = dual(current@solution[free]-current_rhs)
        require(residual <= 1e-8, prefix+'true reference-scaled current residual failed')
        info = data[prefix+'solve']
        require(np.isfinite(info['residual']) and 0 <= info['residual'] <= 1e-8,
                prefix+'native residual report failed')
        require(info['iterations'] > 0 and info['rebuilds'] == (0 if changed else 1),
                prefix+'solve did not exercise intended retained-LU policy')
        require(all(np.isfinite(info[k]) and info[k] >= 0 for k in ('factor_seconds', 'solve_seconds')),
                prefix+'invalid native solver timing')
        if changed:
            require(info['factor_seconds'] == 0., '20% change unexpectedly rebuilt preconditioner')
        solves['changed' if changed else 'initial'] = dict(
            independent_current_dual_residual=residual, native_reported_residual=info['residual'],
            direct_solution_relative_error=agreement, iterations=info['iterations'],
            rebuilds=info['rebuilds'])
    solution = vector(data, 'solution', ns); changed_solution = vector(data, 'changed_solution', ns)
    comparison(checks, 'initial_target', solution, target, 1e-8)
    relative_change = relative(changed_solution, solution)
    require(relative_change > 1e-5, '20% changed-operator witness has no detectable scalar change')
    changed_full = data['alpha']*ref['mass']+1.2*ref['K']
    changed_rhs = (load-changed_full@lift)[free]
    frozen_direct = la.solve(system, changed_rhs, assume_a='pos')
    frozen_residual = dual(changed_full[np.ix_(free, free)]@frozen_direct-changed_rhs)
    require(frozen_residual > 1e-6, 'frozen old operator mutation is not detectable')
    one_adv = float(ref['one']@vector(data, 'advection', ns))
    work = float(target@matrix(data['K'], ns)@target)
    source_work = float(target@vector(data, 'source_load', ns))
    natural_work = float(target@vector(data, 'natural_flux_rhs', ns))
    balances = dict(
        constant_test_advection=dict(native=one_adv, exterior_integral=ref['exterior_advective_flux']),
        diffusion_reaction_work=dict(native=work, independent_integral=ref['work']),
        source_work=dict(native=source_work, independent_integral=ref['source_work']),
        outward_natural_work=dict(native=natural_work, independent_integral=ref['natural_work']))
    for name, values in balances.items():
        reference = values.get('independent_integral', values.get('exterior_integral'))
        values['relative_error'] = abs(values['native']-reference)/abs(reference)
        require(values['relative_error'] <= 1e-10, name+' sign/conservation/work identity failed')
    mutations = dict(
        missing_q_divergence_relative_error=relative(ref['advective_mutation'], ref['advection']),
        omitted_boundary_advection_relative_error=relative(
            ref['advection']-ref['boundary_advection'], ref['advection']),
        flipped_natural_flux_relative_error=relative(-ref['natural_flux_rhs'], ref['natural_flux_rhs']),
        frozen_original_operator_changed_current_dual_residual=frozen_residual)
    require(min(mutations['missing_q_divergence_relative_error'],
                mutations['omitted_boundary_advection_relative_error']) > 1e-3,
            'conservative-advection negative controls are ineffective')
    report.update(solves=solves, identities=balances, mutations=mutations,
                  changed_solution_relative_change=relative_change,
                  dependent_gamma=gamma_checks(data, case, faces, essential, free, max(nx, nz)),
                  geometry=dict(d0=data['d0'], H=data['H'], h=data['h'], volume=area,
                                scalar_geometries=[e['geometry'] for e in cells],
                                boundary_face_count=len(faces)),
                  native_dirichlet_count=data['native_dirichlet_count'],
                  essential_count=len(essential), free_count=nf, scalar_size=ns,
                  residual_norm='sqrt(r_free.T @ inverse(integral(reference_scale^2 phi_free phi_free.T)) @ r_free / volume)',
                  residual_reference_scale='Kin' if case == 'full' else 'Ain/r',
                  alpha_normalisation=False,
                  status='PASS_TINY_SCALAR_FOUNDATION_ONLY')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', required=True, type=Path)
    parser.add_argument('--case', required=True, choices=('full', 'outer'))
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--nx', type=int, default=20)
    parser.add_argument('--nz', type=int, default=24)
    args = parser.parse_args()
    report = dict(schema='WJH27IndependentScalarOracleV1', status='FAIL_SCALAR_FOUNDATION_ORACLE',
                  foundation_only=True, full_coupled_gate_pass=False, case=args.case,
                  input=str(args.native), quadrature=[args.nx, args.nz],
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  independent_basis_sha256=hashlib.sha256(FIXTURE_SOURCE.read_bytes()).hexdigest(),
                  platform=sys.platform,
                  scope='Independent tiny scalar algebra only; no evolving coupled MMS, physical pilot or Task3 qualification.')
    failure = None
    try:
        require(sys.platform.startswith('linux'), 'authoritative scalar verification requires Linux')
        require(args.nx >= 20 and args.nz >= 24, 'independent quadrature must be at least Gauss20x24')
        report['input_sha256'] = hashlib.sha256(args.native.read_bytes()).hexdigest()
        analyze(args.native, args.case, args.nx, args.nz, report)
    except Exception as error:
        failure = error
        report.update(exception=type(error).__name__, reason=str(error))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if failure is not None:
        raise failure
    print(json.dumps(dict(status=report['status'], case=args.case,
                          max_comparison_error=max(v['relative_error'] for v in report['checks'].values()),
                          initial_dual=report['solves']['initial']['independent_current_dual_residual'],
                          changed_dual=report['solves']['changed']['independent_current_dual_residual']), sort_keys=True))


if __name__ == '__main__':
    main()
