#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent tiny mixed-BC/Q3-scalar-force assembly oracle.

Uses the H26 independent Jacobi definition, not native sampled basis/operators.
Default integration is physical Gauss20x20, split at roughness and H with outer
trace ownership. This is a foundation gate only, never a full coupled-model PASS.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'h26'))
import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
from fixtures import modified


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sparse(packet):
    result = np.zeros((packet['rows'], packet['cols']))
    for i, j, value in packet['entries']:
        result[i, j] += value
    require(np.all(np.isfinite(result)), 'nonfinite native matrix')
    return result


def basis(degree, xi, eta, dx, dz, signs):
    bx, gx = modified(degree, xi)
    bz, gz = modified(degree, eta)
    value = (bz[:, :, None]*bx[:, None, :]).reshape(len(xi), -1)*signs
    ddx = (bz[:, :, None]*gx[:, None, :]).reshape(len(xi), -1)*signs*2/dx
    ddz = (gz[:, :, None]*bx[:, None, :]).reshape(len(xi), -1)*signs*2/dz
    return value, ddx, ddz


def norm_error(actual, reference):
    scale = la.norm(reference)
    require(scale > 0, 'relative comparison requires a nonzero declared reference')
    return float(la.norm(actual-reference)/scale)


def expected_masks(data, case):
    elements = data['elements']; n = data['velocity_size']
    xmin = min(v['bounds'][0] for v in elements); xmax = max(v['bounds'][1] for v in elements)
    zmin = min(v['bounds'][2] for v in elements); zmax = max(v['bounds'][3] for v in elements)
    sides = [set(), set(), set(), set()]
    for e in elements:
        xl, xr, zl, zr = e['bounds']
        for j, gid in enumerate(e['u']):
            ix, iz = j % 4, j // 4
            for condition, index in [(abs(xl-xmin) < 1e-12 and ix == 0, 0),
                                     (abs(xr-xmax) < 1e-12 and ix == 1, 1),
                                     (abs(zl-zmin) < 1e-12 and iz == 0, 2),
                                     (abs(zr-zmax) < 1e-12 and iz == 1, 3)]:
                if condition:
                    sides[index].add(gid)
    u_sides = (0, 1, 2, 3) if case == 'closed' else (0, 1, 2) if case == 'slip' else (0, 2)
    v_sides = (0, 1, 2, 3) if case != 'physicaloutlet' else (0, 2, 3)
    essential = [set().union(*(sides[j] for j in indices)) for indices in (u_sides, v_sides)]
    declared = data['free_velocity_ids']
    require(len(declared) == 2, 'free_velocity_ids must contain separate u/v canonical lists')
    for component in range(2):
        ids = declared[component]
        require(len(set(ids)) == len(ids), 'duplicate free velocity ids')
        require(set(ids) == set(range(n))-essential[component],
                'component boundary/free space differs from independent modal-edge mask')
    require(bool(data['gauge']) == (case != 'physicaloutlet'), 'physical boundary/gauge mismatch')
    return essential, [np.array(v, dtype=int) for v in declared]


def verify_native_v(data):
    n = data['velocity_size']; target = np.array(data['target'])
    mapping = np.array(data['v_to_u'], dtype=int); signs = np.array(data['v_to_u_sign'], dtype=float)
    require(len(mapping) == n and sorted(mapping.tolist()) == list(range(n)), 'native v map is not bijective')
    require(len(signs) == n and np.all(np.isin(signs, [-1., 1.])), 'invalid native v signs')
    for e in data['elements']:
        for j, native in enumerate(e['v_native']):
            require(mapping[native] == e['u'][j], 'native v local id does not map to canonical u id')
            require(abs(e['vs_native'][j]-e['us'][j]*signs[native]) < 1e-14,
                    'native v local sign does not map to canonical sign')
    native_target = np.array(data.get('nativeVcoefftarget', data.get('native_v_coefftarget')), dtype=float)
    require(native_target.shape == (n,), 'missing nativeVcoefftarget round-trip witness')
    canonical = np.zeros(n); canonical[mapping] = signs*native_target
    error = norm_error(canonical, target[n:2*n])
    require(error <= 1e-12, 'native v target round trip failed')
    return dict(relative_target_error=error,
                permuted_ids=int(np.sum(mapping != np.arange(n))),
                negative_signs=int(np.sum(signs < 0)))


def assemble(data, order):
    n, p = data['velocity_size'], data['pressure_size']
    nk = data['known_scalar_action']['cols']
    require(data['known_scalar_action']['rows'] == 2*n, 'scalar action shape must be 2Nu x Nk')
    h = float(data.get('h', .10854503542535461)); height = float(data.get('H', .15))
    require(float(data.get('curve', 0.)) == 0., 'this oracle supports the declared flat foundation only')
    m = np.zeros((n, n)); qmass = np.zeros((p, p)); kmass = np.zeros((nk, nk))
    cross = np.zeros((p, nk)); div = np.zeros((p, 2*n)); strain = np.zeros((2*n, 2*n))
    component = np.zeros_like(strain); cp = np.zeros((2*n, p)); ck = np.zeros((2*n, nk))
    mean = np.zeros(p); energy = 0.; volume = 0.; target_errors = []
    z, w = leggauss(order); xx, zz = np.meshgrid(z, z)
    xi, eta = xx.ravel(), zz.ravel(); ww = np.outer(w, w).ravel()
    target = np.array(data['target']); scalar = np.array(data['known_scalar_coefficients'])
    require(target.shape == (2*n+p,) and scalar.shape == (nk,), 'target/scalar coefficient shape mismatch')
    trace_owners = []; trace_cubic_coefficients = []
    for e in data['elements']:
        xl, xr, zl, zr = e['bounds']; dx, dz = xr-xl, zr-zl
        require(dx > 0 and dz > 0, 'nonpositive physical element dimensions')
        x = (xl+xr)/2+xi*dx/2; physical_z = (zl+zr)/2+eta*dz/2
        b, gx, gz = basis(3, xi, eta, dx, dz, np.array(e['us']))
        pressure, _, _ = basis(2, xi, eta, dx, dz, np.array(e['ps']))
        scalar_basis, _, _ = basis(3, xi, eta, dx, dz, np.array(e['ks']))
        ui, pi, ki = np.array(e['u']), np.array(e['p']), np.array(e['k'])
        weight = ww*dx*dz/4; weighted_mu = weight*.02*(1+x+physical_z)
        volume += dx*dz
        m[np.ix_(ui, ui)] += b.T@(weight[:, None]*b)
        qmass[np.ix_(pi, pi)] += pressure.T@(weight[:, None]*pressure)
        kmass[np.ix_(ki, ki)] += scalar_basis.T@(weight[:, None]*scalar_basis)
        cross[np.ix_(pi, ki)] += pressure.T@(weight[:, None]*scalar_basis)
        mean[pi] += pressure.T@weight
        for comp, g in enumerate((gx, gz)):
            div[np.ix_(pi, comp*n+ui)] += pressure.T@(weight[:, None]*g)
        grad = (gx, gz)
        lap = gx.T@(weighted_mu[:, None]*gx)+gz.T@(weighted_mu[:, None]*gz)
        for a in range(2):
            for c in range(2):
                block = grad[c].T@(weighted_mu[:, None]*grad[a])
                if a == c:
                    block += lap
                    component[np.ix_(a*n+ui, c*n+ui)] += lap
                strain[np.ix_(a*n+ui, c*n+ui)] += block
        ux, uz = gx@target[ui], gz@target[ui]
        vx, vz = gx@target[n+ui], gz@target[n+ui]
        energy += float(weighted_mu@(2*ux**2+2*vz**2+(uz+vx)**2))
        exact = [x*x*(1-x)*(2*physical_z-3*physical_z**2),
                 -(2*x-3*x*x)*physical_z**2*(1-physical_z),
                 .01*(x*x+physical_z**2-2/3)]
        represented = [b@target[ui], b@target[n+ui], pressure@target[2*n+pi]]
        target_errors.append([float(np.max(abs(a-b))) for a, b in zip(represented, exact)])
        lo, hi = max(zl, 0.), min(zr, h)
        if hi <= lo:
            continue
        candidates = [v for v in data['elements']
                      if abs(v['bounds'][0]-xl) < 1e-12 and abs(v['bounds'][1]-xr) < 1e-12
                      and v['bounds'][2] <= height+1e-12 and v['bounds'][3] > height+1e-12]
        require(len(candidates) == 1, 'H+ trace element owner absent or ambiguous')
        tr = candidates[0]; txl, txr, tzl, tzr = tr['bounds']
        require(tzl >= height-1e-12 or tzl < height < tzr, 'wrong H trace side')
        eta_h = np.full_like(xi, 2*(height-tzl)/(tzr-tzl)-1)
        _, pressure_x, _ = basis(2, xi, eta_h, dx, tzr-tzl, np.array(tr['ps']))
        scalar_trace, scalar_x, _ = basis(3, xi, eta_h, dx, tzr-tzl, np.array(tr['ks']))
        zband = (lo+hi)/2+eta*(hi-lo)/2
        test, _, _ = basis(3, xi, 2*(zband-zl)/dz-1, dx, dz, np.array(e['us']))
        weight_band = ww*dx*(hi-lo)/4
        cp[np.ix_(ui, np.array(tr['p']))] += test.T@(weight_band[:, None]*pressure_x)
        ck[np.ix_(ui, np.array(tr['k']))] += test.T@(weight_band[:, None]*scalar_x)
        trace_owners.append(dict(test_geometry=e['geometry'], trace_geometry=tr['geometry'],
                                 x_interval=[xl, xr], trace_bounds=tr['bounds'], height=height))
        # Legendre coefficient is independent of native scalar mode ordering.
        leg3 = (5*xi**3-3*xi)/2
        trace_value = scalar_trace@scalar[np.array(tr['k'])]
        cubic = float(7/4*np.sum(ww*trace_value*leg3))
        trace_cubic_coefficients.append(dict(geometry=tr['geometry'], coefficient=cubic))
    return dict(mass=m, pressure_mass=qmass, scalar_mass=kmass, pressure_scalar_mass=cross,
                D=div, G=-div.T, C=cp, K=strain, K_component=component,
                known_scalar_action=ck, mean=mean, energy=energy, volume=volume,
                target_pointwise_errors=np.max(target_errors, axis=0).tolist(),
                trace_owners=trace_owners, trace_cubic_coefficients=trace_cubic_coefficients)


def analyze(path, case, order):
    data = json.loads(path.read_text()); n, p = data['velocity_size'], data['pressure_size']
    essential, free = expected_masks(data, case)
    mapping = verify_native_v(data)
    if case != 'closed':
        require(mapping['permuted_ids']+mapping['negative_signs'] > 0, 'mixed fixture lacks a genuinely different native v map')
    ref = assemble(data, order); checks = {}
    for name in ('mass', 'pressure_mass', 'D', 'G', 'C', 'K', 'known_scalar_action'):
        error = norm_error(sparse(data[name]), ref[name])
        checks[name] = dict(relative_error=error)
        require(error <= 1e-10, name+' independent operator mismatch')
    require(max(ref['target_pointwise_errors']) <= 1e-11, 'native polynomial target differs from independent analytic field')
    require(norm_error(np.array(data['mean']), ref['mean']) <= 1e-11, 'pressure mean action mismatch')
    one = np.zeros(p); seen = {}
    for e in data['elements']:
        for j, gid in enumerate(e['p']):
            value = (1/e['ps'][j]) if j % 3 < 2 and j // 3 < 2 else 0.
            require(gid not in seen or abs(seen[gid]-value) < 1e-14, 'constant pressure shared sign conflict')
            one[gid] = value; seen[gid] = value
    require(norm_error(np.array(data['pressure_one']), one) <= 1e-12, 'native constant pressure representation mismatch')
    keep = np.r_[free[0], n+free[1]]; nv = len(keep)
    a = ref['K'].copy(); alpha, awf = float(data['alpha']), float(data['awf'])
    require(alpha == 2. and awf == 1., 'unexpected frozen tiny fixture alpha/AWF')
    a[:n, :n] += alpha*ref['mass']; a[n:, n:] += alpha*ref['mass']
    gp = ref['G']-awf*ref['C']; gauge = bool(data['gauge'])
    system = np.zeros((nv+p+int(gauge), nv+p+int(gauge)))
    system[:nv, :nv] = a[np.ix_(keep, keep)]
    system[:nv, nv:nv+p] = gp[keep]
    system[nv:nv+p, :nv] = ref['D'][:, keep]
    if gauge:
        system[nv:nv+p, -1] = ref['mean']; system[-1, nv:nv+p] = ref['mean']
    error = norm_error(sparse(data['system']), system)
    checks['system'] = dict(relative_error=error)
    require(error <= 1e-10, 'mixed free-system mismatch')
    constant_c = float(la.norm(ref['C']@one))
    constant_free = float(la.norm(gp[keep]@one))
    require(constant_c <= 1e-12, 'AWF gradient acts on constant pressure')
    require(constant_free <= 1e-11 if gauge else constant_free > 1e-6,
            'restricted constant-pressure action contradicts physical gauge')
    unaugmented = system[:nv+p, :nv+p]
    singular = la.svdvals(unaugmented)
    null_count = int(np.sum(singular <= singular[0]*1e-12))
    require(null_count == int(gauge), 'wrong unaugmented pressure null count')
    target, lift, solution = (np.array(data[name]) for name in ('target', 'lift', 'solution'))
    for component in range(2):
        ids = np.array(sorted(essential[component]), dtype=int)+component*n
        require(np.max(abs(lift[ids]-target[ids]), initial=0.) <= 1e-12, 'essential lift does not match target')
        require(np.max(abs(solution[ids]-lift[ids]), initial=0.) <= 1e-12, 'solution changed an essential coefficient')
    require(np.max(abs(lift[keep]), initial=0.) <= 1e-14 and la.norm(lift[2*n:]) <= 1e-14,
            'lift includes unsupported free velocity or pressure values')
    target_free = np.r_[target[keep], target[2*n:]]
    if gauge:
        target_free = np.r_[target_free, 0.]
    rhs = np.array(data['rhs'])
    checks['rhs'] = dict(relative_error=norm_error(rhs, system@target_free))
    require(checks['rhs']['relative_error'] <= 1e-10, 'algebraic target fixture RHS mismatch')
    direct = la.solve(system, rhs)
    actual = np.r_[solution[keep], solution[2*n:]]
    if gauge:
        actual = np.r_[actual, float(data.get('gauge_multiplier', 0.))]
    agreement = norm_error(actual, direct)
    require(agreement <= 1e-8, 'independent direct solution mismatch')
    residual = system@actual-rhs
    velocity_mass = la.block_diag(ref['mass'][np.ix_(free[0], free[0])], ref['mass'][np.ix_(free[1], free[1])])
    momentum = math.sqrt(max(0., residual[:nv]@la.solve(velocity_mass, residual[:nv], assume_a='pos'))/ref['volume'])
    continuity_vector = ref['D']@solution[:2*n]
    continuity = math.sqrt(max(0., continuity_vector@la.solve(ref['pressure_mass'], continuity_vector, assume_a='pos'))/ref['volume'])
    require(max(momentum, continuity) <= 1e-9, 'true dual residual/weak continuity failed')
    work = float(solution[:2*n]@ref['K']@solution[:2*n])
    target_work = float(target[:2*n]@sparse(data['K'])@target[:2*n])
    work_error = abs(target_work-ref['energy'])/ref['energy']
    require(work_error <= 1e-10, 'full symmetric-strain work identity mismatch')
    scalar = np.array(data['known_scalar_coefficients'])
    known = (-2/3)*ref['known_scalar_action']@scalar
    native_known = np.array(data['known_scalar_load'])
    require(native_known.shape == (2*n,), 'known scalar load must have 2Nu entries')
    known_error = norm_error(native_known, known)
    require(known_error <= 1e-10, 'actual Q3 mechanical-pressure force/sign mismatch')
    projected = la.solve(ref['pressure_mass'], ref['pressure_scalar_mass']@scalar, assume_a='pos')
    projected_load = (-2/3)*ref['C']@projected
    projection_loss = norm_error(projected_load, known)
    old_loss = norm_error(np.array(data['oldQ2k_load']), known)
    max_cubic = max(abs(v['coefficient']) for v in ref['trace_cubic_coefficients'])
    require(max_cubic > 1e-8, 'fixture lacks an independently resolved cubic scalar trace')
    require(projection_loss > 1e-4 and old_loss > 1e-4, 'Q2 replacement negative control is ineffective')
    missing_transpose = norm_error(ref['K_component'], ref['K'])
    flipped_awf = norm_error((ref['G']+awf*ref['C'])[keep], gp[keep])
    flipped_scalar_sign = norm_error(-known, known)
    require(min(missing_transpose, flipped_awf) > 1e-8, 'operator mutations not detectable')
    return dict(schema='WJH27IndependentMixedQ3OracleV1',
                status='PASS_TINY_MIXED_BC_Q3_FOUNDATION_ONLY', full_coupled_gate_pass=False,
                case=case, quadrature=[order, order], input=str(path),
                input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                checks=checks, native_v_map=mapping,
                component_free_counts=[len(v) for v in free], gauge=gauge,
                unaugmented_null_modes=null_count, smallest_singular_values=singular[-5:].tolist(),
                constant_pressure_C_action=constant_c, restricted_constant_pressure_action=constant_free,
                direct_solution_relative_error=agreement, momentum_dual_residual=momentum,
                weak_continuity_dual_norm=continuity,
                norm_scales='Fixed Uref=Lref=force_ref=1 on tiny dimensionless fixture; no alpha/rhs normalization.',
                polynomial_target_max_errors=ref['target_pointwise_errors'],
                full_strain_work=dict(target_discrete=target_work, independent_integral=ref['energy'],
                                      relative_error=work_error, solution_work=work),
                known_scalar=dict(relative_force_error=known_error, independent_Q2_L2_projection_loss=projection_loss,
                                  native_oldQ2_negative_control_loss=old_loss,
                                  maximum_cubic_trace_Legendre_coefficient=max_cubic,
                                  outer_trace_owners=ref['trace_owners']),
                mutations=dict(missing_transpose=missing_transpose, flipped_awf=flipped_awf,
                               flipped_scalar_sign=flipped_scalar_sign),
                scope='Independent tiny algebraic/direct foundation only; no continuous MMS, evolving scalars, physical pilot or Task3 qualification.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', required=True, type=Path)
    parser.add_argument('--case', required=True, choices=['closed', 'slip', 'physicaloutlet'])
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--order', type=int, default=20)
    args = parser.parse_args()
    require(args.order >= 20, 'independent quadrature must be at least Gauss20x20')
    args.out.parent.mkdir(parents=True, exist_ok=True)
    try:
        result = analyze(args.native, args.case, args.order)
    except Exception as error:
        result = dict(schema='WJH27IndependentMixedQ3OracleV1', status='FAIL_FOUNDATION_ORACLE',
                      case=args.case, input=str(args.native), exception=type(error).__name__, reason=str(error),
                      full_coupled_gate_pass=False)
        with args.out.open('x') as stream:
            stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
        raise
    with args.out.open('x') as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=result['status'], case=args.case,
                          max_operator_error=max(v['relative_error'] for v in result['checks'].values()),
                          direct_error=result['direct_solution_relative_error'],
                          Q2_projection_loss=result['known_scalar']['independent_Q2_L2_projection_loss']), sort_keys=True))


if __name__ == '__main__':
    main()
