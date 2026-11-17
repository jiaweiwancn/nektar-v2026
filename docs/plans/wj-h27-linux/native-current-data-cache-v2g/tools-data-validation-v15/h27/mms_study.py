#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent Linux-only completed-endpoint continuous MMS study.

No solver launches. All files are existing completed native endpoints; output
is create-only. Sampled MAX, physical owner-gradient accounts and algebraic
essential reactions are distinct evidence, never interchangeable certificates.
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
import numpy as np
import scipy.linalg as la
from numpy.polynomial.legendre import leggauss
import coupled_mms_check
import coupled_oracle as core

FIELDS = ('U', 'k', 'epsilon', 'Pi', 'gradU', 'gradk', 'gradepsilon', 'gradPi')
TEMPORAL_FIELDS = ('U', 'k', 'epsilon')
RULES = ((20, 24), (24, 32))
GRIDS = (17, 33)
TIME = .4
ORDER_MINIMUM = 1.8
REFERENCE_ROUNDOFF = 1e-12
MAX_REFINEMENT_FRACTION = .01
ERROR_QUADRATURE_FRACTION = .001
NONTEMPORAL_FRACTION = .1
OWNER_RESOLVED_BAND = 1e-10


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clock(deadline):
    if time.monotonic() >= deadline:
        raise TimeoutError('MMS study reached its explicitly admitted analysis deadline')


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def mesh_identity(data):
    return digest(sorted(e['bounds'] for e in data['elements']))


def fields(case, e, state, x, r):
    """Actual signed native fields and one-element-owner physical gradients."""
    x, r = np.asarray(x, float), np.asarray(r, float)
    b, dx, dr = core.basis(e, 'u', x, r); ids = np.asarray(e['u'])
    n = case['data']['sizes'][0]
    U = np.vstack((b@state['flow'][ids], b@state['flow'][n+ids]))
    G = np.stack((dx@state['flow'][ids], dr@state['flow'][ids],
                  dx@state['flow'][n+ids], dr@state['flow'][n+ids]))
    kb, kx, kr = core.basis(e, 'k', x, r); kc = state['k'][e['k']]
    pb, px, pr = core.basis(e, 'p', x, r); pc = state['flow'][2*n+np.asarray(e['p'])]
    result = dict(U=U, gradU=G, k=kb@kc, gradk=np.vstack((kx@kc, kr@kc)),
                  Pi=pb@pc, gradPi=np.vstack((px@pc, pr@pc)))
    if e['epsilon']:
        eb, ex, er = core.basis(e, 'epsilon', x, r); ec = state['epsilon'][e['epsilon']]
        result.update(epsilon=eb@ec, gradepsilon=np.vstack((ex@ec, er@ec)))
    return result


def exact(case, x, r, outer):
    m = case['m']; when = case['data']['time']
    j = {key: np.vstack([np.broadcast_to(np.asarray(v, float), x.shape) for v in f(x, r, when)])
         for key, f in m.functions.items() if key != 'epsilon' or outer}
    result = dict(U=np.vstack((j['u'][0], j['w'][0])),
        gradU=np.vstack((j['u'][1], j['u'][2], j['w'][1], j['w'][2])),
        k=j['k'][0], gradk=j['k'][1:3], Pi=j['Pi'][0], gradPi=j['Pi'][1:3])
    if outer: result.update(epsilon=j['epsilon'][0], gradepsilon=j['epsilon'][1:3])
    return result


def scaled_differences(m, first, second, r, physical=False):
    result = {}
    for key in FIELDS:
        if key not in first or key not in second: continue
        delta = first[key]-second[key]
        value = np.sqrt(np.sum(delta*delta, axis=0)) if delta.ndim == 2 else abs(delta)
        if not physical:
            if key in ('k', 'gradk'): value = value/m.Kin
            elif key in ('epsilon', 'gradepsilon'): value = value*r/m.Ain
            elif key in ('Pi', 'gradPi'): value = value/m.tau0
        result[key] = value
    return result


def rectangle_rule(bounds, nx, nr, grid=False):
    xl, xr, rl, rr = bounds
    if grid:
        x, r = np.meshgrid(np.linspace(xl, xr, nx), np.linspace(rl, rr, nr))
        return x.ravel(), r.ravel(), None
    a, wa = leggauss(nx); b, wb = leggauss(nr)
    x, r = np.meshgrid((xl+xr)/2+(xr-xl)*a/2, (rl+rr)/2+(rr-rl)*b/2)
    w = (wb[:, None]*wa[None, :]*(xr-xl)*(rr-rl)/4).ravel()
    return x.ravel(), r.ravel(), w


def area(m, key):
    return m.R-m.H if key in ('epsilon', 'gradepsilon') else m.R


def endpoint_errors(case, o, grid_size):
    energy = {key: 0. for key in FIELDS}; maxima = {key: 0. for key in FIELDS}
    physical_energy = energy.copy(); physical_maxima = maxima.copy()
    locations = {}
    for e in case['data']['elements']:
        clock(o.deadline)
        x, r, w = o.rule(e, case['state'], case['data']['time'])
        current, target = fields(case, e, case['state'], x, r), exact(case, x, r, bool(e['epsilon']))
        error = scaled_differences(case['m'], current, target, r)
        for key, values in error.items(): energy[key] += float(w@(values*values))
        for key, values in scaled_differences(case['m'], current, target, r, True).items():
            physical_energy[key] += float(w@(values*values))
        gx, gr, _ = rectangle_rule(e['bounds'], grid_size, grid_size, True)
        # Include independent Gauss nodes as well as endpoints/owner grids.
        x, r = np.r_[x, gx], np.r_[r, gr]
        current, target = fields(case, e, case['state'], x, r), exact(case, x, r, bool(e['epsilon']))
        error = scaled_differences(case['m'], current, target, r)
        for key, values in scaled_differences(case['m'], current, target, r, True).items():
            physical_maxima[key] = max(physical_maxima[key], float(np.max(values)))
        for key, values in error.items():
            index = int(np.argmax(values))
            if key not in locations or values[index] > maxima[key]:
                maxima[key] = float(values[index]); locations[key] = dict(
                    x=float(x[index]), r=float(r[index]), z=float(r[index]+case['m'].d0),
                    owner_geometry=e['geometry'])
    return {key: dict(RMS=math.sqrt(energy[key]/area(case['m'], key)),
                      sampled_MAX=maxima[key], sampled_MAX_location=locations[key],
                      physical_RMS=math.sqrt(physical_energy[key]/area(case['m'], key)),
                      physical_sampled_MAX=physical_maxima[key]) for key in FIELDS}


def face_rule(o, e, side, state, order, when, include_exact=True):
    a, w = leggauss(order); xl, xr, rl, rr = e['bounds']
    if side >= 2:
        return (xl+xr)/2+(xr-xl)*a/2, np.full(order, rl if side == 2 else rr, dtype=float), w*(xr-xl)/2
    xf = xl if side == 0 else xr
    K = o.trace(e, 'k', np.array([xf], dtype=float), state)[0]
    cuts = [rl, rr, float(o.cutoff(K)), .4*o.m.H, .75*o.m.H, o.m.h, o.m.H]
    if include_exact:
        cuts.append(float(o.cutoff(o.m.functions['K'](xf, o.m.H, when)[0])))
    tolerance = 64*np.finfo(float).eps*max(1., abs(rl), abs(rr))
    cuts = sorted(set([rl, rr]+[v for v in cuts if rl+tolerance < v < rr-tolerance]))
    radii, weights = [], []
    for lo, hi in zip(cuts[:-1], cuts[1:]):
        radii.extend((lo+hi)/2+(hi-lo)*a/2); weights.extend(w*(hi-lo)/2)
    return np.full(len(radii), xf, dtype=float), np.asarray(radii, float), np.asarray(weights, float)


def exact_scalar(o, x, r, outer, when):
    j = o.exact(x, r, when, outer); k, K = j['k'][0], j['K'][0]
    E = j['epsilon'][0] if outer else K**1.5/o.m.cl/np.maximum(r, o.cutoff(K))
    require(np.all(np.isfinite(E)) and np.min(E) > 0, 'invalid exact prescribed scalar destruction')
    return dict(j=j, Dk=o.m.nu+o.m.cd*k*k/(o.m.sk*E),
                De=o.m.nu+o.m.cd*k*k/(o.m.se*E) if outer else None)


def assemble_current(case, nx, nr, deadline):
    data, m, s, h = case['data'], case['m'], case['state'], case['history']
    o = core.Oracle(data, m, nx, nr, deadline)
    when = data['time']
    # Local study override only: all exterior constitutive integrations resolve
    # the current and analytic cutoffs. The imported frozen oracle is unchanged.
    o.face = lambda e, side: face_rule(o, e, side, s, max(nx, nr), when)
    prior, older = case['prior'], case['older']; dt = data['dt']
    a1, conv1, rough1 = o.advection(prior, h['time'])
    a2, conv2, rough2 = o.advection(older, h['time']-dt)
    explicit = core.combine(a1, 2., a2, -1.)
    history = core.combine(prior, 2/dt, older, -.5/dt)
    forcing = o.forcing(s, when)
    base = core.combine(core.combine(forcing, 1., o.mass(history), 1.), 1., explicit, -1.)
    terms = o.current(s, when)
    residual, norms = o.residual(s, when, 1.5/dt, base, terms)
    errors = {name: o.dual(core.combine(case[name], 1., actual, -1.)) for name, actual in
              (('stored_advection', a1), ('stored_previous_advection', a2))}
    current_advection_error = None
    if 'current_advection' in data:
        actual = core.state(data['current_advection'], data['sizes'])
        require(core.state_hash(actual) == data['advection_sha256'], 'endpoint current advection hash mismatch')
        current_advection_error = o.dual(core.combine(actual, 1., o.advection(s, when)[0], -1.))
    trace = {name: o.trace_error(value, when) for name, value, when in
             (('endpoint', s, when), ('prior', prior, h['time']), ('older', older, h['time']-dt))}
    return dict(o=o, forcing=forcing, terms=terms, residual=residual, norms=norms,
                explicit=explicit, history_errors=errors, current_advection_error=current_advection_error,
                trace=trace, pressure_mean=o.mean.copy())


def budgets(case, assembled, order):
    o = assembled['o']; m = case['m']; s, prior, older = case['state'], case['prior'], case['older']
    dt = case['data']['dt']; when = case['data']['time']; volume = o.volume
    scalar = {key: dict(time_BDF2=0., conservative_advection_EX2=0., production=0., destruction=0.,
                       owner_outward_diffusive_flux=0., exact_prescribed_natural_flux=0.) for key in ('k', 'epsilon')}
    mechanical = dict(velocity_time_BDF2_work=0., advection_EX2_work=0., stress_volume_work=0.,
        pressure_volume_work=0., owner_complete_traction_power=0., sampled_AWF_convection_work=0.,
        sampled_AWF_pressure_work=0., sampled_AWF_known_k_work=0., sampled_AWF_body_work=0., bulk_body_work=0.)
    kinetic = dict(current=0., prior=0., older=0., current_extrapolation=0., prior_extrapolation=0., second_difference=0.)
    gamma_error_integral = 0.; gamma_error_max = 0.
    diffusion_constant = dict(k=0., epsilon=0.)
    for e in o.elements:
        clock(o.deadline); x, r, w = o.rule(e, s, when)
        value, v1, v2 = [fields(case, e, state, x, r) for state in (s, prior, older)]
        a = o.sampled(e, x, r, s)
        U, U1, U2 = value['U'], v1['U'], v2['U']
        derivative = (1.5*U-2*U1+.5*U2)/dt
        mechanical['velocity_time_BDF2_work'] += float(w@np.sum(U*derivative, axis=0))
        for key, vv in (('current', U), ('prior', U1), ('older', U2),
                        ('current_extrapolation', 2*U-U1), ('prior_extrapolation', 2*U1-U2),
                        ('second_difference', U-2*U1+U2)):
            kinetic[key] += float(w@np.sum(vv*vv, axis=0))
        conv = []
        for vv in (v1, v2):
            g = vv['gradU']; conv.append(np.vstack((vv['U'][0]*g[0]+vv['U'][1]*g[1],
                                                     vv['U'][0]*g[2]+vv['U'][1]*g[3])))
        mechanical['advection_EX2_work'] += float(w@np.sum(U*(2*conv[0]-conv[1]), axis=0))
        mechanical['stress_volume_work'] += float(w@np.einsum('qij,qij->q', a['G'], a['V']))
        mechanical['pressure_volume_work'] -= float(w@(value['Pi']*(value['gradU'][0]+value['gradU'][3])))
        b = m.body(x, r, when)
        mechanical['bulk_body_work'] += float(w@(U[0]*b[0]+U[1]*b[2]))
        mask = r < m.h
        if np.any(mask):
            xx = x[mask]; owner = o.owner[e['geometry']]; hh = np.full(len(xx), m.H, dtype=float)
            at1, at2 = [fields(case, owner, state, xx, hh) for state in (prior, older)]
            conv_H = 2*(2*at1['U'][0]*at1['gradU'][0])-2*at2['U'][0]*at2['gradU'][0]
            work = w[mask]*U[0, mask]
            mechanical['sampled_AWF_convection_work'] += float(work@conv_H)
            mechanical['sampled_AWF_pressure_work'] += float(work@o.trace(e, 'p', xx, s, True))
            mechanical['sampled_AWF_known_k_work'] += float(work@(-(2./3)*o.trace(e, 'k', xx, s, True)))
            mechanical['sampled_AWF_body_work'] += float(work@(-m.body(xx, hh, when)[0]))
        for key in ('k', 'epsilon'):
            if key not in value: continue
            _, bx, br = core.basis(e, key, x, r)
            one = o.one[key][e[key]]; grad = value['gradk' if key == 'k' else 'gradepsilon']
            diffusion_constant[key] += float(w@(a['Dk' if key == 'k' else 'De']*
                ((bx@one)*grad[0]+(br@one)*grad[1])))
            scalar[key]['time_BDF2'] += float(w@((1.5*value[key]-2*v1[key]+.5*v2[key])/dt))
            if key == 'k':
                scalar[key]['production'] += float(w@a['production'])
                scalar[key]['destruction'] += float(w@a['E'])
            else:
                scalar[key]['production'] += float(w@(1.44*a['production']*a['E']/a['k']))
                scalar[key]['destruction'] += float(w@(1.82*a['E']*a['E']/a['k']))
    face_rows = []
    for key in ('k', 'epsilon'):
        for e, side in o.faces[key]:
            clock(o.deadline); x, r, w = face_rule(o, e, side, s, order, when)
            value, v1, v2 = [fields(case, e, state, x, r) for state in (s, prior, older)]
            a = o.sampled(e, x, r, s); nx, nz = core.NORMALS[side]
            grad = value['gradk' if key == 'k' else 'gradepsilon']
            D = a['Dk' if key == 'k' else 'De']
            current_flux = -D*(nx*grad[0]+nz*grad[1])
            adv1 = (nx*v1['U'][0]+nz*v1['U'][1])*v1[key]
            adv2 = (nx*v2['U'][0]+nz*v2['U'][1])*v2[key]
            adv = float(w@(2*adv1-adv2)); physical = float(w@current_flux)
            scalar[key]['conservative_advection_EX2'] += adv
            scalar[key]['owner_outward_diffusive_flux'] += physical
            prescribed = 0.; natural_callback = key == 'k' or side != 2
            if natural_callback:
                analytical = exact_scalar(o, x, r, bool(e['epsilon']), when)
                exact_gradient = analytical['j'][key][1:3]
                exact_D = analytical['Dk' if key == 'k' else 'De']
                prescribed = float(w@(-exact_D*(nx*exact_gradient[0]+nz*exact_gradient[1])))
                scalar[key]['exact_prescribed_natural_flux'] += prescribed
            if key == 'epsilon' and side == 2:
                target = o.trace(e, 'k', x, s)**1.5/(m.cl*m.H)
                mismatch = (value['epsilon']-target)/(m.Ain/m.H)
                gamma_error_integral += float(w@(mismatch*mismatch))
                gamma_error_max = max(gamma_error_max, float(np.max(abs(mismatch), initial=0.)))
            face_rows.append(dict(field=key, geometry=e['geometry'], side=side,
                outward_diffusive_owner_flux=physical, prescribed_natural_flux=prescribed,
                conservative_advective_EX2_flux=adv, natural_callback=natural_callback,
                essential_face=(side == 0 or (key == 'epsilon' and side == 2))))
    for e, side in o.faces['u']:
        clock(o.deadline); x, r, w = face_rule(o, e, side, s, order, when)
        value = fields(case, e, s, x, r); a = o.sampled(e, x, r, s)
        nx, nz = core.NORMALS[side]
        tx = nx*(a['V'][:, 0, 0]-value['Pi'])+nz*a['V'][:, 0, 2]
        tz = nx*a['V'][:, 2, 0]+nz*(a['V'][:, 2, 2]-value['Pi'])
        mechanical['owner_complete_traction_power'] += float(w@(value['U'][0]*tx+value['U'][1]*tz))
    algebraic = {}; force = assembled['forcing']; residual = assembled['residual']
    for key in ('k', 'epsilon'):
        record = scalar[key]
        record['manufactured_volume_source'] = float(o.one[key]@force[key])+record['exact_prescribed_natural_flux']
        record['physical_owner_balance'] = (record['time_BDF2']+record['conservative_advection_EX2']+
            record['owner_outward_diffusive_flux']-record['production']+record['destruction']-record['manufactured_volume_source'])
        scale = m.Kin*volume if key == 'k' else m.Ain*math.log(m.R/m.H)
        record['fixed_reference_rate_scale'] = scale
        record['physical_owner_balance_normalized'] = record['physical_owner_balance']/scale
        ids = o.essential[key]; free = o.kfree if key == 'k' else o.efree
        reaction = float(o.one[key][ids]@residual[key][ids])
        free_defect = float(o.one[key][free]@residual[key][free])
        algebraic[key] = dict(essential_constant_test_reaction=reaction,
            free_constant_test_defect=free_defect,
            full_constant_test_residual=float(o.one[key]@residual[key]),
            reaction_plus_free_identity_error=(reaction+free_defect-float(o.one[key]@residual[key]))/scale,
            diffusion_constant_test_cancellation_normalized=diffusion_constant[key]/scale,
            physical_minus_algebraic_balance=(record['physical_owner_balance']-reaction-free_defect)/scale,
            expected_owner_minus_prescribed_flux=(record['owner_outward_diffusive_flux']-record['exact_prescribed_natural_flux'])/scale)
    Ucoeff = s['flow'][:2*o.nu]
    mechanical['manufactured_momentum_work'] = (float(Ucoeff@force['flow'][:2*o.nu])-
        mechanical['bulk_body_work']-mechanical['sampled_AWF_body_work'])
    mechanical['physical_owner_balance'] = (mechanical['velocity_time_BDF2_work']+mechanical['advection_EX2_work']+
        mechanical['stress_volume_work']+mechanical['pressure_volume_work']-mechanical['owner_complete_traction_power']-
        mechanical['manufactured_momentum_work']-mechanical['bulk_body_work']-
        sum(mechanical[key] for key in ('sampled_AWF_convection_work', 'sampled_AWF_pressure_work',
                                        'sampled_AWF_known_k_work', 'sampled_AWF_body_work')))
    mechanical['fixed_reference_rate_scale'] = volume
    mechanical['physical_owner_balance_normalized'] = mechanical['physical_owner_balance']/volume
    bdf_identity = (kinetic['current']+kinetic['current_extrapolation']-kinetic['prior']-
                    kinetic['prior_extrapolation']+kinetic['second_difference'])/(4*dt)
    kinetic['BDF2_G_energy_identity_work'] = bdf_identity
    kinetic['BDF2_G_energy_identity_error_normalized'] = (bdf_identity-mechanical['velocity_time_BDF2_work'])/volume
    free = o.vfree; essential = np.setdiff1d(np.arange(2*o.nu), free)
    mechanical_algebraic = dict(essential_velocity_reaction_power=float(Ucoeff[essential]@residual['flow'][essential]),
        free_velocity_defect_power=float(Ucoeff[free]@residual['flow'][free]),
        full_velocity_test_residual=float(Ucoeff@residual['flow'][:2*o.nu]))
    mechanical_algebraic['physical_minus_algebraic_balance'] = (
        mechanical['physical_owner_balance']-mechanical_algebraic['full_velocity_test_residual'])/volume
    mechanical_algebraic['expected_negative_owner_traction_power'] = -mechanical['owner_complete_traction_power']/volume
    return dict(k=scalar['k'], epsilon=scalar['epsilon'], mechanical=mechanical,
        kinetic_energy=kinetic, scalar_algebraic_reactions=algebraic,
        mechanical_algebraic_reactions=mechanical_algebraic, physical_boundary_faces=face_rows,
        physical_Gamma_target_error=dict(reference_scaled_RMS=math.sqrt(gamma_error_integral),
                                        sampled_MAX=gamma_error_max),
        distinction='Physical balances use owner gradients and full boundary fluxes, including epsilon H and polynomial ground k flux. Essential reactions/shared assembled cancellation are separate algebraic diagnostics, not physical-flux proofs.')


def interior_jumps(case, o, order):
    """Two physical owner fluxes on each shared face; no gradient averaging."""
    total = {key: dict(signed_integral=0., absolute_integral=0.) for key in ('k', 'epsilon', 'traction_power')}
    count = 0; s = case['state']; when = case['data']['time']
    for e in o.elements:
        xl, xr, rl, rr = e['bounds']
        for side in (1, 3):
            if (side == 1 and abs(xr-1) < 1e-12) or (side == 3 and abs(rr-o.m.R) < 1e-12): continue
            matches = [f for f in o.elements if f is not e and
                ((side == 1 and abs(f['bounds'][0]-xr) < 1e-12 and
                  abs(f['bounds'][2]-rl) < 1e-12 and abs(f['bounds'][3]-rr) < 1e-12) or
                 (side == 3 and abs(f['bounds'][2]-rr) < 1e-12 and
                  abs(f['bounds'][0]-xl) < 1e-12 and abs(f['bounds'][1]-xr) < 1e-12))]
            require(len(matches) == 1, 'study requires conforming rectangular whole-face neighbors')
            f = matches[0]; count += 1; clock(o.deadline)
            x, r, w = face_rule(o, e, side, s, order, when)
            va, vb = [fields(case, g, s, x, r) for g in (e, f)]
            a, b = [o.sampled(g, x, r, s) for g in (e, f)]
            nx, nz = core.NORMALS[side]
            for key in ('k', 'epsilon'):
                if key not in va or key not in vb: continue
                gradient = 'gradk' if key == 'k' else 'gradepsilon'; D = 'Dk' if key == 'k' else 'De'
                jump = -a[D]*(nx*va[gradient][0]+nz*va[gradient][1])+b[D]*(nx*vb[gradient][0]+nz*vb[gradient][1])
                total[key]['signed_integral'] += float(w@jump)
                total[key]['absolute_integral'] += float(w@abs(jump))
            power = np.zeros(len(x))
            for law, value, sign in ((a, va, 1.), (b, vb, -1.)):
                tx = nx*(law['V'][:, 0, 0]-value['Pi'])+nz*law['V'][:, 0, 2]
                tz = nx*law['V'][:, 2, 0]+nz*(law['V'][:, 2, 2]-value['Pi'])
                power += sign*(value['U'][0]*tx+value['U'][1]*tz)
            total['traction_power']['signed_integral'] += float(w@power)
            total['traction_power']['absolute_integral'] += float(w@abs(power))
    for key, scale in (('k', o.m.Kin*o.volume), ('epsilon', o.m.Ain*math.log(o.m.R/o.m.H)),
                       ('traction_power', o.volume)):
        total[key]['fixed_reference_rate_scale'] = scale
        total[key]['signed_normalized'] = total[key]['signed_integral']/scale
        total[key]['absolute_normalized'] = total[key]['absolute_integral']/scale
    return dict(shared_face_count=count, fields=total,
                interpretation='Physical two-owner jumps. Conforming CG assembly shares test functions; it does not impose equality of these pointwise owner fluxes.')


def load_case(row, folder, table, table_sha, terminal):
    path = Path(row['native']); path = path if path.is_absolute() else folder/path
    packet = json.loads(path.read_text()); data = packet['state']; h = data['input_history']
    require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True,
            'study input is not a completed registered normal-main endpoint')
    require(packet.get('steady_mode', False) is False, 'advancing MMS study cannot use a steady artifact')
    require(data['schema'] == 'WJH27NativeCoupledStateV2' and data['manufactured'] is True,
            'study requires current native MMS state ABI')
    require(table['status'] == 'PASS' and data['reference_identity'] == table_sha,
            'qualified reference identity differs')
    require(data['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and len(data['source_identity']) == 64,
            'wrong or unbound model/source')
    variant = data['mms_variant']
    require(variant in ('Historical', 'PolynomialScalarsV1'), 'unknown explicit MMS variant')
    require(0 < data['dt'] and data['step'] >= 2 and int(data['step']) == data['step'],
            'endpoint must have completed at least two time steps')
    require(abs(data['time']-terminal) <= 1e-12 and abs(data['step']*data['dt']-terminal) <= 1e-12,
            'completed endpoint does not have the requested physical time')
    require(h['step'] == data['step']-1 and abs(h['time']+data['dt']-terminal) <= 1e-12 and
            abs(h['dt']-data['dt']) <= 1e-14*data['dt'], 'input history is not the pre-final BDF2 history')
    require(len(data['elements']) <= 128 and max(data['sizes']) <= 4096,
            'dense independent study exceeds its prospectively bounded 128-cell/4096-coefficient scope')
    for key in ('quadratureX', 'quadratureR', 'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations'):
        require(key in data, 'new study requires explicit native configuration: '+key)
    m = coupled_mms_check.IndependentMMS(table, body_amplitude=data.get('mms_body_amplitude', 0.), variant=variant)
    s = core.state(data['coefficients'], data['sizes'])
    prior, older = core.state(h['current'], data['sizes']), core.state(h['previous'], data['sizes'])
    adv, advold = core.state(h['advection'], data['sizes']), core.state(h['previous_advection'], data['sizes'])
    require(core.state_hash(s) == data['state_sha256'] and core.state_hash(prior) == data['previous_sha256'] and
            core.state_hash(adv) == data['previous_advection_sha256'], 'state or BDF2 history identity mismatch')
    require(all(np.max(abs(a['flow'][2*data['sizes'][0]:]), initial=0.) == 0 for a in (adv, advold)),
            'stored advection has a pressure load')
    return dict(label=row['label'], role=row.get('role', 'diagnostic'), path=str(path.resolve()),
                native_sha256=sha(path), data=data, state=s, history=h, prior=prior, older=older,
                stored_advection=adv, stored_previous_advection=advold, m=m)


def analysis_case(case, deadline):
    summaries = []; assemblies = []
    for (nx, nr), grid in zip(RULES, GRIDS):
        clock(deadline); a = assemble_current(case, nx, nr, deadline); assemblies.append(a)
        summaries.append(dict(quadrature=[nx, nr], owner_sampling_grid=[grid, grid],
            errors=endpoint_errors(case, a['o'], grid), current_residual=a['norms'],
            history_errors=a['history_errors'], current_advection_error=a['current_advection_error'],
            essential_trace=a['trace'], budgets=budgets(case, a, max(nx, nr)),
            interior_owner_jumps=interior_jumps(case, a['o'], max(nx, nr))))
    coarse, fine = summaries; ac, af = assemblies; o = af['o']
    diff = core.combine(af['residual'], 1., ac['residual'], -1.)
    qerror = o.dual(diff); cont = diff['flow'][2*o.nu:]
    qerror['continuity'] = math.sqrt(max(0., float(cont@la.cho_solve(o.riesz_p, cont)))/o.volume)
    qerror['gauge'] = abs(float((af['pressure_mean']-ac['pressure_mean'])@case['state']['flow'][2*o.nu:]))/(.01*o.volume)
    resolution = {key: dict(RMS_quadrature_difference=abs(fine['errors'][key]['RMS']-coarse['errors'][key]['RMS']),
                           sampled_MAX_grid_difference=abs(fine['errors'][key]['sampled_MAX']-coarse['errors'][key]['sampled_MAX'])) for key in FIELDS}
    for key, row in resolution.items():
        row['RMS_resolved'] = row['RMS_quadrature_difference'] <= max(REFERENCE_ROUNDOFF, ERROR_QUADRATURE_FRACTION*fine['errors'][key]['RMS'])
        row['sampled_MAX_resolved_diagnostic'] = row['sampled_MAX_grid_difference'] <= max(REFERENCE_ROUNDOFF, MAX_REFINEMENT_FRACTION*fine['errors'][key]['sampled_MAX'])
    budget_uncertainty = {key: abs(fine['budgets'][key]['physical_owner_balance_normalized']-
                                 coarse['budgets'][key]['physical_owner_balance_normalized']) for key in ('k', 'epsilon', 'mechanical')}
    b = fine['budgets']; identities = {}
    for key in ('k', 'epsilon'):
        a = b['scalar_algebraic_reactions'][key]
        identities[key+'_owner_vs_reaction_identity'] = abs(a['physical_minus_algebraic_balance']-a['expected_owner_minus_prescribed_flux'])
        identities[key+'_reaction_plus_free'] = abs(a['reaction_plus_free_identity_error'])
        identities[key+'_constant_test_diffusion'] = abs(a['diffusion_constant_test_cancellation_normalized'])
    a = b['mechanical_algebraic_reactions']
    identities['mechanical_owner_vs_reaction_identity'] = abs(a['physical_minus_algebraic_balance']-a['expected_negative_owner_traction_power'])
    identities['BDF2_kinetic_identity'] = abs(b['kinetic_energy']['BDF2_G_energy_identity_error_normalized'])
    current_ok = all(fine['current_residual'][key] <= limit for key, limit in core.LIMITS.items())
    quadrature_ok = all(qerror[key] <= .2*limit for key, limit in core.LIMITS.items())
    history_ok = all(v <= 1e-10 for row in fine['history_errors'].values() for v in row.values())
    if fine['current_advection_error'] is not None:
        history_ok = history_ok and all(v <= 1e-10 for v in fine['current_advection_error'].values())
    traces_ok = all(v <= 1e-10 for row in fine['essential_trace'].values() for v in row.values())
    identity_ok = all(v <= 1e-10 for v in identities.values())
    d = case['data']; native_coverage = {key: d[key] for key in ('minimum_k', 'minimum_epsilon',
        'minimum_covariance_over_k', 'minimum_principal_over_mu')}
    native_coverage_ok = (native_coverage['minimum_k'] > 0 and native_coverage['minimum_epsilon'] > 0 and
        native_coverage['minimum_covariance_over_k'] >= 0 and native_coverage['minimum_principal_over_mu'] > 0)
    return dict(label=case['label'], role=case['role'], native=case['path'], native_sha256=case['native_sha256'],
        source_identity=d['source_identity'], reference_identity=d['reference_identity'], signature=d['signature'],
        variant=d['mms_variant'], body_amplitude=d.get('mms_body_amplitude', 0.), time=d['time'], dt=d['dt'], step=d['step'],
        mesh_identity=mesh_identity(d), cells=len(d['elements']),
        native_configuration={key: d[key] for key in ('quadratureX', 'quadratureR', 'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations')},
        native_reported_coverage=native_coverage, native_reported_coverage_positive=native_coverage_ok,
        independent_principal_certificate=False, assemblies=summaries, errors=fine['errors'],
        error_resolution=resolution, budget_quadrature_uncertainty=budget_uncertainty,
        current_residual_pass=current_ok, current_residual_quadrature_pass=quadrature_ok,
        history_pass=history_ok, trace_pass=traces_ok, algebraic_identity_pass=identity_ok,
        algebraic_identity_errors=identities, residual_quadrature_difference=qerror,
        endpoint_evidence_pass=bool(current_ok and quadrature_ok and history_ok and traces_ok and identity_ok and native_coverage_ok))


def owner(case, x, r):
    found = [e for e in case['data']['elements'] if e['bounds'][0] < x < e['bounds'][1] and e['bounds'][2] < r < e['bounds'][3]]
    require(len(found) == 1, 'common-patch midpoint lacks a unique owner')
    return found[0]


def merged_edges(values):
    result = []
    for v in sorted(values):
        if not result or abs(v-result[-1]) > 1e-12: result.append(v)
    return result


def common_difference(first, second, nx, nr, grid, deadline):
    all_elements = first['data']['elements']+second['data']['elements']
    xs = merged_edges([e['bounds'][i] for e in all_elements for i in (0, 1)])
    rs = merged_edges([e['bounds'][i] for e in all_elements for i in (2, 3)])
    energy = {key: 0. for key in FIELDS}; maxima = energy.copy(); locations = {}
    for xl, xr in zip(xs[:-1], xs[1:]):
        for rl, rr in zip(rs[:-1], rs[1:]):
            clock(deadline); ea, eb = [owner(c, (xl+xr)/2, (rl+rr)/2) for c in (first, second)]
            x, r, w = rectangle_rule((xl, xr, rl, rr), nx, nr)
            error = scaled_differences(first['m'], fields(first, ea, first['state'], x, r),
                                       fields(second, eb, second['state'], x, r), r)
            for key, values in error.items(): energy[key] += float(w@(values*values))
            gx, gr, _ = rectangle_rule((xl, xr, rl, rr), grid, grid, True)
            x, r = np.r_[x, gx], np.r_[r, gr]
            error = scaled_differences(first['m'], fields(first, ea, first['state'], x, r),
                                       fields(second, eb, second['state'], x, r), r)
            for key, values in error.items():
                index = int(np.argmax(values))
                if key not in locations or values[index] > maxima[key]:
                    maxima[key] = float(values[index]); locations[key] = dict(x=float(x[index]), r=float(r[index]),
                        first_owner=ea['geometry'], second_owner=eb['geometry'])
    return {key: dict(RMS=math.sqrt(energy[key]/area(first['m'], key)), sampled_MAX=maxima[key],
                      sampled_MAX_location=locations[key]) for key in FIELDS}


def temporal_analysis(cases, results, deadline):
    selected = sorted((c for c in cases.values() if c['role'] == 'temporal'), key=lambda c: -c['data']['dt'])
    if not selected:
        return dict(status='NOT_EVALUATED', gate_pass=None, reason='No temporal-role endpoints supplied.'), None
    require(len(selected) >= 3, 'temporal assessment requires at least three completed step sizes')
    first = selected[0]; base = first['data']
    require(base['mms_variant'] == 'PolynomialScalarsV1' and abs(base['time']-TIME) <= 1e-12,
            'prospective temporal schedule requires PolynomialScalarsV1 at t=.4')
    for i, case in enumerate(selected):
        d = case['data']
        require(abs(d['dt']-.08/2**i) <= 1e-13, 'temporal schedule must begin .08,.04,.02 and halve thereafter')
        require(mesh_identity(d) == mesh_identity(base) and all(d[key] == base[key] for key in
            ('quadratureX', 'quadratureR', 'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations')),
            'temporal endpoints changed mesh, quadrature or solver tolerances')
    rows = []; passed = True
    for ca, cb in zip(selected[:-1], selected[1:]):
        a, b = results[ca['label']], results[cb['label']]; ratio = a['dt']/b['dt']; fields_order = {}
        for key in FIELDS:
            fields_order[key] = {}
            for norm in ('RMS', 'sampled_MAX'):
                ea, eb = a['errors'][key][norm], b['errors'][key][norm]
                resolved = min(ea, eb) > REFERENCE_ROUNDOFF
                observed = math.log(ea/eb)/math.log(ratio) if resolved else None
                fields_order[key][norm] = dict(order=observed, resolved_above_reference_band=resolved)
                if key in TEMPORAL_FIELDS and norm == 'RMS':
                    good = resolved and observed >= ORDER_MINIMUM and a['error_resolution'][key]['RMS_resolved'] and b['error_resolution'][key]['RMS_resolved']
                    fields_order[key][norm]['gate_pass'] = bool(good); passed = passed and good
        common = [common_difference(ca, cb, nx, nr, grid, deadline) for (nx, nr), grid in zip(RULES, GRIDS)]
        rows.append(dict(coarse=ca['label'], fine=cb['label'], step_ratio=ratio, observed_orders=fields_order,
                         common_point_differences=common))
    return dict(status='PASS' if passed else 'NOT_MET', gate_pass=bool(passed),
        gate_fields=list(TEMPORAL_FIELDS), gate_norm='reference-scaled RMS', minimum_order=ORDER_MINIMUM,
        sampled_MAX_and_pressure_gradient_orders='diagnostic only', pairs=rows,
        startup='Actual BE startup retained in each native history; subsequent steps BDF2/EX2.'), selected[-1]


def max_diameter(case):
    return max(math.hypot(e['bounds'][1]-e['bounds'][0], e['bounds'][3]-e['bounds'][2]) for e in case['data']['elements'])


def compare(row, cases, results, temporal_finest, deadline):
    first, second = cases[row['coarse']], cases[row['fine']]
    a, b = first['data'], second['data']; kind = row['kind']
    require(kind in ('quadrature', 'mesh', 'tolerance', 'time'), 'unknown controlled comparison kind')
    if kind != 'time': require(abs(a['dt']-b['dt']) <= 1e-14*a['dt'], 'controlled comparison changed actual dt')
    else:
        require(temporal_finest is None and abs(a['dt']/b['dt']-2.) <= 1e-12,
                'short-time control requires halved dt and no formal temporal-role rows')
    tolerance_keys = ('linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations')
    if kind != 'mesh': require(mesh_identity(a) == mesh_identity(b), 'non-mesh comparison changed physical spaces')
    if kind != 'quadrature': require(all(a[k] == b[k] for k in ('quadratureX', 'quadratureR')), 'non-quadrature comparison changed native integration')
    if kind != 'tolerance': require(all(a[k] == b[k] for k in tolerance_keys), 'non-tolerance comparison changed solver settings')
    if kind == 'quadrature':
        require(all(b[k] >= a[k] for k in ('quadratureX', 'quadratureR')) and
                any(b[k] > a[k] for k in ('quadratureX', 'quadratureR')), 'quadrature comparison is not a refinement')
    if kind == 'tolerance':
        require(all(b[k] <= a[k] for k in tolerance_keys[:2]) and any(b[k] < a[k] for k in tolerance_keys[:2]),
                'tolerance comparison does not tighten a convergence tolerance')
    if kind == 'mesh':
        require(len(b['elements']) > len(a['elements']) and max_diameter(second) < max_diameter(first),
                'mesh comparison does not refine cell count and maximum diameter')
        for axis in ((0, 1), (2, 3)):
            old = merged_edges([e['bounds'][i] for e in a['elements'] for i in axis])
            new = merged_edges([e['bounds'][i] for e in b['elements'] for i in axis])
            require(all(any(abs(v-w) < 1e-12 for w in new) for v in old), 'mesh study requires nested rectangular edges')
    differences = [common_difference(first, second, nx, nr, grid, deadline) for (nx, nr), grid in zip(RULES, GRIDS)]
    fine = differences[-1]; subordinate = {}; passed = None
    if temporal_finest is not None:
        require(abs(a['dt']-temporal_finest['data']['dt']) <= 1e-14*a['dt'],
                'temporal subordinate checks must use the finest temporal dt')
        target = results[temporal_finest['label']]['errors']; passed = True
        for key in FIELDS:
            threshold = max(REFERENCE_ROUNDOFF, NONTEMPORAL_FRACTION*target[key]['RMS'])
            good = fine[key]['RMS'] <= threshold
            subordinate[key] = dict(RMS_difference=fine[key]['RMS'], finest_temporal_RMS_error=target[key]['RMS'],
                threshold=threshold, pass_=bool(good), sampled_MAX_difference_diagnostic=fine[key]['sampled_MAX'],
                sampled_MAX_fraction_diagnostic=(fine[key]['sampled_MAX']/target[key]['sampled_MAX'] if target[key]['sampled_MAX'] > REFERENCE_ROUNDOFF else None))
            # The inherited subordinate criterion concerns the temporal gate
            # fields. Other fields/gradients are fully quantified diagnostics.
            if key in TEMPORAL_FIELDS: passed = passed and good
    short_time_control = {}; short_time_ok = None
    if kind == 'time':
        short_time_ok = True
        for key in FIELDS:
            target = results[first['label']]['errors'][key]['RMS']
            threshold = max(REFERENCE_ROUNDOFF, NONTEMPORAL_FRACTION*target)
            good = fine[key]['RMS'] <= threshold
            short_time_control[key] = dict(common_point_RMS_difference=fine[key]['RMS'],
                baseline_total_RMS_error=target, threshold=threshold, pass_=bool(good))
            if key in TEMPORAL_FIELDS: short_time_ok = short_time_ok and good
    balance = {}; spatial = {}; owner_pass = None
    ra, rb = results[first['label']], results[second['label']]
    if kind == 'mesh':
        owner_pass = True; ratio = max_diameter(first)/max_diameter(second)
        for key in FIELDS:
            ea, eb = ra['errors'][key]['RMS'], rb['errors'][key]['RMS']
            spatial[key] = dict(maximum_diameter_ratio=ratio,
                RMS_observed_rate=math.log(ea/eb)/math.log(ratio) if min(ea, eb) > REFERENCE_ROUNDOFF else None,
                interpretation='Total matched-time field error; a spatial-rate claim also needs a subordinate temporal-error study.')
        for key in ('k', 'epsilon', 'mechanical'):
            va, vb = [abs(rr['assemblies'][-1]['budgets'][key]['physical_owner_balance_normalized']) for rr in (ra, rb)]
            ua, ub = [rr['budget_quadrature_uncertainty'][key] for rr in (ra, rb)]
            decreased = vb+ub < va-ua
            in_band = max(va+ua, vb+ub) <= OWNER_RESOLVED_BAND
            good = decreased or in_band; owner_pass = owner_pass and good
            balance[key] = dict(coarse_absolute_normalized=va, fine_absolute_normalized=vb,
                coarse_quadrature_uncertainty=ua, fine_quadrature_uncertainty=ub,
                resolved_decrease=bool(decreased), within_fixed_reference_band=bool(in_band), pass_=bool(good))
    else:
        for key in ('k', 'epsilon', 'mechanical'):
            va, vb = [rr['assemblies'][-1]['budgets'][key]['physical_owner_balance_normalized'] for rr in (ra, rb)]
            balance[key] = dict(coarse_normalized=va, fine_normalized=vb, difference=vb-va,
                               interpretation='Quadrature/tolerance/time stability diagnostic; not spatial convergence.')
    return dict(label=row['label'], kind=kind, coarse=first['label'], fine=second['label'],
        common_point_differences=differences, difference_RMS_quadrature_change={key: abs(differences[1][key]['RMS']-differences[0][key]['RMS']) for key in FIELDS},
        difference_sampled_MAX_grid_change={key: abs(differences[1][key]['sampled_MAX']-differences[0][key]['sampled_MAX']) for key in FIELDS},
        temporal_subordinate_check=subordinate, temporal_subordinate_pass=None if passed is None else bool(passed),
        short_time_subordinate_check=short_time_control,
        short_time_subordinate_pass=None if short_time_ok is None else bool(short_time_ok),
        owner_gradient_balance_comparison=balance, owner_gradient_convergence_pass=None if owner_pass is None else bool(owner_pass),
        matched_time_spatial_rates=spatial)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True, help='Explicit coordinator-admitted wall-clock analysis cap; never launches a solver.')
    args = parser.parse_args(); started = time.monotonic()
    report = dict(schema='WJH27IndependentMMSStudyV1', status='FAIL_MMS_STUDY', full_coupled_gate_pass=False,
        temporal_gate=dict(status='NOT_EVALUATED', gate_pass=None), study_gate_pass=False,
        script_sha256=sha(__file__), continuous_mms_sha256=sha(coupled_mms_check.__file__),
        coupled_oracle_sha256=sha(core.__file__), reference_oracle_sha256=sha(core.reference_oracle.__file__),
        independent_basis_sha256=sha(Path(__file__).resolve().parents[1]/'h26'/'fixtures.py'),
        manifest_sha256=sha(args.manifest), table_sha256=sha(args.table),
        scope='Completed advancing MMS endpoints only. No solver launch, physical development, full H27 or Task3 qualification.',
        prospective=dict(independent_rules=[list(v) for v in RULES], owner_sampling_grids=list(GRIDS),
            temporal_RMS_order_minimum=ORDER_MINIMUM, temporal_schedule=[.08, .04, .02], terminal_time=TIME,
            current_residual_limits=core.LIMITS, residual_quadrature_reserve=.2, history_trace_limit=1e-10,
            RMS_quadrature_fraction=ERROR_QUADRATURE_FRACTION, sampled_MAX_resolution_fraction_diagnostic=MAX_REFINEMENT_FRACTION,
            reference_roundoff_band=REFERENCE_ROUNDOFF, temporal_subordinate_fraction=NONTEMPORAL_FRACTION,
            owner_balance_fixed_band=OWNER_RESOLVED_BAND, algebraic_identity_limit=1e-10,
            sample_MAX_is_global_bound=False, temporal_gate_norm='reference-scaled RMS',
            scales=dict(U=1., gradU=1., k='Kin', gradk='Kin/Lref; Lref=1', epsilon='Ain/r',
                        gradepsilon='Ain/(r Lref); Lref=1', Pi='tau0', gradPi='tau0/Lref; Lref=1',
                        k_balance='Kin*domain_volume/Tref', epsilon_balance='Ain*log(R/H)/Tref',
                        mechanical_balance='Uref^2*domain_volume/Tref; Uref=Lref=Tref=1')))
    failure = None
    try:
        require(sys.platform.startswith('linux'), 'numerical study requires authoritative Linux environment')
        require(math.isfinite(args.seconds) and args.seconds > 0, 'analysis requires an explicit finite positive elapsed cap')
        require(not args.out.exists(), 'immutable output already exists')
        deadline = started+args.seconds
        manifest = json.loads(args.manifest.read_text()); table = json.loads(args.table.read_text())
        require(manifest['schema'] == 'WJH27MMSStudyManifestV1', 'wrong study manifest schema')
        terminal = float(manifest['terminal_time']); require(math.isfinite(terminal) and terminal > 0, 'invalid matched physical time')
        rows = manifest['runs']; require(rows and len({v['label'] for v in rows}) == len(rows), 'empty or duplicate labeled runs')
        cases = {v['label']: load_case(v, args.manifest.resolve().parent, table, report['table_sha256'], terminal) for v in rows}
        first = next(iter(cases.values()))['data']
        for case in cases.values():
            d = case['data']
            require(all(d[key] == first[key] for key in ('source_identity', 'reference_identity', 'mms_variant')) and
                    d.get('mms_body_amplitude', 0.) == first.get('mms_body_amplitude', 0.),
                    'study changes model source, reference, MMS variant or prescribed body')
        report.update(terminal_time=terminal, variant=first['mms_variant'], source_identity=first['source_identity'],
            body_amplitude=first.get('mms_body_amplitude', 0.),
            target_representability=('U Q3xQ2, W Q2xQ3, Pi Q2, k Q2xQ2, outer epsilon Q3xQ2 exactly in native spaces; nonlinear constitutive, dependent trace and quadrature errors remain.'
                if first['mms_variant'] == 'PolynomialScalarsV1' else 'Historical k contains a quartic r factor and outer epsilon a rational 1/r factor; exact target is not in the Q3 scalar spaces.'))
        results = {}; report['runs'] = results
        for label, case in cases.items(): results[label] = analysis_case(case, deadline)
        temporal, finest = temporal_analysis(cases, results, deadline); report['temporal_gate'] = temporal
        comparisons = manifest.get('comparisons', [])
        require(len({row['label'] for row in comparisons}) == len(comparisons), 'duplicate comparison labels')
        report['comparisons'] = [compare(row, cases, results, finest, deadline) for row in comparisons]
        endpoint_ok = all(r['endpoint_evidence_pass'] for r in results.values())
        resolution_ok = all(r['error_resolution'][key]['RMS_resolved'] for r in results.values() for key in TEMPORAL_FIELDS)
        mesh_rows = [r for r in report['comparisons'] if r['kind'] == 'mesh']
        owner_ok = bool(mesh_rows) and all(r['owner_gradient_convergence_pass'] for r in mesh_rows)
        controls_complete = {r['kind'] for r in report['comparisons']} >= {'mesh', 'quadrature', 'tolerance'}
        subordinate_ok = finest is not None and controls_complete and all(r['temporal_subordinate_pass'] for r in report['comparisons'])
        report.update(endpoint_evidence_pass=bool(endpoint_ok), error_RMS_resolution_pass=bool(resolution_ok),
            owner_gradient_convergence_status=('PASS' if owner_ok else 'NOT_MET') if mesh_rows else 'NOT_EVALUATED',
            temporal_subordinate_controls_complete=controls_complete,
            temporal_subordinate_gate_pass=None if finest is None else bool(subordinate_ok))
        if finest is None:
            report['study_gate_pass'] = None
            report['status'] = 'PASS_ENDPOINT_ANALYSIS_TEMPORAL_NOT_EVALUATED' if endpoint_ok and resolution_ok else 'NOT_MET_ENDPOINT_ANALYSIS'
        else:
            passed = endpoint_ok and resolution_ok and temporal['gate_pass'] and subordinate_ok and owner_ok
            report['study_gate_pass'] = bool(passed)
            report['status'] = 'PASS_MMS_STUDY_ONLY' if passed else 'NOT_MET_MMS_STUDY'
    except Exception as error:
        failure = error; report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream: stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if failure is not None: raise failure
    print(json.dumps(dict(status=report['status'], study_gate_pass=report['study_gate_pass'],
        temporal_status=report['temporal_gate']['status'], elapsed_seconds=report['elapsed_seconds']), sort_keys=True))
    if not report['endpoint_evidence_pass'] or not report['error_RMS_resolution_pass'] or report['study_gate_pass'] is False: raise SystemExit(2)


if __name__ == '__main__':
    main()
