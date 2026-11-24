#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only analytic, non-solution controls for the physical transient audit.

Exact spatial polynomials and quadratic time factors independently generate
mass and every signed advection-cache coefficient using unsplit Gauss8 rules
(only the sampled wall term is integrated on its exact 0<h support). No native
kernel or manufactured body force is used. These deliberately nonphysical
boundary values must not be reported as a qualified physical solution.
"""
import argparse
import copy
import json
import math
from pathlib import Path
import sys
import time

import physical_oracle as po
import physical_oracle_test as steady_test
import physical_transient_oracle as transient
import numpy as np
from numpy.polynomial.legendre import leggauss


def factors(t):
    return dict(u=1+.03*t+.2*t*t, k=1+.003*t+.002*t*t,
                epsilon=1+.005*t+.001*t*t, p=1+.7*t*t)


def derivatives(t):
    return dict(u=.03+.4*t, k=.003+.004*t, epsilon=.005+.002*t)


def coefficients(base, t):
    result = transient.state(base['coefficients'], base['sizes'])
    n = base['sizes'][0]; f = factors(t)
    result['flow'][:2*n] *= f['u']; result['flow'][2*n:] *= f['p']
    for key in ('k', 'epsilon'): result[key] *= f[key]
    return result


def empty(data):
    n, p, k, e = data['sizes']
    return dict(flow=np.zeros(2*n+p), k=np.zeros(k), epsilon=np.zeros(e))


def tensor_rule(e, radial=None):
    a, w = leggauss(8)
    xl, xr, rl, rr = e['bounds']
    if radial is not None: rl, rr = radial
    xx, zz = np.meshgrid((xl+xr)/2+(xr-xl)*a/2, (rl+rr)/2+(rr-rl)*a/2)
    return xx.ravel(), zz.ravel(), np.outer(w, w).ravel()*(xr-xl)*(rr-rl)/4


def analytic_advection(data, ref, field, t):
    result = empty(data); n = data['sizes'][0]; f = factors(t)
    for e in data['elements']:
        x, r, weight = tensor_rule(e)
        u, w = f['u']*field('u', x, r), f['u']*field('w', x, r)
        conv = np.column_stack((f['u']*(u*(.004+.002*r/ref.top)+w*(.1+.002*(x+10))/ref.top),
                                f['u']*(u*.001-w*.003/ref.top)))
        b = po.basis(e, 'u', x, r)[0]; ids = np.asarray(e['u'])
        for j in range(2): result['flow'][j*n+ids] += b.T@(weight*conv[:, j])
        xl, xr, rl, rr = e['bounds']
        if rl < ref.h:
            xx, rad, ww = tensor_rule(e, (rl, min(rr, ref.h)))
            force = 2*f['u']**2*field('u', xx, ref.H)*(.004+.002*ref.H/ref.top)
            result['flow'][ids] -= po.basis(e, 'u', xx, rad)[0].T@(ww*force)
        for key in ('k', 'epsilon'):
            if not e[key]: continue
            _, dx, dr = po.basis(e, key, x, r)
            q = f[key]*field(key, x, r)
            result[key][e[key]] -= dx.T@(weight*u*q)+dr.T@(weight*w*q)
            # Native scalar conservative exterior flux includes outer epsilon H.
            for side, present in enumerate((xl == -10., xr == 15.,
                    rl == (ref.H if key == 'epsilon' else 0.), rr == ref.top)):
                if not present: continue
                a, gauss_weights = leggauss(8)
                if side < 2:
                    fx, fr = np.full(8, xl if side == 0 else xr), (rl+rr)/2+(rr-rl)*a/2
                    fw = gauss_weights*(rr-rl)/2
                else:
                    fx, fr = (xl+xr)/2+(xr-xl)*a/2, np.full(8, rl if side == 2 else rr)
                    fw = gauss_weights*(xr-xl)/2
                nx, nr = po.NORMALS[side]
                normal = f['u']*(nx*field('u', fx, fr)+nr*field('w', fx, fr))
                fq = f[key]*field(key, fx, fr)
                result[key][e[key]] += po.basis(e, key, fx, fr)[0].T@(fw*normal*fq)
    return result


def analytic_mass(data, ref, field, rate):
    result = empty(data); n = data['sizes'][0]
    integral = dict(k=0., epsilon=0., velocity_square=0.)
    for e in data['elements']:
        x, r, weight = tensor_rule(e)
        b = po.basis(e, 'u', x, r)[0]; ids = np.asarray(e['u'])
        for j, key in enumerate(('u', 'w')):
            result['flow'][j*n+ids] += b.T@(weight*rate['u']*field(key, x, r))
        integral['velocity_square'] += float(weight@(field('u', x, r)**2+field('w', x, r)**2))
        for key in ('k', 'epsilon'):
            if not e[key]: continue
            result[key][e[key]] += po.basis(e, key, x, r)[0].T@(weight*rate[key]*field(key, x, r))
            integral[key] += float(weight@field(key, x, r))
    return result, integral


def packet(ref, step, dt=.04):
    base, field = steady_test.fixture(ref)
    # Each coefficient factor multiplies a separate fresh copy of base arrays.
    now, oldtime = step*dt, (step-1)*dt
    pretime = max(0, step-2)*dt
    new, old, previous = (coefficients(base, t) for t in (now, oldtime, pretime))
    newadv, oldadv, prevadv = (analytic_advection(base, ref, field, t) for t in (now, oldtime, pretime))
    data = dict(base, basis='Modified_A', step=step, time=now, dt=dt,
        coefficients=new, state_sha256=transient.state_hash(new),
        previous_sha256=transient.state_hash(old), advection_sha256=transient.state_hash(newadv),
        previous_advection_sha256=transient.state_hash(oldadv),
        current_advection=newadv, previous_advection=oldadv,
        input_history=dict(step=step-1, time=oldtime, dt=dt, current=old,
                           previous=previous, advection=oldadv, previous_advection=prevadv))
    return data, field


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    args = parser.parse_args()
    po.require(sys.platform.startswith('linux'), 'numerical controls are Linux-only')
    po.require(math.isfinite(args.seconds) and args.seconds > 0 and not args.out.exists(),
               'finite cap and immutable output required')
    began = time.monotonic()
    report = dict(schema='WJH27PhysicalTransientOracleAnalyticControlsV1', status='FAIL',
        source=po.identity(__file__), transient_source=po.identity(transient.__file__),
        steady_source=po.identity(po.__file__), fixture_source=po.identity(steady_test.__file__),
        scope='Non-solution analytic controls; no manufactured force and no physical qualification')
    try:
        po.require(po.identity(args.table)['sha256'] == po.TABLE_SHA, 'wrong immutable table')
        ref = po.Reference(json.loads(args.table.read_text()))
        cases = {'Explicit': [], 'Implicit': []}
        for mode, step in (('Implicit', 1), ('Implicit', 2), ('Explicit', 1), ('Explicit', 2)):
            data, field = packet(ref, step)
            if mode == 'Implicit':
                data['picard_advection'] = mode
            oracle = transient.TransientOracle(data, ref, 20, began+args.seconds)
            core = oracle.assemble()
            errors = {}
            def close(name, actual, expected, scale, tolerance=2e-11):
                error = float(np.max(np.abs(np.asarray(actual)-expected), initial=0.))/scale
                errors[name] = error
                po.require(error <= tolerance, f'{name}: {error} > {tolerance}')
            po.require(core['history_cache_checks']['status'] == 'MET', 'independent analytic history-cache moments differ')
            po.require(core['accounting_identity_status'] == 'MET', 'transient accounting identities failed')
            po.require(core['transient_residual_status'] == 'NOT_MET', 'non-solution fixture falsely qualifies')
            t, dt = data['time'], data['dt']
            rate = (derivatives(t) if step > 1 else
                    {key: (factors(t)[key]-factors(0)[key])/dt for key in ('u', 'k', 'epsilon')})
            expected_mass, integral = analytic_mass(data, ref, field, rate)
            for key in transient.KEYS:
                close('every_'+key+'_BDF_mass_coefficient', oracle.mass[key], expected_mass[key], np.max(np.abs(expected_mass[key])))
            f, g = oracle.scheme['explicit_coefficients']
            expected_transport = (data['current_advection'] if mode == 'Implicit' else
                                  transient.combine(data['input_history']['advection'], f,
                                                    data['input_history']['previous_advection'], g))
            for key in transient.KEYS:
                close('every_'+key+'_advection_coefficient', oracle.advective_action[key], expected_transport[key],
                      np.max(np.abs(expected_transport[key])))
            for key in ('k', 'epsilon'):
                close(key+'_storage_inventory', core['storage_by_region']['whole'][key],
                      rate[key]*integral[key], abs(rate[key]*integral[key]))
            close('current_velocity_times_BDF_storage', core['budgets']['mechanical']['storage_velocity_test_work'],
                  factors(t)['u']*rate['u']*integral['velocity_square'],
                  abs(factors(t)['u']*rate['u']*integral['velocity_square']))
            po.require(np.all(oracle.mass['flow'][2*oracle.current.n:] == 0.), 'pressure received a time mass')
            label = mode.lower()
            alternate = ('explicit_instead_of_current_advection' if mode == 'Implicit' else
                         'current_instead_of_explicit_advection')
            expected_mutations = ['missing_time_mass', alternate, 'missing_'+label+'_AWF_convection',
                'missing_'+label+'_domain_convection', 'nonconservative_'+label+'_scalar_advection']
            if step > 1:
                expected_mutations.append('BE_instead_of_BDF2_mass')
                if mode == 'Explicit':
                    expected_mutations.append('EX1_instead_of_EX2')
            for name in expected_mutations:
                po.require(core['mutations'][name]['status'] == 'DETECTABLE_AT_DECLARED_TOLERANCE', 'inactive nonzero mutation '+name)
            cases[mode].append(dict(scheme=oracle.scheme, analytic_relative_errors=errors,
                accounting_identity_max_error=core['accounting_identity_max_error'],
                cache_checks=core['history_cache_checks'], non_solution_residuals=core['residuals'],
                mutations={name: core['mutations'][name] for name in expected_mutations}))
        report['analytic_cases'] = cases['Explicit']
        report['implicit_advection_cases'] = cases['Implicit']
        # Corrupt a free stored cache with all corresponding byte hashes updated.
        # The independently reconstructed equation must be unchanged while the
        # cache evidence fails, so there is no accidental dependence on caches.
        bad = copy.deepcopy(data)
        gid = bad['free_velocity_ids'][0][0]
        bad['input_history']['advection']['flow'][gid] += .01
        bad['previous_advection'] = copy.deepcopy(bad['input_history']['advection'])
        bad['previous_advection_sha256'] = transient.state_hash(bad['previous_advection'])
        altered = transient.TransientOracle(bad, ref, 20, began+args.seconds)
        altered_core = altered.assemble()
        po.require(altered_core['history_cache_checks']['status'] == 'NOT_MET', 'corrupted accepted cache was not detected')
        residual_change = oracle.current.dual(transient.combine(oracle.total, 1., altered.total, -1.))
        po.require(max(residual_change.values()) == 0., 'stored cache contaminated the independent equations')
        report['corrupted_cache_control'] = dict(status='PASS', equation_change_dual=residual_change,
                                                cache_checks=altered_core['history_cache_checks'])
        rejections = []
        mutations = {
            'skipped_accepted_step': lambda d: d.update(step=d['step']+1),
            'changed_timestep': lambda d: d['input_history'].update(dt=2*d['dt']),
            'wrong_scientific_time': lambda d: d.update(time=d['time']+.01),
            'nonfinite_timestep': lambda d: d.update(dt=float('nan')),
            'wrong_previous_hash': lambda d: d.update(previous_sha256='0'*64),
            'wrong_new_cache_hash': lambda d: d.update(advection_sha256='0'*64),
            'wrong_basis': lambda d: d.update(basis='Ortho_A'),
            'unknown_advection_mode': lambda d: d.update(picard_advection='Unknown'),
            'truncated_history': lambda d: d['input_history']['previous'].update(k=d['input_history']['previous']['k'][:-1]),
        }
        for name, mutate in mutations.items():
            corrupted = copy.deepcopy(data); mutate(corrupted)
            try: transient.TransientOracle(corrupted, ref, 20, began+args.seconds)
            except ValueError as error: rejections.append(dict(name=name, reason=str(error)))
            else: raise ValueError('malformed history accepted: '+name)
        pressure = copy.deepcopy(data)
        pressure['current_advection']['flow'][2*data['sizes'][0]] = .01
        pressure['advection_sha256'] = transient.state_hash(pressure['current_advection'])
        try: transient.TransientOracle(pressure, ref, 20, began+args.seconds)
        except ValueError as error: rejections.append(dict(name='pressure_advection_cache', reason=str(error)))
        else: raise ValueError('pressure advection cache accepted')
        startup, _ = packet(ref, 1)
        startup['input_history']['previous']['k'][0] += 1e-6
        try: transient.TransientOracle(startup, ref, 20, began+args.seconds)
        except ValueError as error: rejections.append(dict(name='startup_history_not_duplicated', reason=str(error)))
        else: raise ValueError('nondoubled startup history accepted')
        report.update(status='PASS', malformed_history_rejections=rejections)
    except Exception as error:
        report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-began
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'reason', 'elapsed_seconds') if key in report}))
    if report['status'] != 'PASS': raise SystemExit(2)


if __name__ == '__main__':
    main()
