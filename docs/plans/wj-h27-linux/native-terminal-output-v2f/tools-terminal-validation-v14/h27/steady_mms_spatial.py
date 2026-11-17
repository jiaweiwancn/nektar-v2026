#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only frozen t=0 Historical MMS spatial study; no solver launches.

Uses the qualified independent steady forcing/current equations, never a
synthetic advancing history. Inputs and create-only output bind exact hashes.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

sys.dont_write_bytecode = True
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import scipy.linalg as la
import coupled_oracle as core
import coupled_mms_check as mms
import mms_study as study

TABLE_SHA = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'
RULES = ((20, 24), (24, 32))
GRIDS = (17, 33)
BAND = 1e-12
FIELDS = ('U', 'k', 'epsilon', 'gradU', 'gradk', 'gradepsilon')
require = core.require


def identity(path):
    path = Path(path)
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def is_hash(value):
    return isinstance(value, str) and len(value) == 64 and all(c in '0123456789abcdef' for c in value)


def load_case(row, folder, table, source):
    path = Path(row['native']); path = path if path.is_absolute() else folder/path
    packet = json.loads(path.read_text()); data = packet['state']
    require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True and
            packet.get('factory') == 'WJNativeCoupledV2' and packet.get('steady_mode') is True,
            'input must be a registered normal-main steady run')
    solver = packet['steady_solver']
    require(solver['converged'] is True and solver['status'] == 'CONVERGED' and
            solver['physical_time_advanced'] == 0, 'native steady solve did not converge without time advance')
    require(data['schema'] == 'WJH27NativeCoupledStateV2' and data['basis'] == 'Modified_A' and
            data['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and data['manufactured'] is True and
            data['mms_variant'] == 'Historical' and data['gauge'] is True, 'wrong steady MMS model/basis/variant')
    require(data['source_identity'] == source and data['reference_identity'] == TABLE_SHA,
            'native scientific source/reference identity mismatch')
    require(data['step'] == data['time'] == 0 and data['input_history']['step'] ==
            data['input_history']['time'] == 0, 'frozen study requires genuine step-zero/time-zero state')
    s = core.state(data['coefficients'], data['sizes'])
    seed = core.state(data['input_history']['current'], data['sizes'])
    cached = core.state(data['current_advection'], data['sizes'])
    previous = core.state(data['previous_advection'], data['sizes'])
    require(core.state_hash(s) == data['state_sha256'] == data['previous_sha256'], 'solved-state hash mismatch')
    require(core.state_hash(cached) == core.state_hash(previous) == data['advection_sha256'] ==
            data['previous_advection_sha256'], 'solved-state advection hash mismatch')
    require(np.all(cached['flow'][2*data['sizes'][0]:] == 0), 'advection cache contains pressure loads')
    m = mms.IndependentMMS(table, data.get('mms_body_amplitude', 0.), 'Historical')
    return dict(label=row['label'], role=row['role'], path=path, native=identity(path), packet=packet,
                data=data, m=m, state=s, seed=seed, cached=cached)


def check_mesh(case, n):
    elements = case['data']['elements']; H = case['m'].H
    xs = np.linspace(0., 1., n+1); rs = np.r_[np.linspace(0., H, n+1), np.linspace(H, .5, n+1)[1:]]
    expected = [(xs[i], xs[i+1], rs[j], rs[j+1]) for i in range(n) for j in range(2*n)]
    require(len(elements) == len(expected), 'mesh is not the prescribed 2*n*n inventory')
    remaining = list(expected)
    for e in elements:
        matches = [i for i, bounds in enumerate(remaining) if np.max(abs(np.asarray(e['bounds'])-bounds)) <= 1e-12]
        require(len(matches) == 1, 'mesh does not have the prescribed nested rectangular cover')
        remaining.pop(matches[0])
    require(not remaining and max(case['data']['sizes']) <= 4096, 'incomplete or oversized bounded MMS inventory')


def weak_budgets(case, o, forcing, adv, terms, residual):
    """Direct volume integrals versus full assembled tests, retaining reactions."""
    s, m = case['state'], case['m']
    scalar = {key: dict(conservative_advection=0., production=0., destruction=0., constant_test_diffusion=0.)
              for key in ('k', 'epsilon')}
    mechanical = dict(stress=0., pressure=0., convection=0., sampled_AWF_convection=0.,
                      sampled_AWF_pressure=0., sampled_AWF_known_k=0.)
    minima = dict(k=math.inf, E=math.inf, covariance_over_k=math.inf)
    for e in o.elements:
        o.clock(); x, r, w = o.rule(e, s, 0.)
        a = o.sampled(e, x, r, s); v = study.fields(case, e, s, x, r)
        G, U = v['gradU'], v['U']; divergence = G[0]+G[3]
        R = m.nu*a['M']-a['V']+(2*a['k']/3)[:, None, None]*np.eye(3)
        minima['k'] = min(minima['k'], float(np.min(a['k'])))
        minima['E'] = min(minima['E'], float(np.min(a['E'])))
        minima['covariance_over_k'] = min(minima['covariance_over_k'], float(np.min(np.linalg.eigvalsh(R)/a['k'][:, None])))
        mechanical['stress'] += float(w@np.einsum('qij,qij->q', a['G'], a['V']))
        mechanical['pressure'] -= float(w@(v['Pi']*divergence))
        conv = np.vstack((U[0]*G[0]+U[1]*G[1], U[0]*G[2]+U[1]*G[3]))
        mechanical['convection'] += float(w@np.sum(U*conv, axis=0))
        mask = r < m.h
        if np.any(mask):
            xx = x[mask]; hh = np.full(len(xx), m.H)
            atH = study.fields(case, o.owner[e['geometry']], s, xx, hh); work = w[mask]*U[0, mask]
            mechanical['sampled_AWF_convection'] -= float(work@(2*atH['U'][0]*atH['gradU'][0]))
            mechanical['sampled_AWF_pressure'] -= float(work@atH['gradPi'][0])
            mechanical['sampled_AWF_known_k'] += float(work@((2./3)*atH['gradk'][0]))
        for key in scalar:
            if key not in v: continue
            grad = v['gradk' if key == 'k' else 'gradepsilon']; row = scalar[key]
            row['conservative_advection'] += float(w@(U[0]*grad[0]+U[1]*grad[1]+v[key]*divergence))
            row['production'] += float(w@(a['production'] if key == 'k' else 1.44*a['production']*a['E']/a['k']))
            row['destruction'] += float(w@(a['E'] if key == 'k' else 1.82*a['E']**2/a['k']))
            _, bx, br = core.basis(e, key, x, r); one = o.one[key][e[key]]
            row['constant_test_diffusion'] += float(w@(a['Dk' if key == 'k' else 'De']*
                ((bx@one)*grad[0]+(br@one)*grad[1])))
    identities = {}
    for key, row in scalar.items():
        one = o.one[key]; free = o.kfree if key == 'k' else o.efree
        essential = o.essential[key]; scale = m.Kin*o.volume if key == 'k' else m.Ain*math.log(m.R/m.H)
        row.update(continuous_frozen_manufactured_weak_load=float(one@forcing[key]),
            full_constant_test_residual=float(one@residual[key]),
            essential_reaction=float(one[essential]@residual[key][essential]),
            free_test_defect=float(one[free]@residual[key][free]), fixed_reference_scale=scale)
        row['independent_weak_balance'] = (row['constant_test_diffusion']+row['conservative_advection']-
            row['production']+row['destruction']-row['continuous_frozen_manufactured_weak_load'])
        identities[key+'_weak_balance'] = abs(row['independent_weak_balance']-row['full_constant_test_residual'])/scale
        identities[key+'_reaction_partition'] = abs(row['essential_reaction']+row['free_test_defect']-row['full_constant_test_residual'])/scale
        identities[key+'_constant_diffusion'] = abs(row['constant_test_diffusion'])/scale
        identities[key+'_conservative_advection'] = abs(row['conservative_advection']-float(one@adv[key]))/scale
    u = s['flow'][:2*o.nu]; full = float(u@residual['flow'][:2*o.nu])
    weak = sum(mechanical.values())-float(u@forcing['flow'][:2*o.nu])
    essential = np.setdiff1d(np.arange(2*o.nu), o.vfree)
    mechanical.update(continuous_frozen_manufactured_weak_work=float(u@forcing['flow'][:2*o.nu]),
        independent_velocity_test_balance=weak, full_velocity_test_residual=full,
        essential_reaction_power=float(u[essential]@residual['flow'][essential]),
        free_velocity_defect_power=float(u[o.vfree]@residual['flow'][o.vfree]))
    identities['mechanical_weak_balance'] = abs(weak-full)/o.volume
    identities['mechanical_reaction_partition'] = abs(mechanical['essential_reaction_power']+
        mechanical['free_velocity_defect_power']-full)/o.volume
    return dict(k=scalar['k'], epsilon=scalar['epsilon'], mechanical=mechanical, identity_errors=identities,
        independent_sampled_admissibility=minima,
        interpretation='Steady weak accounts retain prescribed natural terms and essential reactions. No reaction is forced to zero; physical two-owner flux jumps are reported separately.')


def analyze(case, deadline, row):
    data, s = case['data'], case['state']; assemblies = []; previous = None
    row.update(label=case['label'], role=case['role'], native=case['native'], state_sha256=data['state_sha256'],
        source_identity=data['source_identity'], reference_identity=data['reference_identity'],
        signature=data['signature'], mesh_identity=study.mesh_identity(data), cells=len(data['elements']),
        seed_state_sha256=core.state_hash(case['seed']), steady_solver=case['packet']['steady_solver'],
        configuration={key: data[key] for key in ('quadratureX', 'quadratureR', 'linear_tolerance', 'nonlinear_tolerance')},
        assemblies=assemblies, gate_pass=False)
    for (nx, nr), grid in zip(RULES, GRIDS):
        study.clock(deadline); o = core.Oracle(data, case['m'], nx, nr, deadline)
        o.face = lambda e, side: study.face_rule(o, e, side, s, max(nx, nr), 0.)
        forcing = o.forcing(s, 0., steady=True)
        adv = o.advection(s, 0.)[0]; terms = o.current(s, 0.)
        residual, norms = o.residual(s, 0., 0., core.combine(forcing, 1., adv, -1.), terms)
        budgets = weak_budgets(case, o, forcing, adv, terms, residual)
        current = dict(quadrature=[nx, nr], errors=study.endpoint_errors(case, o, grid),
            current_residual=norms, current_trace=o.trace_error(s, 0.), seed_trace=o.trace_error(case['seed'], 0.),
            stored_current_advection_error=o.dual(core.combine(case['cached'], 1., adv, -1.)),
            removed_analytic_time_load_dual=o.dual(o.removed_analytic_time_load),
            forcing_sha256=core.state_hash(forcing), residual_sha256=core.state_hash(residual),
            weak_budgets=budgets, interior_owner_jumps=study.interior_jumps(case, o, max(nx, nr)))
        assemblies.append(current)
        if previous is not None:
            delta = core.combine(residual, 1., previous[0], -1.); qerror = o.dual(delta)
            c = delta['flow'][2*o.nu:]
            qerror['continuity'] = math.sqrt(max(0., float(c@la.cho_solve(o.riesz_p, c)))/o.volume)
            qerror['gauge'] = abs(float((o.mean-previous[1])@s['flow'][2*o.nu:]))/(.01*o.volume)
            row['residual_quadrature_difference'] = qerror
        previous = (residual, o.mean.copy())
        del o, terms
    fine = assemblies[-1]; row['errors'] = fine['errors']
    row['error_quadrature_difference'] = {key: abs(fine['errors'][key]['RMS']-assemblies[0]['errors'][key]['RMS']) for key in study.FIELDS}
    coverage = {key: float(data[key]) for key in ('minimum_k', 'minimum_epsilon', 'minimum_covariance_over_k', 'minimum_principal_over_mu')}
    row['native_reported_admissibility'] = coverage
    gates = dict(current_equations=all(a['current_residual'][k] <= v for a in assemblies for k, v in core.LIMITS.items()),
        current_equation_quadrature=all(qerror[k] <= .2*v for k, v in core.LIMITS.items()),
        traces=all(v <= 1e-10 for a in assemblies for label in ('current_trace', 'seed_trace') for v in a[label].values()),
        cache=all(v <= 1e-10 for a in assemblies for v in a['stored_current_advection_error'].values()),
        weak_budgets=all(v <= 1e-10 for a in assemblies for v in a['weak_budgets']['identity_errors'].values()),
        native_admissibility=all(math.isfinite(v) and (v >= 0 if k == 'minimum_covariance_over_k' else v > 0) for k, v in coverage.items()),
        independent_sampled_admissibility=all(v >= 0 if k == 'covariance_over_k' else v > 0
            for a in assemblies for k, v in a['weak_budgets']['independent_sampled_admissibility'].items()),
        error_quadrature=all(row['error_quadrature_difference'][k] <= max(BAND, .001*fine['errors'][k]['RMS']) for k in FIELDS))
    row.update(gates=gates, gate_pass=bool(all(gates.values())), independent_principal_certificate=False)


def compare(first, second, results, kind, deadline):
    values = [study.common_difference(first, second, nx, nr, grid, deadline) for (nx, nr), grid in zip(RULES, GRIDS)]
    result = dict(kind=kind, coarse=first['label'], fine=second['label'], common_owner_differences=values)
    result['difference_quadrature_pass'] = all(abs(values[1][k]['RMS']-values[0][k]['RMS']) <=
        max(BAND, .001*values[1][k]['RMS']) for k in FIELDS)
    if kind == 'h':
        rates = {}; ratio = study.max_diameter(first)/study.max_diameter(second)
        require(abs(ratio-2.) <= 1e-12, 'spatial meshes must halve physical diameter')
        for key in study.FIELDS:
            a, b = [results[c['label']]['errors'][key]['RMS'] for c in (first, second)]
            resolved = min(a, b) > BAND
            rate = math.log(a/b)/math.log(ratio) if resolved else None
            limit = 2.5 if key.startswith('grad') else 3.5
            expected = (2 if key.startswith('grad') else 3) if key in ('Pi', 'gradPi') else (3 if key.startswith('grad') else 4)
            rates[key] = dict(coarse_RMS=a, fine_RMS=b, expected_smooth_order=expected,
                observed_order=rate, above_roundoff_band=resolved, prospective_minimum=limit if key in FIELDS else None,
                gate_pass=(rate >= limit if resolved else b <= BAND) if key in FIELDS else None,
                status='RESOLVED_RATE' if resolved else 'ROUND_OFF_BAND_NO_RATE_CLAIM' if b <= BAND else 'UNRESOLVED_RATE')
        result['rates'] = rates
        result['gate_pass'] = bool(result['difference_quadrature_pass'] and all(rates[k]['gate_pass'] for k in FIELDS))
    else:
        checks = {}
        for key in study.FIELDS:
            target = results['h8']['errors'][key]['RMS']; threshold = max(BAND, .1*target)
            checks[key] = dict(common_owner_RMS_difference=values[-1][key]['RMS'], fine_baseline_RMS_error=target,
                threshold=threshold, gate_pass=bool(values[-1][key]['RMS'] <= threshold),
                pressure_diagnostic_only=key in ('Pi', 'gradPi'))
        result.update(subordinate_checks=checks, gate_pass=bool(result['difference_quadrature_pass'] and
                      all(checks[k]['gate_pass'] for k in FIELDS)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('manifest', 'table', 'out'): parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    args = parser.parse_args()
    require(sys.platform.startswith('linux'), 'numerical spatial MMS analysis is Linux-only')
    require(math.isfinite(args.seconds) and args.seconds > 0 and not args.out.exists(), 'finite cap and create-only output required')
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    started = time.monotonic(); deadline = started+args.seconds
    report = dict(schema='WJH27SteadyMMSSpatialReportV1', status='NOT_MET', spatial_MMS_gate_pass=False,
        full_H27_qualification=False, temporal_order_qualified=False, physical_reference_gate_pass=False,
        alpha=0., physical_time_advanced=0., ex2_exercised=False, workers=1, memory_cap_bytes=cap,
        wall_cap_seconds=args.seconds, manifest=identity(args.manifest), table=identity(args.table),
        sources={key: identity(path) for key, path in dict(checker=__file__, current_oracle=core.__file__,
            continuous_MMS=mms.__file__, field_metrics=study.__file__, independent_WJ=core.reference_oracle.__file__,
            independent_basis=Path(core.__file__).resolve().parents[1]/'h26'/'fixtures.py').items()},
        frozen_gates=dict(L2_order=3.5, H1_seminorm_order=2.5, roundoff_band=BAND,
            quadrature_and_tolerance_fraction=.1, current_residual=core.LIMITS,
            weak_identity_and_trace=1e-10, independent_residual_quadrature_fraction=.2),
        scope='Only frozen t=0 Historical continuous MMS spatial qualification; no advancing histories, temporal or physical-reference claim.',
        forcing='core.Oracle.forcing(state,0,steady=True): analytically remove time derivatives; current advection, alpha=0; no projected target forcing.',
        runs={}, comparisons=[])
    try:
        manifest = json.loads(args.manifest.read_text()); table = json.loads(args.table.read_text())
        require(manifest['schema'] == 'WJH27SteadyMMSSpatialManifestV1' and is_hash(manifest['source_identity']), 'manifest schema/source binding mismatch')
        require(report['table']['sha256'] == TABLE_SHA and table['schema'] == 'WJH27V2ReferenceTablesV1' and table['status'] == 'PASS', 'reference table identity/status mismatch')
        rows = manifest['runs']; require(len(rows) == 5 and len({r['label'] for r in rows}) == 5, 'study requires five distinct prescribed runs')
        require({r['label'] for r in rows if r['role'] == 'h'} == {'h2', 'h4', 'h8'} and
                sum(r['role'] == 'quadrature' for r in rows) == sum(r['role'] == 'tolerance' for r in rows) == 1,
                'roles must be h2/h4/h8 plus one quadrature and one tolerance control')
        cases = {r['label']: load_case(r, args.manifest.parent, table, manifest['source_identity']) for r in rows}
        require(len({str(c['path'].resolve()) for c in cases.values()}) == 5, 'native artifacts must be distinct')
        baseline = cases['h8']['data']; controls = {}
        for c in cases.values():
            d = c['data']; n = int(c['label'][1:]) if c['role'] == 'h' else 8; check_mesh(c, n)
            require((d['quadratureX'], d['quadratureR']) == ((16, 24) if c['role'] == 'quadrature' else (12, 18)), 'changed prospective native quadrature')
            target = 3e-13 if c['role'] == 'tolerance' else 3e-12
            require(abs(d['nonlinear_tolerance']/target-1) <= 1e-12, 'changed prospective nonlinear tolerance')
            require(math.isfinite(d['linear_tolerance']) and d['linear_tolerance'] > 0 and
                    (d['linear_tolerance'] <= baseline['linear_tolerance'] if c['role'] == 'tolerance' else d['linear_tolerance'] == baseline['linear_tolerance']), 'unmatched linear tolerance')
            for key in ('constitutive_backend', 'mms_body_amplitude', 'stress_split', 'steady_velocity_preconditioner', 'dt'):
                require(d[key] == baseline[key], 'spatial comparison changed '+key)
            if c['role'] != 'h': controls[c['role']] = c['label']
        expected = [dict(kind='h', coarse='h2', fine='h4'), dict(kind='h', coarse='h4', fine='h8')]+[
            dict(kind=k, coarse='h8', fine=controls[k]) for k in ('quadrature', 'tolerance')]
        if 'comparisons' in manifest:
            require([{k: r[k] for k in ('kind', 'coarse', 'fine')} for r in manifest['comparisons']] == expected,
                    'optional comparisons differ from the prescribed complete matrix')
        for c in sorted(cases.values(), key=lambda c: (len(c['data']['elements']), c['label'])):
            row = {}; report['runs'][c['label']] = row; analyze(c, deadline, row)
        for pair in expected:
            report['comparisons'].append(compare(cases[pair['coarse']], cases[pair['fine']], report['runs'], pair['kind'], deadline))
        # h2/h4 rate is retained diagnostically; only the prescribed fine pair
        # determines the asymptotic order gate. Every run still needs equations.
        passed = all(r['gate_pass'] for r in report['runs'].values()) and all(r['gate_pass'] for r in report['comparisons'][1:])
        report.update(spatial_MMS_gate_pass=bool(passed), status='PASS_STEADY_MMS_SPATIAL_ONLY' if passed else 'NOT_MET',
            pressure='Pi/gradPi errors, rates and controls are diagnostics; pressure is never shifted to improve errors.',
            owner_fluxes='Actual one-sided k/epsilon/traction-power jumps reported at both independent rules; no assembled-cancellation claim of pointwise flux continuity.')
    except Exception as error:
        report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-started
    report['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'reason', 'elapsed_seconds', 'spatial_MMS_gate_pass') if key in report}))
    if not report['spatial_MMS_gate_pass']: raise SystemExit(2)


if __name__ == '__main__':
    main()
