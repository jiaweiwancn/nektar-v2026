#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only implicit-algorithm and Newton-cache output equivalence check.

This compares coefficients, owner-resolved fields, dependent wall quantities
and complete exported histories. Complete PDE residuals are separate audits.
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
from numpy.polynomial.legendre import leggauss
from scipy.sparse.linalg import splu
import physical_oracle as po
import physical_transient_oracle as transient

LIMIT = 1e-8
FIELDS = ('U', 'Pi', 'k', 'epsilon')
STATE_NAMES = ('accepted', 'input_current', 'input_previous')
CACHE_NAMES = ('current_advection', 'previous_advection', 'input_advection', 'input_previous_advection')
REFRESH_REASONS = ('cold', 'alpha', 'drift', 'age', 'control', 'retry', 'injected')
EXACT_KEYS = ('model', 'source_identity', 'reference_identity', 'manufactured', 'basis', 'sizes',
              'geometry', 'elements', 'free_velocity_ids', 'k_free_ids', 'epsilon_free_ids',
              'gauge', 'step', 'time', 'dt', 'quadratureX', 'quadratureR', 'linear_tolerance',
              'nonlinear_tolerance', 'maximum_iterations', 'mms_variant', 'mms_body_amplitude',
              'constitutive_backend', 'stress_split', 'steady_velocity_preconditioner', 'anderson_depth')
require = po.require


def identity(path):
    path = Path(path)
    return dict(path=str(path.resolve()), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def clock(deadline):
    if time.monotonic() >= deadline:
        raise TimeoutError('implicit-algorithm comparison exhausted its admitted cap')


def native_json(text):
    # Native ostream may print IEEE negative zero as the integer-looking -0.
    # Preserve that bit for byte-state hashes; ordinary metadata/map integers
    # remain integers. Default parse_float already preserves -0.0/-0e0.
    return json.loads(text, parse_int=lambda token: -0. if token == '-0' else int(token))


def load(path, algorithm):
    require(path.name.endswith('.result.json'), 'only completed .result.json artifacts are accepted')
    packet = native_json(path.read_text())
    require(packet['schema'] == 'WJH27RegisteredCoupledRunV2' and packet['normal_main'] is True and
            packet['factory'] == 'WJNativeCoupledV2' and packet.get('steady_mode') is False,
            'comparison requires registered normal-main advancing results, not failed/steady states')
    d = packet['state']; h = d['input_history']
    require(d['schema'] == 'WJH27NativeCoupledStateV2' and d['basis'] == 'Modified_A' and
            d['model'] == 'FiniteNuIsotropicWallFullWJCollarV2' and d['implicit_solver'] == algorithm,
            'wrong state/model/basis/implicit solver identity')
    require(type(d['gauge']) is bool and type(d['manufactured']) is bool, 'invalid gauge/manufactured flag')
    require(len(d['source_identity']) == 64 and all(c in '0123456789abcdef' for c in d['source_identity']), 'missing source SHA256')
    scheme = transient.scheme(d)
    raw = dict(accepted=d['coefficients'], input_current=h['current'], input_previous=h['previous'],
               current_advection=d['current_advection'], previous_advection=d['previous_advection'],
               input_advection=h['advection'], input_previous_advection=h['previous_advection'])
    states = {key: transient.state(value, d['sizes']) for key, value in raw.items()}
    hashes = {key: transient.state_hash(value) for key, value in states.items()}
    require(hashes['accepted'] == d['state_sha256'] and hashes['input_current'] == d['previous_sha256'] and
            hashes['current_advection'] == d['advection_sha256'] and
            hashes['previous_advection'] == hashes['input_advection'] == d['previous_advection_sha256'],
            'stored current/previous state or explicit-history hash/shift mismatch')
    for key in CACHE_NAMES:
        require(np.all(states[key]['flow'][2*d['sizes'][0]:] == 0.), 'explicit load contains pressure/continuity entries')
    if scheme['startup']:
        require(hashes['input_current'] == hashes['input_previous'] and
                hashes['input_advection'] == hashes['input_previous_advection'], 'invalid duplicated BE startup history')
    stats = d['coupled_newton']
    for key in ('krylov_iterations', 'jacobian_actions', 'residual_evaluations', 'preconditioner_applications', 'rejected_trials'):
        require(type(stats[key]) is int and stats[key] >= 0, 'malformed Newton diagnostics')
        if algorithm == 'PicardV1': require(stats[key] == 0, 'Picard artifact contains Newton work')
    require(isinstance(stats['uncommitted_guard'], str) and isinstance(stats['physical_acceptance_guard'], str), 'missing guard provenance')
    return dict(data=d, states=states, hashes=hashes, native=identity(path), scheme=scheme)


def cache_policy_statistics(data, expected):
    """Per-step metadata, with first use age1 and exclusive refresh reasons."""
    s = data['coupled_newton']
    require(expected in ('RebuildEveryNewtonV1', 'OperatorDriftV1', 'CurrentDataDriftV1'), 'unknown expected cache policy')
    require(s['preconditioner_policy'] == expected, 'wrong Newton preconditioner policy')
    metric = 'ScaledCurrentQuadratureDataL2V1' if expected == 'CurrentDataDriftV1' else 'CurrentInverseOperatorFrobeniusV1'
    if expected == 'CurrentDataDriftV1' or 'drift_measure' in s:
        require(s['drift_measure'] == metric, 'cache drift measure is missing or differs')
    require(isinstance(s['relative_drift_limit'], (int, float)) and not isinstance(s['relative_drift_limit'], bool) and
            math.isfinite(s['relative_drift_limit']) and 0 <= s['relative_drift_limit'] <= .25,
            'relative drift limit must be finite and lie in [0,.25]')
    require(type(s['maximum_uses']) is int and s['maximum_uses'] > 0, 'maximum uses must be a positive integer')
    require(isinstance(s['step_seconds_cap'], (int, float)) and not isinstance(s['step_seconds_cap'], bool) and
            math.isfinite(s['step_seconds_cap']) and s['step_seconds_cap'] > 0, 'invalid Newton step wall cap')
    for key in ('refreshes', 'reuses', 'maximum_used_age'):
        require(type(s[key]) is int and s[key] >= 0, 'invalid per-step preconditioner counter: '+key)
    reasons = s['refresh_reasons']
    require(set(reasons) == set(REFRESH_REASONS) and
            all(type(v) is int and v >= 0 for v in reasons.values()), 'malformed exclusive refresh reasons')
    require(sum(reasons.values()) == s['refreshes'], 'exclusive refresh reasons do not sum to refreshes')
    used = s['refreshes']+s['reuses']
    require((1 <= s['maximum_used_age'] <= s['maximum_uses']) if used else s['maximum_used_age'] == 0,
            'preconditioner use age is inconsistent or exceeds maximum_uses')
    if s['reuses']: require(s['maximum_used_age'] >= 2, 'a reused factorization cannot have first-use age1')
    if expected == 'RebuildEveryNewtonV1':
        require(s['reuses'] == 0 and s['maximum_used_age'] <= 1,
                'rebuild control reused a preconditioner')
        require(all(reasons[key] == 0 for key in ('cold', 'alpha', 'drift', 'age')),
                'rebuild control violates declared refresh-reason precedence')
    else:
        require(reasons['control'] == 0, 'operator-drift policy reports a rebuild-control refresh')
    return s


def cache_log(path, final, expected, deadline):
    """Bind consecutive per-step evidence once; final JSON duplicates last row."""
    log = path.with_name(path.name[:-len('.result.json')]+'.steps.jsonl')
    require(log.is_file(), 'NewtonCacheV1 requires the sibling complete per-step .steps.jsonl')
    bound_keys = ('schema', 'model', 'source_identity', 'reference_identity', 'signature', 'sizes', 'gauge',
                  'basis', 'manufactured', 'implicit_solver', 'dt', 'quadratureX', 'quadratureR',
                  'linear_tolerance', 'nonlinear_tolerance', 'maximum_iterations', 'mms_variant',
                  'mms_body_amplitude', 'constitutive_backend', 'stress_split', 'steady_velocity_preconditioner',
                  'anderson_depth')
    policy_keys = ('relative_drift_limit', 'maximum_uses', 'step_seconds_cap', 'uncommitted_guard',
                   'physical_acceptance_guard')
    total = dict(refreshes=0, reuses=0, maximum_used_age=0,
                 refresh_reasons={key: 0 for key in REFRESH_REASONS})
    rows = []; previous = None
    with log.open() as stream:
        for line in stream:
            clock(deadline); require(line.strip(), 'blank native step-log record')
            d = native_json(line)
            require(all(d[key] == final[key] for key in bound_keys), 'step-log scientific/configuration binding differs')
            require(type(d['step']) is int and 0 < d['step'] <= final['step'] and
                    math.isfinite(d['time']) and abs(d['time']-d['step']*d['dt']) <= 2e-13*max(1., d['time']),
                    'step-log scientific clock is invalid')
            if previous is not None:
                require(d['step'] == previous['step']+1, 'step log is not consecutive; use InfoSteps=1')
                require(d['previous_sha256'] == previous['state_sha256'] and
                        d['previous_advection_sha256'] == previous['advection_sha256'],
                        'step-log accepted state/advection history chain differs')
            s = cache_policy_statistics(d, expected)
            require(all(s[key] == final['coupled_newton'][key] for key in policy_keys), 'step-log policy settings/guards changed')
            for key in ('refreshes', 'reuses'): total[key] += s[key]
            total['maximum_used_age'] = max(total['maximum_used_age'], s['maximum_used_age'])
            for key in REFRESH_REASONS: total['refresh_reasons'][key] += s['refresh_reasons'][key]
            rows.append(dict(step=d['step'], time=d['time'], refreshes=s['refreshes'], reuses=s['reuses'],
                             maximum_used_age=s['maximum_used_age'], refresh_reasons=s['refresh_reasons']))
            previous = d
    require(previous is not None and previous['step'] == final['step'] and previous['time'] == final['time'],
            'step log does not end at the completed result')
    for key in ('state_sha256', 'previous_sha256', 'advection_sha256', 'previous_advection_sha256', 'coupled_newton'):
        require(previous[key] == final[key], 'final result differs from duplicate final step-log record: '+key)
    require(total['refreshes'] > 0, 'logged trajectory does not exercise preconditioner refresh')
    if expected in ('OperatorDriftV1', 'CurrentDataDriftV1'):
        require(total['reuses'] > 0, 'logged trajectory does not exercise preconditioner reuse')
    return dict(log=identity(log), counts_are='PER_PHYSICAL_STEP', final_result_counted_again=False,
                first_logged_step=rows[0]['step'], last_logged_step=rows[-1]['step'],
                complete_logged_interval=True, totals=total, steps=rows)


def newton_cache_evidence(old, new, old_path, new_path, deadline, new_policy='OperatorDriftV1'):
    a, b = old['data'], new['data']
    sa = cache_policy_statistics(a, 'RebuildEveryNewtonV1')
    sb = cache_policy_statistics(b, new_policy)
    for key in ('relative_drift_limit', 'maximum_uses', 'step_seconds_cap', 'uncommitted_guard', 'physical_acceptance_guard'):
        require(sa[key] == sb[key], 'Newton-cache comparison changed common policy setting/guard: '+key)
    first = cache_log(old_path, a, 'RebuildEveryNewtonV1', deadline)
    second = cache_log(new_path, b, new_policy, deadline)
    require([(r['step'], r['time']) for r in first['steps']] == [(r['step'], r['time']) for r in second['steps']],
            'old/new step-log intervals or scientific clocks differ')
    return dict(old=first, new=second, gate_pass=True,
                scope='Recorded refresh/reuse coverage and bounded ages; fresh residual/Jv/guard equations need independent validation.')


def inventory(data, ref):
    """Validate full signed rectangle inventory without changing any bound."""
    n, p, nk, ne = data['sizes']; g = data['geometry']; elements = data['elements']
    require(all(type(v) is int and v > 0 for v in data['sizes']), 'invalid coefficient sizes')
    require(elements and len({e['geometry'] for e in elements}) == len(elements), 'empty/duplicate geometry IDs')
    require(all(math.isfinite(float(v)) for v in g.values()) and g['left'] < g['right'] and g['top'] > ref.H,
            'invalid geometry')
    for key, value in dict(d0=ref.d0, H=ref.H, h=ref.h).items():
        require(g[key] == value, 'geometry/reference binding differs: '+key)
    tolerance = 64*np.finfo(float).eps*max(1., abs(g['left']), abs(g['right']), g['top'])
    close = lambda a, b: abs(a-b) <= tolerance
    seen = {key: set() for key in ('u', 'p', 'k', 'epsilon')}; one = np.zeros(p); signs_seen = {}
    columns = {}; owners = {}
    for e in elements:
        xl, xr, rl, rr = e['bounds']
        require(all(math.isfinite(v) for v in e['bounds']) and xr > xl and rr > rl and
                xl >= g['left']-tolerance and xr <= g['right']+tolerance and
                rl >= -tolerance and rr <= g['top']+tolerance and not (rl < ref.H-tolerance and rr > ref.H+tolerance),
                'invalid rectangle or element crossing H')
        outer = rl >= ref.H-tolerance
        require(bool(e['epsilon']) == outer, 'epsilon inventory is not outer-only')
        for key, size in (('u', n), ('p', p), ('k', nk), ('epsilon', ne)):
            ids, signs = e[key], e[po.SIGNS[key]]
            if key == 'epsilon' and not outer:
                require(not ids and not signs, 'inner epsilon coefficients present'); continue
            count = 9 if key == 'p' else 16
            require(len(ids) == len(signs) == count and len(set(ids)) == count and
                    all(type(i) is int and 0 <= i < size for i in ids) and all(s in (-1., 1.) for s in signs), 'invalid signed '+key+' map')
            seen[key].update(ids)
        for j, gid in enumerate(e['p']):
            value = 1/e['ps'][j] if j % 3 < 2 and j//3 < 2 else 0.
            require(gid not in signs_seen or signs_seen[gid] == value, 'constant-pressure signed map conflict')
            signs_seen[gid] = value; one[gid] = value
        columns.setdefault((xl, xr), []).append(e)
    ordered = sorted(columns)
    require(close(ordered[0][0], g['left']) and close(ordered[-1][1], g['right']) and
            all(close(a[1], b[0]) for a, b in zip(ordered[:-1], ordered[1:])), 'x cover has gaps/overlaps')
    for column in columns.values():
        column.sort(key=lambda e: e['bounds'][2])
        require(close(column[0]['bounds'][2], 0.) and close(column[-1]['bounds'][3], g['top']) and
                all(close(a['bounds'][3], b['bounds'][2]) for a, b in zip(column[:-1], column[1:])), 'radial cover has gaps/overlaps')
        h = [e for e in column if close(e['bounds'][2], ref.H) and e['epsilon']]
        require(len(h) == 1, 'missing/ambiguous outer H owner')
        for e in column: owners[e['geometry']] = h[0]
    for key, size in (('u', n), ('p', p), ('k', nk), ('epsilon', ne)):
        require(seen[key] == set(range(size)), 'incomplete global '+key+' map')
    for ids, size in ((data['free_velocity_ids'][0], n), (data['free_velocity_ids'][1], n),
                      (data['k_free_ids'], nk), (data['epsilon_free_ids'], ne)):
        require(len(ids) == len(set(ids)) and all(type(i) is int and 0 <= i < size for i in ids), 'invalid free mask')
    return owners, one


def blocks(s, sizes):
    n, p, _, _ = sizes
    return dict(U=s['flow'][:2*n], Pi=s['flow'][2*n:2*n+p], k=s['k'], epsilon=s['epsilon'])


def coefficient_delta(a, b, sizes, ref):
    scales = dict(U=1., Pi=ref.tau, k=ref.Kin, epsilon=ref.A/ref.H)
    aa, bb = blocks(a, sizes), blocks(b, sizes)
    return {key: dict(scaled_MAX=float(np.max(abs(aa[key]-bb[key]))/scale),
                      scaled_coefficient_RMS=float(np.sqrt(np.mean((aa[key]-bb[key])**2))/scale)) for key, scale in scales.items()}


def rule(e, order, grid=False):
    xl, xr, rl, rr = e['bounds']
    if grid: a = np.linspace(-1., 1., order); w = np.zeros(order)
    else: a, w = leggauss(order)
    x, r = np.meshgrid((xl+xr)/2+(xr-xl)*a/2, (rl+rr)/2+(rr-rl)*a/2)
    return x.ravel(), r.ravel(), (w[:, None]*w[None, :]*(xr-xl)*(rr-rl)/4).ravel()


def fields(s, e, sizes, x, r):
    n = sizes[0]; bu, dx, dr = po.basis(e, 'u', x, r); ids = np.asarray(e['u'])
    out = dict(U=np.column_stack((bu@s['flow'][ids], bu@s['flow'][n+ids])),
               gradU=np.column_stack((dx@s['flow'][ids], dr@s['flow'][ids], dx@s['flow'][n+ids], dr@s['flow'][n+ids])))
    for key, name, coefficients in (('p', 'Pi', s['flow'][2*n:]), ('k', 'k', s['k']), ('epsilon', 'epsilon', s['epsilon'])):
        if not e[key]: continue
        b, dx, dr = po.basis(e, key, x, r); c = coefficients[e[key]]
        out[name] = b@c; out['grad'+name] = np.column_stack((dx@c, dr@c))
    return out


def physical_delta(a, b, data, ref, order, grid, deadline):
    delta = {key: a[key]-b[key] for key in a}; norms = po.Norms(); means = [0., 0.]
    for e in data['elements']:
        clock(deadline); x, r, w = rule(e, order); gx, gr, gw = rule(e, grid, True)
        x, r, w = np.r_[x, gx], np.r_[r, gr], np.r_[w, gw]
        value = fields(delta, e, data['sizes'], x, r)
        for key, v in value.items():
            base = key[4:] if key.startswith('grad') else key
            scale = ref.Kin if base == 'k' else ref.tau if base == 'Pi' else ref.A/r if base == 'epsilon' else 1.
            scaled = v/scale if np.ndim(v) == 1 or np.ndim(scale) == 0 else v/scale[:, None]
            # Vector magnitude gives combined physical velocity/gradient norms.
            if scaled.ndim == 2: scaled = np.linalg.norm(scaled, axis=1)
            norms.add('whole' if base != 'epsilon' else 'outer', key, scaled, w)
        pb = po.basis(e, 'p', x, r)[0]
        for j, s in enumerate((a, b)): means[j] += float(w@(pb@s['flow'][2*data['sizes'][0]+np.asarray(e['p'])]))
    volume = (data['geometry']['right']-data['geometry']['left'])*data['geometry']['top']
    return dict(norms=norms.report({key: LIMIT for key in FIELDS}),
                pressure_mean_over_tau=[v/(volume*ref.tau) for v in means], quadrature=order, owner_grid=grid,
                gradient_gate=False, pressure_alignment='NONE; native pressure level is retained')


def dependent_delta(a, b, data, owners, ref, deadline):
    result = dict(K_over_Kin=0., Gamma_over_epsilon_H=0., cutoff_over_H=0., Ew_over_local_reference=0.,
                  native_epsilon_H_difference=0., native_Gamma_pointwise_gap_old=0., native_Gamma_pointwise_gap_new=0.)
    for e in data['elements']:
        clock(deadline); x, r, _ = rule(e, 9); gx, gr, _ = rule(e, 17, True)
        x, r = np.r_[x, gx], np.r_[r, gr]; owner = owners[e['geometry']]
        atH = np.full(len(x), ref.H); basis = po.basis(owner, 'k', x, atH)[0]
        K = [basis@s['k'][owner['k']] for s in (a, b)]
        c = [ref.cutoff(v) for v in K]; A = [v**1.5/ref.cl for v in K]
        result['K_over_Kin'] = max(result['K_over_Kin'], float(np.max(abs(K[0]-K[1]))/ref.Kin))
        result['Gamma_over_epsilon_H'] = max(result['Gamma_over_epsilon_H'], float(np.max(abs(A[0]-A[1]))/ref.A))
        result['cutoff_over_H'] = max(result['cutoff_over_H'], float(np.max(abs(c[0]-c[1]))/ref.H))
        if not e['epsilon']:
            E = [v/np.maximum(r, cutoff) for v, cutoff in zip(A, c)]
            scale = ref.A/np.maximum(r, float(ref.cutoff(ref.Kin)))
            result['Ew_over_local_reference'] = max(result['Ew_over_local_reference'], float(np.max(abs(E[0]-E[1])/scale)))
        if e['geometry'] == owner['geometry']:
            eb = po.basis(owner, 'epsilon', x, atH)[0]
            E = [eb@s['epsilon'][owner['epsilon']] for s in (a, b)]
            result['native_epsilon_H_difference'] = max(result['native_epsilon_H_difference'], float(np.max(abs(E[0]-E[1]))/(ref.A/ref.H)))
            for j, name in enumerate(('old', 'new')):
                key = 'native_Gamma_pointwise_gap_'+name
                result[key] = max(result[key], float(np.max(abs(E[j]-A[j]/ref.H))/(ref.A/ref.H)))
    result['gate_pass'] = all(v <= LIMIT for key, v in result.items() if not key.startswith('native_Gamma_pointwise_gap'))
    result['pointwise_gap_note'] = 'Diagnostic only: endpoint-constrained native Gamma projection need not equal nonlinear Gamma pointwise; exact projected trace/current equations are separate oracle gates.'
    return result


def cache_duals(old, new, data, ref, deadline, *, cache_names=CACHE_NAMES):
    sizes = dict(u=data['sizes'][0], k=data['sizes'][2], epsilon=data['sizes'][3])
    masses = {key: po.SparseMass(size) for key, size in sizes.items()}
    for e in data['elements']:
        clock(deadline); x, r, w = rule(e, 16)
        for key in sizes:
            if not e[key]: continue
            b = po.basis(e, key, x, r)[0]
            scale = ref.Kin if key == 'k' else ref.A/r if key == 'epsilon' else 1.
            masses[key].add(np.asarray(e[key]), b.T@((w*scale**2)[:, None]*b))
    factors = {key: splu(value.finish()) for key, value in masses.items()}
    volume = (data['geometry']['right']-data['geometry']['left'])*data['geometry']['top']
    evolume = (data['geometry']['right']-data['geometry']['left'])*(data['geometry']['top']-ref.H)
    result = {}; n = data['sizes'][0]
    for label in cache_names:
        clock(deadline); d = {key: old[label][key]-new[label][key] for key in old[label]}
        energy = sum(float(v@factors['u'].solve(v)) for v in (d['flow'][:n], d['flow'][n:2*n]))
        result[label] = dict(U=math.sqrt(max(0., energy)/volume),
            k=math.sqrt(max(0., float(d['k']@factors['k'].solve(d['k'])))/volume),
            epsilon=math.sqrt(max(0., float(d['epsilon']@factors['epsilon'].solve(d['epsilon'])))/evolume))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('old', 'new', 'table', 'out'): parser.add_argument('--'+key, type=Path, required=True)
    parser.add_argument('--comparison', choices=('PicardNewtonV1', 'NewtonCacheV1'), default='PicardNewtonV1')
    parser.add_argument('--new-preconditioner', choices=('OperatorDriftV1', 'CurrentDataDriftV1'), default='OperatorDriftV1')
    parser.add_argument('--seconds', type=float, required=True); args = parser.parse_args()
    require(sys.platform.startswith('linux') and math.isfinite(args.seconds) and 0 < args.seconds <= 180,
            'Linux-only analysis requires explicit cap in (0,180] seconds')
    require(not args.out.exists(), 'comparison output is immutable')
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([8*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY]); resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    began = time.monotonic(); deadline = began+args.seconds
    report = dict(schema='WJH27ImplicitAlgorithmComparisonV1', status='NOT_MET', comparison_gate_pass=False,
        comparison=args.comparison,
        full_H27_qualification=False, complete_current_equations_verified_here=False,
        native_json_parser='PreserveLiteralNegativeZeroV1; metadata/map integers retain integer type',
        old=identity(args.old), new=identity(args.new), table=identity(args.table), limit=LIMIT,
        sources={key: identity(path) for key, path in dict(checker=__file__, basis_and_reference=po.__file__,
            history_ABI=transient.__file__, Jacobi_basis=Path(po.__file__).resolve().parents[1]/'h26'/'fixtures.py').items()},
        workers=1, memory_cap_bytes=cap, wall_cap_seconds=args.seconds,
        scope='Algorithm equivalence at matching accepted endpoint and complete exported history only; current PDE equations require the separate independent oracle.')
    try:
        require(args.comparison == 'NewtonCacheV1' or args.new_preconditioner == 'OperatorDriftV1',
                '--new-preconditioner applies only to NewtonCacheV1')
        require(report['table']['sha256'] == po.TABLE_SHA, 'table is not the pinned qualified v2 reference')
        ref = po.Reference(json.loads(args.table.read_text()))
        old_solver = 'PicardV1' if args.comparison == 'PicardNewtonV1' else 'CoupledJFNKV1'
        old, new = load(args.old, old_solver), load(args.new, 'CoupledJFNKV1'); a, b = old['data'], new['data']
        for key in EXACT_KEYS: require(a[key] == b[key], 'comparison changed non-algorithm identity/configuration: '+key)
        require(a['reference_identity'] == po.TABLE_SHA, 'native reference identity mismatch')
        for key in ('step', 'time', 'dt'): require(a['input_history'][key] == b['input_history'][key], 'input-history clock differs')
        for key in ('scales', 'reference_scales'):
            require((key in a) == (key in b) and a.get(key) == b.get(key), 'reported scales differ')
        if args.comparison == 'NewtonCacheV1':
            report['new_preconditioner'] = args.new_preconditioner
            report['preconditioner_reuse_evidence'] = newton_cache_evidence(old, new, args.old, args.new, deadline, args.new_preconditioner)
        owners, _ = inventory(a, ref)
        report.update(source_identity=a['source_identity'], reference_identity=a['reference_identity'],
            exact_metadata_sha256=digest({key: a[key] for key in EXACT_KEYS}),
            algorithm_signatures=dict(old=a['signature'], new=b['signature']),
            implicit_solvers=dict(old=a['implicit_solver'], new=b['implicit_solver']),
            coupled_newton=dict(old=a['coupled_newton'], new=b['coupled_newton']),
            clocks=dict(step=a['step'], time=a['time'], dt=a['dt'], scheme=old['scheme']),
            state_and_history_hashes=dict(old=old['hashes'], new=new['hashes']),
            scales=dict(U=1., Pi=ref.tau, k=ref.Kin, epsilon_coefficient=ref.A/ref.H,
                        epsilon_field='positive local Ain/r', cutoff=ref.H, cache='all-row reference-weighted mass dual'),
            gauge=a['gauge'], pressure_policy='No alignment or constant shift, including gauge cases; compare the native canonical pressure level.',
            coefficient_differences={}, fields={}, dependent={})
        for key in STATE_NAMES+CACHE_NAMES:
            report['coefficient_differences'][key] = coefficient_delta(old['states'][key], new['states'][key], a['sizes'], ref)
        lifts = {}
        for name, scale in (('essential_flow_lift', 1.), ('essential_k_lift', ref.Kin), ('essential_epsilon_lift', ref.A/ref.H)):
            x, y = np.asarray(a[name], float), np.asarray(b[name], float)
            require(x.shape == y.shape and x.ndim == 1 and np.all(np.isfinite(x)) and np.all(np.isfinite(y)), 'invalid stored essential lift')
            lifts[name] = float(np.max(abs(x-y), initial=0.)/scale)
        report['essential_lift_scaled_MAX'] = lifts
        for key in STATE_NAMES:
            report['fields'][key] = [physical_delta(old['states'][key], new['states'][key], a, ref, q, grid, deadline)
                                     for q, grid in ((8, 9), (12, 17))]
            report['dependent'][key] = dependent_delta(old['states'][key], new['states'][key], a, owners, ref, deadline)
        report['explicit_history_full_mass_dual'] = cache_duals(old['states'], new['states'], a, ref, deadline)
        gates = dict(coefficients=all(v['scaled_MAX'] <= LIMIT for row in report['coefficient_differences'].values() for v in row.values()),
            fields=all(v.get('status', 'MET') == 'MET' for rows in report['fields'].values() for row in rows for v in row['norms'].values()),
            dependent=all(v['gate_pass'] for v in report['dependent'].values()),
            lifts=all(v <= LIMIT for v in lifts.values()),
            explicit_history_duals=all(v <= LIMIT for row in report['explicit_history_full_mass_dual'].values() for v in row.values()))
        if args.comparison == 'NewtonCacheV1': gates['preconditioner_reuse_coverage'] = report['preconditioner_reuse_evidence']['gate_pass']
        report.update(gates=gates, comparison_gate_pass=bool(all(gates.values())),
                      status='PASS_IMPLICIT_ALGORITHM_COMPARISON_ONLY' if all(gates.values()) else 'NOT_MET')
    except Exception as error:
        report.update(exception=type(error).__name__, reason=str(error))
    report['elapsed_seconds'] = time.monotonic()-began
    report['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream: json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'reason', 'elapsed_seconds') if key in report}))
    if not report['comparison_gate_pass']: raise SystemExit(2)


if __name__ == '__main__':
    main()
