# SPDX-License-Identifier: MIT
"""Configuration-specific coverage and spatial gates, independent of execution."""
import math
from pathlib import Path

TIME_TOL = 1e-13
SHORT_TIMES = (.0016, .016)
MILESTONES = (.016, .1, .5)
MEMBERS = ('B5', 'B11', 'B5-dt2', 'B11-dt2')


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def prefix(model, policy, wall):
    return ('E' if model == 'E94E95' else model) + '-' + ('N' if policy == 'Natural' else 'P') + '-' + wall


def identity(request, config):
    for key in ('model', 'policy', 'wall'):
        require(request.get(key) == config[key], 'configuration_mismatch:' + key)
    if 'descriptor' in config:
        require(request.get('descriptor') == config['descriptor'], 'configuration_mismatch:descriptor')


def required_times(rows, expected, key='time'):
    observed = [float(x[key]) for x in rows]
    require(all(math.isfinite(t) for t in observed), 'nonfinite_observation_time')
    for t in expected:
        require(sum(abs(s-t) <= TIME_TOL for s in observed) == 1,
                'missing_required_time:' + str(t))
    return observed


def residual_reduction(coarse, fine, model, uncertainty=1e-10):
    """Same physical normalization, every native equation retained alongside gates."""
    checks = []
    equations = ('u', 'k', 'w') if model in ('A93', 'A95') else ('u', 'k')
    for space in ('full_native', 'fixed_common'):
        for name in equations:
            def value(row):
                if space == 'fixed_common':
                    return row['common_test']['fields'][name]['fine_restricted_common_mass_dual']
                return (row['spatial']['native_k_partition_mass_dual']['inner'] if name == 'k'
                        else row['spatial']['native_full_mass_dual'][name])
            c, f = value(coarse), value(fine)
            require(math.isfinite(c) and math.isfinite(f) and c >= 0 and f >= 0,
                    'invalid_residual_norm')
            if max(c, f) <= 1e-8:
                status = 'unresolved_floor'
            elif c-f <= 5*uncertainty:
                status = 'fail' if f > .85*c + 5*uncertainty else 'unresolved_uncertainty'
            else:
                status = 'pass' if f <= .85*c else 'fail'
            checks.append(dict(test_space=space, equation=name, coarse=c, fine=f,
                               ratio=f/c if c else None, status=status,
                               residual_floor=1e-8, uncertainty=uncertainty, maximum_ratio=.85))
    return dict(status='pass' if all(x['status'] in ('pass', 'unresolved_floor') for x in checks) else 'failed_or_unresolved',
                checks=checks, all_equations={label:dict(native=row['spatial']['native_full_mass_dual'],
                k_partition=row['spatial']['native_k_partition_mass_dual'],
                common={k:v['fine_restricted_common_mass_dual'] for k,v in row['common_test']['fields'].items()})
                for label,row in (('coarse',coarse),('fine',fine))},
                scope='spatial reduction only; reconstruction, weak balance, stage solve and field/time gates are separate')


def scalar_coverage(observer, reports, config, horizon):
    """Check the observations against sealed reports, including restart duplicates."""
    if 'configuration' in observer:
        identity(observer['configuration'], config)
    require(observer['status'] == 'pass', 'scalar_observer_failed')
    accepted = [x for x in observer['observations'] if x['stage'] == 'accepted']
    expected = []
    previous_endpoint = None
    for report in reports:
        ts, states = report['audit_times'], report['accepted_states']
        require(len(ts) == len(states) and len(ts) > 0, 'native_accepted_shape')
        if previous_endpoint is None:
            require(abs(ts[0]) <= TIME_TOL, 'scalar_initial_time_missing')
        else:
            require(abs(ts[0]-previous_endpoint) <= TIME_TOL, 'scalar_restart_interval_gap')
        require(abs(ts[-1]-report['time']) <= TIME_TOL, 'scalar_native_terminal_time_mismatch')
        attempts=report.get('pimple_attempts',[])
        require(len(attempts)==len(states)-1, 'scalar_native_accepted_stage_count')
        for first,last,attempt in zip(ts[:-1],ts[1:],attempts,strict=True):
            require(attempt['alpha']>0 and abs(last-first-1/attempt['alpha'])<=TIME_TOL,
                    'scalar_native_accepted_time_grid')
        previous_endpoint=ts[-1]
        expected.extend(zip(ts, states, strict=True))
    require(expected and abs(max(x[0] for x in expected)-horizon) <= TIME_TOL,
            'scalar_native_horizon_incomplete')
    require(accepted and abs(max(x['time'] for x in accepted)-horizon) <= TIME_TOL,
            'scalar_horizon_incomplete')
    require(len(accepted) == len(expected), 'scalar_accepted_coverage_count')
    from rans_qualification import exact_vector_hash
    unique = {}
    duplicates = 0
    for obs, (t,state) in zip(accepted, expected, strict=True):
        require(abs(obs['time']-t) <= TIME_TOL, 'scalar_accepted_time_mismatch')
        digest = exact_vector_hash(state)
        require('state_sha256' in obs, 'scalar_state_identity_missing')
        require(obs['state_sha256'] == digest, 'scalar_accepted_state_mismatch')
        key = next((s for s in unique if abs(s-t) <= TIME_TOL), None)
        if key is None:
            unique[t] = digest
        else:
            require(unique[key] == digest, 'restart_boundary_state_mismatch')
            duplicates += 1
    return dict(status='pass', first_time=min(unique), last_time=max(unique),
                native_accepted_observations=len(expected), unique_accepted_states=len(unique),
                duplicated_restart_boundary_observations=duplicates,
                native_report_intervals=[[r['audit_times'][0],r['audit_times'][-1]] for r in reports],
                actual_state_hashes_checked='state_sha256' in accepted[0],
                scope='exact accepted-state time coverage; native stage and scalar extrema checks remain separate')


def aggregate(contract, receipt):
    """Reject pass-labelled receipts with incomplete or foreign evidence."""
    errors = []
    def check(fn):
        try: fn()
        except (ValueError,KeyError,TypeError) as e: errors.append(str(e))
    check(lambda: identity(receipt['configuration'], contract))
    required = contract['audit_times']
    for member in MEMBERS:
        def member_check(member=member):
            require(member in receipt['companions'], 'missing_companion:' + member)
            c = receipt['companions'][member]
            identity(c['configuration'], contract)
            require(c['status'] == 'pass', 'companion_failed:' + member)
            require(c['scalar_coverage']['last_time'] >= max(required)-TIME_TOL,
                    'scalar_horizon_incomplete:' + member)
        check(member_check)
    for member in MEMBERS[:2]:
        def audit_check(member=member):
            a = receipt['spatial_audits'][member]
            identity(a['configuration'], contract)
            require(a['status'] == 'pass', 'spatial_reconstruction_failed:' + member)
            required_times(a['checks'], required)
        check(audit_check)
    for t in required:
        def reduction_check(t=t):
            c,f = [next(x for x in receipt['spatial_audits'][s]['checks'] if abs(x['time']-t)<=TIME_TOL) for s in MEMBERS[:2]]
            require(residual_reduction(c,f,contract['model'])['status']=='pass', 'evolved_reduction_failed:' + str(t))
        try: check(reduction_check)
        except StopIteration: errors.append('missing_required_time:'+str(t))
    return dict(status='pass' if not errors else 'fail', failures=errors,
                coverage_status='pass' if not any('missing_' in e or 'horizon_' in e or 'configuration_' in e for e in errors) else 'fail',
                evolved_reduction_status='pass' if not any('evolved_reduction_' in e for e in errors) else 'fail')
