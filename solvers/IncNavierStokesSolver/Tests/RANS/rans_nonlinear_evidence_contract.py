# SPDX-License-Identifier: MIT
"""Required experiment coverage, separate from numerical and execution results."""
from rans_duration_contract import identity, require, required_times, MEMBERS


def classify(bundle):
    failures = []
    scopes = {}

    def check(name, action):
        try:
            action()
            scopes[name] = 'pass'
        except (ValueError, KeyError, TypeError) as error:
            failures.append(name + ':' + str(error))
            scopes[name] = 'failed_unresolved_or_incomplete'

    for name in ('aggregate_negative_controls', 'actual_initial_binding',
                 'targeted_regressions', 'historical_receipt_bytes',
                 'current_solver_libraries', 'current_diagnostic_libraries',
                 'forward_interface_tolerance', 'physical_budget_controls',
                 'physical_budget_accounting'):
        check(name, lambda name=name: require(
            bundle['verification'][name]['status'] == 'pass',
            'required_verification_missing_or_failed'))

    admitted_response = []
    for model in ('A93', 'A95'):
        cfg = dict(model=model, policy='Natural', wall='Accepted')

        def duration(model=model, cfg=cfg):
            row = bundle['neutral'][model]
            identity(row['configuration'], cfg)
            require(row['status'] == 'pass', 'neutral_duration_not_qualified')
            required_times(row['comparisons'], (.016, .1, .5))
            require(all(v['status'] == 'pass' for v in row['comparisons']),
                    'required_duration_comparison_failed')
            require(row['final_audit_status'] == 'pass', 'final_audit_missing_or_failed')
            require(set(row['companions']) == set(MEMBERS), 'missing_duration_companion')
            for member in MEMBERS:
                require(row['companions'][member]['last_time'] >= .5-1e-13,
                        'scalar_horizon_incomplete:' + member)

        check(model+'-duration', duration)
        if scopes[model+'-duration'] != 'pass':
            scopes[model+'-response'] = 'not_admitted'
            continue
        admitted_response.append(model)

        def response(model=model, cfg=cfg):
            row = bundle['responses'][model]
            identity(row['configuration'], cfg)
            required_times(row['comparisons'], (.016,))
            require(row['status'] == 'pass', 'response_failed_or_unresolved')
            require(all(v['status'] == 'pass' for v in row['comparisons']),
                    'required_response_comparison_failed')
            for value in row['comparisons']:
                require(value['audit_status'] == 'pass', 'response_final_audit_missing')
                require(value['physical_budget_status'] in ('pass', 'reconstructed'),
                        'physical_flux_reconstruction_missing')
                require(set(value['companions']) == set(MEMBERS),
                        'missing_response_companion')
            if row.get('extended_to_point_one'):
                required_times(row['comparisons'], (.016, .1))

        check(model+'-response', response)

    if admitted_response:
        check('half-amplitude', lambda: require(
            bundle['half_amplitude']['status'] == 'pass',
            'half_amplitude_missing_or_failed'))

        def ep():
            row = bundle['ep_control']
            require(row['configuration']['model'] in admitted_response,
                    'ep_model_neutral_not_qualified')
            require(row['configuration']['wall'] == 'EP', 'wrong_wall_for_ep_control')
            require(row['configuration']['policy'] == 'Natural', 'wrong_policy_for_ep_control')
            require(row['status'] == 'pass', 'ep_response_missing_or_failed')
            required_times(row['comparisons'], (.016,))
            require(all(v['status'] == 'pass' and v['audit_status'] == 'pass'
                        for v in row['comparisons']), 'ep_comparison_or_audit_failed')
            require(row['wall_comparison_status'] == 'pass',
                    'common_grid_wall_comparison_missing_or_failed')
        check('selected-EP-control', ep)
    else:
        scopes['half-amplitude'] = scopes['selected-EP-control'] = 'not_admitted'
    return dict(status='pass' if not failures else 'failed_unresolved_or_incomplete',
                scopes=scopes, failures=failures)


def main():
    import argparse, json
    from pathlib import Path
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    result = classify(json.loads(a.input.read_text()))
    a.output.write_text(json.dumps(result, indent=2)+'\n')
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
