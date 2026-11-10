#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Command-line negative controls for the complete duration/response export."""
import argparse
import copy
import json
from pathlib import Path
import subprocess
import sys
from rans_duration_contract import MEMBERS


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    companions = {m: dict(last_time=.5) for m in MEMBERS}
    bundle = dict(verification={n: dict(status='pass') for n in (
        'aggregate_negative_controls', 'actual_initial_binding', 'targeted_regressions',
        'historical_receipt_bytes', 'current_solver_libraries', 'current_diagnostic_libraries',
        'forward_interface_tolerance', 'physical_budget_controls', 'physical_budget_accounting')},
        neutral={}, responses={}, half_amplitude=dict(status='pass'))
    for model in ('A93', 'A95'):
        cfg = dict(model=model, policy='Natural', wall='Accepted')
        bundle['neutral'][model] = dict(configuration=cfg, status='pass',
            comparisons=[dict(time=t, status='pass') for t in (.016, .1, .5)],
            final_audit_status='pass', companions=companions)
        bundle['responses'][model] = dict(configuration=cfg, status='pass',
            comparisons=[dict(time=.016, status='pass', audit_status='pass',
                              physical_budget_status='reconstructed', companions=companions)])
    bundle['ep_control'] = dict(configuration=dict(model='A93', policy='Natural', wall='EP'),
        status='pass', comparisons=[dict(time=.016, status='pass', audit_status='pass')],
        wall_comparison_status='pass')
    variants = [('positive', bundle, None)]

    def changed(name, reason, fn):
        value = copy.deepcopy(bundle)
        fn(value)
        variants.append((name, value, reason))

    changed('missing-time', 'missing_required_time',
            lambda v: v['neutral']['A93']['comparisons'].pop(1))
    changed('foreign-model', 'configuration_mismatch',
            lambda v: v['neutral']['A95']['configuration'].update(model='E94E95'))
    changed('missing-response', 'A95-response', lambda v: v['responses'].pop('A95'))
    changed('missing-response-audit', 'response_final_audit_missing',
            lambda v: v['responses']['A93']['comparisons'][0].update(audit_status='not_completed'))
    changed('missing-physical-budget', 'physical_flux_reconstruction_missing',
            lambda v: v['responses']['A93']['comparisons'][0].update(physical_budget_status='not_completed'))
    changed('old-library-map', 'current_diagnostic_libraries',
            lambda v: v['verification'].pop('current_diagnostic_libraries'))
    changed('no-wall-comparison', 'common_grid_wall_comparison',
            lambda v: v['ep_control'].update(wall_comparison_status='not_completed'))
    results = []
    for name, value, reason in variants:
        request = a.output/(name+'-input.json')
        output = a.output/(name+'-result.json')
        request.write_text(json.dumps(value)+'\n')
        cmd = [sys.executable, '-B', str(Path(__file__).with_name('rans_nonlinear_evidence_contract.py')),
               '--input', str(request), '--output', str(output)]
        code = subprocess.run(cmd).returncode
        result = json.loads(output.read_text())
        ok = code == 0 and result['status'] == 'pass' if reason is None else (
            code != 0 and any(reason in f for f in result['failures']))
        results.append(dict(name=name, status='pass' if ok else 'fail', exit_code=code,
                            expected_reason=reason, result=result, command=cmd))
    result = dict(status='pass' if all(x['status']=='pass' for x in results) else 'fail',
                  checks=results)
    (a.output/'results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'])
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
