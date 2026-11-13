# SPDX-License-Identifier: MIT
"""Recompute the immutable CT6 evolution gates from finite primitives.

Collection, a successful trajectory, and a passing enclosing label do not
admit a transfer. This module uses only the Python standard library.
"""
import argparse
import json
import math
from pathlib import Path
from ct_accepted_budget_gate import physical_budgets

FIELDS = ('U', 'w', 'k', 'outer_epsilon', 'Pi')
PRIMARY = {(3, .002), (3, .001), (4, .004), (4, .002), (4, .001)}


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def degree(x, y):
    if not all(finite(v) and v >= 0 for v in (x, y)):
        return dict(status='fail', reason='missing_or_nonfinite_degree_primitive')
    margin = .85*x-y
    status = ('unresolved_floor' if max(x, y) <= 1e-7 else
              'pass' if margin > 2e-8 else 'fail' if margin < -2e-8 else
              'unresolved_uncertainty')
    return dict(status=status, coarse_error=x, fine_error=y,
                reduction_margin=margin, ratio_limit=.85, floor=1e-7,
                absolute_margin=2e-8)


def temporal(x, y):
    if not all(finite(v) and v >= 0 for v in (x, y)):
        return dict(status='fail', reason='missing_or_nonfinite_temporal_primitive')
    resolved = min(x, y) > 1e-6
    order = math.log2(x/y) if resolved else None
    return dict(status=('pass' if .8 <= order <= 1.2 else 'fail') if resolved else
                'unresolved_floor', order=order, D004_002=x, D002_001=y,
                signal_floor=1e-6, band=[.8, 1.2])


def comparisons(receipt):
    reasons = []
    if receipt.get('model') != 'A95' or not finite(receipt.get('time')) or abs(receipt['time']-.016) > 1e-13:
        reasons.append('foreign_model_or_incomplete_horizon')
    rows = receipt.get('reference_errors', {})
    expected = {str(v) for v in PRIMARY}
    if set(rows) != expected:
        reasons.append('missing_or_foreign_degree_timestep_member')
    accuracy = {}
    for member, errors in rows.items():
        if set(errors) != set(FIELDS):
            reasons.append('missing_or_foreign_field:'+member)
        for field in FIELDS:
            q = errors.get(field, {})
            v, scale = q.get('scaled_weighted_RMS'), q.get('scale')
            raw = q.get('weighted_RMS')
            valid = (finite(v) and v >= 0 and finite(scale) and scale > 0 and
                     finite(raw) and raw >= 0 and abs(raw/scale-v) <= 2e-14*max(1., abs(v)))
            limit = .01 if field == 'Pi' else .005
            okay = valid and v <= limit
            accuracy[member+':'+field] = dict(status='pass' if okay else 'fail',
                                             scaled_RMS=v, limit=limit)
            if not okay:
                reasons.append('field_accuracy_or_invalid_scale:'+member+':'+field)
    degrees = {}
    orders = {}
    for field in FIELDS:
        c = rows.get(str((3, .001)), {}).get(field, {}).get('scaled_weighted_RMS')
        f = rows.get(str((4, .001)), {}).get(field, {}).get('scaled_weighted_RMS')
        degrees[field] = degree(c, f)
        q = receipt.get('temporal', {}).get('4', {}).get('orders', {}).get(field, {})
        orders[field] = temporal(q.get('D004_002'), q.get('D002_001'))
    degree_ok = (any(q['status'] == 'pass' for q in degrees.values()) and
                 all(q['status'] in ('pass', 'unresolved_floor') for q in degrees.values()))
    order_ok = (any(q['status'] == 'pass' for q in orders.values()) and
                all(q['status'] in ('pass', 'unresolved_floor') for q in orders.values()))
    if not degree_ok:
        reasons.append('required_degree_reduction')
    if not order_ok:
        reasons.append('required_temporal_order_or_unresolved_signal')
    return dict(status='pass' if not reasons else 'fail_or_unresolved', reasons=reasons,
                field_accuracy=accuracy, degree=degrees, temporal_order=orders,
                interpretation='two fixed degrees, not a measured spatial order; below-floor orders stay unresolved')


def audit(receipt):
    result = physical_budgets(receipt)
    reasons = list(result['reasons'])
    mass = receipt.get('mass', {})
    flux, scale = mass.get('outward_net'), mass.get('scale')
    if not (finite(flux) and finite(scale) and scale > 0 and abs(flux)/scale <= 1e-8):
        reasons.append('missing_or_failed_physical_mass')
    matrix = receipt.get('captured_k_matrix', {})
    error, scale = matrix.get('maximum_defect'), matrix.get('scale')
    if not (finite(error) and error >= 0 and finite(scale) and scale > 0 and
            error/scale <= 2e-10 and matrix.get('captured_rows', 0) > 0):
        reasons.append('missing_or_failed_actual_captured_stage')
    for label, row in receipt.get('budgets', {}).items():
        v = row.get('stage_coefficient_validation', {})
        for key in ('actual_flow_stage_maximum_difference', 'actual_pre_k_stage_maximum_difference'):
            value = v.get(key)
            if not finite(value) or value < 0 or value > 2e-13:
                reasons.append('missing_or_failed_stage_binding:'+label+':'+key)
    return dict(status='pass' if not reasons else 'fail', reasons=reasons,
                complete_physical_budget=result)


def primary_coverage(members):
    reasons, keys = [], []
    for row in members:
        key = (row.get('degree'), row.get('dt'))
        keys.append(key)
        if row.get('model') != 'A95' or key not in PRIMARY:
            reasons.append('foreign_model_space_or_timestep')
        t = row.get('time')
        if not finite(t) or abs(t-.016) > 1e-13 or row.get('failure'):
            reasons.append('incomplete_or_rejected_endpoint')
        dt = row.get('dt')
        if finite(dt) and dt > 0:
            steps = round(.016/dt)
            times = row.get('accepted_state_times', [])
            if (row.get('steps') != steps or len(times) != steps+1 or
                any(not finite(v) or abs(v-i*dt) > 1e-13 for i, v in enumerate(times))):
                reasons.append('missing_or_foreign_accepted_state_coverage')
        if row.get('corrections_and_additional_solves') != 0:
            reasons.append('missing_or_nonzero_governing_corrections')
        if not row.get('actual_history_sha256') or not row.get('actual_report_sha256'):
            reasons.append('missing_actual_history_or_report_identity')
        expected_count = {3:1273, 4:2283}.get(row.get('degree'))
        if row.get('stored_slots') != expected_count:
            reasons.append('incorrect_native_count')
        metrics = row.get('independent_endpoint_metrics', {})
        for field in FIELDS:
            v = metrics.get('fields', {}).get(field, {}).get('scaled_weighted_RMS')
            if not finite(v) or v < 0 or v > (.01 if field == 'Pi' else .005):
                reasons.append('missing_or_failed_independent_endpoint_field:'+field)
        for name in ('ground_scaled_RMS', 'full_interface_jump_scaled_RMS',
                     'inner_interface_reference_error_scaled_RMS', 'outer_interface_reference_error_scaled_RMS'):
            v = metrics.get('traction', {}).get(name)
            if not finite(v) or v < 0 or v > .005:
                reasons.append('missing_or_failed_independent_endpoint_traction:'+name)
    if len(keys) != len(PRIMARY) or set(keys) != PRIMARY:
        reasons.append('missing_duplicate_or_foreign_primary_member')
    for p in (3, 4):
        ids = [v.get('actual_initial_sha256') for v in members if v.get('degree') == p]
        if not ids or not all(ids) or len(set(ids)) != 1:
            reasons.append('same_space_actual_initial_coefficients_differ_or_missing')
    return dict(status='pass' if not reasons else 'fail', reasons=sorted(set(reasons)))


def classify(value):
    coverage = primary_coverage(value.get('members', []))
    comparison = comparisons(value.get('comparison', {}))
    members = value.get('members', [])
    audits = value.get('audits', {})
    results = {key:audit(v) for key, v in audits.items()}
    expected = {v.get('actual_report_sha256') for v in members}
    actual = {v.get('budget_binding', {}).get('native_report_sha256') for v in audits.values()}
    identities = {v.get('actual_report_sha256'):v for v in members}
    consistent = True
    for row in audits.values():
        b = row.get('budget_binding', {})
        member = identities.get(b.get('native_report_sha256'), {})
        if (b.get('model') != 'A95' or b.get('time') != member.get('time') or
            b.get('dt') != member.get('dt') or not member):
            consistent = False
    bound = (len(results) == len(PRIMARY) and None not in expected and
             len(actual) == len(PRIMARY) and actual == expected and
             consistent and
             all(v['status'] == 'pass' for v in results.values()))
    # These prerequisite receipts are independently recomputed by the compact
    # exporter. Their complete body and artifact hashes remain separate inputs.
    required = {'continuous_oracle', 'curved_kernel', 'geometry_operator',
                'frozen_pressure_action', 'static_transmission', 'zero_curvature',
                'negative_controls', 'native_restart_rollback', 'resources'}
    prerequisites = value.get('prerequisites', {})
    shared = set(prerequisites) == required and all(v == 'pass' for v in prerequisites.values())
    admitted = shared and coverage['status'] == 'pass' and comparison['status'] == 'pass' and bound
    return dict(schema='CT6EvolutionAdmissionV1', status='pass' if admitted else 'not_admitted',
                primary_coverage=coverage, primary_comparison=comparison, physical_audits=results,
                complete_bound_audit_coverage='pass' if bound else 'fail_or_incomplete',
                shared_prerequisites='pass' if shared else 'fail_or_incomplete',
                A93_E94E95_transfer_eligible=admitted, curved_production_eligible=False,
                physical_experimental_applicability='not_assessed',
                scope='A95 admission only; all-model evolving qualification additionally requires the admitted transfers')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    result = classify(json.loads(a.input.read_text()))
    def serializable(value):
        # Reject before serialization. Keep the failed report writable even
        # when a deliberately invalid input contains NaN/Inf; the original
        # primitive file remains unchanged at --input.
        if isinstance(value, float) and not math.isfinite(value):return None
        if isinstance(value, dict):return {k:serializable(v) for k,v in value.items()}
        if isinstance(value, list):return [serializable(v) for v in value]
        return value
    a.output.write_text(json.dumps(serializable(result), indent=2, allow_nan=False)+'\n')
    return 0 if result['status'] == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
