#!/usr/bin/env python3
"""Read-only aggregation of saved R2 preflight audits; no field reconstruction.

Run on the allocated Linux summary worker. Supports stdin execution with an
explicit hash of this source; writes only compact JSON to standard output.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pin(path):
    path = Path(path)
    return dict(path=str(path), sha256=sha(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preflight', type=Path, required=True)
    parser.add_argument('--curator-sha256', required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Numerical-result curation belongs on the allocated Linux worker')
    if len(args.curator_sha256) != 64 or any(char not in '0123456789abcdef' for char in args.curator_sha256):
        raise ValueError('A complete curator source hash is required')
    root = args.preflight.resolve()
    receipt_path = root/'receipt.json'
    receipt = json.loads(receipt_path.read_text())
    result = dict(schema=2, milestone='M8-R2', kind='finite initial representation coverage',
        created_utc=datetime.now(timezone.utc).isoformat(), curator_sha256=args.curator_sha256,
        reconstruction_or_solver_calls=0, affinity=sorted(os.sched_getaffinity(0)),
        receipt=pin(receipt_path), runner_status=receipt['status'],
        producing_source=receipt['source_commit'], analysis_source=receipt['analysis_source'],
        producing_solver=receipt['inputs']['solver'], original_control=receipt['inputs']['original'],
        producing_runtime_pin=receipt['inputs']['runtime_pin'],
        runtime_before=pin(root/'runtime-before.json') if (root/'runtime-before.json').exists() else None,
        runtime_after=pin(root/'runtime-after.json') if (root/'runtime-after.json').exists() else None,
        passed_operator_gates={name: receipt['inputs'][name] for name in ('inner_gate', 'pressure_gate')},
        native_generator=receipt['inputs']['generator'], independent_oracle=receipt['inputs']['oracle'],
        R1_comparison_report=receipt['inputs']['r1_report'],
        actual_new_initializations=receipt['actual_zero_step_initializations'],
        reused_initializations=receipt['reused_zero_step_initializations'],
        native_stage_calls=receipt['native_stage_calls'], native_evolution_calls=receipt['native_evolution_calls'],
        analysis_Cxx_equivalence_to_producer=receipt.get('analysis_Cxx_equivalence_to_producer'),
        source_runtime_inputs_unchanged=receipt.get('full_runtime_source_inputs_unchanged'),
        input_hashes_unchanged=receipt.get('input_hashes_unchanged'),
        matched_initial_state_identity=receipt.get('cross_law_actual_initial_state_identity'), cases=[])
    for entry in receipt['cases']:
        row = {key: entry.get(key) for key in ('name', 'inner_stress', 'quadrature', 'status',
            'native_execution', 'new_solver_calls', 'xml_sha256', 'report_sha256',
            'initialization_wall_seconds', 'native_separate_initial_evaluation_seconds',
            'initial_state_float64_LE_sha256', 'native_history_canonical_sha256', 'native_history_times')}
        if not entry.get('audit'):
            result['cases'].append(row)
            continue
        path = Path(entry['audit'])
        audited = json.loads(path.read_text())
        row.update(audit=pin(path), audit_status=audited['status'], audit_failure=audited.get('failure'))
        if audited['status'] != 'passed':
            result['cases'].append(row)
            continue
        value = audited['result']; ranges = value['finite_ranges']
        source, trace, domain = value['independent_source_energy'], value['independent_trace'], value['native_scalar_domain']
        errors = value['errors']; error_key = max(errors, key=errors.get)
        row.update(unknowns=value['unknowns'], field_sizes=value['field_sizes'],
            independent_checks=value['independent_checks'],
            maximum_reported_independent_error=dict(name=error_key, scaled_error=errors[error_key]),
            native_scalar_inventory_counts=domain['counts'],
            native_scalar_observations=domain['observations'],
            finite_joint_coverage=value['finite_observed_joint_preflight'],
            gradient_entry_departure=ranges['outer_gradient_entry_departure'],
            gradient_worst=ranges['outer_gradient_worst'], fixed_outer_box_margins=ranges['fixed_outer_box_margins'],
            old_r8_over3_one_sided_witness=value['old_r8_over3_one_sided_derivative_witness'],
            outer_k_over_K0=ranges['outer_k_over_K0'], outer_epsilon_over_epsilon0=ranges['outer_epsilon_over_epsilon0'],
            sampled_K_over_K0=ranges['sampled_K_over_K0'],
            stricter_sampled_K_implementation_margin=value['stricter_implementation_sampled_K_box_margin'],
            scalar_minima=dict(inner_k=ranges['inner_minimum_k'], outer_k=ranges['outer_minimum_k'],
                outer_epsilon=ranges['outer_minimum_epsilon'], sampled_K=ranges['minimum_sampled_K']),
            covariance_minima=dict(inner=ranges['covariance_inner_minimum'], outer=ranges['covariance_outer_minimum'],
                projected_trace=domain['projected_trace_covariance_minimum']),
            inner_a_over_c=ranges['inner_a_over_c'],
            inner_source_stress=dict(zip(('active_volume', 'positive_k_volume', 'minimum_a_over_c',
                'maximum_a_over_c', 'zero_k_volume', 'minimum_radial_tangent'), source['inner_stress'])),
            inner_constitutive_mean_work_production_change=source['inner_constitutive_modification'],
            source_budget=dict(zip(('pre_E', 'E', 'D', 'C', 'limiter_active_volume',
                'cap_active_volume', 'minimum_lambda', 'maximum_destruction_rate'), source['robust_source'])),
            wall_sources=source['wall_sources'],
            trace=dict(actual_GL_sample_epsilon_over_epsilon0=trace['actual_epsilon_ratio_range'],
                projected_epsilon_over_epsilon0=domain['projected_trace_over_epsilon0'],
                projected_box_margin=domain['projected_trace_fixed_box_margin'],
                relative_coefficient_residual=trace['relative_to_E0'],
                chi_zero_epsilon_range=trace['chi_zero_epsilon_ratio_range']),
            normalization_rays=ranges['normalization_rays'],
            raw_unanchored_block_norms=value['raw_unanchored_residual']['actual_initial_block_norms'],
            independent_all_component_traction=value['independent_total_traction'],
            same_state_quadrature_sensitivity=value['same_state_quadrature_sensitivity'])
        if entry.get('transition_comparison'):
            transition_path = Path(entry['transition_comparison'])
            transition = json.loads(transition_path.read_text())['comparison']
            row['initial_transition_change'] = dict(receipt=pin(transition_path),
                common_physical_anchors=transition['common_physical_anchor_count'],
                maximum_absolute_chi_change=transition['maximum_absolute_chi_change'],
                maximum_sensor_change_by_group=transition['maximum_sensor_change_by_group'])
        result['cases'].append(row)
    final_provenance = (receipt.get('full_runtime_source_inputs_unchanged') is True and
        receipt.get('input_hashes_unchanged') is True and receipt.get('source_clean_end') is True and
        receipt.get('analysis_source_clean_end') is True and
        receipt.get('source_commit_end') == receipt['source_commit'] and
        receipt.get('analysis_source_commit_end') == receipt['analysis_source']['commit'])
    result['final_provenance_passed'] = final_provenance
    complete = receipt['status'] == 'passed' and final_provenance and len(result['cases']) == 4 and all(
        row.get('audit_status') == 'passed' for row in result['cases'])
    covered = complete and all(row['finite_joint_coverage'] and
        row['stricter_sampled_K_implementation_margin'] >= 0 for row in result['cases'])
    result['outer6_coverage_decision'] = ('covered on all declared finite sets at both quadratures'
        if covered else 'coverage failed' if complete else 'pending complete independent verification')
    result['outer9_needed'] = not covered if complete else None
    result['new_initializations_expected_after_reuse'] = 3
    result['scope'] = 'Finite initial coverage, not stationarity, polynomial-everywhere positivity, a converged stage, trajectory invariance or physical validation. M8 remains open.'
    result['Jacobian_cost'] = dict(measured=False, authorized_measurement='first actual auxiliary stage within bounded F family',
        planning_estimate_minutes=[10, 20], estimate_is_runtime_conclusion=False)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
