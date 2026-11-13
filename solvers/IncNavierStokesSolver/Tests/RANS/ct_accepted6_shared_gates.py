# SPDX-License-Identifier: MIT
"""Recompute CT6 shared scopes from compact receipts, not enclosing labels."""
from ct_accepted6_gate import finite

MODELS = ('A93', 'A95', 'E94E95')


def bounded(x, limit):
    return finite(x) and 0 <= x <= limit


def positive(x):
    return finite(x) and x > 0


def classify(data):
    checks = {}

    def scope(name, action):
        try:
            conditions = action()
            failed = [k for k, value in conditions.items() if not value]
            checks[name] = dict(status='pass' if not failed else 'fail', failed_primitives=failed)
        except (KeyError, TypeError, ValueError, IndexError) as error:
            checks[name] = dict(status='fail', missing_or_invalid_primitive=repr(error))

    def continuous():
        q = data['continuous_oracle']
        samples, interfaces = q['samples'], q['interfaces']
        return dict(model_region_coverage=(len(samples) == 18 and
            {(r['model'], r['outer']) for r in samples} == {(m, s) for m in MODELS for s in (False, True)}),
            interface_coverage=len(interfaces) == 9 and {(r['model'], r['q']) for r in interfaces} ==
                {(m, x) for m in MODELS for x in (-.9, 0., .9)},
            sources_computed=q['source_evaluation'] is True,
            independent_derivatives=bounded(q['independent_base_gradient_defect'], 2e-11),
            divergence=all(finite(r['divergence']) and abs(r['divergence']) <= 2e-11 for r in samples),
            transmission=all(bounded(v, 2e-10) for r in interfaces for v in r['defects'].values()),
            actual_sampled_domain=all(positive(v) for v in q['sampled_domain'].values()),
            finite_forcing=all(all(finite(v) for v in r['required_body']) and
                finite(r['k_source']) and finite(r['epsilon_source']) for r in samples))

    def kernels():
        q = data['curved_kernel'];rows = q['cases'];negative = q['negative_controls']
        kinds = {'boundary-admissible', 'boundary-inadmissible', 'unapproved-curvature', 'endpoint-mismatch', 'invalid-metric'}
        okay = dict(coverage=len(rows) == 9 and {(r['model'], r['gamma']) for r in rows} ==
            {(m, g) for m in MODELS for g in (-.25, 0., .25)},
            independent_kernel=all(bounded(r['scaled_defect'], 2e-10) and
                bounded(r['independent_quadrature_uncertainty'], 2e-10) for r in rows),
            finite_positive_geometry=all(positive(r['reference']['minimum_stretch']) and
                0 < r['reference']['cutoff'] <= r['reference']['H'] for r in rows),
            column_positivity=all(positive(r['reference']['reference_k']['minimum_over_K']) for r in rows),
            weak_moments=all(abs(v['defect']) <= 2e-10*r['reference']['counterfactual']['epsilon_amplitude']
                for r in rows for v in r['reference']['weak_moments']),
            negative_coverage=len(negative) == 15 and {(r['model'], r['kind']) for r in negative} ==
                {(m, k) for m in MODELS for k in kinds})
        words = {'boundary-inadmissible':'cutoff_budget_inadmissible',
                 'unapproved-curvature':'AcceptedVerificationV1', 'endpoint-mismatch':'endpoint',
                 'invalid-metric':'metric'}
        okay['intended_negative_outcomes'] = all(
            (r['execution']['exit_code'] == 0 and r['native']['status'] == 'pass' and r['margin'] >= 0)
            if r['kind'] == 'boundary-admissible' else
            (r['execution']['exit_code'] != 0 and r['native']['status'] == 'rejected' and
             words[r['kind']].lower() in r['native']['reason'].lower()) for r in negative)
        return okay

    def operator():
        q = data['geometry_operator']
        rows = [q[k] for k in ('3', '4')]
        initial = data['model_initialization']['cases']
        return dict(degree_coverage=set(q) == {'3', '4'},
            model_coverage=set(initial) == {'A95-p3', 'A95-p4', 'A93-p4', 'E94E95-p4'},
            all_model_initial_weak=all(bounded(r['all_row']['maximum_scaled_defect'], 2e-10) and
                bounded(r['continuous_weak']['maximum_scaled_all_row_defect'], 2e-10) for r in initial.values()),
            native_counts=[r['all_row']['actual_stored_slots'] for r in rows] == [188, 319],
            all_rows=all(all(bounded(r['all_row'][k], 2e-10) for k in
                ('maximum_scaled_defect', 'constrained_scaled_defect', 'native_term_sum_scaled_defect')) for r in rows),
            continuous_native_tests=all(bounded(r['continuous_weak']['maximum_scaled_all_row_defect'], 2e-10) and
                bounded(r['continuous_weak']['shared_interface_cancellation'], 2e-10) for r in rows),
            physical_geometry=all(all(bounded(r['geometry'][k], 2e-11) for k in
                ('basis_transform', 'partition_and_constant_gradient', 'cofactor', 'all_element_all_mode_Stokes'))
                and all(bounded(v, 2e-11) for v in r['geometry']['physical_volume_defects'].values()) for r in rows))

    def pressure():
        q = data['frozen_pressure_action']
        return dict(degree_coverage=set(q) == {'3', '4'},
            all_frozen_blocks=all(set(r['blocks']) == {'0', '2', '3'} and
                all(bounded(v, 2e-10) for b in r['blocks'].values() for v in b.values()) for r in q.values()),
            consistent_free_mass=all(bounded(r['free_mass_symmetry'], 2e-10) and positive(r['minimum_free_mass_eigenvalue']) for r in q.values()),
            pressure_rank=all(r['additional_pressure_null_modes'] == 0 and
                positive(r['normalized_consistent_mass_pressure_minimum_singular_ratio']) and
                r['constant_pressure_action'] > 1e-8 for r in q.values()),
            chart_basis=all(bounded(r['native_chart_basis_maximum_defect'], 2e-10) for r in q.values()))

    def static():
        q = data['static_transmission']['checks']
        expected = {m+'-static-p'+str(p) for m in MODELS for p in (3, 4)}
        return dict(coverage=set(q) == expected,
            transmission=all(bounded(r['maximum_primitive_defect'], 2e-10) for r in q.values()),
            nontrivial_flux=all(abs(r['nonzero_conservative_flux']) > 1e-8 for r in q.values()),
            covariance=all(positive(r['minimum_covariance']) for r in q.values()),
            A95_domain=all(positive(r['minimum_A95_radicand']) for k, r in q.items() if k.startswith('A95-')))

    def zero():
        q = data['zero_curvature']['checks'];keys = {k for k in q if '-zero-' in k}
        return dict(coverage=keys == {m+'-zero-p3' for m in MODELS},
            same_operator=all(bounded(q[k]['packet'][v], 2e-10) for k in keys for v in
                ('raw_scaled_defect', 'frozen_matrix_defect', 'profile_defect', 'constitutive_defect')))

    def negatives():
        q = data['negative_controls']['checks']
        expected = {'omit_metric', 'wrong_outlet_orientation', 'omit_normal_motion_traction',
                    'omit_curvature_connection', 'ordinary_curved_production', 'curved_EP',
                    'variable_height_evolution', 'unqualified_cache'}
        expected.update(m+'-'+s for m in MODELS for s in ('omit_tilt', 'fixed_chart_derivative', 'flat_cutoff_on_curved_ray'))
        conditions = dict(coverage=set(q) == expected)
        for k in ('omit_metric', 'wrong_outlet_orientation'):
            conditions[k] = q[k]['observed']['maximum_scaled_defect'] > 2e-10
        for m in MODELS:
            for k in ('omit_tilt', 'fixed_chart_derivative'):
                conditions[m+'-'+k] = q[m+'-'+k]['primitive_signal'] > 1e-8
            conditions[m+'-flat_cutoff'] = abs(q[m+'-flat_cutoff_on_curved_ray']['scaled_wrong_cutoff_budget_defect']) > 1e-8
        conditions['normal_motion'] = q['omit_normal_motion_traction']['full_versus_parallel_traction_signal'] > 1e-8
        conditions['curvature'] = q['omit_curvature_connection']['inner_traction_signal'] > 1e-8
        for k in ('ordinary_curved_production', 'curved_EP', 'variable_height_evolution', 'unqualified_cache'):
            row = q[k];conditions[k] = row['execution']['exit_code'] != 0 and row['intended_reason'].lower() in row['observed_reason'].lower()
        return conditions

    def regression():
        q = data['native_restart_rollback']['checks']
        return dict(coverage=set(q) == {'fresh_CT5_A95_two_steps', 'unchanged_restart', 'deliberate_rejection_rollback', 'incompatible_geometry_history'},
            flat_bitwise=q['fresh_CT5_A95_two_steps']['accepted_states_bitwise_identical'] and q['fresh_CT5_A95_two_steps']['history_payload_bitwise_identical'] and q['fresh_CT5_A95_two_steps']['final_time'] == .004,
            restart_bitwise=q['unchanged_restart']['accepted_state_bitwise_identical'] and q['unchanged_restart']['history_payload_bitwise_identical'] and q['unchanged_restart']['time'] == .004,
            rollback=all(q['deliberate_rejection_rollback'][k] for k in ('last_accepted_state_bitwise_identical', 'rejected_candidate_equals_uninterrupted_next_state', 'initial_state_equals_uninterrupted_initial', 'native_history_bitwise_rollback')) and
                q['deliberate_rejection_rollback']['time'] == .002 and 'deliberate' in q['deliberate_rejection_rollback']['failure'],
            wrong_history=q['incompatible_geometry_history']['execution']['exit_code'] != 0 and
                any(k in q['incompatible_geometry_history']['log_tail'].lower() for k in ('physics', 'signature')))

    def resources():
        q = data['resources']
        return dict(nonempty_closed_receipts=bool(q),
            aggregate=all(bounded(r['final']['maximum_aggregate'], 72) for r in q.values()),
            project=all(bounded(r['final']['maximum_project'], 16) for r in q.values()),
            declared_limits=all(r['final']['aggregate_limit'] == 72 and r['final']['project_limit'] == 16 for r in q.values()))

    for name, action in [('continuous_oracle', continuous), ('curved_kernel', kernels),
        ('geometry_operator', operator), ('frozen_pressure_action', pressure),
        ('static_transmission', static), ('zero_curvature', zero), ('negative_controls', negatives),
        ('native_restart_rollback', regression), ('resources', resources)]:
        scope(name, action)
    return checks
