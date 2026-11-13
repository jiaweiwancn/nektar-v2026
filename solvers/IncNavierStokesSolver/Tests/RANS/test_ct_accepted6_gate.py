#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Admission negatives use synthetic receipts, never scientific evidence."""
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from ct_accepted6_gate import PRIMARY, FIELDS, classify


def fixture():
    members, audits, errors = [], {}, {}
    for p, dt in sorted(PRIMARY):
        report = str((p, dt))
        member = dict(model='A95', degree=p, dt=dt, time=.016, failure='',
                      corrections_and_additional_solves=0, actual_report_sha256=report,
                      actual_history_sha256='synthetic-history-'+report,
                      actual_initial_sha256='synthetic-initial-p'+str(p),
                      stored_slots={3:1273, 4:2283}[p], steps=round(.016/dt),
                      accepted_state_times=[i*dt for i in range(round(.016/dt)+1)])
        members.append(member)
        member['independent_endpoint_metrics'] = dict(
            fields={k:dict(scaled_weighted_RMS=5e-5) for k in FIELDS},
            traction={k:5e-5 for k in ('ground_scaled_RMS', 'full_interface_jump_scaled_RMS',
                'inner_interface_reference_error_scaled_RMS', 'outer_interface_reference_error_scaled_RMS')})
        binding = dict(model='A95', time=.016, dt=dt, native_report_sha256=report)
        body = dict(budget_binding=binding, mass=dict(outward_net=1e-14, scale=.2),
                    captured_k_matrix=dict(maximum_defect=1e-14, scale=1., captured_rows=576),
                    budgets={})
        for label in ('native', 'q48f16', 'q64f16', 'q64f32'):
            body['budgets'][label] = dict(binding=binding,
                equations={k:dict(physical_scale=1., storage_backward_difference=.1,
                    stage_source=.06, stage_physical_boundary_flux=.04,
                    frozen_stage_physical_defect=0.) for k in
                    ('inner_k', 'outer_k', 'outer_epsilon', 'whole_k')},
                reconstruction={k:dict(normalization=1., maximum_defect=1e-15,
                    scaled_defect=1e-15) for k in ('k', 'outer_epsilon')},
                stage_coefficient_validation=dict(actual_flow_stage_maximum_difference=0.,
                                                  actual_pre_k_stage_maximum_difference=0.))
        audits[report] = body
        error = 1e-4 if p == 3 else 5e-5
        errors[report] = {k:dict(scale=1., weighted_RMS=error, scaled_weighted_RMS=error)
                          for k in FIELDS}
    comparison = dict(model='A95', time=.016, reference_errors=errors,
                     temporal={'4':dict(orders={k:dict(D004_002=4e-5, D002_001=2e-5)
                                                for k in FIELDS})})
    return dict(members=members, audits=audits, comparison=comparison,
        prerequisites={k:'pass' for k in ('continuous_oracle', 'curved_kernel', 'geometry_operator',
            'frozen_pressure_action', 'static_transmission', 'zero_curvature', 'negative_controls',
            'native_restart_rollback', 'resources')})


class GateTests(unittest.TestCase):
    def test_complete_synthetic_contract(self):
        self.assertEqual(classify(fixture())['status'], 'pass')

    def denied(self, value):
        result = classify(value)
        self.assertEqual(result['status'], 'not_admitted')
        self.assertFalse(result['A93_E94E95_transfer_eligible'])
        return result

    def test_false_passing_degree_label(self):
        q = fixture();row = q['comparison']['reference_errors'][str((4, .001))]['w']
        row['weighted_RMS'] = row['scaled_weighted_RMS'] = 2e-4
        q['comparison']['status'] = 'pass'
        self.assertEqual(self.denied(q)['primary_comparison']['degree']['w']['status'], 'fail')

    def test_false_passing_order_label(self):
        q = fixture();q['comparison']['temporal']['4']['orders']['Pi'].update(
            status='pass', order=1., D004_002=1e-5, D002_001=2e-5)
        self.assertEqual(self.denied(q)['primary_comparison']['temporal_order']['Pi']['status'], 'fail')

    def test_unresolved_order_is_not_invented(self):
        q = fixture();q['comparison']['temporal']['4']['orders']['k'].update(D004_002=1e-8, D002_001=1e-9)
        result = classify(q)
        self.assertEqual(result['status'], 'pass')
        self.assertIsNone(result['primary_comparison']['temporal_order']['k']['order'])

    def test_missing_endpoint_member(self):
        q = fixture();q['members'].pop()
        self.denied(q)

    def test_incomplete_endpoint(self):
        q = fixture();q['members'][0]['time'] = .008
        self.denied(q)

    def test_foreign_model(self):
        q = fixture();q['members'][0]['model'] = 'A93'
        self.denied(q)

    def test_same_space_initialization(self):
        q = fixture();q['members'][0]['actual_initial_sha256'] = 'different'
        self.denied(q)

    def test_missing_accepted_state(self):
        q = fixture();q['members'][0]['accepted_state_times'].pop(1)
        self.denied(q)

    def test_foreign_budget_binding(self):
        q = fixture();next(iter(q['audits'].values()))['budget_binding']['model'] = 'A93'
        self.denied(q)

    def test_native_only_is_not_complete_quadrature(self):
        q = fixture();first = next(iter(q['audits'].values()));first['status'] = 'pass'
        del first['budgets']['q64f32']
        self.denied(q)

    def test_missing_equation_and_nonfinite(self):
        q = fixture();first = next(iter(q['audits'].values()))
        del first['budgets']['native']['equations']['outer_epsilon']
        self.denied(q)
        q = fixture();next(iter(q['audits'].values()))['mass']['outward_net'] = float('nan')
        self.denied(q)

    def test_reconstruction_does_not_replace_actual_stage_binding(self):
        q = fixture();first = next(iter(q['audits'].values()));first['status'] = 'pass'
        first['captured_k_matrix']['maximum_defect'] = 3.306903700930181e-6
        self.denied(q)

    def test_shared_prerequisites_cannot_be_bypassed(self):
        for key in fixture()['prerequisites']:
            q = fixture();q['prerequisites'][key] = 'fail'
            with self.subTest(key=key):self.denied(q)

    def test_governing_corrections_forbidden(self):
        q = fixture();q['members'][0]['corrections_and_additional_solves'] = 1
        self.denied(q)

    def test_independent_traction_cannot_be_replaced_by_field_accuracy(self):
        q = fixture();q['members'][0]['independent_endpoint_metrics']['traction']['inner_interface_reference_error_scaled_RMS'] = .006
        self.denied(q)

    def test_command_exit_matches_primitives(self):
        script = Path(__file__).with_name('ct_accepted6_gate.py')
        with tempfile.TemporaryDirectory() as d:
            p = Path(d);value = fixture()
            for passed in (True, False):
                if not passed:value['members'].pop()
                (p/'input.json').write_text(json.dumps(value))
                run = subprocess.run([sys.executable, '-B', str(script), '--input', str(p/'input.json'),
                                      '--output', str(p/'result.json')], check=False)
                self.assertEqual(run.returncode, 0 if passed else 1)
                self.assertTrue((p/'result.json').is_file())
            value = fixture()
            next(iter(value['audits'].values()))['budgets']['native']['equations']['outer_epsilon']['storage_backward_difference'] = float('nan')
            (p/'input.json').write_text(json.dumps(value))
            run = subprocess.run([sys.executable, '-B', str(script), '--input', str(p/'input.json'),
                                  '--output', str(p/'nonfinite.json')], check=False)
            self.assertEqual(run.returncode, 1)
            result = json.loads((p/'nonfinite.json').read_text())
            self.assertEqual(result['status'], 'not_admitted')


if __name__ == '__main__':
    unittest.main()
