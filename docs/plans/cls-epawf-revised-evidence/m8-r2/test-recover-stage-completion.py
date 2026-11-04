"""Synthetic completion-recovery metadata checks; no native/scientific calls."""
import copy
import importlib.util
import json
from pathlib import Path
import struct
import subprocess
import sys
import unittest

sys.dont_write_bytecode = True
DIRECTORY = Path(__file__).parent


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, DIRECTORY/filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


base = module('stage_metadata_fixture', 'test-assess-stage-evidence.py')
recovery = module('stage_recovery', 'recover-stage-completion.py')
evidence = recovery.evidence_module()


class RecoveryMetadataTests(unittest.TestCase):
    def setUp(self):
        self.fixture = base.EvidenceMetadataTests()
        self.fixture.setUp()
        self.addCleanup(self.fixture.tearDown)
        f = self.fixture
        self.matched_assessment = f.stored_assessment()
        gate = f.gate(self.matched_assessment)
        self.run = f.root/'frozen-scalar'
        self.case = self.run/'cases/outer6-saturatedr2-frozenscalarr2'
        self.case.mkdir(parents=True)
        self.check_directory = self.run/'completion-check-001'
        self.check_directory.mkdir()
        self.xml = self.case/'case.xml'; self.xml.write_text(f.xml('SaturatedR2', 'FrozenScalarR2'))
        matrix = self.case/'case.feasible.jsonl.solve1.J0.bin'
        matrix.write_bytes(struct.pack('<4d', 1., 0., 0., 1.))
        record = f.record('SaturatedR2', 400., False, [dict(jacobian_binary_file=matrix.name,
            jacobian_binary_shape=[2, 2])])
        record.update(frozen_velocity_pressure_diagnostic=True, terminal='feasible_residual_stagnation',
            scaled_block_norms=[8., 15.])
        scalar_domain = dict(admissible=True, phase='complete', inventory_counts={'metadata_fixture':2},
            minimum_dimensionless_margin=.1, minimum_b=1., minimum_lambda=1., profile_branch_pattern='BP')
        record['domain'] = scalar_domain
        record['iterations'][0].update(iteration=0, state=f.initial, domain=scalar_domain,
            columns=[dict(stencil='metadata_fixture')],
            trials=[dict(fraction=.5, classification='accepted', domain=scalar_domain)])
        self.trace = self.case/'case.feasible.jsonl'
        self.trace.write_text(json.dumps(record)+'\n')
        self.report = self.case/'case.evolution.json'
        self.native = dict(initial=f.initial, state=f.initial, fields=[dict(size=2)], inner_stress_model='SaturatedR2',
            signature='synthetic metadata signature',
            stage_test='FrozenScalarR2', stage_test_target_converged=False, stage_test_target_physical_range=False,
            time=0, steps=0, accepted_states=[f.initial], native_history=[f.initial], native_history_times=[0],
            stage_test_initial_native_history=dict(time=0, steps=0, state=f.initial, history=[f.initial], times=[0]),
            stage_test_levels=[dict(alpha=400., time=.0025, attempted=True, state=f.initial, history=f.initial,
                native_history_advanced=False, converged=False, true_coupled_residual_converged=False,
                physical_range_passed=False, scaled_block_norms=[8., 15.], qualification_audit_failure='',
                iterations=0, residual=[0., 0.], constraints=[], scalar_domain=scalar_domain,
                audit={key:[] for key in ('energy', 'interface', 'strong_continuity', 'r2_qualified',
                    'r2_outer_range', 'r2_covariance', 'r2_scalar_domain')})])
        f.write(self.report, self.native)
        self.log = self.run/'native-SaturatedR2.log'
        self.log.write_text('EP-AWF ScalarFeasibleR2 fresh Jacobian iteration 0 alpha 400 residual 15\n'
            'EP-AWF R2 stage test alpha 400 converged 0 physical range 0 history unchanged\n')
        self.before = self.run/'runtime-before.json'
        f.write(self.before, f.pin); f.write(self.run/'runtime-after.json', f.pin)
        fp = f.source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
        fp.parent.mkdir(parents=True); fp.write_text('# metadata fixture; never executed\n')
        self.post = self.check_directory/'runtime-postcompletion.json'
        self.post_data = dict(copy.deepcopy(f.pin), captured_utc='2026-11-04T01:00:01+00:00',
            compared_fingerprint=str(self.before), compared_fingerprint_sha256=evidence.sha(self.before),
            identical_qualification_inputs_and_runtime=True, differences={})
        f.write(self.post, self.post_data)
        plan = dict(copy.deepcopy(f.plan), alpha_levels=[400.], cases=[dict(inner_stress='SaturatedR2',
            xml=str(self.xml), xml_sha256=evidence.sha(self.xml), representation='outer6', cost=dict(unknowns=2))])
        plan_path = self.run/'cases/stage-plan.json'; f.write(plan_path, plan)
        self.raw = copy.deepcopy(f.receipt)
        for key in ('finished_utc', 'source_commit_end', 'source_clean_end', 'runtime_unchanged'):
            self.raw.pop(key)
        self.raw.update(status='failed', mode='frozen-scalar', alpha_levels=[400.],
            failure=dict(type='BrokenPipeError', message='[Errno 32] Broken pipe'),
            matched_receipt=dict(path=str(f.receipt_path), sha256=evidence.sha(f.receipt_path)), matched_evidence_gate=gate,
            stage_plan_sha256=evidence.sha(plan_path),
            cases=[dict(name=self.case.name, inner_stress='SaturatedR2', xml_path=str(self.xml),
                xml_sha256=evidence.sha(self.xml), status='started')], commands=[])
        for name in ('fingerprint-before', 'prepare', 'native-SaturatedR2', 'fingerprint-after'):
            if name.startswith('native-'):
                self.raw['commands'].append(dict(name=name, argv=['taskset', '-c', '0',
                    str(f.build/'IncNavierStokesSolver'), str(self.xml)], cwd=str(self.case), timeout_seconds=14400, pid=99999999))
            else:
                log = self.run/(name+'.log'); log.write_text('completed metadata fixture\n')
                self.raw['commands'].append(dict(name=name, exit_code=0, log_path=str(log), log_sha256=evidence.sha(log)))
        self.raw_path = self.run/'receipt.json'; f.write(self.raw_path, self.raw)
        fp_log = self.check_directory/'fingerprint.log'; fp_log.write_text('successful fresh fingerprint fixture\n')
        self.check = dict(schema=1, success=True, native_calls=0,
            started_utc='2026-11-04T01:00:00+00:00', finished_utc='2026-11-04T01:00:02+00:00',
            numerical_threads={key:'1' for key in evidence.THREAD_KEYS}, cpu_affinity=[50],
            source_commit_start='native-pin', source_commit_end='native-pin', source_clean_start=True, source_clean_end=True,
            completed_report_present_before_fingerprint=True, no_producing_solver_process_observed=True,
            artifacts_unchanged=True, native_process_exit_code=None, native_process_elapsed_seconds=None,
            preexisting_completed_artifacts=[evidence.identity(path) for path in
                (self.raw_path, self.before, self.log, self.report, self.trace, self.xml)],
            postcompletion_runtime=evidence.identity(self.post),
            fingerprint_command=dict(argv=['taskset', '-c', '50', 'python3', '-B', str(fp), '--root', str(f.root),
                '--output', str(self.post), '--compare', str(self.before)], exit_code=0, log=evidence.identity(fp_log)))
        self.check_path = self.check_directory/'receipt.json'; self.save_check()

    def save_check(self):
        self.fixture.write(self.check_path, self.check)

    def refresh_artifact(self, path):
        self.check['preexisting_completed_artifacts'] = [evidence.identity(path) if pin['path'] == str(path) else pin
            for pin in self.check['preexisting_completed_artifacts']]
        self.save_check()

    def collect(self):
        return recovery.collect(self.raw_path, self.check_path, 'native-pin')

    def stored_recovery(self):
        path = self.run/'additive-recovery.json'
        self.fixture.write(path, dict(schema=1, kind=recovery.KIND, success=True,
            collector=evidence.identity(DIRECTORY/'recover-stage-completion.py'),
            assessor=evidence.identity(DIRECTORY/'assess-stage-evidence.py'),
            native_calls=0, scientific_reconstructions=0, failures=[], recovery=self.collect()))
        return path

    def test_complete_recovery_preserves_failed_supervisor_and_unknown_process_status(self):
        original = self.raw_path.read_bytes()
        result = self.collect()
        self.assertEqual(self.raw_path.read_bytes(), original)
        self.assertEqual(result['original_supervisor_status'], 'failed')
        self.assertIsNone(result['native_process_exit_code'])
        self.assertIsNone(result['native_process_elapsed_seconds'])
        self.assertFalse(result['native_wall_cap_enforcement_after_supervisor_failure_verified'])
        self.assertEqual(len(result['matrices']), 1)
        self.assertEqual(result['complete_trace']['complete_records'][0]['iteration_count'], 1)

    def test_assessor_recovery_keeps_scalar_noncoupled_failure_separate(self):
        answer = evidence.assess(self.raw_path, 'native-pin', self.stored_recovery())
        self.assertEqual(answer['family_status'], 'recovered_native_artifacts')
        self.assertEqual(answer['original_supervisor_status'], 'failed')
        self.assertFalse(answer['cases'][0]['scalar_solve_convergence_evidence'])
        self.assertFalse(answer['cases'][0]['full_coupled_convergence_evidence'])
        self.assertIsNone(answer['cases'][0]['wall_cap_stop_observed'])
        self.assertIsNone(answer['cases'][0]['native_command']['exit_code'])

    def test_normal_assessor_still_rejects_unfinished_raw_receipt(self):
        with self.assertRaisesRegex(ValueError, 'terminal stages-family'):
            evidence.assess(self.raw_path, 'native-pin')

    def test_missing_fresh_post_completion_pin_rejected(self):
        self.post.unlink()
        with self.assertRaises(OSError): self.collect()

    def test_no_completed_artifact_observation_rejected(self):
        self.check['completed_report_present_before_fingerprint'] = False; self.save_check()
        with self.assertRaisesRegex(ValueError, 'fresh post-completion'): self.collect()

    def test_changed_raw_receipt_rejected(self):
        self.raw['extra'] = 'changed after completion observation'; self.fixture.write(self.raw_path, self.raw)
        with self.assertRaisesRegex(ValueError, 'artifact changed'): self.collect()

    def test_changed_native_artifact_rejected(self):
        self.trace.write_text(self.trace.read_text()+'{}\n')
        with self.assertRaisesRegex(ValueError, 'artifact changed'): self.collect()

    def test_incomplete_jsonl_target_rejected_even_with_current_hash(self):
        self.trace.write_text(self.trace.read_text().rstrip('\n')); self.refresh_artifact(self.trace)
        with self.assertRaisesRegex(ValueError, 'complete frozen target'): self.collect()

    def test_mismatched_complete_final_stage_rejected(self):
        self.native['stage_test_levels'][0]['state'] = [9., 2.]
        self.fixture.write(self.report, self.native); self.refresh_artifact(self.report)
        with self.assertRaisesRegex(ValueError, 'terminal record differ'): self.collect()

    def test_fake_native_exit_in_completion_check_rejected(self):
        self.check['native_process_exit_code'] = 0; self.save_check()
        with self.assertRaisesRegex(ValueError, 'remain null'): self.collect()

    def test_fake_native_exit_in_raw_command_rejected(self):
        self.raw['commands'][2]['exit_code'] = 0; self.fixture.write(self.raw_path, self.raw)
        with self.assertRaisesRegex(ValueError, 'exit status/timing'): self.collect()

    def test_fake_native_exit_in_additive_receipt_rejected(self):
        path = self.stored_recovery(); stored = evidence.load(path)
        stored['recovery']['native_command']['exit_code'] = 0; self.fixture.write(path, stored)
        with self.assertRaisesRegex(ValueError, 'unknown process status changed'):
            evidence.assess(self.raw_path, 'native-pin', path)

    def test_fingerprint_captured_before_completed_artifact_check_rejected(self):
        self.post_data['captured_utc'] = '2026-11-04T00:00:00+00:00'
        self.fixture.write(self.post, self.post_data)
        self.check['postcompletion_runtime'] = evidence.identity(self.post); self.save_check()
        with self.assertRaisesRegex(ValueError, 'fresh completed-artifact check'): self.collect()

    def test_curator_preserves_failed_supervisor_fresh_runtime_and_unknown_status(self):
        stored = self.stored_recovery(); output = self.run/'compact-native-summary.json'
        completed = subprocess.run([sys.executable, '-B', str(DIRECTORY/'curate-stages.py'),
            '--run', str(self.run), '--recovery', str(stored), '--output', str(output)],
            capture_output=True, text=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)
        answer = evidence.load(output)
        self.assertEqual(answer['status'], 'recovered_native_artifacts')
        self.assertEqual(answer['original_supervisor_status'], 'failed')
        self.assertEqual(answer['runtime_after'], evidence.identity(self.post))
        self.assertEqual(answer['producer_source_end'], 'native-pin')
        self.assertIsNone(answer['cases'][0]['native_exit_code'])
        self.assertIsNone(answer['cases'][0]['native_seconds'])
        self.assertFalse(answer['native_wall_cap_enforcement_after_supervisor_failure_verified'])
        self.assertTrue(answer['strict_recovered_assessment']['gates_recomputed'])
        self.assertEqual(answer['cases'][0]['nonlinear_solves'][0]['native_iteration_entries'], 1)
        self.assertEqual(answer['cases'][0]['nonlinear_solves'][0]['trial_classifications'], {'accepted':1})
        self.assertEqual(evidence.load(self.raw_path)['status'], 'failed')

    def test_curator_recomputes_recovery_and_rejects_fake_process_status(self):
        stored = self.stored_recovery(); value = evidence.load(stored)
        value['recovery']['native_process_exit_code'] = 0; self.fixture.write(stored, value)
        output = self.run/'must-not-be-written.json'
        completed = subprocess.run([sys.executable, '-B', str(DIRECTORY/'curate-stages.py'),
            '--run', str(self.run), '--recovery', str(stored), '--output', str(output)],
            capture_output=True, text=True, timeout=30)
        self.assertNotEqual(completed.returncode, 0)
        self.assertFalse(output.exists())
        self.assertIn('unknown process status changed', completed.stderr)


if __name__ == '__main__':
    unittest.main()
