#!/usr/bin/env python3
"""Focused stdlib tests for R2 ledger metadata semantics; no numerical calls."""
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import hashlib
import importlib.util
from io import StringIO
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.dont_write_bytecode = True
SCRIPT = Path(__file__).with_name('curate-decision.py').resolve()
SPEC = importlib.util.spec_from_file_location('r2_decision_curator_metadata_tests', SCRIPT)
CURATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CURATOR)


def unicode_hash(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False,
        separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def pin(path):
    return dict(path=path, sha256='a'*64)


def fixture():
    """Explicit terminal incomplete evidence with unavailable root status."""
    runtime = {key: {'unicode_path': '新加卷'} for key in CURATOR.RUNTIME_SECTIONS}
    runtime.update(success=True, clean=True, source_commit=CURATOR.PRODUCER,
        source_root='/producer/source', build_root='/producer/build',
        runtime_artifacts={'build/IncNavierStokesSolver': 'solver-sha'})

    def case(law):
        return dict(inner_stress=law, preflight_h0_sha256='f'*64,
            target_attempted=True, target_start_observed_in_log=True, complete_target_solve_record=False,
            converged_target_evidence=None, full_coupled_convergence_evidence=None,
            scalar_solve_convergence_evidence=None, wall_cap_stop_observed=True,
            issues=[], xml=pin('/'+law+'/case.xml'), log=pin('/'+law+'/native.log'),
            trace_evidence=dict(complete_records=[dict(alpha=3200., converged=False)],
                complete_prefix_bytes=100, complete_prefix_sha256='b'*64, unpublished_tail=None),
            native_postqualification_report_available=False,
            native_history_report_available=False, native_history_unchanged_verified=None)

    wrapper = dict(success=True, failures=[], native_calls=0, assessor=pin('/assessor.py'),
        assessment=dict(source_commit=CURATOR.PRODUCER, mode='matched', family_status='incomplete',
            original_receipt=pin('/receipt.json'), preflight=pin('/preflight.json'),
            original=pin('/original.xml'), stage_plan=pin('/stage-plan.json'),
            runtime=dict(source_commit=CURATOR.PRODUCER, source_root='/producer/source',
                build_root='/producer/build', sections={key: unicode_hash(runtime[key])
                    for key in CURATOR.RUNTIME_SECTIONS}, solver=dict(sha256='solver-sha')),
            cases=[case('Linear'), case('SaturatedR2')],
            bound_files={'/matrix.bin': pin('/matrix.bin')}, conditional_frozen_scalar_allowed=True))
    return runtime, wrapper


class AssessmentEvidence:
    """Supply already pinned metadata; failed binding behavior has a real-file test."""
    def __init__(self, wrapper, receipt_status='incomplete'):
        self.issues = []
        self.wrapper = wrapper
        self.receipt_status = receipt_status

    def add(self, *args, **kwargs):
        return self.wrapper

    def reference(self, role, *args, **kwargs):
        if role == 'strict_stage_original_receipt':
            return dict(status=self.receipt_status, source_commit=CURATOR.PRODUCER,
                finished_utc='terminal-metadata-fixture')
        return {}


class AssessmentTests(unittest.TestCase):
    def setUp(self):
        self.runtime, self.wrapper = fixture()

    def assess(self, wrapper=None, receipt_status='incomplete'):
        evidence = AssessmentEvidence(self.wrapper if wrapper is None else wrapper, receipt_status)
        result, unknowns = CURATOR.stage_assessment(evidence, Path('/assessment.json'), self.runtime)
        return result, unknowns, evidence.issues

    def test_terminal_incomplete_observation_keeps_root_status_unknown(self):
        result, unknowns, issues = self.assess()
        self.assertTrue(result['bounded_evidence_ready'])
        self.assertFalse(issues)
        self.assertTrue(unknowns)
        self.assertIn('completion unavailable', result['bounded_task_outcome'])
        for case in result['cases']:
            self.assertIsNone(case['converged_target_evidence'])
            self.assertFalse(case['complete_target_solve_record'])
            self.assertFalse(case['independent_complete_target_analysis_supplied_by_assessment'])
        self.assertIsNot(result.get('physical_trajectory_qualification'), True)
        self.assertIsNot(result.get('final_stage_evidence_ready'), True)

    def test_positive_target_evidence_survives_missing_final_report(self):
        item = self.wrapper['assessment']['cases'][0]
        item.update(converged_target_evidence=True, full_coupled_convergence_evidence=True,
            complete_target_solve_record=True)
        item['trace_evidence']['complete_records'].append(dict(alpha=400., converged=True))
        self.wrapper['assessment']['conditional_frozen_scalar_allowed'] = False
        result, _, issues = self.assess()
        self.assertTrue(result['bounded_evidence_ready'])
        self.assertFalse(issues)
        self.assertIs(result['cases'][0]['converged_target_evidence'], True)
        self.assertIn('positive target-convergence evidence', result['bounded_task_outcome'])
        self.assertFalse(result['conditional_frozen_scalar_allowed'])
        self.assertFalse(result['fallback_execution_inferred'])

    def test_runtime_section_mismatch_rejects_bounded_ready(self):
        self.wrapper['assessment']['runtime']['sections']['source_inputs'] = 'wrong'
        result, _, issues = self.assess()
        self.assertFalse(result['bounded_evidence_ready'])
        self.assertTrue(issues)

    def test_active_family_is_not_a_terminal_bounded_observation(self):
        self.wrapper['assessment']['family_status'] = 'running'
        result, _, issues = self.assess(receipt_status='running')
        self.assertFalse(result['bounded_evidence_ready'])
        self.assertTrue(issues)

    def test_missing_actual_target_attempt_is_not_promoted(self):
        self.wrapper['assessment']['cases'][1]['target_attempted'] = None
        result, _, issues = self.assess()
        self.assertFalse(result['bounded_evidence_ready'])
        self.assertIsNone(result['cases'][1]['target_attempted'])
        self.assertTrue(issues)

    def test_reported_case_issue_blocks_bounded_ready(self):
        self.wrapper['assessment']['cases'][1]['issues'] = ['unpublished JSONL tail']
        result, _, issues = self.assess()
        self.assertFalse(result['bounded_evidence_ready'])
        self.assertTrue(issues)

    def test_assessor_hash_is_required(self):
        self.wrapper['assessor'].pop('sha256')
        result, _, issues = self.assess()
        self.assertFalse(result['bounded_evidence_ready'])
        self.assertTrue(any(item['kind'] == 'strict_stage_required_reference_hash_missing' for item in issues))


class LedgerInterfaceTests(unittest.TestCase):
    def frozen_family(self, scalar_converged):
        runtime, _ = fixture()
        summary = dict(mode='frozen-scalar', producer_source=CURATOR.PRODUCER,
            producer_source_end=CURATOR.PRODUCER, source_clean=True, runtime_unchanged=True,
            cases=[dict(inner_stress='SaturatedR2', target_attempted=True,
                target_true_coupled_convergence=False, target_native_physical_range_flag=False,
                native_exit_code=None, native_seconds=None,
                levels=[dict(alpha=400., attempted=True, converged=scalar_converged,
                    true_coupled_residual_converged=False, frozen_velocity_pressure_diagnostic=True,
                    scalar_domain=dict(admissible=True))])])

        class SummaryEvidence:
            def add(self, *args, **kwargs):
                return summary

            def reference(self, *args, **kwargs):
                return {}

        with mock.patch.object(CURATOR, 'runtime_binding', return_value=True):
            return CURATOR.stage_family(SummaryEvidence(), Path('/scalar-summary.json'), [], 'scalar', runtime)

    def test_complete_failed_scalar_record_is_reported_as_diagnostic_without_path_audit(self):
        result = self.frozen_family(False)
        row = result['cases'][0]
        self.assertTrue(row['target_record_present'])
        self.assertFalse(row['target']['converged'])
        self.assertIn('frozen-scalar diagnostic only', row['result_from_available_records'])
        self.assertEqual(row['report_bound_analysis_count'], 0)
        self.assertIsNone(row['native_exit_code'])
        self.assertIsNone(row['native_seconds'])
        self.assertFalse(row['all_supplied_metadata_qualification_gates'])

    def test_scalar_convergence_without_coupled_solution_is_never_promoted(self):
        row = self.frozen_family(True)['cases'][0]
        self.assertTrue(row['target']['converged'])
        self.assertFalse(row['native_target_true_coupled_convergence'])
        self.assertTrue(row['frozen_scalar_is_not_a_coupled_root'])
        self.assertFalse(row['all_supplied_metadata_qualification_gates'])
        self.assertIn('remain separate', row['result_from_available_records'])

    def test_actual_reference_hash_mismatch_is_recorded(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'receipt.json'
            path.write_text(json.dumps(dict(status='incomplete'))+'\n')
            evidence = CURATOR.Evidence()
            evidence.reference('original_receipt', dict(path=str(path), sha256='0'*64), required=True)
            self.assertTrue(any(item['kind'] == 'pinned_evidence_hash_mismatch' for item in evidence.issues))

    def test_bounded_final_payload_does_not_supply_complete_or_physical_root(self):
        runtime, _ = fixture()
        assessment = dict(supplied=True, bounded_evidence_ready=True,
            bounded_task_outcome='bounded target attempts observed; completion unavailable', cases=[])

        class MetadataOnlyEvidence:
            """Isolate output wiring; assessor validation is tested above."""
            def __init__(self):
                self.files, self.roles, self.data, self.issues = {}, {}, {}, []

            def add(self, role, path, *args, **kwargs):
                path = str(Path(path).resolve())
                value = runtime if role == 'runtime' else (
                    dict(reviewed_source=CURATOR.REVIEWED) if Path(path).name == 'candidate-contract-v2.json' else {})
                self.roles[role] = path
                self.data[path] = value
                return value

        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); output = root/'ledger.json'
            argv = [str(SCRIPT), '--root', str(root), '--documents', str(root/'documents'),
                '--stage-assessment', str(root/'assessment.json'), '--bounded-final', '--output', str(output)]
            with mock.patch.object(sys, 'argv', argv), mock.patch.object(CURATOR, 'Evidence', MetadataOnlyEvidence), \
                    mock.patch.object(CURATOR, 'stage_family', return_value=dict(supplied=False, cases=[])), \
                    mock.patch.object(CURATOR, 'stage_assessment', return_value=(deepcopy(assessment), ['Target completion unavailable'])), \
                    redirect_stdout(StringIO()):
                self.assertEqual(CURATOR.main(), 0)
            value = json.loads(output.read_text())
            self.assertTrue(value['bounded_task_evidence_ready'])
            self.assertFalse(value['final_stage_evidence_ready'])
            self.assertFalse(value['physical_trajectory_qualification'])
            self.assertFalse(value['M8_closed'])
            self.assertFalse(value['fallback_inferred'])
            self.assertTrue(value['remaining_unknowns'])

    def test_final_and_bounded_final_are_mutually_exclusive(self):
        argv = [str(SCRIPT), '--root', '/unused/root', '--documents', '/unused/documents',
            '--output', '/unused/output.json', '--final', '--bounded-final']
        captured = StringIO()
        with mock.patch.object(sys, 'argv', argv), redirect_stderr(captured):
            with self.assertRaises(SystemExit) as stopped:
                CURATOR.main()
        self.assertEqual(stopped.exception.code, 2)
        self.assertIn('not allowed with argument', captured.getvalue())


if __name__ == '__main__':
    unittest.main()
