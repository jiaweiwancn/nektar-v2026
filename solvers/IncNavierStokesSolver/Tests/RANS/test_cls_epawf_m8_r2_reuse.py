"""Metadata-only rejection tests; no solver/numerical library is imported."""
import copy
import ast
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch

import cls_epawf_m8_r2_reuse as reuse


class ReuseMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='rans-r2-reuse-metadata-')
        self.root = Path(self.temporary.name).resolve()
        self.old, self.current, self.case, self.aux = (self.root/name for name in ('old', 'current', 'case', 'auxiliary-analysis'))
        for path in (self.old, self.current, self.case, self.aux):
            path.mkdir()
        self.code = 'import math\nVALUE=1\nclass Analyzer:\n def value(self): return VALUE\ndef utility(): return 1\ndef verify(): return "old"\ndef main(): pass\n'
        for source in (self.old, self.current):
            (source/reuse.HELPER).parent.mkdir(parents=True)
            (source/reuse.HELPER).write_text(self.code if source == self.old else self.code.replace('"old"', '"new"'))
            (source/'oracle.py').write_text('COEFFICIENT = 3\n')
        self.wrapper = self.old/'readonly-wrapper.py'; self.wrapper.write_text('# historical read-only wrapper\n')
        self.xml = self.case/'case.xml'; self.xml.write_text('<NEKTAR/>\n')
        self.preflight = self.case/'preflight.json'
        self.raw = dict(initial=[1., 2.], fields=[dict(size=2)])
        self.preflight.write_text(json.dumps(self.raw))
        self.matrix = self.case/'case.feasible.jsonl.solve1.J0.bin'; self.matrix.write_bytes(struct.pack('<4d', 1., 0., 0., 1.))
        self.record = dict(alpha=3200., time=.0025, history=[1., 2.], state=[1., 2.],
            inner_stress_model='Linear', policy='ScalarFeasibleR2:v1', frozen_velocity_pressure_diagnostic=False,
            iterations=[dict(jacobian_binary_file=self.matrix.name, jacobian_binary_shape=[2, 2])])
        self.trace = self.case/'case.feasible.jsonl'
        self.line = (json.dumps(self.record)+'\n').encode(); self.trace.write_bytes(self.line)
        self.snapshot = self.aux/'native-record-0.raw.jsonl'; self.snapshot.write_bytes(self.line)
        self.pin = dict(success=True, clean=True, source_commit='producer',
            **{key:{'identity':'unchanged'} for key in reuse.RUNTIME_KEYS})
        self.pin_path = self.case/'runtime-before.json'; self.pin_path.write_text(json.dumps(self.pin))
        self.solver = self.case/'solver'; self.solver.write_bytes(b'pinned test image')
        self.old_imports = {str(path):reuse.digest(path) for path in
            (self.old/reuse.HELPER, self.old/'oracle.py', self.wrapper)}
        files = {str(path):reuse.file_identity(path) for path in
            (self.preflight, self.xml, self.matrix, self.pin_path, self.snapshot, self.old/reuse.HELPER)}
        self.result = dict(success=True, analysis_success=True, failures=[], native_calls=0,
            only_one_completed_native_record_analyzed=True, selected_alpha=3200., physical_stage_time=.0025,
            solved_record=dict(alpha=3200., terminal='bounded_failure'), checks=42, errors={'retained_error':1e-6},
            source_bound_h0_float64_LE_sha256=reuse.vector_digest(self.raw['initial']),
            native_field_descriptors_sha256=reuse.canonical(self.raw['fields']), protected_physics_sha256='physics',
            dense_matrix_files={str(self.matrix):dict(files[str(self.matrix)], shape=[2, 2], dtype='float64', byteorder=sys.byteorder)},
            current_action_gates=[dict(available=True, gate_met=False)],
            exact_affine_weak_continuity_gates=[dict(gate_met=True)])
        self.receipt = dict(success=True, failures=[], native_calls=0, source_clean_start=True,
            source_clean_end=True, source_commit='old', source_commit_end='old', record_index=0,
            command_argv=[str(self.wrapper), '--source', str(self.old), '--preflight', str(self.preflight),
                '--stage-session', str(self.xml), '--runtime-pin', str(self.pin_path)],
            wrapper_sha256=reuse.digest(self.wrapper), inputs=files, imported_source_files_start=self.old_imports,
            active_trace_prefix=dict(path=str(self.trace), prefix_bytes=len(self.line),
                prefix_sha256=reuse.digest(self.trace), prefix_sha256_end=reuse.digest(self.trace),
                snapshot=str(self.snapshot), selected_bytes=len(self.line), selected_sha256=reuse.digest(self.snapshot)),
            producing_source_commit='producer', producing_runtime_sections=reuse.runtime_identity(self.pin),
            producing_solver=dict(path=str(self.solver), sha256=reuse.digest(self.solver)))
        self.write_old()
        self.git_patch = patch.object(reuse, 'git', side_effect=lambda source, *args:
            ('old' if Path(source) == self.old else 'current') if args == ('rev-parse', 'HEAD') else '')
        self.git_patch.start()

    def tearDown(self):
        self.git_patch.stop(); self.temporary.cleanup()

    def write_old(self):
        (self.aux/'receipt.json').write_text(json.dumps(self.receipt))
        (self.aux/'validation.json').write_text(json.dumps(self.result))

    def seal(self):
        return reuse.seal(self.aux, 'current')

    def context(self):
        imported = {str(path):reuse.digest(path) for path in (self.current/reuse.HELPER, self.current/'oracle.py')}
        return reuse.context(self.trace, [copy.deepcopy(self.record)], self.raw, self.xml,
            self.current, self.pin, imported, 'physics', self.case)

    def rejects(self, change):
        self.seal(); current = self.context(); change(current)
        with self.assertRaises(ValueError):
            reuse.load_reuse([self.aux], current)

    def test_reuses_failed_scientific_gate_without_promoting_it(self):
        self.seal(); result = reuse.load_reuse([self.aux], self.context())[0]
        self.assertFalse(result['result']['current_action_gates'][0]['gate_met'])
        self.assertEqual(result['result']['errors']['retained_error'], 1e-6)
        self.assertEqual(result['provenance']['original_checks'], 42)

    def test_missing_seal_is_not_a_fresh_compute_fallback(self):
        with self.assertRaises(ValueError):
            reuse.load_reuse([self.aux], self.context())

    def test_unsuccessful_receipt_cannot_be_sealed(self):
        self.receipt['success'] = False; self.write_old()
        with self.assertRaises(ValueError):
            self.seal()

    def test_unsuccessful_analysis_cannot_be_sealed(self):
        self.result['analysis_success'] = False; self.write_old()
        with self.assertRaises(ValueError):
            self.seal()

    def test_existing_seal_is_not_overwritten(self):
        self.seal()
        with self.assertRaises(ValueError):
            self.seal()

    def test_changed_raw_line_rejected(self):
        self.rejects(lambda c:c['raw_line_sha256'].__setitem__(0, 'different'))

    def test_changed_canonical_record_rejected_even_with_old_raw_hash(self):
        self.rejects(lambda c:c['records'][0].__setitem__('alpha', 400.))

    def test_changed_xml_rejected(self):
        self.rejects(lambda c:c.__setitem__('stage_xml_sha256', 'different'))

    def test_changed_protected_physics_rejected(self):
        self.rejects(lambda c:c.__setitem__('protected_physics_sha256', 'different'))

    def test_changed_fields_rejected(self):
        self.rejects(lambda c:c.__setitem__('fields_sha256', 'different'))

    def test_changed_h0_rejected(self):
        self.rejects(lambda c:c.__setitem__('h0_sha256', 'different'))

    def test_changed_runtime_sections_rejected(self):
        self.rejects(lambda c:c['runtime_sections'].__setitem__('runtime_bindings', 'different'))

    def test_changed_imported_oracle_bytes_rejected(self):
        self.rejects(lambda c:(self.current/'oracle.py').write_text('COEFFICIENT = 4\n'))

    def test_changed_Analyzer_ast_rejected(self):
        self.rejects(lambda c:(self.current/reuse.HELPER).write_text(self.code.replace('return VALUE', 'return VALUE+1')))

    def test_changed_existing_utility_ast_rejected(self):
        self.rejects(lambda c:(self.current/reuse.HELPER).write_text(self.code.replace('return 1', 'return 2')))

    def test_changed_constant_rejected(self):
        self.rejects(lambda c:(self.current/reuse.HELPER).write_text(self.code.replace('VALUE=1', 'VALUE=2')))

    def test_changed_existing_import_ast_rejected(self):
        self.rejects(lambda c:(self.current/reuse.HELPER).write_text(self.code.replace('import math', 'import cmath as math')))

    def test_missing_imported_oracle_rejected(self):
        self.rejects(lambda c:c['imported'].pop(str(self.current/'oracle.py')))

    def test_changed_matrix_rejected(self):
        self.rejects(lambda c:self.matrix.write_bytes(b'corrupt'))

    def test_incompatible_matrix_shape_rejected(self):
        self.result['dense_matrix_files'][str(self.matrix)]['shape'] = [3, 3]; self.write_old(); self.seal()
        with self.assertRaises(ValueError):
            reuse.load_reuse([self.aux], self.context())

    def test_incompatible_matrix_byteorder_rejected(self):
        self.result['dense_matrix_files'][str(self.matrix)]['byteorder'] = 'big' if sys.byteorder == 'little' else 'little'
        self.write_old(); self.seal()
        with self.assertRaises(ValueError):
            reuse.load_reuse([self.aux], self.context())

    def test_changed_validation_output_rejected(self):
        self.rejects(lambda c:(self.aux/'validation.json').write_text('{}'))

    def test_duplicate_record_cache_rejected(self):
        self.seal()
        with self.assertRaises(ValueError):
            reuse.load_reuse([self.aux, self.aux], self.context())

    def test_orchestration_changes_are_allowed_while_scientific_ast_matches(self):
        self.assertEqual(reuse.scientific_core(self.old/reuse.HELPER), reuse.scientific_core(self.current/reuse.HELPER))

    def test_reused_solve_skips_only_solve_and_keeps_final_protocol_audits(self):
        # Execute only the real read-only orchestration AST with metadata
        # stubs. No scientific module, mass matrix or native code is loaded.
        path = Path(reuse.__file__).resolve().parents[4]/reuse.HELPER
        node = next(node for node in ast.parse(path.read_text()).body if isinstance(node, ast.FunctionDef) and node.name == 'verify')
        calls = []
        class AuditStub:
            checks = 3
            errors = {'fresh_error':0.}
            def require(self, condition, message):
                if not condition:
                    raise ValueError(message)
        class AnalyzerStub:
            def __init__(self, *arguments):
                self.audit = AuditStub(); self.failures = []; self.gates = []; self.continuity_gates = []
                self.matrix_files = {}; self.static_coverage = {'fresh':True}
            def protocol(self, records):
                calls.append('protocol'); return {'fresh':True}
            def attempt(self, label, function):
                return function()
            def domain(self, *arguments, **keywords):
                calls.append('initial_domain'); return {'fresh':True}
            def solve(self, record):
                raise AssertionError('An explicitly reused solve was recomputed')
            def levels(self, records):
                calls.append('actual_final_levels'); return [{'fresh_postqualification':True}]
        namespace = dict(Analyzer=AnalyzerStub, ACTION_TOLERANCE=1e-7)
        exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])), str(path), 'exec'), namespace)
        raw = dict(initial=[], initial_scalar_domain={}, stage_test_target_converged=False, stage_test_target_physical_range=False)
        bundle = dict(result=self.result, provenance={'original_checks':42})
        result = namespace['verify'](raw, [{}], 'unused.xml', 'unused-source', reuse_solves={0:bundle})
        self.assertEqual(calls, ['protocol', 'initial_domain', 'actual_final_levels'])
        self.assertFalse(result['current_action_gate_passed'])
        self.assertTrue(result['full_final_case_protocol_and_postqualification_recomputed'])
        self.assertEqual(result['checks'], 3)
        self.assertEqual(result['solves'][0]['terminal'], 'bounded_failure')


if __name__ == '__main__':
    unittest.main()
