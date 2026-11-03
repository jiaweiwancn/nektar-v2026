"""Synthetic metadata-only checks; no solver, array package or oracle imports."""
import copy
import importlib.util
import json
from pathlib import Path
import struct
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
MODULE = Path(__file__).with_name('assess-stage-evidence.py')
SPEC = importlib.util.spec_from_file_location('stage_evidence', MODULE)
evidence = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evidence)


class EvidenceMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='rans-r2-stage-metadata-')
        self.root = Path(self.temporary.name).resolve()
        self.source, self.build = self.root/'source', self.root/'build'
        self.run, self.prior = self.root/'matched-stage', self.root/'preflight'
        for path in (self.source, self.build, self.run/'cases', self.prior):
            path.mkdir(parents=True)
        (self.source/'kernel.cpp').write_text('// unchanged test source\n')
        (self.build/'IncNavierStokesSolver').write_bytes(b'metadata fixture; never executed')
        tool = self.root/'python-tool'; tool.write_bytes(b'metadata-only tool fixture')
        self.pin = dict(success=True, clean=True, source_commit='native-pin',
            source_root=str(self.source), build_root=str(self.build),
            source_inputs={'kernel.cpp':evidence.sha(self.source/'kernel.cpp')},
            runtime_artifacts={'build/IncNavierStokesSolver':evidence.sha(self.build/'IncNavierStokesSolver')},
            runtime_resolved_paths={'build/IncNavierStokesSolver':str(self.build/'IncNavierStokesSolver')},
            runtime_bindings={'fixture':{}}, runtime_environment={'OMP_NUM_THREADS':'1'},
            python_packages={'fixture':{'version':'metadata-only'}},
            tool_images={'python3':dict(resolved=str(tool), sha256=evidence.sha(tool))})
        for name in ('runtime-before.json', 'runtime-after.json'):
            self.write(self.run/name, self.pin)
        self.runtime_pin = self.root/'runtime-fingerprint.json'; self.write(self.runtime_pin, self.pin)
        self.original = self.root/'original.xml'; self.original.write_text('<NEKTAR/>\n')
        self.initial = [1., 2.]
        self.preflight = dict(status='passed', full_runtime_source_inputs_unchanged=True,
            source_commit='native-pin', matched_initial_state_identity_all_quadratures=True, cases=[])
        self.plan = dict(alpha_levels=[3200., 1600., 800., 400.], target_alpha=400.,
            maximum_Newton_updates_per_level=12, time_and_native_history_advanced=False,
            original=dict(path=str(self.original), sha256=evidence.sha(self.original)), cases=[])
        self.receipt = dict(schema=2, status='incomplete', mode='matched', source_commit='native-pin',
            source_commit_end='native-pin', source_clean_start=True, source_clean_end=True,
            runtime_unchanged=True, finished_utc='2026-11-04T20:00:00+00:00',
            target_dt=.0025, target_time=.0025, maximum_Newton_updates_per_level=12,
            maximum_seconds_per_native_case=14400, ranks=1, maximum_concurrent_workers=1,
            numerical_threads={key:'1' for key in evidence.THREAD_KEYS}, representation='outer6',
            alpha_levels=[3200., 1600., 800., 400.],
            matrix_encoding=dict(dtype='IEEE754 binary64', byte_order=sys.byteorder, layout='row-major'),
            original_path=str(self.original), original_sha256=evidence.sha(self.original),
            commands=[], cases=[])
        for name in ('fingerprint-before', 'fingerprint-after', 'prepare'):
            log = self.run/(name+'.log'); log.write_text('metadata fixture\n')
            self.receipt['commands'].append(dict(name=name, exit_code=0, log_path=str(log), log_sha256=evidence.sha(log)))
        self.paths, self.records = {}, {}
        for law in evidence.LAWS:
            path = self.run/'cases'/law.lower(); path.mkdir()
            xml = path/'case.xml'; xml.write_text(self.xml(law, 'MatchedR2'))
            prior_xml = self.prior/(law.lower()+'.xml'); prior_xml.write_text(self.xml(law, 'Off'))
            report = self.prior/(law.lower()+'.json')
            self.write(report, dict(initial=self.initial, fields=[dict(size=2)], time=0, steps=0,
                accepted_states=[self.initial]))
            self.preflight['cases'].append(dict(inner_stress=law, quadrature='base', xml=str(prior_xml),
                xml_sha256=evidence.sha(prior_xml), report=str(report), report_sha256=evidence.sha(report),
                initial_state_float64_LE_sha256=evidence.vector_sha(self.initial)))
            self.plan['cases'].append(dict(inner_stress=law, xml=str(xml), xml_sha256=evidence.sha(xml),
                representation='outer6', cost=dict(unknowns=2)))
            matrix = path/'case.feasible.jsonl.solve1.J0.bin'
            matrix.write_bytes(struct.pack('<4d', 1., 0., 0., 1.))
            # A target matrix may be fully written before its solve record is.
            (path/'case.feasible.jsonl.solve2.J0.bin').write_bytes(matrix.read_bytes())
            self.records[law] = [self.record(law, 3200., False, [dict(jacobian_binary_file=matrix.name,
                jacobian_binary_shape=[2, 2])])]
            trace = path/'case.feasible.jsonl'
            self.write_trace(law, trace)
            log = self.run/('native-'+law+'.log')
            log.write_text('EP-AWF ScalarFeasibleR2 fresh Jacobian iteration 0 alpha 400 residual 13.1962\n')
            command = dict(name='native-'+law, argv=['taskset', '-c', '0', str(self.build/'IncNavierStokesSolver'), str(xml)],
                cwd=str(path), timeout_seconds=14400, exit_code=124, seconds=14401., log_path=str(log), log_sha256=evidence.sha(log))
            self.receipt['commands'].append(command)
            row = dict(inner_stress=law, xml_path=str(xml), xml_sha256=evidence.sha(xml), native_exit_code=124,
                native_seconds=14401., trace_path=str(trace), trace_sha256=evidence.sha(trace), status='native_incomplete')
            self.receipt['cases'].append(row)
            self.paths[law] = dict(directory=path, xml=xml, log=log, trace=trace, row=row, command=command)
        self.preflight_path = self.prior/'receipt.json'
        self.receipt_path = self.run/'receipt.json'
        self.save()

    def tearDown(self):
        self.temporary.cleanup()

    @staticmethod
    def write(path, value):
        Path(path).write_text(json.dumps(value, allow_nan=False)+'\n')

    @staticmethod
    def xml(law, stage):
        return ('<NEKTAR><CONDITIONS><SOLVERINFO><I PROPERTY="RANSWallInnerStress" VALUE="'+law+'"/>'
                '<I PROPERTY="RANSWallStageTest" VALUE="'+stage+'"/></SOLVERINFO></CONDITIONS></NEKTAR>\n')

    def record(self, law, alpha, converged, iterations=None):
        return dict(schema=1, policy='ScalarFeasibleR2:v1', inner_stress_model=law, time=.0025, alpha=alpha,
            tolerance=2e-12, maximum_iterations=12, converged=converged, terminal='recorded_terminal',
            frozen_velocity_pressure_diagnostic=False, history_advanced=False, residual_unchanged=True,
            state=self.initial, history=self.initial, residual=[0., 0.], scaled_block_norms=[0., 0.],
            iterations=iterations or [])

    def write_trace(self, law, path=None):
        path = path or self.paths[law]['trace']
        path.write_text(''.join(json.dumps(record)+'\n' for record in self.records[law]))
        if law in self.paths:
            self.paths[law]['row']['trace_sha256'] = evidence.sha(path)

    def save(self):
        self.write(self.preflight_path, self.preflight)
        self.receipt.update(preflight_path=str(self.preflight_path), preflight_sha256=evidence.sha(self.preflight_path))
        self.write(self.run/'cases/stage-plan.json', self.plan)
        self.receipt['stage_plan_sha256'] = evidence.sha(self.run/'cases/stage-plan.json')
        self.write(self.receipt_path, self.receipt)

    def answer(self):
        self.save()
        return evidence.assess(self.receipt_path, 'native-pin')

    def mutate_log(self, law, text):
        item = self.paths[law]
        item['log'].write_text(text)
        item['command']['log_sha256'] = evidence.sha(item['log'])

    def final(self, law, converged=False):
        item = self.paths[law]
        level = dict(alpha=400., time=.0025, attempted=True, history=self.initial,
            native_history_advanced=False, converged=converged, true_coupled_residual_converged=converged,
            qualification_audit_failure='', physical_range_passed=False,
            audit={key:[] for key in ('energy', 'interface', 'strong_continuity', 'r2_qualified',
                'r2_outer_range', 'r2_covariance', 'r2_scalar_domain')})
        raw = dict(initial=self.initial, fields=[dict(size=2)], inner_stress_model=law,
            stage_test='MatchedR2', stage_test_levels=[level], stage_test_target_converged=converged,
            stage_test_target_physical_range=False, state=self.initial, time=0, steps=0,
            accepted_states=[self.initial], native_history=[self.initial], native_history_times=[0],
            stage_test_initial_native_history=dict(time=0, steps=0, state=self.initial, history=[self.initial], times=[0]))
        path = item['directory']/'case.evolution.json'
        self.write(path, raw)
        item['row'].update(report_path=str(path), report_sha256=evidence.sha(path))
        return path, raw

    def stored_assessment(self):
        path = self.root/'assessment.json'
        self.write(path, dict(success=True, assessor=evidence.identity(MODULE), assessment=self.answer()))
        return path

    def gate(self, assessment=None, **changes):
        kwargs = dict(expected_source='native-pin', preflight=self.preflight_path, original=self.original,
            representation='outer6', runtime_pin=self.runtime_pin)
        kwargs.update(changes)
        return evidence.check_frozen_scalar(assessment, self.receipt_path, **kwargs)

    def test_cap_retains_unknown_completion_and_history(self):
        answer = self.answer()
        self.assertTrue(answer['conditional_frozen_scalar_allowed'])
        for case in answer['cases']:
            self.assertTrue(case['target_start_observed_in_log'])
            self.assertTrue(case['target_attempted'])
            self.assertFalse(case['complete_target_solve_record'])
            self.assertIsNone(case['converged_target_evidence'])
            self.assertIsNone(case['native_history_unchanged_verified'])
            self.assertFalse(case['native_postqualification_report_available'])
            self.assertTrue(case['wall_cap_stop_observed'])
            self.assertEqual(len(case['matrices']), 2)
            self.assertFalse(case['matrices'][1]['referenced_by_complete_record'])

    def test_absent_log_start_is_unknown_and_blocks_conditional(self):
        self.mutate_log('Linear', 'unrelated output only\n')
        answer = self.answer()
        self.assertIsNone(answer['cases'][0]['target_attempted'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_complete_converged_target_cannot_be_discarded_without_final_report(self):
        self.records['Linear'].append(self.record('Linear', 400., True)); self.write_trace('Linear')
        answer = self.answer()
        self.assertTrue(answer['cases'][0]['complete_target_solve_record'])
        self.assertTrue(answer['cases'][0]['converged_target_evidence'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_complete_target_record_establishes_attempt_without_log_start(self):
        self.mutate_log('Linear', 'no fresh-J line: root could converge at iteration zero\n')
        self.records['Linear'].append(self.record('Linear', 400., True)); self.write_trace('Linear')
        answer = self.answer()
        self.assertTrue(answer['cases'][0]['target_attempted'])
        self.assertFalse(answer['cases'][0]['target_start_observed_in_log'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_log_positive_summary_alone_blocks_fallback(self):
        self.mutate_log('Linear', 'EP-AWF R2 stage test alpha 400 converged 1 physical range 0 history unchanged\n')
        answer = self.answer()
        self.assertTrue(answer['cases'][0]['converged_target_evidence'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_conflicting_convergence_preserves_positive_evidence(self):
        self.records['Linear'].append(self.record('Linear', 400., True)); self.write_trace('Linear')
        self.final('Linear', False)
        answer = self.answer()
        self.assertTrue(answer['cases'][0]['converged_target_evidence'])
        self.assertTrue(answer['cases'][0]['issues'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_complete_nonconverged_cases_are_distinct_from_cap_unknowns(self):
        for law in evidence.LAWS:
            self.records[law].append(self.record(law, 400., False)); self.write_trace(law)
            self.final(law, False)
        self.receipt['status'] = 'completed'
        answer = self.answer()
        self.assertTrue(answer['conditional_frozen_scalar_allowed'])
        for case in answer['cases']:
            self.assertFalse(case['converged_target_evidence'])
            self.assertTrue(case['native_postqualification_report_available'])
            self.assertTrue(case['native_history_unchanged_verified'])

    def test_frozen_scalar_success_and_full_coupled_false_are_not_conflicting(self):
        law = 'SaturatedR2'; item = self.paths[law]
        self.receipt.update(mode='frozen-scalar', status='completed', alpha_levels=[400.])
        self.receipt['commands'] = [row for row in self.receipt['commands'] if row['name'] != 'native-Linear']
        self.receipt['cases'] = [item['row']]
        self.plan.update(alpha_levels=[400.])
        self.plan['cases'] = [row for row in self.plan['cases'] if row['inner_stress'] == law]
        item['xml'].write_text(self.xml(law, 'FrozenScalarR2'))
        item['row']['xml_sha256'] = self.plan['cases'][0]['xml_sha256'] = evidence.sha(item['xml'])
        item['row'].update(target_attempted=True, target_converged=False, target_physical_range=False)
        record = self.record(law, 400., True)
        record['frozen_velocity_pressure_diagnostic'] = True
        self.records[law] = [record]; self.write_trace(law)
        self.mutate_log(law, 'EP-AWF R2 stage test alpha 400 converged 0 physical range 0 history unchanged\n')
        path, raw = self.final(law, False)
        raw['stage_test'] = 'FrozenScalarR2'
        raw['stage_test_levels'][0]['converged'] = True
        self.write(path, raw); item['row']['report_sha256'] = evidence.sha(path)
        answer = self.answer(); case = answer['cases'][0]
        self.assertTrue(case['scalar_solve_convergence_evidence'])
        self.assertFalse(case['full_coupled_convergence_evidence'])
        self.assertTrue(case['converged_target_evidence'])
        self.assertTrue(case['full_coupled_root_unresolved'])
        self.assertEqual(case['evidence_status'], 'frozen_scalar_convergence_evidence_available')
        self.assertEqual(case['issues'], [])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_incomplete_tail_is_not_a_complete_record_and_blocks_conditional(self):
        path = self.paths['Linear']['trace']
        with path.open('a') as stream:
            stream.write(json.dumps(self.record('Linear', 400., True)))
        self.paths['Linear']['row']['trace_sha256'] = evidence.sha(path)
        answer = self.answer()
        self.assertFalse(answer['cases'][0]['complete_target_solve_record'])
        self.assertIsNotNone(answer['cases'][0]['trace_evidence']['unpublished_tail'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_malformed_complete_record_rejected(self):
        path = self.paths['Linear']['trace']; path.write_text('{broken}\n')
        self.paths['Linear']['row']['trace_sha256'] = evidence.sha(path)
        with self.assertRaises(ValueError): self.answer()

    def test_active_supervisor_rejected(self):
        self.receipt['status'] = 'running'
        with self.assertRaises(ValueError): self.answer()

    def test_changed_bound_log_rejected(self):
        self.paths['Linear']['log'].write_text('modified after receipt\n')
        with self.assertRaises(ValueError): self.answer()

    def test_changed_trace_rejected(self):
        with self.paths['Linear']['trace'].open('a') as stream: stream.write('{}\n')
        with self.assertRaises(ValueError): self.answer()

    def test_repeated_target_record_rejected(self):
        self.records['Linear'] += [self.record('Linear', 400., False)]*2; self.write_trace('Linear')
        with self.assertRaises(ValueError): self.answer()

    def test_changed_h0_rejected(self):
        self.records['Linear'][0]['history'] = [2., 2.]; self.write_trace('Linear')
        with self.assertRaises(ValueError): self.answer()

    def test_changed_source_rejected(self):
        (self.source/'kernel.cpp').write_text('// source mutated\n')
        with self.assertRaises(ValueError): self.answer()

    def test_changed_solver_image_rejected(self):
        (self.build/'IncNavierStokesSolver').write_bytes(b'other image')
        with self.assertRaises(ValueError): self.answer()

    def test_runtime_section_mismatch_rejected(self):
        after = copy.deepcopy(self.pin); after['python_packages']['fixture']['version'] = 'changed'
        self.write(self.run/'runtime-after.json', after)
        with self.assertRaises(ValueError): self.answer()

    def test_unbound_final_report_rejected(self):
        self.write(self.paths['Linear']['directory']/'case.evolution.json', {})
        with self.assertRaises(ValueError): self.answer()

    def test_changed_xml_beyond_stage_selector_rejected(self):
        xml = self.paths['Linear']['xml']; xml.write_text(self.xml('Linear', 'MatchedR2').replace('</NEKTAR>', '<PARAMETERS>changed</PARAMETERS></NEKTAR>'))
        digest = evidence.sha(xml)
        self.paths['Linear']['row']['xml_sha256'] = digest
        self.plan['cases'][0]['xml_sha256'] = digest
        with self.assertRaises(ValueError): self.answer()

    def test_duplicate_native_command_is_a_retry_and_rejected(self):
        self.receipt['commands'].append(copy.deepcopy(self.paths['Linear']['command']))
        with self.assertRaises(ValueError): self.answer()

    def test_unregistered_extra_native_command_rejected(self):
        extra = copy.deepcopy(self.paths['Linear']['command']); extra['name'] = 'native-unregistered-retry'
        self.receipt['commands'].append(extra)
        with self.assertRaises(ValueError): self.answer()

    def test_positive_receipt_evidence_cannot_be_ignored(self):
        self.paths['Linear']['row']['target_converged'] = True
        answer = self.answer()
        self.assertTrue(answer['cases'][0]['converged_target_evidence'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_matrix_length_rejected_even_when_not_in_complete_record(self):
        (self.paths['Linear']['directory']/'case.feasible.jsonl.solve2.J0.bin').write_bytes(b'bad')
        with self.assertRaises(ValueError): self.answer()

    def test_completed_history_difference_blocks_fallback(self):
        path, raw = self.final('Linear', False)
        raw['native_history'] = [[9., 9.]]; self.write(path, raw)
        self.paths['Linear']['row']['report_sha256'] = evidence.sha(path)
        answer = self.answer()
        self.assertFalse(answer['cases'][0]['native_history_unchanged_verified'])
        self.assertFalse(answer['conditional_frozen_scalar_allowed'])

    def test_exact_assessment_permits_only_conditional_scalar_scope(self):
        path = self.stored_assessment()
        gate = self.gate(path)
        self.assertTrue(gate['both_original_targets_observed'])
        self.assertFalse(gate['coupled_target_retry_permitted'])

    def test_manual_assessment_flag_does_not_override_actual_evidence(self):
        path = self.stored_assessment(); value = evidence.load(path)
        value['assessment']['cases'][0]['converged_target_evidence'] = False
        self.write(path, value)
        with self.assertRaises(ValueError): self.gate(path)

    def test_changed_matrix_after_assessment_rejected(self):
        path = self.stored_assessment()
        matrix = self.paths['Linear']['directory']/'case.feasible.jsonl.solve2.J0.bin'
        matrix.write_bytes(struct.pack('<4d', 2., 0., 0., 1.))
        with self.assertRaises(ValueError): self.gate(path)

    def test_changed_preflight_in_conditional_recipe_rejected(self):
        path = self.stored_assessment(); other = self.root/'other-preflight.json'
        other.write_bytes(self.preflight_path.read_bytes())
        with self.assertRaises(ValueError): self.gate(path, preflight=other)

    def test_changed_representation_in_conditional_recipe_rejected(self):
        path = self.stored_assessment()
        with self.assertRaises(ValueError): self.gate(path, representation='outer9')

    def test_changed_runtime_in_conditional_recipe_rejected(self):
        path = self.stored_assessment(); pin = copy.deepcopy(self.pin)
        pin['runtime_bindings'] = {'changed':{}}; self.write(self.runtime_pin, pin)
        with self.assertRaises(ValueError): self.gate(path)

    def test_complete_positive_target_blocks_direct_conditional_guard(self):
        self.records['SaturatedR2'].append(self.record('SaturatedR2', 400., True)); self.write_trace('SaturatedR2')
        self.save()
        with self.assertRaises(ValueError): self.gate()


if __name__ == '__main__':
    unittest.main()
