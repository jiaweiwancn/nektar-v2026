"""Focused stdlib predictor-contract metadata tests; no solver/operator calls."""
import copy
import importlib.util
import json
from pathlib import Path
import struct
import sys
import tempfile
from types import SimpleNamespace
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location('r3_predictor_metadata', HERE/'prepare-predictor.py')
predictor = importlib.util.module_from_spec(spec); spec.loader.exec_module(predictor)


class PredictorMetadataTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='r3-predictor-metadata-')
        self.addCleanup(self.temporary.cleanup); self.root = Path(self.temporary.name)
        self.state, self.history = [1.]*2165, [2.]*2165
        self.old_h0 = predictor.H0; predictor.H0 = predictor.vector(self.history)
        self.addCleanup(setattr, predictor, 'H0', self.old_h0)
        self.packet = dict(schema=1, inner_stress_model='SaturatedR2', time=.0025, alpha=3200.,
            state=self.state, history=self.history, residual=[.1]*2165, scaled_block_norms=[5e-12]*8, constraints=[],
            source_record=dict(record_index=0, converged=False, terminal='feasible_residual_stagnation'))
        for name in ('state', 'history', 'residual'):
            path = self.root/(name+'.bin'); values = self.packet[name]
            path.write_bytes(struct.pack('<'+str(len(values))+'d', *values))
            self.packet[name+'_binary'] = predictor.identity(path)
        trace = self.root/'historical-trace.jsonl'; trace.write_text('immutable extracted-record fixture\n')
        self.packet['source_record']['trace'] = predictor.identity(trace)
        self.packet_path = self.root/'packet.json'; self.write(self.packet_path, self.packet)
        self.ledger = self.root/'input-ledger.json'; self.update_ledger()
        self.executable = self.root/'solver'; self.executable.write_bytes(b'never executed metadata fixture')
        self.runtime = self.root/'runtime.json'
        self.write(self.runtime, dict(success=True, clean=True, source_commit='R3-source',
            runtime_resolved_paths={'build/IncNavierStokesSolver':str(self.executable)},
            runtime_artifacts={'build/IncNavierStokesSolver':predictor.sha(self.executable)}))
        self.reference, self.current_xml, self.target_xml = (self.root/name for name in ('reference.xml', 'current.xml', 'target.xml'))
        self.reference.write_text(self.xml('ScalarFeasibleR2', 'MatchedR2'))
        self.current_xml.write_text(self.xml('ScalarFeasibleR3', 'Off'))
        self.target_xml.write_text(self.xml('ScalarFeasibleR3', 'ExactTargetR3'))
        self.current = self.root/'initial.json'
        self.write(self.current, dict(initial=self.history, time=0., steps=0, eta=0., dt=.0025,
            inner_stress_model='SaturatedR2', nonlinear_policy='ScalarFeasibleR3', signature='same native signature'))
        self.args = SimpleNamespace(packet=self.packet_path, input_ledger=self.ledger, reference_session=self.reference,
            current_preflight=self.current, current_preflight_session=self.current_xml, target_session=self.target_xml,
            runtime_pin=self.runtime, expected_source='R3-source')

    @staticmethod
    def write(path, value):
        Path(path).write_text(json.dumps(value)+'\n')

    @staticmethod
    def xml(policy, stage):
        action = '<I PROPERTY="RANSWallActionPolicy" VALUE="VerifiedR3Action"/>' if policy.endswith('R3') else ''
        return ('<NEKTAR><GEOMETRY DIM="2"/><CONDITIONS><SOLVERINFO>'
            '<I PROPERTY="RANSWallInnerStress" VALUE="SaturatedR2"/>'
            '<I PROPERTY="RANSWallNonlinearPolicy" VALUE="'+policy+'"/>'
            '<I PROPERTY="RANSWallStageTest" VALUE="'+stage+'"/>'+action+
            '</SOLVERINFO><PARAMETERS><P>H=2</P><P>RANSWallMaxIterations=12</P></PARAMETERS></CONDITIONS></NEKTAR>')

    def update_ledger(self):
        self.write(self.ledger, dict(states=[dict(packet=predictor.identity(self.packet_path))]))

    def test_valid_numerical_predictor_preserves_failed_origin_and_h0(self):
        result = predictor.prepare(self.args)
        self.assertFalse(result['origin']['converged'])
        self.assertEqual(result['target']['alpha'], 400.)
        self.assertFalse(result['target']['history_advanced'])
        self.assertEqual(result['expected_h0_sha256'], predictor.H0)
        self.assertTrue(result['coefficients_are_not_lifted_or_clipped'])

    def test_converged_auxiliary_reclassification_rejected(self):
        self.packet['source_record']['converged'] = True; self.write(self.packet_path, self.packet); self.update_ledger()
        with self.assertRaisesRegex(ValueError, 'saved failed'): predictor.prepare(self.args)

    def test_historical_binary_change_rejected(self):
        path = Path(self.packet['state_binary']['path']); path.write_bytes(b'!'+path.read_bytes()[1:])
        with self.assertRaisesRegex(ValueError, 'binding changed'): predictor.prepare(self.args)

    def test_changed_protected_physics_rejected(self):
        self.target_xml.write_text(self.target_xml.read_text().replace('H=2', 'H=3'))
        with self.assertRaisesRegex(ValueError, 'protected physics'): predictor.prepare(self.args)

    def test_missing_packet_extraction_ledger_binding_rejected(self):
        self.write(self.ledger, {})
        with self.assertRaisesRegex(ValueError, 'input ledger'): predictor.prepare(self.args)

    def test_changed_current_executable_rejected(self):
        self.executable.write_bytes(b'changed executable metadata fixture')
        with self.assertRaisesRegex(ValueError, 'executable changed'): predictor.prepare(self.args)

    def test_changed_actual_h0_rejected(self):
        raw = json.loads(self.current.read_text()); raw['initial'][0] = 3.; self.write(self.current, raw)
        with self.assertRaisesRegex(ValueError, 'original h0'): predictor.prepare(self.args)

    def test_current_signature_newton_budget_mismatch_rejected(self):
        self.target_xml.write_text(self.target_xml.read_text().replace('RANSWallMaxIterations=12', 'RANSWallMaxIterations=16'))
        with self.assertRaisesRegex(ValueError, 'different Newton budget'): predictor.prepare(self.args)


if __name__ == '__main__':
    unittest.main()
