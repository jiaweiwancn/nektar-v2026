# SPDX-License-Identifier: MIT
import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from study_wj_h16 import CASES,checked_exit,same_cp
from study_wj_h15 import qualified_step
from wj_h16_observer import active_source_bindings,packet_fields
from wj_h15_observer import atomic,identity,vector_hash
from run_wj_h16 import reserve_strict,halt,CAPS

class RepairedMonitor(unittest.TestCase):
    def test_child_nonzero_exit_cannot_be_reported_as_success(self):
        p=subprocess.run([sys.executable,'-c','raise SystemExit(7)'])
        with self.assertRaisesRegex(RuntimeError,'child exit 7'):checked_exit(p.returncode)
        checked_exit(1,1)  # an explicitly expected fault is distinct

    def test_six_actual_attempts_and_no_replacement_case(self):
        self.assertEqual(set(CASES),set('ABCDE'))
        self.assertEqual(sum(v[0]for v in CASES.values()),7)  # C pending opportunity
        self.assertEqual(sum(1 if k=='C'else v[0]for k,v in CASES.items()),6)
        self.assertEqual(CAPS['target'],1800);self.assertEqual(CAPS['target_readback'],120)

    def test_durable_cursor_does_not_follow_prepared_or_unpublished_step_commit(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d);m=p/'monitor';m.mkdir();payload=m/'payload';payload.write_bytes(b'bound')
            def commit(n):return dict(schema='WJQualificationCommitV1',binding={'absolute_step':str(n)},**{k:identity(payload)for k in ('checkpoint','callback','report','consumed_ack','prepared_cursor')})
            atomic(m/'step-1.qualified.json',dict(status='PREPARED_PASS'))
            self.assertEqual(qualified_step(p),0)
            atomic(m/'step-1.qualification-commit.json',commit(1));self.assertEqual(qualified_step(p),0)
            atomic(m/'qualified-cursor.json',commit(1));self.assertEqual(qualified_step(p),1)
            atomic(m/'step-2.qualification-commit.json',commit(2));self.assertEqual(qualified_step(p),1)

    def test_strict_slot_once_and_first_halt_preserved(self):
        with tempfile.TemporaryDirectory()as d:
            p=Path(d);(p/'accounting.json').write_text(json.dumps(dict(closed=False)))
            reserve_strict(p,'target')
            with self.assertRaises(AssertionError):reserve_strict(p,'target')
            halt(p,'target_halt',dict(reason='first'));halt(p,'target_halt',dict(reason='later'))
            self.assertEqual(json.loads((p/'accounting.json').read_text())['target_halt']['reason'],'first')

    def test_active_cache_modules_are_verified_not_historical_aliases(self):
        import wj_h15_observer
        rec=identity(wj_h15_observer.__file__);r=active_source_bindings({'source_files':[rec]})
        self.assertEqual(r[0]['active'],rec)
        wrong=dict(rec,sha256='0'*64)
        with self.assertRaises(AssertionError):active_source_bindings({'source_files':[wrong]})

    def test_completed_packet_requires_genuine_stage_and_explicit_context(self):
        r=dict(schema='WJCompletedStepV1',time=.02,dt=.01,last_stage_alpha=100.,stage_audits=[[.02,.01,100.]],step_audits=[[.02]],explicit_audits=[[0.,1.],[.01,2.]],first_retained_explicit_call=2,actual_explicit_calls=3,native_hard_audit={'admissible':[1.]},accepted_states=[[1.]],wj_audit_packets=[{}])
        self.assertEqual(packet_fields(r)['explicit_rows'],2)
        for key,value in (('schema','WJMonitoringStartupV1'),('last_stage_alpha',0.),('first_retained_explicit_call',1)):
            bad=copy.deepcopy(r);bad[key]=value
            with self.assertRaises(AssertionError):packet_fields(bad)

    def test_native_history_and_signed_zero_are_not_reconstructed(self):
        a=dict(state=[0.,1.],fixed_reference=[-0.,2.],fixed_defect=[0.],native_history=[[0.,1.]],native_history_times=[0.],time=.01,steps=1,dt=.01,signature='frozen')
        same_cp(a,copy.deepcopy(a))
        for k in ('state','fixed_reference'):
            b=copy.deepcopy(a);b[k][0]=-b[k][0]
            with self.assertRaises(AssertionError):same_cp(a,b)
        b=copy.deepcopy(a);b['native_history_times']=[.005]
        with self.assertRaises(AssertionError):same_cp(a,b)

if __name__=='__main__':unittest.main()
