# SPDX-License-Identifier: MIT
import copy,json,tempfile,unittest
from pathlib import Path
from wj_h15_observer import identity
from wj_h22_audit_roles import native_full_map_parity,select_auxiliary_audits

def fixture():
    old=dict(state=[1.,-0.],fixed_reference=[2.],fixed_defect=[0.],dt=.0025,time=.005,signature='own-dt',
        original_qh0_state_sha256='original',native_full_audit_status='NOT_SCHEDULED',native_hard_audit={},last_stage_alpha=400.,last_stage_history=[3.])
    audit=copy.deepcopy(old);audit.update(native_full_audit_status='PASS',native_hard_audit=dict(admissible=[1.],margin=[.004]),
                                       last_stage_alpha=0.,last_stage_history=[99.])
    current=copy.deepcopy(old);current.update(native_full_audit_status='PASS',native_hard_audit=copy.deepcopy(audit['native_hard_audit']))
    return current,old,audit

class AuditRoleParity(unittest.TestCase):
    def test_auxiliary_selection_requires_the_sealed_index_and_bytes(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);audit=root/'B-startup.json';audit.write_text('{}')
            commit=root/'B-startup-commit.json';commit.write_text('{}')
            selected=root/'controls/inputs--selected.json';selected.parent.mkdir()
            selected.write_text(json.dumps(dict(states={'M':{'2':dict(audit_report=identity(audit),qualified_commit=identity(commit))}})))
            index=root/'manifest.json';index.write_text(json.dumps(dict(files=[identity(selected)])))
            inputs=dict(prior_indices={'20':dict(publication=identity(index))})
            result=select_auxiliary_audits(inputs)
            self.assertEqual(result['M']['2']['audit_report'],identity(audit))
            audit.write_text('{"changed":true}')
            with self.assertRaisesRegex(AssertionError,'changed pinned input'):select_auxiliary_audits(inputs)

    def test_separate_audit_does_not_replace_producing_callback(self):
        current,old,audit=fixture();saved=copy.deepcopy(old)
        p=native_full_map_parity(current,old,audit)
        self.assertEqual(p['status'],'PASS');self.assertFalse(p['history_or_alpha_replaced'])
        self.assertEqual(old,saved);self.assertEqual(old['last_stage_alpha'],400.)
        self.assertEqual(old['last_stage_history'],[3.])

    def test_empty_producer_map_never_becomes_a_full_map_or_false_PASS(self):
        current,old,audit=fixture()
        with self.assertRaisesRegex(ValueError,'separately pinned'):native_full_map_parity(current,old)
        missing=copy.deepcopy(current);missing['native_full_audit_status']='NOT_SCHEDULED'
        with self.assertRaises(AssertionError):native_full_map_parity(missing,old,audit)

    def test_wrong_state_signed_zero_signature_or_map_is_rejected(self):
        current,old,audit=fixture()
        for key,value in (('state',[1.,0.]),('signature','other-dt'),('time',.01),('fixed_defect',[1.])):
            bad=copy.deepcopy(audit);bad[key]=value
            with self.assertRaises(AssertionError):native_full_map_parity(current,old,bad)
        bad=copy.deepcopy(audit);bad['native_hard_audit']['margin']=[.005]
        with self.assertRaises(AssertionError):native_full_map_parity(current,old,bad)

    def test_historical_completed_full_audit_path_stays_exact(self):
        current,old,audit=fixture();old.update(native_full_audit_status='PASS',native_hard_audit=copy.deepcopy(current['native_hard_audit']))
        p=native_full_map_parity(current,old)
        self.assertEqual(p['map_source_role'],'genuine producing report full audit')
        self.assertEqual(old['last_stage_alpha'],400.)

if __name__=='__main__':unittest.main()
