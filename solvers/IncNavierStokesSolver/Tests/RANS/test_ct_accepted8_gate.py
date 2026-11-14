# SPDX-License-Identifier: MIT
import copy,unittest
from ct_accepted8_gate import pair

class Checks(unittest.TestCase):
    def test_empty_pass_label_rejected(self):
        self.assertEqual(pair(dict(status='pass',qualification=dict(status='pass')))['status'],'fail')
    def test_missing_companion_reason(self):
        r=pair(dict(target=.008,members={},temporal={}))
        self.assertIn('missing_timestep_member',r['reasons'])
    def test_foreign_model(self):
        r=pair(dict(target=.008,members={'0.002':dict(time=.008,failure='',signature='A93',accepted_metrics=[])},temporal={}))
        self.assertIn('0.002:foreign_identity',r['reasons'])
    def test_nonfinite_time(self):
        r=pair(dict(target=.008,members={'0.002':dict(time=float('nan'),signature='',accepted_metrics=[])},temporal={}))
        self.assertIn('0.002:incomplete_or_rejected',r['reasons'])
    def test_passing_label_cannot_hide_temporal_defect(self):
        r=pair(dict(target=.008,members={},temporal={'U':dict(status='pass',measured_drift=1e-5,difference=dict(weighted_RMS=1e-3,scale=1.))}))
        self.assertIn('temporal:U',r['reasons'])
if __name__=='__main__':unittest.main()
