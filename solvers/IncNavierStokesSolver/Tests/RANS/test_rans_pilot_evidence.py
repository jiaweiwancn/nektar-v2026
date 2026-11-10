#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import tempfile,unittest
from pathlib import Path
from rans_pilot_evidence import epsilon_budget_gate,required_comparisons
from rans_preservation_oracle import Reference

class EvidenceTests(unittest.TestCase):
 def test_epsilon_budget_classifies_failure(self):
  ref=Reference('E94E95','Natural','Accepted');scale=2*(ref.Z-ref.H)*(ref.A/ref.H)*ref.us/ref.H
  row=dict(time=.016,alpha=2500.,free_constant_test_defect=0.)
  self.assertEqual(epsilon_budget_gate([row],ref)['status'],'pass')
  self.assertEqual(epsilon_budget_gate([dict(row,free_constant_test_defect=2e-10*scale)],ref)['status'],'fail')
  self.assertEqual(epsilon_budget_gate([],ref)['status'],'fail')
 def test_all_required_comparison_bodies_are_required(self):
  with tempfile.TemporaryDirectory() as d:
   rows=required_comparisons(Path(d))
   self.assertEqual(len(rows),12)
   self.assertTrue(all(x['status']=='fail' for x in rows.values()))
   self.assertIn('A93-N-EP-space-0.016',rows)
if __name__=='__main__':unittest.main()
