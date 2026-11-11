#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Negative controls for the corrective-package accuracy and operator gates."""
import copy
import unittest
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from rans_ownership_contract import budget_gate, classify_scope, OPERATOR, SCHEMA
from classify_rans_physical_budget_accuracy import classify_row
from rans_duration_contract import MEMBERS


class OwnershipContract(unittest.TestCase):
    def fixture(self):
        cfg=dict(model='A95',policy='Natural',wall='Accepted')
        rows=[dict(kind=kind,equation=eq,quadrature_uncertainty=1e-12,
                   signed_values=dict(zip(MEMBERS,[1e-4,1e-5,1e-4,1e-5],strict=True)))
              for kind in ('baseline','perturbed') for eq in ('inner_k','outer_k','outer_epsilon','whole_patch_k')]
        a=dict(configuration=cfg,time=.016,execution_status='completed',reconstruction_status='pass',
               finite_defects=rows,quadrature=[dict(limit=1e-9,changes=dict(volume=1e-12,face=1e-12))])
        accuracy=dict(configuration=cfg,time=.016,source_sha256='analysis',status='pass',checks=list(map(classify_row,rows)))
        members={m:dict(perturbed='p-'+m,baseline='b-'+m) for m in MEMBERS}
        binding=dict(schema='NativeOwnershipBudgetBindingV1',configuration=cfg,numerical_operator=OPERATOR,
            analysis_sha256='analysis',accuracy_sha256='accuracy',cohort='A95-full',amplitude=2e-5,time=.016,
            members=members,controls=dict(volume_orders=[48,64],physical_face_orders=[16,32]))
        expected=dict(configuration=cfg,analysis_sha256='analysis',accuracy_sha256='accuracy',cohort='A95-full',time=.016,members=copy.deepcopy(members))
        return a,accuracy,binding,expected
    def test_passing(self):self.assertEqual(budget_gate(*self.fixture())['status'],'pass')
    def test_reconstruction_not_accuracy(self):
        a,c,b,e=self.fixture();a['finite_defects'][0]['signed_values']['B11']=2e-4
        c.update(status='failed_or_unresolved',checks=list(map(classify_row,a['finite_defects'])))
        self.assertEqual(budget_gate(a,c,b,e)['accuracy_status'],'failed_or_unresolved')
        self.assertEqual(a['reconstruction_status'],'pass')
    def test_missing_receipt(self):
        a,_,b,e=self.fixture()
        with self.assertRaises(KeyError):budget_gate(a,{},b,e)
    def test_foreign_receipt(self):
        a,c,b,e=self.fixture();c['configuration']=dict(c['configuration'],wall='EP')
        with self.assertRaisesRegex(ValueError,'configuration_mismatch'):budget_gate(a,c,b,e)
    def test_stale_receipt(self):
        a,c,b,e=self.fixture();b['accuracy_sha256']='stale'
        with self.assertRaisesRegex(ValueError,'stale_budget_accuracy'):budget_gate(a,c,b,e)
    def test_horizon(self):
        a,c,b,e=self.fixture();a['time']=.1
        with self.assertRaisesRegex(ValueError,'horizon_mismatch'):budget_gate(a,c,b,e)
    def test_companion(self):
        a,c,b,e=self.fixture();del b['members']['B11-dt2']
        with self.assertRaisesRegex(ValueError,'companion_missing'):budget_gate(a,c,b,e)
    def test_pass_label_recomputed(self):
        a,c,b,e=self.fixture();a['finite_defects'][0]['signed_values']['B11']=2e-4
        with self.assertRaisesRegex(ValueError,'primitive_recomputation'):budget_gate(a,c,b,e)
        c['checks']=list(map(classify_row,a['finite_defects']))
        with self.assertRaisesRegex(ValueError,'pass_label_disagrees'):budget_gate(a,c,b,e)
    def test_old_operator(self):
        a,c,b,e=self.fixture();b['numerical_operator']='LegacyCoordinateToleranceV1'
        with self.assertRaisesRegex(ValueError,'old_operator_budget'):budget_gate(a,c,b,e)
    def test_scope_requires_accuracy(self):
        keys=['ownership_geometry','history_cache_CLS','natural_short_matrix','selected_prescribed','A93_neutral',
              'A95_neutral','A95_response_fields','A95_response_native_weak','physical_budget_reconstruction','physical_budget_accuracy']
        r=dict(schema=SCHEMA,numerical_operator=OPERATOR,gates={k:dict(status='pass') for k in keys},
               out_of_scope=dict(A93_response=True),response_max_time=.016)
        self.assertEqual(classify_scope(r)['status'],'pass')
        r['gates']['physical_budget_accuracy']['status']='failed_or_unresolved'
        q=classify_scope(r);self.assertNotEqual(q['status'],'pass');self.assertFalse(q['response_extension_eligible'])
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as d:
            inp=Path(d)/'input.json';out=Path(d)/'output.json';inp.write_text(json.dumps(r))
            result=subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('classify_rans_ownership_scope.py')),
                                   '--input',str(inp),'--output',str(out)],capture_output=True,text=True)
            self.assertEqual(result.returncode,1);self.assertFalse(json.loads(out.read_text())['response_extension_eligible'])


if __name__=='__main__':unittest.main()
