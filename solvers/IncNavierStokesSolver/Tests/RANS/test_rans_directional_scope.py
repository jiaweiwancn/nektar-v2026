#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy,json,os,subprocess,sys,tempfile,unittest
from pathlib import Path
from rans_directional_scope import classify,SCHEMA,SHARED,REQUIRED


class Scope(unittest.TestCase):
    def fixture(self):
        return dict(schema=SCHEMA,numerical_operator='GeneratedVolumeOwnerV1',response_max_time=.016,
            A93_response_out_of_scope=True,layout=dict(fields=[6,5],pressure=[4,3]),unknown_limit=3000,
            gates={k:dict(status='pass') for k in REQUIRED})
    def test_shared_prerequisites_and_A93_independence(self):
        r=self.fixture();r['gates']['A93_mesh_candidate']['status']='failed_or_unresolved'
        self.assertTrue(classify(r)['A95_future_prerequisites_satisfied'])
        for key in SHARED:
            c=copy.deepcopy(r);c['gates'][key]['status']='fail'
            self.assertFalse(classify(c)['A95_future_prerequisites_satisfied'],key)
    def test_reconstruction_is_not_accuracy(self):
        r=self.fixture();r['gates']['physical_budget_accuracy']['status']='failed_or_unresolved'
        self.assertFalse(classify(r)['A95_future_prerequisites_satisfied'])
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR')) as folder:
            p=Path(folder);(p/'in.json').write_text(json.dumps(r))
            process=subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('classify_rans_directional_scope.py')),
                '--input',str(p/'in.json'),'--output',str(p/'out.json')],capture_output=True,text=True)
            self.assertEqual(process.returncode,1)
            self.assertFalse(json.loads((p/'out.json').read_text())['response_extension_eligible'])
    def test_missing_and_mislabelled_prerequisite(self):
        r=self.fixture();del r['gates']['geometry'];self.assertFalse(classify(r)['A95_future_prerequisites_satisfied'])
        r=self.fixture();r['gates']['geometry']['primitive_status']='fail'
        with self.assertRaisesRegex(ValueError,'pass_label_disagrees'):classify(r)
    def test_scope_guards(self):
        for key,value in [('response_max_time',.1),('numerical_operator','LegacyCoordinateToleranceV1'),('unknown_limit',3001)]:
            r=self.fixture();r[key]=value
            with self.assertRaises(ValueError):classify(r)
        self.assertFalse(classify(self.fixture())['response_extension_eligible'])


if __name__=='__main__':unittest.main()
