# SPDX-License-Identifier: MIT
"""Physical comparison conventions and fail-fast WP2 gate contracts."""
import copy, unittest,tempfile,json,sys
from pathlib import Path
from unittest.mock import patch
from wp2_qualification import common_differences,quarter_budget_differences,physical_gates


class PhysicalComparison(unittest.TestCase):
    def samples(self):
        r=[1.,0.,.02,.01,0.,-.013,.002,.005,.003,.001,0.,2.,0.,0.]
        return dict(columns='U,W,k,epsilon,pi,mechanical_pressure,nut,total_shear,total_normal_stress,production,U_x,U_z,W_x,W_z',
                    values=[r[:],r[:]],reference_values=[r[:],r[:]],weights=[.2,.8],regions=[0,1])
    def test_absolute_regional_pressure_and_complete_outer_epsilon(self):
        a=self.samples();b=copy.deepcopy(a);b['values'][0][3]+=100
        for row in b['values']:row[4]+=.003
        d=common_differences(a,b,.1,.02,.0065)
        self.assertEqual(d['epsilon']['dimensional_L2'],0.)
        self.assertAlmostEqual(d['pi']['dimensional_L2'],.003)
        self.assertFalse(quarter_budget_differences(a,b,.1,.02,.0065)['passed'])
    def test_nonzero_small_difference_and_full_gradient(self):
        a=self.samples();b=copy.deepcopy(a);b['values'][1][0]+=1e-5
        self.assertTrue(quarter_budget_differences(a,b,.1,.02,.0065)['passed'])
        b['values'][0][12]=.1
        self.assertFalse(quarter_budget_differences(a,b,.1,.02,.0065)['passed'])
    def test_changed_grid_rejected(self):
        a=self.samples();b=copy.deepcopy(a);b['weights'][0]*=2
        with self.assertRaises(Exception):common_differences(a,b,.1,.02,.0065)


class ContinuationOwnership(unittest.TestCase):
    def test_changed_runtime_and_case_do_not_resume(self):
        from wp2_native import Run
        from run_cls_epawf_picard import sha
        with tempfile.TemporaryDirectory() as temp:
            r=Path(temp);s=r/'solver';s.write_text('immutable fixture');pre=r/'pre.json';pre.write_text('{}');out=r/'run';out.mkdir()
            d=dict(status='completed_bounded_trajectory',model='standard',solver_sha256=sha(s),case_definition=dict(level=4,order=16,projection='RadialEnergyV1'))
            (out/'progress.json').write_text(json.dumps(d))
            Run(r,s,r,pre,out,'standard',4,16,'RadialEnergyV1',True)
            with self.assertRaises(Exception):Run(r,s,r,pre,out,'standard',5,16,'RadialEnergyV1',True)
            s.write_text('different runtime')
            with self.assertRaises(Exception):Run(r,s,r,pre,out,'standard',4,16,'RadialEnergyV1',True)
    def test_unexpected_failure_stops_all_later_blocks(self):
        import run_wp2_preservation as driver
        made=[]
        class Failed:
            def __init__(self,*args):self.result={};self.calls=0;made.append(self)
            def save(self):pass
            def native(self,*args,**kwargs):self.calls+=1;raise RuntimeError('unexpected native rejection')
        args=['wp2','--root','/unused','--solver','/unused','--h1','/unused','--preflight','/unused','--output','/unused','--model','standard','--mode','preserve','--dt','.02']
        with patch.object(driver,'Run',Failed),patch.object(sys,'argv',args):
            with self.assertRaisesRegex(RuntimeError,'unexpected native rejection'):driver.main()
        self.assertEqual(made[0].calls,1)
        self.assertEqual(made[0].result['status'],'terminal_stop')


if __name__=='__main__':unittest.main()
