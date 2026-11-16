# SPDX-License-Identifier: MIT
"""Physical comparison conventions and fail-fast WP2 gate contracts."""
import copy, unittest
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


if __name__=='__main__':unittest.main()
