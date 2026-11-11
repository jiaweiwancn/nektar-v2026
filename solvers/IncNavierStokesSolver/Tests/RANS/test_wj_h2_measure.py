#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
import unittest
import numpy as np
from study_wj_h2 import observation_rule,comparison
from wj_h2_profiles import D0,TOP

class Measurements(unittest.TestCase):
    def packet(self):
        row=dict(physical=[1.,0.,1.],inner=False,U=[1.,0,0],k=.01,epsilon=.001,pi=0.)
        return dict(open_boundary=dict(time=.04,stations=[row],interface=[dict(K=.01,cutoff=.02)],
            boundary_rows=[dict(role='ground',actual_traction=[.0026,0,0])]))
    def test_common_physical_polynomials(self):
        points,w=observation_rule([dict(chart_z=[0.,.25,.5,2.,20.],reference_cutoff=.02)],8,4)
        p=np.array(points);w=np.array(w)
        self.assertAlmostEqual(float(w.sum()),8*(TOP-D0),places=11)
        self.assertAlmostEqual(float(w@p[:,0]**2),(8**3/3)*(TOP-D0),places=9)
        self.assertAlmostEqual(float(w@p[:,2]**2),8*(TOP**3-D0**3)/3,places=8)
    def test_nonfinite_pressure_cannot_pass(self):
        a=self.packet();b=copy.deepcopy(a);b['open_boundary']['stations'][0]['pi']=float('nan')
        with self.assertRaisesRegex(AssertionError,'nonfinite'):comparison(a,b,[1.],.01,.001)
    def test_wrong_time_refused(self):
        a=self.packet();b=copy.deepcopy(a);b['open_boundary']['time']=.08
        with self.assertRaisesRegex(AssertionError,'times differ'):comparison(a,b,[1.],.01,.001)
    def test_positive_field_difference_not_hidden(self):
        a=self.packet();b=copy.deepcopy(a);b['open_boundary']['stations'][0]['k']+=.001
        self.assertEqual(comparison(a,b,[1.],.01,.001)['status'],'UNRESOLVED')
    def test_zero_difference_uses_fixed_pressure_level(self):
        a=self.packet();self.assertEqual(comparison(a,a,[1.],.01,.001)['status'],'PASS')
        b=copy.deepcopy(a);b['open_boundary']['stations'][0]['pi']=.001
        self.assertGreater(comparison(a,b,[1.],.01,.001)['fields']['pi']['normalized_L2'],0)

if __name__=='__main__':unittest.main()
