# SPDX-License-Identifier: MIT
"""Regression for nested axial cuts and owner-preserving comparison."""

import unittest

from c1_paired_grid import subdivisions
from c1_paired_profiles import FIELDS, SCALED, approached, compare


def row(owner, sides, Gxz, traction):
    a, r = sides
    return dict({k:0. for k in FIELDS},owner=owner, x=.5, z=.4,
                a_left=-.5 if a<0 else 0., a_right=0. if a<0 else .5,
                r_low=0. if r<0 else .1, r_high=.1 if r<0 else .15,
                axial_side=a,radial_side=r,Gxz=Gxz,traction_axial_x=traction,
                angular_difference=1e-8,mu=.01)


class PairedProfiles(unittest.TestCase):
    def test_nested_subdivision_preserves_other_intervals(self):
        edges=[-3.,-1.,1.,3.,8.]
        inserted,mids=subdivisions(edges,(0.,3.))
        self.assertEqual(inserted,[0.])
        self.assertEqual(mids,[.5,2.])
        _,more=subdivisions(sorted(edges+inserted+mids),(0.,3.))
        self.assertEqual(more,[.25,.75,1.5,2.5])

    def test_extrema_and_jump_keep_owner_and_direction(self):
        point=(0.,.1)
        a=approached({point:[row(2,(-1,0),.25,.3),row(3,(1,0),.5,.4)]},.15,(-.5,.5,.2))
        b=approached({point:[row(9,(0,0),.75,.5)]},.15,(-.5,.5,.2))
        comp=compare(a,b,{v:1. for v in SCALED.values()},{"margin":1e-9})
        self.assertEqual(comp["approach_samples"],4)
        for label in ("Gxz","gradient_frobenius"):
            result=comp["maxima"][label]
            self.assertEqual(result["approach"],[-1,-1])
            self.assertEqual(result["a"]["owner"],2)
            self.assertEqual(result["b"]["owner"],9)
            self.assertEqual(result["difference"],.5)
            self.assertEqual(result["a"]["cartesian"],[.5,.4])
        self.assertAlmostEqual(comp["maxima"]["jump_traction_axial_x"]["difference"],.1)
        self.assertAlmostEqual(comp["maxima"]["upper_margin"]["upper"],2.1e-8)


if __name__=="__main__":
    unittest.main()
