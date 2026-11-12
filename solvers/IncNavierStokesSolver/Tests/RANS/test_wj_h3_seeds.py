#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import unittest
import tempfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from wj_h2_profiles import Completion,NU,D0,TOP,DELTA,KAPPA,QE
from wj_h2_cases import expressions,continuous_outer,case_xml


class PressureSeeds(unittest.TestCase):
    def test_independent_normal_stress_balance_and_old_failure(self):
        c=Completion();inner,outer,*_=expressions(c,'physical');G,R=continuous_outer(c,D0+c.H)
        # Physical tensors derived independently from p/rho=pi-2k/3.
        Tw=-(float(inner['p'])-2*c.K/3)*np.eye(3)-2*c.K/3*np.eye(3)
        To=-float(outer['p'])*np.eye(3)+2*c.K/3*np.eye(3)-R
        assert_allclose(Tw[2,2],To[2,2],atol=2e-18,rtol=0)
        old_pi=2*c.K/3-R[2,2];old_Tw=-old_pi
        self.assertGreater(abs(old_Tw-To[2,2]),.00559)
        self.assertAlmostEqual(float(inner['p']),-.002795906118110805,places=16)
        self.assertEqual(G[2,2],0.)

    def test_vertical_balance_is_a_distinct_continuous_seed(self):
        c=Completion();C=continuous_outer(c,D0+c.H)[1][2,2]-2*c.K/3
        for z in np.linspace(D0+c.H,TOP,101):
            _,R=continuous_outer(c,z);k=c.k(z);pi=2*k/3-R[2,2]+C
            assert_allclose(-R[2,2]+2*k/3-pi,-C,atol=2e-18,rtol=0)
        v=expressions(c,'physical','VerticalBalanceV1','FrozenSeedTractionV1')
        self.assertNotEqual(v[1]['p'],'0');self.assertEqual(v[0]['p'],expressions(c,'physical')[0]['p'])
        self.assertLess(v[-1]['expression_errors']['outer_pi'],2e-11)
        self.assertLess(v[-1]['expression_errors']['outlet_Txx'],2e-11)
        self.assertLess(v[-1]['expression_errors']['outlet_Tzx'],2e-11)
        self.assertAlmostEqual(float(v[2]['u']),-C)

    def test_analytic_neutral_sign_is_not_reversed(self):
        c=Completion();v=expressions(c,'analytic')
        self.assertAlmostEqual(float(v[0]['p']),-QE*.051**2/(81/20))
        self.assertEqual(v[1]['p'],'0')

    def test_length_label_and_unchanged_formula(self):
        c=Completion();s=c.summary();self.assertEqual(s['blackadar_parameter_fraction_delta'],.09)
        self.assertAlmostEqual(s['asymptotic_length_fraction_delta'],.12055342908373864,places=15)
        for z in np.linspace(D0+c.H,TOP,51):
            L=.09*DELTA;r=z-D0;ell=KAPPA*r*(1+KAPPA*c.H/L)/(1+KAPPA*r/L)
            assert_allclose(c.epsilon(z),(572/6561)**.75*c.k(z)**1.5/ell,rtol=0,atol=0)

    def test_mesh_and_configuration_validation(self):
        with tempfile.TemporaryDirectory() as d:
            p=case_xml(Path(d)/'three','physical',0,xknots=[0.,8/3,16/3,8.],tangential_order=12,
                pressure_seed='VerticalBalanceV1',outlet='FrozenSeedTractionV1')
            r=ET.parse(p).getroot();self.assertIsNone(r.find('./CONDITIONS/FUNCTION[@NAME="BoundaryValuesOuter"]/E[@VAR="p"]'))
            self.assertIsNotNone(r.find('./CONDITIONS/FUNCTION[@NAME="OutletTractionOuter"]/E[@VAR="v"]'))
            for knots in ([0,8,4],[0,float('nan'),8],[-1,8]):
                with self.assertRaises(ValueError):case_xml(Path(d)/'invalid',xknots=knots)
            with self.assertRaises(ValueError):case_xml(Path(d)/'invalid',tangential_order=9)
            with self.assertRaises(ValueError):expressions(Completion(),'physical','SignCorrectedV1','FrozenSeedTractionV1')

if __name__=='__main__':unittest.main()
