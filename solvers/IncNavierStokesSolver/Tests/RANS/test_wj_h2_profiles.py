#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small Linux fixtures for measured-data semantics and boundary session intent."""
import json
import math
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from wj_h2_profiles import Completion,wall,QE,KAPPA,US,D0,NU,TOP,audit_inlet
from wj_h2_cases import case_xml

class Profiles(unittest.TestCase):
    def test_digitized_origins(self):
        c=Completion();self.assertEqual(len(c.rows),60)
        for panel in 'abcd':self.assertEqual(sum(r['panel']==panel for r in c.rows),15)
        r=next(r for r in c.rows if r['panel']=='b')
        self.assertAlmostEqual(float(r['y_mm']),float(r['raw_height'])*250)
        r=next(r for r in c.rows if r['panel']=='c')
        self.assertAlmostEqual(float(r['y_mm']),float(r['raw_height'])*250+3.5)
        self.assertEqual(c.H,.5)

    def test_completion_corner_and_family(self):
        for span in (0.,.5,1.):
            c=Completion(span);assert_allclose(c.epsilon(D0+c.H),(c.K/QE)**1.5/(KAPPA*c.H),rtol=2e-15)
            self.assertGreater(c.K,QE*US**2)
            self.assertGreater(c.k(TOP),0);self.assertGreater(c.epsilon(TOP),0)
        self.assertAlmostEqual((2.15**2+2*1.13**2)/2,3.58815)
        self.assertAlmostEqual((2*2.15**2+1.13**2)/2,5.26095)

    def test_flat_only_calibration(self):
        c=Completion();s=c.summary()
        self.assertTrue(0<c.h<c.H);self.assertEqual(len(s['fit']['residual']),3)
        self.assertLess(s['fit']['rms'],.02)
        self.assertGreater(s['cutoff']['reference_k_min'],0)
        self.assertGreater(s['cutoff']['Lambda'],1)
        self.assertGreater(s['cutoff']['molecular_fraction_H'],0)
        assert_allclose(s['cutoff']['Q_endpoints'],0,atol=2e-15)

    def test_negative_principal_is_retained(self):
        # This checks full tensor/action values, never a universal q threshold.
        c=Completion(0.);a=audit_inlet(c,count=41)
        self.assertEqual(a['status'],'OBSTRUCTED')
        self.assertLess(a['minima']['principal_in_plane']['value'],0)
        self.assertGreater(a['minima']['covariance']['value'],0)

    def test_production_session_roles(self):
        with tempfile.TemporaryDirectory() as d:
            for kind in ('physical','analytic'):
                p=case_xml(Path(d)/kind,kind,0);r=ET.parse(p).getroot();design=json.loads(p.with_name('case-design.json').read_text())
                info={e.get('PROPERTY'):e.get('VALUE') for e in r.findall('./CONDITIONS/SOLVERINFO/I')}
                self.assertEqual(info['SolverType'],'RANSEquilibriumWallEvolution')
                self.assertEqual(info['RANSWallBoundaryMode'],'OpenTunnelV1')
                self.assertEqual(info['RANSWallFunctionCoordinates'],'Physical')
                self.assertIsNone(r.find('./CONDITIONS/FUNCTION[@NAME="BoundaryValuesOuter"]/E[@VAR="p"]'))
                self.assertIn(.5,design['chart_z'])
                if kind=='analytic':
                    self.assertEqual(info['RANSWallOpenTop'],'ReferenceTractionV1')
                    self.assertIsNotNone(r.find('./CONDITIONS/FUNCTION[@NAME="TopDiffusiveFluxOuter"]/E[@VAR="epsilon"]'))
                    self.assertGreater(design['reference_k_min'],0)
                    # Distinct transverse outlet traction is required.
                    v=r.find('./CONDITIONS/FUNCTION[@NAME="OutletTractionOuter"]/E[@VAR="v"]')
                    self.assertAlmostEqual(float(v.get('VALUE')),US**2)
                else:self.assertIsNone(r.find('./CONDITIONS/FUNCTION[@NAME="TopTractionOuter"]'))

if __name__=='__main__':unittest.main()
