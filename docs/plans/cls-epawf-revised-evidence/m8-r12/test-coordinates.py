#!/usr/bin/env python3
"""Focused independent checks of frozen residual coordinates and mixed rows."""
import importlib.util
from pathlib import Path
import os
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np

class Coordinates(unittest.TestCase):
    def test_feasible_seed_preserves_affine_origin_and_objective(self):
        J=np.array([[2.,1.,0.],[0.,3.,1.],[1.,0.,1.]])
        r=np.array([4.,-2.,.3]);s=.125
        rf=r.copy();rf[:2]/=s
        seed=np.array([-.1,0.,-.2])
        B=np.diag([2.,3.,0.]);T=np.linalg.solve(J,B)
        newton=np.linalg.solve(J,-rf)
        c=np.zeros(3);c[:2]=(rf+J@seed)[:2]/np.diag(B)[:2]
        self.assertLess(np.linalg.norm(seed-(newton+T@c)),1e-14)
        for y in (np.zeros(3),np.array([.1,-.2,0.])):
            d=seed+T@y
            self.assertAlmostEqual((rf+J@d)[2],0.)
            native=(rf+J@d)[:2]/np.diag(B)[:2]
            self.assertAlmostEqual(float(native@native),float((y+c)@(y+c)))
            physical=r+s*J@d
            self.assertAlmostEqual(physical[2],(1.-s)*r[2])
        # Replacing d_N by the seed, without the objective shift, is different.
        self.assertGreater(np.linalg.norm(seed-newton),1.)

    def test_numerical_budget_does_not_change_physical_xml(self):
        path=(Path(os.environ['RANS_SOURCE_ROOT'])/'docs/plans/cls-epawf-revised-evidence/m8-r3/prepare-predictor.py'
              if 'RANS_SOURCE_ROOT' in os.environ else Path(__file__).parents[1]/'m8-r3/prepare-predictor.py')
        spec=importlib.util.spec_from_file_location('r3',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        root=ET.fromstring('<NEKTAR><CONDITIONS><PARAMETERS><P>beta=.4</P></PARAMETERS><SOLVERINFO/></CONDITIONS></NEKTAR>')
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'case.xml';ET.ElementTree(root).write(p)
            original=module.protected_xml(p)
            ET.SubElement(root.find('.//SOLVERINFO'),'I',PROPERTY='RANSWallR12CutLimit',VALUE='256')
            ET.SubElement(root.find('.//SOLVERINFO'),'I',PROPERTY='RANSWallR12Experiment',VALUE='Complete')
            ET.ElementTree(root).write(p)
            self.assertEqual(original,module.protected_xml(p))
            root.find('.//PARAMETERS/P').text='beta=.41';ET.ElementTree(root).write(p)
            self.assertNotEqual(original,module.protected_xml(p))

if __name__=='__main__':unittest.main()
