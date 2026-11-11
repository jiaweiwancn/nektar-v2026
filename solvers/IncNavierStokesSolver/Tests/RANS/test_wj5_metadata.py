#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Custom-normal physical rules and immutable requested-intent fixtures."""
import copy
from pathlib import Path
import tempfile
import unittest
import numpy as np
from numpy.testing import assert_allclose
from wj2_cases import local_normal_knots,validated_knots,normal_knots,case_xml
from wj4_measure import rule,boundary_rule,owned_rule
from wj_immutable import write_once,xml_identity

class Meshes(unittest.TestCase):
    def test_local_mesh_is_refinement_only(self):
        l1=normal_knots(1);m1=local_normal_knots('M1');m2=local_normal_knots('M2')
        self.assertEqual((len(l1),len(m1),len(m2)),(31,35,43))
        self.assertTrue(set(l1)<set(m1)<set(m2))
        self.assertEqual([v for v in m1 if v not in l1],[.085,.095,.125,.175])
        self.assertEqual([v for v in m2 if v not in m1],[.0825,.0875,.0925,.0975,.1125,.1375,.1625,.1875])
    def test_invalid_knots(self):
        for knots in ([0,.1,1,2,2],[0,1,2],[0,.1,float('nan'),1,2],[0,.1,1,1.5],[-1,.1,1,2]):
            with self.assertRaises(ValueError):validated_knots(knots,0.,2.,(.1,1.))
    def test_actual_composite_and_surface_measures(self):
        meshes=[local_normal_knots('M1'),local_normal_knots('M2')]
        for order,xorder in ((6,4),(8,6)):
            rows=rule(zmeshes=meshes,order=order,xorder=xorder,all_faces=True)
            q=np.array([r for r in rows if r[3]>0]);x,z,w=q[:,0],q[:,2],q[:,3]
            assert_allclose([sum(w),w@x**6,w@z**6],[2,2/7,128/7],rtol=3e-14,atol=3e-14)
            surface=boundary_rule(zmeshes=meshes,order=order,xorder=xorder)
            assert_allclose(sum(surface.values()),5.,rtol=0,atol=3e-14)
            self.assertTrue(all(r[3]==0 for r in rows if r[6]))
        physical=rule(zmeshes=meshes,all_faces=True)
        rows=owned_rule(physical,1,[0,.5,1],meshes[0])
        for z in meshes[0][1:-1]:
            corner=[r for r in rows if r[:3]==[.5,0.,z]]
            self.assertEqual(len({r[5] for r in corner}),4)
            self.assertEqual({(r[6],r[7]) for r in corner},{(-1,-1),(-1,1),(1,-1),(1,1)})
    def test_knots_and_comparison_role_are_immutable(self):
        with tempfile.TemporaryDirectory() as t:
            path=Path(t)/'intent.json';a=dict(level=1,zknots=local_normal_knots('M1'),role='horizontal',dt=3.2e-5,target=128)
            write_once(path,a);write_once(path,a)
            for change in (dict(zknots=local_normal_knots('M2')),dict(role='timestep'),dict(dt=1.6e-5),dict(target=32)):
                b=copy.deepcopy(a);b.update(change)
                with self.assertRaisesRegex(ValueError,'immutable intent'):write_once(path,b)
    def test_case_defaults_and_custom_identity(self):
        with tempfile.TemporaryDirectory() as t:
            base=case_xml(Path(t)/'a',level=1)
            same=case_xml(Path(t)/'b',level=1,zknots=normal_knots(1))
            custom=case_xml(Path(t)/'c',level=1,zknots=local_normal_knots('M1'))
            self.assertNotEqual(xml_identity(same)['semantic_sha256'],xml_identity(custom)['semantic_sha256'])
            # Checkpoint destinations are intentionally case-specific.
            import xml.etree.ElementTree as ET
            for xml in (base,same):
                root=ET.parse(xml).getroot();self.assertEqual(len(root.findall('./GEOMETRY/ELEMENT/Q')),30)

if __name__=='__main__':unittest.main()
