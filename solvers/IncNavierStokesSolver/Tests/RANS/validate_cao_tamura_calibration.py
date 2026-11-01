#!/usr/bin/env python3
"""Algebra/data-isolation checks for CV3; these are not CFD calibration results."""
import math
import json
import os
from pathlib import Path
import sys
import tempfile
import xml.etree.ElementTree as ET
import unittest

sys.dont_write_bytecode=True
from cao_tamura_experiment import config, reference, inlet_values, frozen_data, mean_derivative, mean_ratio, generate, environment, wall_epsilon_expression
from cao_tamura_calibration import design, freeze, probes, thickness, select


class CalibrationGuards(unittest.TestCase):
    def test_shared_inlet_and_policy_traces(self):
        # Different closure Q/kappa must not silently retune physical inlet data.
        for surface in ('S','R'):
            for interpolation,epsilon_profile in ((b,e) for b in ('Geometric','Arithmetic') for e in ('NeutralLog','LocalShear')):
                inputs=[]
                for model in ('standard_original','cls_original'):
                    c=config(dict(surface=surface,model=model,ks=0 if surface=='S' else .004,
                        inlet_delta=.21,outer_blend=interpolation,epsilon_profile=epsilon_profile,velocity_scaling='incident_ReH',
                        inlet_reference_q=10/3,inlet_reference_kappa=math.sqrt(.48*1.3/(10/3))))
                    r=reference(c);inputs.append((c,r))
                    self.assertAlmostEqual(inlet_values(c['hm'],c,r)['epsilon'],r['wall_epsilon_near'],places=14)
                for i in range(201):
                    z=.022+(.8-.022)*i/200
                    a,b=[inlet_values(z,c,r) for c,r in inputs]
                    for key in a:self.assertAlmostEqual(a[key],b[key],places=13)
                    self.assertGreater(min(a['k'],a['epsilon']),0.)
                self.assertNotEqual(inputs[0][1]['q'],inputs[1][1]['q'])
                self.assertNotEqual(inputs[0][1]['kappa'],inputs[1][1]['kappa'])
        for model in ('standard','standard_fixed','cls_fixed'):
            c=config(dict(model=model));r=reference(c)
            self.assertAlmostEqual(inlet_values(c['hm'],c,r)['epsilon'],r['wall_epsilon_near'],places=14)

    def test_startup_boundary_matches_wall_policy(self):
        with tempfile.TemporaryDirectory() as root:
            for model in ('standard','standard_original','cls_original','standard_fixed','cls_fixed'):
                folder=Path(root)/model
                m=generate(dict(model=model,k_scale=1.5,inlet_reference_q=10/3,
                    inlet_reference_kappa=math.sqrt(.48*1.3/(10/3))),folder)
                tree=ET.parse(folder/'case.xml')
                entry=next(n for n in tree.getroot().find("CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='0']")
                           if n.attrib['VAR'] in ('eps','eps_tilde'))
                self.assertEqual(entry.attrib['VALUE'],wall_epsilon_expression(m['config']))
                c=m['config'];r=m['reference']
                self.assertAlmostEqual(inlet_values(c['hm'],c,r)['epsilon'],r['wall_epsilon_near'],places=13)

    def test_measured_rms_shear_family(self):
        for surface in ('S','R'):
            common=dict(surface=surface,k_profile='MeasuredRms',epsilon_profile='ShearLength',
                inlet_wall_collar=True,inlet_delta=.21,velocity_scaling='incident_ReH',
                inlet_reference_q=10/3,inlet_reference_kappa=math.sqrt(.48*1.3/(10/3)))
            cs=[config({**common,'model':m}) for m in ('standard_original','cls_original')]
            rs=[reference(c) for c in cs]
            for i in range(401):
                z=cs[0]['hm']*(.8/cs[0]['hm'])**(i/400)
                vs=[inlet_values(z,c,r) for c,r in zip(cs,rs)]
                for v in vs:self.assertGreater(min(v['k'],v['epsilon']),0.)
                if z>=cs[0]['point_n']:
                    for key in vs[0]:self.assertAlmostEqual(vs[0][key],vs[1][key],places=12)
                if z>cs[0]['hm']*1.001:
                    h=z*1e-5
                    fd=(mean_ratio(z+h,cs[0],rs[0]['coefficients'])-mean_ratio(z-h,cs[0],rs[0]['coefficients']))/(2*h)
                    self.assertAlmostEqual(fd,mean_derivative(z,cs[0],rs[0]['coefficients']),delta=3e-6)
            for c,r in zip(cs,rs):
                h=c['hm'];k=r['k_near'];u=math.sqrt(k/r['inlet_q'])
                self.assertAlmostEqual(r['Uinf']*mean_derivative(h,c,r['coefficients']),u/(r['inlet_kappa']*h),places=10)
                self.assertAlmostEqual(inlet_values(h,c,r)['epsilon'],r['wall_epsilon_near'],places=13)

    def test_actual_mean_gradient(self):
        for cap in (False,True):
            for surface in ('S','R'):
                c=config(dict(surface=surface,mean_outer_cap=cap));r=reference(c)
                for z in (.012,.04,.19,.225,.27,.32,.7):
                    h=1e-6*z
                    fd=(mean_ratio(z+h,c,r['coefficients'])-mean_ratio(z-h,c,r['coefficients']))/(2*h)
                    self.assertAlmostEqual(fd,mean_derivative(z,c,r['coefficients']),delta=2e-6)

    def test_design_and_holdout_barrier(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root)
            for surface,count in (('S',7),('R',9)):
                self.assertEqual(len(design(surface)),count)
                rows=probes(surface,root/(surface+'.csv'))
                self.assertFalse(any(r['id'].startswith('f3d') for r in rows))
                self.assertFalse(any(r['id'].startswith(('f4','f5','f6','f7','f8','f9')) for r in rows))
                for case in design(surface):
                    if surface=='S':self.assertEqual(case['parameters']['ks'],0.)
            with self.assertRaises(ValueError):freeze(root/'empty',root/'absent.json',root/'cannot-freeze','test')
            self.assertFalse((root/'cannot-freeze').exists())
            # Scenarios shift the observation coordinate by 1mm only.
            rough=probes('R',root/'rough.csv');indexed={r['id']:r for r in rough}
            for r in rough:
                if r['id'].endswith(':d0.0025'):
                    self.assertAlmostEqual(r['z'],indexed[r['id'].rsplit(':',1)[0]]['z']-.001)
        frozen_data()  # no source data edits allowed

    def test_edge_requires_sustained_crossing(self):
        rows=[dict(id=f'edge:0:{i}',z=str(i*.1),u=str(v))
              for i,v in enumerate((.7,1.,.98,1.,1.))]
        self.assertAlmostEqual(thickness(rows,1.),.25)
        self.assertIsNone(thickness([dict(id='edge:0:0',z='.022',u='1')],1.))
        self.assertIsNone(thickness([dict(id='edge:0:0',z='.022',u='.9')],1.))

    def test_design_must_conclude_and_match_its_inputs(self):
        with tempfile.TemporaryDirectory() as root:
            root=Path(root);paths=[]
            for surface in ('S','R'):
                for candidate in design(surface):
                    for model in ('standard_original','cls_original'):
                        path=root/surface/candidate['id']/model/'development.json';path.parent.mkdir(parents=True)
                        path.write_text(json.dumps(dict(status='developing',complete=False,config=candidate['parameters'])))
                        paths.append(path)
            self.assertFalse(select(root)['design_executed'])
            for path in paths:
                d=json.loads(path.read_text());d['status']='failed_or_timed_out';path.write_text(json.dumps(d))
            result=select(root);self.assertTrue(result['design_executed']);self.assertEqual(result['selected'],dict(S=None,R=None))
            d=json.loads(paths[0].read_text());d['config']['k_scale']=2.;paths[0].write_text(json.dumps(d))
            with self.assertRaises(ValueError):select(root)

    def test_repeated_environment_setup_does_not_grow(self):
        old=os.environ.get('LD_LIBRARY_PATH')
        try:
            with tempfile.TemporaryDirectory() as root:
                root=Path(root);(root/'geometry-lib').mkdir();(root/'geometry-lib/libSpatialDomains.so').touch()
                environment(root);first=os.environ['LD_LIBRARY_PATH']
                for _ in range(1000):environment(root)
                self.assertEqual(os.environ['LD_LIBRARY_PATH'],first)
        finally:
            if old is None:os.environ.pop('LD_LIBRARY_PATH',None)
            else:os.environ['LD_LIBRARY_PATH']=old


if __name__=='__main__':unittest.main()
