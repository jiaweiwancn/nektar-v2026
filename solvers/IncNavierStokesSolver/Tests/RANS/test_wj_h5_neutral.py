#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from wj_h5_neutral import *
from wj_h5_oracle import check_reference
from wj_h5_cases import case_xml,mesh_knots
from wj_h2_cases import analytic_expressions,inner_normal_knots


class NeutralInlet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.fit=mean_fit()
    def config(self,H=5.5,**kw):return NeutralReference(H/SCALE_MM,self.fit['z0_e'],**kw)

    def test_constants_and_physical_viscosity(self):
        assert_allclose([QE,KAPPA,DEFAULT_NU/DEFAULT_US*SCALE_MM],[3.386780140528718,.381917863546733,.0503778337531],rtol=1e-12)
        c=self.config();self.assertAlmostEqual(c.K,.008809015145515192,places=17)
        assert_allclose(np.diag(c.Rc)/c.u_tau**2,[3.094095437,2.257853427,1.421611417],atol=5e-10,rtol=0)
        self.assertEqual(c.Rc[0,2],-c.u_tau**2)
        self.assertNotAlmostEqual(c.nu,c.u_tau*c.z0_e/3.97,places=8)

    def test_data_origins_and_common_fit(self):
        f=self.fit;self.assertEqual(len(f['points']),3)
        assert_allclose([f['z0_e_mm'],f['rms']],[.227050812848672,.00269823112168],rtol=3e-12,atol=1e-14)
        self.assertEqual(len(mean_fit(omit_lowest=True)['points']),2)
        e=mean_fit(window_mm=50);self.assertEqual(len(e['points']),6);self.assertLess(e['points'][-1]['z_mm'],50.)
        for low,high in f['reading_bounds_mm'].values():self.assertLess(low,f['z0_e_mm']);self.assertGreater(high,f['z0_e_mm'])
        with self.assertRaisesRegex(ValueError,'empty'):mean_fit(window_mm=1)
        from wj_h2_profiles import Completion
        with self.assertRaisesRegex(ValueError,'empty'):Completion(H=.125)

    def test_observation_and_feasibility_conflict(self):
        f=feasibility(self.fit)
        assert_allclose([f['minima'][p]['z_mm']for p in 'bcd'],[10.0889788525,10.1539229009,10.2222203003],atol=5e-11,rtol=0)
        assert_allclose([f['fig3_conservative_ceiling_mm'],f['primary_h_equals_H_mm'],f['printed_h_equals_H_mm']],
            [6.2987042140,4.41351554172,3.80908930644],atol=5e-10,rtol=0)
        self.assertEqual(f['primary_two_percent_interval'],'EMPTY');self.assertEqual(f['all_7mm_outside']['primary'],'EMPTY')
        assert_allclose(f['all_7mm_outside']['primary_required_h_at_3p5mm'],4.486214437,atol=6e-9,rtol=0)
        for H in (5.,5.5,6.):
            regions=observation_regions(self.config(H));self.assertTrue(all(r['region']=='outer WJ'for r in regions if r['panel']in 'abcd'))
            self.assertEqual(regions[-1]['region'],'inner AWF');self.assertEqual(regions[-1]['z_mm'],7.)

    def test_three_reference_profiles_and_d_lt_2nu_extension(self):
        expected=[(4.381467179,1.047926912,2.562624069,1.907080936,1.283037),
                  (4.359758967,1.052746607,2.653345528,1.933718675,1.280847),
                  (4.341801417,1.056795932,2.736517844,1.956862238,1.279006)]
        for H,e in zip((5.,5.5,6.),expected):
            c=self.config(H);v=c.summary()
            assert_allclose([v['h_s_mm'],v['rc_mm'],v['Lambda'],v['k_ground_over_K']],e[:4],atol=8e-9,rtol=0)
            assert_allclose(v['reference_covariance_min_over_utau2'],e[4],atol=6e-7,rtol=0)
            self.assertLess(v['d_over_nu'],2);self.assertTrue(all(v['positivity_extension'].values()))
            self.assertTrue(all(c.w['Q'](r)<0 for r in np.linspace(0,c.H,31)[1:-1]))
            self.assertTrue(all(c.k_inner(r)>=c.K-1e-15 for r in np.linspace(0,c.H,31)))
            self.assertEqual(check_reference(c)['status'],'PASS')

    def test_reject_invalid_or_infeasible_configuration(self):
        for H in (0.,float('nan'),3.5/SCALE_MM):
            with self.assertRaises(ValueError):NeutralReference(H,self.fit['z0_e'])
        with self.assertRaisesRegex(ValueError,'branch'):NeutralReference(.1375,.001)
        with self.assertRaises(ValueError):self.config(nu=0.)
        with self.assertRaises(ValueError):self.config(u_tau=float('inf'))

    def test_full_traction_gradients_pressure_and_epsilon_trace(self):
        c=self.config();H=c.H;G=np.zeros((3,3));G[0,2]=c.u_tau**2/c.D(H)
        piw=-c.K/NE;Tw=c.D(H)*(G+G.T)-piw*np.eye(3);To=-c.Rc+2*c.K/3*np.eye(3)
        assert_allclose(Tw[:,2],To[:,2],atol=3e-18,rtol=0)
        assert_allclose(Tw[:,0],[c.K/NE,0,c.u_tau**2],atol=3e-18,rtol=0)
        assert_allclose(To[:,0],[-c.K/NE,0,c.u_tau**2],atol=3e-18,rtol=0)
        assert_allclose(c.U_inner(H),c.U_outer(H),atol=2e-15,rtol=0)
        assert_allclose(G[0,2]/(c.u_tau/KAPPA/H),1-c.z0_e/H,atol=1e-15,rtol=0)
        assert_allclose(c.epsilon(H),c.destruction(H),atol=1e-17,rtol=0)
        self.assertNotAlmostEqual(c.nt(H)*G[0,2],c.u_tau**2,places=7)

    def test_explicit_parameters_reach_wall_functions_and_observers(self):
        c=self.config(nu=DEFAULT_NU*1.1,u_tau=.052,fit_identity='explicit parameter test')
        v=analytic_expressions(c.H,c.h,u_tau=c.u_tau,nu=c.nu,d0=.1,top=19.,z0_e=c.z0_e)
        self.assertEqual(float(v[1]['k']),c.K);self.assertEqual(float(v[2]['v']),c.u_tau**2)
        self.assertIn('0.10000000000000001',v[1]['u']);self.assertGreater(abs(c.h-self.config().h),1e-4)
        from wj_h4_points import evaluate
        owner=dict(region='inner',chart=[0.,0.,c.H/2]);r=dict(parameters=[c.nu,c.H,c.h],geometry_signature='fixture:CompositeTunnelV1:0.1:19')
        def field(r,s,b,o):return (c.K,np.zeros(3))if b==2 else (0.,np.array([0.,0.,1.])if b==0 else np.zeros(3))
        with patch('wj_h4_points.field',side_effect=field),patch('wj_h4_points.matching',return_value=owner):
            p=evaluate(r,[],[0,0,.1+c.H/2],owner)
        assert_allclose(p['nu_d'],c.D(c.H/2),atol=1e-16,rtol=0)
        self.assertIsNone(p['epsilon'])

    def test_low_H_generation_has_one_interface_and_no_old_completion(self):
        with tempfile.TemporaryDirectory()as folder:
            for H in (5.,5.5,6.):
                c=self.config(H);ik,outer,z=mesh_knots(c)
                self.assertTrue(all(0<r<c.H for r in ik[1:-1]));self.assertEqual(z.count(c.H),1)
                self.assertTrue(all(a<b for a,b in zip(z[:-1],z[1:])))
                with patch('wj_h2_cases.Completion',side_effect=AssertionError('old inlet family must not be used')):
                    path=case_xml(Path(folder)/str(H),c)
                root=ET.parse(path).getroot();design=json.loads(path.with_name('case-design.json').read_text())
                self.assertLess(max(design['expression_errors'].values()),2e-11)
                self.assertIsNone(root.find('./CONDITIONS/SOLVERINFO/I[@PROPERTY="RANSWallRestart"]'))
                self.assertEqual(root.find('./CONDITIONS/SOLVERINFO/I[@PROPERTY="RANSWallReference"]').get('VALUE'),'Functions')
                self.assertIsNone(root.find('./CONDITIONS/FUNCTION[@NAME="TopTractionOuter"]/E[@VAR="v"]'))
                self.assertIsNone(root.find('./CONDITIONS/FUNCTION[@NAME="BoundaryValuesInner"]/E[@VAR="p"]'))
                self.assertFalse('Blackadar'in path.read_text());self.assertIn('NumSteps = 0',path.read_text())
                J=(c.top-c.d0-c.H)/(c.top-c.H)
                assert_allclose([c.d0+c.H+(q-c.H)*J for q in z[len(ik):]],[c.d0+q for q in outer[1:]],atol=4e-15,rtol=0)
        old_w=dict(rc=.02,zv=.00008);self.assertEqual(inner_normal_knots(.5,old_w,.04),[0.,.001,.005,.02,.04,.05,.1,.25,.5])

    def test_exterior_weak_load_uses_projected_test_trace(self):
        from wj_h2_weak import independent_weak
        c=self.config();H=c.H
        el=dict(geometry=3,degree=3,lower=[0,0,0],upper=[1,0,H])
        outer=dict(geometry=4,degree=3,lower=[0,0,H],upper=[1,0,H+.1])
        rule=np.polynomial.legendre.leggauss(24)[0];samples=[]
        test=lambda t:np.array([[1]*4,[.2]*4,[.2+.3*t-.4*t*t+.1*t**3]*4])
        for x in rule:samples.append(dict(role='outlet',geometry=3,chart=[1,0,H*(x+1)/2],modes=test((x+1)/2).tolist()))
        t=.3;flux=np.array([4.,0,2*c.K,0,0,0]);expected=np.zeros((3,6));expected[:,:3]=test(t)[:,:3]*flux[:3]
        row=dict(role='outlet',geometry=3,inner=True,chart=[1,0,H*t],weight=1.,normal=[1,0,0],
            U=[2,0,0],k=c.K,epsilon=1.,G=np.zeros((3,3)).tolist(),K=c.K,
            grad_k=[0,0,0],grad_epsilon=[0,0,0],applied_traction=[0,0,0])
        report=dict(parameters=[c.nu,H,c.h],geometry_signature='fixture:CompositeTunnelV1:0.0875:20',
            fields=[dict(mesh=[el,outer])],signature='test',open_boundary=dict(boundary_rows=[row]))
        packet=dict(alpha_mass=0,samples=[],volume=1,boundary_mode_samples=samples,
            native_inner=expected.tolist(),native_outer=np.zeros((3,6)).tolist(),native_terms=dict(exterior_boundary=expected.tolist()))
        v=independent_weak(packet,report)
        self.assertEqual(v['status'],'PASS');self.assertIn('projected',v['exterior_test_scope'])
        # An unprojected hat vanishes here; it gives the wrong exterior load.
        self.assertNotEqual(expected[2,0],0.)


if __name__=='__main__':unittest.main()
