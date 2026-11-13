#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from numpy.testing import assert_allclose
from wj_h13_domain import X,H,D0,TOP,topology,cg_keys,dof_counts,storage,mapped,certificate,geometry_checks,source_capacity,preflight
from wj_h13_cases import primary,boundary,protocol,prepare
from wj_h5_cases import case_xml,mesh_knots
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_parameter,write

SOURCE=Path(__file__).resolve().parents[4]

def tiny(nx=2):
    r=ET.fromstring('<NEKTAR><GEOMETRY/><EXPANSIONS/><CONDITIONS/></NEKTAR>')
    xs=np.linspace(-10,15,nx+1).tolist();tensor_mesh(r,[0.,.05,H,1.,TOP],3,2,25.,1.,points=[8,8],xknots=xs)
    return r,topology(r)

class FullHeightPairedDomain(unittest.TestCase):
    def test_frozen_reference_full_height_and_molecular_anchor(self):
        c,v=primary(SOURCE);self.assertEqual(c.K,.008809015145515192);self.assertEqual(c.nu,.051*.005/3.97)
        self.assertEqual(c.h,v['frozen_primary']['h_s']);self.assertGreater(c.U_outer(TOP-D0),1.)
        self.assertEqual(len(v['historical_heights']),3);self.assertTrue(v['no_new_fit'])
        self.assertNotEqual(c.z0_e,.2/40)
    def test_complete_traction_pressure_and_unspecialized_epsilon_budget(self):
        c,_=primary(SOURCE);v=boundary(c);B=v['top']['B'];rT=TOP-D0
        self.assertAlmostEqual(B,v['top']['epsilon_equation_coefficient'],places=20)
        for factor in (.7,1.,1.4):
            cd=factor*572/6561;self.assertAlmostEqual(cd*c.K*c.K/1.3,factor*c.u_tau**4/1.3,places=20)
        self.assertGreater(v['top']['epsilon_outward_diffusive_flux'],0.)
        self.assertEqual(v['top']['epsilon_signed_gradient_flux'],-v['top']['epsilon_outward_diffusive_flux'])
        self.assertEqual(v['outlet']['inner'],[c.K/(81/20),0.,c.u_tau**2])
        self.assertEqual(v['outlet']['outer'],[-c.K/(81/20),0.,c.u_tau**2])
        # Derive normal/shear traction independently from full stress tensors.
        G=np.zeros((3,3));G[0,2]=c.u_tau**2/c.D(H)
        Tw=c.D(H)*(G+G.T)-v['pressure']['inner']*np.eye(3)
        To=2*c.K/3*np.eye(3)-c.Rc
        assert_allclose(Tw[:,2],To[:,2],rtol=0,atol=2e-18)
        assert_allclose(To[:,0],v['outlet']['outer'],rtol=0,atol=2e-18)
        gradient=-c.A/rT**2;current_D=(572/6561)*(1.1*c.K)**2/(1.3*(.9*c.A/rT))
        self.assertGreater(abs(-current_D*gradient-v['top']['epsilon_outward_diffusive_flux']),1e-8)
    def test_paired_composite_foot_reach_and_interface(self):
        _,t=tiny();self.assertEqual(geometry_checks(t,False)['status'],'PASS')
        # Small geometry fixture includes explicit feet, separate from target.
        r=ET.fromstring('<NEKTAR><GEOMETRY/><EXPANSIONS/><CONDITIONS/></NEKTAR>');tensor_mesh(r,[0.,H,1.,TOP],3,2,25.,1.,xknots=[-10,-2.5,0,2.5,15]);t=topology(r)
        check=geometry_checks(t,True);self.assertLess(check['opposite_cofactor_error'],1e-12)
        a=mapped(1.,0.,True);self.assertGreater(abs(a['physical'][2]-(math.cos(math.pi/5)**2+D0)),1e-3)
        for foot in check['foot_limits']:
            assert_allclose(foot['owners'][0]['physical'],foot['owners'][1]['physical'],atol=1e-15)
            self.assertNotEqual(foot['owners'][0]['kappa_phys'],foot['owners'][1]['kappa_phys'])
        with self.assertRaises(ValueError):certificate(True,H=2.)
        with self.assertRaises(ValueError):mapped(0.,2.,True,region='inner')
    def test_target_count_entities_and_epsilon_identity_rows(self):
        for nx,expected in ((1,2328),(2,4048),(14,24688)):
            c,_=primary(SOURCE);_,_,z=mesh_knots(c,bisect_once=True);r=ET.fromstring('<NEKTAR><GEOMETRY/><EXPANSIONS/><CONDITIONS/></NEKTAR>')
            tensor_mesh(r,z,3,2,25.,1.,xknots=np.linspace(-10,15,nx+1));n=dof_counts(topology(r));self.assertEqual(n['state'],expected)
            self.assertEqual(n['free_solve_blocks']['epsilon']-n['independent']['epsilon'],3*nx+1)
            self.assertEqual(n['fields']['pi_inner'],(nx+1)*17);self.assertEqual(n['fields']['pi_outer'],(nx+1)*35)
            if nx==1:self.assertEqual((n['free_solve_blocks']['momentum_predictor'],n['free_solve_blocks']['pressure']),(897,104))
        self.assertEqual(n['free_solve_blocks']['epsilon'],4327)
        self.assertGreater(storage(n)['simultaneous_Pimple_payload_lower_bound'],4*1024**3)
    def test_multi_column_trace_identity_without_projection(self):
        # Pure coefficient injection on three columns, orientation signs
        # explicit; no reference reevaluation or radial/mass projection.
        from validate_cls_epawf_m3 import basis
        radial=np.array([.2,.4,.07,-.013]);sign=np.array([1.,-1.,1.,-1.]);lifted=radial*sign
        recovered=sign*lifted
        for column in range(3):
            coeff=np.zeros((4,4));coeff[:,0]=recovered;coeff[:,1]=recovered
            px,dx=basis(3,[-1.,-.3,.2,1.]);pr,dr=basis(3,[-1.,0.,1.])
            for j in range(3):
                for i in range(4):
                    self.assertAlmostEqual(np.outer(pr[j],px[i]).ravel()@coeff.ravel(),pr[j]@radial,places=15)
                    self.assertAlmostEqual(np.outer(pr[j],dx[i]).ravel()@coeff.ravel(),0.,places=15)
            assert_allclose(coeff[:,0]*sign,lifted,rtol=0,atol=0)
        self.assertEqual(len(cg_keys(3,2,3)),(3*3+1)*(3*2+1))
    def test_generator_preserves_defaults_and_primary_physical_functions(self):
        c,_=primary(SOURCE)
        with tempfile.TemporaryDirectory()as d:
            d=Path(d);a=case_xml(d/'legacy',c,bisect_once=True);b=case_xml(d/'paired',c,bisect_once=True,xknots=X)
            fa=lambda p:[ET.tostring(v)for v in ET.parse(p).findall('./CONDITIONS/FUNCTION')]
            self.assertEqual(fa(a),fa(b));self.assertEqual(json.loads((a.parent/'case-design.json').read_text())['predicted_state_size'],2328)
            self.assertEqual(json.loads((b.parent/'case-design.json').read_text())['predicted_state_size'],24688)
    def test_observation_scope_and_duration_not_a_preservation_pass(self):
        c,_=primary(SOURCE);_,t=tiny();p=protocol(c,t);stations=p['stations'];self.assertIn(7.5,[v['x_code']for v in stations])
        self.assertEqual(p['status'],'PREPARED_NOT_ADMITTED');self.assertIsNone(p['stages']['B']['application_error_budgets'])
        self.assertGreater(p['stages']['B']['clocks']['proposed_preservation_horizon'],p['stages']['B']['clocks']['T_slowest_outer_transit'])
        self.assertFalse(p['old_checkpoint_conversion'])
    def test_external_preflight_refuses_before_any_native_allocation(self):
        from unittest.mock import patch
        with tempfile.TemporaryDirectory()as d:
            m=prepare(Path(d)/'case',SOURCE);p=m['flat_session']['path']
            with patch('subprocess.Popen',side_effect=AssertionError('native launch prohibited')):
                v=preflight(p,m,SOURCE)
            self.assertEqual(v['status'],'OBSTRUCTED');self.assertFalse(v['native_launch_allowed']);self.assertTrue(v['NumSteps_zero_does_not_bypass_capacity'])
            r=ET.parse(p).getroot();set_parameter(r,'NumSteps',1);write(r,p)
            with self.assertRaises(ValueError):preflight(p,m,SOURCE)
    def test_no_native_or_build_command_in_offline_adapter(self):
        from run_wj_h13 import offline_intent
        for command in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','cmake','ninja'):
            with self.assertRaises(AssertionError):offline_intent([command])
        with self.assertRaises(AssertionError):offline_intent(['python3','study_wj_h13.py','--execute'])
        self.assertIsNone(offline_intent(['python3','-B','study_wj_h13.py']))

if __name__=='__main__':unittest.main()
