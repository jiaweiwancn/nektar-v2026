#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import tempfile
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h9 import admit_command
from wj_h9_local import MATRIX,phase_bindings,delta_B,telescope,safe_ratio,same,MomentAudit,actions,DIRECTION
from wj_h9_weak import reaction,pressure_terms,ProjectedTests,shapes
from wj_h7_analysis import Rows
from wj_h4_points import vector_hash
from wj_h5_neutral import NeutralReference,QE,KAPPA
from wj_h6_projection import project_field
from test_wj_h6_projection import fixture


def phases():
    P=np.array([1.,2.,3.,4.,5.,6.]);C=P.copy();F=C.copy();F[[0,1,4,5]]+=[.1,.2,.3,.4];K=F.copy();K[2]+=.1;A=K.copy();A[3]+=.2;history=P.copy()
    packet=dict(time=.01,alpha=100.,committed=C.tolist(),post_flow=F.tolist(),pre_k=F.tolist(),raw=K.tolist(),corrected=K.tolist(),history=history.tolist(),fixed_reference=P.tolist(),k_bound=dict(activations=0,mass_norm_correction=0))
    audit=lambda s,t,a:dict(state_sha256=vector_hash(s),stage_history_sha256=vector_hash(history),callback_history=history.tolist(),alpha_mass=a,time=t)
    report=dict(offsets=list(range(7)),dt=.01,steps=1,accepted_states=[P.tolist(),A.tolist()],wj_audit_packets=[audit(P,0,0),audit(A,.01,100.)],fixed_reference=P.tolist())
    return report,packet

class StartupAuditTests(unittest.TestCase):
    def test_fixed_seven_attempt_zero_native_scope(self):
        self.assertEqual(len(MATRIX),7);self.assertEqual([n for b,dt,n in MATRIX if b=='quarter'],[1,2,3,4])
        self.assertTrue(all(dt*n<=.01 for b,dt,n in MATRIX));admit_command(['python3','-B','study_wj_h9.py'])
        for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','cmake','make','WJKernelProbe'):
            with self.assertRaises(ValueError):admit_command([name,'case.xml'])
    def test_real_phase_masks_and_aliases(self):
        r,p=phases();states,b=phase_bindings(r,p,1);self.assertTrue(b['P_C_bitwise']);self.assertEqual(b['field_masks'],{'C->F':[0,1,4,5],'F->K':[2],'K->A':[3]})
        p['committed'][0]+=.01;states,b=phase_bindings(r,p,1);self.assertFalse(b['P_C_bitwise']);self.assertAlmostEqual(b['lift_max'],.01)
        r,p=phases();p['post_flow'][2]+=.1;p['pre_k']=p['post_flow'].copy()
        with self.assertRaises(AssertionError):phase_bindings(r,p,1)
    def test_immutable_history_reference_and_rejected_lineage(self):
        r,p=phases();p['history'][0]+=.1
        with self.assertRaises(AssertionError):phase_bindings(r,p,1)
        r,p=phases();p['fixed_reference'][0]=-0.
        with self.assertRaises(AssertionError):phase_bindings(r,p,1)
        r,p=phases();r['steps']=0
        with self.assertRaises(ValueError):phase_bindings(r,p,1)
        r,p=phases();p['corrected'][2]+=.01
        with self.assertRaises(AssertionError):phase_bindings(r,p,1)
    def test_intermediate_trace_is_not_repaired(self):
        H=.15;K=.0088;newK=K*1.001;old_e=(K/QE)**1.5/(KAPPA*H);new_e=(newK/QE)**1.5/(KAPPA*H)
        role=dict(stage='K',accepted=False,epsilon=old_e,target=new_e)
        self.assertGreater(abs(role['epsilon']-role['target']),1e-7);self.assertEqual(role['epsilon'],old_e)
        roleA=dict(stage='A',accepted=True,epsilon=new_e,target=new_e);self.assertEqual(roleA['epsilon']-roleA['target'],0)
    def test_all_Q1_moment_linearity_and_pressure_independence(self):
        r=fixture();r['geometry_signature']='fixture:CompositeTunnelV1:0.2:2';c=NeutralReference(.5,.0056762703212168,d0=.2,top=2)
        rows=Rows(r,[[8,0,.7]],c);audit=MomentAudit(rows);C=np.zeros(len(r['state']));co,_=project_field(r,0,lambda x,z:x*x);C[:len(co)]=co;F=np.zeros_like(C)
        result=audit.identity(C,F)
        for v in result.values():self.assertLess(v['linear_identity_max'],1e-12);self.assertLess(v['repair_max'],1e-12);self.assertGreater(max(abs(np.array(v['C']))),.1)
        shifted=C.copy();shifted[r['offsets'][4]:]+=3.
        for key,v in audit.moments(C,False).items():assert_allclose(v['moments'],audit.moments(shifted,False)[key]['moments'],atol=0,rtol=0)
    def test_signed_impulses_and_no_regauging(self):
        pattern=np.array([-2.,1.,4.]);seed=np.array([.3,-.2,.1]);impulses=[]
        for dt in (.01,.005,.0025):impulses.append(dt*((seed+pattern/dt)-seed))
        for v in impulses:assert_allclose(v,pattern,rtol=0,atol=1e-15)
        self.assertGreater(abs(np.mean(pattern)),.1);self.assertIsNone(safe_ratio(0,0));self.assertEqual(safe_ratio(2,4),.5)
    def test_exact_delta_B_has_two_nonzero_terms(self):
        G=np.array([[.1,.2,.7],[-.1,.3,.2],[.4,.2,-.1]]);a=dict(k=.2,epsilon=.1,G=G);b=dict(k=.3,epsilon=.12,G=G+.03*np.arange(9).reshape(3,3));v=delta_B(a,b)
        self.assertGreater(np.linalg.norm(v['gradient_term']),.01);self.assertGreater(np.linalg.norm(v['timescale_term']),.01);assert_allclose(v['residual'],0,atol=5e-16,rtol=0)
    def test_telescoping_nonzero_functional_increments(self):
        values=np.array([[1.,2.],[2.,4.],[3.,1.],[2.,5.]]);functional=values[:,0]**2+values[:,1]**3;pieces,error=telescope(functional)
        self.assertTrue(all(abs(pieces)>0));self.assertEqual(error,0)
    def test_pressure_roughness_and_reaction_signs(self):
        pi=np.array([2.]);grad=np.array([[[.3,0,.4]]]);phi=np.array([[.7]]);v,f=pressure_terms(pi,grad,phi,np.array([.5]),np.array([True]))
        assert_allclose(v,[[[-.6,0,-.8]]],atol=0,rtol=0);assert_allclose(f,[[[-.35,0,0]]],atol=0,rtol=0)
        assert_allclose(pressure_terms(pi,grad,phi,np.array([.5]),np.array([False]))[1],0,atol=0,rtol=0)
        T=np.array([np.diag([1.,2.,3.])]);shift=T-2*np.eye(3);n=np.array([0.,0,1.]);given=np.array([[.1,0,0]])
        assert_allclose(-(reaction('top',shift,n,given)-reaction('top',T,n,given)),[2*n],atol=0,rtol=0)
        assert_allclose(reaction('outlet',shift,n,given),reaction('outlet',T,n,given),atol=0,rtol=0)
    def test_final_traction_reuse_would_fail_intermediate(self):
        F=np.array([np.diag([2.,1.,3.])]);A=np.array([np.diag([4.,1.,7.])]);n=np.array([-1.,0,0]);given=np.zeros((1,3))
        expected=np.array([[-2.,0,0]]);assert_allclose(reaction('inlet',F,n,given),expected,atol=0,rtol=0)
        self.assertGreater(np.max(abs(reaction('inlet',A,n,given)-expected)),1)
    def test_projected_polynomial_tests_and_gradients(self):
        r=fixture();H=.5;J=(20-.0875-H)/(20-H);samples=[];boundary=[];g=np.polynomial.legendre.leggauss(4)[0]
        for el in r['fields'][0]['mesh']:
            geom=el['geometry'];inner=el['upper'][2]<=H;lo,hi=el['lower'],el['upper'];xx,rr=np.meshgrid(4+4*g,(lo[2]+hi[2])/2+(hi[2]-lo[2])/2*g)
            for x,z in zip(xx.ravel(),rr.ravel()):
                modes=[]
                for m in range(3):
                    v=1. if m==0 else x/5 if m==1 else x*z+.3*z*z+.02*x**3
                    gv=[0.,0,0]if m==0 else [.2,0,0]if m==1 else [z+.06*x*x,0,(x+.6*z)/(1 if inner else J)]
                    ev=1. if m==0 else x/5 if m==1 else .5*x*z+.1*z**3
                    eg=[0.,0,0]if m==0 else [.2,0,0]if m==1 else [.5*z,0,(.5*x+.3*z*z)/J]
                    modes.append(dict(v=v,grad_v=gv,p=1. if m==0 else x/5 if m==1 else x*z,e=0 if inner else ev,grad_e=[0,0,0]if inner else eg))
                samples.append(dict(geometry=geom,chart=[x,0,z],inner=inner,modes=modes))
        t=ProjectedTests(r,dict(samples=samples,boundary_mode_samples=[]));v,d=t.values(2,'v',np.array([3.]),np.array([1.2]));self.assertAlmostEqual(v[2,0],3*1.2+.3*1.2**2+.02*27,places=12);self.assertGreater(v[2,0],1)
    def test_saved_k_matvec_without_a_solve(self):
        from wj_h9_local import matrix_check
        M=np.array([[2.,.1],[.2,3.]]);delta=np.array([1e-6,-2e-6]);F=np.array([.01,.02]);K=F+delta
        packet=dict(transport_rows=[0,1],raw_increment=delta.tolist(),transport_residual=(-M@delta).tolist(),raw_linear_residual=[0.,0.],raw=K.tolist(),corrected=K.tolist(),pre_k=F.tolist())
        with tempfile.TemporaryDirectory()as folder:
            p=Path(folder)/'matrix.bin';M.tofile(p)
            with patch('numpy.linalg.solve',side_effect=AssertionError('no governing solve')):
                v=matrix_check(packet,p)
            self.assertEqual(v['matvecs'],3);self.assertEqual(v['solves'],0)
            (M+np.eye(2)).tofile(p)
            with self.assertRaises(AssertionError):matrix_check(packet,p)

    def test_local_only_interface_packet_does_not_invent_wall_samples(self):
        from types import SimpleNamespace
        from study_wj_h6 import configuration
        c=configuration(6.);v=Rows.__new__(Rows);v.config=c;v.rows=[([4,0,c.d0+c.H],dict(chart=[4,0,c.H],interface_role='matching_inner',region='inner')),([4,0,c.d0+c.H],dict(chart=[4,0,c.H],interface_role='matching_outer',region='outer'))];v.point_count=2;v.outer=np.array([False,True]);v.weights=np.zeros(2)
        d=dict(K=np.full(2,c.K),U=np.zeros((2,3)),k=np.full(2,c.K),epsilon=np.full(2,c.A/c.H),T=np.zeros((2,3,3)),nu_d=np.ones(2),grad_k=np.zeros((2,3)),pi=np.zeros(2),turbulent_production=np.zeros(2),production=np.zeros(2),reference_production=np.zeros(2),destruction=np.zeros(2))
        p=v.interface_history(d);self.assertIsNone(p['wall_shear_range']);self.assertIn('not_evaluated',p['wall_shear_scope']);self.assertEqual(len(p['rows']),1)

    def test_anchor_categories_keep_physical_and_derived_loads_distinct(self):
        from wj_h9_weak import compare_anchor
        value=np.zeros((2,3,6));value[0,0,0]=1.;physical=dict(components={'time_mass':value.sum(axis=0).tolist()},regional=value.tolist())
        native=dict(components={'time_mass':value.sum(axis=0).tolist(),'total_unmasked':value.sum(axis=0).tolist(),'spatial_sum':np.zeros((3,6)).tolist(),'total_constrained':np.full((3,6),99.).tolist()},regional=value.tolist())
        result=compare_anchor(physical,native);self.assertEqual(result['status'],'PASS');self.assertIn('total_constrained',result['separate_native_diagnostic_roles'])
        native['components']['time_mass'][0][0]=2.
        self.assertEqual(compare_anchor(physical,native)['status'],'UNRESOLVED')

    def test_actual_root_derivative_and_fixed_direction(self):
        self.assertLess(abs(np.trace(DIRECTION)),1e-16)
        G=np.array([[[.01,0,.8],[0,0,0],[.001,0,-.002]]]);v=actions(np.array([.009]),np.array([.002]),G)
        self.assertLess(abs(v['implicit_derivative_error'][0]),2e-11);self.assertGreater(v['F_N'][0],0);self.assertAlmostEqual(v['dimensional_fixed_action'][0]/v['nu_eff'][0],v['normalized_fixed_action'][0],places=14)

if __name__=='__main__':unittest.main()
