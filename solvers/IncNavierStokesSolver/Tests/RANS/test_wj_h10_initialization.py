#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
from wj_h10_analysis import MATRIX,QUALITY,gates,residual_check,arithmetic,remaining_matrix
from wj_h9_local import same
from run_wj_h10 import native_intent,runtime_snapshot
from cls_epawf_m3_cases import set_info,set_parameter,write
import xml.etree.ElementTree as ET

class FlatInletInitialization(unittest.TestCase):
    def test_three_first_steps_only(self):
        self.assertEqual(MATRIX,(('coarse',.01),('fine',.005),('quarter',.0025)))
        self.assertEqual(len(MATRIX),3)
        from study_wj_h10 import small_fixture_tree
        tree=small_fixture_tree(ET.parse(Path(__import__('wj_h5_cases').__file__).with_name('EquilibriumWallFlat.xml')).getroot())
        variables={v.text for v in tree.findall('./CONDITIONS/VARIABLES/V')}
        self.assertIn('epsilon',variables)
        for e in tree.findall('./EXPANSIONS/E'):self.assertTrue(set(e.get('FIELDS').split(','))<=variables)
        self.assertEqual(sum(1 for _ in MATRIX),3)
    def test_quality_and_compatibility_are_separate(self):
        s=dict(status='PASS',min_in_plane=.003,divergence_scaled_max=1e-16,min_covariance=.001,min_transverse=.1)
        m={'4':dict(moments=[0.]*34),'5':dict(moments=[0.]*70)}
        self.assertEqual(gates(s,m,startup=True),'PASS');m['4']['moments'][3]=3e-15
        self.assertEqual(gates(s,m,startup=True),'OBSTRUCTED');m['4']['moments'][3]=0
        s['divergence_scaled_max']=2e-12;self.assertEqual(gates(s,m,startup=True),'OBSTRUCTED');self.assertEqual(gates(s,m,startup=False),'PASS')
        s['min_in_plane']=QUALITY-1e-9;self.assertEqual(gates(s,m,startup=False),'QUALITY_LOST')
        s['status']='UNRESOLVED';self.assertEqual(gates(s,m,startup=False),'UNRESOLVED')
    def test_nonfinite_cannot_pass(self):
        s=dict(status='PASS',min_in_plane=.003,divergence_scaled_max=float('nan'),min_covariance=.1,min_transverse=.1)
        self.assertEqual(gates(s,{'4':dict(moments=[0.])},startup=True),'UNRESOLVED')
    def test_all104_moments_required(self):
        s=dict(status='PASS',min_in_plane=.003,divergence_scaled_max=0.,min_covariance=.1,min_transverse=.1)
        self.assertEqual(gates(s,{'4':dict(moments=[0.]*34)},startup=True),'UNRESOLVED')
        self.assertEqual(gates(s,{'4':dict(moments=[0.]*34),'5':dict(moments=[float('nan')]*70)},startup=True),'UNRESOLVED')
    def test_signed_zero_preservation_is_bitwise(self):
        self.assertFalse(same([0.,1.],[-0.,1.]));self.assertTrue(same([-0.,1.],[-0.,1.]))
    def test_saved_predictor_correction_nonzero(self):
        BG=np.array([[2.,-1.],[.5,3.]]);dp=np.array([.2,-.4]);star=np.array([.6,.7]);du=star-dp@BG
        check=residual_check(du,star-dp@BG,abs(star)+abs(dp@BG));self.assertEqual(check['status'],'PASS')
        with self.assertRaises(AssertionError):residual_check(du,star+dp@BG,abs(star)+abs(dp@BG))
    def test_no_pressure_regauge_in_signed_impulse(self):
        a=np.array([.2,-.3,.7]);b=np.array([.1,-.2,.5]);j=.01*(b-a)
        self.assertNotEqual(float(np.mean(j)),0.);np.testing.assert_allclose(j,.01*b-.01*a,rtol=1e-15,atol=1e-18)
    def test_immutable_dt_and_step_admission(self):
        from run_wj_h4 import sha
        with tempfile.TemporaryDirectory()as t:
            r=Path(t);folder=r/'science/cases/coarse/cp0';folder.mkdir(parents=True);x=ET.fromstring('<NEKTAR><CONDITIONS><SOLVERINFO/><PARAMETERS/></CONDITIONS></NEKTAR>')
            set_info(x,'RANSWallInitialVelocity','FlatInletTraceV1');set_info(x,'RANSWallRestart','None');set_parameter(x,'TimeStep',.01);set_parameter(x,'NumSteps',0);p=folder/'case.xml';write(x,p)
            v=dict(branch='coarse',dt=.01,cp0=dict(path=str(p),sha256=sha(p)),run=dict(path='unused',sha256='unused'));(r/'science/input-intent.json').write_text(json.dumps(dict(cases=[v])))
            self.assertEqual(native_intent(r,[str(r/'runtime/IncNavierStokesSolver'),str(p)],'cp0'),'cp0:coarse:cp0')
            set_parameter(x,'TimeStep',.005);write(x,p)
            with self.assertRaises(AssertionError):native_intent(r,[str(r/'runtime/IncNavierStokesSolver'),str(p)],'cp0')
    def test_zero_relaunch_resume_scope(self):
        b=dict(cp0_reservations=1,governing_attempts=1,strict_readbacks=1)
        self.assertEqual(remaining_matrix(b,True),MATRIX[1:])
        for changed in (dict(b,strict_readbacks=2),dict(b,h10_halt={'status':'UNRESOLVED'}),dict(b,governing_attempts=0)):
            with self.assertRaises(AssertionError):remaining_matrix(changed,True)
    def test_actual_arithmetic_accepts_decoded_list_rhs(self):
        from wj_h4_points import vector_hash
        c=np.zeros(108);c[2:4]=1;f=c.copy();f[4:6]=[.2,-.4];k=f.copy();k[2]=1.2;a=k.copy();a[3]=1.3
        h=vector_hash(c);packets=[dict(state_sha256=h,time=0.,alpha_mass=0),dict(state_sha256=vector_hash(a),time=.01,alpha_mass=100.,callback_history=c.tolist(),stage_history_sha256=h)]
        counters={x:1 for x in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')};counters['bound_activations']=0
        r=dict(offsets=[0,1,2,3,4,38,108],dt=.01,steps=1,accepted_states=[c.tolist(),a.tolist()],wj_audit_packets=packets,fixed_reference=c.tolist(),initialization_lift_calls_iterations=[0,0],initial_velocity_calls=0,pimple_attempts=[counters])
        snap=dict(time=.01,alpha=100.,history=c.tolist(),fixed_reference=c.tolist(),committed=c.tolist(),post_flow=f.tolist(),pre_k=f.tolist(),raw=k.tolist(),corrected=k.tolist(),k_bound=dict(activations=0,mass_norm_correction=0))
        dp=np.zeros(104);dp[:2]=[.2,-.4];BG=np.zeros((104,2));BG[:2]=np.eye(2);S=np.zeros((104,104));S[:2,:2]=np.eye(2)
        mom=dict(rows=[0,1],solution=[.2,-.4],rhs=[20.,-40.],linear_residual=[0.,0.]);p=dict(rows=list(range(4,108)),solution=dp.tolist(),rhs=dp.tolist(),linear_residual=[0.]*104)
        class Moments:
            def moments(self,x,solve=True):return {'4':dict(moments=[float(x[0]),float(x[1])]+[0.]*32),'5':dict(moments=[0.]*70)}
        result=arithmetic(r,snap,mom,p,np.eye(2),BG,S,100*np.eye(2),Moments());self.assertEqual(result['status'],'PASS');self.assertEqual(result['total_velocity_increment'],[0.,0.])
    def test_native_runtime_snapshot_excludes_preexec_union(self):
        root=Path('/h10');host='/usr/lib/libc.so.6';private='/h10/runtime/glibc/libc.so.6';exe='/h10/runtime/RANSWallResolutionProbe.elf'
        self.assertIsNone(runtime_snapshot(root,[host]))
        self.assertIsNone(runtime_snapshot(root,[exe]))
        self.assertEqual(runtime_snapshot(root,[exe,private])['status'],'PASS')
        self.assertEqual(runtime_snapshot(root,[exe,host,private])['status'],'OBSTRUCTED')
    def test_rejected_phase_is_not_a_strict_packet(self):
        with tempfile.TemporaryDirectory()as t:
            r=Path(t);(r/'science').mkdir();(r/'science/input-intent.json').write_text(json.dumps(dict(cases=[])));p=r/'bad.json';p.write_text(json.dumps(dict(operation='candidate')))
            with self.assertRaises(AssertionError):native_intent(r,[str(r/'runtime/RANSWallResolutionProbe'),'missing.xml',str(p)],'strict')

if __name__=='__main__':unittest.main()
