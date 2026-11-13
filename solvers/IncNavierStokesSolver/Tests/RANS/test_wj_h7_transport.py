#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import tempfile
from pathlib import Path
import unittest
import numpy as np
from numpy.testing import assert_allclose
from wj_h3_analysis import principal
from wj_h7_analysis import principal_batch,classify,context
from run_wj_h7 import MATRIX,admit_xml
from study_wj_h7 import clone,compatible,accounting
from wj_h4_points import vector_hash,aligned
from wj_h5_cases import case_xml
from study_wj_h6 import configuration

class TransportTests(unittest.TestCase):
    def test_matrix_limits_and_clone(self):
        self.assertEqual(sum(n for h,dt,n in MATRIX),36);self.assertEqual(len(MATRIX),6)
        with tempfile.TemporaryDirectory()as temp:
            root=Path(temp);src=case_xml(root/'parent',configuration(5.5),bisect_once=True)
            cp=root/'fine-cp0';cp.mkdir();xml=clone(src,cp/'case.xml',.005,0)
            self.assertEqual(admit_xml(xml,True)[1],0)
            run=root/'run';run.mkdir();xml=clone(src,run/'case.xml',.01,4,root/'CP')
            self.assertEqual(admit_xml(xml)[1],4)
            clone(src,xml,.01,5,root/'CP')
            with self.assertRaises(AssertionError):admit_xml(xml)
    def test_cross_dt_and_state_signature_rejection(self):
        a=dict(dt=.01,state=[0.,1.],fixed_reference=[0.,1.],signature='a',steps=0,time=0)
        compatible(a,a,.01)
        with self.assertRaisesRegex(ValueError,'timestep'):compatible(a,a,.005)
        with self.assertRaisesRegex(ValueError,'state'):compatible(a,dict(a,state=[-0.,1.]),.01)
        with self.assertRaisesRegex(ValueError,'signature'):compatible(a,dict(a,signature='b'),.01)
        compatible(a,dict(a,dt=.005,signature='fine'),.005,same_dt=False)
    def test_integer_time_and_odd_steps(self):
        self.assertEqual(aligned(6,.005,.03+2e-17),aligned(3,.01,.03))
        self.assertNotEqual(aligned(5,.005,.025),aligned(3,.01,.03))
        self.assertEqual(len([s for s in range(9)if s%2]),4)
        with self.assertRaises(ValueError):aligned(2.5,.005,.0125)
    def test_genuine_context_and_rejected_binding(self):
        s=[1.,2.];h=[.9,1.8];p=dict(state_sha256=vector_hash(s),stage_history_sha256=vector_hash(h),callback_history=h,alpha_mass=100.,time=.01)
        r=dict(dt=.01,steps=1,accepted_states=[s],wj_audit_packets=[p],physical_range_events=[dict(state=s,history=h,alpha=100.,time=.02,wj_packet=dict(p,time=.02))])
        self.assertEqual(context(r)['role'],'accepted');self.assertEqual(context(r,rejected=True)['step'],2)
        with self.assertRaises(ValueError):context(r,0)
        with self.assertRaisesRegex(ValueError,'binding'):context(dict(r,accepted_states=[[2.,2.]]))
    def test_batched_full_action_matches_existing(self):
        G=np.array([[[.01,.03,.9],[.02,-.03,.1],[-.02,.01,.02]],[[.02,0,.7],[0,.01,0],[.01,0,-.04]],np.zeros((3,3))])
        k=np.array([.02,.009,.01]);e=np.array([.01,.002,.001]);_,_,_,m,t,a=principal_batch(k,e,G)
        for i in range(3):
            ref=principal(k[i],e[i],G[i]);assert_allclose([m[i],t[i]],[ref['in_plane'],ref['transverse']],atol=2e-11,rtol=0)
    def test_fixed_scale_zero_safe_norms(self):
        from wj_h7_analysis import Rows
        from types import SimpleNamespace
        r=Rows.__new__(Rows);r.rows=[0,1,2];r.config=SimpleNamespace();r.weights=np.array([.5,.5,1.]);r.outer=np.array([False,True,True])
        a={k:np.zeros((3,3))if k=='U'else np.zeros((3,3,3))if k in ('G','R')else np.zeros(3)for k in ('U','G','R','k','epsilon','pi','mechanical_pressure')};b={k:v.copy()for k,v in a.items()};a['U'][:,0]=2
        self.assertAlmostEqual(r.norms(a,b)['full']['U']['L2'],np.sqrt(8));self.assertEqual(r.norms(a,b)['full']['epsilon']['L2'],0)
    def test_bound_calls_are_not_activations(self):
        a=dict(scheduled_complete=True,momentum_predictor_solves=1,pressure_system_solves=1,k_transport_solves=1,epsilon_transport_solves=1,bound_calls=1,bound_activations=0)
        r=dict(pimple_attempts=[a],steps=1,failure='');self.assertEqual(accounting(r)['bound_calls'],1)
        with self.assertRaises(AssertionError):accounting(dict(r,pimple_attempts=[dict(a,bound_activations=1)]))
        partial=dict(a,scheduled_complete=False,k_transport_solves=0,bound_calls=0)
        self.assertEqual(accounting(dict(r,pimple_attempts=[partial],steps=0,failure='original'))['failure'],'original')
    def test_quality_and_physical_status_are_distinct(self):
        self.assertEqual(classify(.003,1,.001,.001),'PASS');self.assertEqual(classify(.001,1,.001,.001),'QUALITY_LOST')
        self.assertEqual(classify(-1e-3,1,.001,.001),'OBSTRUCTED');self.assertEqual(classify(-1e-12,1,.001,.001),'UNRESOLVED')
        self.assertEqual(classify(float('nan'),1,1,1),'UNRESOLVED')

if __name__=='__main__':unittest.main()
