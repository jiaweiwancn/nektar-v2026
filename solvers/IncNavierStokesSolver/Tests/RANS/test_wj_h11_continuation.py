#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
import json
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET
import numpy as np
from wj_h11_analysis import MATRIX,VIEWS,admit_xml,selected,checkpoint_match,ratio,native_reserve,common_step
from wj_h10_analysis import arithmetic
from wj_h4_points import vector_hash
from run_wj_h7 import complete_operation,record_coverage_halt
from cls_epawf_m3_cases import set_parameter,write


def report():
    c=np.zeros(108);c[2:4]=1;f=c.copy();f[4:6]=[.2,-.4];k=f.copy();k[2]=1.2;a=k.copy();a[3]=1.3;h=vector_hash(c)
    states=[c.tolist(),c.tolist(),c.tolist(),a.tolist()];packets=[]
    for i,s in enumerate(states):packets.append(dict(state_sha256=vector_hash(s),time=(i+1)*.01,alpha_mass=100. if i else 0.,callback_history=c.tolist(),stage_history_sha256=h))
    counters={x:1 for x in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls')};counters['bound_activations']=0
    r=dict(offsets=[0,1,2,3,4,38,108],dt=.01,steps=4,accepted_states=states,wj_audit_packets=packets,fixed_reference=c.tolist(),initialization_lift_calls_iterations=[0,0],initial_velocity_calls=0,pimple_attempts=[None,None,counters],signature='active')
    snap=dict(time=.04,alpha=100.,history=c.tolist(),fixed_reference=c.tolist(),committed=c.tolist(),post_flow=f.tolist(),pre_k=f.tolist(),raw=k.tolist(),corrected=k.tolist(),k_bound=dict(activations=0,mass_norm_correction=0))
    return r,snap


class CompatibleContinuation(unittest.TestCase):
    def test_fixed_additional_matrix(self):
        self.assertEqual([n for b,dt,n,target in MATRIX],[3,7,15]);self.assertEqual(sum(n for b,dt,n,target in MATRIX),25)
        self.assertEqual(sum(len(v)for v in VIEWS.values()),11)
        self.assertTrue(all(abs(target*dt-.04)<1e-15 for b,dt,n,target in MATRIX))
    def test_passive_cannot_bypass_native_or_build_admission(self):
        from run_wj_h11 import intent
        for name in ('IncNavierStokesSolver','RANSWallResolutionProbe','cmake'):
            with self.assertRaises(AssertionError):intent(Path('/unused'),[name],'passive')
        self.assertEqual(intent(Path('/unused'),['python3','-B'],'passive'),(None,0))
    def test_restart_absolute_index_not_local_one(self):
        r,s=report();c,a,i=selected(r,4);self.assertEqual(i,2);self.assertEqual(c['step'],4)
        with self.assertRaises(ValueError):selected(r,1)
        with self.assertRaises(ValueError):selected(r,5)
    def test_selected_native_arithmetic_later_attempt(self):
        r,s=report();dp=np.zeros(104);dp[:2]=[.2,-.4];BG=np.zeros((104,2));BG[:2]=np.eye(2);S=np.zeros((104,104));S[:2,:2]=np.eye(2)
        m=dict(rows=[0,1],solution=[.2,-.4],rhs=[20.,-40.],linear_residual=[0.,0.]);p=dict(rows=list(range(4,108)),solution=dp.tolist(),rhs=dp.tolist(),linear_residual=[0.]*104)
        class Moments:
            def moments(self,x,solve=True):return {'4':dict(moments=[float(x[0]),float(x[1])]+[0.]*32),'5':dict(moments=[0.]*70)}
        self.assertEqual(arithmetic(r,s,m,p,np.eye(2),BG,S,100*np.eye(2),Moments(),step=4)['binding']['step'],4)
    def test_wrong_producer_alpha_and_history_reject(self):
        r,s=report();r['wj_audit_packets'][-1]['alpha_mass']=200.
        with self.assertRaises(ValueError):selected(r,4)
        r,s=report();r['wj_audit_packets'][-1]['stage_history_sha256']='bad'
        with self.assertRaises(ValueError):selected(r,4)
        r,s=report();r['accepted_states'][-1][0]=1
        with self.assertRaises(ValueError):selected(r,4)
    def test_checkpoint_suffix_has_no_step_authority(self):
        r,s=report();p=dict(steps=4,dt=.01,time=.04,signature='active',state=r['accepted_states'][-1],checkpoint='CP.1')
        self.assertTrue(checkpoint_match(p,r,4));p['steps']=1;self.assertFalse(checkpoint_match(p,r,4))
        p.update(steps=4,dt=.005);self.assertFalse(checkpoint_match(p,r,4))
    def test_exact_TIME_dt_and_additional_count(self):
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);x=ET.fromstring('<NEKTAR><GEOMETRY/><EXPANSIONS/><CONDITIONS><SOLVERINFO/><PARAMETERS/></CONDITIONS></NEKTAR>')
            for k,v in [('Time',0),('TimeStep',.01),('NumSteps',1),('IO_CheckSteps',1)]:set_parameter(x,k,v)
            old=root/'old.xml';write(x,old);set_parameter(x,'NumSteps',3);p=root/'new.xml';write(x,p);self.assertEqual(admit_xml(p,'coarse',old),3)
            for k,v in [('NumSteps',4),('Time',.01),('TimeStep',.005)]:
                y=ET.fromstring(p.read_text());set_parameter(y,k,v);bad=root/'bad.xml';write(y,bad)
                with self.assertRaises(AssertionError):admit_xml(bad,'coarse',old)
    def test_atomic_capacity_no_relaunch_or_borrowing(self):
        b=dict(reserved_attempts=0);native_reserve(b,'block:coarse',3)
        with self.assertRaises(AssertionError):native_reserve(b,'block:coarse',3)
        native_reserve(b,'block:fine',7);native_reserve(b,'block:quarter',15);self.assertEqual(b['reserved_attempts'],25)
        with self.assertRaises(AssertionError):native_reserve(b,'extra',3)
    def test_halt_merge_first_witness_and_earliest_time(self):
        with tempfile.TemporaryDirectory()as t:
            p=Path(t)/'journal.json';p.write_text(json.dumps(dict(entries={'op':dict(status='RESERVED',units=3)},reserved_attempts=3)))
            e=dict(status='OBSTRUCTED',role='accepted',trajectory='fine',time=.03)
            record_coverage_halt(p,e);record_coverage_halt(p,dict(e,time=.01));complete_operation(p,'op',dict(status='PASS',seconds=1))
            b=json.loads(p.read_text());self.assertEqual(b['coverage_halt']['time'],.03);self.assertEqual(b['coverage_earliest']['fine']['time'],.01);self.assertEqual(b['reserved_attempts'],3)
            with self.assertRaises(AssertionError):native_reserve(b,'block:fine',7)
    def test_common_integer_times_and_fine_off_common(self):
        self.assertEqual([common_step(dt,.03+2e-17)for dt in(.01,.005,.0025)],[3,6,12])
        self.assertEqual(common_step(.0025,.0175),7)
        with self.assertRaises(ValueError):common_step(.005,.0125)
    def test_floor_has_no_precise_zero_or_rate(self):
        self.assertEqual(ratio(0,0,4e-11)['status'],'UNRESOLVED')
        v=ratio(1e-6,1e-12,4e-11);self.assertEqual(v['status'],'INDISTINGUISHABLE');self.assertIsNone(v['rho']);self.assertAlmostEqual(v['rho_upper_bound'],4e-5)
        self.assertEqual(ratio(1e-6,2e-6,4e-11)['status'],'DOES_NOT_CONTRACT');self.assertEqual(ratio(1e-6,.5e-6,4e-11)['status'],'CONTRACTS')
        with self.assertRaises(ValueError):ratio(float('nan'),0,1e-11)

if __name__=='__main__':unittest.main()
