#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
from pathlib import Path
import unittest
import numpy as np
from wj_h12_phase import (ATTEMPTS,WINDOWS,BRANCHES,OWNER_KEYS,PHASES,SparseRows,checked_context,bind,
                         local,decompositions,window,paired_window,cancellation,identity,kernel_input,
                         predictor_eligible,scatter,validate_owner,DIRECTION,role_context)
from wj_h4_points import chart,owner_point,vector_hash
from wj_h9_local import same
from test_wj_h11_continuation import report

def fixture():
    r=dict(parameters=[6.423173803526447e-5,.15,.1085450354253546],
           geometry_signature='fixture:CompositeTunnelV1:0.0875:20',offsets=[i*48 for i in range(7)])
    faces=[.2375,OWNER_KEYS[0][2],OWNER_KEYS[2][2],.4];mesh=[]
    for j,(a,b)in enumerate(zip(faces[:-1],faces[1:])):
        lo=chart(r,[0,0,a]).tolist();hi=chart(r,[8,0,b]).tolist()
        mesh.append(dict(geometry=17+j,lower=lo,upper=hi,degree=3,ids=list(range(j*16,(j+1)*16)),sign=[-1 if i%3==0 else 1 for i in range(16)]))
    r['fields']=[dict(mesh=copy.deepcopy(mesh))for _ in range(6)];state=np.zeros(288)
    for b in range(6):
        for el,a,z in zip(mesh,faces[:-1],faces[1:]):
            # Modified_A endpoint modes exactly represent a physical affine field.
            c=np.zeros((4,4));c[0,:2]=[b+3*a,b+16+3*a];c[1,:2]=[b+3*z,b+16+3*z]
            state[np.array(el['ids'])+r['offsets'][b]]=c.ravel()*el['sign']
    owners=[]
    for geo,side,z in OWNER_KEYS:
        p=[8.,0.,z];el=next(e for e in mesh if e['geometry']==geo)
        owners.append(dict(physical=p,owner=owner_point(r,el,p)))
    return r,state,owners

class SavedPhaseDiagnosis(unittest.TestCase):
    def test_exact_fourteen_attempt_equal_windows(self):
        self.assertEqual(len(ATTEMPTS),14);self.assertEqual(len(ATTEMPTS)*4*4,224)
        for w,(lo,hi)in enumerate(WINDOWS):
            for b,dt in BRANCHES:
                steps=[s for ww,bb,d,s in ATTEMPTS if (ww,bb)==(w,b)]
                self.assertAlmostEqual(len(steps)*dt,hi-lo)
                self.assertAlmostEqual((steps[0]-1)*dt,lo);self.assertAlmostEqual(steps[-1]*dt,hi)
        self.assertEqual([s for w,b,dt,s in ATTEMPTS if w==1 and b=='quarter'],[13,14,15,16])
    def test_wrong_integer_time_and_wrong_producer(self):
        r,s=report();r['wj_audit_packets'][-1]['time']=.03
        with self.assertRaises(ValueError):checked_context(r,4)
        with self.assertRaises(ValueError):bind(r,s,4,{'sha256':'wrong'},{'sha256':'right'})
        # A predecessor at a wrong but internally aligned time must also fail.
        r,s=report();r['wj_audit_packets'][-2]['time']=.02
        with self.assertRaises(ValueError):checked_context(r,3)
    def test_sparse_owner_signs_and_physical_gradient(self):
        r,s,o=fixture();cache=SparseRows(r,o);v=cache.linear(s)
        self.assertLess(cache.reference_check(s),1e-13)
        for row,rec in zip(v,o):
            np.testing.assert_allclose(row['U'],[16+3*rec['physical'][2],0,17+3*rec['physical'][2]],atol=1e-13)
            np.testing.assert_allclose(row['G'][0],[2,0,3],atol=1e-13)
            np.testing.assert_allclose(row['G'][2],[2,0,3],atol=1e-13)
        bad=copy.deepcopy(o);bad[0]['owner']['geometry']=18
        with self.assertRaises(ValueError):SparseRows(r,bad)
        bad=copy.deepcopy(o);bad[0]['owner']['roundoff_snap'][2]+=1e-10
        with self.assertRaises(ValueError):validate_owner(r,bad[0],0)
    def test_masks_aliases_and_actual_stale_epsilon(self):
        r,s=report();r['fixed_defect']=[0.]*108
        a=r['pimple_attempts'][2];a.update(scalar_order='k_then_epsilon',history_sha256=vector_hash(s['history']),incoming_sha256=vector_hash(s['committed']))
        states,v=bind(r,s,4,{'sha256':'a'},{'sha256':'a'})
        self.assertTrue(v['P_C_bitwise']);self.assertTrue(same(states['K'][3:4],states['F'][3:4]))
        self.assertFalse(same(states['K'][3:4],states['A'][3:4]))
        bad=copy.deepcopy(s);bad['raw'][0]=1;bad['corrected'][0]=1
        with self.assertRaises(AssertionError):bind(r,bad,4,{}, {})
    def test_predecessor_role_keeps_its_own_time_and_step(self):
        b=dict(time=.04,step=4,predecessor_time=.03,predecessor_step=3,stage_history_sha256='saved',alpha=100.)
        p=role_context(b,'P');a=role_context(b,'A');c=role_context(b,'C')
        self.assertEqual((p['time'],p['state_absolute_step']),(.03,3))
        self.assertEqual((a['time'],a['state_absolute_step']),(.04,4))
        self.assertFalse(c['accepted_role']);self.assertIsNone(c['state_absolute_step'])
        self.assertEqual(p['attempt_time'],a['attempt_time'])
    def test_fixed_action_normalization_gap_and_B_identities(self):
        G=np.array([[.003,0,.8],[0,0,0],[.002,0,-.001]])
        a=local(dict(G=G,k=.009,epsilon=.002));b=local(dict(G=G*1.002,k=.0091,epsilon=.00203))
        self.assertEqual(a['constitutive_status'],'DEFINED');d=decompositions(a,b)
        self.assertGreater(np.linalg.norm(d['delta_B']['gradient_term']),0)
        self.assertGreater(np.linalg.norm(d['delta_B']['timescale_term']),0)
        self.assertLess(d['optimizer_identity']['maximum_error'],1e-15)
        self.assertGreater(abs(d['normalization_term']),0)
        self.assertGreater(abs(d['dimensional_term']),0)
        self.assertAlmostEqual(np.linalg.norm(DIRECTION),1.)
        bad=local(dict(G=G,k=-1.,epsilon=.002));self.assertEqual(decompositions(a,bad)['status'],'NOT_EVALUABLE')
    def test_equal_window_nonzero_incoming_and_signed_cancellation(self):
        def phases(start,changes):
            vals=np.r_[start,start+np.cumsum(changes)]
            return {p:dict(mhat=x)for p,x in zip(PHASES,vals)}
        a=window([phases(2.,[0.,-3.,2.,1.2])],'mhat')
        b=window([phases(1.,[0.,.5,-.2,.1]),phases(1.4,[0.,.1,-.1,.3])],'mhat')
        p=paired_window(a,b,'mhat');self.assertEqual(float(p['incoming']),1.)
        self.assertAlmostEqual(float(p['endpoint']),.5)
        self.assertGreater(float(p['cancellation']['cancellation']),.8)
        c=cancellation(0.,[1e-12,-2e-12],-1e-12,1.)
        self.assertFalse(c['resolved']);self.assertTrue(np.isnan(c['cancellation']))
        with self.assertRaises(ValueError):identity(1e-10,[1e-10])
    def test_exact_signed_aliases_and_missing_predictor_scope(self):
        v=dict(k=.009,epsilon=.002,G=np.zeros((3,3)));u=copy.deepcopy(v);u['G'][0,0]=-0.
        self.assertNotEqual(kernel_input(v),kernel_input(u))
        self.assertEqual(sum(predictor_eligible(b,s)for w,b,dt,s in ATTEMPTS),3)
        self.assertFalse(any(predictor_eligible(b,s)for w,b,dt,s in ATTEMPTS if w==1))
        with self.assertRaises(ValueError):scatter([0,0],[1,2],4)
    def test_reject_all_native_strict_mass_and_build_admission(self):
        from run_wj_h12 import passive_intent
        for command in ('IncNavierStokesSolver','RANSWallResolutionProbe','WJKernelProbe','cmake','ninja'):
            with self.assertRaises(AssertionError):passive_intent([command])
        for flag in ('--governing','--strict','--cp0'):
            with self.assertRaises(AssertionError):passive_intent(['python3','study_wj_h12.py',flag])
        self.assertIsNone(passive_intent(['python3','-B','study_wj_h12.py','--mode','science']))

if __name__=='__main__':unittest.main()
