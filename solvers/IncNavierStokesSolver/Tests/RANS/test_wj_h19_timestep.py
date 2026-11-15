# SPDX-License-Identifier: MIT
import copy,sys,unittest
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import numpy as np
from wj_h19_analysis import physical_pair,spatial_xml,fixed_action,principal_screen,classification,metrics,GATE,FLOOR
from wj_h7_analysis import principal_batch
from wj_h15_observer import vector_hash
from run_wj_h19 import native_selection

class TimestepSensitivity(unittest.TestCase):
    def report(self,dt,step):
        state=[float(step),-0.];h=[1.,0.]
        return dict(steps=step,dt=dt,accepted_states=[state],wj_audit_packets=[dict(state_sha256=vector_hash(state),stage_history_sha256=vector_hash(h),
            callback_history=h,time=step*dt,alpha_mass=1/dt)])

    def test_physical_pairs_use_fine_two_four_not_equal_step_or_suffix(self):
        for n in (1,2):
            a,b=physical_pair(self.report(.005,n),self.report(.0025,2*n),n)
            self.assertEqual((a['step'],b['step']),(n,2*n));self.assertNotEqual(a['alpha'],b['alpha'])
            with self.assertRaises((AssertionError,ValueError,IndexError)):physical_pair(self.report(.005,n),self.report(.0025,n),n)
        bad=self.report(.0025,2);bad['wj_audit_packets'][0]['time']=.01
        with self.assertRaises(AssertionError):physical_pair(self.report(.005,1),bad,1)
        bad=self.report(.0025,2);bad['wj_audit_packets'][0]['stage_history_sha256']='wrong producer'
        with self.assertRaises(ValueError):physical_pair(self.report(.005,1),bad,1)

    def test_xml_dt_is_the_only_new_parameter_allowance(self):
        x=ET.fromstring('<NEKTAR><CONDITIONS><PARAMETERS><P>Time = 0</P><P>TimeStep = .0025</P><P>Kinvis = .0001</P><P>NumSteps = 0</P></PARAMETERS><SOLVERINFO/><FUNCTION NAME="Initial"><E VALUE="log(y)"/></FUNCTION></CONDITIONS></NEKTAR>')
        y=copy.deepcopy(x);y.find('./CONDITIONS/PARAMETERS/P[2]').text='TimeStep = .005'
        self.assertEqual(spatial_xml(x),spatial_xml(y))
        for path,text in (('./CONDITIONS/PARAMETERS/P[1]','Time = .005'),('./CONDITIONS/PARAMETERS/P[3]','Kinvis = .1')):
            z=copy.deepcopy(y);z.find(path).text=text;self.assertNotEqual(spatial_xml(x),spatial_xml(z))
        y.find('./CONDITIONS/FUNCTION/E').set('VALUE','sqrt(y)');self.assertNotEqual(spatial_xml(x),spatial_xml(y))

    def test_fixed_directions_reproduce_minimum_and_detect_signed_change(self):
        k=np.array([.0088]);e=np.array([.002]);G=np.zeros((1,3,3));G[0,0,2]=.6;G[0,2,0]=.001
        R,nt,N,m,t,angle=principal_batch(k,e,G)
        a=fixed_action(k,e,G,angle)
        self.assertLess(abs(a['normalized'][0]-m[0]),2e-11)
        self.assertLess(abs(np.sum(a['a']*a['xi'])),1e-15)
        d=dict(k=k,epsilon=e,G=G,margin=m,angle=angle)
        screen,arrays=principal_screen(d,d,np.array([True]));self.assertEqual(screen['status'],'PASS')
        self.assertTrue(np.all(arrays['minimized_signed_difference']==0))
        changed=dict(d,margin=m+1.)
        self.assertEqual(principal_screen(changed,d,np.array([True]))[0]['status'],'EXCEEDED')

    def test_principal_budget_floor_is_unresolved_not_positive_allowance(self):
        k=np.array([.0088]);e=np.array([.002]);G=np.zeros((1,3,3));G[0,0,2]=.6
        d=dict(k=k,epsilon=e,G=G,margin=np.array([GATE+FLOOR]),angle=np.array([.3]))
        screen,_=principal_screen(d,d,np.array([True]));self.assertEqual(screen['status'],'UNRESOLVED')
        self.assertIsNone(screen['checks']['minimized']['max_utilization'])
        self.assertFalse(classification(1e-12,1e-3)['resolved']);self.assertEqual(classification(2e-3,1e-3)['status'],'EXCEEDED')

    def test_rms_is_vector_scaled_and_height_relative_epsilon_is_not_hidden(self):
        exact=dict(U=np.zeros((2,3)),G=np.zeros((2,3,3)),k=np.ones(2),R=np.zeros((2,3,3)),pi=np.zeros(2),epsilon=np.array([0.,1e-5]))
        a=copy.deepcopy(exact);a['U'][1]=[3e-5,0.,4e-5];a['epsilon'][1]*=1.02
        r=SimpleNamespace(rows=[None,None],outer=np.array([False,True]),weights=np.array([3.75,494.0625]),exact=exact,config=SimpleNamespace(K=1.,u_tau=.051))
        v=metrics(r,a,exact)
        self.assertAlmostEqual(v['outer']['U']['maximum']['value'],5e-5)
        self.assertAlmostEqual(v['outer']['epsilon']['rms']['value'],.02)
        self.assertEqual(v['outer']['epsilon']['rms']['status'],'EXCEEDED');self.assertNotIn('epsilon',v['inner'])

    def test_one_bootstrap_and_two_attempt_native_command_ceiling(self):
        root=Path('/home/jay/build/wj-h19-target-timestep-20261115')
        script=str(Path(__file__).with_name('study_wj_h19.py'))
        for kind,n in (('bootstrap',0),('target',2),('component',0)):
            argv=[sys.executable,'-B',script,'--root',str(root),'--mode',kind];native_selection(root,argv,kind,n)
            with self.assertRaises(AssertionError):native_selection(root,argv,kind,n+1)
            with self.assertRaises(AssertionError):native_selection(root,argv+['--dt','.01'],kind,n)

    def test_review_stop_occurs_after_commit_and_before_next_advance(self):
        source=(Path(__file__).parents[2]/'EquationSystems/RANSEquilibriumWallEvolution.cpp').read_text()
        start=source.index('const auto qualified=J::Object')
        self.assertLess(source.index('M::Atomic(prefix+".qualification-commit.json"',start),source.index('Require(next!="REVIEW_STOP"',start))
        self.assertIn('TEMPORAL_REVIEW_STOP',source)

if __name__=='__main__':unittest.main()
