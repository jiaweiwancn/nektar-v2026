# SPDX-License-Identifier: MIT
import copy,sys,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from wj_h15_observer import vector_hash,atomic,identity
from wj_h20_analysis import DT,ADDITIONAL,START,END,STRIDE,TIMES,absolute_step,selected,availability,ratio,continuation,all_pair_screens,shared_direction_actions,FLOOR
from wj_h20_observer import resolve
from run_wj_h20 import native_selection,reserve_strict,usage
from study_wj_h20 import validate_monitor_intent,ACK_TIMEOUT,admission
from wj_h15_observer import bits

class ThreeLevel(unittest.TestCase):
    def report(self,level,step):
        state=[float(step),-0.];h=[1.,0.];dt=DT[level]
        return dict(steps=step,dt=dt,accepted_states=[state],wj_audit_packets=[dict(state_sha256=vector_hash(state),stage_history_sha256=vector_hash(h),callback_history=h,time=step*dt,alpha_mass=1/dt)])

    def test_exact_map_order_attempts_and_CP0_command(self):
        self.assertEqual(list(ADDITIONAL),['M','C','F']);self.assertEqual(sum(ADDITIONAL.values()),22)
        self.assertEqual([[absolute_step(z,t)for z in ('C','M','F')]for t in TIMES],[[1,2,4],[2,4,8],[3,6,12],[4,8,16]])
        self.assertEqual({z:START[z]+ADDITIONAL[z]for z in DT},END)
        self.assertEqual([DT[z]*STRIDE[z]for z in DT],[.005]*3)
        self.assertEqual(absolute_step('M',np.nextafter(.015,np.inf)),6)
        root=Path('/home/jay/build/wj-h20-three-level-20261115');base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h20.py')),'--root',str(root),'--mode']
        self.assertEqual(native_selection(root,base+['bootstrap'],'bootstrap',0),'F')
        for z,n in ADDITIONAL.items():
            self.assertEqual(native_selection(root,base+['target','--level',z],'target',n),z)
            with self.assertRaises(AssertionError):native_selection(root,base+['target','--level',z],'target',n+1)
        with self.assertRaises(AssertionError):native_selection(root,base+['target','--level','F','--dt','.01'],'target',16)

    def test_wrong_integer_time_producer_and_alpha_rejected(self):
        for z,n in (('C',3),('M',6),('F',12)):
            r=self.report(z,n);self.assertEqual(selected(r,z,.015)['step'],n)
            bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['time']=.02
            with self.assertRaises(AssertionError):selected(bad,z,.015)
            bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['stage_history_sha256']='wrong producer'
            with self.assertRaises(ValueError):selected(bad,z,.015)
            bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['alpha_mass']=0
            with self.assertRaises(AssertionError):selected(bad,z,.015)

    def test_restored_alpha_zero_is_not_original_completed_callback(self):
        r=self.report('M',4);r['wj_audit_packets'][0]['alpha_mass']=0
        with self.assertRaises(AssertionError):selected(r,'M',.01)
        self.assertEqual(selected(self.report('M',4),'M',.01)['alpha'],400)

    def test_planned_absence_is_not_PASS_or_future_wait(self):
        self.assertEqual(availability('M'),'NOT_YET_AVAILABLE')
        self.assertEqual(continuation('PASS','NOT_YET_AVAILABLE'),'PROCEED')
        with tempfile.TemporaryDirectory()as t:
            with self.assertRaisesRegex(AssertionError,'required computed partner absent'):resolve(Path(t),dict(states={}),'M',.015)

    def test_absolute_stop_is_separate_from_contraction_and_physics(self):
        for diagnosis in ('NONCONTRACTING','UNRESOLVED'):
            self.assertEqual(continuation('PASS','PASS',diagnosis),'PROCEED')
        self.assertEqual(continuation('PASS','EXCEEDED'),'REVIEW_STOP')
        self.assertEqual(continuation('PASS','UNRESOLVED'),'REVIEW_STOP')
        self.assertEqual(continuation('OBSTRUCTED','PASS'),'PHYSICAL_STOP')

    def test_ratio_floor_keeps_signed_data_and_no_false_zero_order(self):
        self.assertEqual(ratio(FLOOR,1e-12)['status'],'UNRESOLVED')
        r=ratio(1e-6,1e-12);self.assertIsNone(r['rho']);self.assertIsNone(r['p_obs']);self.assertAlmostEqual(r['upper_bound'],FLOOR/1e-6)
        self.assertEqual(ratio(1e-6,2e-6)['status'],'NONCONTRACTING')
        self.assertAlmostEqual(ratio(1e-6,5e-7)['p_obs'],1.)

    def test_direct_CF_exceedance_is_not_hidden_by_adjacent_PASS(self):
        data={z:dict(value=v)for z,v in (('C',0.),('M',.75),('F',1.5))}
        seen=[]
        def metric(rows,a,b):
            d=abs(a['value']-b['value']);seen.append(d)
            return {'full':{'U':{k:dict(status='EXCEEDED'if d>1 else'PASS')for k in ('rms','maximum')}}}
        def principal(a,b,outer):return dict(status='PASS'),{}
        with patch('wj_h20_analysis.metrics',metric),patch('wj_h20_analysis.principal_screen',principal):
            r,_=all_pair_screens(SimpleNamespace(outer=np.array([True])),data)
        self.assertEqual(seen,[.75,.75,1.5]);self.assertEqual(r['status'],'EXCEEDED');self.assertEqual(r['pairs']['C-M']['status'],'PASS');self.assertEqual(r['pairs']['C-F']['status'],'EXCEEDED')

    def test_three_common_directions_evaluate_all_three_actual_states(self):
        data={z:dict(k=np.array([1.+i]),epsilon=np.ones(1),G=np.full((1,3,3),i),angle=np.array([.1+.2*i]))for i,z in enumerate(('C','M','F'))};calls=[]
        def action(k,e,G,angle):calls.append((float(k[0]),float(angle[0])));return dict(normalized=k+angle)
        with patch('wj_h20_analysis.fixed_action',action):a=shared_direction_actions(data,np.array([0]))
        self.assertEqual(len(calls),9)
        for z in ('C','M','F'):
            for level in ('C','M','F'):self.assertAlmostEqual(a[z+level]['normalized'][0],data[level]['k'][0]+data[z]['angle'][0])

    def test_final_strict_scope_is_unique_and_never_fourth(self):
        with tempfile.TemporaryDirectory()as t:
            p=Path(t);atomic(p/'accounting.json',dict(closed=False))
            for z in ('M','C','F'):reserve_strict(p,z)
            with self.assertRaises(AssertionError):reserve_strict(p,'M')
            with self.assertRaises(AssertionError):reserve_strict(p,'extra')

    def test_H19_commit_before_review_stop_and_H18_default_remain(self):
        source=(Path(__file__).parents[2]/'EquationSystems/RANSEquilibriumWallEvolution.cpp').read_text()
        at=source.index('const auto qualified=J::Object')
        self.assertLess(source.index('M::Atomic(prefix+".qualification-commit.json"',at),source.index('Require(next!="REVIEW_STOP"',at))
        self.assertIn('CheckpointAuditV1',source)

    def test_native_monitor_timeout_preflight_rejects_failed_180s_intent(self):
        source=(Path(__file__).parents[2]/'EquationSystems/RANSEquilibriumWallEvolution.cpp').read_text()
        self.assertIn('get<double>("ack_timeout_seconds")<=120',source)
        for level in DT:
            intent=dict(schema='WJMonitorIntentV1',nonce='a'*32,level=level,maximum_absolute_step=END[level],restored_absolute_step=START[level],recovery_stride=1,scientific_stride=STRIDE[level],dt_bits=bits(DT[level]),ack_timeout_seconds=ACK_TIMEOUT)
            validate_monitor_intent(intent)
            for timeout in (0,180):
                with self.assertRaisesRegex(AssertionError,'native bounded ACK timeout'):validate_monitor_intent(dict(intent,ack_timeout_seconds=timeout))
            with self.assertRaises(AssertionError):validate_monitor_intent(dict(intent,maximum_absolute_step=END[level]+1))

    def test_first_native_halt_keeps_later_branches_closed(self):
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);atomic(root/'accounting.json',dict(target_halt={'operation':'original failed M launch'}))
            with patch('study_wj_h20.frozen',lambda root:None):
                for level in ('C','F'):
                    with self.assertRaises(AssertionError):admission(root,level)

    def test_resumption_carries_all_prior_work_without_reclaiming_charges(self):
        book=dict(carryforward_seconds=714.0562545959547,carryforward_native_seconds=42.09456458600471,entries={'failed':dict(phase='target',reserved_seconds=10,seconds=3.5)})
        before=usage(book)
        self.assertAlmostEqual(before['total'],714.0562545959547+60+300+3.5)
        self.assertAlmostEqual(before['native'],42.09456458600471+3.5)
        book['entries']['unknown']=dict(phase='target',reserved_seconds=30)
        self.assertEqual(usage(book)['total']-before['total'],30)
        self.assertEqual(usage(book)['native']-before['native'],30)

if __name__=='__main__':unittest.main()
