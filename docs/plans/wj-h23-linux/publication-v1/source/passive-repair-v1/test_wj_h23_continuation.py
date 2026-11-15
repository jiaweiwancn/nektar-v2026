# SPDX-License-Identifier: MIT
import copy,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from wj_h15_observer import identity,bits,vector_hash
from wj_h23_continuation import *
from run_wj_h23 import command_kind,reserve_native,usage,NATIVE
from wj_h20_analysis import continuation,ratio,FLOOR

def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return identity(path)

class H23Continuation(unittest.TestCase):
    def test_restage_paths_fail_before_admission_until_explicit_rebinding(self):
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);historical=put(root/'history.json',{})
            old=dict(active=dict(path='source-v2/numerics.py',sha256='identical',bytes=1))
            current=dict(active=dict(path='source-v3/numerics.py',sha256='identical',bytes=1))
            pins={k:put(root/(k+'.json'),{})for k in ('active_adapter','new_native_identity','runtime_dependencies')}
            c=dict(status='PASS',actual_imports_verified=True,historical_manifest=historical,active_sources=[old],**pins)
            pin=put(root/'old-compatibility.json',c)
            with patch('wj_h23_continuation.active_source_bindings',return_value=[current]):
                with self.assertRaisesRegex(AssertionError,'rebind before source freeze'):active_execution_bindings(pin)
                new_pin=put(root/'new-compatibility.json',dict(c,active_sources=[current]))
                self.assertEqual(active_execution_bindings(new_pin)['status'],'PASS')
            self.assertEqual(load(root/'old-compatibility.json')['active_sources'],[old])

    def report(self,z,n):
        state=[float(n),-0.];h=[2.,-0.];t=n*DT[z]
        return dict(steps=n,time=t,dt=DT[z],state=state,last_stage_history=h,accepted_states=[state],
            wj_audit_packets=[dict(state_sha256=vector_hash(state),stage_history_sha256=vector_hash(h),callback_history=h,time=t,alpha_mass=1/DT[z])])

    def test_exact_nonzero_matrix_and_common_time_map(self):
        self.assertEqual(sum(ADDITIONAL.values()),28);self.assertEqual(START,dict(C=4,M=8,F=16))
        self.assertEqual([[absolute_step(z,t)for z in ('C','M','F')]for t in TIMES],[[5,10,20],[6,12,24],[7,14,28],[8,16,32]])
        self.assertEqual([absolute_step(z,.02)for z in ('C','M','F')],[4,8,16])
        self.assertEqual([branch_forecast(z)for z in ORDER],[920,780,1400])
        self.assertLess(1.25*sum(branch_forecast(z)for z in ORDER),4200)
        for t in (.005,.01,.045):
            with self.assertRaises(AssertionError):absolute_step('F',t)

    def test_actual_integer_time_producer_and_genuine_startup_alpha(self):
        for z in DT:
            for n in (START[z],END[z]):
                r=self.report(z,n);self.assertEqual(selected(r,z,n)['step'],n)
                for key,value in (('time',r['time']+DT[z]),('alpha_mass',0.),('stage_history_sha256','wrong producer')):
                    bad=copy.deepcopy(r)
                    if key=='time':bad['wj_audit_packets'][0][key]=value
                    else:bad['wj_audit_packets'][0][key]=value
                    with self.assertRaises((AssertionError,ValueError)):selected(bad,z,n)
        self.assertEqual(absolute_step('F',np.nextafter(.035,np.inf)),28)

    def test_nonzero_F_intent_captures_off_and_no_callback_override(self):
        with tempfile.TemporaryDirectory()as t:
            pin=put(Path(t)/'completion.json',dict(binding=dict(absolute_step=16),callback=dict(path='authenticated')))
            v=dict(schema='WJMonitorIntentV1',h23_role='ContinuationFrom02V1',level='F',maximum_absolute_step=32,
                restored_absolute_step=16,recovery_stride=1,scientific_stride=4,dt_bits=bits(.00125),ack_timeout_seconds=120.,
                substep_evidence_steps=[],flow_capture='Off',snapshot_capture='Off',restart_completion=pin)
            validate_intent(v,16)
            for key,value in (('restored_absolute_step',0),('restored_absolute_step','16'),('maximum_absolute_step',16),
                ('flow_capture','MappedIncrementV1'),('substep_evidence_steps',[17]),('snapshot_capture','On'),
                ('restart_callback',dict(path='alternate')),('ack_timeout_seconds',121)):
                with self.assertRaises(AssertionError):validate_intent(dict(v,**{key:value}),16)
            with self.assertRaises(AssertionError):validate_intent(v,32)

    def test_partner_availability_is_phase_specific_and_missing_fails(self):
        self.assertEqual(partners('M'),());self.assertEqual(partners('C'),('M',));self.assertEqual(partners('F'),('C','M'))
        inputs=dict(anchors=dict(M='sealed M8'))
        with tempfile.TemporaryDirectory()as t:
            self.assertEqual(resolve(Path(t),inputs,'M',.02),'sealed M8')
            with self.assertRaisesRegex(AssertionError,'required computed H23 partner missing'):resolve(Path(t),inputs,'M',.025)
        self.assertEqual(continuation('PASS','NOT_YET_AVAILABLE'),'PROCEED')

    def test_governor_deadline_does_not_include_future_strict(self):
        remaining=local_forecast('M',14)
        self.assertEqual(remaining,185)
        self.assertTrue(deadline_ok(1.25*remaining+.01,0,'M',14))
        self.assertFalse(deadline_ok(1.25*remaining-.01,0,'M',14))
        self.assertLess(1.25*remaining,1.25*(remaining+120))

    def test_campaign_forecast_uses_elapsed_not_active_cap_twice(self):
        book=dict(entries={'old':dict(phase='passive',kind='passive',seconds=30,reserved_seconds=100),
            'active':dict(phase='native',kind='M',reserved_seconds=1050)})
        f=progress_forecast(book,'M',10,250)
        self.assertEqual(f['actual_native'],250);self.assertEqual(f['actual_total_with_closure'],640)
        self.assertEqual(f['protected_strict_seconds'],360);self.assertFalse(f['active_reservation_counted_as_spent'])
        self.assertEqual(f['status'],'PASS')

    def test_native_command_counts_order_no_retry_and_protected_stop_read(self):
        root=Path('/home/jay/build/wj-h23-flat-continuation-20261115')
        base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h23.py')),'--root',str(root),'--mode']
        for k in NATIVE:self.assertEqual(command_kind(root,base+[k]),k)
        for k in ('bootstrap','verify_M2','cp0','build','probe','M_retry'):
            with self.assertRaises(AssertionError):command_kind(root,base+[k])
        b=dict(entries={},native_reservations=[],reserved_attempts=0)
        reserve_native(b,'M')
        with self.assertRaises(AssertionError):reserve_native(b,'M')
        b['coverage_halt']={'first':'retained'};reserve_native(b,'strict_M')
        with self.assertRaises(AssertionError):reserve_native(b,'C')
        b.pop('coverage_halt');reserve_native(b,'C');reserve_native(b,'strict_C');reserve_native(b,'F');reserve_native(b,'strict_F')
        self.assertEqual(b['reserved_attempts'],28)

    def test_startup_strict_keeps_original_callback_and_no_terminal_fallback(self):
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);directory=root/'monitor';directory.mkdir()
            a={k:put(root/(k+'.json'),dict(role=k))for k in ('report','callback','checkpoint','completion')}
            q=dict(schema='WJStartupVerificationCommitV1',binding=dict(absolute_step='16'),completion=a['completion'],checkpoint=a['checkpoint'])
            put(directory/'startup.verification-commit.json',q);put(directory/'latest-fully-verified-checkpoint.json',q)
            r=dict(time=.02,last_stage_alpha=800.)
            with patch('wj_h23_continuation.qualified',return_value=(r,{}))as guard:
                request,selection=strict_request(directory,a,'F',False)
                self.assertEqual(request['native_report'],a['report']['path']);self.assertEqual(request['completed_sidecar'],a['callback'])
                self.assertEqual(selection['step'],16);self.assertEqual(request['import_role'],'QualifiedV1');guard.assert_called_once()
            with self.assertRaisesRegex(AssertionError,'no fallback'):strict_request(directory,a,'F',True)
            q['schema']='PREPARED_NOT_NATIVE_COMMIT';put(directory/'latest-fully-verified-checkpoint.json',q)
            with self.assertRaises(AssertionError):strict_request(directory,a,'F',False)

    def test_noncontraction_is_descriptive_but_absolute_allowance_stops(self):
        self.assertEqual(continuation('PASS','PASS',ratio(1e-6,2e-6)['status']),'PROCEED')
        self.assertEqual(continuation('PASS','PASS',ratio(FLOOR,1e-6)['status']),'PROCEED')
        self.assertEqual(continuation('PASS','UNRESOLVED'),'REVIEW_STOP')
        self.assertEqual(continuation('PASS','EXCEEDED'),'REVIEW_STOP')
        self.assertEqual(continuation('OBSTRUCTED','PASS'),'PHYSICAL_STOP')

    def test_completion_integrity_rejects_signature_history_and_callback_corruption(self):
        # Exercise the actual completed-reader checks with a tiny serialization fixture.
        import hashlib
        z='F';n=16;r=self.report(z,n)
        r.update(signature='own-F',run_configuration_sha256='unused',fixed_reference=[2.],fixed_defect=[0.],original_qh0_state_sha256='qh0')
        with tempfile.TemporaryDirectory()as t:
            root=Path(t);configuration=put(root/'config.json',{})
            r['run_configuration_sha256']=configuration['sha256'];report=put(root/'report.json',r)
            side=put(root/'side.json',dict(schema='WJCompletedCallbackV1',state_sha256=vector_hash(r['state']),
                callback_history_sha256=vector_hash(r['last_stage_history']),callback_history=r['last_stage_history'],time=.02,alpha_mass=800.))
            a=dict(report=report,checkpoint=put(root/'State.bin',{}),callback=side)
            b=dict(absolute_step='16',physical_signature_sha256=hashlib.sha256(b'own-F').hexdigest(),checkpoint_sha256=a['checkpoint']['sha256'],
                report_sha256=report['sha256'],callback_sidecar_sha256=side['sha256'],state_sha256=vector_hash(r['state']),
                callback_history_sha256=vector_hash(r['last_stage_history']),time_bits=bits(.02),alpha_bits=bits(800.),
                fixed_reference_sha256=vector_hash(r['fixed_reference']),fixed_defect_sha256=vector_hash(r['fixed_defect']),original_qh0_sha256='qh0')
            completion=dict(schema='WJNativeCompletionV1',binding=b,configuration=configuration,**a)
            a['completion']=put(root/'completion.json',completion)
            with patch('wj_h23_continuation.checkpoint_payload',return_value={}),patch('wj_h23_continuation.same_cp'):
                read_completed(a,z,n)
                for k,v in (('alpha_bits',bits(0.)),('physical_signature_sha256','wrong'),('callback_history_sha256','wrong'),('original_qh0_sha256','wrong')):
                    bad=copy.deepcopy(completion);bad['binding'][k]=v;pin=put(root/'bad-completion.json',bad)
                    with self.assertRaises(AssertionError):read_completed(dict(a,completion=pin),z,n)

if __name__=='__main__':unittest.main()
