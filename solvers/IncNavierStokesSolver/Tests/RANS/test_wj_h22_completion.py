# SPDX-License-Identifier: MIT
import copy,json,sys,tempfile,unittest
from pathlib import Path
from unittest.mock import patch,Mock
from wj_h15_observer import identity,bits,vector_hash
from wj_h22_completion import validate_role,typed_strict,restored_parity,require_M2_gate
from run_wj_h22_completion import command_kind,reserve_native,usage,NATIVE
from wj_h22_replay import DT,flow_vectors,phases
from test_wj_h22_phase import actual_fixture,vector_fixture

def put(path,value):
    path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value));return identity(path)

def records():
    r=dict(schema='WJCompletedStepV1',state=[1.,-0.],fixed_reference=[2.],fixed_defect=[0.],dt=.0025,time=.005,steps=2,
        signature='own-M',native_history=[[1.,-0.]],native_history_times=[.005],last_stage_history=[3.,-0.],
        last_stage_alpha=400.,last_stage_time=.005,original_qh0_state_sha256='qh0',native_full_audit_status='NOT_SCHEDULED',native_hard_audit={})
    audited=copy.deepcopy(r);audited.update(native_full_audit_status='PASS',native_hard_audit=dict(admissible=[1.],margin=[.004]),last_stage_alpha=0.,last_stage_history=[9.])
    current=copy.deepcopy(r);current.update(schema='WJMonitoringStartupV1',native_full_audit_status='PASS',native_hard_audit=copy.deepcopy(audited['native_hard_audit']))
    return r,audited,current,dict(eta=0.)

def startup_fixture(root):
    r,_,_,cp=records();directory=root/'monitor';directory.mkdir()
    anchor={k:put(root/(k+'.json'),dict(role=k))for k in ('report','callback','checkpoint')}
    b=dict(absolute_step=2,state_sha256=vector_hash(r['state']),time_bits=bits(.005),alpha_bits=bits(400.),callback_history_sha256=vector_hash(r['last_stage_history']))
    anchor['completion']=put(root/'completion.json',dict(schema='WJNativeCompletionV1',binding=b,**anchor))
    result=put(directory/'startup.result.json',dict(status='PASS',binding=b))
    verified=put(directory/'startup.verified.json',dict(status='PASS',binding=b,result=result))
    q=dict(schema='WJStartupVerificationCommitV1',binding=b,checkpoint=anchor['checkpoint'],completion=anchor['completion'],result=verified)
    put(directory/'startup.verification-commit.json',q);put(directory/'latest-fully-verified-checkpoint.json',q)
    return directory,anchor,r,cp,q

class CompletionRoles(unittest.TestCase):
    def test_only_zero_M_verification_four_F_advances_and_two_reads(self):
        root=Path('/home/jay/build/wj-h22-completion-20261115')
        base=[sys.executable,'-B',str(Path(__file__).with_name('study_wj_h22_completion.py')),'--root',str(root),'--mode']
        for kind in NATIVE:self.assertEqual(command_kind(root,base+[kind]),kind)
        for kind in ('C','M','bootstrap','cp0','build','strict_C','F_retry'):
            with self.assertRaises(AssertionError):command_kind(root,base+[kind])
        b=dict(carry=dict(total=1169.9461947359669,native=465.6705588850309,build=135.61016353601008,analysis=11.737643021973781),
               entries={},new_native_reservations=[],new_attempt_capacity=0)
        self.assertLess(usage(b)['native']+1.25*sum(NATIVE.values()),1800)
        with self.assertRaises(AssertionError):reserve_native(b,'F')
        reserve_native(b,'verify_M2')
        with self.assertRaises(AssertionError):reserve_native(b,'verify_M2')
        reserve_native(b,'strict_M2');reserve_native(b,'F');self.assertEqual(b['new_attempt_capacity'],4)
        reserve_native(b,'strict_F')
        with self.assertRaises(AssertionError):reserve_native(b,'F')
        bad=copy.deepcopy(b);bad['new_native_reservations']=[];bad['integrity_halt']={'first':'retained'}
        with self.assertRaises(AssertionError):reserve_native(bad,'verify_M2')

    def test_retained_role_cannot_advance_reset_or_override_callback(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);completion=put(root/'completion.json',dict(binding=dict(absolute_step=2),callback=dict(path='authenticated')))
            v=dict(schema='WJMonitorIntentV1',completion_role='RetainedM2StartupV1',level='M',maximum_absolute_step=2,restored_absolute_step=2,
                scientific_stride=2,recovery_stride=1,dt_bits=bits(.0025),substep_evidence_steps=[],flow_capture='Off',ack_timeout_seconds=120.,restart_completion=completion)
            validate_role(v,0)
            with self.assertRaises(AssertionError):validate_role(v,1)
            for key,value in (('maximum_absolute_step',4),('restored_absolute_step',0),('dt_bits',bits(.00125)),
                ('maximum_absolute_step',2.),('restored_absolute_step','2'),
                ('flow_capture','MappedIncrementV1'),('substep_evidence_steps',[2]),('ack_timeout_seconds',0),('ack_timeout_seconds',121),
                ('restart_callback',dict(path='substituted'))):
                with self.assertRaises(AssertionError):validate_role(dict(v,**{key:value}),0)
            put(root/'cp0-completion.json',dict(binding=dict(absolute_step=0),callback=dict(path='cp0')))
            with self.assertRaises(AssertionError):validate_role(dict(v,restart_completion=identity(root/'cp0-completion.json')),0)

    def test_F_role_requires_own_initial_history_and_exact_capture_selection(self):
        v=dict(schema='WJMonitorIntentV1',completion_role='MissingFReplayV1',level='F',maximum_absolute_step=4,restored_absolute_step=0,
            scientific_stride=4,recovery_stride=1,dt_bits=bits(.00125),substep_evidence_steps=[1,2,3,4],flow_capture='MappedIncrementV1',ack_timeout_seconds=120.)
        validate_role(v,4)
        for key,value in (('maximum_absolute_step',5),('restored_absolute_step',2),('substep_evidence_steps',[4]),('restart_callback',{}),('level','M')):
            with self.assertRaises(AssertionError):validate_role(dict(v,**{key:value}),4)

    def test_restored_producer_history_parity_precedes_auxiliary_map(self):
        r,a,current,cp=records()
        with tempfile.TemporaryDirectory()as tmp:
            audit=put(Path(tmp)/'B-audit.json',a);descriptor=dict(completion='actual-completion',callback='actual-callback')
            with patch('wj_h22_completion.read_completed',side_effect=[(r,cp),(r,cp)]):
                result=restored_parity(current,descriptor,{},dict(audit_report=audit),cp)
            self.assertEqual(result['status'],'PASS');self.assertFalse(result['startup_is_producer'])
            self.assertFalse(result['full_map_parity']['history_or_alpha_replaced'])
            bad=copy.deepcopy(current);bad['last_stage_history']=[9.]
            with patch('wj_h22_completion.read_completed',side_effect=[(r,cp),(r,cp)]),patch('wj_h22_completion.native_full_map_parity')as maps:
                with self.assertRaises(AssertionError):restored_parity(bad,descriptor,{},dict(audit_report=audit),cp)
                maps.assert_not_called()
            bad=copy.deepcopy(current);bad['schema']='WJCompletedStepV1'
            with self.assertRaises(AssertionError):restored_parity(bad,descriptor,{},dict(audit_report=audit),cp)

    def test_strict_startup_commit_keeps_original_completed_producer(self):
        with tempfile.TemporaryDirectory()as tmp:
            directory,anchor,r,cp,q=startup_fixture(Path(tmp))
            with patch('wj_h22_completion.read_completed',return_value=(r,cp)):
                request=typed_strict(directory,anchor,'M',2)
            self.assertEqual(request['import_role'],'QualifiedV1');self.assertEqual(request['native_report'],anchor['report']['path'])
            self.assertEqual(request['completed_sidecar'],anchor['callback']);self.assertIn('startup.verification-commit',request['qualified_cursor']['path'])
            for schema,step in (('PREPARED_NOT_NATIVE_COMMIT',2),('WJStartupVerificationCommitV1',0),('PendingNativeV1',2)):
                bad=copy.deepcopy(q);bad['schema']=schema;bad['binding']['absolute_step']=step
                put(directory/'latest-fully-verified-checkpoint.json',bad)
                with self.assertRaises(AssertionError):typed_strict(directory,anchor,'M',2)
            put(directory/'latest-fully-verified-checkpoint.json',q)
            wrong=copy.deepcopy(anchor);wrong['callback']=put(Path(tmp)/'other-callback.json',{})
            with self.assertRaises(AssertionError):typed_strict(directory,wrong,'M',2)
            (directory/'latest-fully-verified-checkpoint.json').unlink()
            with self.assertRaises(FileNotFoundError):typed_strict(directory,anchor,'M',2)

    def test_flags_or_prepared_PASS_cannot_admit_F(self):
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);put(root/'inputs/selected.json',dict(retained_attempts={'M':{'2':dict(anchor={})}}))
            put(root/'science/M2/verification.json',dict(status='PASS'))
            put(root/'science/M2/final-strict/result.json',dict(status='PASS'))
            with self.assertRaises(FileNotFoundError):require_M2_gate(root)

    def test_strict_completed_F_commit_is_a_distinct_supported_role(self):
        with tempfile.TemporaryDirectory()as tmp:
            directory,anchor,r,cp,q=startup_fixture(Path(tmp))
            r.update(dt=.00125,steps=4,last_stage_alpha=800.)
            q.update(schema='WJQualificationCommitV1',**{k:anchor[k]for k in ('report','callback','checkpoint')})
            q['binding'].update(absolute_step=4,alpha_bits=bits(800.))
            put(directory/'step-4.qualification-commit.json',q);put(directory/'latest-fully-verified-checkpoint.json',q)
            with patch('wj_h22_completion.read_completed',return_value=(r,cp)):
                request=typed_strict(directory,anchor,'F',4)
            self.assertEqual(request['import_role'],'QualifiedV1')
            self.assertTrue(request['qualified_cursor']['path'].endswith('step-4.qualification-commit.json'))

    def test_old_and_new_capture_sidecars_are_distinct_exact_bindings(self):
        offsets,r,p,P=actual_fixture();states,_=phases(r,p,P,offsets)
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);one=root/'old';two=root/'new';one.mkdir();two.mkdir()
            old,op,_=vector_fixture(one,states,r);new,np,_=vector_fixture(two,states,r)
            for path,packet in ((old,op),(new,np)):
                flow_vectors(identity(path),states,r,offsets,packet['execution_identity'])
            self.assertNotEqual(op['execution_identity'],np['execution_identity'])
            with self.assertRaises(AssertionError):flow_vectors(identity(old),states,r,offsets,np['execution_identity'])

    def test_startup_only_observer_assesses_step2_and_stays_alive(self):
        from wj_h22_completion_observer import daemon
        with tempfile.TemporaryDirectory()as tmp:
            root=Path(tmp);directory=root/'monitor';directory.mkdir()
            completion=put(root/'completion.json',dict(binding=dict(absolute_step=2),callback={}))
            selected=put(root/'selected.json',dict(retained_attempts={'M':{'2':dict(anchor={'checkpoint':put(root/'state.bin',{})})}},
                oracles={'M':{'2':{}}},auxiliary_audits={'M':{'2':{}}},common_support=put(root/'support.json',dict(points=[]))))
            intent=dict(schema='WJMonitorIntentV1',completion_role='RetainedM2StartupV1',level='M',maximum_absolute_step=2,restored_absolute_step=2,
                scientific_stride=2,recovery_stride=1,dt_bits=bits(.0025),substep_evidence_steps=[],flow_capture='Off',ack_timeout_seconds=120.,
                restart_completion=completion,observer_pid=__import__('os').getpid(),directory=str(directory),nonce='fixture-nonce-123456789',selected_inputs=selected)
            ip=root/'intent.json';put(ip,intent);_,_,current,_=records()
            ready=dict(report=put(root/'startup.json',current),structure=put(root/'structure.json',dict(fields=[],offsets=[])),binding={'absolute_step':2})
            put(directory/'startup.ready.json',ready)
            put(directory/'startup.verification-commit.json',dict(schema='WJStartupVerificationCommitV1',completion=completion))
            rows=Mock(rows=[],point_count=0);cache=Mock(report={'fields':[],'offsets':[]});cache.view.return_value=(rows,{}, {})
            with patch('wj_h22_completion_observer.load_compatible',return_value=(cache,{'cache':{}})),\
                 patch('wj_h22_completion_observer.checkpoint_payload',return_value={}),\
                 patch('wj_h22_completion_observer.restored_parity',return_value={'status':'PASS'})as parity,\
                 patch('wj_h22_completion_observer.check_current')as current_check,\
                 patch('wj_h22_completion_observer.common_points',return_value=[]),\
                 patch('wj_h22_completion_observer.assess',return_value=({'summary':{'status':'PASS'}},{})),\
                 patch('wj_h22_completion_observer.prefix')as prefix,patch('wj_h22_completion_observer.wait_stop')as wait:
                daemon(root,ip)
                parity.assert_called_once();current_check.assert_called_once_with(current,cache,2)
                prefix.assert_not_called();wait.assert_called_once_with(directory)
            self.assertTrue((directory/'startup-observer-complete.json').exists())
            self.assertFalse((directory/'step-2.ack.json').exists())

if __name__=='__main__':unittest.main()
