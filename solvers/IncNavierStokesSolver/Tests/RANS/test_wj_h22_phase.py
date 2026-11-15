# SPDX-License-Identifier: MIT
import copy,hashlib,json,sys,tempfile,unittest
import xml.etree.ElementTree as ET
from pathlib import Path
from unittest.mock import patch
import numpy as np
from wj_h15_observer import identity,bits,vector_hash
from wj_h22_replay import DT,END,MATRIX,selected,phases,flow_vectors,integer,read_completed
from wj_h22_phase import branch_window,cross_window
from run_wj_h22 import command_level
from study_wj_h22 import validate_intent,protected_case,xml_case
from wj_h22_observer import prepared_cursor
import test_wj_h20_three_level as h20_fixture

def actual_fixture():
    offsets=[0,2,4,6,8,9,10];P=np.array([1.,-0.,0.,0.,.009,.009,.002,.002,.001,.001])
    F=P.copy();F[[0,2]]+=[.2,.2];F[[8,9]]+=[.0001,-.0002]
    K=F.copy();K[4]+=1e-6;A=K.copy();A[6]+=2e-6
    report=dict(state=A.tolist(),fixed_reference=P.tolist(),time=.005,steps=1,last_stage_alpha=200.,last_stage_history=P.tolist())
    packet=dict(committed=P.tolist(),post_flow=F.tolist(),pre_k=F.tolist(),raw=K.tolist(),corrected=K.tolist(),history=P.tolist(),
        time=.005,alpha=200.,fixed_reference=P.tolist(),k_bound=dict(activations=0,mass_norm_correction=0))
    return offsets,report,packet,dict(state=P.tolist())

def vector_fixture(folder,states,report):
    arrays=dict(predictor=np.array([.3,-0.]),bGz=np.array([.1,-.2]),du=np.array([.2,.2]),dz=np.array([.0001,-.0002]))
    data=b'';records=[]
    for name,v in arrays.items():
        raw=np.asarray(v,dtype='<f8').tobytes();records.append(dict(name=name,byte_offset=len(data),count=len(v),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()));data+=raw
    payload=folder/'flow.bin';payload.write_bytes(data);execution=folder/'execution.json';execution.write_text('{"source":"fixture","ELF":"fixture"}')
    packet=dict(schema='WJMappedFlowVectorsV1',status='FLOW_PROVISIONAL_NOT_TIME_COMMIT',absolute_step=1,time_bits=bits(.005),alpha_bits=bits(200.),
        base_sha256=vector_hash(states['C']),flow_state_sha256=vector_hash(states['F']),callback_history_sha256=vector_hash(report['last_stage_history']),
        execution_identity=identity(execution),additional_governing_operations=0,pressure_velocity_sign='actual bGz; signed correction is -bGz',
        payload=identity(payload),arrays=records,velocity_rows=[0,2],pressure_rows=[8,9])
    path=folder/'flow.json';path.write_text(json.dumps(packet));return path,packet,arrays

class NativePhaseReplay(unittest.TestCase):
    def test_exact_seven_attempt_scope_without_CP0_or_bootstrap(self):
        self.assertEqual(len(MATRIX),7);self.assertEqual([z for z,n in MATRIX],['C','M','M','F','F','F','F'])
        root=Path('/home/jay/build/wj-h22-early-phase-replay-20261115');script=Path(__file__).with_name('study_wj_h22.py')
        base=[sys.executable,'-B',str(script),'--root',str(root),'--mode']
        for z,n in END.items():
            self.assertEqual(command_level(root,base+['target','--level',z],'target',n),z)
            with self.assertRaises(AssertionError):command_level(root,base+['target','--level',z],'target',n+1)
        with self.assertRaises(AssertionError):command_level(root,base+['bootstrap'],'bootstrap',0)

    def test_own_integer_step_time_callback_and_alpha(self):
        for z,n in MATRIX:
            r=h20_fixture.ThreeLevel().report(z,n);r.update(state=r['accepted_states'][0],time=n*DT[z])
            self.assertEqual(selected(r,z,n)['step'],n)
            bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['time']+=DT[z]
            with self.assertRaises((AssertionError,ValueError)):selected(bad,z,n)
            bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['stage_history_sha256']='other producer'
            with self.assertRaises(ValueError):selected(bad,z,n)
            bad=copy.deepcopy(r);bad['wj_audit_packets'][0]['alpha_mass']=0
            with self.assertRaises(AssertionError):selected(bad,z,n)
        self.assertEqual(integer('4'),4)
        for value in (1.5,'1.5','01',True,float('nan')):
            with self.assertRaises((AssertionError,ValueError)):integer(value)

    def test_actual_masks_aliases_and_provisional_K_epsilon(self):
        offsets,r,p,P=actual_fixture();states,b=phases(r,p,P,offsets)
        self.assertEqual(b['field_masks'],{'C->F':[0,1,4,5],'F->K':[2],'K->A':[3]})
        np.testing.assert_array_equal(states['K'][6:8],states['F'][6:8])
        self.assertNotEqual(states['K'][6],states['A'][6])
        bad=copy.deepcopy(p);bad['raw'][0]+=1;bad['corrected'][0]+=1
        with self.assertRaises(AssertionError):phases(r,bad,P,offsets)
        bad=copy.deepcopy(P);bad['state'][1]=0.
        with self.assertRaises(AssertionError):phases(r,p,bad,offsets)

    def test_actual_bGz_sign_vectors_and_roundoff_not_bitwise_subtraction(self):
        offsets,r,p,P=actual_fixture();states,_=phases(r,p,P,offsets)
        with tempfile.TemporaryDirectory()as d:
            path,packet,v=vector_fixture(Path(d),states,r)
            increments,proof=flow_vectors(identity(path),states,r,offsets,packet['execution_identity'])
            self.assertEqual(increments['pressure'][0],-.1);self.assertEqual(increments['pressure'][2],.2)
            self.assertTrue(np.signbit(v['predictor'][1]))
            self.assertGreater(proof['checks']['native_du']['maximum_error'],0.)
            self.assertLess(proof['checks']['native_du']['maximum_error'],proof['checks']['native_du']['maximum_bound'])
            raw=(Path(d)/'flow.bin').read_bytes();wrong=bytearray(raw);wrong[16:32]=np.array([-.1,.2],dtype='<f8').tobytes()
            (Path(d)/'flow.bin').write_bytes(wrong);packet['payload']=identity(Path(d)/'flow.bin')
            packet['arrays'][1]['sha256']=hashlib.sha256(wrong[16:32]).hexdigest();path.write_text(json.dumps(packet))
            with self.assertRaises(ValueError):flow_vectors(identity(path),states,r,offsets,packet['execution_identity'])

    def test_vector_source_history_row_and_corruption_refusals(self):
        offsets,r,p,P=actual_fixture();states,_=phases(r,p,P,offsets)
        with tempfile.TemporaryDirectory()as d:
            path,original,v=vector_fixture(Path(d),states,r)
            for key,value in [('callback_history_sha256','wrong'),('absolute_step',1.5),('velocity_rows',[0,8]),('velocity_rows',[0,2.5]),('base_sha256','wrong'),
                              ('execution_identity',dict(path='other-source',bytes=1,sha256='wrong'))]:
                packet=copy.deepcopy(original);packet[key]=value;path.write_text(json.dumps(packet))
                with self.assertRaises((AssertionError,ValueError)):flow_vectors(identity(path),states,r,offsets,original['execution_identity'])
            path.write_text(json.dumps(original));(Path(d)/'flow.bin').write_bytes(b'truncated')
            with self.assertRaises(AssertionError):flow_vectors(identity(path),states,r,offsets,original['execution_identity'])

    def test_complete_window_zero_incoming_and_all_split_terms(self):
        def row(start,flow,k,epsilon):
            values=[start,start,start+flow,start+flow+k,start+flow+k+epsilon]
            return {p:dict(mhat=v)for p,v in zip(('P','C','F','K','A'),values)}
        coarse=branch_window([row(1.,.4,-.3,.1)],'mhat')
        middle=branch_window([row(1.,.1,-.05,.04),row(1.09,.15,-.1,.03)],'mhat')
        pair=cross_window(coarse,middle,'mhat')
        self.assertEqual(float(pair['incoming']),0.)
        self.assertGreater(abs(float(pair['contributions']['flow'])),0)
        self.assertGreater(abs(float(pair['contributions']['k'])),0)
        self.assertGreater(abs(float(pair['contributions']['epsilon'])),0)
        self.assertLess(pair['identity']['maximum_error'],1e-15)
        changed=copy.deepcopy(middle);changed['start']=np.asarray(1.01)
        with self.assertRaises(AssertionError):cross_window(coarse,changed,'mhat')
        undefined=branch_window([{'P':{},'C':{},'F':{},'K':{},'A':{}}],'mhat')
        self.assertEqual(undefined['status'],'NOT_EVALUABLE')

    def test_intent_timeout_and_all_selected_capture_steps(self):
        for z in DT:
            v=dict(schema='WJMonitorIntentV1',level=z,maximum_absolute_step=END[z],restored_absolute_step=0,
                recovery_stride=1,scientific_stride=END[z],dt_bits=bits(DT[z]),substep_evidence_steps=list(range(1,END[z]+1)),ack_timeout_seconds=120.)
            validate_intent(v)
            for timeout in (0.,120.001):
                bad=dict(v,ack_timeout_seconds=timeout)
                with self.assertRaises(AssertionError):validate_intent(bad)
            with self.assertRaises(AssertionError):validate_intent(dict(v,substep_evidence_steps=[]))

    def test_prepared_sessions_preserve_own_dt_and_functions(self):
        with tempfile.TemporaryDirectory()as d:
            for z,dt in DT.items():
                source=Path(d)/(z+'-source.xml')
                source.write_text('<NEKTAR><CONDITIONS><PARAMETERS><P>Time = 0</P><P>TimeStep = '+str(dt)+'</P><P>NumSteps = 0</P></PARAMETERS>'
                    '<SOLVERINFO><I PROPERTY="RANSWallPimpleVerify" VALUE="Off"/></SOLVERINFO>'
                    '<FUNCTION NAME="InitialConditionsOuter"><E VAR="u" VALUE="log(z-.0875)"/></FUNCTION></CONDITIONS></NEKTAR>')
                path=xml_case(source,Path(d)/z,Path(d)/'existing-own-CP0',z)
                self.assertEqual(protected_case(ET.parse(source).getroot()),protected_case(ET.parse(path).getroot()))
                changed=ET.parse(path).getroot()
                for p in changed.findall('./CONDITIONS/PARAMETERS/P'):
                    if p.text.split('=')[0].strip()=='TimeStep':p.text='TimeStep = .123'
                self.assertNotEqual(protected_case(changed),protected_case(ET.parse(source).getroot()))

    def test_native_writer_is_copy_only_and_dense_rejection_survives(self):
        root=Path(__file__).resolve().parents[2]/'EquationSystems'
        helper=(root/'RANSWallFlowVectorCapture.h').read_text();sweep=(root/'RANSWallEvolutionPimple.cpp').read_text()
        block=sweep.split('if(m_pimpleFlowVectors && m_pimpleSnapshots)',1)[1].split('auto budget=',1)[0]
        for forbidden in ('Apply(', 'Multiply(', 'Factor(', 'Assemble'):
            self.assertNotIn(forbidden,helper);self.assertNotIn(forbidden,block)
        self.assertIn('predictor,bGz,du,dz',block)
        self.assertLess(sweep.index('const auto flowState=s;'),sweep.index('RANSWallFlowCapture::Write'))
        self.assertIn('MappedBlockV1 rejects legacy dense arithmetic/verification before allocation',(root/'RANSWallEvolutionPatch.cpp').read_text())

    def test_prepared_cursor_binds_actual_inventory_for_native_wait(self):
        binding=dict(state_sha256='actual-state',absolute_step='1')
        p=prepared_cursor(binding,dict(sha256='actual-ready'),dict(sha256='actual-inventory'),dict(sha256='actual-prefix-proof'))
        self.assertEqual(p['binding'],binding)
        self.assertEqual(p['ready_sha256'],'actual-ready')
        self.assertEqual(p['inventory_sha256'],'actual-inventory')
        self.assertEqual(p['prefix_parity']['sha256'],'actual-prefix-proof')
        self.assertEqual(p['role'],'PREPARED_NOT_NATIVE_COMMIT')

if __name__=='__main__':unittest.main()
