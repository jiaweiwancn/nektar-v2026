#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
import json
import struct
import sys
from types import SimpleNamespace
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from wj_qualification import wj3_gates,aggregate
from wj_immutable import write_once,execution_intent,artifacts,verify_artifacts,Ledger,xml_identity,read

class Gates(unittest.TestCase):
    def setUp(self):
        self.s=json.loads((Path(__file__).resolve().parents[4]/'docs/plans/wj3-linux/EP.json').read_text())
    def test_existing_evidence(self):
        for policy in ('EP','Budget'):
            s=json.loads((Path(__file__).resolve().parents[4]/('docs/plans/wj3-linux/'+policy+'.json')).read_text())
            self.assertEqual(wj3_gates(s)['status'],'PASS')
    def test_false_pass_fixtures(self):
        for key,value,status in [('cross_obstruction',{'reason':'native hard failure'},'OBSTRUCTED'),('transactions',{},'NOT_RUN'),('main',[],'NOT_RUN')]:
            s=copy.deepcopy(self.s);s[key]=value;self.assertEqual(wj3_gates(s)['status'],status)
        for a,b,status in [(1e-14,1e-14,'UNRESOLVED'),(1e-8,2e-8,'FAIL'),(float('nan'),1e-8,'FAIL'),(float('inf'),1e-8,'FAIL')]:
            s=copy.deepcopy(self.s);s['temporal']['rates']['velocity'].update(A=a,B=b)
            self.assertEqual(wj3_gates(s)['status'],status)
        s=copy.deepcopy(self.s);s['spatial'][2]['raw_reference_norms'][0]=1.
        self.assertEqual(wj3_gates(s)['status'],'FAIL')
        s=copy.deepcopy(self.s);s['spatial']=s['spatial'][:2]
        self.assertEqual(wj3_gates(s)['status'],'NOT_RUN')
    def test_incomplete_blocks_and_negative_norms(self):
        s=copy.deepcopy(self.s);s['spatial'][0]['raw_reference_norms']=[]
        self.assertNotEqual(wj3_gates(s)['status'],'PASS')
        s=copy.deepcopy(self.s);s['spatial'][0]['raw_reference_norms'][0]=-1.
        self.assertEqual(wj3_gates(s)['status'],'FAIL')
    def test_incomplete_native_targets(self):
        s=copy.deepcopy(self.s);s['main'][0]['steps']=1
        self.assertEqual(wj3_gates(s)['status'],'NOT_RUN')
        s=copy.deepcopy(self.s);s['main'][0]['comparison']=None
        self.assertEqual(wj3_gates(s)['status'],'NOT_RUN')
    def test_unknown_component_cannot_pass(self):
        for value in ({},{'status':'UNKNOWN'},None):
            self.assertEqual(aggregate({'required':value})['status'],'FAIL')
    def test_constant_is_distinct(self):
        self.s['temporal']['rates']['k'].update(A=1e-16,B=2e-16,classification='roundoff_limited_constant')
        self.assertEqual(wj3_gates(self.s)['status'],'PASS')

class Immutable(unittest.TestCase):
    def test_legacy_control_budget_units(self):
        from wj4_validate import legacy_control_budget
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'case.xml'
            for n,expected in ((0,1),(2,2)):
                p.write_text('<N><PARAMETERS><P>NumSteps = '+str(n)+'</P></PARAMETERS></N>')
                tests=[dict(name='control',command=['/bin/IncNavierStokesSolver',str(p)]),dict(name='kernel',command=['/bin/WJKernelProbe'])]
                rows=legacy_control_budget(tests);self.assertEqual(len(rows),1);self.assertEqual(rows[0]['attempts'],expected)
                self.assertEqual(rows[0]['policy'],'SharedLegacyControl')

    def test_transaction_intermediate_packets_are_required(self):
        from verify_wj4_transactions import verify
        with tempfile.TemporaryDirectory() as temp:
            attempts=[]
            for i in range(32):
                path=Path(temp)/(str(i)+'.json');path.write_text(json.dumps(dict(raw=[i],corrected=[i])))
                attempts.append(dict(history_sha256=str(i),incoming_sha256=str(i),time=i+1,alpha=1.,governing_solves_complete=True,
                    scheduled_complete=True,momentum_predictor_solves=1,pressure_system_solves=1,k_transport_solves=1,
                    epsilon_transport_solves=1,bound_calls=1,bound_activations=0,additional_governing_equation_solves=0,
                    additional_correction_or_restoration_solves=0,k_snapshot=str(path)))
            whole=dict(accepted_states=[[i] for i in range(33)],audit_times=list(range(33)),
                state_audits=[dict(k=i) for i in range(33)],wj_audit_packets=[dict(k=i) for i in range(33)],pimple_attempts=attempts,
                state=[32.],native_history=[[32.]],native_history_times=[32.],time=32.,steps=32,physical_samples=dict(k=[32.]))
            first=copy.deepcopy(whole);second=copy.deepcopy(whole)
            for key in ('accepted_states','audit_times','state_audits','wj_audit_packets'):
                first[key]=first[key][:17];second[key]=second[key][16:]
            first['pimple_attempts']=first['pimple_attempts'][:16];second['pimple_attempts']=second['pimple_attempts'][16:]
            self.assertEqual(verify(whole,first,second)['status'],'PASS')
            for key in ('accepted_states','audit_times','state_audits','wj_audit_packets','pimple_attempts'):
                changed=copy.deepcopy(second)
                if key=='pimple_attempts':changed[key][4]['history_sha256']='changed'
                else:changed[key][4]='changed'
                with self.assertRaisesRegex(ValueError,'mismatch'):verify(whole,first,changed)

    def test_orphan_command_log_is_retained(self):
        from wj3_receipts import command
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);log=root/'interrupted.log';log.write_bytes(b'failed attempt witness\n')
            args=SimpleNamespace(root=root,label='interrupted',command=[sys.executable,'-c','print("replacement")'])
            with self.assertRaises(FileExistsError):command(args)
            self.assertEqual(log.read_bytes(),b'failed attempt witness\n')
            self.assertFalse((root/'command-receipts.json').exists())

    def test_changed_intents_and_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);files=[]
            for name in ('case.xml','State.bin','solver','report.json','log','arithmetic.json','matrix.bin','metric.bin'):
                p=root/name;p.write_bytes(name.encode());files.append(p)
            base=execution_intent(['solver','case.xml'],files,{'code':'v1'},{'runtime':'v1'},{'dt':1e-5,'mesh':[0,.5,1]})
            write_once(root/'intent.json',base);write_once(root/'intent.json',base)
            for key in ('argv','runtime','code','settings'):
                changed=copy.deepcopy(base);changed[key]='different'
                with self.assertRaisesRegex(ValueError,'intent/design mismatch'):write_once(root/'intent.json',changed)
            for file in files:
                index=artifacts(files);old=file.read_bytes();file.write_bytes(old+b'stale')
                with self.assertRaisesRegex(ValueError,'artifact mismatch'):verify_artifacts(index)
                file.write_bytes(old)
            write_once(root/'design.json',{'T':.008192})
            with self.assertRaisesRegex(ValueError,'intent/design mismatch'):write_once(root/'design.json',{'T':.004096})
    def test_native_negative_zero(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'state.json';p.write_text('{"state":[-0,0,-0.0]}')
            values=read(p)['state']
            self.assertEqual(struct.pack('<d',values[0]),struct.pack('<d',-0.0))
            self.assertEqual(struct.pack('<d',values[1]),struct.pack('<d',0.0))
            self.assertEqual(struct.pack('<d',values[2]),struct.pack('<d',-0.0))
    def test_xml_bytes_and_semantics(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)/'case.xml';p.write_text('<N><P>TimeStep=1e-5</P></N>')
            first=xml_identity(p);p.write_text('<N>\n<P>TimeStep=1e-5</P></N>')
            second=xml_identity(p)
            self.assertNotEqual(first['sha256'],second['sha256'])
            self.assertEqual(first['semantic_sha256'],second['semantic_sha256'])
            p.write_text('<N><P>TimeStep=2e-5</P></N>')
            self.assertNotEqual(first['semantic_sha256'],xml_identity(p)['semantic_sha256'])
    def test_session_wall_time_is_not_reset_or_double_counted(self):
        with tempfile.TemporaryDirectory() as temp:
            l=Ledger(temp,'EP',cap=3,wall_cap=200)
            with patch('wj_immutable.time.time',return_value=1000):l.begin_session('first')
            with patch('wj_immutable.time.time',return_value=1010):
                l.reserve('command',{},1,100);l.finish('command',1,5,'receipt')
                self.assertEqual(l.totals()['seconds'],10)
            with patch('wj_immutable.time.time',return_value=1020):l.finish_session('first')
            l=Ledger(temp,'EP',cap=3,wall_cap=200)
            self.assertEqual(l.totals()['seconds'],20)
            with patch('wj_immutable.time.time',return_value=2000):l.begin_session('second')
            with patch('wj_immutable.time.time',return_value=2030):
                self.assertEqual(l.totals()['seconds'],50)
                l.finish_session('second')
            self.assertEqual(l.totals()['seconds'],50)
    def test_persistent_reservations(self):
        with tempfile.TemporaryDirectory() as temp:
            l=Ledger(temp,'EP',cap=3,wall_cap=200);l.reserve('a',{},2,100)
            l=Ledger(temp,'EP',cap=3,wall_cap=200)
            with self.assertRaisesRegex(RuntimeError,'RESOURCE_LIMIT'):l.reserve('b',{},2,1)
            with self.assertRaisesRegex(RuntimeError,'RESOURCE_LIMIT'):l.reserve('b',{},1,101)
            l.finish('a',1,50,'receipt');l.reserve('b',{},2,100)
            self.assertEqual(l.totals()['attempts'],3)
            with self.assertRaisesRegex(ValueError,'already reserved'):l.reserve('b',{},0,1)
            with self.assertRaisesRegex(ValueError,'limit mismatch'):Ledger(temp,'EP',cap=4,wall_cap=200)

if __name__=='__main__':unittest.main()
