#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import copy
import json
from pathlib import Path
import tempfile
import unittest
from wj_qualification import wj3_gates
from wj_immutable import write_once,execution_intent,artifacts,verify_artifacts,Ledger

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
    def test_constant_is_distinct(self):
        self.s['temporal']['rates']['k'].update(A=1e-16,B=2e-16,classification='roundoff_limited_constant')
        self.assertEqual(wj3_gates(self.s)['status'],'PASS')

class Immutable(unittest.TestCase):
    def test_changed_intents_and_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);files=[]
            for name in ('case.xml','State.bin','solver','report.json','log'):
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
