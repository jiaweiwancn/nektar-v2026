#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused semantic-producer, import and joint both-order coverage fixtures."""
import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from publish_wj6 import execution_stage_binding
from wj7_lineage import verify_import,bind_producer
from wj7_qualification import expected_tuples,tuple_gate,joint_ready,final_gate,load_measurement
from wj6_qualification import sum_envelopes
from wj_immutable import file_identity


class Metadata(unittest.TestCase):
    def test_both_order_comparisons_integrate_the_complete_boundary(self):
        from wj4_measure import rule,comparison
        from wj7_qualification import consistent_comparison
        for order,xorder in ((8,4),(10,6)):
            settings=dict(order=order,xorder=xorder)
            rows=[]
            for x,y,z,w,side,*extra in rule(**settings):
                rows.append(dict(position=[x,y,z],weight=w,side=side,inner=z<1 or (z==1 and side<0),
                    value=[0.]*6,represented_initial=[0.]*6,geometry=0,element=0,requested_x_side=0,
                    face_sides_xz=[0,0],provenance='fixture'))
            m=dict(observations=rows)
            result=consistent_comparison(m,m,'EP',settings)
            self.assertAlmostEqual(result['constrained_boundary']['fields']['velocity']['measure'],5.)
            self.assertAlmostEqual(result['constrained_boundary']['fields']['outer_epsilon']['measure'],3.)
            if order==10:
                old=comparison(m,m,'EP',settings)
                self.assertAlmostEqual(old['constrained_boundary']['fields']['velocity']['measure'],4.)
    def test_real_publisher_stage_report_distinct_from_arithmetic_matrix(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);matrix=root/'epsilon.metric.bin';matrix.write_bytes(bytes(32))
            descriptor=root/'epsilon.json';f=dict(file_identity(matrix),shape=[2,2])
            descriptor.write_text(json.dumps(dict(matrix=f,mass_or_correction_map=f)))
            report=root/'report.json'
            for arithmetic in (False,True):
                attempt=dict(time=1.)
                if arithmetic:attempt['epsilon_arithmetic']=file_identity(descriptor)
                report.write_text(json.dumps(dict(pimple_attempts=[attempt],state=[1.],accepted_states=[[0.],[1.]],
                    checkpoint='CP1',last_stage_history=[.5],last_stage_time=1.,last_stage_alpha=2.,wj_audit_packets=[])))
                producer=file_identity(report);d,stage,allocations=execution_stage_binding(producer)
                stages={producer['path']:stage}
                self.assertEqual(list(stages),[str(report)])
                self.assertEqual(stage['report'],producer)
                self.assertNotEqual(stage['report'],file_identity(matrix))
                self.assertEqual(len(allocations),2 if arithmetic else 0)
            with self.assertRaises(ValueError):execution_stage_binding(file_identity(descriptor))

    def test_import_retains_original_reference_and_signed_history(self):
        original=dict(dt=3.2e-5,signature='physical',steps=32,time=.001024,state=[1.,-0.],
            native_history=[[1.,-0.]],native_history_times=[.001024],fixed_reference=[.5,0.],
            last_stage_history=[2.,-0.],last_stage_alpha=31250.,last_stage_time=.001024)
        actual=copy.deepcopy(original);actual.update(last_stage_alpha=0,pimple_attempts=[],implicit_callback_count=0,
            initialization_lift_calls_iterations=[0,0])
        self.assertEqual(verify_import(original,actual,original)['status'],'PASS')
        for key,value in [('dt',1.6e-5),('state',[1.,0.]),('native_history',[[0.,-0.]]),
                          ('fixed_reference',[.6,0.]),('steps',0),('time',0.),('signature','wrong')]:
            bad=copy.deepcopy(actual);bad[key]=value
            with self.assertRaises(ValueError):verify_import(original,bad,original)
        bad=copy.deepcopy(actual);bad['initialization_lift_calls_iterations']=[1,1]
        with self.assertRaises(ValueError):verify_import(original,bad,original)

    def test_wrong_producer_is_not_accepted_by_matching_path(self):
        with tempfile.TemporaryDirectory() as temp:
            f=Path(temp)/'report.json';f.write_text('{}');identity=file_identity(f)
            c=dict(report=str(f),report_sha256=identity['sha256'])
            with self.assertRaisesRegex(ValueError,'producer_mismatch'):bind_producer(c,dict(entries={}))
            with self.assertRaisesRegex(ValueError,'producer_role'):bind_producer(c,dict(entries={str(f):dict(report=dict(identity,sha256='0'*64))}))

    def test_tuple_coverage_all_orders_times_fields_scopes_and_no_duplicates(self):
        rows=[dict(key=list(k)) for k in expected_tuples(128)]
        self.assertEqual(len(rows),1008);self.assertEqual(tuple_gate(rows,128)['status'],'PASS')
        for match in (lambda k:k[0]==10,lambda k:k[2]==96,lambda k:k[4]=='k',lambda k:k[3]=='central'):
            filtered=[r for r in rows if not match(r['key'])]
            self.assertNotEqual(tuple_gate(filtered,128)['status'],'PASS')
        self.assertNotEqual(tuple_gate(rows+[rows[0]],128)['status'],'PASS')

    def test_joint_barrier_and_final_target_cannot_be_skipped(self):
        good=dict(end=64,gate=dict(status='PASS'))
        self.assertTrue(joint_ready(dict(EP=good,Budget=good),64))
        for bad in (dict(end=32,gate=dict(status='PASS')),dict(end=64,gate=dict(status='NOT_RUN')),dict(end=64,gate=dict(status='UNRESOLVED'))):
            self.assertFalse(joint_ready(dict(EP=good,Budget=bad),64))
        self.assertNotEqual(final_gate(dict(EP=good,Budget=good))['status'],'PASS')
        final=dict(end=128,gate=dict(status='PASS'))
        self.assertEqual(final_gate(dict(EP=final,Budget=final))['status'],'PASS')

    def test_asymmetric_pair_uncertainties_are_both_retained(self):
        self.assertEqual(sum_envelopes(dict(k=dict(L2=2.,sampled_Linf=3.)),dict(k=dict(L2=11.,sampled_Linf=17.))),
            dict(k=dict(L2=13.,sampled_Linf=20.)))

    def test_absolute_steps_are_not_reset_at_package_import(self):
        with tempfile.TemporaryDirectory() as temp:
            f=Path(temp)/'producer.json';f.write_text(json.dumps(dict(steps=96,time=.003072)))
            c=dict(path=str(f),sha256=file_identity(f)['sha256'],report=str(f),report_sha256=file_identity(f)['sha256'],steps=64,time=.003072)
            with self.assertRaisesRegex(ValueError,'absolute native step'):load_measurement(c)


if __name__=='__main__':unittest.main()
