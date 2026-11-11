#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused execution-resource and complete transaction metadata fixtures."""
import copy
from pathlib import Path
import tempfile
import unittest

from wj6_resources import DEFAULT, validate
from run_wj5_native import executed_m1_obstruction
from verify_wj4_transactions import verify
from wj_immutable import write_once,artifacts
from wj6_qualification import sum_envelopes,common_extension_authorized,required_accuracy,fresh_startup_gate


class Metadata(unittest.TestCase):
    def test_observation_resume_cannot_reuse_failed_or_governing_receipts(self):
        from qualify_wj6_observations import reusable_observation,failure_status
        good=dict(exit_code=0,actual_attempts=0,artifacts=dict(output=dict(sha256='a'*64)))
        self.assertTrue(reusable_observation(good))
        for bad in (dict(good,exit_code=125),dict(good,actual_attempts=1),dict(good,artifacts={})):
            self.assertFalse(reusable_observation(bad))
        self.assertEqual(failure_status('physical audit failure',dict(exit_code=1,failure=None)),'OBSTRUCTED')
        self.assertEqual(failure_status('readback failed',dict(failure='RESOURCE_LIMIT: worker denial')),'RESOURCE_LIMIT')
    def test_resource_denial_alias_requires_no_launched_attempt(self):
        from resume_wj6_transaction import prelaunch_alias
        good=dict(exit_code=125,actual_attempts=0,resource=dict(not_launched=True))
        self.assertEqual(prelaunch_alias(good,'readback'),'readback-worker-resume1')
        for changed in (dict(good,actual_attempts=1),dict(good,exit_code=134),dict(good,resource={}),None):
            self.assertEqual(prelaunch_alias(changed,'readback'),'readback')
    def test_independent_order_row_gate_cannot_be_missing_or_failed(self):
        from qualify_wj6_observations import complete_orders
        self.assertEqual(complete_orders(dict(status='PASS'), None)['status'], 'NOT_RUN')
        self.assertEqual(complete_orders(dict(status='PASS'), dict(status='UNRESOLVED'))['status'], 'UNRESOLVED')
        self.assertEqual(complete_orders(dict(status='OBSTRUCTED'), dict(status='PASS'))['status'], 'OBSTRUCTED')
        self.assertEqual(complete_orders(dict(status='PASS'), dict(status='PASS'))['status'], 'PASS')
    def test_fresh_startup_requires_both_immutable_inputs(self):
        self.assertEqual(fresh_startup_gate(None)['status'],'NOT_RUN')
        with tempfile.TemporaryDirectory() as temp:
            paths=[Path(temp)/n for n in ('coarse','fine')]
            for p in paths:p.write_bytes(p.name.encode())
            packet=dict(exact_physical_primary_gauge_history=True,native_labels_differ=True,checkpoints=artifacts(paths))
            self.assertEqual(fresh_startup_gate(packet)['status'],'PASS')
            bad=copy.deepcopy(packet);bad['native_labels_differ']=False
            self.assertEqual(fresh_startup_gate(bad)['status'],'FAIL')
            paths[0].write_bytes(b'stale')
            with self.assertRaises(ValueError):fresh_startup_gate(packet)
    def test_both_asymmetric_uncertainty_contributions(self):
        self.assertEqual(sum_envelopes(dict(u=dict(L2=2.,sampled_Linf=3.)),dict(u=dict(L2=7.,sampled_Linf=11.))),
                         dict(u=dict(L2=9.,sampled_Linf=14.)))
        with self.assertRaises(ValueError):sum_envelopes(dict(u=dict(L2=2.)),{})

    def test_no_missing_or_failed_primary_extension(self):
        good=dict(primary=dict(status='PASS'),transactions=dict(status='PASS'),extension_resource=dict(status='PASS'))
        self.assertTrue(common_extension_authorized(dict(EP=good,Budget=good)))
        self.assertFalse(common_extension_authorized(dict(EP=good)))
        for key in good:
            bad=copy.deepcopy(good);bad[key]['status']='UNRESOLVED'
            self.assertFalse(common_extension_authorized(dict(EP=good,Budget=bad)))
        records={'0':dict(errors=dict(accuracy_gate=dict(status='PASS'))),'4':dict(errors=dict(accuracy_gate=dict(status='OBSTRUCTED'))),'32':dict(errors=dict(accuracy_gate=dict(status='PASS')))}
        self.assertEqual(required_accuracy(records,['0','4','32'])['status'],'OBSTRUCTED')
        del records['4']
        self.assertNotEqual(required_accuracy(records,['0','4','32'])['status'],'PASS')
    def test_finite_resource_design(self):
        validate(DEFAULT)
        for key,value in [('rss_bytes',float('nan')),('native_seconds',float('inf')),
                          ('rss_bytes',8.5),('server_workers',73),('native_seconds',False)]:
            data=dict(DEFAULT);data[key]=value
            with self.assertRaises(ValueError):validate(data)
        with self.assertRaises(ValueError):validate({})

    def test_missing_m1_never_authorizes_m2(self):
        self.assertFalse(executed_m1_obstruction([{},{}]))
        self.assertFalse(executed_m1_obstruction([dict(screens=dict(M1=dict(status='FAIL')))]))
        entry=dict(status='OBSTRUCTED',terminal=dict(steps=32,sha256='a'*64))
        self.assertTrue(executed_m1_obstruction([dict(screens=dict(M1=entry)),{}]))
        entry['status']='PASS'
        self.assertFalse(executed_m1_obstruction([dict(screens=dict(M1=entry))]))

    def test_eight_and_legacy_thirtytwo_step_transactions(self):
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp)
            for n in (8,32):
                attempts=[]
                for i in range(n):
                    path=root/('k-'+str(n)+'-'+str(i)+'.json');write_once(path,dict(raw=[float(i)],corrected=[float(i)]))
                    attempts.append(dict(history_sha256=str(i),incoming_sha256=str(i),time=float(i+1),alpha=1.,
                        governing_solves_complete=True,scheduled_complete=True,momentum_predictor_solves=1,
                        pressure_system_solves=1,k_transport_solves=1,epsilon_transport_solves=1,bound_calls=1,
                        bound_activations=0,additional_governing_equation_solves=0,
                        additional_correction_or_restoration_solves=0,k_snapshot=str(path)))
                whole=dict(accepted_states=[[float(i)] for i in range(n+1)],audit_times=list(map(float,range(n+1))),
                    state_audits=[dict(test=i) for i in range(n+1)],wj_audit_packets=[dict(test=i) for i in range(n+1)],
                    pimple_attempts=attempts,state=[float(n)],native_history=[[float(n)]],native_history_times=[float(n)],
                    time=float(n),steps=n,physical_samples=dict(value=1.))
                first=copy.deepcopy(whole);second=copy.deepcopy(whole)
                for key in ('accepted_states','audit_times','state_audits','wj_audit_packets'):
                    first[key]=first[key][:n//2+1];second[key]=second[key][n//2:]
                first['pimple_attempts']=first['pimple_attempts'][:n//2];second['pimple_attempts']=second['pimple_attempts'][n//2:]
                result=verify(whole,first,second,n) if n==8 else verify(whole,first,second)
                self.assertEqual(result['accepted_states_including_startup'],n+1)
                bad=copy.deepcopy(second);bad['wj_audit_packets'][1]['test']=-1
                with self.assertRaises(ValueError):verify(whole,first,bad,n)
                with self.assertRaises(ValueError):verify(whole,first,second,float(n))


if __name__=='__main__':unittest.main()
