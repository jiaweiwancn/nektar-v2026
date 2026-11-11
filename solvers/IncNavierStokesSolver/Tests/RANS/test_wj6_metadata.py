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
from wj_immutable import write_once
from wj6_qualification import sum_envelopes,common_extension_authorized,required_accuracy


class Metadata(unittest.TestCase):
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
