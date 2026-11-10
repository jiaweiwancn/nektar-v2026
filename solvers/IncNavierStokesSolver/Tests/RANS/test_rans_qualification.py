#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Failure-contract checks for scientific qualification receipts."""
import copy
import unittest
from rans_qualification import coverage,compare_samples

class FailureContract(unittest.TestCase):
    def setUp(self):
        self.a=dict(time=.0016,grid=[[0,0,1],[0,0,3]],grid_regions=[0,1],grid_values=[[1,0,1,0,2],[1,0,1,.1,-2]])
    def test_coverage_missing_duplicate_unexpected(self):
        want={'a','b'}
        for names in (['a'],['a','b','b'],['a','b','c']):self.assertEqual(coverage([dict(name=x) for x in names],want,key=lambda x:x['name'])['status'],'fail')
        self.assertEqual(coverage([dict(name=x) for x in ('a','b')],want,key=lambda x:x['name'])['status'],'pass')
    def test_reject_comparison_mismatch(self):
        for key,value in [('time',.0015),('grid',[[0,0,1]]),('grid_values',[[1,0,1,0,2]]),('grid_regions',[1,0]),('grid_weights',[1,2])]:
            b=copy.deepcopy(self.a);b[key]=value
            with self.assertRaises(ValueError):compare_samples(self.a,b,[1.,1.])
    def test_common_pressure_gauge_preserves_jump(self):
        b=copy.deepcopy(self.a)
        for v in b['grid_values']:v[4]+=7
        self.assertEqual(compare_samples(self.a,b,[1.,1.])['fields']['Pi']['weighted_RMS'],0.)
        b['grid_values'][0][4]+=2
        self.assertEqual(compare_samples(self.a,b,[1.,1.])['fields']['Pi']['weighted_RMS'],1.)
    def test_L2_vs_RMS(self):
        b=copy.deepcopy(self.a);b['grid_values'][0][0]+=1;b['grid_values'][1][0]+=1
        f=compare_samples(self.a,b,[2.,2.])['fields']['u'];self.assertEqual(f,dict(weighted_RMS=1.,L2=2.))

if __name__=='__main__':unittest.main()
