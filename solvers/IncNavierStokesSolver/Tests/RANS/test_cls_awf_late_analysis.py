#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Scientific estimator checks against known first-order and unresolved data."""
import unittest
from analyse_cls_awf_late_ivp import assess, estimate


class LateEstimates(unittest.TestCase):
    policy=dict(near_zero_normalized_difference=1e-11,coherent_order_interval=[.8,1.2])

    def test_known_linear_error_labels(self):
        # q_h=q_exact+C*h at h=4,2,1: the original coarse error is 4C.
        q=[4e-5,2e-5,1e-5];dc=q[0]-q[1];df=q[1]-q[2]
        d=estimate(dc,df,dc,df,1e-4,self.policy)
        self.assertAlmostEqual(d['p'],1.)
        for key,value in [('normalized_fine_estimate',q[2]),('normalized_middle_estimate',q[1]),('normalized_original_dt_estimate',q[0])]:
            self.assertAlmostEqual(d[key],value,places=15)
        self.assertTrue(d['original_dt_within_target'])

    def test_small_pair_difference_is_not_small_original_error(self):
        d=estimate(1e-4,5e-5,1e-4,5e-5,1e-4,self.policy)
        self.assertTrue(d['observed_pairwise_sensitivity_below_target'])
        self.assertFalse(d['original_dt_within_target'])

    def test_no_order_for_unresolved_signal(self):
        d=estimate(1e-12,6e-13,1e-12,6e-13,1e-4,self.policy)
        self.assertIsNone(d['p']);self.assertIsNone(d['normalized_original_dt_estimate'])
        self.assertTrue(d['observed_pairwise_sensitivity_below_target'])

    def test_unresolved_component_does_not_hide_resolved_error(self):
        large=estimate(1e-4,5e-5,1e-4,5e-5,1e-4,self.policy)
        tiny=estimate(1e-12,6e-13,1e-12,6e-13,1e-4,self.policy)
        for row in (large,tiny):row['observed_span_below_target']=True
        label,adequate=assess([large,tiny])
        self.assertFalse(adequate)
        self.assertEqual(label,'late_sensitivity_requires_further_diagnosis')


if __name__=='__main__':unittest.main()
