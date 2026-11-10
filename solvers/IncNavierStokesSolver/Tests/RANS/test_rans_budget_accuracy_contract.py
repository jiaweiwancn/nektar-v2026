#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import unittest
from classify_rans_physical_budget_accuracy import classify_row


class BudgetAccuracy(unittest.TestCase):
    def row(self,values):
        return dict(kind='perturbed',equation='whole_patch_k',quadrature_uncertainty=1e-12,
            signed_values=dict(zip(('B5','B11','B5-dt2','B11-dt2'),values)),
            spatial_status='pass',temporal_status='pass')
    def test_accuracy_requires_both_reductions(self):
        result=classify_row(self.row([1e-4,1e-5,1e-4+1e-9,1e-5+1e-9]))
        self.assertEqual(result['spatial_status'],'pass');self.assertEqual(result['temporal_status'],'pass')
    def test_pass_labels_cannot_hide_failed_spatial_check(self):
        self.assertEqual(classify_row(self.row([1e-4,1.1e-4,1e-4,1.1e-4]))['spatial_status'],'failed_or_unresolved')
    def test_temporal_uncertainty_is_separate(self):
        self.assertEqual(classify_row(self.row([1e-4,1e-5,2e-4,1e-5]))['temporal_status'],'unresolved')
    def test_floor_does_not_establish_a_rate(self):
        self.assertEqual(classify_row(self.row([1e-11,2e-11,1e-11,2e-11]))['spatial_status'],'unresolved_floor')


if __name__=='__main__':unittest.main()
