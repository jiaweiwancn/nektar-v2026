"""Focused metadata fixtures using the actual Audit.near implementation.

No Analyzer/field construction, source model, residual or native call occurs.
Run on the authoritative Linux host with one numerical thread.
"""
import math
import unittest

from cls_epawf_m8_r1_nonlinear import Audit
from cls_epawf_m8_r2_stage_analysis import Analyzer, ray_coordinate_correspondence


class RayCorrespondenceTests(unittest.TestCase):
    def setUp(self):
        self.audit = Audit.__new__(Audit)
        self.audit.checks = 0
        self.audit.errors = {}

    @staticmethod
    def row(radius=0.191552638660500013, weight=.1):
        return dict(inventory='volume_normalization_ray', geometry=17, block=2,
            side=0, strict=False, dependent_trace=False,
            coordinate=[-2.5, 0., radius], weight=weight)

    def check(self, native, independent):
        return ray_coordinate_correspondence(sorted(native, key=Analyzer.point_key),
            sorted(independent, key=Analyzer.point_key), self.audit, Analyzer.point_key)

    def test_one_ulp_crosses_rounding_bin_and_keeps_original_checks(self):
        native = self.row()
        independent = self.row(math.nextafter(native['coordinate'][2], -math.inf))
        self.assertNotEqual(Analyzer.point_key(native), Analyzer.point_key(independent))
        result = self.check([native], [independent])
        self.assertEqual(result['version'], 'v2')
        self.assertFalse(result['old_rounded_key_equality'])
        self.assertEqual(result['maximum_raw_coordinate_difference'], 2.7755575615628914e-17)
        self.assertTrue(result['additional_absolute_coordinate_safeguard_applicable'])
        self.assertEqual(result['additional_absolute_coordinate_safeguard'], 1e-13)
        self.assertEqual(result['original_coordinate_tolerance'], 2e-11)
        self.assertEqual(result['original_weight_tolerance'], 2e-11)
        self.assertIn('independent_complete_dynamic_ray_coordinates', self.audit.errors)
        self.assertIn('independent_complete_dynamic_ray_weights', self.audit.errors)

    def test_each_structural_field_must_remain_exact(self):
        changes = dict(inventory='other_inventory', geometry=18, block=3,
            side=1, strict=True, dependent_trace=True)
        for key, value in changes.items():
            with self.subTest(field=key):
                other = self.row(); other[key] = value
                with self.assertRaisesRegex(ValueError, 'structural inventory differs'):
                    self.check([self.row()], [other])

    def test_wrong_count_rejects(self):
        with self.assertRaisesRegex(ValueError, 'inventory lengths differ'):
            self.check([self.row(), self.row()], [self.row()])

    def test_changed_keys_outside_strict_raw_safeguard_reject(self):
        native, independent = self.row(.1915526386606), self.row(.1915526386604)
        self.assertNotEqual(Analyzer.point_key(native), Analyzer.point_key(independent))
        self.assertGreater(abs(native['coordinate'][2]-independent['coordinate'][2]), 1e-13)
        with self.assertRaisesRegex(ValueError, 'absolute-coordinate safeguard'):
            self.check([native], [independent])

    def test_weight_above_original_tolerance_rejects(self):
        with self.assertRaisesRegex(ValueError, 'independent_complete_dynamic_ray_weights'):
            self.check([self.row(weight=.1)], [self.row(weight=.1+3e-11)])

    def test_changed_rounded_key_still_runs_original_weight_check(self):
        native = self.row()
        independent = self.row(math.nextafter(native['coordinate'][2], -math.inf), weight=.1+3e-11)
        with self.assertRaisesRegex(ValueError, 'independent_complete_dynamic_ray_weights'):
            self.check([native], [independent])

    def test_unchanged_key_preserves_ordinary_near_path(self):
        native, independent = self.row(3.0000000000004), self.row(3.0000000000001)
        self.assertEqual(Analyzer.point_key(native), Analyzer.point_key(independent))
        result = self.check([native], [independent])
        self.assertTrue(result['old_rounded_key_equality'])
        self.assertGreater(result['maximum_raw_coordinate_difference'], 1e-13)
        self.assertFalse(result['additional_absolute_coordinate_safeguard_applicable'])
        self.assertIsNone(result['additional_absolute_coordinate_safeguard'])
        self.assertIn('independent_complete_dynamic_ray_coordinates', self.audit.errors)
        self.assertIn('independent_complete_dynamic_ray_weights', self.audit.errors)

    def test_any_changed_key_applies_safeguard_to_every_pair(self):
        native, independent = self.row(), self.row(math.nextafter(self.row()['coordinate'][2], -math.inf))
        other_native, other_independent = self.row(3.0000000000004), self.row(3.0000000000001)
        self.assertEqual(Analyzer.point_key(other_native), Analyzer.point_key(other_independent))
        with self.assertRaisesRegex(ValueError, 'absolute-coordinate safeguard'):
            self.check([native, other_native], [independent, other_independent])


if __name__ == '__main__':
    unittest.main()
