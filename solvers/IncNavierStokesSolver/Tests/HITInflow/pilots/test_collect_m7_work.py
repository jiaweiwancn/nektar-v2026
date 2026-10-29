#!/usr/bin/env python3
"""Focused native-work arithmetic/preservation tests; execute only on Linux."""

import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import collect_m7_work as work


def linear_field_fixture():
    """Exact slab moments of u=U+a*x, v=b*x, w=c*x over [0,9]."""
    speed, viscosity, length, area = 12.7, 1.5e-5, 9.0, .36
    a, b, c = .01, .02, -.015
    slopes_squared = a*a+b*b+c*c
    molecular = viscosity*(2*a*a+b*b+c*c)
    rows = []
    for index in range(90):
        lower, upper = index/10, (index+1)/10
        center = (lower+upper)/2
        mean_x_squared = (lower*lower+lower*upper+upper*upper)/3
        rows.append({"xmin": lower, "xmax": upper, "xmid": center, "volume": (upper-lower)*area,
            "mean_u": speed+a*center, "mean_v": b*center, "mean_w": c*center,
            "instantaneous_K_about_prescribed_mean": .5*slopes_squared*mean_x_squared,
            "molecular_strain_dissipation": molecular, "svv_xy": 0.0, "svv_z": 0.0,
            "divergence_rms": abs(a)})
    volume = area*length
    flow = {"volume": volume, "viscosity": viscosity,
        "instantaneous_resolved_strain_dissipation": molecular,
        "instantaneous_svv_work_xy": 0.0, "instantaneous_svv_work_z": 0.0, "divergence_rms": abs(a),
        "velocity_l2": math.sqrt(volume*(speed*speed+speed*a*length+slopes_squared*length*length/3))}
    return rows, flow, speed, viscosity


class NativeWorkTests(unittest.TestCase):
    def test_90_slab_closure_against_analytic_linear_field(self):
        rows, flow, speed, viscosity = linear_field_fixture()
        self.assertEqual(len(rows), 90)
        checks = work.native_closure(rows, flow, speed, viscosity)
        for key in ("volume", "molecular_strain_dissipation", "divergence_rms", "velocity_l2"):
            self.assertAlmostEqual(checks[key]["from_slabs"], checks[key]["native_global"], places=12)

    def test_work_and_volume_mismatch_rejected(self):
        rows, flow, speed, viscosity = linear_field_fixture()
        for key, factor, message in (("volume", 1.001, "volumes"),
                ("instantaneous_resolved_strain_dissipation", 1.1, "work")):
            with self.subTest(key=key):
                changed = dict(flow)
                changed[key] *= factor
                with self.assertRaisesRegex(ValueError, message):
                    work.native_closure(rows, changed, speed, viscosity)

    def test_first_center_edge_and_interior_interpolation(self):
        rows = [work.add_work({"xmin": 0., "xmax": 1., "xmid": .5,
                    "molecular_strain_dissipation": 2., "svv_xy": .1, "svv_z": .1, "divergence_rms": .2}),
                work.add_work({"xmin": 1., "xmax": 2., "xmid": 1.5,
                    "molecular_strain_dissipation": 4., "svv_xy": .8, "svv_z": .4, "divergence_rms": .4})]
        for x, center, edge in ((0., .5, True), (.5, .5, False), (2., 1.5, True)):
            with self.subTest(x=x):
                values, metadata = work.interpolate(rows, x)
                self.assertEqual(metadata["source_slab_centers_m"], [center])
                self.assertFalse(metadata["exact_inlet_value"])
                self.assertIn("nearest" if edge else "matching", metadata["method"])
                self.assertEqual(values["molecular_strain_dissipation"], 2. if center == .5 else 4.)
        values, metadata = work.interpolate(rows, 1.)
        self.assertEqual(metadata["fraction"], .5)
        self.assertAlmostEqual(values["total_svv_over_molecular"], .7/3)
        self.assertNotAlmostEqual(values["total_svv_over_molecular"], .2)  # Ratio after work interpolation.
        with self.assertRaisesRegex(ValueError, "outside"):
            work.interpolate(rows, -.1)

    def test_scale_closed_form_and_undefined_zero_work(self):
        result = work.scales(.04, .06, 1e-5)
        self.assertAlmostEqual(result["lambda_surrogate_m"], .01, places=14)
        self.assertAlmostEqual(result["Re_lambda_surrogate"], 200., places=11)
        self.assertAlmostEqual(result["eta_surrogate_m"]**4*.06, (1e-5)**3, delta=1e-27)
        missing = work.scales(.04, 0., 1e-5)
        self.assertIsNone(missing["lambda_surrogate_m"])
        self.assertIsNone(missing["eta_surrogate_m"])
        self.assertIsNone(missing["Re_lambda_surrogate"])

    def test_decay_secants_preserve_decay_and_growth_sign(self):
        planes = [{"station_x_m": x, "variance_m2_s2": [variance]*3, "K_m2_s2": 1.5*variance}
                  for x, variance in ((0., 2.), (2., 1.), (4., 1.5))]
        intervals = work.decay_intervals(planes, 2.)
        for key in ("minus_U_delta_K_over_delta_x_m2_s3", "minus_1p5_U_delta_varu_over_delta_x_m2_s3"):
            self.assertEqual([row[key] for row in intervals], [1.5, -.75])
        with self.assertRaisesRegex(ValueError, "unique"):
            work.decay_intervals([planes[0], planes[0]], 2.)

    def test_output_preservation_before_analysis_read(self):
        with tempfile.TemporaryDirectory(prefix="m7-work-preservation-") as directory:
            root = Path(directory)
            case, source = root/"case", root/"source"
            case.mkdir()
            source.mkdir()
            existing = case/"native.his"
            existing.write_text("preserve native history\n")
            with self.assertRaisesRegex(ValueError, "fresh"):
                work.collect(SimpleNamespace(output=existing, analysis=root/"missing-analysis.json"))
            self.assertEqual(existing.read_text(), "preserve native history\n")
            for output in (case/"new.json", source/"new.json"):
                with self.subTest(output=output), self.assertRaisesRegex(ValueError, "outside"):
                    work.validate_output(output, case, source)
            self.assertEqual(work.validate_output(root/"work.json", case, source), (root/"work.json").resolve())

    def test_external_copy_requires_explicit_source_root(self):
        with tempfile.TemporaryDirectory(prefix="m7-work-copy-") as directory:
            output = Path(directory)/"work.json"
            with self.assertRaisesRegex(ValueError, "--source-root"):
                work.validate_output(output, source_root=None)


if __name__ == "__main__":
    unittest.main()
