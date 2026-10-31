#!/usr/bin/env python3
"""Small deterministic tests; execute on the authoritative Linux host."""

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest

import numpy as np

import analyze_band100 as analysis


def comparison_fixture():
    manifest = {
        "physics": {"mean_velocity_m_s": 12.7, "kinvis_m2_s": 1.5e-5,
            "normalization": "Preserve", "interpolation": "LogLinear",
            "low_extension": "K4", "high_extension": "Gaussian"},
        "synthesis": {"seed_uint64": "20261028", "counts": [2548, 154, 154],
            "periods_m": [10, .6, .6]},
        "sampling": {"stations_x_m": [0.0], "probe_shape_yz": [2, 2],
            "probes": [{"xyz": [0.0, y*.3, z*.3]} for y in range(2) for z in range(2)]},
        "mesh": {"period_y_m": .6, "period_z_m": .6, "hom_modes_z": 192},
        "knobs": {"length": 1.0}, "initial_conditions": {"kind": "mean_only"}}
    window = {"requested_start_s": .2, "requested_end_s": .6,
        "first_last_samples_s": [.2, .6], "samples": 1601, "sample_dt_s": .00025}
    bands = [{"k1_limits_rad_m": list(band), "power_m2_s2": i + 1.0}
             for i, band in enumerate(analysis.BANDS)]
    return {"schema": analysis.SCHEMA,
        "control_identity": analysis.control_identity(manifest, window, ["spectrum hash"]),
        "synthesis_band": {"cutoff_rad_m": 600, "taper_start_rad_m": 480},
        "numerical_settings": {"mesh": manifest["mesh"], "time_step_s": 1.25e-5,
            "svv": {"model_xy": "dg", "model_z": "ExpKernel"}},
        "statistics": {"planes": [{"station_x_m": 0.0,
            "variance_m2_s2": [1.0, 2.0, 3.0], "K_m2_s2": 3.0}]},
        "spectra": [{"station_x_m": 0.0, "components": {
            c: {"welch": {"band_powers": copy.deepcopy(bands), "F_component_100_m3_s2": .01}}
            for c in analysis.COMPONENTS}}]}


class Band100Tests(unittest.TestCase):
    def test_integral_inserts_endpoints_on_nonuniform_grid(self):
        # F(k)=2k+1, integrated from1 to10, independent of bin alignment.
        self.assertAlmostEqual(analysis.integrate_density([0, 3, 8, 15],
            [1, 7, 17, 31], 1, 10), 108.0, places=12)
        self.assertEqual(analysis.integrate_density([0, 3, 8, 15], [0]*4, 1, 10), 0)

    def test_invalid_density_or_outside_band_is_rejected(self):
        for x, y, low, high in (([0, 1], [1, 1], -.1, 1),
                ([0, 1], [1, 1], 0, 2), ([0, 1], [1, -1], 0, 1),
                ([0, 0], [1, 1], 0, 1), ([0, 1], [1, math.nan], 0, 1),
                ([0, 1], [1], 0, 1)):
            with self.subTest(x=x, y=y), self.assertRaises(ValueError):
                analysis.integrate_density(x, y, low, high)

    def test_three_component_power_and_parseval_on_sparse_grid(self):
        n, dt = 1024, .00025
        t = np.arange(n)*dt
        values = np.zeros((n, 2, 2, 3))
        for y in range(2):
            for z in range(2):
                phase = .1*y+.2*z
                values[:, y, z, 0] = 12.7 + 2*np.sin(2*np.pi*125*t+phase)
                values[:, y, z, 1] = 3*np.sin(2*np.pi*62.5*t+phase)
                values[:, y, z, 2] = .5*np.cos(2*np.pi*250*t+phase)
        result = analysis.component_spectra(values, dt, 12.7)
        for c, expected in zip(analysis.COMPONENTS, (2.0, 4.5, .125)):
            with self.subTest(component=c):
                for kind in ("boxcar", "welch"):
                    self.assertAlmostEqual(result[c][kind]["F_component_integral_m2_s2"],
                                           expected, places=12)
                    self.assertLess(result[c][kind]["parseval_relative_residual"], 1e-12)
                    self.assertNotIn("E11_m3_s2", result[c][kind])
                self.assertEqual(len(result[c]["welch"]["band_powers"]), 4)
        self.assertEqual(result["u"]["welch"]["component_density_name"], "F11")
        self.assertEqual(result["v"]["welch"]["component_density_name"], "F22")
        self.assertEqual(result["w"]["welch"]["component_density_name"], "F33")
        self.assertAlmostEqual(result["u"]["welch"]["band_powers"][2]["power_m2_s2"], 2, places=12)

    def test_matched_numerics_and_bandwidth_are_labelled_separately(self):
        reference = comparison_fixture()
        current = copy.deepcopy(reference)
        current["numerical_settings"]["time_step_s"] *= .5
        current["statistics"]["planes"][0]["variance_m2_s2"][0] *= 1.01
        current["spectra"][0]["components"]["u"]["welch"]["band_powers"][0]["power_m2_s2"] *= 1.03
        result = analysis.compare_analyses(current, reference)
        self.assertTrue(result["matched_controls_eligible"])
        self.assertTrue(result["eligible_for_numerical_refinement_screen"])
        self.assertEqual(result["kind"], "matched_numerical_control")
        self.assertAlmostEqual(result["maximum_absolute_variance_K_change_percent"], 1, places=12)
        self.assertAlmostEqual(result["maximum_absolute_band_power_change_percent"], 3, places=12)
        current["synthesis_band"] = {"cutoff_rad_m": 800, "taper_start_rad_m": 640}
        result = analysis.compare_analyses(current, reference)
        self.assertTrue(result["matched_controls_eligible"])
        self.assertTrue(result["bandwidth_changed"])
        self.assertFalse(result["eligible_for_numerical_refinement_screen"])
        self.assertEqual(result["kind"], "input_bandwidth_model_sensitivity")

    def test_control_identity_mismatches_block_comparison(self):
        reference = comparison_fixture()
        changes = {
            "counts": [3000, 154, 154], "seed_uint64": "42", "periods_m": [32, .6, .6],
            "stations_x_m": [.1], "probe_shape_yz": [4, 1],
            "probe_coordinates_m": [[0, 0, 0]], "spectrum_sha256": ["different spectrum"],
            "domain_lengths_m": [9, .6, .6]}
        for key, value in changes.items():
            current = copy.deepcopy(reference)
            current["control_identity"][key] = value
            result = analysis.compare_analyses(current, reference)
            self.assertFalse(result["matched_controls_eligible"], key)
            self.assertEqual(result["stations"], [], key)
        for key, value in (("sample_dt_s", .0005), ("samples", 801), ("requested_end_s", .5)):
            current = copy.deepcopy(reference)
            current["control_identity"]["window"][key] = value
            self.assertFalse(analysis.compare_analyses(current, reference)["matched_controls_eligible"])
        current = copy.deepcopy(reference)
        current["control_identity"]["window"]["first_last_samples_s"][1] += 1e-14
        self.assertTrue(analysis.compare_analyses(current, reference)["matched_controls_eligible"])

    def test_zero_reference_and_malformed_station_counts(self):
        reference = comparison_fixture()
        current = copy.deepcopy(reference)
        reference["statistics"]["planes"][0]["variance_m2_s2"][0] = 0
        result = analysis.compare_analyses(current, reference)
        self.assertIsNone(result["maximum_absolute_variance_K_change_percent"])
        self.assertTrue(result["stations"][0]["variances_and_K"]["u_variance"]["zero_reference"])
        current["statistics"]["planes"] = []
        with self.assertRaisesRegex(ValueError, "station count"):
            analysis.compare_analyses(current, reference)

    def test_output_is_fresh_and_outside_case_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, case = root/"source", root/"case"
            source.mkdir()
            case.mkdir()
            manifest = case/"manifest.json"
            manifest.write_text("{}")
            for target in (source/"new", case/"new", case, source):
                with self.subTest(target=target), self.assertRaises(ValueError):
                    analysis.validate_output(manifest, target, source)
            self.assertEqual(analysis.validate_output(manifest, root/"fresh", source), root/"fresh")

    def test_end_to_end_sparse_temporal_analysis_and_immutable_receipt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, case = root/"source", root/"case"
            (source/"library/SolverUtils").mkdir(parents=True)
            (source/"solvers/IncNavierStokesSolver").mkdir(parents=True)
            case.mkdir()
            probes = [{"id": 2*y+z, "station_index": 0, "iy": y, "iz": z,
                       "xyz": [0.0, y*.3, z*.3]} for y in range(2) for z in range(2)]
            spectrum = case/"spectrum.dat"
            spectrum.write_text("fixture bound by hash; no CFD or spectrum fit\n")
            session = case/"pilot.xml"
            session.write_text("<NEKTAR/>\n")
            base = comparison_fixture()
            manifest = {"schema": "NEKTAR_HIT_PILOT_V1", "case": session.name,
                "history_file": "native.his", "knobs": {"length": 1.0},
                "physics": base["control_identity"]["physics"],
                "mesh": base["numerical_settings"]["mesh"],
                "svv": base["numerical_settings"]["svv"],
                "initial_conditions": {"kind": "mean_only"},
                "sampling": {"stations_x_m": [0.0], "probe_shape_yz": [2, 2], "probes": probes},
                "time": {"dt_s": 1.25e-5, "initial_time_s": 0.0, "record_start_s": 0.0,
                         "end_s": .16, "sample_dt_s": .00025},
                "synthesis": {"seed_uint64": "20261028", "counts": [2548, 154, 154],
                    "periods_m": [10, .6, .6], "cutoff_rad_m": 600, "taper_start_rad_m": 480},
                "sha256": {spectrum.name: hashlib.sha256(spectrum.read_bytes()).hexdigest()}}
            manifest_path = case/"manifest.json"
            manifest_path.write_text(json.dumps(manifest))
            history = case/"native.his"
            with history.open("w") as stream:
                stream.write("# History data for variables (:u,v,w,p,) at points:\n")
                for p in probes:
                    stream.write("# %d %.19e %.19e %.19e\n" % (p["id"], *p["xyz"]))
                for i in range(641):
                    time = i*.00025
                    for p in probes:
                        wave = math.sin(2*math.pi*125*time+.1*p["id"])
                        stream.write("%.19e %.19e %.19e %.19e 0\n" %
                                     (time, 12.7+wave, .5*wave, 2*wave))
            experimental = root/"experimental.json"
            experimental.write_text(json.dumps({"schema": "nektar.hit.cbc1971.source.v1",
                "table2a": {"rows": [[".05", "570", None, None], ["1.0", "114", None, None]]}}))
            original = hashlib.sha256(history.read_bytes()).hexdigest()
            args = argparse.Namespace(manifest=manifest_path, output=root/"analysis", source_root=source,
                reference_analysis=None, experimental_source=experimental, start=None, end=None,
                probe_batch=4, welch_segment_samples=512, no_plot=True)
            result = analysis.analyze(args)
            self.assertEqual(result["sampling"]["probe_shape_yz"], [2, 2])
            self.assertEqual(result["numerical_settings"]["mesh"]["hom_modes_z"], 192)
            self.assertEqual(result["window"]["samples"], 641)
            receipt = analysis.read_json(root/"analysis/receipt.json")
            self.assertEqual(receipt["status"], "completed")
            self.assertEqual(receipt["analysis_sha256"], analysis.m7.sha(root/"analysis/analysis.json"))
            self.assertEqual(hashlib.sha256(history.read_bytes()).hexdigest(), original)
            args.reference_analysis = args.output
            args.output = root/"matched"
            compared = analysis.analyze(args)
            self.assertTrue(compared["reference_comparison"]["matched_controls_eligible"])
            self.assertEqual(compared["reference_comparison"]["maximum_absolute_band_power_change_percent"], 0)
            with self.assertRaisesRegex(ValueError, "fresh"):
                analysis.analyze(args)


if __name__ == "__main__":
    unittest.main()
