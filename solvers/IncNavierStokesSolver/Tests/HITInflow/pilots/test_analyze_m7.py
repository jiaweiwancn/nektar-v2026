#!/usr/bin/env python3
"""Small estimator/normalization tests. Execute only on authoritative Linux."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import analyze_m7 as analysis


class M7AnalysisTests(unittest.TestCase):
    def test_one_sided_even_nyquist_odd_last_bin_and_jacobian(self):
        for n in (128, 129):
            with self.subTest(n=n):
                t = np.arange(n)
                if n % 2 == 0:
                    signal = np.sin(2*np.pi*7*t/n) + .25*(-1.0)**t
                    expected = .5+.25**2
                else:
                    signal = np.cos(2*np.pi*(n//2)*t/n)
                    expected = .5
                data = np.stack((12.7+signal, -2+signal), axis=1)
                result = analysis.one_sided_psd(data, .01)
                self.assertAlmostEqual(result["integrated_PSD_m2_s2"], expected, places=12)
                self.assertLess(result["parseval_relative_residual"], 1e-12)
                last_variance = result["Sf_m2_s2_per_Hz"][-1]*result["frequency_resolution_Hz"]
                self.assertAlmostEqual(last_variance, .25**2 if n % 2 == 0 else .5, places=12)
                converted = analysis.longitudinal_spectrum(result, 12.7)
                self.assertAlmostEqual(converted["E11_integral_m2_s2"], expected, places=12)
                self.assertLess(converted["Jacobian_parseval_relative_residual"], 1e-12)

    def test_periodic_hann_welch_window_energy_and_normalization(self):
        t = np.arange(128)
        data = np.stack((12.7+2*np.sin(2*np.pi*8*t/128),
                         3+np.cos(2*np.pi*20*t/128)), axis=1)
        result = analysis.one_sided_psd(data, .01, segment=64, window="hann", overlap=32)
        self.assertEqual(result["segments"], 3)
        self.assertEqual(result["window"], "periodic_Hann")
        self.assertAlmostEqual(result["window_energy_sum_squared_weights"], 24, places=12)
        self.assertAlmostEqual(result["equivalent_noise_bandwidth_Hz"], 1.5/(64*.01), places=12)
        self.assertAlmostEqual(result["integrated_PSD_m2_s2"], 1.25, places=12)

    def test_frozen_advected_pair_delay_sign_overlap_and_two_normalizations(self):
        times = np.arange(601)*.01
        delay = .3
        velocity = np.empty((len(times), 2, 2, 4, 3))
        for y in range(2):
            for z in range(4):
                phase = .2*y+.1*z
                for c in range(3):
                    def wave(t):
                        return np.sin(2*np.pi*3*t+phase+c*.3)+.4*np.cos(2*np.pi*.7*t+.5*phase)
                    velocity[:, 0, y, z, c] = 12+c+y+wave(times)
                    velocity[:, 1, y, z, c] = -4+c+z+2*wave(times-delay)
        selected = (times >= 1) & (times <= 5)
        moment = analysis.moments(velocity[selected])
        mean = moment.mean.reshape(2, 2, 4, 3)
        variance = np.mean((moment.m2/moment.n).reshape(2, 2, 4, 3), axis=(1, 2))
        pair = analysis.delayed_pair(times, velocity, 0, 1, delay, 1, 5, mean, variance, True)
        self.assertEqual(pair["pairs_per_probe"], 371)
        self.assertAlmostEqual(pair["actual_paired_sample_span_s"], 3.7, places=12)
        self.assertEqual(pair["upstream_queries_before_record_start"], 0)
        np.testing.assert_allclose(pair["overlap_Pearson_components"], 1, atol=2e-12, rtol=0)
        extension = analysis.delayed_pair(times, velocity, 0, 1, delay, 1, 5, mean, variance, False)
        self.assertEqual(extension["pairs_per_probe"], 401)
        self.assertEqual(extension["upstream_queries_before_record_start"], 30)
        self.assertTrue(extension["full_downstream_record_covered"])
        np.testing.assert_allclose(extension["overlap_Pearson_components"], 1, atol=2e-12, rtol=0)
        self.assertIn("need not lie", pair["normalizations"])

    def test_sampling_flags_and_exact_constant_signal_zero_handling(self):
        times = np.arange(64)*.01
        data = np.zeros((64, 1, 2, 2, 3))
        signal = np.sin(2*np.pi*np.arange(64)/8)
        signal[32:] *= 3
        for c in range(3):
            data[:, 0, :, :, c] = ((c+1)*signal)[:, None, None]
        result, _ = analysis.station_statistics(times, data, [0.0], 4)
        self.assertEqual(result["qualification"][0]["engineering_screen"],
                         "not_qualified_at_2_percent_sampling_precision")
        self.assertGreater(result["qualification"][0]["full_to_final_half_change_percent"][0], 2)
        self.assertEqual(set(result["blocks"]), {"4", "8", "16"})
        data[:] = 12.7
        result, _ = analysis.station_statistics(times, data, [0.0], 4)
        self.assertEqual(result["planes"][0]["variance_m2_s2"], [0, 0, 0])
        self.assertEqual(result["qualification"][0]["engineering_screen"], "screens_passed_only")

    def test_periodic_zero_mode_is_retained_and_removed_only_in_labelled_diagnostic(self):
        t = np.arange(64)
        amplitude = np.sin(2*np.pi*4*t/64)
        values = np.zeros((64, 2, 32, 3))
        values[..., 0] = amplitude[:, None, None]*(2+np.cos(2*np.pi*np.arange(32)/32))[None, None, :]
        result = analysis.z_correlations(values, .6)
        self.assertAlmostEqual(result["zero_mode_retained"]["finite_half_period_integrals_m"][0], 4*.6/9, places=12)
        self.assertAlmostEqual(result["zero_mode_removed"]["finite_half_period_integrals_m"][0], 0, places=12)
        self.assertEqual(result["integration_limit_m"], .3)
        self.assertEqual(result["zero_mode_retained"]["first_zero_positive_lobe"][0]["status"],
                         "no_zero_within_half_period_truncated")
        self.assertAlmostEqual(result["zero_mode_removed"]["first_zero_positive_lobe"][0]["first_zero_m"], .15, places=12)
        self.assertIn("mathematically zero", result["zero_mode_removed"]["full_half_period_interpretation"])

    def test_first_zero_positive_lobe_linear_crossing_and_truncation(self):
        radius = [0, .1, .2, .3]
        crossing = analysis.first_zero_positive_lobe(radius, [1, .5, -.5, -.7])
        self.assertAlmostEqual(crossing["first_zero_m"], .15, places=12)
        self.assertAlmostEqual(crossing["positive_lobe_integral_m"], .0875, places=12)
        self.assertFalse(crossing["truncated_at_half_period"])
        truncated = analysis.first_zero_positive_lobe(radius, [1, .8, .6, .4])
        self.assertAlmostEqual(truncated["positive_lobe_integral_m"], .21, places=12)
        self.assertTrue(truncated["truncated_at_half_period"])
        self.assertIsNone(truncated["first_zero_m"])
        self.assertIsNone(analysis.first_zero_positive_lobe(radius, [0, 0, 0, 0])["positive_lobe_integral_m"])

    def test_analysis_destinations_are_fresh_and_outside_case_source_and_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, source = root/"case", root/"source"
            case.mkdir()
            (source/"library/SolverUtils").mkdir(parents=True)
            (source/"solvers/IncNavierStokesSolver").mkdir(parents=True)
            manifest, history, receipt = case/"manifest.json", case/"native.his", case/"run.receipt.json"
            for path in (manifest, history, receipt):
                path.write_text("preserve input\n")
            cache, output = root/"cache/history.npz", root/"analysis/result.json"
            analysis.validate_paths(manifest, history, cache, output, source)
            for bad in (history, receipt, case/"new.json", source/"new.json", cache, cache.with_suffix(".npz.json")):
                with self.subTest(output=str(bad)), self.assertRaises(ValueError):
                    analysis.validate_paths(manifest, history, cache, bad, source)
            output.parent.mkdir()
            output.write_text("preserve old analysis\n")
            with self.assertRaisesRegex(ValueError, "fresh"):
                analysis.validate_paths(manifest, history, cache, output, source)
            for bad in (history, case/"cache.npz", source/"cache.npz"):
                with self.subTest(cache=str(bad)), self.assertRaises(ValueError):
                    analysis.validate_paths(manifest, history, bad, root/"fresh.json", source)
            external_history = root/"external-native.his"
            external_history.write_text("external input\n")
            with self.assertRaisesRegex(ValueError, "native history"):
                analysis.validate_paths(manifest, external_history, external_history, root/"fresh.json", source)
            self.assertEqual(history.read_text(), "preserve input\n")
            self.assertEqual(receipt.read_text(), "preserve input\n")
            self.assertEqual(output.read_text(), "preserve old analysis\n")

    def test_native_ascii_cache_identity_and_manifest_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, source = root/"case", root/"source"
            case.mkdir()
            (source/"library/SolverUtils").mkdir(parents=True)
            (source/"solvers/IncNavierStokesSolver").mkdir(parents=True)
            probes = [{"id": 2*y+z, "station_index": 0, "iy": y, "iz": z,
                       "xyz": [0.0, y*.5, z*.5]} for y in range(2) for z in range(2)]
            manifest = {"schema": "NEKTAR_HIT_PILOT_V1", "sampling": {"stations_x_m": [0.0],
                "probe_shape_yz": [2, 2], "probes": probes}, "mesh": {"period_y_m": 1.0, "period_z_m": 1.0},
                "time": {"initial_time_s": 0.0, "end_s": .63, "sample_dt_s": .01, "record_start_s": .32}}
            manifest_path = case/"manifest.json"
            manifest_path.write_text(json.dumps(manifest))
            history = case/"native.his"
            with history.open("w") as stream:
                stream.write("# History data for variables (:u,v,w,p,) at points:\n")
                for p in probes:
                    stream.write("# %d %.19e %.19e %.19e\n" % (p["id"], *p["xyz"]))
                for i in range(64):
                    for p in probes:
                        stream.write("%.19e %.19e 0 0 0\n" % (i*.01, 12.7+i*.001+p["id"]*.1))
            cache = root/"cache/history.npz"
            analysis.create_cache(cache, history, manifest_path, manifest, None, source)
            times, velocity, receipt = analysis.load_cache(cache, manifest_path, manifest)
            self.assertEqual(velocity.shape, (64, 1, 2, 2, 3))
            self.assertEqual(len(times), 64)
            self.assertFalse(velocity.flags.writeable)
            self.assertEqual(receipt["metadata"]["native_history_sha256"], hashlib.sha256(history.read_bytes()).hexdigest())
            manifest_path.write_text(manifest_path.read_text()+"\n")
            with self.assertRaisesRegex(ValueError, "manifest"):
                analysis.load_cache(cache, manifest_path, manifest)


if __name__ == "__main__":
    unittest.main()
