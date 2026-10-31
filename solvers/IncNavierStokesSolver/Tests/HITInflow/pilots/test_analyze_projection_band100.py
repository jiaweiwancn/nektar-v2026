#!/usr/bin/env python3
"""Small temporal-projection tests; execute on authoritative Linux."""

import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

import analyze_projection_band100 as analysis


def metadata(samples=641, counts=(2, 1, 0), dt=.00025):
    return {"ranks": len(counts), "temporal_samples": samples,
        "time_grid": {"mode": "uniform", "regular_samples": samples,
            "evaluation_calls": samples + 1, "start_s": .2, "sample_dt_s": dt,
            "end_s": .2 + (samples - 1)*dt, "repeat_time_s": .2,
            "csv_excludes_repeat_check": True}, "repeat_check_performed": True,
        "repeat_identical": True, "csv_exported_points_per_rank": list(counts),
        "csv_exported_points_per_time": sum(counts), "csv_exported_rows": sum(counts)*samples,
        "native_points_per_time": sum(counts) + 4,
        "projection_max_error": [.5, .5, .5]}


def arrays(samples=641, points=3, dt=.00025):
    times = .2 + np.arange(samples)*dt
    raw = np.zeros((samples, points, 3))
    for p in range(points):
        signal = sum(amplitude*np.sin(2*np.pi*(index/(512*dt))*times + .1*p)
            for index, amplitude in zip((4, 8, 14, 22, 26), (1.0, .8, .6, .4, .2)))
        raw[:, p, 0] = 12.7 + signal
        raw[:, p, 1] = .7*signal
        raw[:, p, 2] = .5*signal
    projected = .9*raw
    projected[:, :, 0] = 12.7 + .9*(raw[:, :, 0] - 12.7)
    return times, raw, projected


def write_csv(path, times, raw, projected, coordinates):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(analysis.HEADER)
        for i, time in enumerate(times):
            for p, xyz in enumerate(coordinates):
                writer.writerow([format(float(time), ".17g"), *xyz,
                    *[format(float(v), ".17g") for v in raw[i, p]],
                    *[format(float(v), ".17g") for v in projected[i, p]]])


class ProjectionBand100Tests(unittest.TestCase):
    def test_factor_point_nine_has_point_eight_one_power_in_all_bands(self):
        _, raw, projected = arrays()
        result = analysis.compare_projection(raw, projected, .00025, 12.7)
        for component in ("u", "v", "w"):
            row = result["spectral_projection_loss"][component]
            for band in row["welch_band_power_loss"]:
                self.assertAlmostEqual(band["projected_over_raw_ratio"], .81, places=12)
                self.assertAlmostEqual(band["signed_loss_percent"], 19, places=11)
            self.assertAlmostEqual(row["welch_F100_loss"]["projected_over_raw_ratio"], .81, places=12)
        self.assertFalse(result["informational_projection_screen"]["within_5_percent_all_welch_bands"])
        np.testing.assert_allclose(result["errors"]["sampled_rms_projection_error_m_s"],
            np.sqrt(np.mean((projected - raw)**2, axis=(0, 1))), rtol=1e-14, atol=0)

    def test_exact_identity_and_zero_power_are_handled_explicitly(self):
        _, raw, _ = arrays()
        result = analysis.compare_projection(raw, raw.copy(), .00025, 12.7)
        self.assertEqual(result["errors"]["sampled_max_abs_projection_error_m_s"], [0, 0, 0])
        self.assertTrue(result["informational_projection_screen"]["within_5_percent_all_welch_bands"])
        constant = np.full((641, 2, 3), 12.7)
        zero = analysis.compare_projection(constant, constant, .00025, 12.7)
        self.assertFalse(zero["informational_projection_screen"]["all_band_relative_losses_defined"])
        self.assertFalse(zero["informational_projection_screen"]["within_5_percent_all_welch_bands"])
        self.assertTrue(zero["spectral_projection_loss"]["u"]["welch_band_power_loss"][0]["both_zero"])

    def test_native_metadata_regular_grid_and_counts_are_required(self):
        original = metadata()
        self.assertEqual(analysis.temporal_grid(original)["exported_points"], 3)
        for section, key, value in ((None, "temporal_samples", 0),
                ("time_grid", "mode", "legacy_fractions"),
                ("time_grid", "end_s", .9),
                (None, "csv_exported_rows", 12),
                (None, "repeat_identical", False),
                (None, "csv_exported_points_per_rank", [1, 1, 0])):
            changed = copy.deepcopy(original)
            (changed if section is None else changed[section])[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                analysis.temporal_grid(changed)

    def test_missing_frame_changed_coordinates_nonfinite_and_extra_rows_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"rank.csv"
            times, raw, projected = arrays(samples=3, points=2)
            xyz = [[0, 0, 0], [0, .1, .2]]
            grid = analysis.temporal_grid(metadata(samples=3, counts=(2,)))
            write_csv(path, times, raw, projected, xyz)
            self.assertEqual(analysis.load_rank_csv(path, 0, grid)["raw"].shape, (3, 2, 3))
            original = path.read_text().splitlines()
            cases = []
            cases.append(original[:-2])
            changed = list(original)
            parts = changed[3].split(",")
            parts[2] = ".001"
            changed[3] = ",".join(parts)
            cases.append(changed)
            changed = list(original)
            parts = changed[3].split(",")
            parts[0] = str(times[1] + .00001)
            changed[3] = ",".join(parts)
            cases.append(changed)
            changed = list(original)
            parts = changed[1].split(",")
            parts[4] = "nan"
            changed[1] = ",".join(parts)
            cases.append(changed)
            cases.append(original + [original[-1]])
            for changed in cases:
                path.write_text("\n".join(changed) + "\n")
                with self.assertRaises(ValueError):
                    analysis.load_rank_csv(path, 0, grid)

    def test_complete_native_export_receipt_and_empty_rank(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            case, source = root/"case", root/"source"
            case.mkdir()
            (source/"library/SolverUtils").mkdir(parents=True)
            (source/"solvers/IncNavierStokesSolver").mkdir(parents=True)
            stem = "projection-temporal"
            summary = case/(stem + ".probe.json")
            summary.write_text(json.dumps(metadata()))
            (case/(stem + ".xml")).write_text("<NEKTAR><CONDITIONS><PARAMETERS>"
                "<P>Ubar = 12.7</P></PARAMETERS></CONDITIONS></NEKTAR>")
            times, raw, projected = arrays()
            write_csv(case/(stem + ".probe-rank0.csv"), times, raw[:, :2], projected[:, :2],
                      [[0, 0, 0], [0, .1, .2]])
            write_csv(case/(stem + ".probe-rank1.csv"), times, raw[:, 2:], projected[:, 2:],
                      [[0, .2, .3]])
            empty = case/(stem + ".probe-rank2.csv")
            empty.write_text(",".join(analysis.HEADER) + "\n")
            before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in case.iterdir()}
            args = argparse.Namespace(case=case, stem=stem, output=root/"analysis", source_root=source,
                probe_batch=4, welch_segment_samples=512, no_plot=True)
            result = analysis.analyze(args)
            self.assertEqual(result["sampled_points"]["empty_ranks"], [2])
            self.assertEqual(result["sampled_points"]["count"], 3)
            receipt = analysis.band100.read_json(root/"analysis/receipt.json")
            self.assertEqual(receipt["status"], "completed")
            self.assertEqual(receipt["analysis_sha256"], analysis.m7.sha(root/"analysis/analysis.json"))
            self.assertEqual({p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in case.iterdir()}, before)
            with self.assertRaisesRegex(ValueError, "fresh"):
                analysis.analyze(args)
            args.output = case/"inside-case"
            with self.assertRaisesRegex(ValueError, "outside"):
                analysis.analyze(args)


if __name__ == "__main__":
    unittest.main()
