#!/usr/bin/env python3
"""Measure sampled temporal inlet-projection loss through100 rad/m on Linux.

Inputs are the uniform-time CSV/JSON exports of HITBoundaryProbe. This tests
the boundary projection only, before CFD. Exported points are neither an exact
plane quadrature nor independent realizations; shared-coordinate duplicates
are retained with their native export multiplicities.
"""

import argparse
import csv
import datetime
import math
import os
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True

import numpy as np

import analyze_band100 as band100
import analyze_m7 as m7


SCHEMA = "NEKTAR_HIT_BAND100_PROJECTION_ANALYSIS_V1"
RECEIPT_SCHEMA = "NEKTAR_HIT_BAND100_PROJECTION_RECEIPT_V1"
HEADER = ["time", "x", "y", "z", "u", "v", "w", "projected_u", "projected_v", "projected_w"]
require = band100.require


def integer(value, label, minimum=0):
    require(isinstance(value, int) and not isinstance(value, bool) and value >= minimum,
            label + " must be an integer >= " + str(minimum))
    return value


def temporal_grid(metadata):
    ranks = integer(metadata["ranks"], "ranks", 1)
    samples = integer(metadata["temporal_samples"], "temporal_samples", 2)
    grid = metadata["time_grid"]
    require(grid["mode"] == "uniform", "projection analysis requires uniform temporal exports")
    require(integer(grid["regular_samples"], "regular_samples", 2) == samples,
            "native temporal sample counts disagree")
    require(integer(grid["evaluation_calls"], "evaluation_calls", samples) == samples + 1,
            "native temporal export needs one separate repeat evaluation")
    start, dt, end = (float(grid[key]) for key in ("start_s", "sample_dt_s", "end_s"))
    require(all(math.isfinite(v) for v in (start, dt, end)) and dt > 0,
            "native temporal grid must be finite with positive dt")
    expected_end = start + (samples - 1)*dt
    tolerance = 64*sys.float_info.epsilon*max(abs(start), abs(end), dt) + 2e-12*dt
    require(abs(end - expected_end) <= tolerance, "native temporal end does not match its count/dt")
    require(metadata.get("repeat_check_performed") is True
            and metadata.get("repeat_identical") is True
            and grid.get("csv_excludes_repeat_check") is True,
            "native temporal export must have a successful separate repeat check")
    repeat_time = float(grid["repeat_time_s"])
    require(math.isfinite(repeat_time) and start - tolerance <= repeat_time <= end + tolerance,
            "native repeat check must use a time inside the exported interval")
    counts = metadata["csv_exported_points_per_rank"]
    require(isinstance(counts, list) and len(counts) == ranks,
            "native per-rank exported point counts disagree with ranks")
    counts = [integer(value, "rank exported point count") for value in counts]
    total = integer(metadata["csv_exported_points_per_time"], "exported points per time", 1)
    require(sum(counts) == total, "native global/per-rank exported point counts disagree")
    require(integer(metadata["csv_exported_rows"], "exported rows", 1) == samples*total,
            "native exported row count disagrees with its temporal layout")
    native = integer(metadata["native_points_per_time"], "native points per time", total)
    return {"samples": samples, "ranks": ranks, "start_s": start, "dt_s": dt,
            "end_s": end, "tolerance_s": tolerance, "points_per_rank": counts,
            "exported_points": total, "native_points": native}


def load_rank_csv(path, rank, grid):
    """Read exact ordered point frames and reject missing or unstable frames."""
    expected_points = grid["points_per_rank"][rank]
    times, coordinates, raw_frames, projected_frames = [], None, [], []
    frame_time, frame_coords, frame_raw, frame_projected = None, [], [], []

    def finish():
        nonlocal coordinates, frame_coords, frame_raw, frame_projected
        require(len(frame_coords) == expected_points,
                f"rank{rank}: incomplete/extra point frame")
        actual = np.asarray(frame_coords, dtype=np.float64)
        if coordinates is None:
            coordinates = actual
        else:
            require(np.array_equal(actual, coordinates),
                    f"rank{rank}: ordered native coordinates changed across frames")
        index = len(times)
        require(index < grid["samples"], f"rank{rank}: extra temporal frame")
        expected_time = grid["start_s"] + index*grid["dt_s"]
        require(abs(frame_time - expected_time) <= grid["tolerance_s"],
                f"rank{rank}: nonuniform/incomplete native time grid")
        times.append(frame_time)
        raw_frames.append(frame_raw)
        projected_frames.append(frame_projected)
        frame_coords, frame_raw, frame_projected = [], [], []

    with path.open("r", encoding="ascii", newline="") as source:
        reader = csv.reader(source)
        require(next(reader, None) == HEADER, f"rank{rank}: unexpected projection CSV header")
        for line_number, row in enumerate(reader, 2):
            require(expected_points > 0, f"rank{rank}: data present for declared empty rank")
            require(len(row) == len(HEADER), f"rank{rank} line{line_number}: expected ten columns")
            try:
                values = [float(value) for value in row]
            except ValueError as error:
                raise ValueError(f"rank{rank} line{line_number}: nonnumeric CSV value") from error
            require(all(math.isfinite(value) for value in values),
                    f"rank{rank} line{line_number}: nonfinite CSV value")
            time = values[0]
            if frame_time is None:
                frame_time = time
            elif time != frame_time:
                finish()
                require(time > frame_time, f"rank{rank}: reversed/duplicate native time frame")
                frame_time = time
            frame_coords.append(values[1:4])
            frame_raw.append(values[4:7])
            frame_projected.append(values[7:10])
            require(len(frame_coords) <= expected_points, f"rank{rank}: too many rows in a native frame")
    if expected_points == 0:
        return None
    require(frame_time is not None, f"rank{rank}: missing temporal rows")
    finish()
    require(len(times) == grid["samples"], f"rank{rank}: incomplete temporal export")
    return {"rank": rank, "times": np.asarray(times, dtype=np.float64),
            "coordinates": coordinates, "raw": np.asarray(raw_frames, dtype=np.float64),
            "projected": np.asarray(projected_frames, dtype=np.float64)}


def loss(raw, projected):
    raw, projected = float(raw), float(projected)
    require(math.isfinite(raw) and math.isfinite(projected) and min(raw, projected) >= 0,
            "spectral powers must be finite/nonnegative")
    ratio = projected/raw if raw else None
    require(ratio is None or math.isfinite(ratio), "projection power ratio overflow")
    return {"raw": raw, "projected": projected, "projected_over_raw_ratio": ratio,
            "signed_loss_percent": 100*(1-ratio) if ratio is not None else None,
            "zero_raw_reference": raw == 0, "both_zero": raw == projected == 0,
            "sign_convention": "Positive loss is attenuation; negative loss is amplification."}


def compare_projection(raw, projected, dt, speed, segment=512, batch=256):
    require(raw.shape == projected.shape and raw.ndim == 3 and raw.shape[2] == 3
            and raw.shape[1] > 0 and np.isfinite(raw).all() and np.isfinite(projected).all(),
            "projection arrays must be finite matching time,point,component data")
    spectra = {"raw": band100.component_spectra(raw[:, :, None, :], dt, speed, segment, batch),
               "projected": band100.component_spectra(projected[:, :, None, :], dt, speed, segment, batch)}
    errors = projected - raw
    _, raw_fluctuations = m7.centered(raw)
    _, projected_fluctuations = m7.centered(projected)
    raw_variance = np.mean(raw_fluctuations*raw_fluctuations, axis=(0, 1))
    projected_variance = np.mean(projected_fluctuations*projected_fluctuations, axis=(0, 1))
    rms_error = np.sqrt(np.mean(errors*errors, axis=(0, 1)))
    raw_rms = np.sqrt(raw_variance)
    comparisons, changes = {}, []
    for c, component in enumerate(band100.COMPONENTS):
        own, mapped = spectra["raw"][component], spectra["projected"][component]
        for kind in ("boxcar", "welch"):
            if kind == "boxcar":
                for values, variance in ((own, raw_variance[c]), (mapped, projected_variance[c])):
                    integral = values[kind]["F_component_integral_m2_s2"]
                    error = abs(integral - variance)/variance if variance else abs(integral)
                    require(error <= 5e-10, "projection boxcar PSD does not recover sampled variance")
        bands = []
        for expected, a, b in zip(band100.BANDS, own["welch"]["band_powers"],
                                   mapped["welch"]["band_powers"]):
            require(a["k1_limits_rad_m"] == b["k1_limits_rad_m"] == list(expected),
                    "projection physical bands disagree")
            difference = loss(a["power_m2_s2"], b["power_m2_s2"])
            changes.append(difference["signed_loss_percent"])
            bands.append({"k1_limits_rad_m": list(expected), **difference})
        comparisons[component] = {"welch_band_power_loss": bands,
            "welch_F100_loss": loss(own["welch"]["F_component_100_m3_s2"],
                                      mapped["welch"]["F_component_100_m3_s2"])}
    defined = all(value is not None for value in changes)
    maximum = max((abs(value) for value in changes if value is not None), default=None)
    return {"spectra": spectra, "spectral_projection_loss": comparisons,
        "errors": {"sampled_rms_projection_error_m_s": rms_error.tolist(),
            "sampled_max_abs_projection_error_m_s": np.max(np.abs(errors), axis=(0, 1)).tolist(),
            "sampled_mean_projection_error_m_s": np.mean(errors, axis=(0, 1)).tolist(),
            "raw_temporal_rms_m_s": raw_rms.tolist(),
            "projected_temporal_rms_m_s": np.sqrt(projected_variance).tolist(),
            "rms_error_percent_of_raw_temporal_rms": [float(100*error/rms) if rms else
                (0.0 if error == 0 else None) for error, rms in zip(rms_error, raw_rms)]},
        "informational_projection_screen": {"threshold_percent": 5.0,
            "all_band_relative_losses_defined": defined,
            "maximum_absolute_welch_band_loss_percent": maximum,
            "within_5_percent_all_welch_bands": defined and maximum is not None and maximum <= 5,
            "scope": "Boundary projection at exported probes only, all3components and four10–100rad/m bands. This is not a CFD-convergence, experimental-fidelity, exact-plane or statistical-confidence verdict. Endpoint100 is reported separately."}}


def mean_speed(session):
    parameters = ET.parse(session).findall("./CONDITIONS/PARAMETERS/P")
    found = [node.text.split("=", 1)[1].strip() for node in parameters
             if node.text and "=" in node.text and node.text.split("=", 1)[0].strip().lower() == "ubar"]
    require(len(found) == 1, "projection session must declare one numeric Ubar parameter")
    try:
        value = float(found[0])
    except ValueError as error:
        raise ValueError("projection Ubar must be a numeric SI speed") from error
    require(math.isfinite(value) and value > 0, "projection Ubar must be finite and positive")
    return value


def plot_result(result, output):
    try:
        import matplotlib
    except ImportError:
        return {"available": False, "reason": "Matplotlib is unavailable"}
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    raw = result["spectra"]["raw"]["u"]["welch"]
    projected = result["spectra"]["projected"]["u"]["welch"]
    x, a, b = (np.asarray(value) for value in (raw["k1_rad_m"], raw["F_component_m3_s2"],
                                             projected["F_component_m3_s2"]))
    keep = (x > 0) & (x <= 110) & (a > 0) & (b > 0)
    if np.any(keep):
        axes[0].loglog(x[keep], a[keep], label="Raw native boundary values")
        axes[0].loglog(x[keep], b[keep], "--", label="Constrained projection")
        axes[1].plot(x[keep], b[keep]/a[keep], label="Projected/raw power density")
    axes[0].set(xlabel="Longitudinal k1 (rad/m)", ylabel="One-sided F11 (m³/s²)", xlim=(3, 110))
    axes[1].axhline(1, color="gray", linewidth=.8)
    axes[1].axhspan(.95, 1.05, color="gray", alpha=.15)
    axes[1].set(xlabel="Longitudinal k1 (rad/m)", ylabel="Projected/raw F11", xlim=(0, 110))
    for axis in axes:
        axis.axvline(100, color="gray", linestyle=":", linewidth=.8)
        axis.grid(True, which="both", alpha=.2)
        axis.legend(fontsize=7)
    fig.suptitle("CBC100: sampled temporal boundary projection")
    fig.text(.01, .015, "Projection only; sampled native points are neither exact plane quadrature nor independent realizations.", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .93))
    path = output / "projection-band100.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return {"available": True, "path": path.name, "sha256": m7.sha(path)}


def analyze(args):
    case, stem = args.case.resolve(), args.stem
    require(case.is_dir() and stem and Path(stem).name == stem and stem not in (".", ".."),
            "case must be a directory and stem a plain basename")
    summary, session = case/(stem + ".probe.json"), case/(stem + ".xml")
    metadata_hash = m7.sha(summary)
    metadata = band100.read_json(summary)
    grid = temporal_grid(metadata)
    require(args.probe_batch > 0 and args.welch_segment_samples >= 4
            and args.welch_segment_samples % 2 == 0,
            "require positive probe batch and an even Welch segment>=4")
    expected_csv = [case/(stem + f".probe-rank{rank}.csv") for rank in range(grid["ranks"])]
    actual_csv = sorted(case.glob(stem + ".probe-rank*.csv"))
    require(set(actual_csv) == set(expected_csv), "native rank CSV set is missing or contains unexpected files")
    source = m7.source_root(args.source_root)
    protected = [summary, session, *expected_csv, Path(__file__).resolve(), Path(band100.__file__).resolve(),
                 Path(m7.__file__).resolve(), Path(m7.ah.__file__).resolve()]
    if (case/"manifest.json").is_file():
        protected.append(case/"manifest.json")
    output = band100.validate_output(summary, args.output, source, protected)
    # Snapshot before parsing to bind the complete export, not only selected rows.
    hashes = {str(path): m7.sha(path) for path in protected}
    require(hashes[str(summary)] == metadata_hash, "native summary changed while inputs were inspected")
    speed = mean_speed(session)
    output.mkdir(parents=True)
    receipt_path = output/"receipt.json"
    receipt = {"schema": RECEIPT_SCHEMA, "status": "running", "command": sys.argv,
        "case": str(case), "stem": stem, "output": str(output),
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat()}
    band100.write_json(receipt_path, receipt, atomic=True)
    try:
        ranks = [load_rank_csv(path, rank, grid) for rank, path in enumerate(expected_csv)]
        selected = [rank for rank in ranks if rank is not None]
        require(selected, "native projection export has no sampled points")
        times = selected[0]["times"]
        require(all(np.array_equal(rank["times"], times) for rank in selected),
                "native rank temporal grids disagree")
        raw = np.concatenate([rank["raw"] for rank in selected], axis=1)
        projected = np.concatenate([rank["projected"] for rank in selected], axis=1)
        require(raw.shape[1] == grid["exported_points"], "assembled projection point count disagrees")
        result = compare_projection(raw, projected, grid["dt_s"], speed,
                                    args.welch_segment_samples, args.probe_batch)
        native_maximum = metadata["projection_max_error"]
        sampled_maximum = result["errors"]["sampled_max_abs_projection_error_m_s"]
        require(isinstance(native_maximum, list) and len(native_maximum) == 3
                and all(isinstance(v, (int, float)) and math.isfinite(v) and v >= 0 for v in native_maximum),
                "native projection maximum must contain three finite nonnegative components")
        require(all(sample <= native + 2e-12*max(1.0, sample, native)
                    for sample, native in zip(sampled_maximum, native_maximum)),
                "sampled projection error exceeds the native all-point maximum")
        result.update(schema=SCHEMA, verification_status="ANALYZED_TEMPORAL_BOUNDARY_PROJECTION",
            case_directory=str(case), native_stem=stem, command=sys.argv,
            window={"samples": len(times), "sample_dt_s": grid["dt_s"],
                "first_last_samples_s": [float(times[0]), float(times[-1])], "mean_speed_m_s": speed},
            sampled_points={"count": grid["exported_points"], "native_points_per_time": grid["native_points"],
                "per_rank": [{"rank": rank["rank"], "coordinates_m": rank["coordinates"].tolist()}
                             for rank in selected], "empty_ranks": [i for i, rank in enumerate(ranks) if rank is None],
                "interpretation": "Native exported index subsets, retaining repeated/shared coordinates. Arithmetic probe weights; no exact-plane quadrature or independent-probe claim."},
            native_metadata=metadata,
            band_integration={"bands_rad_m": [list(b) for b in band100.BANDS],
                "method": "Linear interpolation at exact physical band endpoints followed by trapezoidal density integration."},
            claim_limits=["Boundary projection only; no CFD evolution, full-distance comparison or experimental-fidelity verdict.",
                "Both PSDs remove each signal/probe's own full-record temporal mean.",
                "Overlapping Welch segments and native point subsets are not independent replicates.",
                "Native all-point maximum errors and sampled temporal-band losses have different coverage.",
                "The fixed mean-convection conversion and finite temporal sample band remain explicit assumptions."],
            inputs={"sha256": hashes}, software={"python": sys.version, "numpy": np.__version__})
        if args.no_plot:
            result["figure"] = {"available": False, "reason": "--no-plot"}
        else:
            os.environ.setdefault("MPLCONFIGDIR", str(output/"mplconfig"))
            result["figure"] = plot_result(result, output)
        for path in protected:
            require(m7.sha(path) == hashes[str(path)], "projection input changed during analysis: " + str(path))
        analysis_path = output/"analysis.json"
        band100.write_json(analysis_path, result)
        receipt.update(status="completed", analysis_sha256=m7.sha(analysis_path), figure=result["figure"],
            finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        band100.write_json(receipt_path, receipt, atomic=True)
        return result
    except BaseException as error:
        receipt.update(status="failed", error_type=type(error).__name__, error=str(error),
            finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        band100.write_json(receipt_path, receipt, atomic=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--stem", default="projection-temporal")
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory outside native case/source")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--probe-batch", type=int, default=256)
    parser.add_argument("--welch-segment-samples", type=int, default=512)
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            analyze(args)
    except (ValueError, OSError, KeyError, TypeError, OverflowError, FloatingPointError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
