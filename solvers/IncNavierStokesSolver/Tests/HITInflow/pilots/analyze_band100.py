#!/usr/bin/env python3
"""Temporal-only diagnostics for the CBC100 feasibility pilot; run on Linux.

The transverse grid is an explicitly sampled set of probes, not an exact plane
quadrature or a collection of independent realizations. No periodic-z transform,
full-distance validation, or experimental-fidelity verdict is produced here.
"""

import argparse
import datetime
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True

import numpy as np

import analyze_m7 as m7


HERE = Path(__file__).resolve().parent
SCHEMA = "NEKTAR_HIT_BAND100_ANALYSIS_V1"
RECEIPT_SCHEMA = "NEKTAR_HIT_BAND100_RECEIPT_V1"
BANDS = ((10.0, 20.0), (20.0, 40.0), (40.0, 70.0), (70.0, 100.0))
COMPONENTS = ("u", "v", "w")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite JSON number")
    elif isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            finite_tree(item)


def read_json(path):
    def invalid(token):
        raise ValueError("nonfinite JSON token: " + token)
    result = json.loads(Path(path).read_text(encoding="utf-8"), parse_constant=invalid)
    finite_tree(result)
    return result


def write_json(path, value, atomic=False):
    finite_tree(value)
    destination = path.with_name(path.name + ".tmp") if atomic else path
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + "\n")
    if atomic:
        os.replace(destination, path)


def integrate_density(k, density, lower, upper):
    """Integrate the piecewise-linear density after inserting exact endpoints."""
    x, y = np.asarray(k, dtype=np.float64), np.asarray(density, dtype=np.float64)
    require(x.ndim == y.ndim == 1 and len(x) == len(y) and len(x) >= 2,
            "density needs matching one-dimensional coordinate/value arrays")
    require(np.isfinite(x).all() and np.isfinite(y).all() and np.all(y >= 0)
            and np.all(np.diff(x) > 0), "density must be finite, nonnegative and ordered")
    require(math.isfinite(lower) and math.isfinite(upper)
            and x[0] <= lower < upper <= x[-1], "band is outside the available density")
    inside = (x > lower) & (x < upper)
    xx = np.concatenate(([lower], x[inside], [upper]))
    yy = np.concatenate(([np.interp(lower, x, y)], y[inside], [np.interp(upper, x, y)]))
    result = float(np.sum((xx[1:] - xx[:-1]) * (yy[1:] + yy[:-1]) * .5))
    require(math.isfinite(result) and result >= 0, "integrated density is not finite/nonnegative")
    return result


def band_powers(spectrum):
    return [{"k1_limits_rad_m": list(band),
             "power_m2_s2": integrate_density(spectrum["k1_rad_m"],
                 spectrum["F_component_m3_s2"], *band)} for band in BANDS]


def component_spectra(values, dt, speed, segment=512, batch=256):
    """All three temporal powers, averaged after each probe's own centering."""
    require(values.ndim == 4 and values.shape[-1] == 3,
            "station velocity must have time,y,z,component dimensions")
    require(segment >= 4 and segment % 2 == 0 and len(values) >= segment,
            "record must cover an even Welch segment of at least four samples")
    result = {}
    for c, name in enumerate(COMPONENTS):
        signal = values[..., c].reshape(len(values), -1)
        estimates = {}
        for key, length, window, overlap in (("boxcar", None, "boxcar", 0),
                ("welch", segment, "hann", segment // 2)):
            estimate = m7.longitudinal_spectrum(m7.one_sided_psd(
                signal, dt, length, window, overlap, batch), speed)
            estimate["F_component_m3_s2"] = estimate.pop("E11_m3_s2")
            estimate["F_component_integral_m2_s2"] = estimate.pop("E11_integral_m2_s2")
            estimate["component_density_name"] = ("F11", "F22", "F33")[c]
            estimate["component_alignment"] = ("longitudinal" if c == 0 else "transverse")
            estimate["band_powers"] = band_powers(estimate)
            estimate["F_component_100_m3_s2"] = float(np.interp(
                100, estimate["k1_rad_m"], estimate["F_component_m3_s2"]))
            estimates[key] = estimate
        result[name] = estimates
    return result


def validate_output(manifest, output, source, protected=()):
    manifest, output, source = (Path(p).resolve() for p in (manifest, output, source))
    require(not output.exists(), "analysis directory must be fresh; never overwrite prior results")
    require(m7.outside(output, manifest.parent) and m7.outside(output, source),
            "analysis directory must be outside the case and source directories")
    for path in protected:
        path = Path(path).resolve()
        require(output != path and m7.outside(path, output),
                "analysis directory cannot contain or overwrite an input")
    return output


def control_identity(manifest, window, spectrum_hashes):
    """Physical/random identities that must agree before matched comparisons."""
    synthesis, mesh, physics = (manifest[key] for key in ("synthesis", "mesh", "physics"))
    return {
        "window": window,
        "stations_x_m": manifest["sampling"]["stations_x_m"],
        "probe_coordinates_m": [p["xyz"] for p in manifest["sampling"]["probes"]],
        "probe_shape_yz": manifest["sampling"]["probe_shape_yz"],
        "seed_uint64": str(synthesis["seed_uint64"]),
        "counts": synthesis["counts"], "periods_m": synthesis["periods_m"],
        "physics": {key: physics.get(key) for key in ("mean_velocity_m_s", "kinvis_m2_s",
            "normalization", "interpolation", "low_extension", "high_extension")},
        "spectrum_sha256": sorted(spectrum_hashes),
        "domain_lengths_m": [manifest["knobs"]["length"], mesh["period_y_m"], mesh["period_z_m"]],
        "initial_conditions": {key: manifest.get("initial_conditions", {}).get(key)
            for key in ("kind", "files_receipt_sha256")},
    }


def clock_equal(left, right):
    if isinstance(left, dict) and isinstance(right, dict):
        return left.keys() == right.keys() and all(clock_equal(left[k], right[k]) for k in left)
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(clock_equal(a, b) for a, b in zip(left, right))
    if isinstance(left, (float, int)) and isinstance(right, (float, int)):
        return math.isclose(left, right, rel_tol=2e-12, abs_tol=2e-14)
    return left == right


def difference(value, reference):
    value, reference = float(value), float(reference)
    require(math.isfinite(value) and math.isfinite(reference) and min(value, reference) >= 0,
            "comparison quantities must be finite/nonnegative")
    relative = 100 * (value - reference) / reference if reference else (0.0 if value == 0 else None)
    if relative is not None and not math.isfinite(relative):
        relative = None
    return {"case": value, "reference": reference,
            "absolute_difference": value - reference,
            "signed_relative_change_percent": relative,
            "zero_reference": reference == 0}


def compare_analyses(current, reference):
    require(current.get("schema") == reference.get("schema") == SCHEMA,
            "reference analysis has an incompatible schema")
    finite_tree(current)
    finite_tree(reference)
    own, other = current["control_identity"], reference["control_identity"]
    reasons = []
    for key in sorted(set(own) | set(other)):
        equal = (clock_equal(own.get(key), other.get(key)) if key == "window"
                 else own.get(key) == other.get(key))
        if not equal:
            reasons.append(key + " differs")
    bandwidth_changed = current["synthesis_band"] != reference["synthesis_band"]
    numerical_changed = current["numerical_settings"] != reference["numerical_settings"]
    result = {"matched_controls_eligible": not reasons, "ineligible_reasons": reasons,
        "bandwidth_changed": bandwidth_changed, "numerical_settings_changed": numerical_changed,
        "kind": ("input_bandwidth_model_sensitivity" if bandwidth_changed else "matched_numerical_control"),
        "eligible_for_numerical_refinement_screen": not reasons and not bandwidth_changed,
        "case_numerical_settings": current["numerical_settings"],
        "reference_numerical_settings": reference["numerical_settings"],
        "interpretation": "Aligned sampled-probe statistics, not independent replicas. Bandwidth changes alter the physical input. Spatial/Fourier refinement with SVV also changes the stabilization operator.",
        "stations": []}
    if reasons:
        return result
    rows, moment_changes, spectral_changes = [], [], []
    count = len(own["stations_x_m"])
    require(count > 0 and all(len(rows) == count for rows in (current["spectra"], reference["spectra"],
            current["statistics"]["planes"], reference["statistics"]["planes"])),
            "reference spectral/statistical station count differs")
    for own_stats, ref_stats, own_psd, ref_psd in zip(current["statistics"]["planes"],
            reference["statistics"]["planes"], current["spectra"], reference["spectra"]):
        require(own_stats["station_x_m"] == ref_stats["station_x_m"]
                == own_psd["station_x_m"] == ref_psd["station_x_m"], "station order differs")
        quantities = {name: difference(a, b) for name, a, b in zip(
            ("u_variance", "v_variance", "w_variance", "K"),
            own_stats["variance_m2_s2"] + [own_stats["K_m2_s2"]],
            ref_stats["variance_m2_s2"] + [ref_stats["K_m2_s2"]])}
        moment_changes.extend(v["signed_relative_change_percent"] for v in quantities.values())
        components, endpoints = {}, {}
        for c in COMPONENTS:
            own_bands = own_psd["components"][c]["welch"]["band_powers"]
            ref_bands = ref_psd["components"][c]["welch"]["band_powers"]
            require(len(own_bands) == len(ref_bands) == len(BANDS), "reference band count differs")
            bands = []
            for a, b, expected in zip(own_bands, ref_bands, BANDS):
                require(a["k1_limits_rad_m"] == b["k1_limits_rad_m"] == list(expected),
                        "reference physical bands differ")
                change = difference(a["power_m2_s2"], b["power_m2_s2"])
                spectral_changes.append(change["signed_relative_change_percent"])
                bands.append({"k1_limits_rad_m": list(expected), **change})
            components[c] = bands
            endpoints[c] = difference(own_psd["components"][c]["welch"]["F_component_100_m3_s2"],
                                      ref_psd["components"][c]["welch"]["F_component_100_m3_s2"])
        rows.append({"station_x_m": own_stats["station_x_m"], "variances_and_K": quantities,
                     "welch_band_power_changes": components, "welch_endpoint_100_changes": endpoints})
    def maximum(changes):
        return None if not changes or any(v is None for v in changes) else max(abs(v) for v in changes)
    mv, mp = maximum(moment_changes), maximum(spectral_changes)
    result.update(stations=rows, maximum_absolute_variance_K_change_percent=mv,
        maximum_absolute_band_power_change_percent=mp,
        informational_within_2_percent_variances_K=(mv <= 2 if mv is not None else None),
        informational_within_5_percent_band_power=(mp <= 5 if mp is not None else None),
        screen_caveat="Engineering screens only; no sampling-precision or experimental-fidelity verdict. A bandwidth comparison is model sensitivity even when thresholds pass.")
    return result


def upstream_overlay(source, stations):
    require(source.get("schema") == "nektar.hit.cbc1971.source.v1", "unexpected CBC source schema")
    rows = [[float(row[0]) * 100, float(row[1]) * 1e-6]
            for row in source["table2a"]["rows"] if row[1] is not None and 0 < float(row[0]) * 100 <= 100]
    return {"available": any(x == 0 for x in stations), "station_x_m": 0.0,
        "theta": 42, "Table2a_F11_SI": rows,
        "interpretation": "Raw upstream one-sided longitudinal measurements, descriptive overlay at x=0 only. No experimental match is assigned to other pilot stations; no fitting or fidelity score."}


def plot_spectra(result, output):
    try:
        import matplotlib
    except ImportError:
        return {"available": False, "reason": "Matplotlib is unavailable"}
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    stations = result["spectra"]
    fig, axes = plt.subplots(1, len(stations), figsize=(4 * len(stations), 4.3), squeeze=False)
    for axis, station in zip(axes[0], stations):
        estimate = station["components"]["u"]["welch"]
        x, y = np.asarray(estimate["k1_rad_m"]), np.asarray(estimate["F_component_m3_s2"])
        keep = (x > 0) & (x <= 110) & (y > 0)
        if np.any(keep):
            axis.loglog(x[keep], y[keep], label="Sampled-probe Welch u")
        if station["station_x_m"] == 0:
            reference = result["upstream_experimental_overlay"]["Table2a_F11_SI"]
            axis.loglog([p[0] for p in reference], [p[1] for p in reference], "kD:",
                        markersize=3, label="CBC theta42, raw full-band")
        axis.axvline(100, color="gray", linestyle="--", linewidth=.8)
        axis.set(xlabel="Longitudinal k1 (rad/m)", title=f"x = {station['station_x_m']:g} m",
                 xlim=(3, 110))
        axis.grid(True, which="both", alpha=.2)
        axis.legend(fontsize=7)
    axes[0][0].set_ylabel("One-sided F11 (m³/s²)")
    fig.suptitle("CBC100 pilot: temporal spectra at sampled probes")
    fig.text(.01, .015, "No independent-probe, exact-plane, full-distance or experimental-fidelity verdict.", fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .93))
    path = output / "band100-spectra.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return {"available": True, "path": path.name, "sha256": m7.sha(path)}


def analyze(args):
    manifest_path = args.manifest.resolve()
    manifest = read_json(manifest_path)
    require(manifest.get("schema") == "NEKTAR_HIT_PILOT_V1", "unexpected pilot manifest schema")
    stations, ny, nz, _, _ = m7.probe_layout(manifest)
    require(manifest["physics"]["normalization"] == "Preserve", "CBC100 pilot requires Preserve normalization")
    require(args.probe_batch > 0, "probe-batch must be positive")
    start = manifest["time"]["record_start_s"] if args.start is None else args.start
    end = manifest["time"]["end_s"] if args.end is None else args.end
    dt, speed = manifest["time"]["sample_dt_s"], manifest["physics"]["mean_velocity_m_s"]
    require(all(math.isfinite(v) for v in (start, end, dt, speed))
            and manifest["time"]["initial_time_s"] <= start < end <= manifest["time"]["end_s"]
            and dt > 0 and speed > 0, "invalid requested recording window, sample interval or speed")
    history = (manifest_path.parent / manifest["history_file"]).resolve()
    case = (manifest_path.parent / manifest["case"]).resolve()
    experimental = args.experimental_source.resolve()
    source = m7.source_root(args.source_root)
    reference_path = args.reference_analysis.resolve() if args.reference_analysis else None
    if reference_path and reference_path.is_dir():
        reference_path /= "analysis.json"
    protected = [manifest_path, history, case, experimental, Path(__file__).resolve(),
                 Path(m7.__file__).resolve(), Path(m7.ah.__file__).resolve()]
    reference = None
    if reference_path:
        reference = read_json(reference_path)
        reference_receipt = reference_path.parent / "receipt.json"
        ref_receipt = read_json(reference_receipt)
        require(ref_receipt.get("schema") == RECEIPT_SCHEMA and ref_receipt.get("status") == "completed"
                and ref_receipt["analysis_sha256"] == m7.sha(reference_path),
                "reference analysis must have its matching completed receipt")
        protected.extend((reference_path, reference_receipt))
    spectrum_hashes = []
    for name, checksum in manifest.get("sha256", {}).items():
        if Path(name).suffix == ".dat" and "spectrum" in Path(name).name.lower():
            path = (manifest_path.parent / name).resolve()
            require(m7.sha(path) == checksum, "spectrum changed since preparation")
            protected.append(path)
            spectrum_hashes.append(checksum)
    require(spectrum_hashes, "CBC pilot needs a bound tabulated spectrum")
    output = validate_output(manifest_path, args.output, source, protected)
    cache, analysis_path = output / "history.npz", output / "analysis.json"
    m7.validate_paths(manifest_path, history, cache, analysis_path, source, protected)
    output.mkdir(parents=True)
    receipt_path = output / "receipt.json"
    receipt = {"schema": RECEIPT_SCHEMA, "status": "running", "command": sys.argv,
        "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "manifest": str(manifest_path), "output": str(output)}
    write_json(receipt_path, receipt, atomic=True)
    try:
        hashes = {str(path): m7.sha(path) for path in protected if path != history}
        native_stat = m7.stable_stat(history)
        m7.create_cache(cache, history, manifest_path, manifest,
                        min(start, manifest["time"]["record_start_s"]), source)
        times, velocity, cache_receipt = m7.load_cache(cache, manifest_path, manifest)
        tolerance = 2e-12 * max(abs(start), abs(end), dt)
        selected_ids = np.flatnonzero((times >= start - tolerance) & (times <= end + tolerance))
        require(len(selected_ids) >= 32, "record needs at least 32 selected samples")
        selected_times = times[selected_ids]
        require(m7.ah.close(float(selected_times[0]), math.ceil((start-tolerance)/dt)*dt, dt)
                and m7.ah.close(float(selected_times[-1]), math.floor((end+tolerance)/dt)*dt, dt),
                "history does not cover the complete requested sampled window")
        selected = velocity[selected_ids[0]:selected_ids[-1] + 1]
        stats, _ = m7.station_statistics(selected_times, selected, stations, args.probe_batch)
        spectra = []
        for index, station in enumerate(stations):
            components = component_spectra(selected[:, index], dt, speed,
                                            args.welch_segment_samples, args.probe_batch)
            for c, component in enumerate(COMPONENTS):
                variance = stats["planes"][index]["variance_m2_s2"][c]
                integral = components[component]["boxcar"]["F_component_integral_m2_s2"]
                residual = abs(integral - variance) / variance if variance else abs(integral)
                require(residual <= 5e-10, "boxcar PSD does not recover sampled-probe variance")
                components[component]["boxcar"]["relative_difference_from_full_window_variance"] = residual
            spectra.append({"station_x_m": station, "components": components})
        window = {"requested_start_s": start, "requested_end_s": end,
            "first_last_samples_s": [float(selected_times[0]), float(selected_times[-1])],
            "samples": len(selected_times), "sample_dt_s": dt}
        hashes[str(history)] = cache_receipt["metadata"]["native_history_sha256"]
        numerical = {"mesh": manifest["mesh"], "time_step_s": manifest["time"]["dt_s"],
                     "svv": manifest["svv"]}
        result = {"schema": SCHEMA, "command": sys.argv,
            "verification_status": "ANALYZED_TEMPORAL_SAMPLED_PROBES",
            "case_directory": str(manifest_path.parent), "window": window,
            "sampling": {"probe_shape_yz": [ny, nz], "probes_per_station": ny*nz,
                "interpretation": "Arithmetic average of own-probe temporal moments/powers at fixed sampled locations. This grid need not resolve all transverse modes; no exact full-plane quadrature, z spectrum or independent-probe claim."},
            "control_identity": control_identity(manifest, window, spectrum_hashes),
            "synthesis_band": {key: manifest["synthesis"][key]
                               for key in ("cutoff_rad_m", "taper_start_rad_m")},
            "numerical_settings": numerical, "statistics": stats, "spectra": spectra,
            "band_integration": {"bands_rad_m": [list(b) for b in BANDS],
                "method": "Piecewise-linear density versus angular wavenumber; interpolate exact endpoints, then trapezoidal integration. Identical physical bands in every control."},
            "upstream_experimental_overlay": upstream_overlay(read_json(experimental), stations),
            "claim_limits": ["No full-distance or experimental-fidelity verdict.",
                "One realization and correlated probes/overlapping Welch segments are not independent replicas.",
                "Raw upstream Table2a is descriptive; shell truncation changes F11 even below its cutoff.",
                "No experimental reference is assigned to interior short-pilot stations.",
                "Temporal-to-spatial conversion uses prescribed mean speed; intrinsic evolution is retained.",
                "Native histories have no explicit antialias prefilter; temporal aliasing must be checked separately.",
                "Conditional block screens do not prove independent blocks, stationarity or confidence coverage."],
            "software": {"python": sys.version, "numpy": np.__version__},
            "inputs": {"sha256": hashes, "cache_sha256": cache_receipt["cache_sha256"],
                "native_history_hash_scope": "Hashed twice during immutable-cache creation; unchanged file stat required through analysis."}}
        if reference is not None:
            result["reference_comparison"] = {"reference_path": str(reference_path),
                "reference_sha256": m7.sha(reference_path), **compare_analyses(result, reference)}
        if not args.no_plot:
            os.environ.setdefault("MPLCONFIGDIR", str(output / "mplconfig"))
            result["figure"] = plot_spectra(result, output)
        else:
            result["figure"] = {"available": False, "reason": "--no-plot"}
        require(m7.stable_stat(history) == native_stat, "native history changed during analysis")
        for path in protected:
            if path != history:
                require(m7.sha(path) == hashes[str(path)], "input changed during analysis: " + str(path))
        write_json(analysis_path, result)
        receipt.update(status="completed", analysis_sha256=m7.sha(analysis_path),
            cache_sha256=cache_receipt["cache_sha256"], figure=result["figure"],
            finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(receipt_path, receipt, atomic=True)
        return result
    except BaseException as error:
        receipt.update(status="failed", error_type=type(error).__name__, error=str(error),
            finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(receipt_path, receipt, atomic=True)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Fresh directory outside case/source")
    parser.add_argument("--reference-analysis", type=Path,
                        help="Completed analysis.json or its directory, with a matching receipt")
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--experimental-source", type=Path, default=HERE.parent / "cbc1971-source.json")
    parser.add_argument("--start", type=float)
    parser.add_argument("--end", type=float)
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
