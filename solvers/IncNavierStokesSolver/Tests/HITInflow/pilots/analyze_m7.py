#!/usr/bin/env python3
"""Reproducible, descriptive CBC1971 comparisons for the limited-band HIT case.

Run on Linux with NumPy. Cache validated native ASCII once, retaining earlier
upstream samples for explicitly labelled delayed correlations. No downstream
fitting, independent-probe assumption, universal experimental error, or full-DNS
accuracy claim is made. Existing analyze_history estimators supply the scalar
moment summaries and conditional block intervals.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np

import analyze_history as ah

HERE = Path(__file__).resolve().parent
CACHE_SCHEMA = "NEKTAR_HIT_M7_HISTORY_CACHE_V1"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_stat(path):
    value = path.stat()
    return value.st_size, value.st_mtime_ns


def source_root(explicit=None):
    candidates = [Path(explicit).resolve()] if explicit is not None else []
    if explicit is None:
        for script in (Path(__file__).resolve(), Path(ah.__file__).resolve()):
            candidates.extend(script.parents)
    for candidate in candidates:
        if ((candidate/"library/SolverUtils").is_dir()
                and (candidate/"solvers/IncNavierStokesSolver").is_dir()):
            return candidate
    raise ValueError("cannot locate source checkout; pass --source-root when using external script copies")


def outside(path, directory):
    try:
        path.relative_to(directory)
    except ValueError:
        return True
    return False


def validate_paths(manifest_path, history, cache, output=None, source=None, protected=()):
    """Validate destinations before allocating arrays or reading huge histories."""
    manifest_path, history, cache = (Path(p).resolve() for p in (manifest_path, history, cache))
    checkout = source_root(source)
    sidecar = cache.with_suffix(cache.suffix+".json")
    targets = [cache, sidecar] + ([Path(output).resolve()] if output is not None else [])
    require(len(set(targets)) == len(targets), "analysis output/cache/sidecar paths must be distinct")
    inputs = {manifest_path, history, Path(__file__).resolve(), Path(ah.__file__).resolve()}
    inputs.update(Path(p).resolve() for p in protected)
    for target in targets:
        require(target not in inputs, "output/cache/sidecar cannot overwrite an input or native history")
        require(outside(target, manifest_path.parent), "output/cache/sidecar must be outside the case directory")
        require(outside(target, checkout), "output/cache/sidecar must be outside the source checkout")
    if output is not None:
        require(not Path(output).resolve().exists(), "analysis output must be a fresh path; never overwrite an existing file")
    return checkout


def probe_layout(manifest):
    sampling = manifest["sampling"]
    stations = sampling["stations_x_m"]
    ny, nz = sampling["probe_shape_yz"]
    require(stations and len(set(stations)) == len(stations) and
            all(math.isfinite(x) for x in stations), "invalid station coordinates")
    require(all(isinstance(n, int) and not isinstance(n, bool) and n >= 2 for n in (ny, nz)),
            "probe-grid counts must be integers >=2")
    probes = sampling["probes"]
    expected = [(s, y, z) for s in range(len(stations)) for y in range(ny) for z in range(nz)]
    require(len(probes) == len(expected) and [p["id"] for p in probes] == list(range(len(probes)))
            and [(p["station_index"], p["iy"], p["iz"]) for p in probes] == expected,
            "cache requires complete canonical station,y,z probe order")
    y0, z0 = probes[0]["xyz"][1:]
    width_y, width_z = manifest["mesh"]["period_y_m"], manifest["mesh"]["period_z_m"]
    for point, (s, y, z) in zip(probes, expected):
        xyz = point["xyz"]
        require(len(xyz) == 3 and all(math.isfinite(v) for v in xyz) and
                all(ah.close(a, b, 1e-3) for a, b in zip(xyz,
                    (stations[s], y0+y*width_y/ny, z0+z*width_z/nz))),
                "probes must share the same uniform half-open yz grid at every station")
    signature = hashlib.sha256(json.dumps(probes, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    return stations, ny, nz, probes, signature


def create_cache(cache, history, manifest_path, manifest, cache_start, source=None):
    cache, history, manifest_path = (Path(p).resolve() for p in (cache, history, manifest_path))
    validate_paths(manifest_path, history, cache, source=source)
    sidecar = cache.with_suffix(cache.suffix+".json")
    require(not cache.exists() and not sidecar.exists(), "cache/sidecar already exists; never overwrite an immutable cache")
    stations, ny, nz, probes, signature = probe_layout(manifest)
    initial = manifest["time"]["initial_time_s"]
    end, dt = manifest["time"]["end_s"], manifest["time"]["sample_dt_s"]
    start = initial if cache_start is None else cache_start
    require(math.isfinite(start) and initial <= start <= manifest["time"]["record_start_s"],
            "cache-start must lie between initial time and recording start")
    tolerance = 2e-12 * max(abs(start), abs(end), dt)
    first_step = math.ceil((start - initial - tolerance) / dt)
    last_step = math.floor((end - initial + tolerance) / dt)
    count = last_step - first_step + 1
    require(count >= 2, "cache window needs >=2 samples")
    times = np.empty(count, dtype=np.float64)
    velocity = np.empty((count, len(stations), ny, nz, 3), dtype=np.float64)
    before = stable_stat(history)
    manifest_hash = sha(manifest_path)
    parsed_hash = hashlib.sha256()
    selected = 0
    for time, flat in ah.history_frames(history, probes, dt, parsed_hash):
        if start - tolerance <= time <= end + tolerance:
            require(selected < count, "native history contains unexpected extra samples")
            times[selected] = time
            velocity[selected] = np.asarray(flat, dtype=np.float64).reshape(len(stations), ny, nz, 3)
            selected += 1
    require(selected == count and ah.close(times[0], initial + first_step * dt, dt)
            and ah.close(times[-1], initial + last_step * dt, dt), "native cache window is incomplete")
    metadata = {"schema": CACHE_SCHEMA, "manifest_sha256": manifest_hash,
                "probe_order_sha256": signature, "native_history": str(history),
                "native_history_sha256": parsed_hash.hexdigest(),
                "shape_T_station_y_z_component": list(velocity.shape),
                "dtype": "float64", "sample_dt_s": dt,
                "first_last_times_s": [float(times[0]), float(times[-1])],
                "velocity_components": ["u", "v", "w"], "pressure": "validated by native reader, not cached"}
    cache.parent.mkdir(parents=True, exist_ok=True)
    created = False
    try:
        with cache.open("xb") as destination:
            created = True
            np.savez(destination, time=times, velocity=velocity,
                     metadata_json=np.asarray(json.dumps(metadata, sort_keys=True)))
        after_hash = sha(history)
        require(before == stable_stat(history) and parsed_hash.hexdigest() == after_hash
                and sha(manifest_path) == manifest_hash,
                "native history changed during cache creation")
        receipt = {"schema": CACHE_SCHEMA, "cache_sha256": sha(cache), "cache_bytes": cache.stat().st_size,
                   "native_history_sha256_before_cache_write": parsed_hash.hexdigest(),
                   "native_history_sha256_after_cache_write": after_hash, "metadata": metadata}
        with sidecar.open("x", encoding="utf-8") as destination:
            json.dump(receipt, destination, indent=2, allow_nan=False)
            destination.write("\n")
    except Exception:
        if created:
            cache.unlink(missing_ok=True)
        raise


def load_cache(cache, manifest_path, manifest):
    _, ny, nz, _, signature = probe_layout(manifest)
    sidecar = cache.with_suffix(cache.suffix + ".json")
    receipt = json.loads(sidecar.read_text())
    require(receipt["schema"] == CACHE_SCHEMA and receipt["cache_sha256"] == sha(cache), "cache hash/schema mismatch")
    with np.load(cache, allow_pickle=False) as archive:
        require(set(archive.files) == {"time", "velocity", "metadata_json"}, "unexpected NPZ keys")
        metadata = json.loads(archive["metadata_json"].item())
        times, velocity = archive["time"], archive["velocity"]
    require(metadata == receipt["metadata"] and metadata["manifest_sha256"] == sha(manifest_path)
            and metadata["probe_order_sha256"] == signature, "cache metadata does not match its manifest/receipt")
    require(receipt["native_history_sha256_before_cache_write"] == receipt["native_history_sha256_after_cache_write"]
            == metadata["native_history_sha256"], "native cache source hashes disagree")
    shape = (len(times), len(manifest["sampling"]["stations_x_m"]), ny, nz, 3)
    require(times.dtype == np.float64 and velocity.dtype == np.float64 and times.ndim == 1
            and velocity.shape == shape and list(shape) == metadata["shape_T_station_y_z_component"], "invalid cache shape/dtype")
    require(len(times) >= 2 and np.isfinite(times).all() and np.isfinite(velocity).all(), "nonfinite/empty cache")
    dt = manifest["time"]["sample_dt_s"]
    require(metadata["sample_dt_s"] == dt and all(ah.close(a, b, dt) for a, b in
            zip((float(times[0]), float(times[-1])), metadata["first_last_times_s"])),
            "cache time metadata does not match its arrays")
    tolerance = 64 * sys.float_info.epsilon * max(float(np.max(np.abs(times))), dt) + 2e-12 * dt
    require(np.all(np.diff(times) > 0) and np.all(np.abs(np.diff(times) - dt) <= tolerance), "cache times are nonuniform")
    times.flags.writeable = velocity.flags.writeable = False
    return times, velocity, receipt


def centered(values):
    """Stable temporal centering, retaining exactly zero constant signals."""
    shifted = values-values[0]
    mean_shift = np.mean(shifted, axis=0, dtype=np.float64)
    return values[0]+mean_shift, shifted-mean_shift


def moments(values, batch=256):
    flat = values.reshape(len(values), -1, 3)
    result = ah.Moments(flat.shape[1] * 3)
    result.n = len(flat)
    mean = np.empty(flat.shape[1:], dtype=np.float64)
    m2 = np.empty_like(mean)
    for first in range(0, flat.shape[1], batch):
        mean[first:first+batch], delta = centered(flat[:, first:first+batch])
        m2[first:first + batch] = np.einsum("tpc,tpc->pc", delta, delta)
    result.mean, result.m2 = mean.reshape(-1), m2.reshape(-1)
    return result


def scalar_values(plane):
    return list(plane["variance_m2_s2"]) + [plane["K_m2_s2"]]


def station_statistics(times, velocity, stations, batch):
    points_per_station = velocity.shape[2] * velocity.shape[3]
    ids = [list(range(s * points_per_station, (s + 1) * points_per_station)) for s in range(len(stations))]
    full = moments(velocity, batch)
    planes = ah.plane_summary(full, ids)
    windows = {}
    for name, first, last in (("full", 0, len(times)), ("first_half", 0, len(times)//2),
                              ("final_half", len(times)//2, len(times)), ("last_75_percent", len(times)//4, len(times))):
        summary = planes if name == "full" else ah.plane_summary(moments(velocity[first:last], batch), ids)
        windows[name] = {"samples": last-first, "first_last_times_s": [float(times[first]), float(times[last-1])],
                         "planes": summary}
    block_sets = {}
    mean_full = full.mean.reshape(-1, 3)
    for count in (4, 8, 16):
        per = len(times) // count
        require(per >= 2, "record needs at least two samples in each of16 blocks")
        blocks = [moments(velocity[b*per:(b+1)*per], batch) for b in range(count)]
        results = []
        for s, point_ids in enumerate(ids):
            energy, own = [], []
            for block in blocks:
                central = block.m2.reshape(-1, 3)[point_ids] / block.n
                offset = block.mean.reshape(-1, 3)[point_ids] - mean_full[point_ids]
                energy.append(np.mean(central + offset*offset, axis=0).tolist())
                own.append(np.mean(central, axis=0).tolist())
            intervals = [ah.block_interval([v[c] for v in energy]) for c in range(3)]
            intervals.append(ah.block_interval([0.5*sum(v) for v in energy]))
            results.append({"station_x_m": stations[s], "component_variance_and_K": intervals,
                            "own_block_centered_component_variances": own})
        block_sets[str(count)] = {"count": count, "samples_per_block": per,
            "sample_span_s": float(times[per-1]-times[0]), "excluded_tail_samples": len(times)-per*count,
            "definition": "Point energies centered on each probe's full-window mean; equal-count nonoverlapping blocks. Intervals concern mean block energy, not the distinct full-window variance.",
            "interval_assumptions": "Conditional Student-t approximation: stationarity/weak block dependence are not established by nonoverlap; spatial probes are not independent samples.", "planes": results}
    qualification = []
    for s, x in enumerate(stations):
        values = scalar_values(planes[s])
        half = scalar_values(windows["final_half"]["planes"][s])
        changes = [100*abs(ah.relative_difference(a, b)) for a, b in zip(values, half)]
        widths = {}
        for count, block in block_sets.items():
            widths[count] = []
            for value, interval in zip(values, block["planes"][s]["component_variance_and_K"]):
                low, high = interval["conditional_95_percent_t_interval"]
                width = (high-low)/2
                widths[count].append(100*width/value if value else (0.0 if width == 0 else None))
        failed = any(change > 2 for change in changes) or any(
            value is None or value > 2 for row in widths.values() for value in row)
        qualification.append({"station_x_m": x, "order": ["u_variance", "v_variance", "w_variance", "K"],
            "full_to_final_half_change_percent": changes, "conditional_interval_halfwidth_percent_of_full": widths,
            "engineering_screen": "not_qualified_at_2_percent_sampling_precision" if failed else "screens_passed_only",
            "caveat": "Engineering flags, not hypothesis tests, independence proof or asymptotic confidence coverage."})
        covariance = np.zeros((3, 3))
        data = velocity[:, s].reshape(len(times), -1, 3)
        for first in range(0, data.shape[1], batch):
            _, delta = centered(data[:, first:first+batch])
            covariance += np.einsum("tpi,tpj->ij", delta, delta)/(len(times)*data.shape[1])
        planes[s]["reynolds_stress_covariance_m2_s2"] = covariance.tolist()
    return {"planes": [{"station_x_m": x, **p} for x, p in zip(stations, planes)],
            "windows": windows, "blocks": block_sets, "qualification": qualification}, full


def one_sided_psd(data, dt, segment=None, window="boxcar", overlap=0, batch=256):
    """Mean of probe powers, after each probe's full-record mean removal."""
    require(data.ndim == 2 and len(data) >= 2 and dt > 0, "PSD expects time,probe data and positive dt")
    length = len(data) if segment is None else segment
    require(4 <= length <= len(data) and 0 <= overlap < length, "invalid PSD segment/overlap")
    weight = np.ones(length) if window == "boxcar" else 0.5 - 0.5*np.cos(2*np.pi*np.arange(length)/length)
    norm = float(np.sum(weight*weight))
    starts = list(range(0, len(data)-length+1, length-overlap))
    power = np.zeros(length//2+1)
    target, full_variance = 0.0, 0.0
    for first in range(0, data.shape[1], batch):
        _, fluctuation = centered(data[:, first:first+batch])
        full_variance += float(np.sum(fluctuation*fluctuation)) / len(data)
        for start in starts:
            values = fluctuation[start:start+length] * weight[:, None]
            transform = np.fft.rfft(values, axis=0)
            power += np.sum(np.abs(transform)**2, axis=1) * dt/norm
            target += float(np.sum(values*values)) / norm
    divisor = len(starts) * data.shape[1]
    power /= divisor
    target /= divisor
    power[1:-1 if length % 2 == 0 else None] *= 2
    frequency = np.fft.rfftfreq(length, dt)
    integral = float(np.sum(power) / (length*dt))
    residual = abs(integral-target)/target if target else abs(integral)
    require(residual <= 5e-11, "one-sided PSD Parseval check failed")
    return {"frequency_Hz": frequency.tolist(), "Sf_m2_s2_per_Hz": power.tolist(),
            "window": "periodic_Hann" if window == "hann" else "boxcar", "segment_samples": length,
            "overlap_samples": overlap, "segments": len(starts), "excluded_tail_samples": len(data)-(starts[-1]+length),
            "detrending": "one temporal mean per probe over full recording window; no segment detrending or plane centering",
            "frequency_resolution_Hz": 1/(length*dt), "FFT_duration_s": length*dt,
            "window_energy_sum_squared_weights": norm,
            "equivalent_noise_bandwidth_Hz": norm/(dt*float(np.sum(weight))**2),
            "integrated_PSD_m2_s2": integral, "weighted_segment_variance_m2_s2": target,
            "full_window_variance_m2_s2": full_variance/data.shape[1], "parseval_relative_residual": residual,
            "uncertainty": "Spatial and overlapping-segment averages are not independent replicates; no PSD confidence interval inferred."}


def longitudinal_spectrum(spectrum, speed):
    frequency = np.asarray(spectrum["frequency_Hz"])
    k = 2*math.pi*frequency/speed
    energy = speed/(2*math.pi)*np.asarray(spectrum["Sf_m2_s2_per_Hz"])
    integral = float(np.sum(energy)*(k[1]-k[0]))
    target = spectrum["integrated_PSD_m2_s2"]
    residual = abs(integral-target)/target if target else abs(integral)
    require(residual <= 5e-11, "E11 frequency/wavenumber Jacobian check failed")
    spectrum.update(k1_rad_m=k.tolist(), E11_m3_s2=energy.tolist(), E11_integral_m2_s2=integral,
        Jacobian_parseval_relative_residual=residual,
        Taylor_hypothesis="Uses prescribed mean convection; intrinsic evolution/pressure/SVV are not removed")
    return spectrum


def plane_acf(values, dt, max_lag, batch=256):
    data = values.reshape(len(values), -1, 3)
    n = len(data)
    maximum = min(n-1, max_lag)
    nfft = 1 << (2*n-2).bit_length()
    sums = np.zeros((maximum+1, 3))
    for first in range(0, data.shape[1], batch):
        _, fluctuation = centered(data[:, first:first+batch])
        transform = np.fft.rfft(fluctuation, n=nfft, axis=0)
        products = np.fft.irfft(transform*transform.conjugate(), n=nfft, axis=0)[:maximum+1]
        sums += np.sum(products, axis=1)
    covariance = sums / (data.shape[1] * (n-np.arange(maximum+1))[:, None])
    correlation = np.divide(covariance, covariance[0], out=np.zeros_like(covariance), where=covariance[0] != 0)
    return {"lag_s": (dt*np.arange(maximum+1)).tolist(), "component_correlations": correlation.T.tolist(),
        "definition": "Noncircular lag covariance with N-lag pairs, about each probe's full-window mean; normalized by full-window component variance. This is not overlap-Pearson normalization.",
        "zero_variance_components": (covariance[0] == 0).tolist()}


def delayed_pair(times, velocity, source, target, delay, start, end, full_mean, full_variance, primary):
    tolerance = 2e-12 * max(abs(end), abs(start), delay, 1.0)
    queries = times-delay
    mask = (times >= start-tolerance) & (times <= end+tolerance)
    mask &= (queries >= times[0]-tolerance) & (queries <= times[-1]+tolerance)
    if primary:
        mask &= queries >= start-tolerance
    selected = np.flatnonzero(mask)
    if len(selected) < 2:
        return {"available": False, "reason": "fewer than two valid paired samples", "pairs_per_probe": len(selected)}
    query = np.clip(queries[selected], times[0], times[-1])
    hi = np.clip(np.searchsorted(times, query, side="left"), 1, len(times)-1)
    lo = hi-1
    alpha = (query-times[lo])/(times[hi]-times[lo])
    shape = (-1, 1, 1, 1)
    upstream = (1-alpha.reshape(shape))*velocity[lo, source] + alpha.reshape(shape)*velocity[hi, source]
    downstream = velocity[selected, target]
    means_a, da = centered(upstream)
    means_b, db = centered(downstream)
    variance_a = np.mean(da*da, axis=(0, 1, 2))
    variance_b = np.mean(db*db, axis=(0, 1, 2))
    covariance = np.mean(da*db, axis=(0, 1, 2))
    offset = np.mean((means_a-full_mean[source])*(means_b-full_mean[target]), axis=(0, 1))
    denominator_pair = np.sqrt(variance_a*variance_b)
    denominator_full = np.sqrt(full_variance[source]*full_variance[target])
    pearson = [float(c/d) if d else None for c, d in zip(covariance, denominator_pair)]
    normalized = [float(c/d) if d else None for c, d in zip(covariance+offset, denominator_full)]
    require(all(value is None or abs(value) <= 1+2e-11 for value in pearson), "paired Cauchy check failed")
    return {"available": True, "pairs_per_probe": len(selected), "delay_s": delay,
        "downstream_first_last_times_s": [float(times[selected[0]]), float(times[selected[-1]])],
        "upstream_first_last_query_times_s": [float(query[0]), float(query[-1])],
        "actual_paired_sample_span_s": float(times[selected[-1]]-times[selected[0]]),
        "ideal_joint_poststart_overlap_s": max(0.0, end-start-delay),
        "full_downstream_record_covered": len(selected) == int(np.sum((times >= start-tolerance) & (times <= end+tolerance))),
        "upstream_queries_before_record_start": int(np.sum(query < start-tolerance)),
        "overlap_Pearson_components": pearson, "full_window_normalized_covariance_components": normalized,
        "paired_component_variances_upstream": variance_a.tolist(), "paired_component_variances_downstream": variance_b.tolist(),
        "normalizations": "Pearson: pointwise paired-overlap means/variances, then plane covariance normalization. Other estimate: full-record pointwise means and full-record plane variances; it need not lie in[-1,1].",
        "interpolation": "Fixed exact mean-convection delay; linear interpolation in time, no fitted lag. No antialias/interpolation-error correction.",
        "qualification": "Both queries >=record_start" if primary else "Stationary-forced-inlet extension only; earlier upstream samples are not startup-qualified interior samples."}


def first_zero_positive_lobe(separation, correlation):
    radius, value = np.asarray(separation), np.asarray(correlation)
    require(len(radius) == len(value) and len(radius) >= 2 and np.all(np.diff(radius) > 0)
            and np.isfinite(radius).all() and np.isfinite(value).all(), "invalid sampled autocorrelation")
    if value[0] <= 0:
        return {"status": "zero_variance_undefined", "positive_lobe_integral_m": None,
                "first_zero_m": None, "integration_limit_m": None, "truncated_at_half_period": None}
    crossing = np.flatnonzero(value[1:] <= 0)
    if not len(crossing):
        return {"status": "no_zero_within_half_period_truncated",
                "positive_lobe_integral_m": float(np.trapz(value, radius)),
                "first_zero_m": None, "integration_limit_m": float(radius[-1]),
                "truncated_at_half_period": True,
                "caveat": "Positive lobe remains unresolved beyond the finite half-period; not an unbounded integral scale."}
    right = int(crossing[0])+1
    zero = radius[right-1] + value[right-1]*(radius[right]-radius[right-1])/(value[right-1]-value[right])
    integral = float(np.trapz(value[:right], radius[:right])) + .5*value[right-1]*(zero-radius[right-1])
    return {"status": "first_zero_linearly_interpolated", "positive_lobe_integral_m": float(integral),
            "first_zero_m": float(zero), "integration_limit_m": float(zero),
            "truncated_at_half_period": False,
            "caveat": "Trapezoidal sampled positive lobe with a linearly located first zero; a finite periodic diagnostic, not unbounded Lg/Lf."}


def z_correlations(values, width, batch=64):
    n, ny, nz, _ = values.shape
    require(nz % 2 == 0, "half-period z integral requires even probe Nz")
    power = np.zeros((nz//2+1, 3))
    for first in range(0, ny, batch):
        _, fluctuation = centered(values[:, first:first+batch])
        transform = np.fft.rfft(fluctuation, axis=2)
        power += np.sum(np.abs(transform)**2, axis=(0, 1)) / (n*ny*nz)
    result = {"separation_m": (np.arange(nz//2+1)*width/nz).tolist(), "integration_limit_m": width/2,
              "definition": "Periodic z autocovariance of temporal fluctuations, averaged over time/y/z origins. u is transverse; w longitudinal by displacement/component alignment. Finite half-period integrals, not unbounded isotropic Lg/Lf."}
    removed = power.copy()
    removed[0] = 0
    for label, covariance in (("zero_mode_retained", np.fft.irfft(power, n=nz, axis=0)),
                               ("zero_mode_removed", np.fft.irfft(removed, n=nz, axis=0))):
        correlation = np.divide(covariance[:nz//2+1], covariance[0], out=np.zeros((nz//2+1, 3)), where=covariance[0] != 0)
        result[label] = {"component_correlations": correlation.T.tolist(),
            "finite_half_period_integrals_m": np.trapz(correlation, dx=width/nz, axis=0).tolist(),
            "normalization_variances_m2_s2": covariance[0].tolist(),
            "zero_variance_components": (covariance[0] == 0).tolist(),
            "first_zero_positive_lobe": [first_zero_positive_lobe(result["separation_m"], correlation[:, c]) for c in range(3)],
            "full_half_period_interpretation": (
                "For the complete periodic cosine expansion, the full half-period integral contains only the zero-mode contribution; it is not an unbounded integral scale."
                if label == "zero_mode_retained" else
                "For the complete periodic cosine expansion with zero mode removed, the full half-period integral is mathematically zero up to numerical quadrature/roundoff; it is not a measured Lg/Lf.")}
    return result


def table4_targets(source):
    keys = source["table4"]["columns"]
    scales = {"u_rms_cm_s": ("u_rms_m_s", 0.01), "epsilon_cm2_s3": ("epsilon_m2_s3", 1e-4),
              "eta_cm": ("eta_m", .01), "lambda_cm": ("lambda_m", .01), "Lg_cm": ("Lg_m", .01),
              "Lf_cm": ("Lf_m", .01), "Re_lambda": ("Re_lambda", 1.0)}
    result = {}
    for row in source["table4"]["rows"]:
        raw = dict(zip(keys, row))
        result[raw["theta"]] = {new: float(raw[old])*scale if raw[old] is not None else None
                               for old, (new, scale) in scales.items()}
    return result


def linear_reference(x, xp, yp):
    if x < xp[0] or x > xp[-1]:
        return {"value": None, "method": "outside published range; no extrapolation"}
    for i, coordinate in enumerate(xp):
        if abs(x-coordinate) <= 1e-12*max(1.0, abs(x)):
            return {"value": yp[i], "method": "exact published row", "published_coordinate": coordinate}
    hi = int(np.searchsorted(xp, x))
    fraction = (x-xp[hi-1])/(xp[hi]-xp[hi-1])
    return {"value": yp[hi-1]+fraction*(yp[hi]-yp[hi-1]),
            "method": "explicit linear interpolation between published rows; not an exact measurement",
            "bracket_coordinates": [xp[hi-1], xp[hi]], "bracket_values": [yp[hi-1], yp[hi]]}


def comparisons(source, table5, manifest, stations, statistics, spectra, reference_x):
    experiment = source["experiment"]
    mesh, approach, theta0 = (float(experiment[k]) for k in ("grid_mesh_m", "approach_speed_m_s", "reference_theta"))
    speed = manifest["physics"]["mean_velocity_m_s"]
    targets = table4_targets(source)
    result = []
    for column, theta in enumerate(source["stations"], 1):
        x = reference_x + speed*mesh/approach*(float(theta)-theta0)
        nearest = min(range(len(stations)), key=lambda s: abs(stations[s]-x))
        match = nearest if abs(stations[nearest]-x) <= 1e-7*max(1.0, abs(x)) else None
        published = [{"k1_rad_m": float(row[0])*100, "E11_m3_s2": float(row[column])*1e-6,
                      "printed_row": row} for row in source["table2a"]["rows"] if row[column] is not None]
        item = {"theta": theta, "expected_station_x_m": x, "matched_station_x_m": stations[match] if match is not None else None,
                "table4_published_SI": targets[theta], "table2a_published_SI": published}
        if match is not None:
            model = statistics["planes"][match]
            urms = model["rms_m_s"][0]
            item["table4_direct_probe_comparison"] = {"u_rms_m_s": urms,
                "u_rms_relative_difference": urms/targets[theta]["u_rms_m_s"]-1,
                "K_model_m2_s2": model["K_m2_s2"],
                "experimental_K_isotropic_inference_m2_s2": 1.5*targets[theta]["u_rms_m_s"]**2,
                "unavailable_from_probes_alone": ["native molecular dissipation", "eta", "lambda", "Lg", "Lf", "Re_lambda"],
                "caveat": "Only urms is directly compared. K reference is an isotropic inference. Native strain/SVV work and assumption-based surrogates are separate observables."}
            for name, spectrum in spectra[match].items():
                if name not in ("welch_primary", "boxcar_parseval") or not spectrum.get("available", True):
                    continue
                k, e = spectrum["k1_rad_m"], spectrum["E11_m3_s2"]
                mapped = []
                for row in published:
                    resolved = k[1] <= row["k1_rad_m"] <= k[-1]
                    estimate = float(np.interp(row["k1_rad_m"], k, e)) if resolved else None
                    mapped.append({**row, "model_E11_m3_s2": estimate,
                        "model_relative_difference": estimate/row["E11_m3_s2"]-1 if estimate is not None else None,
                        "model_interpolation": "linear interpolation of estimated PSD density; no extrapolation" if resolved else "outside positive FFT-bin range",
                        "k1_below_inlet_radial_cutoff": row["k1_rad_m"] < manifest["synthesis"]["cutoff_rad_m"],
                        "band_caution": "Raw Table2a is not shell-filtered; high radial k contributes even at low k1. No direct band-matched reference is implied."})
                item["table2a_"+name+"_descriptive_comparison"] = mapped
        result.append(item)
    grid = next(g for g in table5["grids"] if g["id"] == "grid_5_08_cm")
    return result, grid, mesh


def analyze(args):
    manifest_path = args.manifest.resolve()
    manifest_hash = sha(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    require(manifest.get("schema") == "NEKTAR_HIT_PILOT_V1", "unsupported manifest schema")
    case_path = manifest_path.parent/manifest["case"]
    history = (args.history or manifest_path.parent/manifest["history_file"]).resolve()
    cache, output = args.cache.resolve(), args.output.resolve()
    checkout = validate_paths(manifest_path, history, cache, output, args.source_root,
        protected=(case_path, args.experimental_source, args.correlation_source))
    spectrum_type = ET.parse(case_path).find("./HITINFLOW/REGION/SPECTRUM").get("TYPE") if case_path.is_file() else None
    family = args.case_family or {"Tabulated": "cbc1971", "VonKarman": "von-karman"}.get(spectrum_type)
    require(family in ("cbc1971", "von-karman"), "cannot identify case family; provide --case-family")
    require(not spectrum_type or family == {"Tabulated": "cbc1971", "VonKarman": "von-karman"}.get(spectrum_type),
            "case-family disagrees with native session spectrum type")
    stations, ny, nz, _, _ = probe_layout(manifest)
    require(nz >= manifest["mesh"]["hom_modes_z"],
            "periodic z diagnostics require probe Nz>=solver HomModesZ; do not silently alias retained modes")
    require(args.probe_batch > 0 and args.welch_segment_samples >= 4 and math.isfinite(args.acf_max_lag_s)
            and args.acf_max_lag_s > 0, "invalid batch/PSD segment/ACF lag size")
    if not cache.exists():
        create_cache(cache, history, manifest_path, manifest, args.cache_start, checkout)
    cache_before = stable_stat(cache)
    times, velocity, cache_receipt = load_cache(cache, manifest_path, manifest)
    require(output != Path(cache_receipt["metadata"]["native_history"]).resolve(),
            "analysis output cannot overwrite cached native history")
    dt, speed, nu = (manifest["time"]["sample_dt_s"], manifest["physics"]["mean_velocity_m_s"], manifest["physics"]["kinvis_m2_s"])
    start = manifest["time"]["record_start_s"] if args.start is None else args.start
    end = manifest["time"]["end_s"] if args.end is None else args.end
    require(all(math.isfinite(v) for v in (start, end, dt, speed, nu)) and 0 <= start < end
            and dt > 0 and speed > 0 and nu > 0, "invalid physical time/speed/viscosity")
    tolerance = 2e-12*max(abs(end), abs(start), dt)
    mask = (times >= start-tolerance) & (times <= end+tolerance)
    indices = np.flatnonzero(mask)
    require(len(indices) >= 32 and ah.close(times[indices[0]], math.ceil((start-tolerance)/dt)*dt, dt)
            and ah.close(times[indices[-1]], math.floor((end+tolerance)/dt)*dt, dt), "cache does not cover the complete requested recording window")
    selected_times = times[indices]
    selected = velocity[indices[0]:indices[-1]+1]
    stats, full = station_statistics(selected_times, selected, stations, args.probe_batch)
    full_mean = full.mean.reshape(len(stations), ny, nz, 3)
    full_variance = np.asarray([p["variance_m2_s2"] for p in stats["planes"]])
    spectra, correlations, transverse, surrogates = [], [], [], []
    max_lag = max(1, round(args.acf_max_lag_s/dt))
    for s, x in enumerate(stations):
        signal = selected[:, s, :, :, 0].reshape(len(selected), -1)
        boxcar = one_sided_psd(signal, dt, batch=args.probe_batch)
        length = args.welch_segment_samples
        welch = one_sided_psd(signal, dt, length, "hann", length//2, args.probe_batch) if length <= len(signal) else {"available": False, "reason": "record shorter than declared Welch segment"}
        for spectrum in (boxcar, welch):
            if spectrum.get("available", True):
                longitudinal_spectrum(spectrum, speed)
                target = stats["planes"][s]["variance_m2_s2"][0]
                difference = spectrum["integrated_PSD_m2_s2"]-target
                spectrum["relative_difference_from_full_window_variance"] = difference/target if target else (0.0 if difference == 0 else None)
        relative = boxcar["relative_difference_from_full_window_variance"]
        require(relative is not None and abs(relative) <= 5e-10,
                "boxcar PSD does not recover temporal-central station variance")
        spectra.append({"station_x_m": x, "welch_primary": welch, "boxcar_parseval": boxcar})
        correlations.append({"station_x_m": x, **plane_acf(selected[:, s], dt, max_lag, args.probe_batch)})
        transverse.append({"station_x_m": x, **z_correlations(selected[:, s], manifest["mesh"]["period_z_m"])})
        k = np.asarray(boxcar["k1_rad_m"])
        e = np.asarray(boxcar["E11_m3_s2"])
        derivative = float(np.sum(k*k*e) * (k[1]-k[0]))
        variance = stats["planes"][s]["variance_m2_s2"][0]
        epsilon = 15*nu*derivative
        lam = math.sqrt(variance/derivative) if derivative else None
        surrogates.append({"station_x_m": x, "resolved_temporal_Taylor_longitudinal_derivative_moment_s2": derivative,
            "epsilon_Taylor_isotropy_surrogate_m2_s3": epsilon, "lambda_surrogate_m": lam,
            "eta_surrogate_m": (nu**3/epsilon)**.25 if epsilon else None,
            "Re_lambda_surrogate": math.sqrt(variance)*lam/nu if lam is not None else None,
            "caveat": "Finite sampled frequency band plus Taylor frozen-convection and local isotropy assumptions. Boxcar endpoint discontinuities can inflate high-frequency derivative moments. Not native molecular/total dissipation, unbounded length or measured microscale."})
    source = json.loads(args.experimental_source.read_text())
    table5 = json.loads(args.correlation_source.read_text())
    require(source["schema"] == "nektar.hit.cbc1971.source.v1"
            and table5["schema"] == "nektar.hit.cbc1971.full_band_correlation.v1", "unexpected experimental source schema")
    comparisons_out, grid, mesh = comparisons(source, table5, manifest, stations, stats, spectra, args.reference_x)
    reference = min(range(len(stations)), key=lambda s: abs(stations[s]-args.reference_x))
    require(abs(stations[reference]-args.reference_x) <= 1e-8, "reference-x must be a sampled station")
    pairs = []
    for s, x in enumerate(stations):
        if x < args.reference_x:
            continue
        delay = (x-args.reference_x)/speed
        primary = delayed_pair(times, velocity, reference, s, delay, start, end, full_mean, full_variance, True)
        secondary = (delayed_pair(times, velocity, reference, s, delay, start, end, full_mean, full_variance, False)
                     if abs(args.reference_x-args.forced_inlet_x) <= 1e-8 else {"available": False, "reason": "reference is not the declared forced inlet"})
        xp, yp = zip(*[(float(row[0]), float(row[1])) for row in grid["rows"]])
        xi = (x-args.reference_x)/mesh
        reference_value = linear_reference(xi, xp, yp)
        brackets = []
        for lag in sorted(set((math.floor(delay/dt), math.ceil(delay/dt)))):
            if lag >= 0:
                brackets.append(delayed_pair(times, velocity, reference, s, lag*dt, start, end, full_mean, full_variance, True))
        pairs.append({"upstream_x_m": args.reference_x, "downstream_x_m": x, "xi_delta_x_over_M": xi,
            "delay_delta_x_over_U_s": delay, "primary_joint_poststart": primary,
            "secondary_stationary_forced_inlet_extension": secondary, "fixed_floor_ceil_delay_sensitivity": brackets,
            "Table5_reference": reference_value, "descriptive_only": "Full-band published curve versus limited-band evolving solution; no fidelity verdict or lag fitting."})
    input_paths = [manifest_path, cache, cache.with_suffix(cache.suffix+".json"), args.experimental_source.resolve(),
                   args.correlation_source.resolve(), Path(__file__).resolve(), Path(ah.__file__).resolve()]
    if case_path.is_file():
        input_paths.append(case_path.resolve())
    require(args.output.resolve() not in input_paths, "output cannot overwrite an input")
    require(stable_stat(cache) == cache_before and sha(cache) == cache_receipt["cache_sha256"], "cache changed during analysis")
    require(sha(manifest_path) == manifest_hash, "manifest changed during analysis")
    result = {"schema": "NEKTAR_HIT_M7_ANALYSIS_V1", "purpose": __doc__, "command": sys.argv,
        "source_root": str(checkout),
        "case_family": family, "case_id": args.case_id or manifest_path.parent.name,
        "case_directory": str(manifest_path.parent), "manifest_path": str(manifest_path),
        "physics": manifest["physics"], "synthesis": manifest["synthesis"],
        "manifest_seed": manifest["synthesis"]["seed_uint64"], "cache_receipt": cache_receipt,
        "window": {"requested_start_s": start, "requested_end_s": end, "first_last_samples_s": [float(selected_times[0]), float(selected_times[-1])],
                   "samples": len(selected_times), "sample_dt_s": dt, "sample_Nyquist_Hz": .5/dt,
                   "prescribed_convection_speed_m_s": speed, "anti_aliasing": "No postprocessing antialias filter imposed; sampling/aliasing qualification remains necessary."},
        "statistics": stats, "longitudinal_spectra": spectra, "temporal_autocorrelations": correlations,
        "periodic_z_correlations": transverse, "assumption_based_temporal_surrogates": surrogates,
        "Table2a_Table4_comparisons": comparisons_out, "mean_convection_delayed_pairs": pairs,
        "published_Table5_selected_grid": grid, "experimental_source": source["source"],
        "claim_limits": ["Own-probe temporal means; instantaneous plane means are never removed.",
            "No independent spatial-probe/segment/seed assumption or universal experimental error bar.",
            "Across CBC/VK families equal seeds intentionally provide common random mode inputs, not four independent realizations.",
            "Block-energy intervals differ from full-window variance; screening flags do not establish confidence coverage.",
            "Far-field timestep/Fourier/span/outlet/SVV-strength errors remain unquantified as recorded by M6.",
            "Raw Table2a at k1 below cutoff is not a radial-shell-filtered experimental reference.",
            "Full-band overlays and Table4 differences are descriptive; native work/adjustment qualification is separate."],
        "software": {"python": sys.version, "numpy": np.__version__}, "sha256": {str(path): sha(path) for path in input_paths}}
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as destination:
        destination.write(json.dumps(result, indent=2, allow_nan=False)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True, help="Immutable NPZ; created from validated native ASCII if absent")
    parser.add_argument("--history", type=Path, help="Cache creation only; defaults to manifest history_file")
    parser.add_argument("--cache-start", type=float, help="Optional retained cache start; default initial time, preserving earlier inlet samples")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-root", type=Path,
                        help="Source checkout protected from output/cache writes; required if neither script is inside a checkout")
    parser.add_argument("--case-family", choices=("cbc1971", "von-karman"), help="Normally inferred from native session spectrum TYPE")
    parser.add_argument("--case-id", help="Defaults to the manifest directory name")
    parser.add_argument("--start", type=float)
    parser.add_argument("--end", type=float)
    parser.add_argument("--probe-batch", type=int, default=256)
    parser.add_argument("--welch-segment-samples", type=int, default=512)
    parser.add_argument("--acf-max-lag-s", type=float, default=.1)
    parser.add_argument("--reference-x", type=float, default=0.0)
    parser.add_argument("--forced-inlet-x", type=float, default=0.0)
    parser.add_argument("--experimental-source", type=Path, default=HERE.parent/"cbc1971-source.json")
    parser.add_argument("--correlation-source", type=Path, default=HERE.parent/"cbc1971_full_band_correlation.json")
    args = parser.parse_args()
    try:
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            analyze(args)
    except (ValueError, OSError, KeyError, TypeError, OverflowError, FloatingPointError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
