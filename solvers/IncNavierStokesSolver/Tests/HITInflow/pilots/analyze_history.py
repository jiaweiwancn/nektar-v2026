#!/usr/bin/env python3
"""Analyze native physical-space HistoryPoints without plane-mean subtraction.

Temporal central moments are computed independently at each fixed probe, then
averaged over each half-open yz plane. Block intervals are conditional sampling
diagnostics, not an experimental uncertainty estimate or proof of independence.
Only the Python standard library is required; the history is streamed twice.
"""

import argparse
from collections import deque
import cmath
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys


class Moments:
    def __init__(self, size):
        self.n = 0
        self.mean = [0.0] * size
        self.m2 = [0.0] * size

    def add(self, values):
        self.n += 1
        for i, value in enumerate(values):
            delta = value - self.mean[i]
            self.mean[i] += delta / self.n
            self.m2[i] += delta * (value - self.mean[i])


def close(a, b, scale=1.0):
    return abs(a - b) <= 2e-12 * max(abs(a), abs(b), scale)


def history_frames(path, probes, sample_dt, hasher=None):
    """Validate native point headers and consecutive equal-time point frames."""
    seen_points = 0
    seen_variables = False
    current = []
    frame_time = previous_time = None
    with path.open("rb") as source:
        for line_number, raw in enumerate(source, 1):
            if hasher is not None:
                hasher.update(raw)
            try:
                line = raw.decode("ascii").strip()
            except UnicodeDecodeError as error:
                raise ValueError(f"line {line_number}: history must be ASCII") from error
            if not line:
                continue
            if line.startswith("#"):
                if current or previous_time is not None:
                    raise ValueError(f"line {line_number}: unexpected header inside history data")
                if "History data for variables" in line:
                    match = re.search(r"\(:(.*?)\)", line)
                    variables = [v for v in match.group(1).split(",") if v] if match else []
                    if seen_variables or variables != ["u", "v", "w", "p"]:
                        raise ValueError("expected one native u,v,w,p history header")
                    seen_variables = True
                elif "Wavespace" in line or "Wavespace" in line.replace(" ", ""):
                    raise ValueError("wave-space histories are not supported")
                else:
                    fields = line.split()
                    if len(fields) == 5 and fields[1].isdigit():
                        point_id = int(fields[1])
                        if point_id != seen_points or point_id >= len(probes):
                            raise ValueError("native history point IDs must match manifest order")
                        xyz = [float(value) for value in fields[2:]]
                        if any(not math.isfinite(v) for v in xyz) or any(
                                not close(a, b, 1e-3) for a, b in zip(xyz, probes[point_id]["xyz"])):
                            raise ValueError(f"native history coordinate mismatch at point {point_id}")
                        seen_points += 1
                continue
            if not seen_variables or seen_points != len(probes):
                raise ValueError("history header is incomplete or probe count differs from manifest")
            fields = line.split()
            if len(fields) != 5:
                raise ValueError(f"line {line_number}: expected time,u,v,w,p")
            try:
                values = [float(value) for value in fields]
            except ValueError as error:
                raise ValueError(f"line {line_number}: invalid numeric history row") from error
            if any(not math.isfinite(value) for value in values):
                raise ValueError(f"line {line_number}: non-finite history value")
            time = values[0]
            if not current:
                frame_time = time
                if previous_time is not None:
                    # A subtraction of absolute solver times incurs roundoff
                    # on their magnitude, not merely on the small interval.
                    interval_tolerance = (64 * sys.float_info.epsilon *
                                          max(abs(time), abs(previous_time), sample_dt)
                                          + 2e-12 * sample_dt)
                    if time <= previous_time or abs(time - previous_time - sample_dt) > interval_tolerance:
                        raise ValueError(f"line {line_number}: nonuniform/duplicate history time")
            elif not close(time, frame_time, sample_dt):
                raise ValueError(f"line {line_number}: incomplete equal-time point frame")
            current.extend(values[1:4])
            if len(current) == 3 * len(probes):
                yield frame_time, current
                previous_time = frame_time
                current = []
    if current:
        raise ValueError("history ends with an incomplete point frame")
    if previous_time is None:
        raise ValueError("history contains no complete frames")


def plane_summary(moment, station_ids):
    if moment.n < 1:
        return None
    summaries = []
    for indices in station_ids:
        means, variances, heterogeneity = [], [], []
        for component in range(3):
            point_means = [moment.mean[3 * i + component] for i in indices]
            point_variances = [max(0.0, moment.m2[3 * i + component] / moment.n) for i in indices]
            means.append(statistics.fmean(point_means))
            variances.append(statistics.fmean(point_variances))
            heterogeneity.append({"min": min(point_variances), "max": max(point_variances),
                                  "spatial_std": statistics.pstdev(point_variances)})
        summaries.append({"mean_velocity_m_s": means, "variance_m2_s2": variances,
                          "rms_m_s": [math.sqrt(v) for v in variances],
                          "K_m2_s2": 0.5 * sum(variances),
                          "point_variance_spatial_distribution_m2_s2": heterogeneity,
                          "spatial_std_temporal_mean_m_s": [statistics.pstdev(
                              [moment.mean[3 * i + c] for i in indices]) for c in range(3)]})
    return summaries


def t975(df):
    # Exact tabulated two-sided 95% t critical values through df=30; expansion
    # beyond that avoids a SciPy dependency (error <~2e-6 by df=31).
    table = [None, 12.7062047364, 4.30265272975, 3.18244630528, 2.77644510520,
             2.57058183564, 2.44691184879, 2.36462425159, 2.30600413503,
             2.26215716285, 2.22813885196, 2.20098516008, 2.17881282966,
             2.16036865646, 2.14478668792, 2.13144954556, 2.11990529922,
             2.10981557783, 2.10092204024, 2.09302405441, 2.08596344727,
             2.07961384473, 2.07387306790, 2.06865761042, 2.06389856163,
             2.05953855275, 2.05552943864, 2.05183051648, 2.04840714180,
             2.04522964213, 2.04227245630]
    if df <= 30:
        return table[df]
    z = statistics.NormalDist().inv_cdf(0.975)
    return (z + (z**3 + z) / (4 * df) + (5*z**5 + 16*z**3 + 3*z) / (96*df**2)
            + (3*z**7 + 19*z**5 + 17*z**3 - 15*z) / (384*df**3)
            + (79*z**9 + 776*z**7 + 1482*z**5 - 1920*z**3 - 945*z) / (92160*df**4))


def lag1(values):
    left, right = values[:-1], values[1:]
    ml, mr = statistics.fmean(left), statistics.fmean(right)
    numerator = sum((a - ml) * (b - mr) for a, b in zip(left, right))
    denominator = math.sqrt(sum((a - ml)**2 for a in left) * sum((b - mr)**2 for b in right))
    return numerator / denominator if denominator else None


def block_interval(values):
    mean = statistics.fmean(values)
    half_width = t975(len(values) - 1) * statistics.stdev(values) / math.sqrt(len(values))
    return {"block_values": values, "mean_of_blocks": mean,
            "conditional_95_percent_t_interval": [mean - half_width, mean + half_width],
            "lag1_block_correlation": lag1(values) if len(values) > 2 else None}


def relative_difference(a, b):
    denominator = max(abs(a), abs(b))
    return (b - a) / denominator if denominator else 0.0


def analyze(args):
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "NEKTAR_HIT_PILOT_V1":
        raise ValueError("manifest must have schema NEKTAR_HIT_PILOT_V1")
    probes = manifest["sampling"]["probes"]
    stations = manifest["sampling"]["stations_x_m"]
    if not probes or [probe["id"] for probe in probes] != list(range(len(probes))):
        raise ValueError("manifest probe IDs must be consecutive")
    station_ids = [[p["id"] for p in probes if p["station_index"] == station]
                   for station in range(len(stations))]
    if any(not ids for ids in station_ids):
        raise ValueError("manifest contains an empty probe plane")
    sample_dt = manifest["time"]["sample_dt_s"]
    start = manifest["time"]["record_start_s"] if args.start is None else args.start
    end = manifest["time"]["end_s"] if args.end is None else args.end
    if not all(math.isfinite(v) for v in (start, end, sample_dt)) or not 0 <= start < end or sample_dt <= 0:
        raise ValueError("require finite 0<=start<end and positive manifest sample interval")
    history = (args.history or manifest_path.parent / manifest["history_file"]).resolve()
    before = history.stat()
    file_hash = hashlib.sha256()
    overall = Moments(3 * len(probes))
    first_time = last_time = raw_first_time = raw_last_time = None
    raw_count = 0
    tolerance = 2e-12 * max(abs(start), abs(end), sample_dt)
    for time, values in history_frames(history, probes, sample_dt, file_hash):
        raw_count += 1
        raw_first_time = time if raw_first_time is None else raw_first_time
        raw_last_time = time
        if start - tolerance <= time <= end + tolerance:
            overall.add(values)
            first_time = time if first_time is None else first_time
            last_time = time
    expected_first = math.ceil((start - tolerance) / sample_dt) * sample_dt
    expected_last = math.floor((end + tolerance) / sample_dt) * sample_dt
    if (first_time is None or not close(first_time, expected_first, sample_dt)
            or not close(last_time, expected_last, sample_dt)):
        raise ValueError("history does not completely cover the requested sampled time window")
    if args.blocks < 2 or overall.n < 2 * args.blocks:
        raise ValueError("need blocks>=2 and at least two selected samples per block")
    per_block = overall.n // args.blocks
    used = per_block * args.blocks
    blocks = [Moments(3 * len(probes)) for _ in range(args.blocks)]
    windows = {name: Moments(3 * len(probes)) for name in
               ("first_half", "second_half", "last_75_percent")}
    window_times = {name: [] for name in windows}
    try:
        lags = sorted(set(int(part.strip()) for part in args.correlation_lags.split(",") if part.strip()))
    except ValueError as error:
        raise ValueError("correlation-lags must contain comma-separated integer sample lags") from error
    if any(lag < 1 or lag >= overall.n for lag in lags):
        raise ValueError("correlation lags must lie in [1,selected_sample_count)")
    ring = deque(maxlen=max(lags) if lags else 0)
    correlations = {lag: [[0.0] * 3 for _ in stations] for lag in lags}
    shape = manifest["sampling"]["probe_shape_yz"]
    ny, nz = shape
    spectra = None
    if args.z_spectra:
        if nz < manifest["mesh"]["hom_modes_z"]:
            raise ValueError("z-spectra requires probe Nz>=solver HomModesZ to avoid unreported spatial aliasing")
        expected = [(s, iy, iz) for s in range(len(stations)) for iy in range(ny) for iz in range(nz)]
        if [(p["station_index"], p["iy"], p["iz"]) for p in probes] != expected:
            raise ValueError("z-spectra requires canonical station,y,z probe order")
        weights = [[cmath.exp(-2j * math.pi * m * iz / nz) / nz for iz in range(nz)]
                   for m in range(nz // 2 + 1)]
        spectra = [[[0.0] * (nz // 2 + 1) for _ in range(3)] for _ in stations]
    selected = 0
    for time, values in history_frames(history, probes, sample_dt):
        if not start - tolerance <= time <= end + tolerance:
            continue
        if selected < used:
            blocks[selected // per_block].add(values)
        names = ["first_half" if selected < overall.n // 2 else "second_half"]
        if selected >= overall.n // 4:
            names.append("last_75_percent")
        for name in names:
            windows[name].add(values)
            if not window_times[name]:
                window_times[name].append(time)
            if len(window_times[name]) == 1:
                window_times[name].append(time)
            window_times[name][-1] = time
        centered = [value - mean for value, mean in zip(values, overall.mean)]
        for lag in lags:
            if selected >= lag:
                previous = ring[-lag]
                for station, ids in enumerate(station_ids):
                    for component in range(3):
                        correlations[lag][station][component] += sum(
                            centered[3*i+component] * previous[3*i+component] for i in ids) / len(ids)
        if lags:
            ring.append(centered)
        if spectra is not None:
            for station in range(len(stations)):
                for iy in range(ny):
                    first = (station * ny + iy) * nz
                    for component in range(3):
                        row = [centered[3*(first + iz)+component] for iz in range(nz)]
                        for mode, weight in enumerate(weights):
                            coefficient = sum(value * w for value, w in zip(row, weight))
                            factor = 1 if mode == 0 or (nz % 2 == 0 and mode == nz // 2) else 2
                            spectra[station][component][mode] += factor * abs(coefficient)**2
        selected += 1
    after = history.stat()
    if selected != overall.n or (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError("history changed while being analyzed; analyze a completed immutable record")
    summary = plane_summary(overall, station_ids)
    window_summaries = {name: {"samples": moment.n, "first_last_times_s": window_times[name],
                                "planes": plane_summary(moment, station_ids)}
                        for name, moment in windows.items()}
    block_results = []
    for station, ids in enumerate(station_ids):
        values = [[], [], []]
        own_center = [[], [], []]
        for block in blocks:
            for component in range(3):
                indices = [3 * i + component for i in ids]
                own_center[component].append(statistics.fmean(block.m2[i] / block.n for i in indices))
                values[component].append(statistics.fmean(
                    block.m2[i] / block.n + (block.mean[i] - overall.mean[i])**2 for i in indices))
        kinetic = [0.5 * sum(values[c][b] for c in range(3)) for b in range(args.blocks)]
        block_results.append({"station_x_m": stations[station],
                              "component_variance_m2_s2": [block_interval(v) for v in values],
                              "K_m2_s2": block_interval(kinetic),
                              "component_own_block_centered_variance_m2_s2": own_center})
    startup = []
    for station, x in enumerate(stations):
        first, second = (window_summaries[name]["planes"][station] for name in ("first_half", "second_half"))
        startup.append({"station_x_m": x,
                        "half_window_K_relative_change": relative_difference(first["K_m2_s2"], second["K_m2_s2"]),
                        "half_window_component_variance_relative_change": [relative_difference(a, b) for a, b in
                            zip(first["variance_m2_s2"], second["variance_m2_s2"])],
                        "last_half_vs_full_K_relative_change": relative_difference(summary[station]["K_m2_s2"], second["K_m2_s2"])})
    result = {
        "schema": "NEKTAR_HIT_PILOT_STATISTICS_V1",
        "inputs": {"history": str(history), "manifest": str(manifest_path),
                   "history_sha256": file_hash.hexdigest(),
                   "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
                   "analyzer_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "estimator": {"definition": "For each fixed probe: temporal mean and sum((u-mean)^2)/N; then arithmetic mean over plane probes. K=0.5*(var_u+var_v+var_w).",
                      "instantaneous_plane_mean_subtracted": False,
                      "population_variance_denominator": "N (matches native cumulative ReynoldsStresses)",
                      "probe_grid_integrates_plane_exactly": False,
                      "pressure": "parsed/validated, excluded from velocity statistics"},
        "window": {"requested_start_s": start, "requested_end_s": end,
                   "first_sample_s": first_time, "last_sample_s": last_time,
                   "sample_dt_s": sample_dt, "samples": overall.n,
                   "raw_history_frames": raw_count, "raw_first_last_times_s": [raw_first_time, raw_last_time]},
        "planes": [{"station_x_m": x, **entry} for x, entry in zip(stations, summary)],
        "blocks": {"count": args.blocks, "samples_per_block": per_block,
                   "sample_span_per_block_s": (per_block - 1) * sample_dt,
                   "block_start_separation_s": per_block * sample_dt,
                   "tail_samples_excluded_from_blocks": overall.n - used,
                   "definition": "Nonoverlapping equal-count contiguous blocks; component energy in each block is centered on that probe's full-window temporal mean. Own-block centered values are provided separately.",
                   "interval_assumptions": "Student-t interval for the mean of approximately independent stationary block energies. Nonoverlap alone does not establish independence; inspect temporal/block correlations and block-length sensitivity. The full-window mean is itself estimated. These are conditional sampling diagnostics, not experimental uncertainty.",
                   "planes": block_results},
        "startup_windows": window_summaries, "startup_comparison": startup,
        "limitations": ["Finite record and probe grid require duration/grid checks", "No experimental downstream data were loaded",
                        "No pass/fail thresholds or stationarity/independence claims are inferred automatically",
                        "At record-start=0 history includes initial state; native ReynoldsStresses begins at the first sampled solver step"]}
    if lags:
        result["temporal_correlations"] = {
            "definition": "Time-pair covariance about each probe's full-window temporal mean, spatially averaged, divided by that component's plane-averaged full-window variance; no instantaneous plane centering.",
            "lags": [{"lag_samples": lag, "lag_s": lag * sample_dt, "pairs_per_point": overall.n - lag,
                      "planes": [{"station_x_m": stations[s], "component_correlations": [
                          correlations[lag][s][c] / ((overall.n-lag) * summary[s]["variance_m2_s2"][c])
                          if summary[s]["variance_m2_s2"][c] else None for c in range(3)]}
                                 for s in range(len(stations))]} for lag in lags]}
    if spectra is not None:
        width = manifest["mesh"]["period_z_m"]
        for station in range(len(stations)):
            for component in range(3):
                spectra[station][component] = [value / (overall.n * ny) for value in spectra[station][component]]
        result["z_spectra"] = {
            "definition": "One-sided discrete variance per z Fourier mode, averaged over time and y; full-window per-point temporal mean removed before DFT. Mode0 is retained. Sum over modes equals sampled component temporal variance by Parseval.",
            "kz_rad_m": [2 * math.pi * m / width for m in range(nz // 2 + 1)],
            "probe_nz": nz,
            "planes": [{"station_x_m": stations[s], "component_variance_per_mode_m2_s2": spectra[s],
                        "parseval_relative_residual": [relative_difference(summary[s]["variance_m2_s2"][c], sum(spectra[s][c])) for c in range(3)]}
                       for s in range(len(stations))]}
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", nargs="?", type=Path, help="Defaults to manifest directory/native.his")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="JSON output; stdout if omitted")
    parser.add_argument("--start", type=float, help="Inclusive start, defaults to manifest record-start")
    parser.add_argument("--end", type=float, help="Inclusive end, defaults to manifest end")
    parser.add_argument("--blocks", type=int, default=8)
    parser.add_argument("--correlation-lags", default="", help="Optional comma-separated positive integer sample lags")
    parser.add_argument("--z-spectra", action="store_true", help="Optional stdlib DFT; requires probe Nz>=solver HomModesZ")
    args = parser.parse_args()
    try:
        analyze(args)
    except (ValueError, OSError, KeyError, TypeError, OverflowError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
