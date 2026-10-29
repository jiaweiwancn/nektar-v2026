#!/usr/bin/env python3
"""Join M7 temporal statistics with hashed native late-checkpoint slab work.

Standard-library postprocessing only. Instantaneous total resolved strain work
is not qualified fluctuation dissipation; three-snapshot ranges are not CIs.
"""

import argparse
from bisect import bisect_left
import csv
import datetime
import hashlib
import io
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET


TIMES = (1.30, 1.45, 1.60)
SLABS = 90
COLUMNS = ("xmin", "xmax", "xmid", "volume", "mean_u", "mean_v", "mean_w",
           "instantaneous_K_about_prescribed_mean", "molecular_strain_dissipation",
           "svv_xy", "svv_z", "divergence_rms")
WORK = ("molecular_strain_dissipation", "svv_xy", "svv_z", "total_svv_work",
        "total_svv_over_molecular", "divergence_rms")
SOURCE_ROOT = next((p for p in Path(__file__).resolve().parents
                    if (p/"CMakeLists.txt").is_file() and (p/"solvers").is_dir()), None)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_tree(value):
    if isinstance(value, float):
        require(math.isfinite(value), "nonfinite JSON number")
    elif isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


class Inputs:
    def __init__(self):
        self.files = {}

    def raw(self, path, expected=None):
        path = path.resolve()
        before = path.stat()
        raw = path.read_bytes()
        after = path.stat()
        require((before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
                f"input changed while reading: {path}")
        record = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
        require(expected is None or expected == record, f"receipt hash/size mismatch: {path}")
        require(str(path) not in self.files or self.files[str(path)] == record,
                f"input changed between reads: {path}")
        self.files[str(path)] = record
        return raw

    def read(self, path, expected=None):
        result = json.loads(self.raw(path, expected))
        finite_tree(result)
        return result


def close(a, b, rel=1e-9, absolute=1e-14):
    return math.isclose(a, b, rel_tol=rel, abs_tol=absolute)


def validate_output(path, case=None, source_root=SOURCE_ROOT):
    original = path.expanduser()
    output = original.resolve()
    require(not original.exists() and not original.is_symlink(), "output must be a fresh file; overwrite is forbidden")
    require(source_root is not None, "cannot discover a source checkout; copied scripts require --source-root")
    roots = [source_root.resolve()] + ([case.resolve()] if case is not None else [])
    require(all(output != root and root not in output.parents for root in roots),
            "output must be outside the source checkout and native case directory")
    require(output.parent.is_dir(), "output parent directory must already exist")
    return output


def range_summary(values):
    defined = [v for v in values if v is not None]
    return {"min": min(defined) if len(defined) == len(values) else None,
            "mean": math.fsum(defined) / len(values) if len(defined) == len(values) else None,
            "max": max(defined) if len(defined) == len(values) else None,
            "defined_snapshots": len(defined)}


def add_work(row):
    row = dict(row)
    row["total_svv_work"] = row["svv_xy"] + row["svv_z"]
    molecular = row["molecular_strain_dissipation"]
    row["total_svv_over_molecular"] = row["total_svv_work"] / molecular if molecular > 0 else None
    return row


def native_closure(rows, flow, speed, viscosity):
    volume = math.fsum(r["volume"] for r in rows)
    require(close(volume, flow["volume"], rel=1e-10), "slab/global volumes do not close")
    require(close(viscosity, flow["viscosity"], rel=2e-12, absolute=1e-18), "native/statistical viscosity mismatch")
    checks = {"volume": {"from_slabs": volume, "native_global": flow["volume"]}}
    for column, key in (("molecular_strain_dissipation", "instantaneous_resolved_strain_dissipation"),
                        ("svv_xy", "instantaneous_svv_work_xy"), ("svv_z", "instantaneous_svv_work_z")):
        value = math.fsum(r["volume"] * r[column] for r in rows) / volume
        require(close(value, flow[key]), f"slab/global work does not close: {column}")
        checks[column] = {"from_slabs": value, "native_global": flow[key], "difference": value-flow[key]}
    divergence = math.sqrt(math.fsum(r["volume"] * r["divergence_rms"]**2 for r in rows) / volume)
    require(close(divergence, flow["divergence_rms"]), "slab/global divergence RMS does not close")
    checks["divergence_rms"] = {"from_slabs": divergence, "native_global": flow["divergence_rms"]}
    norm_squared = math.fsum(r["volume"] * (2*r["instantaneous_K_about_prescribed_mean"]
                             + 2*speed*r["mean_u"] - speed**2) for r in rows)
    require(norm_squared >= 0 and close(math.sqrt(norm_squared), flow["velocity_l2"]),
            "slab fixed-mean energy/global velocity norm does not close")
    checks["velocity_l2"] = {"from_slabs": math.sqrt(norm_squared), "native_global": flow["velocity_l2"]}
    return checks


def load_profile(inputs, profile, manifest, speed, viscosity):
    expected = profile["output_files"]
    csv_paths = [Path(p) for p in expected if p.endswith(".slabs.csv")]
    flow_paths = [Path(p) for p in expected if p.endswith(".flow.json")]
    require(len(csv_paths) == len(flow_paths) == 1, "expected one slab CSV and flow JSON per checkpoint")
    flow_path, csv_path = flow_paths[0], csv_paths[0]
    flow = inputs.read(flow_path, expected[str(flow_path)])
    require(flow == profile["flow"], "embedded flow values differ from the hashed native JSON")
    raw = inputs.raw(csv_path, expected[str(csv_path)])
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
    require(tuple(reader.fieldnames or ()) == COLUMNS, "unexpected native slab CSV columns")
    rows = [{key: float(value) for key, value in row.items()} for row in reader]
    require(len(rows) == profile["slabs"] == manifest["mesh"]["nx"] == SLABS,
            "expected all 90 native x slabs")
    require(all(math.isfinite(v) for row in rows for v in row.values()), "nonfinite native slab CSV")
    mesh = manifest["mesh"]
    span = mesh["nx"] * mesh["hx_m"]
    for index, row in enumerate(rows):
        require(row["xmax"] > row["xmin"] and row["volume"] > 0
                and close(row["xmid"], .5*(row["xmin"]+row["xmax"]), rel=1e-12), "invalid slab geometry")
        require(row["molecular_strain_dissipation"] >= 0 and row["divergence_rms"] >= 0
                and row["instantaneous_K_about_prescribed_mean"] >= 0, "negative native squared quantity")
        require(close(row["volume"], (row["xmax"]-row["xmin"])*mesh["period_y_m"]*mesh["period_z_m"], rel=1e-10),
                "slab is not a complete full-y/full-z column")
        if index:
            require(row["xmin"] > rows[index-1]["xmin"] and close(row["xmin"], rows[index-1]["xmax"], rel=1e-10),
                    "unordered, overlapping or missing native slabs")
    require(close(rows[0]["xmin"], 0) and close(rows[-1]["xmax"], span, rel=1e-10),
            "slabs do not cover the complete streamwise domain")
    checks = native_closure(rows, flow, speed, viscosity)
    checkpoint_files = profile["checkpoint_inputs"]
    info_paths = [Path(p) for p in checkpoint_files if Path(p).name == "Info.xml"]
    require(len(info_paths) == 1, "expected one checkpoint metadata Info.xml")
    info = ET.fromstring(inputs.raw(info_paths[0], checkpoint_files[str(info_paths[0])]))
    require(close(float(info.findtext("Metadata/Time")), profile["time_s"], rel=2e-12),
            "checkpoint metadata/receipt times differ")
    require(info.findtext("Metadata/Provenance/GitSHA1") == profile["field_build_revision"],
            "checkpoint metadata/receipt build revisions differ")
    rows = [add_work(row) for row in rows]
    volume = flow["volume"]
    spatial = {}
    for key in WORK:
        values = [row[key] for row in rows]
        spatial[key] = {"min": min(values) if None not in values else None,
                        "volume_weighted_mean": math.fsum(row["volume"]*row[key] for row in rows)/volume
                            if None not in values else None,
                        "max": max(values) if None not in values else None}
    global_work = add_work({"molecular_strain_dissipation": flow["instantaneous_resolved_strain_dissipation"],
                           "svv_xy": flow["instantaneous_svv_work_xy"], "svv_z": flow["instantaneous_svv_work_z"],
                           "divergence_rms": flow["divergence_rms"]})
    return {"time_s": profile["time_s"], "checkpoint_index": profile["index"],
            "field_build_revision": profile["field_build_revision"], "rows": rows,
            "flow": flow, "global_work": global_work, "closure": checks, "spatial_summary": spatial,
            "checkpoint_input_hashes_from_receipt": checkpoint_files}


def interpolate(rows, x):
    centers = [row["xmid"] for row in rows]
    require(rows[0]["xmin"] <= x <= rows[-1]["xmax"], "station is outside native slab domain")
    hi = bisect_left(centers, x)
    if hi < len(rows) and centers[hi] == x:
        return ({key: rows[hi][key] for key in WORK}, {
            "method": "matching slab center; still a full x-slab average, not plane work",
            "source_slab_centers_m": [x], "exact_inlet_value": False})
    if hi == 0 or hi == len(rows):
        index = 0 if hi == 0 else len(rows)-1
        return ({key: rows[index][key] for key in WORK}, {
            "method": "nearest slab center at domain edge; no extrapolation",
            "source_slab_centers_m": [centers[index]], "exact_inlet_value": False,
            "note": ("The inlet uses the first full x-slab average, not exact inlet-plane work" if index == 0
                     else "The outlet uses the last full x-slab average, not exact outlet-plane work")})
    fraction = (x-centers[hi-1])/(centers[hi]-centers[hi-1])
    result = {key: rows[hi-1][key]+fraction*(rows[hi][key]-rows[hi-1][key]) for key in WORK
              if key != "total_svv_over_molecular"}
    result = add_work(result)
    return result, {"method": "linear interpolation of slab-center work; not exact plane gradients",
                    "source_slab_centers_m": [centers[hi-1], centers[hi]], "fraction": fraction,
                    "exact_inlet_value": False}


def scales(variance_u, molecular, viscosity):
    if molecular <= 0:
        return {"lambda_surrogate_m": None, "eta_surrogate_m": None, "Re_lambda_surrogate": None,
                "unavailable_reason": "native molecular work is not positive"}
    lam = math.sqrt(15*viscosity*variance_u/molecular)
    return {"lambda_surrogate_m": lam, "eta_surrogate_m": (viscosity**3/molecular)**.25,
            "Re_lambda_surrogate": math.sqrt(variance_u)*lam/viscosity}


def decay_intervals(planes, speed):
    result = []
    for first, second in zip(planes[:-1], planes[1:]):
        dx = second["station_x_m"]-first["station_x_m"]
        require(dx > 0, "temporal-statistics station coordinates must be unique and increasing")
        result.append({"x_interval_m": [first["station_x_m"], second["station_x_m"]],
            "interval_midpoint_m": .5*(first["station_x_m"]+second["station_x_m"]),
            "minus_U_delta_K_over_delta_x_m2_s3": -speed*(second["K_m2_s2"]-first["K_m2_s2"])/dx,
            "minus_1p5_U_delta_varu_over_delta_x_m2_s3":
                -1.5*speed*(second["variance_m2_s2"][0]-first["variance_m2_s2"][0])/dx})
    return result


def collect(args):
    source_root = getattr(args, "source_root", None) or SOURCE_ROOT
    validate_output(args.output, source_root=source_root)
    inputs = Inputs()
    analysis = inputs.read(args.analysis)
    require(analysis["schema"] == "NEKTAR_HIT_M7_ANALYSIS_V1", "unsupported M7 analysis schema")
    case = Path(analysis["case_directory"]).resolve()
    validate_output(args.output, case, source_root)
    manifest_path = Path(analysis["manifest_path"]).resolve()
    manifest = inputs.read(manifest_path)
    manifest_hash = inputs.files[str(manifest_path)]["sha256"]
    require(manifest["schema"] == "NEKTAR_HIT_PILOT_V1"
            and analysis["sha256"][str(manifest_path)] == manifest_hash
            and analysis["cache_receipt"]["metadata"]["manifest_sha256"] == manifest_hash,
            "M7 analysis/historical cache receipt/manifest digest mismatch")
    require(case == manifest_path.parent and analysis["physics"] == manifest["physics"]
            and analysis["synthesis"] == manifest["synthesis"]
            and analysis["manifest_seed"] == manifest["synthesis"]["seed_uint64"], "M7 case/physics/synthesis mismatch")
    viscosity, speed = (manifest["physics"][key] for key in ("kinvis_m2_s", "mean_velocity_m_s"))
    require(all(math.isfinite(v) and v > 0 for v in (viscosity, speed)), "positive finite viscosity/speed required")
    snapshots, ignored_times = [], []
    for path in args.profiles_receipts:
        receipt = inputs.read(path)
        require(receipt["schema"] == "NEKTAR_HIT_CHECKPOINT_PROFILES_V1" and receipt["status"] == "completed"
                and Path(receipt["case"]).resolve() == case, "profile receipt is incomplete or belongs to another case")
        source_inputs = receipt["source_inputs"]
        require(source_inputs[str(manifest_path)]["sha256"] == manifest_hash, "profile source/analysis manifest mismatch")
        for source, expected in source_inputs.items():
            inputs.raw(Path(source), expected)
            require(source not in analysis["sha256"] or analysis["sha256"][source] == expected["sha256"],
                    "profile source/analysis source digest mismatch")
        for profile in receipt["profiles"]:
            time = profile["time_s"]
            if not any(close(time, target, rel=2e-12) for target in TIMES):
                ignored_times.append(time)
                continue
            window = analysis["window"]
            require(window["requested_start_s"]-1e-12 <= time <= window["requested_end_s"]+1e-12,
                    "native snapshot is outside the temporal recording window")
            snapshots.append(load_profile(inputs, profile, manifest, speed, viscosity))
    snapshots.sort(key=lambda item: item["time_s"])
    require(len(snapshots) == 3 and all(close(s["time_s"], t, rel=2e-12) for s, t in zip(snapshots, TIMES)),
            "require exactly one completed snapshot at each of 1.30, 1.45 and 1.60 s")
    grid = [[row[key] for row in snapshots[0]["rows"]] for key in ("xmin", "xmax", "xmid", "volume")]
    require(all([[row[key] for row in s["rows"]] for key in ("xmin", "xmax", "xmid", "volume")] == grid
                for s in snapshots), "native checkpoint slab grids differ")
    profile = {"xmin_m": grid[0], "xmax_m": grid[1], "xmid_m": grid[2], "volume_m3": grid[3],
               "columns": {key: [range_summary([s["rows"][i][key] for s in snapshots]) for i in range(SLABS)] for key in WORK}}
    planes = sorted(analysis["statistics"]["planes"], key=lambda p: p["station_x_m"])
    require(planes, "M7 analysis contains no temporal station statistics")
    station_output = []
    published = {p["matched_station_x_m"]: p["table4_published_SI"] for p in analysis.get("Table2a_Table4_comparisons", [])
                 if p["matched_station_x_m"] is not None}
    for plane in planes:
        variance = plane["variance_m2_s2"]
        require(len(variance) == 3 and all(v >= 0 for v in variance) and plane["K_m2_s2"] >= 0
                and close(plane["K_m2_s2"], .5*math.fsum(variance), rel=1e-12), "invalid temporal variance/K")
        x = plane["station_x_m"]
        values = []
        for snapshot in snapshots:
            native, interpolation = interpolate(snapshot["rows"], x)
            values.append({"time_s": snapshot["time_s"], "work": native,
                           "instantaneous_total_gradient_surrogates": scales(variance[0], native["molecular_strain_dissipation"], viscosity)})
        scale_keys = ("lambda_surrogate_m", "eta_surrogate_m", "Re_lambda_surrogate")
        station_output.append({"station_x_m": x, "temporal_variance_u_m2_s2": variance[0],
            "temporal_K_m2_s2": plane["K_m2_s2"], "native_spatial_sampling": interpolation,
            "snapshots": values, "native_work_snapshot_ranges": {key: range_summary([v["work"][key] for v in values]) for key in WORK},
            "instantaneous_total_gradient_surrogate_ranges": {key: range_summary([v["instantaneous_total_gradient_surrogates"][key] for v in values]) for key in scale_keys},
            "published_Table4_SI_from_analysis": published.get(x)})
    require(str(args.output.resolve()) not in inputs.files, "output must not overwrite an input")
    return {"schema": "NEKTAR_HIT_M7_NATIVE_WORK_V1", "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "case_id": analysis["case_id"], "case_family": analysis["case_family"], "case_directory": str(case),
        "manifest_seed": analysis["manifest_seed"], "statistics_window": analysis["window"],
        "kinvis_m2_s": viscosity, "prescribed_mean_speed_m_s": speed, "snapshot_times_s": [s["time_s"] for s in snapshots],
        "script": {"path": str(Path(__file__).resolve()), "sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
        "verified_input_files": inputs.files, "ignored_other_snapshot_times_s": ignored_times,
        "snapshots": [{key: value for key, value in s.items() if key != "rows"} for s in snapshots],
        "slab_profiles_three_snapshot_ranges": profile, "stations": station_output,
        "adjacent_sampled_station_decay_intervals": decay_intervals(planes, speed),
        "published_Table4_matched_station_decay_intervals": decay_intervals([p for p in planes if p["station_x_m"] in published], speed),
        "definitions": {"molecular": "Native instantaneous total resolved 2*nu*S:S volume average; not qualified temporal fluctuation dissipation",
            "work_units": "m2/s3; total_svv_over_molecular is dimensionless; divergence_rms is 1/s",
            "ranges": "Arithmetic mean/min/max over exactly three instantaneous snapshots at each fixed slab/station, not a time integral or confidence interval",
            "spatial_summary": "Per-snapshot spatial min/max and slab-volume-weighted arithmetic mean. Mean of local SVV/molecular ratios differs from global_work's ratio of global means. Global divergence RMS uses the separate square-integral closure",
            "station_work": "Interpolate full x-slab averages by slab centers; domain edges use nearest center. Inlet values are not exact inlet-plane work",
            "lambda_surrogate": "sqrt(15*nu*temporal_var_u/native_instantaneous_molecular_work)",
            "eta_surrogate": "(nu^3/native_instantaneous_molecular_work)^(1/4)",
            "Re_lambda_surrogate": "sqrt(temporal_var_u)*lambda_surrogate/nu",
            "scales": "Instantaneous-total-gradient surrogates mixing full-record temporal var(u) with a single total-strain snapshot; no qualified fluctuation epsilon or Table4 fidelity claim",
            "decay": "-U*delta_K/delta_x and -1.5*U*delta_var_u/delta_x are interval secants; the second invokes Eq.38 isotropy. Unclosed budgets omit pressure, transport, nonuniform-mean and time effects. No error bars added",
            "hash_scope": "Analysis, manifest, source session/table, receipt, native CSV/flow JSON and checkpoint Info.xml rehashed. Cache binding uses the historical analysis/cache receipt; NPZ and native history are not freshly rehashed by this collector. Large field partitions and executable retain recorded receipt hashes and were not reread"},
        "closure_tolerances": {"volume_relative": 1e-10, "work_divergence_velocity_relative": 1e-9, "absolute": 1e-14},
        "claim_limits": ["Three late snapshots are descriptive and do not qualify mean fluctuation dissipation or steady energy-budget closure",
            "Published Table4 values are copied only from the existing analysis for context; no fitting or additional experimental source was loaded",
            "Zero molecular work makes its ratio/scales unavailable; incomplete defined ranges remain null rather than dropping snapshots"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--profiles-receipts", type=Path, nargs="+", required=True)
    parser.add_argument("--source-root", type=Path, help="Protected checkout root; required for a copied script outside a checkout")
    parser.add_argument("--output", type=Path, required=True, help="Fresh JSON file outside the native case and source checkout")
    args = parser.parse_args()
    try:
        result = collect(args)
        finite_tree(result)
        with args.output.expanduser().resolve().open("x", encoding="utf-8") as destination:
            destination.write(json.dumps(result, indent=2, allow_nan=False)+"\n")
    except (ValueError, OSError, KeyError, TypeError, OverflowError, ZeroDivisionError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
