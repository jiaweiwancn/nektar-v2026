#!/usr/bin/env python3
"""Compare completed pilot summaries; all thresholds and block CIs are advisory."""

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import statistics
import sys

sys.dont_write_bytecode = True
ANALYZER = Path(__file__).with_name("analyze_history.py")
SPEC = importlib.util.spec_from_file_location("hit_pilot_analysis", ANALYZER)
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def finite_tree(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError("nonfinite JSON number")
    if isinstance(value, dict):
        for item in value.values():
            finite_tree(item)
    elif isinstance(value, list):
        for item in value:
            finite_tree(item)


def read(path):
    result = json.loads(path.read_text(encoding="utf-8"))
    finite_tree(result)
    return result


def canonical(manifest):
    synthesis = manifest["synthesis"]
    source = sorted(value for key, value in manifest["sha256"].items()
                    if Path(key).suffix == ".dat" and "spectrum" in Path(key).name.lower())
    if not source:
        return None
    return {key: synthesis[key] for key in (
        "periods_m", "counts", "seed_uint64", "cutoff_rad_m", "taper_start_rad_m")} | {
        "source_spectrum_sha256": source,
        "spectrum_definition": {key: manifest["physics"][key] for key in (
            "mean_velocity_m_s", "normalization", "interpolation", "low_extension", "high_extension")}}


def quantities(plane):
    variances, rms = plane["variance_m2_s2"], plane["rms_m_s"]
    if len(variances) != 3 or len(rms) != 3 or min(variances + rms) < 0:
        raise ValueError("invalid component variance/RMS")
    if not all(math.isfinite(value) for value in variances + rms + [plane["K_m2_s2"]]):
        raise ValueError("nonfinite variance/RMS/K")
    if plane["K_m2_s2"] < 0:
        raise ValueError("negative kinetic energy")
    return dict(zip(("u", "v", "w", "K"), variances + [plane["K_m2_s2"]]))


def load_case(directory):
    paths = {name: directory / name for name in (
        "manifest.json", "statistics.json", "run.receipt.json")}
    manifest, stats, run = (read(paths[name]) for name in paths)
    flow = directory / (Path(manifest["case"]).stem + ".flow.json")
    paths[flow.name] = flow
    flow_metrics = read(flow)
    if run["solver"]["exit_code"] or run["solver"].get("native_iteration_limit_exceeded", False):
        raise ValueError("solver failed or exceeded its native iteration limit")
    hashes = {name: digest(path) for name, path in paths.items()}
    if stats["inputs"]["manifest_sha256"] != hashes["manifest.json"]:
        raise ValueError("statistics/manifest digest mismatch")
    planes = {plane["station_x_m"]: quantities(plane) for plane in stats["planes"]}
    blocks = {plane["station_x_m"]: plane for plane in stats["blocks"]["planes"]}
    if len(planes) != len(stats["planes"]) or len(blocks) != len(stats["blocks"]["planes"]):
        raise ValueError("duplicate station coordinates")
    for x in planes:
        if not math.isfinite(x) or x not in blocks:
            raise ValueError("invalid or missing block station")
        if len(blocks[x]["component_variance_m2_s2"]) != 3:
            raise ValueError("block station needs three velocity components")
        for block in blocks[x]["component_variance_m2_s2"] + [blocks[x]["K_m2_s2"]]:
            ci = block["conditional_95_percent_t_interval"]
            if len(ci) != 2 or ci[0] > ci[1] or len(block["block_values"]) != stats["blocks"]["count"]:
                raise ValueError("invalid conditional block interval/count")
    resources = {name: {key: run[name].get(key) for key in (
        "wall_seconds", "sampled_process_tree_peak_rss_kib", "sampled_process_tree_peak_threads")}
                 for name in ("initialization", "solver") if name in run}
    return {"manifest": manifest, "stats": stats, "planes": planes, "blocks": blocks,
            "canonical": canonical(manifest), "metadata": {
                "source_revision": run["revision"], "input_hashes": hashes,
                "statistics_inputs": stats["inputs"], "mesh": manifest["mesh"],
                "knobs": manifest["knobs"], "svv": manifest["svv"], "flow": flow_metrics,
                "window": stats["window"], "resources": resources}}


def eligibility(case, reference):
    reasons = []
    if case["canonical"] is None or reference["canonical"] is None:
        reasons.append("unknown canonical inflow source hash")
    elif case["canonical"] != reference["canonical"]:
        reasons.append("canonical inflow differs (including span/seed/source)")
    keys = ("requested_start_s", "requested_end_s", "first_sample_s", "last_sample_s", "sample_dt_s")
    # Only clock roundoff is tolerated; station-coordinate matching stays exact.
    same_window = all(math.isclose(case["stats"]["window"][key], reference["stats"]["window"][key],
                      rel_tol=2e-12, abs_tol=2e-14) for key in keys)
    same_window &= case["stats"]["window"]["samples"] == reference["stats"]["window"]["samples"]
    if not same_window:
        reasons.append("record/sample window, sample dt or sample count differs")
    for key in ("count", "samples_per_block", "tail_samples_excluded_from_blocks", "definition"):
        if case["stats"]["blocks"][key] != reference["stats"]["blocks"][key]:
            reasons.append("block layout/definition differs")
            break
    if case["stats"]["blocks"]["count"] < 2:
        reasons.append("fewer than two paired blocks")
    return not reasons, same_window, reasons


def compare(case, reference):
    paired, same_window, reasons = eligibility(case, reference)
    svv_changed = case["manifest"]["svv"] != reference["manifest"]["svv"]
    mesh_changed = any(case["manifest"]["mesh"][key] != reference["manifest"]["mesh"][key]
                       for key in ("hx_m", "hy_m", "polynomial_degree"))
    fourier_changed = any(case["manifest"]["mesh"][key] != reference["manifest"]["mesh"][key]
                          for key in ("hom_modes_z", "period_z_m"))
    coupled = (svv_changed or (mesh_changed and case["manifest"]["svv"]["model_xy"] == "dg")
               or (fourier_changed and case["manifest"]["svv"]["model_z"] != "none"))
    rows, changes = [], []
    for x in sorted(case["planes"].keys() & reference["planes"].keys()):
        row = {"station_x_m": x, "quantities": {}}
        for component, key in enumerate(("u", "v", "w", "K")):
            value, base = case["planes"][x][key], reference["planes"][x][key]
            change = 100 * (value - base) / base if base else (0.0 if value == 0 else None)
            if change is not None and not math.isfinite(change):
                change = None
            own_block, ref_block = (item["blocks"][x]["component_variance_m2_s2"][component]
                if component < 3 else item["blocks"][x]["K_m2_s2"] for item in (case, reference))
            own_ci, ref_ci = (block["conditional_95_percent_t_interval"] for block in (own_block, ref_block))
            result = {"case": value, "reference": base, "signed_relative_change_percent": change,
                      "conditional_block_95_intervals": {"case": own_ci, "reference": ref_ci,
                          "overlap": max(own_ci[0], ref_ci[0]) <= min(own_ci[1], ref_ci[1])},
                      "lag1_block_correlation": {"case": own_block.get("lag1_block_correlation"),
                                                 "reference": ref_block.get("lag1_block_correlation")}}
            if paired:
                differences = [a - b for a, b in zip(own_block["block_values"], ref_block["block_values"])]
                mean = statistics.fmean(differences)
                width = ANALYSIS.t975(len(differences) - 1) * statistics.stdev(differences) / math.sqrt(len(differences))
                result["paired_block_difference_case_minus_reference"] = {
                    "block_values": differences, "mean": mean,
                    "conditional_95_percent_t_interval": [mean - width, mean + width],
                    "lag1_block_correlation": ANALYSIS.lag1(differences) if len(differences) > 2 else None}
            row["quantities"][key] = result
            changes.append(abs(change) if change is not None else None)
        rows.append(row)
    maximum = max((value for value in changes if value is not None), default=None)
    return {"paired_blocks_eligible": paired, "unpaired_reasons": reasons,
            "same_record_window": same_window, "resolution_model_coupling": coupled,
            "svv_settings_changed": svv_changed, "stations": rows,
            "maximum_absolute_relative_change_percent": maximum,
            "informational_all_within_2_percent":
                maximum <= 2 if paired and changes and None not in changes else None,
            "no_common_stations": not rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--reference", default="h05")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    cases, skipped = {}, []
    study_path = root / "study.receipt.json"
    try:
        statuses = {item["id"]: item for item in read(study_path)["cases"]} if study_path.exists() else {}
        for directory in sorted(path for path in root.iterdir() if path.is_dir()):
            if not (directory / "manifest.json").exists() and directory.name not in statuses:
                continue
            try:
                if statuses.get(directory.name, {}).get("status") == "failed":
                    raise ValueError(statuses[directory.name].get("error", "study reports failure"))
                cases[directory.name] = load_case(directory)
            except (OSError, ValueError, KeyError, TypeError) as error:
                skipped.append({"id": directory.name, "reason": str(error)})
        seen = set(cases) | {item["id"] for item in skipped}
        skipped += [{"id": name, "reason": item.get("error", "study reports failure")}
                    for name, item in statuses.items() if item.get("status") == "failed" and name not in seen]
        if args.reference not in cases:
            raise ValueError(f"reference {args.reference!r} has no completed valid statistics: {skipped}")
        result = {"schema": "NEKTAR_HIT_PILOT_COMPARISON_V1", "root": str(root),
                  "reference": args.reference, "comparator_sha256": digest(Path(__file__)),
                  "study_receipt_sha256": digest(study_path) if study_path.exists() else None,
                  "t_interval_source_sha256": digest(ANALYZER), "skipped_cases": skipped,
                  "reference_metadata": cases[args.reference]["metadata"],
                  "cases": {name: {**case["metadata"], "canonical_inflow": case["canonical"],
                      "comparison": compare(case, cases[args.reference])} for name, case in cases.items()},
                  "interpretation": "Descriptive sensitivity only; the 2% flag makes no milestone closure decision. Paired differences use aligned time blocks of spatially averaged, full-window-centered energies; mesh/Hom/probe grids may differ. Student-t intervals remain conditional on approximate stationarity and independent blocks; nonoverlap does not establish independence. SVV and DG h/p changes couple resolution/model effects. No downstream experimental data were loaded."}
        args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
