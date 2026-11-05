#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect compact native reduced-search evidence, never infer M8 convergence."""
import argparse
from collections import Counter
import json
from pathlib import Path
import struct

from run_cls_epawf_picard import load, require, sha, state_hash, write_json, H0


def records(directory):
    path = directory / "target/solver.log"
    if not path.exists():
        return []
    return [json.loads(line[15:]) for line in path.read_text().splitlines()
            if line.startswith("EP-AWF Reduced ")]


def compact(row):
    result = dict(row)
    if "active" in result:
        result["active"] = [{k: v for k, v in cut.items() if k != "row"}
                            for cut in result["active"]]
    for key in ("shift_domain", "half_domain"):
        if key in result:
            domain = result.pop(key)
            result[key+"_summary"] = {k: v for k, v in domain.items()
                                       if k not in ("profile_branch_pattern", "limiting_point")}
    return result


def collect(directory):
    manifest = load(directory / "manifest.json")
    rows = records(directory)
    responses = [r for r in rows if r["event"] == "non_k_response"]
    require(all(r["attained"] <= r["required"] for r in responses if r["converged"]),
            "a converged response misses its declared H tolerance")
    histories = [compact(r) for r in rows if r["event"] in
                 ("protocol", "initial", "accepted", "converged", "repeatability",
                  "direction", "trial", "model_retry", "prediction_decomposition",
                  "inaccurate_or_failed_responses", "trial_packet", "model_packet",
                  "ground_cut", "ground_inward_roundoff", "restoration_balance")]
    checks = [compact(r) for r in rows if r["event"] == "response_check"]
    for row in checks:
        if row["iteration"]:
            row.pop("signed_observations", None)
    outcome = dict(path=str(directory), manifest=manifest, history=histories,
                   response_checks=checks,
                   non_k=dict(evaluations=len(responses),
                              converged=sum(r["converged"] for r in responses),
                              failed=Counter(r["failure"] for r in responses if not r["converged"]),
                              flow_solves=sum(r["flow_solves"] for r in responses),
                              epsilon_solves=sum(r["epsilon_solves"] for r in responses),
                              maximum_attained=max((r["attained"] for r in responses if r["converged"]), default=None),
                              maximum_attained_ratio=max((r["attained"]/r["required"] for r in responses
                                                          if r["converged"]), default=None)))
    # Old development controls used a build path subsequently rebuilt. Keep
    # their launch identities and report that distinction; current immutable
    # image runs must still match their launch executable and libraries.
    outcome["runtime_at_collection"] = {}
    for name in ("IncNavierStokesSolver", "RANSWallResolutionProbe"):
        launched = manifest.get(name)
        if not launched:
            continue
        paths = {launched["path"]: launched["sha256"], **launched["libraries"]}
        current = {path: sha(Path(path)) if Path(path).is_file() else None for path in paths}
        differences = [path for path in paths if current[path] != paths[path]]
        immutable = "/binaries-" in launched["path"]
        require(not immutable or not differences,
                "an immutable native image/runtime changed since launch")
        outcome["runtime_at_collection"][name] = dict(sha256=current,
            differences_from_launch=differences, immutable_image=immutable)
    report = directory / "target/case.evolution.json"
    if report.exists():
        state = load(report)
        outcome["native"] = {k: state.get(k) for k in ("steps", "time", "failure")}
        if state["steps"] == 0:
            startup = load(directory / "startup/case.evolution.json")
            unchanged = (state_hash(state["state"]) == H0 and state["time"] == 0 and
                         state["native_history"] == startup["native_history"] and
                         state["native_history_times"] == startup["native_history_times"])
            require(unchanged, "failed native solve changed accepted state or physical history")
            outcome["native"]["accepted_state_and_history_preserved"] = True
        readback = manifest.get("readback", {}).get("original_block_norms")
        accepted = [r for r in rows if r["event"] in ("initial", "accepted", "converged")]
        if readback and accepted:
            last = accepted[-1]["norms"]
            discrepancy = max(abs(a-b) for a, b in zip(last, readback))
            refined = any(r["purpose"] == "refine_base_accuracy" and
                          r["seconds"] > accepted[-1]["seconds"] for r in responses)
            require(refined or discrepancy < 1e-9,
                    "fresh original residual differs from the terminal reduced state")
            outcome["native"]["readback_norm_difference"] = discrepancy
            outcome["native"]["H_refined_after_last_logged_state"] = refined
    else:
        outcome["native"] = dict(status="no_native_report", accepted_physical_steps=None,
                                 note="Do not infer native rejection/history behavior from a killed process.")
    outcome["affine_packets"] = []
    for path in sorted((directory / "target").glob("reduced-model-*.bin")):
        with path.open("rb") as stream:
            n, m = struct.unpack("<QQ", stream.read(16))
        require(path.stat().st_size == 16+8*(n*n+2*n+m*(n+1)), "affine packet shape differs")
        outcome["affine_packets"].append(dict(path=str(path), variables=n, cuts=m,
                                               bytes=path.stat().st_size, sha256=sha(path)))
    checked = directory / "reduced-checks.json"
    if checked.exists():
        outcome["operator_checks"] = load(checked)["status"]
    return outcome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_root", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--case", action="append", dest="cases")
    args = parser.parse_args()
    cases = args.cases or [p.name for p in sorted(args.run_root.iterdir())
                           if p.is_dir() and (p / "manifest.json").exists()]
    result = dict(scope="Native operator and original residual evidence; statuses retain failed stage results.",
                  run_root=str(args.run_root), cases={n: collect(args.run_root / n) for n in cases})
    write_json(args.output, result)


if __name__ == "__main__":
    main()
