#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check a completed run that retains an independent previous-k direction.

This is a focused operator/history regression for a bounded single-stage
fixture. Passing does not certify original M8 convergence or physical evolution.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import test_cls_epawf_reduced_response as operator_checks
from run_cls_epawf_picard import load, require, state_hash


def records(run):
    return [json.loads(line[15:], parse_int=float)
            for line in (run / "target/solver.log").read_text().splitlines()
            if line.startswith("EP-AWF Reduced ")]


def without_seconds(value):
    if isinstance(value, dict):
        return {key: without_seconds(item) for key, item in value.items() if key != "seconds"}
    if isinstance(value, list):
        return [without_seconds(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--control", type=Path,
                        help="completed matched run with the previous direction queued last")
    args = parser.parse_args()
    subprocess.run([sys.executable, operator_checks.__file__, str(args.run)], check=True)
    rows = records(args.run)
    protocol, = [r for r in rows if r["event"] == "protocol"]
    require(protocol.get("previous_direction_priority") == "after_conditional_krylov_and_rate_seeds",
            "run does not exercise the repaired previous-direction priority")
    directions = [r for r in rows if r["event"] == "direction" and r["iteration"] > 0
                  and "previous" in r["basis_names"]]
    require(directions, "fixture did not retain an independent previous-k response")
    selected = []
    for direction in directions:
        iteration = direction["iteration"]
        checks = [r for r in rows if r["event"] == "response_check" and
                  r["iteration"] == iteration and r["name"] == "previous" and
                  r["seconds"] < direction["seconds"]]
        require(checks and checks[-1]["usable"] and checks[-1]["independent_component"] > 0,
                "retained previous direction lacks its independent feasible two-size response")
        require(any(r["event"] == "previous_direction_seed" and r["iteration"] == iteration
                    for r in rows), "previous response was not explicitly queued")
        selected.append(dict(model_seconds=direction["seconds"], **{key: checks[-1][key] for key in
                         ("iteration", "h", "relative_error", "observation_relative_error",
                          "signed_merit_slope", "half_merit_slope", "independent_component")}))
    result = dict(status="passed_previous_direction_operator_checks_not_m8_convergence",
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  run=str(args.run), previous_responses=selected,
                  original_operator_checks=str(args.run / "reduced-checks.json"))
    output = args.run / "previous-direction-checks.json"
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        control = records(args.control)
        require((args.run / "target/case.xml").read_bytes() ==
                (args.control / "target/case.xml").read_bytes(), "matched inputs differ")
        for event in ("initial", "accepted"):
            points = [next(row for row in run if row["event"] == event)
                      for run in (rows, control)]
            require(points[0]["state_sha256"] == points[1]["state_sha256"] and
                    state_hash(points[0]["norms"]) == state_hash(points[1]["norms"]),
                    "unmodified first model changed the matched numerical state or norms")
        prefixes = [[row for row in run if row["event"] == "response_check" and
                     row["iteration"] == 0] for run in (rows, control)]
        serialized = [json.dumps(without_seconds(prefix), sort_keys=True) for prefix in prefixes]
        require(len(prefixes[0]) >= 3 and serialized[0] == serialized[1],
                "first-model response values differ beyond elapsed seconds")
        require(not any(row["event"] == "direction" and "previous" in row["basis_names"]
                        for row in control), "control already retained the previous direction")
        require(all(load(run / "target/case.evolution.json")["steps"] == 0
                    for run in (args.run, args.control)),
                "comparison fixture expects rejected native stages, not physical evolution")
        result["matched_control"] = dict(
            run=str(args.control), identical_target_xml=True,
            initial_and_first_update_state_and_norms_bitwise_identical=True,
            first_model_response_count=len(prefixes[0]),
            first_model_response_values_identical_except_seconds=True,
            first_model_response_digest=hashlib.sha256(serialized[0].encode()).hexdigest(),
            numerical_residual_histories={str(run): [row["reduced_norm"] for row in data
                                                    if row["event"] == "accepted"]
                                          for run, data in ((args.run, rows),
                                                            (args.control, control))})
        result["status"] = "passed_previous_direction_comparison_not_m8_convergence"
        output = args.run / "previous-direction-comparison.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
