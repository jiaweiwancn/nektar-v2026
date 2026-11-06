#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check actual native progress after an initially failed radius-model query.

This is a bounded solver recovery test, not a certificate of M8 convergence.
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
    return [json.loads(line[15:]) for line in (run / "target/solver.log").read_text().splitlines()
            if line.startswith("EP-AWF Reduced ")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--control", type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, operator_checks.__file__, str(args.run)], check=True)
    rows = records(args.run)
    used = []
    for position, row in enumerate(rows):
        if row["event"] != "radius_query_retry" or row["purpose"] != "main_direction" or \
                row["initial_failure"] != "quadratic_iteration_limit" or not row["solved"]:
            continue
        end = next((i for i in range(position+1, len(rows))
                    if rows[i]["event"] in ("model_retry", "accepted", "converged")), len(rows))
        following = rows[position+1:end]
        directions = [r for r in following if r["event"] == "direction" and r["qp_solved"] and
                      r["ground_model_admissible"]]
        signed = [r for r in following if r["event"] == "combined_direction_check" and r["available"] and
                  r["predicted_merit_slope"] < 0 and r["half_merit_slope"] < 0]
        trials = [r for r in following if r["event"] == "trial" and r["accepted"]]
        if not directions or not signed or not trials:
            continue
        before = [r for r in rows[:position] if r["event"] in ("initial", "accepted")][-1]
        require(trials[-1]["full_norm"] < max(before["norms"]),
                "model recovery did not lower the original full residual")
        used.append(dict(retry=row, direction={k:directions[-1][k] for k in
            ("basis_dimension", "response_rank", "displacement", "radius", "active_constraints")},
            signed_check={k:v for k,v in signed[-1].items() if k != "signed_active_observations"},
            actual_trial={k:v for k,v in trials[-1].items() if k != "qualification_inventory"},
            original_full_residual_before=max(before["norms"])))
    require(used, "fixture lacks an independently checked accepted direction after the failed initial query")
    native = load(args.run / "target/case.evolution.json")
    result = dict(status="passed_radius_query_recovery_checks_not_m8_convergence",
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  run=str(args.run), recovered_directions=used,
                  accepted_physical_steps=native["steps"], time=native["time"])
    if args.control:
        first = min(r["retry"]["iteration"] for r in used)
        def prefix(data):
            return [(r["state_sha256"],state_hash(r["norms"])) for r in data
                    if r["event"] in ("initial", "accepted") and r["iteration"] <= first]
        require(prefix(rows) == prefix(records(args.control)),
                "the original-start numerical prefix changes before recovery")
        result["matched_prefix"] = dict(control=str(args.control),
            before_iteration=first, original_states_and_norms_bitwise_identical=True,
            scope="Observed first-stage prefix; run budgets and producing executables are recorded separately")
    (args.run / "radius-query-checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
