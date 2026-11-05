#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check controlled reduced evaluations and the constrained direction diagnostics."""

import argparse
import json
from pathlib import Path

from run_cls_epawf_picard import require, load, state_hash, H0, write_json


def check(directory):
    manifest = load(directory / "manifest.json")
    report = load(directory / "target/case.evolution.json")
    records = []
    for line in (directory / "target/solver.log").read_text().splitlines():
        if line.startswith("EP-AWF Reduced "):
            records.append(json.loads(line[len("EP-AWF Reduced "):]))
    repeat, = [r for r in records if r["event"] == "repeatability"]
    require(repeat["failed_response_rolled_back"] and repeat["k_coefficients_fixed"],
            "controlled response rollback/fixed-k assertion failed")
    require(repeat["same_k_difference"] < 1e-10 and repeat["tightened_k_difference"] < 2e-6,
            "repeated reduced residual is inconsistent with the inner tolerance")
    checks = [r for r in records if r["event"] == "response_check"]
    require(len(checks) >= 3 and sum(c["usable"] for c in checks) >= 3,
            "insufficient independently consistent reduced responses")
    warm_starts = [r for r in records if r["event"] == "half_response_warm_start"]
    if any(r.get("half_response_warm_start") for r in records if r["event"] == "protocol"):
        require(warm_starts, "the half-response warm start was not checked independently")
    require(all(r["attained"] <= r["required"] and
                r["canonical_attained"] <= r["required"] and
                r["same_local_response"] == (r["relative_action_difference"] <= .01 and
                                             r["observation_action_difference"] <= .01)
                for r in warm_starts), "invalid half-response warm-start status")
    directions = [r for r in records if r["event"] == "direction"]
    require(directions, "the reduced direction was not constructed")
    require(any(d["qp_solved"] and d["displacement"] > 0 for d in directions),
            "no nonzero constrained model direction was solved")
    require(all(d["displacement"] <= d["radius"]*(1+1e-12) for d in directions),
            "final correction exceeds its radius")
    require(all(d["orthogonal_probes"] <= d["basis_dimension"] and
                (d["orthogonal_probes"] != d["basis_dimension"] or
                 d["maximum_column_correlation"] < 1e-10)
                for d in directions if "orthogonal_probes" in d),
            "orthogonal field probes retained correlated model columns")
    require(all(d["radius_model_solves"] >= 1 and d["radius_multiplier"] >= 0. and
                d["radius_complementarity"] >= 0. for d in directions
                if "radius_model_solves" in d and d["qp_solved"]),
            "invalid radius-constrained model status")
    for direction in directions:
        for name in direction["basis_names"]:
            matching = [c for c in checks if c["name"] == name and
                        c["iteration"] == direction["iteration"] and
                        c["seconds"] <= direction["seconds"]]
            require(matching and matching[-1]["usable"] and
                    matching[-1]["relative_error"] <= .03 and
                    matching[-1]["observation_relative_error"] <= .05,
                    "a direction used an unchecked or inaccurate response")
            checked_warm_starts = [r for r in warm_starts if r["name"] == name and
                                  r["iteration"] == direction["iteration"] and
                                  r["seconds"] <= direction["seconds"]]
            require(not checked_warm_starts or checked_warm_starts[-1]["same_local_response"],
                    "a direction used an inconsistent warm-start branch")
    trials = [r for r in records if r["event"] == "trial"]
    require(trials, "no constrained direction received an equilibrated native trial")
    corrections = [r for r in records if r["event"] == "second_order_trial"]
    require(all(r["correction"] <= r["correction_radius"]*(1+1e-12) and
                r["displacement"] <= r["radius"]*(1+1e-12) and
                (not r["valid"] or r["non_k_error"] <= r["required"])
                for r in corrections),
            "nonlinear qualification correction violated a field bound or H accuracy")
    require(all(r["violation"] < r["before_violation"] for r in corrections if r["improved"]),
            "a failed nonlinear qualification correction was marked useful")
    combined = [r for r in records if r["event"] == "combined_direction_check"]
    require(all(r["attained"] <= r["required"] and r["half_attained"] <= r["required"] and
                r["half_relative_residual_error"] <= .03 and
                r["half_observation_relative_error"] <= .05
                for r in combined if r["available"]),
            "combined response does not match its checked reduced model")
    inward = [r for r in records if r["event"] == "ground_inward_roundoff"]
    require(all(r["after"]["admissible"] and r["displacement"] <= r["radius"]*(1+1e-12)
                for r in inward), "ground-inward correction violated the ground or radius bound")
    responses = [r for r in records if r["event"] == "non_k_response"]
    require(all(r["attained"] <= r["required"] for r in responses if r["converged"]),
            "a partially converged response was accepted")
    if report["steps"] == 0:
        startup = load(directory / "startup/case.evolution.json")
        require(state_hash(report["state"]) == H0 and report["time"] == 0 and
                [state_hash(slot) for slot in report["native_history"]] ==
                [state_hash(slot) for slot in startup["native_history"]] and
                report["native_history_times"] == startup["native_history_times"],
                "reduced failure advanced the native state or history")
    result = dict(status="passed_operator_checks", manifest_status=manifest["status"],
                  repeatability=repeat, half_response_warm_starts=warm_starts,
                  action_checks=checks, directions=directions,
                  trials=trials, combined_direction_checks=combined,
                  second_order_corrections=corrections)
    write_json(directory / "reduced-checks.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    check(parser.parse_args().directory)
