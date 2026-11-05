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
    directions = [r for r in records if r["event"] == "direction"]
    require(directions, "the reduced direction was not constructed")
    require(any(d["qp_solved"] and d["displacement"] > 0 for d in directions),
            "no nonzero constrained model direction was solved")
    require(all(d["displacement"] <= d["radius"]*(1+1e-12) for d in directions),
            "final correction exceeds its radius")
    for direction in directions:
        for name in direction["basis_names"]:
            matching = [c for c in checks if c["name"] == name and
                        c["iteration"] == direction["iteration"] and
                        c["seconds"] <= direction["seconds"]]
            require(matching and matching[-1]["usable"] and
                    matching[-1]["relative_error"] <= .03 and
                    matching[-1]["observation_relative_error"] <= .05,
                    "a direction used an unchecked or inaccurate response")
    trials = [r for r in records if r["event"] == "trial"]
    require(trials, "no constrained direction received an equilibrated native trial")
    inward = [r for r in records if r["event"] == "ground_inward_roundoff"]
    require(all(r["after"]["admissible"] and r["displacement"] <= r["radius"]*(1+1e-12)
                for r in inward), "ground-inward correction violated the ground or radius bound")
    responses = [r for r in records if r["event"] == "non_k_response"]
    require(all(r["attained"] <= r["required"] for r in responses if r["converged"]),
            "a partially converged response was accepted")
    if report["steps"] == 0:
        startup = load(directory / "startup/case.evolution.json")
        require(state_hash(report["state"]) == H0 and report["time"] == 0 and
                report["native_history"] == startup["native_history"] and
                report["native_history_times"] == startup["native_history_times"],
                "reduced failure advanced the native state or history")
    result = dict(status="passed_operator_checks", manifest_status=manifest["status"],
                  repeatability=repeat, action_checks=checks, directions=directions,
                  trials=trials)
    write_json(directory / "reduced-checks.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    check(parser.parse_args().directory)
