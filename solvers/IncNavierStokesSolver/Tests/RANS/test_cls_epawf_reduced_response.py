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
    require(all(d["displacement"] <= d["radius"]*(1+1e-12) for d in directions),
            "final correction exceeds its radius")
    require(all(any(c["name"] == name and c["iteration"] == d["iteration"] and
                        c["usable"] and c["relative_error"] <= .03 and
                        c["observation_relative_error"] <= .05 for c in checks)
                for d in directions for name in d["basis_names"]),
            "a direction used an unchecked or inaccurate response")
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
                  trials=[r for r in records if r["event"] == "trial"])
    write_json(directory / "reduced-checks.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    check(parser.parse_args().directory)
