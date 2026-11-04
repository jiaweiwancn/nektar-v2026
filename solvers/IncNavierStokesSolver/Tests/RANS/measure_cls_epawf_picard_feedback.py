#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Measure native k responses after sequential flow and epsilon equilibration."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import platform
import xml.etree.ElementTree as ET

import numpy as np

from cls_epawf_m3_cases import set_parameter, write
from run_cls_epawf_picard import H0, load, prepare, require, run, sha, stages, state_hash, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("guess", type=Path)
    parser.add_argument("diagnostic", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    require(platform.system() == "Linux" and 1 <= args.workers <= 4,
            "use Linux and an explicit allocation of at most four ranks")
    args.solver = args.solver.resolve(strict=True)
    args.output = args.output.resolve()
    require(args.output.is_relative_to(Path("/media/jay/新加卷/nektar-v2026-runs")),
            "use the designated run storage")
    args.output.mkdir()
    packet = load(args.guess)
    data = json.loads(args.diagnostic.read_text())
    require(state_hash(packet["state"]) == data["guess_sha256"] and
            state_hash(packet["history"]) == H0, "diagnostic and historical history binding differ")
    rows = data["free_rows"]
    directions = {a["name"]: a["direction"] for a in data["actions"]}
    names = ("ground_inward", "sampled_K", "near_tangent_1")
    tolerance = 2e-9

    def measure(job):
        name, h = job
        case = args.output / (name if not h else f"{name}-{h:g}")
        case.mkdir()
        initial = dict(packet, state=list(packet["state"]))
        if h:
            for row, value in zip(rows, directions[name]):
                initial["state"][row] += h*value
        guess = case / "guess.json"
        write_json(guess, initial)
        prepare(argparse.Namespace(output=case, verify=False, steps=1, sweeps=55,
                                    flow=1, k=.75, epsilon=1, pseudo_alpha=0, perturbation=0,
                                    ground_projection=False, guess_file=guess))
        target = case / "target"
        root = ET.parse(target / "case.xml").getroot()
        set_parameter(root, "RANSWallPicardNonKTolerance", tolerance)
        write(root, target / "case.xml")
        receipt = run([str(args.solver), "case.xml"], target, "solver.log", 480)
        report = load(target / "case.evolution.json")
        records = stages(target)
        result = dict(name=name, h=h, receipt=receipt, failure=report["failure"])
        events = [row for row in records if row["event"] == "non_k_equilibrated"]
        if events:
            final, = events
            require(report["steps"] == report["time"] == 0 and
                    state_hash(report["state"]) == H0 and
                    all(report["rejected_state"][i] == initial["state"][i] for i in rows),
                    "non-k equilibration advanced history or changed a fixed k coefficient")
            require(final["non_k_norm"] <= tolerance, "non-k solve missed the diagnostic tolerance")
            result["final"] = final
        return result

    jobs = [("base", 0.)] + [(name, h) for name in names for h in (1e-4, 1e-5)]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        cases = list(pool.map(measure, jobs))
    summary = dict(solver=str(args.solver), solver_sha256=sha(args.solver),
                   diagnostic_sha256=sha(args.diagnostic), guess_sha256=sha(args.guess),
                   non_k_tolerance=tolerance, workers=args.workers, cases=cases)
    if all("final" in case for case in cases):
        mass = np.array(data["mass"])
        jac = np.array(data["k_jacobian"])
        residual = np.array(cases[0]["final"]["k_residual"])
        dual = np.linalg.solve(mass, residual)

        def dual_norm(value):
            return float(np.sqrt(value @ np.linalg.solve(mass, value)))

        for case in cases[1:]:
            fixed = jac @ np.array(directions[case["name"]])
            response = (np.array(case["final"]["k_residual"]) - residual) / case["h"]
            case["action"] = dict(conditional_norm=dual_norm(fixed),
                                  equilibrated_norm=dual_norm(response),
                                  feedback_norm=dual_norm(response-fixed),
                                  conditional_residual_projection=float(fixed @ dual/dual_norm(residual)),
                                  equilibrated_residual_projection=float(response @ dual/dual_norm(residual)),
                                  outer_gradient_difference=(case["final"]["outer_range"][0] -
                                                             cases[0]["final"]["outer_range"][0]))
        summary["status"] = "measured"
    else:
        summary["status"] = "incomplete"
    write_json(args.output / "feedback.json", summary)
    return 0 if summary["status"] == "measured" else 1


if __name__ == "__main__":
    raise SystemExit(main())
