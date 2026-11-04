#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare two bounded native sweeps with Picard verification on and off."""

import argparse
from pathlib import Path
import platform

from run_cls_epawf_picard import H0, prepare, require, run, stages, load, state_hash, sha, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    require(platform.system() == "Linux", "use the Linux validation build")
    args.solver = args.solver.resolve(strict=True)
    args.output = args.output.resolve()
    require(args.output.is_relative_to(Path("/media/jay/新加卷/nektar-v2026-runs")),
            "use the designated run storage")
    args.output.mkdir()
    result = {"solver": str(args.solver), "solver_sha256": sha(args.solver)}
    trajectories = []
    for verify in (False, True):
        name = "verified" if verify else "unverified"
        params = argparse.Namespace(output=args.output / name, verify=verify, steps=1, sweeps=2,
                                    flow=1, k=.5, epsilon=.5, pseudo_alpha=0, perturbation=0,
                                    ground_projection=False)
        params.output.mkdir()
        prepare(params)
        target = params.output / "target"
        receipt = run([str(args.solver), "case.xml"], target, "solver.log", 240)
        report = load(target / "case.evolution.json")
        require(receipt["exit_code"] != 0 and "Picard sweep limit" in report["failure"],
                "the bounded two-sweep check must reach its deliberate stop")
        require(report["steps"] == report["time"] == 0 and state_hash(report["state"]) == H0,
                "operator check advanced the native trajectory")
        records = stages(target)
        accepted = [row for row in records if row["event"] == "accepted_sweep"]
        require(len(accepted) == 2, "expected two accepted sweeps")
        trajectories.append([(row["state_sha256"], row["norms"]) for row in accepted])
        require(state_hash(report["rejected_state"]) == accepted[-1]["state_sha256"],
                "published trial differs from final logged sweep")
        result[name] = dict(receipt=receipt, states=trajectories[-1])
        if verify:
            solves = {row["block"]: row["dimension"] for row in records
                      if row["event"] == "solve"}
            require(solves == {"velocity_predictor": 980, "pressure_correction": 85,
                               "k": 504, "epsilon_transport_trace": 254}, "incorrect block dimensions")
            check, = [row for row in records if row["event"] == "frozen_saddle_comparison"]
            predicted = check["original_predictor_continuity_gauge"]
            corrected = check["original_corrected_continuity_gauge"]
            require(check["velocity_correction"] > 0 and
                    check["coefficient_difference"] < 1e-10 and
                    max(corrected) < 1e-12 and max(corrected) < max(predicted) and
                    predicted[-1] < 1e-12, "invalid predictor/corrector verification")
            result[name]["saddle"] = check
    require(trajectories[0] == trajectories[1], "verification changed accepted sweep states")
    result["status"] = "passed"
    write_json(args.output / "verification.json", result)


if __name__ == "__main__":
    main()
