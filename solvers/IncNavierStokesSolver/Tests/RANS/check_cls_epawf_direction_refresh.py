#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check measured direction refresh and constrained re-solve; not M8 convergence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import test_cls_epawf_reduced_response as operator_checks
from check_cls_epawf_response_chord import records
from run_cls_epawf_picard import require, state_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--control", type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, operator_checks.__file__, str(args.run)], check=True)
    data = records(args.run)
    refreshed = [r for r in data if r["event"] == "response_direction_refreshed"]
    require(refreshed, "no measured action refreshed a reused constrained direction")
    packets = [r["path"] for r in data if r["event"] == "model_packet"]
    require(len(set(packets)) == len(packets), "an earlier affine model packet was overwritten")
    useful = []
    for trial in [r for r in data if r["event"] == "trial" and r["accepted"]]:
        same = [r for r in refreshed if r["iteration"] == trial["iteration"] and r["seconds"] < trial["seconds"]]
        if not same:
            continue
        refresh = same[-1]
        direction = [r for r in data if r["event"] == "direction" and r["iteration"] == trial["iteration"] and
                     refresh["seconds"] < r["seconds"] < trial["seconds"]][-1]
        require(direction["response_model_kind"] == "historical_secant" and direction["qp_solved"] and
                direction["directional_refreshes"] == refresh["refresh"],
                "the measured action was not used to re-solve the constrained model")
        # Every re-solved model is retained separately: a successful small
        # interpolation defect must not overwrite the original failed model.
        packet = [r for r in data if r["event"] == "model_packet" and r["iteration"] == trial["iteration"] and
                  refresh["seconds"] < r["seconds"] < direction["seconds"]][-1]
        require("-refresh"+str(refresh["refresh"])+"-" in packet["path"] and
                (args.run / "target" / packet["path"]).is_file(), "refreshed affine evidence was overwritten")
        checked = [r for r in data if r["event"] == "response_chord_direction_check" and r["passed"] and
                   r["iteration"] == trial["iteration"] and direction["seconds"] < r["seconds"] < trial["seconds"]]
        require(checked and all(r["bound_scaled_observation_model_error"] <= .1 for r in checked),
                "a refreshed direction lacks an independent bound-scaled actual response check")
        useful.append(dict(refresh=refresh, direction={k:direction[k] for k in [
            "iteration", "basis_dimension", "response_rank", "active_constraints", "radius", "displacement",
            "predicted_norm", "predicted_violation", "directional_refreshes"]},
            checks=checked, trial={k:v for k,v in trial.items() if k != "qualification_inventory"}))
    require(useful, "refreshed actions did not produce an actually accepted original-equation correction")
    result = dict(status="passed_direction_refresh_checks_not_m8_convergence", run=str(args.run),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), corrections=useful)
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        def normalized_xml(run):
            root = ET.parse(run / "target/case.xml").getroot()
            for parent in root.iter():
                for child in list(parent):
                    if child.tag == "P" and (child.text or "").split("=")[0].strip() == "RANSWallReducedDirectionalRefreshUpdates":
                        parent.remove(child)
            return ET.tostring(root)
        require(normalized_xml(args.run) == normalized_xml(args.control),
                "direction refresh control changes another physical or numerical input")
        images = [operator_checks.execution_manifest(run)["IncNavierStokesSolver"] for run in (args.run, args.control)]
        require(images[0]["sha256"] == images[1]["sha256"] and images[0]["libraries"] == images[1]["libraries"],
                "direction refresh comparison changes executable or libraries")
        control = records(args.control)
        first = refreshed[0]["iteration"]
        def prefix(rows):
            return [(r["state_sha256"], state_hash(r["norms"])) for r in rows
                    if r["event"] in ("initial", "accepted") and r["iteration"] <= first]
        require(prefix(data) == prefix(control), "direction refresh changes states before it is used")
        result["matched_control"] = dict(run=str(args.control), identical_input_except_refresh_budget=True,
            identical_executable_and_libraries=True, original_prefix_bitwise_identical=True,
            accepted_original_residual_histories={str(run): [max(r["norms"]) for r in rows if r["event"] == "accepted"]
                for run, rows in ((args.run, data), (args.control, control))},
            controlled_H_evaluations={str(run): sum(r["event"] == "non_k_response" for r in rows)
                for run, rows in ((args.run, data), (args.control, control))})
    (args.run / "direction-refresh-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
