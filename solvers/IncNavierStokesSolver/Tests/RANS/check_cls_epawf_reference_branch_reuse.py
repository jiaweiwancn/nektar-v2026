#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check actual original-start progress across reference-weight switches; not M8 convergence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import test_cls_epawf_reduced_response as operator_checks
from check_cls_epawf_response_chord import records
from run_cls_epawf_picard import H0, require, state_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--control", type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, operator_checks.__file__, str(args.run)], check=True)
    manifest = operator_checks.execution_manifest(args.run)
    stage = manifest["target"]["first_stage_input"]
    require(stage["unknowns"] == 2165 and stage["alpha"] == 400 and stage["time"] == .0025 and
            stage["tolerance"] == 2e-12 and stage["incoming_guess_sha256"] == stage["guess_sha256"] == H0 and
            stage["history_sha256"] == "f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca",
            "reference-branch test did not use the unchanged original h0 target")
    data = records(args.run)
    bridges = [r for r in data if r["event"] == "response_reference_branch_bridge"]
    require(bridges, "no actual accepted endpoint crossed an allowed reference branch")
    used = []
    for trial in [r for r in data if r["event"] == "trial" and r["accepted"]]:
        model = [r for r in data if r["event"] == "direction" and
                 r["iteration"] == trial["iteration"] and r["seconds"] < trial["seconds"]][-1]
        if not model.get("reference_branch_bridged"):
            continue
        previous = [r for r in data if r["event"] in ("initial", "accepted") and
                    r["seconds"] < model["seconds"]][-1]
        require(model["response_model_kind"] == "historical_secant" and
                trial["full_norm"] < max(previous["norms"]) and trial["violation"] < 1e-8,
                "bridged historical model did not produce qualified original-equation progress")
        used.append(dict(direction={k:model[k] for k in [
            "iteration", "basis_dimension", "response_rank", "radius", "directional_refreshes",
            "reference_branch_bridged"]}, trial={k:v for k,v in trial.items() if k != "qualification_inventory"}))
    require(used, "no current-branch checked historical direction produced an actual original-start correction")
    result = dict(status="passed_reference_branch_reuse_checks_not_m8_convergence",
                  run=str(args.run), script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  original_target=stage, branch_switches=bridges, accepted_corrections=used)
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        def normalized_xml(run):
            root = ET.parse(run / "target/case.xml").getroot()
            for parent in root.iter():
                for child in list(parent):
                    if child.tag == "I" and child.attrib.get("PROPERTY") == "RANSWallReducedReferenceBranchReuse":
                        parent.remove(child)
            return ET.tostring(root)
        require(normalized_xml(args.run) == normalized_xml(args.control),
                "reference-branch comparison changes another input or budget")
        images = [operator_checks.execution_manifest(p)["IncNavierStokesSolver"] for p in (args.run, args.control)]
        require(images[0]["sha256"] == images[1]["sha256"] and images[0]["libraries"] == images[1]["libraries"],
                "reference-branch comparison changes executable or libraries")
        control = records(args.control)
        first = used[0]["direction"]["iteration"]
        def prefix(rows):
            return [(r["state_sha256"], state_hash(r["norms"])) for r in rows
                    if r["event"] in ("initial", "accepted") and r["iteration"] <= first]
        require(prefix(data) == prefix(control), "reference reuse changes states before its first direction")
        result["matched_control"] = dict(run=str(args.control), identical_input_except_reference_reuse=True,
            identical_executable_and_libraries=True, original_prefix_bitwise_identical=True,
            accepted_original_residual_histories={str(p): [max(r["norms"]) for r in rows if r["event"] == "accepted"]
                for p, rows in ((args.run, data), (args.control, control))},
            controlled_H_evaluations={str(p): sum(r["event"] == "non_k_response" for r in rows)
                for p, rows in ((args.run, data), (args.control, control))})
    (args.run / "reference-branch-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
