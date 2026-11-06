#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check a bounded native fixture using recovered interior source-model derivatives.

These are direction, original-operator and history checks, not M8 convergence.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import test_cls_epawf_reduced_response as operator_checks
from run_cls_epawf_picard import load, require, state_hash


def records(run):
    return [json.loads(line[15:], parse_int=float)
            for line in (run / "target/solver.log").read_text().splitlines()
            if line.startswith("EP-AWF Reduced ")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--control", type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, operator_checks.__file__, str(args.run)], check=True)
    rows = records(args.run)
    protocol, = [row for row in rows if row["event"] == "protocol"]
    require(protocol.get("interior_source_model_anchors"), "interior source anchor is not enabled")
    anchors = [row for row in rows if row["event"] == "source_linearization_anchor" and
               row["available"] and row["anchor_reductions"] > 0]
    require(anchors, "fixture did not recover blocked derivatives at an interior anchor")
    used = []
    for anchor in anchors:
        iteration = anchor["iteration"]
        refinements = [row for row in rows if row["event"] == "source_direction_refinement" and
                       row["iteration"] == iteration and row["refinement"] == anchor["refinement"] and
                       row["seconds"] > anchor["seconds"] and row["improved"]]
        checks = [row for row in rows if row["event"] == "source_refined_direction_check" and
                  row["iteration"] == iteration and row["available"] and
                  row["predicted_merit_slope"] < 0 and row["half_merit_slope"] < 0]
        trials = [row for row in rows if row["event"] == "trial" and row["iteration"] == iteration and
                  row["accepted"] and row["proposal_origin"] == "nonlinear_source_refined_direction"]
        if refinements and checks and trials:
            used.append(dict(anchor=anchor, refinement=refinements[0],
                             signed_actual_check=checks[0], actual_trial={
                                 key: value for key, value in trials[0].items()
                                 if key != "qualification_inventory"}))
    require(used, "recovered model derivatives did not produce an independently checked accepted direction")
    result = dict(status="passed_source_anchor_direction_checks_not_m8_convergence",
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  run=str(args.run), used_anchors=used)
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        def numerical_xml(run, enabled):
            root = ET.parse(run / "target/case.xml").getroot()
            values = [row.get("VALUE") for row in root.iter()
                      if row.get("PROPERTY") == "RANSWallReducedInteriorSourceAnchor"]
            require(values == ["On" if enabled else "Off"], "incorrect source-anchor control selector")
            for parent in root.iter():
                for child in list(parent):
                    if child.get("PROPERTY") == "RANSWallReducedInteriorSourceAnchor":
                        parent.remove(child)
            return ET.tostring(root)
        require(numerical_xml(args.run, True) == numerical_xml(args.control, False),
                "source-anchor comparison changes another numerical or physical input")
        manifests = [operator_checks.execution_manifest(run) for run in (args.run, args.control)]
        images = [manifest["IncNavierStokesSolver"] for manifest in manifests]
        require(images[0]["sha256"] == images[1]["sha256"] and
                images[0]["libraries"] == images[1]["libraries"],
                "source-anchor comparison uses different executable/runtime identities")
        control = records(args.control)
        first_anchor = min(row["iteration"] for row in anchors)
        points = [[row for row in data if row["event"] in ("initial", "accepted") and
                   row["iteration"] <= first_anchor] for data in (rows, control)]
        require([(row["state_sha256"], state_hash(row["norms"])) for row in points[0]] ==
                [(row["state_sha256"], state_hash(row["norms"])) for row in points[1]],
                "the trajectory differs before the first recovered anchor")
        result["matched_control"] = dict(run=str(args.control), identical_input_except_selector=True,
            identical_executable_and_libraries=True, prior_updates_bitwise_identical=True,
            residual_histories={str(run): [max(row["norms"]) for row in data
                                          if row["event"] == "accepted"]
                               for run, data in ((args.run, rows), (args.control, control))},
            accepted_physical_steps=[load(run / "target/case.evolution.json")["steps"]
                                     for run in (args.run, args.control)])
    (args.run / "source-anchor-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
