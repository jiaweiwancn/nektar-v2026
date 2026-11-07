#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check actual verified reduced trials after optional source work; not M8 convergence."""
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
    here = Path(__file__).resolve().parent
    subprocess.run([sys.executable, str(here / "check_cls_epawf_reference_branch_reuse.py"),
                    str(args.run)], check=True)
    rows = records(args.run)
    deferred = [r for r in rows if r["event"] == "source_refinement_deferred"]
    useful = []
    for trial in [r for r in rows if r["event"] == "trial" and r["accepted"]]:
        model = [r for r in rows if r["event"] == "direction" and
                 r["iteration"] == trial["iteration"] and r["seconds"] < trial["seconds"]][-1]
        events = [r for r in deferred if r["phase"] == "verified_linear_trial" and
                  r["iteration"] == trial["iteration"] and model["seconds"] < r["seconds"] < trial["seconds"]]
        if not events:
            continue
        require(trial["proposal_origin"] == "linear_reduced_direction" and
                model["reference_branch_bridged"] and model["directional_refreshes"] >= 2,
                "actual fallback trial did not use the re-solved measured reduced direction")
        require(not any(r["event"] == "source_direction_refinement" and r["iteration"] == trial["iteration"] and
                        events[-1]["seconds"] < r["seconds"] < trial["seconds"] for r in rows),
                "optional source refinements blocked the declared actual linear trial")
        require(any(r["event"] == "source_ray_selection" and r["iteration"] == trial["iteration"] and
                    model["seconds"] < r["seconds"] < trial["seconds"] for r in rows),
                "verified linear trial discarded the complete nonlinear destruction ray model")
        useful.append(dict(model={k: model[k] for k in ["iteration", "basis_dimension", "response_rank",
                      "directional_refreshes", "radius", "displacement"]}, deferred=events[-1],
                      trial={k: v for k, v in trial.items() if k != "qualification_inventory"}))
    require(useful, "no qualified actual original-start correction exercised verified-linear fallback")
    result = dict(status="passed_verified_linear_trials_not_m8_convergence", run=str(args.run),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope="Actual constrained reduced correction after measured source-response refresh; complete original equations and qualification remain authoritative",
                  corrections=useful)
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        def xml_without_policy(path):
            root = ET.parse(path / "target/case.xml").getroot()
            for parent in root.iter():
                for child in list(parent):
                    if child.tag == "I" and child.attrib.get("PROPERTY") == "RANSWallReducedPreferVerifiedLinear":
                        parent.remove(child)
            return ET.tostring(root)
        require(xml_without_policy(args.run) == xml_without_policy(args.control),
                "verified-linear comparison changes another input")
        images = [operator_checks.execution_manifest(p)["IncNavierStokesSolver"] for p in (args.run, args.control)]
        require(images[0]["sha256"] == images[1]["sha256"] and images[0]["libraries"] == images[1]["libraries"],
                "verified-linear comparison changes executable or libraries")
        first = useful[0]["model"]["iteration"]
        control = records(args.control)
        def prefix(data):
            return [(r["state_sha256"], state_hash(r["norms"])) for r in data
                    if r["event"] in ("initial", "accepted") and r["iteration"] <= first]
        require(prefix(rows) == prefix(control), "preference changed the original-start prefix before its first use")
        result["matched_control"] = dict(run=str(args.control), input_only_policy_differs=True,
            same_executable_and_libraries=True, original_prefix_bitwise_identical=True,
            residual_histories={str(path): [max(r["norms"]) for r in data if r["event"] == "accepted"]
                                for path, data in ((args.run, rows), (args.control, control))},
            H_evaluations={str(path): sum(r["event"] == "non_k_response" for r in data)
                           for path, data in ((args.run, rows), (args.control, control))})
    (args.run / "verified-linear-trial-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
