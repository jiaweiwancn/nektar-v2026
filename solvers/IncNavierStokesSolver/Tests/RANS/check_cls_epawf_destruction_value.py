#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check restricted destruction values and matched native trajectories, not M8 convergence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

import test_cls_epawf_reduced_response as operator_checks
from run_cls_epawf_picard import require, state_hash


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
    checks = [r for r in rows if r["event"] == "restricted_destruction_value_check"]
    require(len(checks) >= 3 and all(r["bitwise_equal"] and r["rows"] == 504 for r in checks),
            "fixture did not check restricted destruction values on three native samples")
    result = dict(status="passed_restricted_destruction_value_checks_not_m8_convergence",
                  run=str(args.run), script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  checked_values=checks)
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        control = records(args.control)
        def normalized_xml(run):
            root = ET.parse(run / "target/case.xml").getroot()
            for parent in root.iter():
                for child in list(parent):
                    if child.get("PROPERTY") == "RANSWallReducedFastDestruction":parent.remove(child)
            return ET.tostring(root)
        require(normalized_xml(args.run) == normalized_xml(args.control),
                "restricted destruction comparison changes another numerical or physical input")
        manifests = [operator_checks.execution_manifest(run) for run in (args.run, args.control)]
        images = [m["IncNavierStokesSolver"] for m in manifests]
        require(images[0]["sha256"] == images[1]["sha256"] and images[0]["libraries"] == images[1]["libraries"],
                "restricted destruction comparison changes executable or libraries")
        protocols = [[r for r in data if r["event"] == "protocol"][0] for data in (rows, control)]
        require(protocols[0]["restricted_destruction_value"] and
                not protocols[1]["restricted_destruction_value"], "incorrect destruction-value selectors")
        def points(data):
            return [(r["state_sha256"], state_hash(r["norms"])) for r in data
                    if r["event"] in ("initial", "accepted", "converged")]
        require(points(rows) == points(control), "restricted destruction changes accepted native states or original norms")
        result["matched_control"] = dict(run=str(args.control), identical_input_except_component_selector=True,
            identical_executable_and_libraries=True, all_original_states_and_norms_bitwise_identical=True,
            target_seconds={str(run): m["target"]["seconds"] for run, m in zip((args.run,args.control),manifests)},
            accepted_original_residual_history=[max(r["norms"]) for r in rows if r["event"] == "accepted"])
    (args.run / "destruction-value-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
