#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check native, freshly safeguarded response reuse; not M8 convergence."""
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
    directions = [r for r in rows if r["event"] == "direction"]
    used = []
    for trial in [r for r in rows if r["event"] == "trial" and r["accepted"]]:
        direction = [r for r in directions if r["iteration"] == trial["iteration"] and
                     r["seconds"] < trial["seconds"]][-1]
        if direction.get("response_model_kind") == "historical_secant":
            used.append(dict(direction={k:direction[k] for k in [
                "iteration", "basis_dimension", "response_model_age", "response_rank", "radius"]},
                actual_trial={k:v for k,v in trial.items() if k != "qualification_inventory"}))
    require(used, "no historical response produced a freshly checked accepted original-equation correction")
    # A successful interpolation/model KKT solve must never be promoted to
    # a measured current reduced first-order certificate.
    output = args.run / "rejected-chord-first-order-audit.json"
    audit = subprocess.run([sys.executable, str(Path(__file__).with_name(
        "audit_cls_epawf_reduced_first_order.py")), str(args.run.parent), str(output),
        "--case", args.run.name+":"+str(used[0]["direction"]["iteration"])],
        capture_output=True, text=True, timeout=600)
    require(audit.returncode != 0 and "not a freshly measured original first-order audit" in audit.stderr and
            not output.exists(), "historical response was mislabelled as fresh first-order information")
    result = dict(status="passed_response_chord_checks_not_m8_convergence", run=str(args.run),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  accepted_chord_corrections=used, fresh_first_order_audit_refused=True)
    if args.control:
        subprocess.run([sys.executable, operator_checks.__file__, str(args.control)], check=True)
        def normalized_xml(run):
            root = ET.parse(run / "target/case.xml").getroot()
            for parent in root.iter():
                for child in list(parent):
                    if child.tag == "P" and (child.text or "").split("=")[0].strip() == "RANSWallReducedResponseReuseUpdates":
                        parent.remove(child)
            return ET.tostring(root)
        require(normalized_xml(args.run) == normalized_xml(args.control),
                "response chord comparison changes another numerical or physical input")
        manifests = [operator_checks.execution_manifest(run) for run in (args.run, args.control)]
        images = [m["IncNavierStokesSolver"] for m in manifests]
        require(images[0]["sha256"] == images[1]["sha256"] and images[0]["libraries"] == images[1]["libraries"],
                "response chord comparison changes executable or libraries")
        control = records(args.control)
        protocols = [[r for r in data if r["event"] == "protocol"][0] for data in (rows, control)]
        require(protocols[0]["response_chord_update_budget"] > 0 and
                protocols[1]["response_chord_update_budget"] == 0, "invalid response reuse control")
        first = min(r["iteration"] for r in rows if r["event"] == "response_chord_reused")
        def prefix(data):
            return [(r["state_sha256"], state_hash(r["norms"])) for r in data
                    if r["event"] in ("initial", "accepted") and r["iteration"] <= first]
        require(prefix(rows) == prefix(control), "response reuse changes states before the first reused model")
        result["matched_control"] = dict(run=str(args.control), identical_input_except_reuse_budget=True,
            identical_executable_and_libraries=True, original_prefix_bitwise_identical=True,
            accepted_original_residual_histories={str(run): [max(r["norms"]) for r in data if r["event"] == "accepted"]
                for run, data in ((args.run, rows), (args.control, control))},
            controlled_H_evaluations={str(run): sum(r["event"] == "non_k_response" for r in data)
                for run, data in ((args.run, rows), (args.control, control))})
    (args.run / "response-chord-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
