#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check retained native current-base actions and useful trials; not M8 convergence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from check_cls_epawf_response_chord import records
from run_cls_epawf_picard import require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--control", type=Path, required=True)
    args = parser.parse_args()
    directory = Path(__file__).resolve().parent
    subprocess.run([sys.executable, str(directory / "check_cls_epawf_reference_branch_reuse.py"),
                    str(args.run), "--control", str(args.control)], check=True)
    rows = records(args.run)
    measured = [r for r in rows if r["event"] == "response_direction_refreshed"]
    protected = [r for r in measured if r["preserved_directions"] > 0]
    require(protected, "no actual correlated current-base measurements exercised protection")
    useful = []
    for trial in [r for r in rows if r["event"] == "trial" and r["accepted"]]:
        direction = [r for r in rows if r["event"] == "direction" and
                     r["iteration"] == trial["iteration"] and r["seconds"] < trial["seconds"]][-1]
        if not direction.get("reference_branch_bridged") or direction["directional_refreshes"] < 2:
            continue
        checks = [r for r in rows if r["event"] == "response_chord_direction_check" and
                  r["iteration"] == trial["iteration"] and direction["seconds"] < r["seconds"] < trial["seconds"]]
        require(checks and all(r["passed"] for r in checks),
                "protected direction lacked accurate independent current-base action checks")
        useful.append(dict(direction={k: direction[k] for k in [
            "iteration", "directional_refreshes", "basis_dimension", "response_rank", "radius", "displacement"]},
            trial={k: v for k, v in trial.items() if k != "qualification_inventory"}))
    require(useful, "protected measurements were never used in an accepted qualified original-start correction")
    result = dict(status="passed_protected_current_actions_not_m8_convergence", run=str(args.run),
                  scope="Preserved current-base measured responses and actual useful native trial, unchanged original h0 target; no nonlinear stationarity or physical acceptance claim",
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  refreshes=measured, accepted_corrections=useful)
    (args.run / "protected-action-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
