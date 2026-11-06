#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check native use of accuracy-refined surrogate actions, not M8 convergence."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import test_cls_epawf_reduced_response as operator_checks
from run_cls_epawf_picard import require


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    args = parser.parse_args()
    subprocess.run([sys.executable, operator_checks.__file__, str(args.run)], check=True)
    rows = [json.loads(line[15:]) for line in
            (args.run / "target/solver.log").read_text().splitlines()
            if line.startswith("EP-AWF Reduced ")]
    used = []
    for precision in rows:
        if (precision["event"] != "source_derivative_precision" or
                not precision["available"] or precision["initial_error"] <= .03):
            continue
        iteration = precision["iteration"]
        later = [row for row in rows if row["iteration"] == iteration and
                 row["seconds"] > precision["seconds"]]
        refinements = [row for row in later if row["event"] == "source_direction_refinement" and
                       row["refinement"] == precision["refinement"] and row["improved"]]
        checks = [row for row in later if row["event"] == "source_refined_direction_check" and
                  row["available"] and row["predicted_merit_slope"] < 0 and
                  row["half_merit_slope"] < 0]
        trials = [row for row in later if row["event"] == "trial" and row["accepted"] and
                  row["proposal_origin"] == "nonlinear_source_refined_direction"]
        points = [row for row in later if row["event"] == "accepted"]
        base = [row for row in rows if row["event"] in ("initial", "accepted") and
                row["seconds"] < precision["seconds"]]
        if not (refinements and checks and trials and points and base):
            continue
        require(max(points[0]["norms"]) < max(base[-1]["norms"]),
                "refined source direction did not reduce the original full residual")
        used.append(dict(precision=precision, refinement=refinements[0],
                         signed_actual_check=checks[0], actual_trial={
                             key: value for key, value in trials[0].items()
                             if key != "qualification_inventory"},
                         original_full_before=max(base[-1]["norms"]),
                         original_full_after=max(points[0]["norms"])))
    require(used, "no accuracy-refined source action produced a checked, accepted residual reduction")
    result = dict(status="passed_source_precision_direction_checks_not_m8_convergence",
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  run=str(args.run), used_precision_refinements=used)
    (args.run / "source-precision-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    print(result["status"])


if __name__ == "__main__":
    main()
