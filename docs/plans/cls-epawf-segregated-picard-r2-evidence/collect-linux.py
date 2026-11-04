#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Summarize the checked Picard experiments while leaving heavy payloads on Linux."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

ATTEMPTS = ("conditional-original", "conditional-endpoint",
            "reduced-original", "reduced-endpoint", "reduced-radius")


def load(path, **kwargs):
    return json.loads(path.read_text(), **kwargs)


def pin(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest()}


def diagnostic(path):
    data = load(path)
    keys = ("alpha", "time", "guess_sha256", "history_sha256", "outer_range",
            "outer_gradient_witness", "k_norm", "ground_extrema")
    summary = {name: data[name] for name in keys}
    summary["artifact"] = pin(path)
    summary["dimension"] = len(data["free_rows"])
    summary["actions"] = [{k: v for k, v in action.items() if k != "direction"}
                          for action in data["actions"]]
    assert summary["dimension"] == 504
    assert len(summary["actions"]) == 6
    assert all(action["samples"] and action["best_relative_error"] < 3e-6
               for action in summary["actions"]), "restricted action check exceeded tolerance"
    for name in ("frozen_projection", "true_projection"):
        p = data[name]
        row = {"max_unprojected": max(map(abs, p["unprojected"])),
               "max_projected": max(map(abs, p["projected"]), default=None),
               "min_affine_ground": min(x["value"] for x in p["affine_ground_extrema"]),
               "cuts": p["cuts"], "trials": p["trials"], "kkt": p["kkt"]}
        if row["kkt"]:
            row["kkt"]["active"] = [{k: v for k, v in active.items() if k != "row"}
                                    for active in row["kkt"]["active"]]
        summary[name] = row
    return summary


def history_identity(path):
    report = load(path, parse_int=float)

    def coefficients(values):
        return struct.pack(f"<{len(values)}d", *values)

    return (report["steps"], report["time"], report["native_history_times"],
            coefficients(report["state"]), [coefficients(v) for v in report["native_history"]])


def attempt(root, name):
    case = root / name
    manifest = load(case / "manifest.json")
    summary = {key: manifest[key] for key in
               ("source", "status", "failure", "thread_limits", "target", "readback", "xml_sha256")
               if key in manifest}
    summary["manifest"] = pin(case / "manifest.json")
    startup = history_identity(case / "startup/case.evolution.json")
    summary["forced_rejection_preserved_history"] = (
        startup == history_identity(case / "rejected/case.evolution.json"))
    summary["failed_stage_preserved_history"] = (
        startup == history_identity(case / "target/case.evolution.json"))
    assert summary["forced_rejection_preserved_history"]
    if manifest["target"]["steps"] == 0:
        assert summary["failed_stage_preserved_history"]
    sweeps = load(case / "target/sweeps.json")
    summary["sweeps"] = [{key: r[key] for key in ("sweep", "norms", "relaxation", "seconds",
                                                 "trace_coefficient_defect", "solve_counts", "state_sha256")} |
                          {"ground_minimum": r["ground"]["minimum"],
                           "finite_scalar_margin": r["scalar_domain"]["minimum_dimensionless_margin"]}
                          for r in sweeps if r["event"] in ("initial", "accepted_sweep")]
    summary["qualification_samples"] = [r for r in sweeps if r["event"] == "qualification"]
    summary["reduced_trials"] = [r for r in sweeps if r["event"] == "reduced_k_trial"]
    summary["conditional_k"] = {"trial_count": sum(r["event"] == "conditional_k_trial" for r in sweeps),
                                "accepted_locally": sum(r["event"] == "conditional_k_trial" and r["accepted"]
                                                       for r in sweeps),
                                "fallback_count": sum(r["event"] == "conditional_k_fallback" for r in sweeps)}
    if "readback" in summary and manifest["target"]["steps"] == 0:
        assert summary["readback"]["original_block_norms"] == summary["sweeps"][-1]["norms"]
    audit = root / "terminal-audit" / name / "audit.json"
    if audit.is_file():
        entry = load(audit)[0]["audit"]
        summary["terminal_audit"] = pin(audit) | {k: entry[k] for k in
            ("r2_qualified", "r2_outer_range", "r2_covariance", "r4_outer_gradient_witness")}
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root
    result = {"verification": load(root / "verification/verification.json"),
              "final_verification": load(root / "final-verification/verification.json"),
              "binaries": (root / "binaries.sha256").read_text().splitlines(),
              "regressions": pin(root / "regression.log"),
              "constitutive": pin(root / "robust-kernel.json"),
              "legacy_matrix_hashes": (root / "legacy-action-regression/matrices.sha256").read_text().splitlines(),
              "diagnostics": {name: diagnostic(root / name / "target/k-diagnostic.json")
                              for name in ATTEMPTS[:2]},
              "attempts": {name: attempt(root, name) for name in ATTEMPTS}}
    feedback = load(root / "feedback/feedback.json")
    for case in feedback["cases"]:
        case.get("final", {}).pop("k_residual", None)
    result["feedback"] = feedback
    result["qualification_direction"] = load(root / "qualification-direction/comparison.json")
    for name in ("conditional-terminal", "best-terminal"):
        terminal = root / f"{name}-diagnostic/target/k-diagnostic.json"
        result["diagnostics"][name] = diagnostic(terminal)
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
