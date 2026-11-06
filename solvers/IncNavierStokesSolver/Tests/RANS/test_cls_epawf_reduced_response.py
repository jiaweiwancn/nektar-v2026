#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Check controlled reduced evaluations and the constrained direction diagnostics."""

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

from run_cls_epawf_picard import require, load, state_hash, H0, write_json, sha


def execution_manifest(directory):
    """Use the actual native launch receipt for ordinary or diagnostic runs."""
    ordinary = directory / "manifest.json"
    if ordinary.exists():
        return load(ordinary)
    continuation = load(directory.parent / "continuation.json")
    matching = [case for case in continuation["cases"]
                if Path(case["path"]).resolve() == directory.resolve()]
    require(len(matching) == 1 and continuation["method"] == "reduced",
            "missing or ambiguous reduced continuation launch receipt")
    case, = matching
    require(case.get("history_preserved") and case["accepted_steps"] == 0 and
            case["accepted_time"] == 0 and case["native"]["command"] ==
            [continuation["solver"], "case.xml"],
            "continuation operator check lacks a completed no-advance receipt")
    return dict(status=case["status"], scope=continuation["scope"],
                IncNavierStokesSolver=dict(path=continuation["solver"],
                    sha256=continuation["solver_sha256"],
                    libraries=continuation["IncNavierStokesSolver"]["libraries"]))


def check(directory):
    manifest = execution_manifest(directory)
    report = load(directory / "target/case.evolution.json")
    records = []
    for line in (directory / "target/solver.log").read_text().splitlines():
        if line.startswith("EP-AWF Reduced "):
            records.append(json.loads(line[len("EP-AWF Reduced "):]))
    repeat, = [r for r in records if r["event"] == "repeatability"]
    require(repeat["failed_response_rolled_back"] and repeat["k_coefficients_fixed"],
            "controlled response rollback/fixed-k assertion failed")
    require(repeat["same_k_difference"] < 1e-10 and repeat["tightened_k_difference"] < 2e-6,
            "repeated reduced residual is inconsistent with the inner tolerance")
    checks = [r for r in records if r["event"] == "response_check"]
    require(len(checks) >= 3 and sum(c["usable"] for c in checks) >= 3,
            "insufficient independently consistent reduced responses")
    warm_starts = [r for r in records if r["event"] == "half_response_warm_start"]
    if any(r.get("half_response_warm_start") for r in records if r["event"] == "protocol"):
        require(warm_starts, "the half-response warm start was not checked independently")
    require(all(r["attained"] <= r["required"] and
                r["canonical_attained"] <= r["required"] and
                r["same_local_response"] == (r["relative_action_difference"] <= .01 and
                                             r["observation_action_difference"] <= .01)
                for r in warm_starts), "invalid half-response warm-start status")
    directions = [r for r in records if r["event"] == "direction"]
    require(directions, "the reduced direction was not constructed")
    protocol, = [r for r in records if r["event"] == "protocol"]
    require(all(d["basis_dimension"] <= protocol["maximum_basis"] <= 512 for d in directions),
            "the direction exceeds its declared adaptive response budget")
    require(0 <= protocol["response_krylov_budget"] <= 128,
            "the response Krylov budget exceeds its explicit limit")
    require(any(d["qp_solved"] and d["displacement"] > 0 for d in directions),
            "no nonzero constrained model direction was solved")
    require(all(d["displacement"] <= d["radius"]*(1+1e-12) for d in directions),
            "final correction exceeds its radius")
    require(all(d["orthogonal_probes"] <= d["basis_dimension"] and
                (d["orthogonal_probes"] != d["basis_dimension"] or
                 d["maximum_column_correlation"] < 1e-10)
                for d in directions if "orthogonal_probes" in d),
            "orthogonal field probes retained correlated model columns")
    require(all(d["radius_model_solves"] >= 1 and d["radius_multiplier"] >= 0. and
                d["radius_complementarity"] >= 0. for d in directions
                if "radius_model_solves" in d and d["qp_solved"]),
            "invalid radius-constrained model status")
    require(all(d["augmented_face_solves"] >= 0 and d["augmented_face_fallbacks"] >= 0
                for d in directions if "augmented_face_solves" in d),
            "invalid response-coordinate face-solve accounting")
    require(all(0 <= d["maximum_ball_redundant_cuts"] <= d["original_model_cuts"]
                for d in directions if "maximum_ball_redundant_cuts" in d),
            "invalid conservative ball-cut filter accounting")
    require(all(0 <= d["warm_face_fallbacks"] <= d["warm_face_solves"] <= d["radius_model_solves"]
                for d in directions if "warm_face_solves" in d),
            "invalid large-model warm-face solve accounting")
    for direction in directions:
        for name in direction["basis_names"]:
            matching = [c for c in checks if c["name"] == name and
                        c["iteration"] == direction["iteration"] and
                        c["seconds"] <= direction["seconds"]]
            require(matching and matching[-1]["usable"] and
                    matching[-1]["relative_error"] <= .03 and
                    matching[-1]["observation_relative_error"] <= .05,
                    "a direction used an unchecked or inaccurate response")
            checked_warm_starts = [r for r in warm_starts if r["name"] == name and
                                  r["iteration"] == direction["iteration"] and
                                  r["seconds"] <= direction["seconds"]]
            require(not checked_warm_starts or checked_warm_starts[-1]["same_local_response"],
                    "a direction used an inconsistent warm-start branch")
    trials = [r for r in records if r["event"] == "trial"]
    require(trials, "no constrained direction received an equilibrated native trial")
    corrections = [r for r in records if r["event"] == "second_order_trial"]
    require(all(r["correction"] <= r["correction_radius"]*(1+1e-12) and
                r["displacement"] <= r["radius"]*(1+1e-12) and
                (not r["valid"] or r["non_k_error"] <= r["required"])
                for r in corrections),
            "nonlinear qualification correction violated a field bound or H accuracy")
    require(all(r["violation"] < r["before_violation"] for r in corrections if r["improved"]),
            "a failed nonlinear qualification correction was marked useful")
    secants = [r for r in records if r["event"] == "constraint_response_secant"]
    require(all(r["field_signal"] > r["minimum_field_signal"] and not r["failure"]
                for r in secants if r["updated"]),
            "a below-signal or failed constraint response was updated")
    combined = [r for r in records if r["event"] == "combined_direction_check"]
    require(all(r["attained"] <= r["required"] and r["half_attained"] <= r["required"] and
                r["half_relative_residual_error"] <= .03 and
                r["half_observation_relative_error"] <= .05
                for r in combined if r["available"]),
            "combined response does not match its checked reduced model")
    inward = [r for r in records if r["event"] == "ground_inward_roundoff"]
    require(all(r["after"]["admissible"] and r["displacement"] <= r["radius"]*(1+1e-12)
                for r in inward), "ground-inward correction violated the ground or radius bound")
    responses = [r for r in records if r["event"] == "non_k_response"]
    require(all(r["attained"] <= r["required"] for r in responses if r["converged"]),
            "a partially converged response was accepted")
    prepared = [r for r in records if r["event"] == "response_preconditioner_prepared"]
    if protocol.get("response_field_preconditioner_reuse"):
        cached = [r for r in responses if r.get("base_field_preconditioner")]
        require(prepared and cached and warm_starts,
                "response preconditioner reuse lacks native applications or independent fresh branch checks")
        require(all(r["velocity_size"] > 0 and r["pressure_gauge_size"] > 0 and
                    r["epsilon_trace_size"] > 0 for r in prepared),
                "response preconditioner has an invalid separate-field shape")
    require(all(0 <= r.get("cached_flow_solves", 0) <= r["flow_solves"] and
                0 <= r.get("cached_epsilon_solves", 0) <= r["epsilon_solves"] and
                (not r.get("base_field_preconditioner") or
                 r["purpose"].startswith(("response_", "half_response_")))
                for r in responses),
            "cached response accounting is invalid or actual trial did not use current operators")
    source_selections = [r for r in records if r["event"] == "source_ray_selection"]
    require(all(0 <= r["fraction"] <= 1 and r["evaluations"] <= 12 and
                r["invalid"] <= r["evaluations"] and
                (not r["fraction"] or r["predicted_norm"] < r["base_norm"]*(1-1e-6))
                for r in source_selections), "invalid surrogate-only source-ray selection")
    source_predictions = [r for r in records if r["event"] == "source_ray_prediction"]
    require(all(r["non_k_error"] <= r["required"] and
                r["actual_displacement"] <= r["radius"]*(1+1e-12)
                for r in source_predictions),
            "a source-ray actual trial lacked controlled H or exceeded its radius")
    refinements = [r for r in records if r["event"] == "source_direction_refinement"]
    anchors = [r for r in records if r["event"] == "source_linearization_anchor"]
    require(protocol.get("interior_source_model_anchors", False) or not anchors,
            "an interior source anchor was used without its selector")
    require(all(1 <= r["attempts"] <= 3 and 0 <= r["anchor_reductions"] <= 2 and
                r["anchor_fraction"] == 2.**(-r["anchor_reductions"]) and
                0 <= r["model_evaluations"] <= r["model_evaluation_budget"] ==
                4*r["dimension"]+64 and
                (not r["available"] or r["available_derivatives"] == r["dimension"] and
                 r["two_scale_error"] <= .03) for r in anchors),
            "source anchor exceeded its retry/evaluation budget or used unresolved derivatives")
    require(all(r["dimension"] == r["available_derivatives"] and r["two_scale_error"] <= .03 and
                r["predicted_norm"] < r["before_norm"] and r["original_merit_slope"] < 0.
                for r in refinements if r["improved"]),
            "a source-model direction used inaccurate derivatives or lacks predicted descent")
    refined_directions = [r for r in records if r["event"] == "source_refined_direction"]
    require(all(r["displacement"] <= r["radius"]*(1+1e-12) and r["original_merit_slope"] < 0.
                for r in refined_directions), "source refinement violates its bound or true linear descent")
    refined_checks = [r for r in records if r["event"] == "source_refined_direction_check"]
    require(all(r["attained"] <= r["required"] and r["half_attained"] <= r["required"] and
                r["half_relative_residual_error"] <= .03
                for r in refined_checks if r["available"]),
            "source-refined direction does not match the independent local reduced action")
    if report["steps"] == 0:
        startup = load(directory / "startup/case.evolution.json")
        require(state_hash(report["state"]) == H0 and report["time"] == 0 and
                [state_hash(slot) for slot in report["native_history"]] ==
                [state_hash(slot) for slot in startup["native_history"]] and
                report["native_history_times"] == startup["native_history_times"],
                "reduced failure advanced the native state or history")
    result = dict(status="passed_operator_checks", manifest_status=manifest["status"],
                  repeatability=repeat, half_response_warm_starts=warm_starts,
                  action_checks=checks, directions=directions,
                  trials=trials, combined_direction_checks=combined,
                  second_order_corrections=corrections, constraint_response_secants=secants,
                  source_ray_selections=source_selections, source_ray_predictions=source_predictions,
                  source_direction_refinements=refinements, source_refined_directions=refined_directions,
                  source_refined_direction_checks=refined_checks)
    result["response_preconditioners"] = prepared
    result["source_linearization_anchors"] = anchors
    write_json(directory / "reduced-checks.json", result)


def compare_verification(verified, unverified):
    """Compare identical native numerical inputs with observation On/Off."""
    def numerical_xml(directory):
        root = ET.parse(directory / "target/case.xml").getroot()
        values = [e.get("VALUE") for e in root.iter()
                  if e.get("PROPERTY") == "RANSWallReducedVerify"]
        require(values == ["On" if directory == verified else "Off"],
                "the comparison requires verified and unverified native inputs")
        for element in root.iter():
            if element.get("PROPERTY") == "RANSWallPicardGuessFile":
                guess = Path(element.get("VALUE"))
                if not guess.is_absolute():
                    guess = directory / "target" / guess
                require(guess.is_file(), "verification guess file is missing")
                # Continuation writes separate copies of the same packet.
                # Compare every input byte, rather than treating those paths
                # as different numerical data or ignoring the guess input.
                element.set("VALUE", "sha256:" + sha(guess))
        for element in root.iter():
            for child in list(element):
                if child.get("PROPERTY") == "RANSWallReducedVerify":
                    element.remove(child)
        return ET.tostring(root)

    require(numerical_xml(verified) == numerical_xml(unverified),
            "verification comparison changes a numerical input or iteration budget")
    a, b = (execution_manifest(p) for p in (verified, unverified))
    require(a["IncNavierStokesSolver"]["sha256"] ==
            b["IncNavierStokesSolver"]["sha256"] and
            a["IncNavierStokesSolver"]["libraries"] ==
            b["IncNavierStokesSolver"]["libraries"],
            "verification comparison uses different solver/runtime identities")
    x, y = (load(p / "target/case.evolution.json") for p in (verified, unverified))
    for key in ("steps", "time", "native_history_times", "last_stage_alpha", "last_stage_time"):
        require(x[key] == y[key], "verification changes native time or stage bookkeeping")
    hashes = {}
    for key in ("state", "rejected_state", "initial", "last_stage_history"):
        hashes[key] = state_hash(x[key])
        require(hashes[key] == state_hash(y[key]),
                "verification changes native coefficients or physical history")
    require([state_hash(q) for q in x["native_history"]] ==
            [state_hash(q) for q in y["native_history"]],
            "verification changes a physical history slot")
    require([state_hash(q) for q in x["accepted_states"]] ==
            [state_hash(q) for q in y["accepted_states"]],
            "verification changes an accepted native state")
    def numerical_updates(directory):
        result = []
        for line in (directory / "target/solver.log").read_text().splitlines():
            if not line.startswith("EP-AWF Reduced "):
                continue
            row = json.loads(line[15:])
            if row["event"] in ("initial", "accepted", "converged", "converged_algebraic_diagnostic"):
                result.append({k: row[k] for k in
                               ("event", "iteration", "alpha", "time", "state_sha256", "norms")})
        require(result, "verification comparison lacks numerical update receipts")
        return result

    updates = numerical_updates(verified)
    require(updates == numerical_updates(unverified),
            "verification changes an intermediate numerical update or its original residual")
    write_json(verified / "verification-comparison.json",
               dict(status="passed_reduced_verification_observability",
                    scope="Bitwise numerical trajectory and history comparison; not M8 convergence",
                    verified=str(verified), unverified=str(unverified),
                    solver_sha256=a["IncNavierStokesSolver"]["sha256"],
                    numerical_updates=updates,
                    state_sha256=hashes, accepted_physical_steps=x["steps"], time=x["time"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--compare-unverified", type=Path,
                        help="Also require bitwise identical numerical states/history with verification Off")
    args = parser.parse_args()
    check(args.directory)
    if args.compare_unverified:
        compare_verification(args.directory, args.compare_unverified)
