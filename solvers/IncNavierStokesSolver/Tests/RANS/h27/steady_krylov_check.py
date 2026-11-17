#!/usr/bin/env python3
"""Independent analytic checks for the generic native steady Krylov probe.

Run on Linux after steady_krylov_probe. This checks manufactured finite
dimensional equations, never claims a coupled PDE solution or solver speed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


SX = np.array([1e-3, 3.0, 2e3, 0.02])
B = np.array([1e4, 0.002, 7.0, 0.3])
ROW = np.array([0.25, 2.0, 0.5, 1.5]) / B
M = np.array([[1, .3, -.2, .1], [.2, 1, .4, -.1],
              [-.1, .3, 1, .2], [.4, -.2, .1, 1]])
INITIAL_Y = np.array([.3, 4, .4, 1.8])


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def finite_tree(value, path="input") -> None:
    if isinstance(value, dict):
        for key, item in value.items():
            finite_tree(item, path + "." + key)
    elif isinstance(value, list):
        for index, item in enumerate(value):
            finite_tree(item, f"{path}[{index}]")
    elif isinstance(value, float):
        require(math.isfinite(value), path + " is nonfinite")


def residual(x):
    y = np.asarray(x) / SX
    require(y[1] > 0, "returned state is outside the log domain")
    f = np.array([y[0] ** 2 - 1, math.log(y[1]), y[2] ** 3 - 1,
                  math.exp(y[3]) - math.e])
    return B * (M @ f)


def diagonal(y):
    return np.array([2 * y[0], 1 / y[1], 3 * y[2] ** 2, math.exp(y[3])])


def action(x, direction):
    return B * (M @ (diagonal(np.asarray(x) / SX) * (np.asarray(direction) / SX)))


def relative_error(a, b, scale):
    a, b, scale = np.asarray(a), np.asarray(b), np.asarray(scale)
    return float(np.max(np.abs((a - b) / scale)) /
                 max(1., float(np.max(np.abs(b / scale)))))


def check(data):
    finite_tree(data)
    require(data["schema"] == "WJH27SteadyKrylovProbeV1", "incorrect probe schema")
    require(data["status"] == "RECORDED", "probe did not record its evidence")
    require(data["solver_identity"] == "WJH27SteadyJFNKV1", "incorrect solver identity")
    for key, expected in (("unknown_scale", SX), ("physical_residual_scale", B),
                          ("residual_row_scale", ROW)):
        require(relative_error(data[key], expected, expected) <= 2e-15,
                "changed declared scaling: " + key)
    require(data["acceptance_norm"] == 1e-10, "changed acceptance threshold")
    solves = {case["name"]: case for case in data["solves"]}
    expected_solves = {"validity_rejection", "typed_residual_rejection", "newton_bound",
                       "linear_bound", "line_search_bound", "deadline_before_base",
                       "deadline_after_base", "best_unaccepted_state"}
    require(set(solves) == expected_solves and len(data["solves"]) == len(solves),
            "missing, duplicate or unexpected solve controls")
    maximum_residual_error = maximum_preconditioner_error = 0.
    solve_summary = {}
    for name, case in solves.items():
        r, system = case["result"], case["system"]
        require(case["maximum_seconds"] > 0, name + ": no explicit finite wall budget")
        if system == "coupled_analytic":
            require(relative_error(case["initial"], SX * INITIAL_Y, SX) < 1e-15,
                    name + ": changed initial state")
        else:
            require(system == "no_real_root" and name == "best_unaccepted_state",
                    "unexpected analytic system")
        if r["has_evaluated_state"]:
            if system == "coupled_analytic":
                f, scale = residual(r["x"]), B
            else:
                x = r["x"][0]
                f, scale = np.array([1 - x + .99 * x * x]), np.ones(1)
            err = relative_error(r["residual"], f, scale)
            maximum_residual_error = max(maximum_residual_error, err)
            require(err <= 2e-13, name + ": returned residual is stale or incorrect")
            norm = float(np.max(np.abs(f / scale)))
            require(abs(r["true_norm"] - norm) <= 2e-13 * max(1., norm),
                    name + ": returned true norm is not the fixed physical norm")
            require(r["true_norm"] <= r["last_accepted_norm"] + 2e-13,
                    name + ": returned state is not best evaluated state")
        else:
            require(name == "deadline_before_base" and not r["residual"],
                    name + ": missing evaluated state without pre-base deadline")
        stats = r["statistics"]
        require(all(value >= 0 for value in stats.values()), name + ": negative statistic")
        require(stats["accepted_newton_steps"] == sum(n["accepted"] for n in r["newton"]),
                name + ": accepted history count disagrees")
        require(stats["accepted_newton_steps"] == sum(n["accepted"] for n in r["line_search"]),
                name + ": line-search and Newton histories disagree")
        for sample in case.get("preconditioner_samples", []):
            # Physical RHS -> physical correction. Both independent diagonal
            # residual units and unknown units must be handled exactly once.
            expected = SX * np.asarray(sample["raw_rhs"]) / B / diagonal(sample["prepared_y"])
            err = relative_error(sample["physical_correction"], expected, SX)
            maximum_preconditioner_error = max(maximum_preconditioner_error, err)
            require(err <= 2e-13, name + ": physical right-preconditioner units disagree")
        if name in {"validity_rejection", "typed_residual_rejection"}:
            require(r["converged"] and r["status"] == "CONVERGED" and r["has_evaluated_state"],
                    name + ": analytic nonlinear solve failed")
            require(r["true_norm"] <= data["acceptance_norm"], name + ": true residual failed")
            require(r["best_from_accepted"], name + ": converged state was not accepted")
            require(np.max(np.abs(np.asarray(r["x"]) / SX - 1)) <= 2e-9,
                    name + ": incorrect analytic root")
            require(stats["preconditioner_applications"] > 0 and
                    stats["preconditioner_preparations"] > 0 and case["preconditioner_samples"],
                    name + ": right preconditioning was not exercised")
            require(any(fd["purpose"] == "true-linear-residual" for fd in r["finite_differences"]),
                    name + ": actual accumulated-direction linear check was not exercised")
            rejected = [ls for ls in r["line_search"] if not ls["admissible"]]
            require(rejected and stats["line_search_rejections"] > 0,
                    name + ": feasible line search was not exercised")
            marker = "positive nonlinear coordinate" if name == "validity_rejection" else "typed residual log-domain"
            require(any(marker in ls["reason"] for ls in rejected), name + ": wrong rejection path")
            require(any(ls["accepted"] and ls["step"] < 1 for ls in r["line_search"]),
                    name + ": rejected full step was not followed by an accepted feasible step")
            bases = {nr["iteration"]: nr for nr in r["newton"]}
            for ls in r["line_search"]:
                if ls["accepted"]:
                    threshold = (1 - 1e-4 * ls["step"]) * bases[ls["newton"]]["base_norm"]
                    require(ls["admissible"] and ls["evaluated"] and ls["sufficient_decrease"] and
                            ls["true_norm"] <= max(data["acceptance_norm"], threshold),
                            name + ": accepted trial violates fixed-norm Armijo")
            for lr in r["linear"]:
                if lr["converged"]:
                    require(lr["true_checked"] and
                            lr["true_norm"] <= .001 * bases[lr["newton"]]["base_norm"],
                            name + ": false inner true-residual convergence")
        else:
            require(not r["converged"] and r["status"] == "NOT_CONVERGED",
                    name + ": bounded failure falsely converged")
            prefix = {"newton_bound": "NEWTON_NOT_CONVERGED", "linear_bound": "LINEAR_NOT_CONVERGED",
                      "line_search_bound": "LINE_SEARCH_NOT_CONVERGED", "deadline_before_base": "DEADLINE",
                      "deadline_after_base": "DEADLINE", "best_unaccepted_state": "LINE_SEARCH_NOT_CONVERGED"}[name]
            require(r["reason"].startswith(prefix + ":"), name + ": wrong bounded-failure reason")
            if name == "newton_bound":
                require(stats["accepted_newton_steps"] == 1, "Newton bound did not preserve its accepted update")
            if name == "linear_bound":
                require(stats["linear_cycles"] == 1 and stats["linear_iterations"] == 1 and
                        stats["accepted_newton_steps"] == 0, "declared linear cycle bound was not applied")
            if name == "line_search_bound":
                require(len(r["line_search"]) == 1 and not r["line_search"][0]["admissible"] and
                        stats["accepted_newton_steps"] == 0, "declared line-search bound was not applied")
            if name.startswith("deadline_"):
                require(stats["total_seconds"] >= case["maximum_seconds"], name + ": wall deadline was not reached")
            if name == "deadline_before_base":
                require(not r["has_evaluated_state"] and stats["residual_evaluations"] == 0,
                        "pre-base wall cap manufactured an evaluated solution")
            if name == "deadline_after_base":
                require(r["has_evaluated_state"] and stats["preconditioner_preparations"] == 1 and
                        stats["preparation_seconds"] >= .05 and stats["accepted_newton_steps"] == 0,
                        "cooperative deadline failed to preserve evaluated initial state")
                require(relative_error(r["x"], case["initial"], SX) == 0,
                        "deadline changed initial state without an accepted/evaluated trial")
            if name == "best_unaccepted_state":
                require(not r["best_from_accepted"] and stats["accepted_newton_steps"] == 0 and
                        abs(r["x"][0] - 1) < 1e-8 and abs(r["true_norm"] - .99) < 1e-8 and
                        r["last_accepted_norm"] == 1, "best unaccepted state provenance was lost")
        solve_summary[name] = {"status": r["status"], "reason": r["reason"],
                               "true_norm": r["true_norm"], "has_evaluated_state": r["has_evaluated_state"],
                               "accepted_newton_steps": stats["accepted_newton_steps"]}

    probes = {case["name"]: case for case in data["derivative_probes"]}
    require(set(probes) == {"central", "forward", "backward", "typed_forward", "reduced_central", "infeasible_stencil"}
            and len(probes) == len(data["derivative_probes"]), "incorrect derivative probe inventory")
    maximum_directional_error = maximum_change_error = 0.
    for name, probe in probes.items():
        exact = action(probe["base"], probe["direction"])
        require(relative_error(probe["analytic_action"], exact, B) <= 2e-14,
                name + ": probe's analytic action is incorrect")
        require(probe["multipliers"] == [1, .5, .25], name + ": changed prescribed FD ladder")
        entries = probe["finite_differences"]
        completed_entries = [fd for fd in entries if fd["method"] != "infeasible-reduction"]
        if name == "infeasible_stencil":
            require(not probe["completed"] and not probe["actions"] and not completed_entries and
                    probe["reason"].startswith("FD_INFEASIBLE:") and len(entries) == 3,
                    name + ": declared stencil exhaustion was hidden")
        else:
            require(probe["completed"] and len(probe["actions"]) == 3 and len(completed_entries) == 3,
                    name + ": missing prescribed directional action")
            method = {"central": "central", "forward": "forward-second-order", "backward": "backward-second-order",
                      "typed_forward": "forward-second-order", "reduced_central": "central"}[name]
            require(all(fd["method"] == method for fd in completed_entries), name + ": wrong differentiation rule")
            for value in probe["actions"]:
                err = relative_error(value, exact, B)
                maximum_directional_error = max(maximum_directional_error, err)
                require(err <= 2e-8, name + ": directional derivative disagrees with exact Jacobian")
            changes = []
            for before, after in zip(probe["actions"], probe["actions"][1:]):
                a, b = ROW * np.asarray(after), ROW * np.asarray(before)
                changes.append(float(np.linalg.norm(a - b) / max(np.linalg.norm(a), np.linalg.norm(b))))
            err = relative_error(probe["successive_scaled_changes"], changes, np.ones(2))
            maximum_change_error = max(maximum_change_error, err)
            require(err <= 2e-14, name + ": halving diagnostic ignores residual row scaling")
        coordinate_norm = max(1., float(np.max(np.abs(np.asarray(probe["base"]) / SX))))
        direction_norm = float(np.max(np.abs(np.asarray(probe["direction"]) / SX)))
        action_index = 0
        for fd in entries:
            expected_relative = 1e-5 * probe["multipliers"][action_index] * .5 ** fd["reduction"]
            require(abs(fd["relative_step"] / expected_relative - 1) <= 2e-14,
                    name + ": FD step changed outside its declared ladder")
            require(abs(fd["scaled_direction_norm"] / direction_norm - 1) <= 2e-14 and
                    abs(fd["scaled_coordinate_step"] / (expected_relative * coordinate_norm) - 1) <= 2e-14 and
                    abs(fd["parameter_step"] * direction_norm / (expected_relative * coordinate_norm) - 1) <= 2e-14,
                    name + ": scaled-coordinate FD displacement is wrong")
            method = fd["method"]
            if method == "central":
                require(fd["plus_admissible"] and fd["minus_admissible"], name + ": central stencil is infeasible")
            elif method == "forward-second-order":
                require(fd["plus_admissible"] and not fd["minus_admissible"] and fd["half_admissible"] and
                        fd["minus_reason"], name + ": forward rule lacks actual rejected backward trial")
                if name == "typed_forward":
                    require("typed lower-bound" in fd["minus_reason"], "typed trial exception was not recorded")
            elif method == "backward-second-order":
                require(not fd["plus_admissible"] and fd["minus_admissible"] and fd["half_admissible"] and
                        fd["plus_reason"], name + ": backward rule lacks actual rejected forward trial")
            else:
                require(method == "infeasible-reduction" and not fd["plus_admissible"] and
                        not fd["minus_admissible"] and fd["plus_reason"] and fd["minus_reason"],
                        name + ": derivative reduction lacks actual infeasible trials")
            if method != "infeasible-reduction":
                action_index += 1
        if name == "reduced_central":
            require(any(fd["method"] == "infeasible-reduction" for fd in entries),
                    name + ": derivative reduction was not exercised")

    fatal = {case["name"]: case["reason"] for case in data["fatal_exceptions"]}
    markers = {"fatal_trial_residual": "fatal-residual-marker", "fatal_trial_validity": "fatal-validity-marker",
               "fatal_preconditioner": "fatal-preconditioner-marker", "typed_base_is_fatal": "inadmissible steady base",
               "malformed_residual": "incorrect size", "invalid_true_norm": "TrueNorm returned",
               "missing_wall_cap": "explicit finite wall budget", "invalid_unknown_scale": "positive finite vectors"}
    require(set(fatal) == set(markers) and len(fatal) == len(data["fatal_exceptions"]), "fatal control inventory mismatch")
    for name, marker in markers.items():
        require(marker in fatal[name], name + ": fatal programming/input error was swallowed or misclassified")
    return {"maximum_returned_residual_error": maximum_residual_error,
            "maximum_physical_preconditioner_error": maximum_preconditioner_error,
            "maximum_analytic_directional_error": maximum_directional_error,
            "maximum_halving_diagnostic_error": maximum_change_error,
            "solve_cases": solve_summary, "derivative_cases": len(probes), "fatal_controls": len(fatal),
            "scope": "Generic finite-dimensional solver only; no coupled PDE convergence, physical qualification or performance claim",
            "wall_cap_scope": "Cooperative deadline; the intentionally sleeping callback cannot be interrupted"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("immutable checker output already exists")
    raw = args.input.read_bytes()
    report = {"schema": "WJH27SteadyKrylovCheckV1", "status": "FAIL",
              "input_sha256": hashlib.sha256(raw).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    try:
        report.update(check(json.loads(raw)))
        report["status"] = "PASS"
    except Exception as error:
        report["error"] = f"{type(error).__name__}: {error}"
    with args.output.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(report, allow_nan=False))
    if report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
