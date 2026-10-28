#!/usr/bin/env python3
"""Check raw solver logs for the small HITInflow Milestone 3 fixtures.

No error is rounded, reference-subtracted, or silently discarded. Run this on
Linux solver logs; generating fixtures alone is not numerical validation.
"""

import argparse
import json
import math
from pathlib import Path
import re
import sys
import xml.etree.ElementTree as ET


VARIABLES = ("u", "v", "w", "p")
ERROR_PATTERN = re.compile(
    r"L\s+(2|inf)\s+error\s*\(variable\s+([^)]+)\)\s*:\s*(\S+)", re.IGNORECASE)


def parse_errors(path, expected_steps, expected_time):
    text = path.read_text(encoding="utf-8", errors="replace")
    if "Total Computation Time" not in text:
        raise ValueError(f"{path}: missing completed-solver marker")
    status = re.findall(r"Steps:\s*(\d+)\s+Time:\s*(\S+)", text)
    if not status:
        raise ValueError(f"{path}: missing final timestep status")
    steps, time = int(status[-1][0]), float(status[-1][1])
    if steps != expected_steps or not math.isclose(time, expected_time, abs_tol=1e-12):
        raise ValueError(f"{path}: final status is {steps} steps, time {time}; "
                         f"expected {expected_steps} steps, time {expected_time}")
    found = {"l2": {}, "linf": {}}
    for metric, var, token in ERROR_PATTERN.findall(text):
        value = float(token)
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"{path}: invalid {metric} error for {var}: {token}")
        found["l2" if metric == "2" else "linf"][var.strip()] = value
    for metric, values in found.items():
        missing = set(VARIABLES) - values.keys()
        if missing:
            raise ValueError(f"{path}: missing {metric} errors for {sorted(missing)}")
    return found


def parameters(path):
    root = ET.parse(path).getroot()
    values = {}
    for element in root.findall("./CONDITIONS/PARAMETERS/P"):
        name, value = element.text.split("=", 1)
        values[name.strip().lower()] = float(value)
    return values


def temporal(args):
    if len(args.run) < 3:
        raise ValueError("Temporal convergence needs at least three --run SESSION LOG pairs")
    runs = []
    for session, log in args.run:
        params = parameters(session)
        errors = parse_errors(log, round(params["numsteps"]),
                              params["timestep"] * params["numsteps"])
        vector_error = math.sqrt(sum(errors["l2"][var]**2 for var in VARIABLES[:3]))
        if vector_error <= 0.0:
            raise ValueError(f"{log}: zero vector error cannot establish temporal convergence")
        runs.append({"session": str(session), "log": str(log),
                     "dt": params["timestep"],
                     "final_time": params["timestep"] * params["numsteps"],
                     "l2": errors["l2"], "linf": errors["linf"],
                     "velocity_vector_l2": vector_error})
    runs.sort(key=lambda item: item["dt"], reverse=True)
    orders = []
    for coarse, fine in zip(runs, runs[1:]):
        if not math.isclose(coarse["final_time"], fine["final_time"], abs_tol=1e-13):
            raise ValueError("All manufactured runs must finish at the same physical time")
        ratio = coarse["dt"] / fine["dt"]
        if not math.isclose(ratio, 2.0, abs_tol=1e-12):
            raise ValueError("Use consecutive timestep halvings for temporal convergence")
        order = math.log(coarse["velocity_vector_l2"] / fine["velocity_vector_l2"]) / math.log(ratio)
        orders.append(order)
    result = {"kind": "temporal", "runs": runs,
              "velocity_vector_l2_orders": orders, "minimum_order": args.min_order,
              "maximum_coarse_velocity_vector_l2": args.max_coarse_error}
    result["passed"] = (all(order >= args.min_order for order in orders)
                        and runs[0]["velocity_vector_l2"] <= args.max_coarse_error)
    return result


def smoke(args):
    errors = parse_errors(args.log, 20, 0.02)
    result = {"kind": args.command, "log": str(args.log), **errors}
    if args.command == "mean-only":
        result["maximum_error"] = args.max_error
        result["passed"] = max(value for row in errors.values() for value in row.values()) <= args.max_error
    else:
        # Stochastic fixture has no ExactSolution: the printed errors are norms
        # relative to zero. These checks establish a finite, nonzero smoke run,
        # not agreement with an exact flow or an experimental spectrum.
        volume = (2.0 * math.pi)**2
        result["norm_reference"] = "zero; no exact stochastic solver solution asserted"
        result["maximum_norm"] = args.max_norm
        result["passed"] = (
            0.5 * math.sqrt(volume) < errors["l2"]["u"] < 2.0 * math.sqrt(volume)
            and errors["l2"]["v"] > 1e-12 and errors["l2"]["w"] > 1e-12
            and max(value for row in errors.values() for value in row.values()) < args.max_norm)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, help="Save the complete unrounded acceptance result")
    commands = parser.add_subparsers(dest="command", required=True)
    convergence = commands.add_parser("temporal", help="Check raw final-time L2 errors")
    convergence.add_argument("--run", nargs=2, type=Path, action="append", required=True,
                             metavar=("SESSION", "LOG"))
    convergence.add_argument("--min-order", type=float, default=1.8)
    convergence.add_argument("--max-coarse-error", type=float, default=5e-5)
    mean = commands.add_parser("mean-only", help="Check zero-RMS exact mean preservation")
    mean.add_argument("log", type=Path)
    mean.add_argument("--max-error", type=float, default=1e-10)
    stochastic = commands.add_parser("stochastic", help="Check a nonzero stable open-outlet smoke run")
    stochastic.add_argument("log", type=Path)
    stochastic.add_argument("--max-norm", type=float, default=20.0)
    args = parser.parse_args()
    try:
        result = temporal(args) if args.command == "temporal" else smoke(args)
    except (OSError, ValueError, KeyError, ET.ParseError) as error:
        parser.error(str(error))
    output = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.json:
        args.json.write_text(output, encoding="utf-8")
    print(output, end="")
    if not result["passed"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
