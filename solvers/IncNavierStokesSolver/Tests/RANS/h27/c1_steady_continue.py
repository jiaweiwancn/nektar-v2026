#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Continue diagnostic steady roots on rebuilt normal charts with audit evidence."""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from resources import identity


HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text())


def call(command, log, seconds, cwd):
    wrapped = ["/usr/bin/time", "-v", "-o", str(log)+".cost", "timeout",
               "--signal=TERM", "--kill-after=10s", f"{seconds:.3f}s", *map(str, command)]
    with log.open("x") as stream:
        return subprocess.run(wrapped, cwd=cwd, stdin=subprocess.DEVNULL, stdout=stream,
                              stderr=subprocess.STDOUT, check=False).returncode


def stage(args, case, seed, height, increment, deadline):
    starts = time.monotonic()
    available = deadline-starts
    cap = min(args.stage_seconds, available-35.)
    prepare = [sys.executable, HERE/"c1_steady_case.py", args.template, case, "--seed", seed,
               "--remap", "--height", height, "--quadrature", *args.quadrature, "--seconds", cap,
               "--eta", args.eta, "--regularization-viscosity", args.viscosity,
               "--regularization-radius", *args.regularization_radius,
               "--pseudo-time", args.pseudo_time, "--beta", args.beta,
               "--scalar-subdivision", args.scalar_subdivision]
    row = dict(height=height, increment=increment, seed=identity(seed), case=str(case),
               eta=args.eta, predicted="geometry predictor followed by the current fixed-height steady residual",
               commands=[list(map(str, prepare))])
    status = call(prepare, args.out/(case.name+".setup.log"), min(30., available), args.out)
    if status:
        row.update(decision="setup failure", status=status, seconds=time.monotonic()-starts)
        return row, None
    native = [args.solver, "run.xml"]
    row["commands"].append(list(map(str, native)))
    row["session"] = identity(case/"run.xml")
    row["status"] = call(native, case/"run.log", min(cap+30., deadline-time.monotonic()), case)
    ledger = case/"state.steady.json"
    candidate = case/"state.diagnostic-steady.json"
    if ledger.is_file():
        nonlinear = read(ledger)
        row.update(reason=nonlinear["reason"], original_or_modified_true_norm=nonlinear["true_norm"],
                   nonlinear_steps=nonlinear["accepted_newton_steps"])
    if candidate.is_file():
        state = read(candidate)
        row.update(candidate=identity(candidate), original_residuals=state["original_residuals"],
                   solved_residuals=state["residuals"], original_margin=state["original_principal_over_mu"],
                   modified_margin=state["modified_principal_over_original_mu"],
                   K_over_Kin=state["state"]["interface_K_over_Kin"])
        for program, name in (("c1_oracle.py", "equations"), ("c1_margins.py", "margins")):
            for regularized in ([False, True] if args.eta > 0. else [False]):
                audit_cap = min(args.audit_seconds, deadline-time.monotonic()-5.)
                if audit_cap <= 5.:
                    break
                kind = name+("-modified" if regularized else "")
                path = case/(kind+".json")
                check = [sys.executable, HERE/program, "--native", candidate, "--diagnostic",
                         "--table", args.table, "--out", path, "--order", args.audit_order,
                         "--seconds", audit_cap]
                if regularized:
                    check.append("--regularized")
                row["commands"].append(list(map(str, check)))
                code = call(check, path.with_suffix(".log"), audit_cap+4., case)
                row[kind] = identity(path) if code == 0 else {"exit_status": code}
    row["seconds"] = time.monotonic()-starts
    equation_report = case/("equations-modified.json" if args.eta > 0. else "equations.json")
    if row["status"] != 0 or not candidate.is_file():
        row["decision"] = "reduce increment after failed steady solve"
        return row, None
    if not equation_report.is_file() or not read(equation_report)["residual_pass"]:
        row["decision"] = "hold for independent quadrature check"
        return row, None
    if not (case/"margins.json").is_file() or (args.eta > 0. and not (case/"margins-modified.json").is_file()):
        row["decision"] = "hold for independent principal-response check"
        return row, None
    row["decision"] = "continue as diagnostic seed"
    return row, candidate


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("template", "out", "solver", "seed", "table"):
        p.add_argument("--"+key, required=True, type=lambda s: Path(s).resolve())
    p.add_argument("--target", type=float, required=True)
    p.add_argument("--increment", type=float, required=True)
    p.add_argument("--minimum", type=float, required=True)
    p.add_argument("--maximum", type=float, required=True)
    p.add_argument("--quadrature", type=int, nargs=2, required=True)
    p.add_argument("--audit-order", type=int, default=20)
    p.add_argument("--stage-seconds", type=float, required=True)
    p.add_argument("--audit-seconds", type=float, default=180.)
    p.add_argument("--campaign-seconds", type=float, required=True)
    p.add_argument("--eta", type=float, default=0.)
    p.add_argument("--viscosity", type=float, default=0.)
    p.add_argument("--regularization-radius", type=float, nargs=2, default=(0.,0.))
    p.add_argument("--pseudo-time", type=float, default=0.)
    p.add_argument("--beta", type=float, default=0.)
    p.add_argument("--scalar-subdivision", type=int, choices=range(7), default=0)
    args = p.parse_args()
    seed, bound = args.seed, read(args.seed)
    if (bound["schema"] != "WJH27DiagnosticCandidateV1" or bound["role"] != "STEADY_SEARCH"
            or max(bound["residuals"]) > 1e-10):
        p.error("a converged diagnostic steady seed is required")
    g = bound["state"]["geometry"]["hill_height"]
    if not (0 < args.minimum <= args.increment <= args.maximum and g < args.target <= 1):
        p.error("invalid continuation heights or increments")
    args.out.mkdir(parents=True)
    deadline = time.monotonic()+args.campaign_seconds
    delta, number = args.increment, 0
    with (args.out/"attempts.jsonl").open("x") as journal:
        while g < args.target and delta >= args.minimum and deadline-time.monotonic()>60.:
            number += 1
            height = min(args.target, g+delta)
            row, candidate = stage(args, args.out/f"{number:03d}-height-{height:.9g}",
                                   seed, height, delta, deadline)
            row["solver"] = identity(args.solver)
            journal.write(json.dumps(row, allow_nan=False)+"\n")
            journal.flush()
            os.fsync(journal.fileno())
            if candidate:
                seed, g = candidate, height
                factor = 1.7 if row["nonlinear_steps"] <= 7 else 1.2
                delta = min(args.maximum, factor*delta)
            elif row["decision"].startswith("hold"):
                break
            else:
                delta *= .5
    print(json.dumps(dict(reached_height=g, last_seed=str(seed), target=args.target)))


if __name__ == "__main__":
    main()
