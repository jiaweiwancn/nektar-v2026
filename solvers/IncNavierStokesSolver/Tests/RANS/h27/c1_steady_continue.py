#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Continue independently admitted steady roots, with fixed-height deadline recovery."""

import argparse
import json
import math
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


def accepted_packet(path, height):
    state = read(path)
    if not (state.get("schema") == "WJH27DiagnosticCandidateV1" and
            state.get("role") == "STEADY_SEARCH" and
            state.get("iterate_custody") == "LAST_ACCEPTED" and
            state.get("physical_history") is False and
            abs(state["state"]["geometry"]["hill_height"]-height) <= 2e-14):
        raise ValueError("missing typed last accepted state at the requested height")
    return state


def termination(code, ledger, candidate):
    """Classify native work before running independent admission checks."""
    if code == 0 and ledger and ledger.get("converged") and candidate:
        return "AUDIT_REQUIRED"
    if not ledger or not candidate:
        return "DEADLINE_WITHOUT_CUSTODY" if code in (124, 137) else "SETUP_OR_CUSTODY_FAILURE"
    reason = ledger["reason"]
    if reason.startswith("DEADLINE"):
        progress = (ledger.get("accepted_newton_steps", 0) > 0 and ledger.get("newton") and
                    ledger["last_accepted_norm"] < ledger["newton"][0]["base_norm"])
        return "BUDGET_EXHAUSTED_WITH_PROGRESS" if progress else "BUDGET_EXHAUSTED"
    if reason.startswith(("LINEAR_", "SPARSE_FACTOR", "SPARSE_SOLVE", "SPARSE_BACKSOLVE",
                          "SPARSE_FD_", "SPARSE_SUPPORT", "FD_", "ARNOLDI_")):
        return "CURRENT_LINEAR_FAILURE"
    if reason.startswith(("LINE_SEARCH_", "NEWTON_DIRECTION_", "NEWTON_ROUNDOFF")):
        return "GLOBALIZATION_FAILURE"
    if reason.startswith("NEWTON_NOT_CONVERGED"):
        return "NEWTON_ITERATION_BUDGET"
    return "SETUP_OR_CUSTODY_FAILURE"


def original_admission(candidate, table, receipts):
    """Apply the original scientific gates to individually digest-bound evidence."""
    state = read(candidate)
    equations, margins, points = (read(receipts[name])
                                  for name in ("equations", "margins", "points"))
    expected_packet = identity(candidate)["sha256"]
    expected_table = identity(table)["sha256"]
    if any(r["native"]["sha256"] != expected_packet or r["table"]["sha256"] != expected_table
           for r in (equations, margins, points)):
        return "AUDIT_IDENTITY_FAILURE"
    native = state["state"]
    if (equations["source_identity"] != native["source_identity"] or
            points["coefficient_sha256"] != native["state_sha256"] or
            any(r["candidate_role"] != "STEADY_SEARCH" or r["regularization"] is not None or
                r["step"] != 0 or r["time"] != 0. or r["order"] <= native["quadratureX"] or
                r["radial_order"] <= native["quadratureR"] for r in (equations, margins))):
        return "AUDIT_IDENTITY_FAILURE"
    if not (state["regularization"]["eta"] == 0. and not state["physical_history"] and
            state.get("iterate_custody") != "UNCOMMITTED_TRIAL" and
            max(state["original_residuals"]) <= 1e-10 and equations["residual_pass"] and
            all(abs(equations["residuals"][key]) <= limit for key, limit in
                (("momentum", 1e-9), ("continuity", 1e-9), ("k", 1e-8), ("epsilon", 1e-8))) and
            equations["epsilon_trace"]["interior_projection_defect"] <= 1e-10 and
            abs(equations["power"]["algebraic_error"]) <= 1e-10 and
            abs(equations["mass_integration_identity_error"]) <= 1e-10 and
            abs(equations["interface"]["physical_face_measure_max_difference"]) <= 1e-10 and
            abs(equations["interface"]["shared_trace_max"]) <= 1e-10):
        return "INDEPENDENT_EQUATION_FAILURE"
    if not (min(state["physical_whole_cell_scalar_lower_bounds"]) > 0. and
            state["minimum_covariance_over_k"] > 1e-10 and
            min(native["interface_K_table_margin_over_Kin"]) >= 0. and
            equations["integration_minimum_jacobian"] > 0. and
            equations["raw_boundary"]["outlet_minimum_normal_velocity"] >= -1e-12):
        return "ORIGINAL_DOMAIN_FAILURE"
    minima = points["local_minima"]
    if (not margins["one_sided_faces"] or not
            {"inner", "outer", "matching_inner", "matching_outer", "faces_inner",
             "faces_outer", "cutoff", "forcing_support"} <= margins["regions"].keys() or
            not {"inner", "outer", "faces_inner", "faces_outer", "forcing_support"} <=
            {record["region"] for record in minima} or
            any(record["receipt"]["sha256"] != identity(receipts["margins"])["sha256"]
                for record in minima)):
        return "AUDIT_IDENTITY_FAILURE"
    if (state["original_principal_over_mu"] <= 1e-8 or
            not all(region["sample_pass"] for region in margins["regions"].values()) or
            any(not record["values"] or any(
                value["differentiation_spread"] < 0. or len(value["minima"]) < 3 or
                any(m["over_mu"]-value["differentiation_spread"] <= 1e-8 for m in value["minima"])
                for value in record["values"]) for record in minima)):
        return "ADVERSE_CONVERGED_ROOT"
    return "QUALIFIED"


def run_audits(args, case, candidate, deadline, row):
    requests = [("equations", "c1_oracle.py", False), ("margins", "c1_margins.py", False)]
    if args.eta > 0.:
        requests.extend((("equations-modified", "c1_oracle.py", True),
                         ("margins-modified", "c1_margins.py", True)))
    if args.require_original_admissible:
        requests.append(("points", "c1_refined_witness.py", False))
    paths = {}
    for kind, tool, modified in requests:
        cap = min(args.audit_seconds, deadline-time.monotonic()-10.)
        if cap <= 5.:
            return "HOLD_FOR_AUDIT"
        path = case/(kind+".json")
        check = [sys.executable, HERE/tool]
        if kind == "points":
            check.append("inspect")
        else:
            check.append("--diagnostic")
        check.extend(("--native", candidate, "--table", args.table, "--out", path,
                      "--order", args.audit_order, "--seconds", cap))
        if modified:
            check.append("--regularized")
        if args.require_original_admissible and kind == "margins":
            check.append("--one-sided-faces")
        if kind == "points":
            for a, r in [(-8.29, .1111), (2.5, 1.107899), (1.875, .11541061001735001)]:
                check.extend(("--point", a, r))
            for region in ("inner", "outer", "faces_inner", "faces_outer", "forcing_support"):
                check.extend(("--minima", region, paths["margins"]))
        row["commands"].append(list(map(str, check)))
        status = call(check, path.with_suffix(".log"), cap+4., case)
        if status or not path.is_file():
            row[kind] = {"exit_status": status}
            return "HOLD_FOR_AUDIT"
        paths[kind] = path
        row[kind] = identity(path)
    if args.require_original_admissible:
        return original_admission(candidate, args.table, paths)
    source = paths["equations-modified" if args.eta > 0. else "equations"]
    return "QUALIFIED" if read(source)["residual_pass"] else "INDEPENDENT_EQUATION_FAILURE"


def stage(args, case, seed, last_qualified, height, increment, deadline, resume=False):
    starts = time.monotonic()
    available = deadline-starts
    cap = min(args.resume_seconds if resume else args.stage_seconds, available-35.)
    prepare = [sys.executable, HERE/"c1_steady_case.py", args.template, case, "--seed", seed,
               "--height", height, "--quadrature", *args.quadrature, "--seconds", cap,
               "--eta", args.eta, "--regularization-viscosity", args.viscosity,
               "--regularization-radius", *args.regularization_radius,
               "--pseudo-time", args.pseudo_time, "--beta", args.beta,
               "--scalar-subdivision", args.scalar_subdivision, "--inverse", args.inverse,
               "--maximum-inverse-age", args.maximum_inverse_age,
               "--reuse-linear-iterations", args.reuse_linear_iterations]
    if not resume:
        prepare.append("--remap")
    row = dict(height=height, increment=increment, seed=identity(seed), case=str(case),
               last_qualified=identity(last_qualified), eta=args.eta, same_height_resume=resume,
               requested_native_seconds=cap, commands=[list(map(str, prepare))])
    status = call(prepare, args.out/(case.name+".setup.log"), min(30., available), args.out)
    if status:
        row.update(decision="SETUP_OR_CUSTODY_FAILURE", status=status, seconds=time.monotonic()-starts)
        return row, None
    native = [args.solver, "run.xml"]
    row["commands"].append(list(map(str, native)))
    row["session"] = identity(case/"run.xml")
    row["status"] = call(native, case/"run.log", min(cap+30., deadline-time.monotonic()), case)
    ledger = case/"state.steady.json"
    candidate = case/"state.diagnostic-steady.json"
    trial = case/"state.diagnostic-uncommitted.json"
    nonlinear = read(ledger) if ledger.is_file() else None
    if nonlinear:
        row.update(ledger=identity(ledger), reason=nonlinear["reason"],
                   best_from_accepted_iterate=nonlinear["best_from_accepted_iterate"],
                   best_evaluated_norm=nonlinear["true_norm"], last_accepted_norm=nonlinear["last_accepted_norm"],
                   nonlinear_steps=nonlinear["accepted_newton_steps"])
    if trial.is_file():
        row["uncommitted_trial"] = identity(trial)
    try:
        state = accepted_packet(candidate, height) if candidate.is_file() else None
    except ValueError as error:
        row.update(decision="SETUP_OR_CUSTODY_FAILURE", reason=str(error), seconds=time.monotonic()-starts)
        return row, None
    if state:
        if (state["regularization"]["eta"] != args.eta or
                (state["state"]["quadratureX"], state["state"]["quadratureR"]) != tuple(args.quadrature)):
            row.update(decision="SETUP_OR_CUSTODY_FAILURE", reason="native equation controls differ",
                       seconds=time.monotonic()-starts)
            return row, None
        row.update(last_accepted_candidate=identity(candidate), original_residuals=state["original_residuals"],
                   solved_residuals=state["residuals"], original_margin=state["original_principal_over_mu"],
                   modified_margin=state["modified_principal_over_original_mu"],
                   K_over_Kin=state["state"]["interface_K_over_Kin"])
    row["decision"] = termination(row["status"], nonlinear, state)
    if row["decision"] == "AUDIT_REQUIRED":
        row["decision"] = run_audits(args, case, candidate, deadline, row)
    row["seconds"] = time.monotonic()-starts
    return row, candidate if state else None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("template", "out", "solver", "seed", "table"):
        p.add_argument("--"+key, required=True, type=lambda s: Path(s).resolve())
    p.add_argument("--seed-receipts", type=Path, nargs=3, metavar=("EQUATIONS", "MARGINS", "POINTS"))
    p.add_argument("--resume-pending", type=lambda s: Path(s).resolve())
    p.add_argument("--target", type=float, required=True)
    p.add_argument("--increment", type=float, required=True)
    p.add_argument("--minimum", type=float, required=True)
    p.add_argument("--maximum", type=float, required=True)
    p.add_argument("--quadrature", type=int, nargs=2, required=True)
    p.add_argument("--audit-order", type=int, default=28)
    p.add_argument("--stage-seconds", type=float, required=True)
    p.add_argument("--resume-seconds", type=float, default=0.)
    p.add_argument("--audit-seconds", type=float, default=180.)
    p.add_argument("--campaign-seconds", type=float, required=True)
    p.add_argument("--eta", type=float, default=0.)
    p.add_argument("--viscosity", type=float, default=0.)
    p.add_argument("--regularization-radius", type=float, nargs=2, default=(0.,0.))
    p.add_argument("--pseudo-time", type=float, default=0.)
    p.add_argument("--beta", type=float, default=0.)
    p.add_argument("--scalar-subdivision", type=int, choices=range(7), default=0)
    p.add_argument("--inverse", choices=("BlockV1", "ColoredSparseV1"), default="BlockV1")
    p.add_argument("--maximum-inverse-age", type=int, default=0)
    p.add_argument("--reuse-linear-iterations", type=int, default=12)
    p.add_argument("--require-original-admissible", action="store_true")
    args = p.parse_args()
    if not (all(math.isfinite(value) and value > 0.
                for value in (args.stage_seconds, args.audit_seconds, args.campaign_seconds)) and
            math.isfinite(args.resume_seconds) and args.resume_seconds >= 0.):
        p.error("positive finite native, audit and campaign budgets are required")
    seed, bound = args.seed, read(args.seed)
    if (bound["schema"] != "WJH27DiagnosticCandidateV1" or bound["role"] != "STEADY_SEARCH"
            or bound.get("iterate_custody") == "UNCOMMITTED_TRIAL" or max(bound["residuals"]) > 1e-10):
        p.error("a converged accepted diagnostic steady seed is required")
    if args.require_original_admissible:
        if args.eta != 0. or bound["regularization"]["eta"] != 0.:
            p.error("original admissible continuation requires literal eta=0")
        if not args.seed_receipts:
            p.error("supply the independent original equation, response and point receipts for the anchor")
        receipts = dict(zip(("equations", "margins", "points"), args.seed_receipts))
        decision = original_admission(seed, args.table, receipts)
        if decision != "QUALIFIED":
            p.error("anchor independent admission failed: "+decision)
    g = bound["state"]["geometry"]["hill_height"]
    if not (0 < args.minimum <= args.increment <= args.maximum and g < args.target <= 1):
        p.error("invalid continuation heights or increments")
    pending = args.resume_pending
    if pending:
        h = read(pending)["state"]["geometry"]["hill_height"]
        packet = accepted_packet(pending, h)
        if not (g < h <= args.target and args.resume_seconds > 0. and
                packet["regularization"]["eta"] == args.eta and
                (packet["state"]["quadratureX"], packet["state"]["quadratureR"]) ==
                tuple(args.quadrature)):
            p.error("explicit same-height resume must lie above the qualified anchor within the target")
    args.out.mkdir(parents=True)
    deadline = time.monotonic()+args.campaign_seconds
    delta, number, decision = args.increment, 0, "NO_STAGE"
    solver = identity(args.solver)
    with (args.out/"attempts.jsonl").open("x") as journal:
        while g < args.target and delta >= args.minimum and deadline-time.monotonic()>60.:
            number += 1
            resume = pending is not None
            height = h if resume else min(args.target, g+delta)
            row, candidate = stage(args, args.out/f"{number:03d}-height-{height:.9g}",
                                   pending if resume else seed, seed, height, height-g, deadline, resume)
            row["solver"] = solver
            journal.write(json.dumps(row, allow_nan=False)+"\n")
            journal.flush()
            os.fsync(journal.fileno())
            decision = row["decision"]
            if decision == "QUALIFIED":
                last_delta = height-g
                seed, g, pending = candidate, height, None
                uses = row["nonlinear_steps"]
                factor = 1.5 if uses <= 4 else 1.15 if uses <= 8 else .8
                delta = max(args.minimum, min(args.maximum, factor*last_delta))
            elif (decision == "BUDGET_EXHAUSTED_WITH_PROGRESS" and not resume and
                  args.resume_seconds > 0. and deadline-time.monotonic()>60.):
                pending, h = candidate, height
            else:
                pending = candidate if candidate else pending
                break
    print(json.dumps(dict(reached_height=g, last_seed=str(seed), last_qualified=identity(seed),
                          pending=identity(pending) if pending else None, decision=decision,
                          target=args.target), allow_nan=False))


if __name__ == "__main__":
    main()
