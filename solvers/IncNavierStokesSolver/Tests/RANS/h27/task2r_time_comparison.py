#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare registered Bf trajectories at common times and verify cold histories."""

import argparse
from pathlib import Path
import time

import physical_oracle as po
from physical_cache_history_check import binary_history
from physical_transient_oracle import state_hash
from resources import atomic, identity
from task2r_advancing import DEFAULT_FD_STEP, DEFAULT_STEP_PC_MAXIMUM_USES, read, require, verified
from task2r_advancing_evidence import HASH_KEYS
from task2s_observations import compare


def accepted(case, step, allow_interrupted=False):
    index_path = case / "state.science-index.json"
    index = read(index_path)
    require((index["run_complete"] or allow_interrupted) and step > index["initial_step"],
            "comparison needs a process-accepted positive-time frame")
    frame, = (f for f in index["frames"] if f["step"] == step)
    path = verified({k:frame[k] for k in ("path","sha256","bytes")})
    data, _ = po.load_native(path)
    require(path == case / f"state.science-{step}.json"
            and frame["origin"] == "current_process_accepted_step"
            and data["step"] == step and data["time"] == frame["time"]
            and data["state_sha256"] == state_hash(data["coefficients"]),
            "scientific frame or coefficient identity differs")
    return read(case/"admission.json"), data, identity(path), identity(index_path)


def common_time(coarse, fine, coarse_step, fine_step, out, seconds, coarse_audit=None, fine_audit=None):
    deadline = time.monotonic()+seconds
    a, sa, ba, ia = accepted(coarse,coarse_step,allow_interrupted=coarse_audit is not None)
    b, sb, bb, ib = accepted(fine,fine_step,allow_interrupted=fine_audit is not None)
    inputs = ("role","source","application","geometry","table","initial_fields","certificate",
              "represented_proof","step_preconditioner")
    require(a["role"] in ("startup","perturbed") and all(a.get(k) == b.get(k) for k in inputs)
            and a.get("fd_step",DEFAULT_FD_STEP) == b.get("fd_step",DEFAULT_FD_STEP)
            and a.get("step_pc_maximum_uses", DEFAULT_STEP_PC_MAXIMUM_USES)
            == b.get("step_pc_maximum_uses", DEFAULT_STEP_PC_MAXIMUM_USES)
            and sa["elements"] == sb["elements"] and a["timestep"] >= b["timestep"]
            and abs(sa["time"]-sb["time"]) <= 2e-13*max(1.,sa["time"]),
            "trajectory or exact common-time binding differs")
    audits = {}
    for role, path, frame, plan in (("coarse",coarse_audit,ba,a), ("fine",fine_audit,bb,b)):
        if path is not None:
            audit = equation_audit(path,frame,plan["table"])
            require(audit["equations"] == "MET", f"{role} common-time frame failed its independent audit")
            audits[role] = audit
    reference = po.load_reference(verified(a["table"]))
    p,q = (po.PhysicalOracle(state,reference,20,deadline) for state in (sa,sb))
    rows, interior, inventory = compare(p,q,interior=(-2.5,8),detailed=True)
    failures = [key for group in (rows,interior) for key,row in group.items() if row["status"] != "MET"]
    result = dict(schema="WJH27R3TimeComparisonV1",status="MET" if not failures else "NOT_MET",
        role=a["role"],time=sa["time"],timesteps=[sa["dt"],sb["dt"]],steps=[coarse_step,fine_step],
        frame_hashes=[ba,bb],indices=[ia,ib],table=a["table"],inventory=inventory,independent_order=20,
        original_rows=rows,fixed_window_rows=interior,failures=failures)
    if audits:
        result["endpoint_equation_audits"] = audits
        result["complete_intervals"] = [read(verified(index))["run_complete"] for index in (ia,ib)]
    atomic(out,result)


def exact_clock(packet, boundary):
    return all(packet[k] == boundary[k] for k in ("step","time","dt","signature"))


def equation_audit(path, native, table):
    result = read(path)
    require(result["schema"] == "WJH27IndependentPhysicalTransientOracleV1"
            and all(result[name][key] == bound[key] for name,bound in (("native",native),("table",table))
                    for key in ("path","sha256")), "independent equation audit has different scientific inputs")
    return dict(artifact=identity(path), equations=result.get("accepted_time_discrete_equations_status","NOT_MET"))


def evolution(case, out, seconds, recovered_prefix=None, accepted_audit=None):
    deadline = time.monotonic()+seconds
    index = read(case/"state.science-index.json")
    interrupted = not index["run_complete"]
    require(not interrupted or (recovered_prefix is not None and accepted_audit is not None),
            "an interrupted prefix needs a verified recovered start and an audited accepted endpoint")
    plan, final, final_file, index_file = accepted(case,index["final_step"],allow_interrupted=interrupted)
    audit = equation_audit(accepted_audit,final_file,plan["table"]) if accepted_audit else None
    require(not interrupted or audit["equations"] == "MET", "interrupted endpoint audit did not pass")
    role = plan["role"]
    if recovered_prefix is None:
        frame, = (row for row in index["frames"] if row["origin"] == "fresh_initial_observation")
        start = verified({k:frame[k] for k in ("path","sha256","bytes")})
    else:
        anchor = read(recovered_prefix)
        checkpoint = anchor["recovery"]
        require(anchor["schema"] == "WJH27R3RecoveredPrefixAnchorV1" and anchor["status"] == "MET"
                and role == "restart" and checkpoint["step"] == index["initial_step"]
                and plan["restart"] == {k:checkpoint[k] for k in ("path","sha256","bytes")},
                "evolution and recovered prefix anchors differ")
        start = verified(anchor["initial"])
        role = "perturbed"
    initial, _ = po.load_native(start)
    require(initial["step"] == initial["time"] == 0 and initial["elements"] == final["elements"]
            and initial["signature"] == final["signature"], "initial and accepted inventories differ")
    ref = po.load_reference(verified(plan["table"]))
    before, after = (po.PhysicalOracle(data,ref,20,deadline) for data in (initial,final))
    rows, window, inventory = compare(after,before,interior=(-2.5,8),detailed=True)
    atomic(out,dict(schema="WJH27R3ObservedEvolutionV1",role=role,
        initial=identity(start),accepted=final_file,scientific_index=index_file,
        interrupted=interrupted,accepted_equation_audit=audit,
        time=final["time"],inventory=inventory,original_rows=rows,fixed_window_rows=window,
        interpretation="Signed final-minus-initial evolution. Comparison tolerances provide scale "
                       "context for measured evolution and do not govern this change."))


def cold(uninterrupted, resumed, split, checkpoint, out, seconds, uninterrupted_audit=None,
         resumed_audit=None):
    deadline = time.monotonic()+seconds
    parent = read(uninterrupted/"state.science-index.json")
    child = read(resumed/"state.science-index.json")
    first = read(split/"state.science-index.json")
    pa, p, pf, pi = accepted(uninterrupted,parent["final_step"])
    ca, c, cf, ci = accepted(resumed,child["final_step"])
    sa, s, sf, si = accepted(split,first["final_step"])
    invariant = ("geometry", "table", "application", "source", "step_preconditioner")
    initial = ("role", "initial_fields", "certificate", "represented_proof")
    require(p["step"] == c["step"] > s["step"] >= 2 and exact_clock(p,c)
            and all(pa[k] == ca[k] == sa[k] for k in invariant)
            and pa.get("fd_step",DEFAULT_FD_STEP) == ca.get("fd_step",DEFAULT_FD_STEP)
            == sa.get("fd_step",DEFAULT_FD_STEP)
            and all(pa.get(k) == sa.get(k) for k in initial)
            and p["signature"] == s["signature"] and ca["restart"] == identity(checkpoint)
            and checkpoint == split/f"state.recovery-{s['step']}.bin"
            and s["step"] % sa["recovery_stride"] == 0,
            "restart misses the fixed clock, run signature or on-cadence recovery")
    ref = po.load_reference(verified(pa["table"]))
    imported_files = (f for f in child["frames"] if f["origin"] == "imported_accepted_history")
    imported, = imported_files
    path = verified({k:imported[k] for k in ("path","bytes","sha256")})
    imported_state, _ = po.load_native(path)
    _, recovery = binary_history(checkpoint,s,s["dt"],deadline)
    histories_ok = exact_clock(imported_state,s) and exact_clock(recovery,s) and all(
        recovery["hashes"][role] == s[digest] == imported_state[digest] for role,digest in HASH_KEYS.items())
    run = read(resumed/"state.result.json")
    require(run["restart"] and run["normal_main"] and not run["steady_mode"], "cold branch did not reload in normal main")
    load = run["recovery_load_seconds"]
    rebuild = run["recovery_cache_rebuild_seconds"]
    initialisation = run["initialisation_seconds"]
    require(load > 0 and rebuild > 0 and load + rebuild <= initialisation,
            "recovery reader and cache timings are absent or overlap")
    endings = []
    for case, state in ((uninterrupted,p),(resumed,c)):
        final = case/"state.final.bin"
        fields, native = binary_history(final,state,state["dt"],deadline)
        require(all(native["hashes"][role] == state[digest] for role,digest in HASH_KEYS.items()),
                "native scientific and binary terminal histories differ")
        endings.append((fields,native))
    old,d0 = endings[0]
    new,d1 = endings[1]
    exact = d0["hashes"] == d1["hashes"]
    if exact:
        differences = {role:0. for role in HASH_KEYS}
        dependent = dict(K_over_Kin=0.,Gamma_over_epsilon_H=0.,cutoff_over_H=0.,
                         Ew_over_local_reference=0.,native_epsilon_H_difference=0.)
        physical = {role:0. for role in ("current","previous")}
        caches = {role:0. for role in ("advection","previous_advection")}
    else:
        from implicit_algorithm_check import coefficient_delta, dependent_delta, inventory, physical_delta, cache_duals
        differences = {role:coefficient_delta(old[role],new[role],p["sizes"],ref) for role in HASH_KEYS}
        owners, _ = inventory(p,ref)
        dependent = dependent_delta(old["current"],new["current"],p,owners,ref,deadline)
        physical = {role:physical_delta(old[role],new[role],p,ref,20,17,deadline)
                    for role in ("current","previous")}
        caches = cache_duals(old, new, p, ref, deadline,
                            cache_names=("advection", "previous_advection"))
    limit = 1e-8
    # The exact route proves all physical and cache differences simultaneously.
    norms_ok = exact
    if not exact:
        norms_ok = (dependent["gate_pass"] and all(
            max(row["RMS"],row["MAX"]) <= limit
            for data in physical.values() for row in data["norms"].values())
            and all(value["scaled_MAX"] <= limit for role in ("current","previous")
                    for value in differences[role].values())
            and all(value <= limit for row in caches.values() for value in row.values()))
    audits = []
    for audit_path, endpoint in ((uninterrupted_audit, pf), (resumed_audit, cf)):
        if audit_path is None:
            continue
        audits.append(equation_audit(audit_path,endpoint,pa["table"]))
    status = "MET"
    if not histories_ok or not norms_ok or any(row["equations"] != "MET" for row in audits):
        status = "NOT_MET"
    elif len(audits) != 2:
        status = "OPEN"
    restart_report = dict(schema="WJH27R3ColdRecoveryV1",status=status,
        state_histories_match_recovery=histories_ok,bitwise_terminal_roles=exact,limit=limit,
        endpoint_time=p["time"],endpoint_step=p["step"],on_cadence_recovery=recovery,
        positive_time_frames=[pf,cf,sf],scientific_indices=[pi,ci,si],imported_frame=identity(path),
        endpoint_histories=[d0,d1],coefficient_changes=differences,dependent_changes=dependent,physical_changes=physical,
        cache_changes=caches,recovery_initialisation_seconds=initialisation,
        recovery_load_seconds=load,recovery_cache_rebuild_seconds=rebuild,
        scheduled_split_write_seconds=read(split/"state.result.json")["recovery_seconds"],
        endpoint_equation_audits=audits)
    atomic(out,restart_report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action",required=True)
    cmd = sub.add_parser("common-time")
    for key,kind in (("coarse",Path),("fine",Path),("coarse_step",int),("fine_step",int),("out",Path)):
        cmd.add_argument(key,type=kind)
    cmd.add_argument("--seconds",required=True,type=float)
    cmd.add_argument("--coarse-audit", type=Path)
    cmd.add_argument("--fine-audit", type=Path)
    cmd = sub.add_parser("evolution")
    cmd.add_argument("case",type=Path)
    cmd.add_argument("out",type=Path)
    cmd.add_argument("--seconds",required=True,type=float)
    cmd.add_argument("--recovered-prefix",type=Path)
    cmd.add_argument("--accepted-audit",type=Path)
    cmd = sub.add_parser("cold")
    for key in ("uninterrupted","resumed","split","checkpoint","out"):
        cmd.add_argument(key,type=Path)
    cmd.add_argument("--seconds",required=True,type=float)
    cmd.add_argument("--uninterrupted-audit", type=Path)
    cmd.add_argument("--resumed-audit", type=Path)
    args=vars(parser.parse_args())
    action=args.pop("action")
    {"common-time":common_time,"cold":cold,"evolution":evolution}[action](**args)


if __name__=="__main__":
    main()
