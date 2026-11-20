#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect digest-bound native frames, independent audits and elapsed costs on Bf."""

import argparse
import json
import math
from pathlib import Path
import statistics

from resources import atomic, identity
from task2r_advancing import BINDING, OUTLET, identity_of_state, read, verified
from task2r_run import GIB, require


ROLES = ("current", "previous", "advection", "previous_advection")
HASH_KEYS = dict(zip(ROLES, ("state_sha256", "previous_sha256", "advection_sha256",
                            "previous_advection_sha256")))


def active_step_bounds(steps, setup, export, operation, initial_step):
    """Locate each native scope up to the full unassigned operation time."""
    gap = operation["seconds"] - setup - sum(steps) - export
    require(gap >= 0, "native timing scopes exceed the charged operation")
    pauses = [(row["start_seconds"],row.get("end_seconds",row.get("settled_seconds",operation["seconds"])))
              for row in operation["worker_pauses"]]
    def overlap(left,right):
        return sum(max(0.,min(right,b)-max(left,a)) for a,b in pauses)

    out = []
    start = setup
    for step, duration in enumerate(steps,initial_step+1):
        end = start+duration
        out.append(dict(step=step,complete_seconds=duration,
            active_lower_seconds=max(0.,duration-overlap(start,end+gap)),
            active_upper_seconds=duration-overlap(start+gap,end)))
        start = end
    return out


def bind_runtime(root, executable, source_root, out):
    from task2_outlet_baseline import source_identity

    base = read(root / "evidence" / BINDING)["inherited"]["binding"]
    source, inventory = source_identity(source_root)
    require(source.encode() in executable.read_bytes(), "R3 solver lacks the archived source identity")
    atomic(out, dict(schema="WJH27R3RebuiltRuntimeV1", application=identity(executable),
                     source=source, source_inventory=inventory, state=base["state"],
                     seed_source=base["source"]))


def bind_parity(root, label, out):
    base = read(root / "evidence" / BINDING)["inherited"]["binding"]
    case = root / "science-r3" / label
    plan = read(case / "admission.json")
    result = case / "state.result.json"
    attempt = root / "operations" / label / "receipt.json"
    data = identity_of_state(result)
    original = identity_of_state(verified(base["state"]))
    require(read(attempt)["result"]["status"] == "PASS"
            and plan["phase"] == "r3-parity" and read(result)["steady_solver"]["converged"]
            and data["source_identity"] == plan["source"] and data["step"] == data["time"] == 0
            and data["initial_field_transfer"]["packet_sha256"] == base["state"]["sha256"]
            and not data["initial_field_transfer"]["projected"] and max(data["residuals"]) <= 1e-10,
            "same-mesh typed preparation did not converge or changed its recipe")
    exact = ("state_sha256", "reference_identity", "outlet_identity", "inner_k_transport",
             "elements", "sizes", "residuals", "quadratureX", "quadratureR")
    require(all(data[key] == original[key] for key in exact),
            "rebuilt application differs from the original Bf science at step zero")
    atomic(out, dict(schema="WJH27R3SameMeshParityV1", status="MET",
                     application=plan["application"], source=plan["source"], state=identity(result),
                     original=base["state"], coefficient_sha256=data["state_sha256"],
                     exact_fields=exact, admission=identity(case / "admission.json"),
                     receipt=identity(attempt)))


def layout_packet(root, parity_path, out):
    from perturbation_certificate import BF_SHAPE, independent_layout
    from task2_outlet_baseline import mesh_axes

    parity = read(parity_path)
    require(parity["status"] == "MET", "perturbation requires verified rebuilt parity")
    data = identity_of_state(verified(parity["state"]))
    geometry = read(root / "evidence" / BINDING)["inherited"]["binding"]["mesh"]
    x, z, _ = mesh_axes(verified(geometry))
    d0, H = data["geometry"]["d0"], data["geometry"]["H"]
    r = [v-d0 for v in z]
    r = [next((edge for edge in (0., H) if abs(value-edge) < 1e-14), value) for value in r]
    shape = independent_layout(x, r, H, BF_SHAPE)
    atomic(out, dict(schema="WJH27PerturbationNativeProbeV1", mode="layout",
        provenance="Native signed geometry from the rebuilt registered Bf steady replay",
        tensor_inventory_verified=True, H=H, domain=[x[0], x[-1], r[0], r[-1]],
        source_identity=data["source_identity"], reference_identity=data["reference_identity"],
        base_state_sha256=data["state_sha256"], mesh_xml_sha256=geometry["sha256"],
        x_edges=x, r_edges=r, elements=[dict(bounds=e["bounds"]) for e in data["elements"]],
        layouts={BF_SHAPE: shape}, same_mesh_parity=identity(parity_path)))


def audit_plan(root, label, native, anchor_label, cap, rss_gib, address_gib, analyzer,
               maximum_pause=None):
    from physical_oracle import load_native

    case = root / "science-r3" / label
    require(3 <= rss_gib <= 20 and 0 < address_gib <= 20
            and (maximum_pause is None or 0 < maximum_pause <= cap)
            and analyzer.name == "physical_transient_oracle.py",
            "invalid independent audit resources")
    anchor = root / "operations" / anchor_label / "receipt.json"
    prior = read(anchor)
    require(prior["request"]["phase"] in ("r3-advance", "r3-restart")
            and native.parent == root / "science-r3" / anchor_label,
            "independent audit needs a recorded Bf run")
    packet, _ = load_native(native)
    require(packet["step"] > 0 and packet["outlet_policy"] == OUTLET and "input_history" in packet,
            "audit target lacks genuine accepted histories")
    base_path = root / "evidence" / BINDING
    bound = read(base_path)
    table = verified(bound["inherited"]["table"])
    case.mkdir()
    command = ["python3", "-B", str(analyzer), "--native", str(native), "--table", str(table),
               "--out", str(case / "audit.json"), "--seconds", str(cap-15),
               "--address-space-gib", str(address_gib)]
    atomic(case / "admission.json", dict(schema="WJH27R3AdmissionV1", label=label,
        phase="r3-audit", role="accepted_time_discrete", seconds=cap, peak_new_bytes=192*1024**2,
        rss_limit_bytes=rss_gib*GIB, address_space_gib=address_gib, cpu=bound["affinity"]["cpu"],
        maximum_paused_seconds=maximum_pause,
        bf_binding=identity(base_path), command=command, anchor=identity(anchor),
        inputs=[identity(p) for p in (native, table, analyzer, anchor)],
        resource_anchors=[bound["x15_parity"], bound["inherited"]["binding"]["operation"]],
        resource_note="Single audit process, sequential independent Gauss orders 20 and 24; "
                      "maximum RSS checked separately from inherited virtual memory limits."))


def perturbation_plan(root, label, initial_label, parity, certificate, cap, rss_gib, analyzer,
                      maximum_pause=None):
    artifact = root / "science-r3" / label
    start = root / "science-r3" / initial_label
    plan = read(start / "admission.json")
    receipt = root / "operations" / initial_label / "receipt.json"
    basis = read(parity)
    require(plan["phase"] == "r3-initial" and plan["source"] == basis["source"]
            and read(receipt)["result"]["status"] == "PASS"
            and plan["initial_fields"] == basis["state"] and plan["certificate"] == identity(certificate)
            and (maximum_pause is None or 0 < maximum_pause <= cap)
            and 0 < rss_gib <= 20 and analyzer.name == "task2r_perturbation_audit.py",
            "represented perturbation audit requires the named step-zero initialized data")
    anchor = verified(plan["bf_binding"])
    bound = read(anchor)
    base, native = verified(basis["state"]), start / "state.result.json"
    table = verified(bound["inherited"]["table"])
    file_inputs = (base, native, table, certificate, receipt, parity, analyzer)
    command = ["python3", "-B", str(analyzer), "--base", str(base), "--native", str(native),
               "--table", str(table), "--certificate", str(certificate),
               "--out", str(artifact / "audit.json"), "--seconds", str(cap-15)]
    artifact.mkdir()
    atomic(artifact / "admission.json", dict(schema="WJH27R3AdmissionV1", label=label,
        phase="r3-audit", role="represented_perturbation", seconds=cap, peak_new_bytes=48*1024**2,
        maximum_paused_seconds=maximum_pause,
        cpu=bound["affinity"]["cpu"], bf_binding=identity(anchor), rss_limit_bytes=rss_gib*GIB,
        command=command, inputs=[identity(p) for p in file_inputs], anchor=identity(receipt),
        resource_note="Exact rational signed basis and independent dependent-trace and principal checks."))


def seal_run(root, label, out):
    import physical_transient_oracle as transient

    case = root / "science-r3" / label
    plan = read(case / "admission.json")
    receipt_path = root / "operations" / label / "receipt.json"
    receipt = read(receipt_path)
    index_path = case / "state.science-index.json"
    evidence = dict(schema="WJH27R3NativeAdvanceEvidenceV1", label=label, role=plan["role"],
                    admission=identity(case / "admission.json"), operation=identity(receipt_path),
                    status=receipt["result"]["status"], exit=receipt["result"]["exit"],
                    reason=receipt["result"]["reason"])
    if not index_path.exists():
        atomic(out, evidence)
        return
    index = read(index_path)
    steps_log = case / "state.steps.jsonl"
    steps = [json.loads(line) for line in steps_log.read_text().splitlines()]
    dt = plan["timestep"]
    initial_step = plan.get("initial_step", 0)
    require(index["source_identity"] == plan["source"] and index["dt"] == dt
            and index["coefficients_included"] and index["initial_step"] == initial_step
            and index["initial_time"] == initial_step*dt
            and len(steps) == index["final_step"]-initial_step, "scientific clock inventory differs")
    by_step = {row["step"]: row for row in steps}
    for i, row in enumerate(steps, initial_step+1):
        require(row["step"] == i and row["time"] == i*dt
                and row["signature"] == index["signature"]
                and row["outlet_policy"] == OUTLET and row["implicit_solver"] == "CoupledJFNKV1"
                and row["inner_k_transport"] == "ScalarCdV2" and max(row["residuals"]) <= 1e-10
                and min(row[k] for k in ("minimum_k", "minimum_epsilon", "minimum_covariance_over_k",
                                         "minimum_principal_over_mu")) > 0,
                f"accepted native step {i} missed the clock, model, tolerance or guards")
    frames = {}
    for frame in index["frames"]:
        n = frame["step"]
        path = verified({k: frame[k] for k in ("path", "bytes", "sha256")})
        packet = identity_of_state(path)
        require(path == case / f"state.science-{n}.json" and n not in frames
                and packet["step"] == n and packet["time"] == frame["time"] == dt*n
                and packet["signature"] == index["signature"], "invalid scientific frame lineage")
        if n > initial_step:
            require(frame["origin"] == "current_process_accepted_step", "invalid accepted frame origin")
            require(all(packet[k] == by_step[n][k] for k in HASH_KEYS.values()), "frame and step log disagree")
            h = packet["input_history"]
            require(h["step"] == n-1 and h["time"] == dt*(n-1), "frame has the wrong input clock")
            digest = {HASH_KEYS[k]: transient.state_hash(h[k]) for k in ROLES}
            previous = by_step.get(n-1, frames.get(n-1, {}).get("hashes"))
            if previous:
                require(all(digest[k] == previous[k] for k in HASH_KEYS.values()), "explicit history roles differ")
        else:
            origin = "imported_accepted_history" if initial_step else "fresh_initial_observation"
            require(n == initial_step and frame["origin"] == origin, "invalid initial observation lineage")
            if not initial_step:
                require("input_history" not in packet and packet["previous_sha256"] == packet["state_sha256"]
                        and packet["previous_advection_sha256"] == packet["advection_sha256"],
                        "fresh initial histories were not duplicated")
        frames[n] = dict(frame, hashes={k: packet[k] for k in HASH_KEYS.values()})
    evidence.update(index=identity(index_path), step_log=identity(steps_log), frames=frames,
                    accepted_steps=len(steps), final_time=index["final_time"],
                    recoveries=[identity(p) for p in sorted(case.glob("state.*.bin"))])
    if receipt["result"]["status"] == "PASS":
        result = case / "state.result.json"
        payload = read(result)
        require(index["run_complete"] and len(steps) == plan["requested_steps"]
                and payload["normal_main"] and payload["steady_mode"] is False
                and payload["state"]["time"] == index["final_time"], "positive-time output is incomplete")
        complete = payload["complete_step_seconds"]
        engine = {k: sum(row["disjoint_seconds"][k] for row in steps)
                  for k in steps[-1]["disjoint_seconds"]}
        setup = sum(payload[k] for k in ("construction_seconds", "initialisation_seconds"))
        engine_total = sum(engine.values())
        frames_seconds = payload["scientific_frame_seconds"]
        recovery = sum(payload["recovery_seconds"])
        export = payload["final_export_seconds"]
        process = receipt["result"].get("child_elapsed_seconds") or receipt["result"]["seconds"]
        other_step = sum(complete)-engine_total-frames_seconds-recovery
        remainder = process-setup-sum(complete)-export
        require(other_step >= -1e-8 and remainder >= 0 and
                abs(engine_total-sum(row["seconds"] for row in steps)) <= 1e-6,
                "whole-process timing scopes overlap")
        equation = sum(engine[k] for k in ("constitutive", "wall", "assembly", "flow", "scalar",
                                        "forcing", "residual"))
        detail = {}
        if "timing_detail_seconds" in steps[0]:
            parts = {key: sum(row["timing_detail_seconds"][key] for row in steps)
                     for key in steps[0]["timing_detail_seconds"]}
            momentum = parts["momentum_advection"]
            transport = parts["scalar_advection"]
            coupling = parts["scalar_momentum_coupling"]
            movement = parts["history_movement"]
            require(0 <= momentum+transport <= engine["assembly"] and 0 <= coupling <= engine["flow"]
                    and 0 <= movement <= engine["other"], "step subscopes overlap")
            detail = dict(engine, **parts)
            detail.update(assembly=engine["assembly"]-momentum-transport,
                          flow=engine["flow"]-coupling, other=engine["other"]-movement)
        bounds = active_step_bounds(complete,setup,export,receipt["result"],initial_step)
        # Exclude BE startup, the BDF2 coefficient change and the first
        # local step after recovery.
        warm = bounds[1 if initial_step else 2:]
        evidence.update(result=identity(result), timing=dict(
            process_seconds=process, charge_seconds=receipt["result"]["seconds"],
            construction_seconds=payload["construction_seconds"],
            initialisation_seconds=payload["initialisation_seconds"],
            complete_step_seconds=complete, engine_seconds=engine, scientific_frame_seconds=frames_seconds,
            detailed_engine_seconds=detail if detail else None,
            active_step_bounds=bounds,
            recovery_seconds=recovery, other_step_seconds=other_step, final_export_seconds=export,
            process_remainder_seconds=remainder,
            attributed_fraction=1-(other_step+(detail or engine)["other"]+remainder)/process,
            equation_share=equation/process,
            warm_step_mean_seconds=statistics.mean(row["complete_seconds"] for row in warm) if warm else None,
            warm_active_mean_bounds_seconds=[
                statistics.mean(row[key] for row in warm)
                for key in ("active_lower_seconds", "active_upper_seconds")] if warm else None),
            observed=[{k: row[k] for k in ("step", "time", "residuals", "rebuilds", "retries",
                                           "coupled_newton", "constitutive_calls", "disjoint_seconds")}
                      for row in steps])
    failure = case / "state.failure.json"
    if failure.exists():
        evidence["failure"] = identity(failure)
    atomic(out, evidence)


def time_matrix(root, startup, perturbation, parity_path, represented, preconditioner, out,
                recovered_prefix=None, coarse_dt=None, selection_evidence=None):
    from task2r_advancing import STEP_PRECONDITIONERS

    parity = read(parity_path)
    proof = read(represented)
    require(preconditioner in STEP_PRECONDITIONERS and parity["status"] == "MET"
            and proof["status"] == "MET" and proof["base"] == parity["state"]
            and all(path.parent == root / "evidence" for path in (startup, perturbation, parity_path, out)),
            "matrix needs a qualified Bf application and represented initial state")
    prefixes = [read(path) for path in (startup, perturbation)]
    levels = []
    for index, (prefix, role) in enumerate(zip(prefixes, ("startup", "perturbed"))):
        require(prefix["status"] == "PASS", "matrix bracket needs both completed evolving prefixes")
        plan = read(verified(prefix["admission"]))
        measured = prefix["timing"]
        step_cost = measured["complete_step_seconds"]
        hashes = {frame["hashes"]["state_sha256"] for frame in prefix["frames"].values()}
        accepted_steps = prefix["accepted_steps"]
        cost_source = "Observed complete native steps."
        if role == "perturbed" and recovered_prefix:
            before = read(recovered_prefix)
            checkpoint = before["recovery"]
            accepted = next(frame for frame in prefix["frames"].values() if frame["step"] == checkpoint["step"])
            require(prefix["role"] == "restart"
                    and plan["restart"] == {k:checkpoint[k] for k in ("path", "sha256", "bytes")}
                    and checkpoint["dt"] == plan["timestep"]
                    and accepted["origin"] == "imported_accepted_history"
                    and before["state_hashes"] == accepted["hashes"], "recovered prefix cannot supply a bracket")
            first_frame = before["initial"]
            early_steps = []
            for step in range(checkpoint["step"], 0, -1):
                require(before["schema"] == "WJH27R3RecoveredPrefixAnchorV1" and before["status"] == "MET"
                        and before["recovery"]["step"] == step and before["initial"] == first_frame,
                        "recovered accepted steps are discontinuous")
                verified(before["accepted"])
                rows = map(json.loads, verified(before["step_log"]).read_text().splitlines())
                earlier, = (row for row in rows if row["step"] == step)
                require(all(earlier[key] == before["state_hashes"][key] for key in HASH_KEYS.values()),
                        "recovered accepted history and step log disagree")
                hashes.add(earlier["state_sha256"])
                early_steps.append(earlier["seconds"])
                if step > 1:
                    before = read(verified(before["preceding_anchor"]))
            hashes.add(identity_of_state(verified(first_frame))["state_sha256"])
            output = (measured["scientific_frame_seconds"] + measured["recovery_seconds"])/accepted_steps
            step_cost = [*(value+output for value in reversed(early_steps)), *step_cost]
            accepted_steps += len(early_steps)
            cost_source = "Earlier measured accepted engines and an I/O allowance from the completed resumed run."
        else:
            require(prefix["role"] == role, "matrix prefix role differs")
        require(accepted_steps >= 3 and len(hashes) >= 4, "matrix prefix has insufficient evolution")
        trial = plan["timestep"]
        if index == 0:
            prefix_dt = trial
            dt = trial if coarse_dt is None else coarse_dt
            require(dt in (trial, trial/2) and (dt == trial or selection_evidence is not None),
                    "a reduced bracket needs a measured prefix result")
        require(trial == prefix_dt, "the evolving prefixes have different bracket anchors")
        operation = read(verified(prefix["operation"]))
        require(operation["result"]["status"] == "PASS", "prefix operation did not pass")
        setup = measured["construction_seconds"] + measured["initialisation_seconds"] + measured["final_export_seconds"]
        pause_fraction = operation["result"]["worker_paused_seconds"] / operation["result"]["seconds"]
        starts = operation["request"]["storage_before"]["r3_science_bytes"]
        peak = operation["result"]["peak_r3_science_bytes"] - starts
        peak_mib = max(550, math.ceil(peak / 1024 ** 2) + 150)
        for denominator in (1, 2, 4):
            n = 4 * denominator
            duration = setup + sum(step_cost[:2]) + (n-2) * max(step_cost[2:])
            cap = math.ceil(1.35 * duration + 60)
            require(cap <= 7200, "measured matrix admission exceeds the per-process advancement limit")
            levels.append(dict(role=role, label=f"r3-time-{role}-{denominator}", timestep=dt/denominator,
                steps=n, common_steps=[2*denominator, n], science_stride=2*denominator,
                recovery_stride=2*denominator, estimated_elapsed_seconds=duration, cost_source=cost_source,
                operation_cap_seconds=cap,
                maximum_pause_seconds=min(cap-1, max(180, min(cap//2, math.ceil(1.2*pause_fraction*cap)))),
                rss_gib=16, peak_mib=peak_mib, native_step_cap_seconds=plan["native_step_cap_seconds"]))
    require(len(levels) == 6 and prefixes[0]["label"] != prefixes[1]["label"], "matrix inventory differs")
    atomic(out, dict(schema="WJH27R3ProspectiveTimeMatrixV1", prefixes=[identity(startup), identity(perturbation)],
        parity=identity(parity_path), represented_perturbation=identity(represented),
        recovered_prefix=identity(recovered_prefix) if recovered_prefix else None,
        prefix_timestep=prefix_dt, bracket_selection_evidence=identity(selection_evidence) if selection_evidence else None,
        preconditioner=preconditioner, timesteps=[dt, dt/2, dt/4], common_times=[2*dt, 4*dt], levels=levels,
        finer_witness="Require D/8 at both common times if D/4 is the first viable timestep.",
        evaluation="Compare the two transients separately at both times across the full owner inventory."))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    modes = parser.add_subparsers(dest="action", required=True)
    for name, args in (
        ("bind-runtime", (("executable", Path), ("source_root", Path), ("out", Path))),
        ("bind-parity", (("label", str), ("out", Path))),
        ("layout", (("parity_path", Path), ("out", Path))),
        ("audit-plan", (("label", str), ("native", Path), ("anchor_label", str), ("cap", int),
                        ("rss_gib", int), ("address_gib", int), ("analyzer", Path))),
        ("perturbation-plan", (("label", str), ("initial_label", str), ("parity", Path),
                               ("certificate", Path), ("cap", int), ("rss_gib", int), ("analyzer", Path))),
        ("time-matrix", (("startup", Path), ("perturbation", Path), ("parity_path", Path),
                         ("represented", Path), ("preconditioner", str), ("out", Path))),
        ("seal", (("label", str), ("out", Path))),
    ):
        mode = modes.add_parser(name)
        for key, datatype in args:
            mode.add_argument(key, type=datatype)
        if name in ("audit-plan", "perturbation-plan"):
            mode.add_argument("--maximum-pause", type=int)
        if name == "time-matrix":
            mode.add_argument("--recovered-prefix", type=Path)
            mode.add_argument("--coarse-dt", type=float)
            mode.add_argument("--selection-evidence", type=Path)
    args = vars(parser.parse_args())
    action = args.pop("action")
    method = {"bind-runtime": bind_runtime, "bind-parity": bind_parity, "layout": layout_packet,
              "audit-plan": audit_plan, "perturbation-plan": perturbation_plan,
              "time-matrix": time_matrix, "seal": seal_run}[action]
    method(**args)


if __name__ == "__main__":
    main()
