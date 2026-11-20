#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect digest-bound native frames, independent audits and elapsed costs on Bf."""

import argparse
import json
from pathlib import Path
import statistics

from resources import atomic, identity
from task2r_advancing import BINDING, OUTLET, identity_of_state, read, verified
from task2r_run import GIB, require


ROLES = ("current", "previous", "advection", "previous_advection")
HASH_KEYS = dict(zip(ROLES, ("state_sha256", "previous_sha256", "advection_sha256",
                            "previous_advection_sha256")))


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


def audit_plan(root, label, native, anchor_label, cap, rss_gib, address_gib, analyzer):
    from physical_oracle import load_native

    case = root / "science-r3" / label
    require(3 <= rss_gib <= 20 and 0 < address_gib <= 20
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
        bf_binding=identity(base_path), command=command, anchor=identity(anchor),
        inputs=[identity(p) for p in (native, table, analyzer, anchor)],
        resource_anchors=[bound["x15_parity"], bound["inherited"]["binding"]["operation"]],
        resource_note="Single audit process, sequential independent Gauss orders 20 and 24; "
                      "maximum RSS checked separately from inherited virtual memory limits."))


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
        evidence.update(result=identity(result), timing=dict(
            process_seconds=process, charge_seconds=receipt["result"]["seconds"],
            construction_seconds=payload["construction_seconds"],
            initialisation_seconds=payload["initialisation_seconds"],
            complete_step_seconds=complete, engine_seconds=engine, scientific_frame_seconds=frames_seconds,
            recovery_seconds=recovery, other_step_seconds=other_step, final_export_seconds=export,
            process_remainder_seconds=remainder,
            attributed_fraction=1-(other_step+engine["other"]+remainder)/process,
            equation_share=equation/process,
            warm_step_mean_seconds=statistics.mean(complete[1:]) if len(complete)>1 else None),
            observed=[{k: row[k] for k in ("step", "time", "residuals", "rebuilds", "retries",
                                           "coupled_newton", "constitutive_calls", "disjoint_seconds")}
                      for row in steps])
    failure = case / "state.failure.json"
    if failure.exists():
        evidence["failure"] = identity(failure)
    atomic(out, evidence)


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
        ("seal", (("label", str), ("out", Path))),
    ):
        mode = modes.add_parser(name)
        for key, datatype in args:
            mode.add_argument(key, type=datatype)
    args = vars(parser.parse_args())
    action = args.pop("action")
    method = dict(zip(("bind-runtime", "bind-parity", "layout", "audit-plan", "seal"),
                      (bind_runtime, bind_parity, layout_packet, audit_plan, seal_run)))[action]
    method(**args)


if __name__ == "__main__":
    main()
