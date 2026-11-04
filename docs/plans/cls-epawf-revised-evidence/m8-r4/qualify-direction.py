#!/usr/bin/env python3
"""Qualify the actual local R4 direction and gate the single pilot."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from prepare import bound, load, pin, require, save


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epoch", type=Path, required=True)
    parser.add_argument("--cases", type=Path, required=True)
    args = parser.parse_args()
    epoch, cases = args.epoch.resolve(), args.cases.resolve()
    source = epoch / "source"
    tests = source / "solvers/IncNavierStokesSolver/Tests/RANS"
    sys.path.insert(0, str(tests))
    import numpy as np
    from validate_cls_epawf_m8_r3_actions import validate_case

    receipt_path = cases / "local-run/receipt.json"
    receipt = load(receipt_path)
    require(receipt["status"] == "completed" and receipt["runtime_unchanged"], "Local native run failed")
    collection = load(bound(receipt["durable_collection"]["file"]))
    require(collection["protocol"]["local_direction_experiment"], "Wrong local run protocol")
    states = collection["completed_iterates"]
    manifest = load(bound(receipt["manifest"]))
    base = states[0]
    require(base["state"]["sha256"] == manifest["state"]["sha256"] and
        base["residual"]["sha256"] == manifest["origin"]["residual"]["sha256"],
        "Local solve changed the saved state or residual")
    assembled = collection["jacobians"][0]["event"]["matrix"]
    require(assembled["sha256"] == manifest["origin"]["matrix"]["sha256"],
            "Current local action did not reproduce the retained R3 matrix")
    bound(assembled)
    events = [json.loads(line) for line in Path(collection["progress"]["path"]).read_text().splitlines()]
    directions = [e for e in events if e["kind"] == "r4_direction"]
    require(len(directions) == 1, "No unique computed R4 direction")
    proposal = directions[0]
    step = bound(proposal["step"])
    report = cases / "action.json"
    probe = epoch / "build/RANSNonpressureActionProbe"
    runtime = load(epoch / "runtime-fingerprint.json")
    require(pin(probe)["sha256"] == runtime["runtime_artifacts"]["build/RANSNonpressureActionProbe"],
            "Current action probe has a different runtime identity")
    command = ["taskset", "-c", "0", str(probe), str(cases/"probe.xml"),
               str(cases/"direction-state.json"), str(report), str(step)]
    with (cases/"action.log").open("x") as log:
        subprocess.run(command, cwd=cases, stdout=log, stderr=subprocess.STDOUT, timeout=120, check=True)
    raw = load(report)
    require(pin(raw["matrix_file"])["sha256"] == assembled["sha256"] and
        hashlib.sha256(np.asarray(raw["state"], dtype="<f8").tobytes()).hexdigest() == base["state"]["sha256"] and
        hashlib.sha256(np.asarray(raw["history"], dtype="<f8").tobytes()).hexdigest() == manifest["history"]["sha256"] and
        hashlib.sha256(np.asarray(raw["directions"][0]["direction"], dtype="<f8").tobytes()).hexdigest() == proposal["step"]["sha256"],
        "Independent action input differs from the local direction/target")
    verification = validate_case(raw, source)
    save(cases/"direction-validation.json", dict(command=command, report=pin(report),
        probe=pin(probe), matrix=pin(raw["matrix_file"]), direction=proposal["step"], result=verification))

    newton = next(e for e in events if e["kind"] == "r4_newton_comparison")
    old_limit = load(cases/"wall-check.json")["r3_terminal_norm"]
    final = states[-1]
    actual = max(final["scaled_block_norms"])
    base_norm = max(base["scaled_block_norms"])
    meaningful = (len(states) == 2 and final["scalar_domain"]["admissible"] and
        final["scalar_domain"]["r4_wall_trace"]["admissible"] and
        actual < .99 * min(old_limit, base_norm, max(newton["actual_norms"])))
    times = collection["terminal"]["profiling"]
    trial = (2*times["scalar_domain"]["seconds"]/times["scalar_domain"]["calls"]
        + times["residual"]["seconds"]/times["residual"]["calls"]
        + times["constraints"]["seconds"]/times["constraints"]["calls"]
        + times["mass_dual_norms"]["seconds"]/times["mass_dual_norms"]["calls"])
    bound_budget = 12*(times["complete_action_assembly"]["seconds"] +
        times["r4_direction"]["seconds"] + times["factorization"]["seconds"] + 32*trial)
    practical = bound_budget+90 < 600
    gate = dict(schema=1, source_commit=receipt["source_commit"], executable_sha256=receipt["executable"]["sha256"],
        absolute_action_tolerance=1e-7, required_state_action_gates_met=verification["success"],
        active_dependencies_covered=verification["success"], useful_feasible_native_descent=meaningful and practical,
        local_evidence=pin(receipt_path), base_state=manifest["state"],
        qualification=pin(cases/"direction-validation.json"), base_norm=base_norm,
        local_norm=actual, newton_norm=max(newton["actual_norms"]), r3_terminal_norm=old_limit,
        threshold="At least 1% below base, safeguarded Newton, and R3 terminal maxima",
        estimated_12_update_32_trial_seconds=bound_budget, setup_and_audit_allowance_seconds=90,
        direction_seconds=times["r4_direction"]["seconds"],
        scope="One measured direction and actual local trial; no root or physical-step qualification")
    save(cases/"r4-gate.json", gate)
    if gate["required_state_action_gates_met"] and gate["useful_feasible_native_descent"]:
        (cases/"budget.txt").write_text(
            f"One 600s R4 pilot from the unchanged saved base. Local direction cost {gate['direction_seconds']:.6g}s. "
            f"12 updates at measured costs and up to 32 trials project {bound_budget:.6g}s, "
            "plus 90s for setup, publication and range/terminal checks.\n")
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
