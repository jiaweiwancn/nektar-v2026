#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect signed bounded-direction results and verify accepted checkpoint lineage."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from c1_action_gap import read_audit
from c1_bounded_trials import SETTINGS, load, packed
from resources import identity


def assemble(case, seed, validate_vectors):
    origin = load(seed)
    previous = origin
    initial = max(origin["original_residuals"])
    start = int((case/"start.epoch").read_text())
    basis = {}
    trials = []
    for path in case.glob("state.action-*.json"):
        m = load(path)
        if m.get("stage") == "bounded-current-basis":
            assert m["correction_mass"] == 0 and m["fixed_k_tangent"]
            basis[m["newton"]] = (identity(path), m)
        if m.get("stage") == "bounded-direction":
            trials.append((path, m))
    trials.sort(key=lambda item: (item[1]["newton"],item[1]["proposal"]))
    rows, acceptances = [], []
    for path, trial in trials:
        index = trial["newton"]
        current = basis[index][1]
        assert (trial["base_coefficient_sha256"], trial["bounded_policy_sha256"],
                trial["equation_identity"], trial["free_id_sha256"]) == (
                current["base_coefficient_sha256"], current["bounded_policy_sha256"],
                current["equation_identity"], current["free_id_sha256"])
        assert trial["correction_mass"] == 0 and trial["fixed_k_tangent"]
        row = {"base":index, "proposal":trial["proposal"], "receipt":identity(path),
               "payload":trial["raw_payload"], "measurements":trial["measurements"],
               "notes":trial["notes"],"base_coefficients":trial["base_coefficient_sha256"],
               "blocks":{part["name"]:part["physical_blocks"] for part in trial["arrays"]
                         if "physical_blocks" in part},
               "factor":trial["factor"], "stencils":trial["stencils"]}
        rows.append(row)
        if trial["notes"]["decision"] == "BOUNDED_ACCEPTED":
            acceptances.append((path,row))
    all_checkpoints = list(sorted(case.glob("state.accepted-*.json")))
    assert len(acceptances) >= len(all_checkpoints)
    checkpoints, products = [], 0
    for path, (trial_path,row) in zip(all_checkpoints,acceptances):
        packet = load(path)
        controller = packet["merit_controller"]
        radius = controller["state"]["radius"]
        m = row["measurements"]
        assert row["base_coefficients"] == previous["state"]["state_sha256"]
        assert all(packet["state"][key] == origin["state"][key] for key in SETTINGS)
        assert packet["iterate_custody"] == "LAST_ACCEPTED" and packet["attempted_time"] == 0
        assert np.isclose(max(packet["original_residuals"]), m["physical_actual"], rtol=1e-12, atol=1e-16)
        assert radius["used"] == m["radius"] and radius["scaled_displacement"] == m["scaled_displacement"]
        if validate_vectors:
            previous_free, _ = packed(previous, previous["state"]["state_sha256"])
            current_free, _ = packed(packet, packet["state"]["state_sha256"])
            _, vectors = read_audit(trial_path)
            dx = np.asarray(controller["accepted_displacement"], dtype="<f8")
            assert controller["accepted_displacement_sha256"] == hashlib.sha256(dx.tobytes()).hexdigest()
            assert (np.array_equal(vectors["base"][:,0],previous_free) and
                    np.array_equal(vectors["trial_free"][:,0],current_free) and
                    np.array_equal(current_free-previous_free,dx)), "accepted displacement differs"
            assert (vectors["raw_remainder"] ==
                    vectors["raw_actual_residual"]-vectors["raw_model_residual"]).all()
        products += int(m["current_products_at_base"])
        elapsed = path.stat().st_mtime - start
        final = max(packet["original_residuals"])
        checkpoints.append({"checkpoint":identity(path), "coefficients":packet["state"]["state_sha256"],
            "displacement_sha256":controller["accepted_displacement_sha256"],
            "policy_sha256":controller["policy_sha256"], "radius":radius, "margins": {
                key:packet[key] for key in ("minimum_covariance_over_k","original_principal_over_mu",
                                            "physical_whole_cell_scalar_lower_bounds")},
            "trace_K_over_Kin":packet["state"]["interface_K_over_Kin"],
            "residuals":packet["original_residuals"], "cumulative_products":products,
            "elapsed_to_checkpoint":elapsed, "actual_decrease_ratio":controller["state"]["actual_to_predicted_decrease"],
            "physical_reduction_per_product":(initial-final)/max(products,1),
            "physical_reduction_per_process_second":(initial-final)/max(elapsed,1.)})
        previous = packet
    outcome = case/"state.steady.json"
    solve = load(outcome) if outcome.exists() else None
    return {"seed":identity(seed), "cooperative_terminal":None if solve is None else {
                "receipt":identity(outcome),"status":solve["status"],"reason":solve["reason"],
                "seconds":solve["total_seconds"],"physical_time_advanced":solve["physical_time_advanced"]},
            "basis":{i:{"receipt":rec,"payload":row["raw_payload"],"measurements":row["measurements"],
                         "base_coefficients":row["base_coefficient_sha256"],"factor":row["factor"]}
                     for i,(rec,row) in basis.items()},
            "trials":rows, "accepted":checkpoints,"initial_max":initial,"best_max":max(previous["original_residuals"])}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    parser.add_argument("seed", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--skip-raw-vectors", action="store_true")
    args = parser.parse_args()
    result = assemble(args.case.resolve(),args.seed.resolve(),not args.skip_raw_vectors)
    with args.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
