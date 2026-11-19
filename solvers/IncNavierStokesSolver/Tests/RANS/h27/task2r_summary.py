#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Extract compact verified Task 2R physical and spatial evidence."""

import argparse
from collections import Counter
import json
from pathlib import Path

from resources import atomic, identity
from task2r_reference import D0, H


FIELDS = ("U", "W", "k", "R", "Pi", "traction", "G", "epsilon")
L0 = .04
TIMING_KEYS = ("construction_seconds", "initialisation_seconds", "complete_advancing_seconds",
               "final_export_seconds")


def read(path):
    return json.loads(path.read_text())


def extrema(norms, fields):
    rows = {}
    for field in fields:
        selected = {k: v for k, v in norms.items() if k.endswith("/"+field)}
        if selected:
            r = max(selected, key=lambda k: selected[k]["RMS"])
            m = max(selected, key=lambda k: selected[k]["MAX"])
            rows[field] = dict(RMS=selected[r]["RMS"], RMS_group=r, MAX=selected[m]["MAX"], MAX_group=m)
    return rows


def steady_timing(data):
    # These callbacks have separate RAII timers. Line-search time spans
    # callbacks already counted here and must not be charged a second time.
    components = dict(data["callback_seconds"], nonlinear_validation=data["nonlinear_validation_seconds"])
    components["other_inside_krylov"] = data["total_seconds"] - sum(components.values())
    return components


def cases(root):
    result = {}
    for location in sorted((root/"science").glob("r3-*")):
        inputs = location/"inputs.json"
        native = location/"state.result.json"
        independent = location/"independent.json"
        steady = location/"state.steady.json"
        if not inputs.is_file():
            continue
        spec = read(inputs)
        current = dict(inputs=identity(inputs), mesh=spec["geometry"], session=spec["session"],
                       cells=spec["cells"], order=spec["modes"], quadrature=spec["quadrature"],
                       convergence="NO_STEADY_RESULT")
        if steady.is_file():
            response = read(steady)
            current.update(native_steady=identity(steady),
                           convergence=response["status"], reason=response["reason"], norm=response["true_norm"],
                           native_krylov_seconds=response["total_seconds"], krylov_attribution=steady_timing(response),
                           accepted_newton_steps=response["accepted_newton_steps"],
                           jacobian_actions=response["jacobian_actions"],
                           physical_time_advanced=response["physical_time_advanced"])
        if native.is_file():
            record = read(native)
            state = record["state"]
            current.update(native_state=identity(native), source_identity=state["source_identity"],
                           reference_identity=state["reference_identity"],
                           process_timers={k:record[k] for k in TIMING_KEYS},
                           native_guards={k:state[k] for k in ("minimum_k", "minimum_epsilon",
                                                              "minimum_covariance_over_k", "minimum_principal_over_mu")})
        if independent.is_file():
            doc = read(independent)
            core = doc["core"]
            current.update(independent=identity(independent), current_equations=doc["equations"],
                           physical=doc["physical"], current_residuals=core["residuals"],
                           accounting_identity_error=core["accounting_identity_max_error"],
                           domain_norms=core["physical_envelope"], stations=doc["observations"],
                           station_errors=extrema(doc["profiles"]["norms"], ("U","W","k","epsilon")),
                           station_status=doc["profiles"]["status"],
                           historical_inner=extrema(doc["historical_inner_diagnostic"], ("U","k")),
                           owner_fluxes={k:v for k,v in doc["interfaces"]["norms"].items() if k.startswith("H/")},
                           ground_vs_new_reference=core["ground_shear_diagnostic"],
                           normalized_budget_defects={k:{
                               f:core["budgets"][k][f]/core["budgets"][k]["fixed_scale"]
                               for f in ("mixed_physical_balance","all_owner_gradient_balance")}
                               for k in ("k","epsilon","mechanical")})
        result[location.name] = current
    return result


def comparisons(root):
    result = {}
    for location in sorted((root/"science").glob("r3-*")):
        for candidate in sorted(location.glob("*comparison.json")):
            data = read(candidate)
            records = data["norms"]
            not_met = [k for k,v in records.items() if v["status"] != "MET"]
            result[location.name+"/"+candidate.name] = dict(status=data["status"], witness=identity(candidate),
                native=data["native"], paired_native=data["paired_native"], control=data["control"],
                extrema=extrema(records,FIELDS), norm_count=len(records), not_met=not_met,
                failed_norms={name:records[name] for name in not_met},
                common_owner_subcells=data["common_owner_subcells"])
    return result


def summarize(root, out):
    acc = read(root/"accounting.json")
    if [v for v in acc["entries"].values() if v["status"]=="RESERVED" and v["phase"]!="evidence"]:
        raise ValueError("scientific operation still reserved")
    projection = {}
    for name in ("r2-graded","r3-axis-redistribution"):
        p = root/"science"/name/"projection.json"
        if p.exists():
            projection[name] = dict(witness=identity(p), **read(p))
    operations = {}
    seconds = Counter()
    for label, event in sorted(acc["entries"].items()):
        if not label.startswith("r3-") or event["status"] == "RESERVED":
            continue
        receipt = root/"operations"/label/"receipt.json"
        operations[label] = dict(
            **{key:event[key] for key in ("phase", "status", "seconds", "peak_tree_rss_bytes",
                                          "peak_host_workers", "reason")},
            receipt=identity(receipt))
        seconds[event["phase"]] += event["seconds"]
    variants = dict(table=identity(root/"science/r1-qualified/tables-exact-inner-v1.json"),
        compiled_artifact=identity(root/"runtime-r2-v1/IncNavierStokesSolver"),
        calibration_policy="ExactInnerFiniteNuV1", inner_transport="ScalarCdV2",
        source_identity="061c48fe6805b1b1bb66a8effffc3bfc28c5f86777c5e669f52afd9f8c0d861f")
    atomic(out, dict(schema="WJH27Task2RPhysicalQualificationEvidenceV2", geometry=dict(length_unit_m=L0,
         matching_wall_radius_code=H, matching_height_m=L0*H,
         matching_plane_elevation_above_original_floor_m=L0*(D0+H)),
         identities=variants, scientific_runs=cases(root), comparisons=comparisons(root), projection=projection,
         localization=[identity(p) for p in sorted((root/"science").glob("r3-*/*witness.json"))]
         + [identity(p) for p in sorted((root/"science").glob("r3-*/localization.json"))],
         proposals=[identity(p) for p in sorted((root/"evidence").glob("r3-*-admission.json"))],
         operations=operations, settled_seconds_before_publication_by_phase=seconds))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", type=Path)
    p.add_argument("out", type=Path)
    a = p.parse_args()
    summarize(a.root, a.out)


if __name__=="__main__":
    main()
