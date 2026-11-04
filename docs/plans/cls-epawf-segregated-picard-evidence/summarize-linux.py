#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect compact sweep histories and decompose the saved terminal k residual."""

import argparse
import hashlib
import json
from pathlib import Path
import struct

import numpy as np


def read(path):
    return json.loads(path.read_text())


def pin(path):
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(part)
    return dict(path=str(path), sha256=digest.hexdigest(), bytes=path.stat().st_size)


def history_identity(path):
    report = json.loads(path.read_text(), parse_int=float)
    def coefficient_bytes(values):
        return struct.pack(f"<{len(values)}d", *values)
    return (report["steps"], report["time"], report["native_history_times"],
            coefficient_bytes(report["state"]),
            [coefficient_bytes(slot) for slot in report["native_history"]])


def k_terms(readback):
    k = 2
    space = readback["fields"][k]
    offset = readback["offsets"][k]
    fixed = {row - offset for row, _ in readback["constraints"]}
    free = np.array([i for i in range(space["size"]) if i not in fixed])
    # This fixture has affine rectangles and cubic Modified_A modes. A
    # four-point Gauss rule integrates every one-dimensional mass product.
    z, w = np.polynomial.legendre.leggauss(4)
    phi = np.array([(1-z)/2, (1+z)/2, (1-z*z)/4, z*(1-z*z)/2])
    local_mass = (phi*w) @ phi.T
    cell_mass = np.kron(local_mass, local_mass)
    mass = np.zeros((space["size"], space["size"]))
    volume = 0.
    for el in space["mesh"]:
        if el["degree"] != 3 or len(el["ids"]) != 16:
            raise ValueError("Expected the historical tensor cubic k space")
        ids = el["ids"]
        area = (el["upper"][0]-el["lower"][0])*(el["upper"][2]-el["lower"][2])
        mass[np.ix_(ids, ids)] += area/4 * np.outer(el["sign"], el["sign"]) * cell_mass
        volume += area
    names = ("time_mass", "conservative_advection", "k_diffusion",
             "inner_k_production", "inner_k_destruction",
             "outer_k_production", "outer_k_destruction", "prescribed_scalar_source")
    columns = [readback["native"]["stage_constrained"]]
    columns.extend(readback["terms"][name] for name in names)
    loads = np.array(columns).T[offset+free, :]
    dual = np.linalg.solve(mass[np.ix_(free, free)], loads)
    parameters = readback["parameters"]
    rate = parameters[5]*readback["ustar"]/parameters[1]
    scale = rate*np.sqrt(volume)
    norms = np.sqrt(np.maximum(0., np.sum(loads*dual, axis=0)))/scale
    stage_norm = readback["norms"]["native_stage"][k]
    if abs(norms[0]-stage_norm) > 1e-10*max(1., stage_norm):
        raise ValueError("Independently assembled k mass does not reproduce the native dual norm")
    # Signed contribution along the residual direction. Individual positive
    # and negative terms can cancel. The sum recovers the full stage norm.
    projections = (loads[:, 1:].T @ dual[:, 0] / (scale*scale*norms[0])
                   if norms[0] else np.zeros(len(names)))
    return {name: dict(norm=float(norm), signed_projection=float(projection))
            for name, norm, projection in zip(names, norms[1:], projections)}


def collect(root, name):
    run = root / name
    manifest = read(run / "manifest.json")
    records = read(run / "target/sweeps.json")
    result = {key: manifest[key] for key in ("source", "status", "ranks", "thread_limits",
                                              "xml_sha256", "startup", "rejected", "target",
                                              "IncNavierStokesSolver", "RANSWallResolutionProbe")
              if key in manifest}
    result["manifest"] = pin(run / "manifest.json")
    result["sweep_log"] = pin(run / "target/solver.log")
    initial = history_identity(run / "startup/case.evolution.json")
    result["deliberate_rejection_history_unchanged"] = (
        history_identity(run / "rejected/case.evolution.json") == initial)
    if manifest["target"]["steps"] == 0:
        result["failed_stage_history_unchanged"] = (
            history_identity(run / "target/case.evolution.json") == initial)
    result["sweeps"] = []
    for row in records:
        event = row["event"]
        if event in ("initial", "accepted_sweep"):
            domain = row["scalar_domain"]
            result["sweeps"].append({key: row[key] for key in
                                     ("event", "sweep", "norms", "relaxation", "trace_coefficient_defect",
                                      "seconds", "solve_counts", "ground")} |
                                    {key: domain[key] for key in
                                     ("minimum_inner_k", "minimum_outer_k", "minimum_outer_epsilon",
                                      "minimum_sampled_K", "minimum_projected_trace_epsilon",
                                      "minimum_b", "minimum_lambda", "minimum_dimensionless_margin")})
        if event in ("fixed_point", "frozen_saddle_comparison"):
            result.setdefault("verification", []).append(row)
        if event == "solve":
            result.setdefault("first_block_solves", {}).setdefault(row["block"], row)
    projection = [row for row in records if row["event"] == "ground_projection"]
    if projection:
        result["projection_attempts"] = len(projection)
        result["projection_successes"] = sum(row["solved"] for row in projection)
        result["projection_max_cuts"] = max(row["cuts"] for row in projection)
    readback = run / "target/readback.json"
    if readback.is_file():
        data = read(readback)
        result["readback"] = dict(file=pin(readback), original_norms=data["norms"]["native_stage"],
                                  term_reassembly_difference=data["norms"]["native_minus_total_constrained"],
                                  k_terms=k_terms(data))
    audit = root / "terminal-audit" / name / "audit.json"
    if audit.is_file():
        values = read(audit)[0]["audit"]
        result["terminal_qualification"] = {key: values[key] for key in
                                           ("r2_qualified", "r2_outer_range", "r2_covariance",
                                            "r4_outer_gradient_witness", "r2_scalar_domain")}
        result["terminal_audit"] = pin(audit)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    result = {name: collect(args.root, name) for name in args.cases}
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
