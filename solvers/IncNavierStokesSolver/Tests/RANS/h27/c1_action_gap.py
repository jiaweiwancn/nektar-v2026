#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Verify and reconstruct the native R8 nonaccepted action diagnostics."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def read_audit(receipt):
    metadata = json.loads(receipt.read_text())
    if (metadata["schema"], metadata["role"], metadata["physical_history"]) != (
            "WJH27SteadyActionAuditV1", "NONACCEPTED_ACTION_VECTOR", False):
        raise ValueError("expected a nonaccepted native action receipt")
    if hashlib.sha256((metadata["signature"]+"|steady").encode()).hexdigest() != metadata["equation_identity"]:
        raise ValueError("original equation identity differs")
    record = metadata["raw_payload"]
    payload = Path(record["path"])
    data = np.memmap(payload, mode="r", dtype=np.uint8)
    if data.size != record["bytes"] or hashlib.sha256(data).hexdigest() != record["sha256"]:
        raise ValueError("raw action payload differs")
    arrays = {}
    offset = 0
    for spec in metadata["arrays"]:
        name = spec["name"]
        count = spec["rows"]*spec["columns"]
        end = offset+8*count
        raw = data[offset:end]
        if name in arrays or spec["offset_bytes"] != offset or end > data.size or (
                hashlib.sha256(raw).hexdigest() != spec["sha256"]):
            raise ValueError(f"array identity or layout differs for {name}")
        arrays[name] = raw.view("<f8").reshape(spec["rows"], spec["columns"], order="F")
        offset = end
    if offset != data.size or sum(metadata["free_counts"]) != arrays["row_scale"].shape[0]:
        raise ValueError("free-row layout or raw payload length differs")
    return metadata, arrays


def summarize(receipt, prospective_forcing):
    meta, arrays = read_audit(receipt)
    vector = lambda name: arrays[name][:, 0]
    norm = lambda x: float(np.linalg.norm(x))
    row_scale, x_scale = vector("row_scale"), vector("unknown_scale")
    base = vector("raw_base_residual")
    forcing = meta["measurements"].get("forcing") or prospective_forcing
    statistics = {row["name"]: row for row in meta["arrays"] if "physical_blocks" in row}
    output = {
        "receipt": str(receipt.resolve()),
        "receipt_sha256": hashlib.sha256(receipt.read_bytes()).hexdigest(),
        "stage": meta["stage"],
        "source": meta["source_identity"],
        "equation": meta["equation_identity"],
        "coefficient": meta["base_coefficient_sha256"],
        "factor_serial": meta["factor_serial"],
        "forcing": forcing,
        "scaled_rhs": norm(row_scale*base),
        "physical_rhs": max(statistics["raw_base_residual"]["physical_blocks"]),
        "raw": meta["raw_payload"],
        "direction_scaled_infinity": float(np.max(np.abs(vector("direction")/x_scale))),
        "rows": statistics,
    }
    if "scaled_stored_actions" in arrays:
        v, z, w, h = (arrays[key] for key in ("scaled_arnoldi_basis", "scaled_preconditioned_directions",
                                            "scaled_stored_actions", "hessenberg"))
        i = int(meta["measurements"]["basis_columns"])
        y, r0 = vector("least_squares_coefficients"), vector("scaled_restart_rhs")
        d0 = vector("scaled_seed_direction")
        stored, arnoldi = r0-w@y, r0-v@(h@y)
        fresh = -row_scale*(base+vector("raw_fresh_action"))
        reconstructed = x_scale*(d0+z@y)
        target = norm(row_scale*base)*forcing
        output["independent_reconstruction"] = dict(
            scaled_stored=norm(stored), scaled_arnoldi=norm(arnoldi), scaled_fresh=norm(fresh),
            stored_minus_arnoldi=norm(stored-arnoldi), stored_minus_fresh=norm(stored-fresh),
            stored_minus_fresh_over_scaled_target=norm(stored-fresh)/target,
            basis_orthogonality_frobenius=norm(v[:, :i].T@v[:, :i]-np.eye(i)),
            arnoldi_relation_frobenius=norm(w-v@h),
            least_squares_reconstruction=abs(norm(arnoldi)-meta["measurements"]["projected_scaled"]),
            direction_reconstruction=norm(reconstructed-vector("direction")),
            scaled_direction_amplification=meta["measurements"]["scaled_direction_sum"]/norm(reconstructed/x_scale),
        )
    if "raw_action_0" in arrays:
        output["fixed_direction"] = dict(
            scaled_change=[norm(row_scale*(vector(f"raw_action_{j}")-vector(f"raw_action_{j-1}")))
                           for j in (1, 2, 3)],
            repeat_action_equal=np.array_equal(vector("raw_action_4"), vector("raw_action_1")),
            repeat_residual_difference_equal=np.array_equal(vector("raw_difference_4"), vector("raw_difference_1")),
        )
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipts", nargs="+", type=Path)
    parser.add_argument("--forcing", type=float, default=.03, help="target for the standalone smoke test")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if not 0 < args.forcing < 1:
        parser.error("forcing must be a fraction")
    reports = [summarize(path, args.forcing) for path in args.receipts]
    with args.out.open("x") as stream:
        json.dump(reports,stream,indent=2,allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
