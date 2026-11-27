#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Solve a signed current-basis model for the opt-in bounded B3 policy."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from c1_action_gap import read_audit

HELPER_SHA = "5818b2c71dea79b11d4c5eb3f431a98682529631f5d4556cd82e2057f2749941"
HELPERS = Path(__file__).resolve().parents[5]/"docs/plans/wj-h27-linux/b3-frozen-size"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def choose(receipt, radius, output):
    if sha(HELPERS/"frozen_size.py") != HELPER_SHA:
        raise ValueError("pinned joint-model solver differs")
    sys.path.insert(0, str(HELPERS))
    from frozen_size import BLOCKS, JointModel, array_spec

    metadata, arrays = read_audit(receipt)
    basis_sha = sha(receipt)
    if (metadata["stage"] != "bounded-current-basis" or metadata["correction_mass"] != 0 or
            not metadata["fixed_k_tangent"] or len(metadata["bounded_policy_sha256"]) != 64):
        raise ValueError("expected a signed unshifted hybrid current basis")
    q, response = (arrays[name] for name in ("scaled_correction_basis", "raw_residual_and_actions"))
    sr = arrays["row_scale"][:,0]
    if (q.shape[1] < 1 or q.shape[1] > 24 or
            np.linalg.norm(q.T @ q - np.eye(q.shape[1])) > 1.e-10 or
            not np.array_equal(response[:,0], arrays["raw_base_residual"][:,0])):
        raise ValueError("current basis is not orthonormal or does not bind the original residual")
    base = array_spec(metadata, "raw_base_residual")
    physical, scaled = max(base["physical_blocks"]), base["scaled_l2"]
    if min(physical, scaled, radius) <= 0 or not np.isfinite([physical, scaled, radius]).all():
        raise ValueError("invalid original model scales or radius")
    normalization = arrays["response_column_normalization"][:,0]
    parts = []
    for i, block in enumerate(BLOCKS):
        gram = arrays[f"normalized_physical_gram_{i}"]
        values, vectors = np.linalg.eigh(gram)
        peak = max(float(values[-1]), 1.e-300)
        if values[0] < -1.e-10*peak:
            raise ValueError(f"non-Hilbert native block {block}")
        parts.append(((np.sqrt(np.maximum(values, 0.))[:,None]*vectors.T)*normalization[None,:])/physical)
    parts.append(np.linalg.qr(sr[:,None]*response/scaled, mode="r"))
    model = JointModel(q, parts)
    _, y, info = model.solve(radius, max_iterations=100)
    observed = float(np.max(np.abs(q @ y)))
    if not np.isfinite(y).all() or observed > radius*(1.+1.e-10):
        raise ValueError("nonfinite or infeasible radius proposal")
    lower, gap, iterations = (info[k] for k in ("model_lower", "model_gap", "model_iterations"))
    if gap is not None and gap <= 2.e-4:
        status = "GAP_TARGET"
    elif iterations == 100:
        status = "ITERATION_CAP"
    else:
        status = "LP_LIMIT"
    info.update(status=status, basis_sha256=basis_sha, actual_scaled_size=observed,
                base_physical=physical, base_scaled=scaled)
    with Path(str(output)+".json").open("x") as f:
        json.dump(info, f, indent=2, allow_nan=False)
        f.write("\n")
    with output.open("x") as f:
        print("WJH27BoundedProposalV1", basis_sha, format(radius, ".17g"), len(y),
              format(info["model_score"], ".17g"), int(lower is not None),
              format(lower or 0., ".17g"), iterations, status, file=f)
        for value in y:
            print(format(value, ".17g"), file=f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("basis", type=Path)
    parser.add_argument("radius", type=float)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    choose(args.basis, args.radius, args.output)
