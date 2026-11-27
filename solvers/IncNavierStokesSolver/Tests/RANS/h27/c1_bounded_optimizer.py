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

HELPER_SHA = "4f63b4952dab67c7e962c101ee1f2e66af3e7dcc9fcfc9b76e9241c60b681f3c"
HELPERS = Path(__file__).resolve().parents[5]/"docs/plans/wj-h27-linux/b3-frozen-size"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def incumbents(receipt, output, metadata, arrays, model, archived):
    q = model.q
    for path in archived:
        binding, data = read_audit(path)
        names = ("base_coefficient_sha256", "free_id_sha256", "equation_identity", "reference_identity")
        if (any(binding[name] != metadata[name] for name in names) or
                any(not np.array_equal(data[name], arrays[name])
                    for name in ("base", "unknown_scale", "row_scale", "raw_base_residual"))):
            raise ValueError("archived incumbent has a different accepted base or scaling")
        coordinates = q.T @ (data["direction"][:,0]/arrays["unknown_scale"][:,0])
        yield path.name, np.linalg.solve(model.transform, coordinates)
    if "scaled_previous_displacement" in arrays:
        coordinates = q.T @ arrays["scaled_previous_displacement"][:,0]
        yield "reprojected accepted displacement", np.linalg.solve(model.transform, coordinates)
    # Only same-basis proposals are portable as coordinates across radii.
    siblings = sorted(receipt.parent.glob(receipt.name+".proposal-*")) if output.parent == receipt.parent else ()
    for path in siblings:
        if path == output or path.suffix == ".json":
            continue
        with path.open() as stream:
            header = stream.readline().split()
            y = np.asarray([float(x) for x in stream.read().split()])
        if (len(header)!=9 or header[:2] != ["WJH27BoundedProposalV2", sha(receipt)] or
                int(header[3])!=q.shape[1] or y.shape!=(q.shape[1],) or not np.isfinite(y).all()):
            raise ValueError("stored proposal has an invalid same-basis binding")
        yield path.name, np.linalg.solve(model.transform, y)


def choose(receipt, radius, output, seed_model=None, *, iterations=None, seconds=None,
           relative_gap=None, archived=()):
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
    gram_minima = []
    for i, block in enumerate(BLOCKS):
        gram = arrays[f"normalized_physical_gram_{i}"]
        values, vectors = np.linalg.eigh(gram)
        peak = max(float(values[-1]), 1.e-300)
        if values[0] < -1.e-10*peak:
            raise ValueError(f"non-Hilbert native block {block}")
        gram_minima.append(float(values[0]/peak))
        parts.append(((np.sqrt(np.maximum(values, 0.))[:,None]*vectors.T)*normalization[None,:])/physical)
    parts.append(np.linalg.qr(sr[:,None]*response/scaled, mode="r"))
    model = JointModel(q, parts)
    previous = seed_model(model, radius) if seed_model is not None else None
    model_ratio = iterations is not None
    if model_ratio:
        if not (1 <= iterations <= 1000 and 0 < seconds <= 120 and 0 < relative_gap < 1):
            raise ValueError("invalid reduced optimization budget")
        seeds = list(incumbents(receipt, output, metadata, arrays, model, archived))
        if previous is not None:
            seeds.append(("current active gradients", previous))
        _, y, info = model.solve(radius, max_iterations=iterations, relative_gap=relative_gap,
                                 seconds=seconds, incumbents=seeds)
    else:
        _, y, info = model.solve(radius, previous=previous, max_iterations=100)
    if seed_model is not None:
        info["incumbent_method"] = "current active gradients"
        info["seed_model_score"] = None if previous is None else float(
            max(np.linalg.norm(b + a @ previous) for b, a in model.parts))
    observed = float(np.max(np.abs(q @ y)))
    if not np.isfinite(y).all() or observed > radius*(1.+1.e-10):
        raise ValueError("nonfinite or infeasible radius proposal")
    lower, used, status = (info[k] for k in ("model_lower", "model_iterations", "status"))
    info.update(basis_sha256=basis_sha, actual_scaled_size=observed,
                base_physical=physical, base_scaled=scaled, gram_minima_relative=gram_minima,
                model_singular_values=model.ranks.tolist(), declared_iterations=iterations,
                declared_seconds=seconds, declared_relative_gap=relative_gap)
    with Path(str(output)+".json").open("x") as f:
        json.dump(info, f, indent=2, allow_nan=False)
        f.write("\n")
    with output.open("x") as f:
        print("WJH27BoundedProposalV2" if model_ratio else "WJH27BoundedProposalV1",
              basis_sha, format(radius, ".17g"), len(y),
              format(info["model_score"], ".17g"), int(lower is not None),
              format(lower or 0., ".17g"), used, status, file=f)
        for value in y:
            print(format(value, ".17g"), file=f)


def parse(description):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("basis", type=Path)
    parser.add_argument("radius", type=float)
    parser.add_argument("output", type=Path)
    parser.add_argument("--iterations",type=int)
    parser.add_argument("--seconds",type=float)
    parser.add_argument("--relative-gap",type=float)
    parser.add_argument("--incumbent",type=Path,action="append",default=[])
    args = parser.parse_args()
    limits = [args.iterations, args.seconds, args.relative_gap]
    if any(x is None for x in limits) and (not all(x is None for x in limits) or args.incumbent):
        parser.error("declare the iteration, time and relative-gap controls together")
    return args


if __name__ == "__main__":
    args = parse(__doc__)
    choose(args.basis, args.radius, args.output, iterations=args.iterations,
           seconds=args.seconds, relative_gap=args.relative_gap, archived=args.incumbent)
