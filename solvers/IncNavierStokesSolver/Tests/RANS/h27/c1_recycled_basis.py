#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare a signed recycled span and its k/epsilon block enrichment."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from c1_action_gap import read_audit
from c1_bounded_trials import load, packed, save
from resources import identity


def digest(data):
    return hashlib.sha256(data).hexdigest()


def prepare(checkpoint, receipt, directory, retained_only=False):
    packet = load(checkpoint)
    state = packet["state"]
    _, ordering = packed(packet, state["state_sha256"])
    metadata, arrays = read_audit(receipt)
    controller = packet["merit_controller"]
    if (metadata["stage"] != "bounded-current-basis" or
            metadata["free_id_sha256"] != ordering or
            metadata["equation_identity"] != packet["equation_identity"] or
            metadata["reference_identity"] != state["reference_identity"] or
            metadata["bounded_policy_sha256"] != controller["policy_sha256"] or
            metadata["quadrature"] != [state["quadratureX"], state["quadratureR"]]):
        raise ValueError("recycled span must come from the same signed original equation and policy")
    dx = np.asarray(controller["accepted_displacement"], dtype="<f8")
    sx, sr = arrays["unknown_scale"][:,0], arrays["row_scale"][:,0]
    if digest(dx.tobytes()) != controller["accepted_displacement_sha256"] or dx.shape != sx.shape:
        raise ValueError("accepted displacement differs")
    v = dx/sx
    q = arrays["scaled_correction_basis"].copy()
    base_dim = q.shape[1]
    if np.linalg.norm(q.T @ q-np.eye(base_dim)) >= 1e-10 or np.any(sx <= 0) or np.any(sr <= 0):
        raise ValueError("invalid recycled span or scales")
    offsets = np.r_[0, np.cumsum(metadata["free_counts"])]
    enriched = q.copy()
    additions = []
    blocks = () if retained_only else ((4, "k"), (5, "outer_epsilon"))
    for block, name in blocks:
        axis = np.zeros_like(v)
        axis[offsets[block]:offsets[block+1]] = v[offsets[block]:offsets[block+1]]
        size = np.linalg.norm(axis)
        for _ in range(2):
            axis -= enriched @ (enriched.T @ axis)
        outside = np.linalg.norm(axis)
        kept = bool(outside > 1e-10*max(size, 1e-30))
        if kept:
            enriched = np.column_stack((enriched,axis/outside))
        additions.append({"block":name, "uncorrected_norm":float(size),
                          "outside_span_norm":float(outside), "kept":kept})
    if not retained_only and not (base_dim < enriched.shape[1] <= 24):
        raise ValueError("block split supplied no bounded new direction")
    directory.mkdir()
    frames = [("recycled",q,[])]
    if not retained_only:
        frames.append(("enriched",enriched,additions))
    for name, matrix, components in frames:
        out = directory/(name+".bin")
        with out.open("xb") as stream:
            stream.write(matrix.astype("<f8",copy=False).tobytes(order="F"))
        payload = identity(out)
        binding = digest("|".join(["WJH27RecycledBasisV1",state["reference_identity"],ordering,
            digest(sx.astype("<f8",copy=False).tobytes()),digest(sr.astype("<f8",copy=False).tobytes()),
            str(state["quadratureX"]),str(state["quadratureR"]),state["geometry"]["normal_map_signature"],
            state["outlet_identity"],payload["sha256"]]).encode())
        save(directory/(name+".json"),{
            "schema":"WJH27RecycledBasisV1", "checkpoint":identity(checkpoint),
            "current_checkpoint_coefficient_sha256":state["state_sha256"], "previous_basis":identity(receipt),
            "previous_basis_coefficients":metadata["base_coefficient_sha256"],"source_equation":packet["equation_identity"],
            "accepted_displacement_sha256":controller["accepted_displacement_sha256"],
            "binding":binding, "input":payload, "columns":matrix.shape[1], "added_axes":components,
            "orthogonality_error":float(np.linalg.norm(matrix.T @ matrix-np.eye(matrix.shape[1])))})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for param in ("checkpoint", "basis", "directory"):
        parser.add_argument(param, type=Path)
    parser.add_argument("--retained-only",action="store_true")
    args = parser.parse_args()
    prepare(args.checkpoint.resolve(),args.basis.resolve(),args.directory.resolve(),args.retained_only)
