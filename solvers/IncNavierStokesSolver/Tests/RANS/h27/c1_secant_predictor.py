#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Extrapolate two constrained remaps using their common free inventory."""

import argparse
import copy
import hashlib
import json
from pathlib import Path
import struct

from resources import identity


def load(path):
    return json.loads(path.read_text(), parse_int=lambda x: -0.0 if x == "-0" else int(x))


def prepare(args):
    roots = [load(path) for path in args.roots]
    moved = [load(path) for path in args.transfers]
    for path, root, remap in zip(args.roots, roots, moved):
        if not (root["role"] == remap["role"] == "STEADY_SEARCH" and
                root.get("iterate_custody") != "UNCOMMITTED_TRIAL" and
                max(root["original_residuals"]) <= 1e-10 and root["regularization"]["eta"] == 0. and
                remap["iterate_custody"] == "LAST_ACCEPTED" and remap["regularization"]["eta"] == 0. and
                remap["initial_packet_sha256"] == identity(path)["sha256"]):
            raise ValueError("secant requires two exact remaps from the converged original anchors")
    old, new = (remap["state"] for remap in moved)
    for name in ("signature", "geometry", "sizes", "elements", "free_velocity_ids",
                 "k_free_ids", "epsilon_free_ids", "quadratureX", "quadratureR"):
        if old[name] != new[name]:
            raise ValueError(f"secant remaps differ in {name}")
    g0, g1 = (root["state"]["geometry"]["hill_height"] for root in roots)
    target = new["geometry"]["hill_height"]
    if not g0 < g1 < target:
        raise ValueError("secant requires increasing anchor and target heights")
    ratio = (target-g1)/(g1-g0)
    velocity, pressure, _, _ = new["sizes"]
    u, v = new["free_velocity_ids"]
    free = dict(flow=[*u, *(velocity+j for j in v),
                      *range(2*velocity, 2*velocity+pressure)],
                k=new["k_free_ids"], epsilon=new["epsilon_free_ids"])
    state = copy.deepcopy(new)
    h = hashlib.sha256()
    for name in ("flow", "k", "epsilon"):
        values = state["coefficients"][name]
        before = old["coefficients"][name]
        for j in free[name]:
            values[j] += ratio*(values[j]-before[j])
        for value in values:
            h.update(struct.pack("<d", value))
    state["state_sha256"] = h.hexdigest()
    for key in ("interface_K_over_Kin", "interface_K_table_margin_over_Kin", "interface_cutoff_range",
                "essential_flow_lift", "essential_k_lift", "essential_epsilon_lift"):
        state[key] = None
    state["interface_derived_status"] = "UNEVALUATED_PREDICTOR"
    packet = dict(schema="WJH27DiagnosticCandidateV1", role="STEADY_SEARCH", physical_history=False,
                  iterate_custody="UNCOMMITTED_TRIAL", reason="constrained secant predictor",
                  residual_kind="UNEVALUATED_SECANT_GUESS", attempted_step=0, attempted_time=0.,
                  attempted_iteration=-1, residuals=None, regularization={"eta":0.}, state=state,
                  predictor=dict(anchor_heights=[g0,g1], target_height=target, ratio=ratio,
                                 roots=[identity(p) for p in args.roots],
                                 transfers=[identity(p) for p in args.transfers]))
    with args.out.open("x") as output:
        json.dump(packet, output, separators=(",",":"), allow_nan=False)
        output.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--roots", nargs=2, type=Path, required=True, metavar=("OLDER", "NEWER"))
    parser.add_argument("--transfers", nargs=2, type=Path, required=True, metavar=("OLDER", "NEWER"))
    parser.add_argument("--out", type=Path, required=True)
    prepare(parser.parse_args())
