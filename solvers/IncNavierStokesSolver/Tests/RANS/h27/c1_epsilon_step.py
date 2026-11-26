#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently replay one accepted diagnostic direction on a selected outer epsilon row."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

from c1_oracle import Audit, coefficients, dot, load_c1_state, state_hash
import physical_oracle as po
from resources import identity


def matching_project(audit, cells, k):
    """Reconstruct the represented matching Dirichlet modes at the trial K."""
    trace = {}
    q, w = leggauss(16)
    for cell in cells:
        a0, a1, r0, _ = cell["bounds"]
        if abs(r0-audit.ref.H) > 1e-12:
            continue
        edge = list(po.EDGES[2])
        r = lambda a: np.full_like(a, audit.ref.H)
        a = (a0+a1)/2+(a1-a0)*q/2
        length, _ = audit.hill.face(cell, 2, a, r(a))
        weights = length*w*(a1-a0)/2
        ae = np.array([a0, a1])
        endpoint_basis = po.basis(cell, "epsilon", ae, r(ae))[0][:, edge[:2]]
        target = lambda p: (po.basis(cell, "k", p, r(p))[0]@k[cell["k"]])**1.5/(audit.ref.cl*audit.ref.H)
        endpoints = np.linalg.solve(endpoint_basis, target(ae))
        b = po.basis(cell, "epsilon", a, r(a))[0][:, edge]
        tests = b[:, 2:]
        modes = np.linalg.solve(tests.T@(weights[:, None]*tests), tests.T@(weights*(target(a)-b[:, :2]@endpoints)))
        ids = np.asarray(cell["epsilon"])[edge]
        for idx, value in zip(ids, [*endpoints, *modes]):
            if idx in trace:
                po.require(abs(value-trace[idx]) < 2e-13, "independent matching owners disagree")
            trace[idx] = value
    return trace


def outer_row(audit, owners, row):
    """Integrate signed advection, diffusion, production and destruction."""
    audit.clock()
    result = np.zeros(4)
    for cell in owners:
        j = cell["epsilon"].index(row)
        a, r, w = audit.rule(cell)
        metric = audit.hill.metrics(cell, a, r)
        basis, u, _, _, _, prod, scalars, diffusion, grad, _ = audit.sample(cell, a, r, metric)
        value, dq = basis["epsilon"]
        E, k = scalars["epsilon"], scalars["k"]
        phi, g = value[:, j], dq[:, j, :]
        result += (w*metric["det"]) @ np.column_stack((
            -E*dot(g, u), diffusion["epsilon"]*dot(g, grad["epsilon"]),
            -1.44*phi*prod*E/k, 1.82*phi*E**2/k))
        xl, xr, rl, rr = cell["bounds"]
        exterior = (xl-audit.hill.left, xr-audit.hill.right, rl-audit.ref.H, rr-audit.top)
        for side, distance in enumerate(exterior):
            if abs(distance) > 1e-12:
                continue
            a, r, w = audit.face(cell, side)
            length, normal = audit.hill.face(cell, side, a, r)
            velocity = audit.velocity(po.basis(cell, "u", a, r)[0], np.array(cell["u"]), audit.s)
            b = po.basis(cell, "epsilon", a, r)[0]
            result[0] += (w*length) @ (b[:, j] * dot(velocity, normal) * (b@audit.s["epsilon"][cell["epsilon"]]))
            # The prescribed top flux is constant and cancels in the remainder.
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("base", "accepted", "ledger", "table", "out"):
        parser.add_argument("--"+name, type=Path, required=True)
    parser.add_argument("--row", type=int, required=True)
    parser.add_argument("--iteration", type=int, required=True)
    args = parser.parse_args()
    i = args.iteration
    x, mx = load_c1_state(args.base, True)
    y, my = load_c1_state(args.accepted, True)
    for meta, sequence in ((mx, i), (my, i+1)):
        po.require(meta["iterate_custody"] == "LAST_ACCEPTED" and meta["regularization"]["eta"] == 0.
                   and meta["accepted_checkpoint"]["sequence"] == sequence, "expected consecutive original accepted states")
    for key in ("signature", "geometry", "elements", "source_identity"):
        po.require(x[key] == y[key], "accepted states have different "+key)
    po.require(mx["equation_identity"] == my["equation_identity"] and
               mx["initial_packet_sha256"] == my["initial_packet_sha256"] and
               mx["accepted_checkpoint"]["session_sha256"] == my["accepted_checkpoint"]["session_sha256"],
               "different run lineage")
    line, = [row for row in json.loads(args.ledger.read_text())["line_search"]
             if row["newton"] == i and row["accepted"]]
    for before, after in zip(line["actual_blocks"], my["original_residuals"]):
        po.require(abs(before-after) <= 5e-12, "ledger does not bind the accepted residual")
    alpha = line["step"]
    audit = Audit(x, po.load_reference(args.table), 28, time.monotonic()+60)
    po.require(not x["manufactured"] and args.row in audit.free["epsilon"], "requires a physical original free outer row")
    owners = [e for e in audit.mesh if args.row in e["epsilon"]]
    base, accepted = audit.s, coefficients(y["coefficients"])
    po.require(state_hash(accepted) == y["state_sha256"], "accepted coefficient digest differs")
    increments = {k: (accepted[k]-base[k])/alpha for k in base}
    lift = matching_project(audit, owners, base["k"])
    po.require(lift, "row must touch the matching lift through an owner")
    ids = np.array(list(lift))
    endpoint_errors = []
    for state in (base, accepted):
        projected = matching_project(audit, owners, state["k"])
        endpoint_errors.append(float(max(abs(v-state["epsilon"][j]) for j,v in projected.items())))
    po.require(max(endpoint_errors) <= 2e-13, "independent matching coefficient reconstruction differs")

    def sample(fraction):
        audit.s = {key: val+fraction*increments[key] for key,val in base.items()}
        updated = matching_project(audit, owners, audit.s["k"])
        audit.s["epsilon"][ids] = [updated[j] for j in ids]
        current = outer_row(audit, owners, args.row)
        audit.s["epsilon"][ids] = base["epsilon"][ids]
        frozen = outer_row(audit, owners, args.row)
        return np.stack((current, current-frozen))

    f0, fa = sample(0.), sample(alpha)
    h = min(4e-6, alpha/8.)
    derivatives = [(sample(step)-sample(-step))/(2*step) for step in (h, h/2)]
    remainder = fa-f0-alpha*derivatives[1]
    with args.out.open("x") as stream:
        json.dump(dict(schema="WJH27AcceptedEpsilonRowV1", iteration=i, row=args.row, accepted_step=alpha,
                       inputs=[identity(p) for p in (args.base, args.accepted, args.ledger, args.table)],
                       coefficients=[x["state_sha256"], y["state_sha256"]], integration_order=28,
                       owners=[dict(id=e["geometry"], bounds=e["bounds"]) for e in owners],
                       reconstruction_error=endpoint_errors, fraction_spacing=[h,h/2],
                       maximum_row_evaluations=12, terms=["advection","diffusion","production","destruction"],
                       signed_row_remainder=remainder[0].tolist(), signed_total=float(remainder[0].sum()),
                       dependent_lift_overlap=remainder[1].tolist(),
                       spacing_contribution=(alpha*abs(derivatives[0]-derivatives[1])).tolist(),
                       native_full_epsilon_remainder=line["nonlinear_remainder_blocks"][4],
                       comparison_scope="One selected unscaled signed weak row."),
                  stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
