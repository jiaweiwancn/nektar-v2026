#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reconstruct rejected WJ fields and compare one-sided constitutive witnesses."""

import argparse
import json
from pathlib import Path
import re
import time

import numpy as np
from scipy.optimize import minimize_scalar

from c1_margins import response
from c1_oracle import Audit, load_c1_state
import physical_oracle as po
from resources import identity


def point(audit, cell, a, radius, side):
    x, r = np.array([a]), np.array([radius])
    outer = bool(cell["epsilon"])
    metric = audit.hill.metrics(cell, x, r, outer=outer)
    _, velocity, gradient, _, _, production, scalars, _, _, _ = audit.sample(cell, x, r, metric)
    k, epsilon = scalars["k"][0], scalars["epsilon"][0]
    K = audit.trace(cell, "k", x)[0]
    nu = audit.ref.nu
    aux = 0. if outer else audit.ref.aux(np.array([K]), r)[0]
    q = min(1., max(0., (radius/audit.ref.H-.4)/.35))
    chi = q**3*(10-15*q+6*q*q)
    cartesian = [0, 2]
    G = np.zeros((3, 3))
    G[np.ix_(cartesian, cartesian)] = gradient[0]
    normal, tangent = np.zeros(3), np.zeros(3)
    normal[cartesian], tangent[cartesian] = metric["n"][0], metric["t"][0]
    wall_projection = np.outer(tangent, normal) + np.outer(normal, tangent)
    I = np.eye(3)

    def covariance(g):
        shear = np.einsum("i,nij,j->n", tangent, g+g.swapaxes(1, 2), normal)
        wall = 2*k*I/3 - aux*shear[:, None, None]*wall_projection
        wj, root, defect, derivative = po.full_wj(
            np.full(len(g), k), np.full(len(g), epsilon), g)
        return (1-chi)*wall + chi*wj, root, defect, derivative

    R, root, root_defect, root_derivative = covariance(G[None, :, :])
    M = G+G.T
    S = (k/epsilon)*(M/2-np.trace(G)*I/3)
    O = (.5*k/epsilon)*(G-G.T)
    inv_s, inv_o, inv_4 = map(float, (np.trace(S@S), np.trace(O@O), np.trace(S@O@O)))
    mu = nu+(1-chi)*aux+chi*.6*k*k/epsilon*root[0]/(root[0]**2-2*inv_o)

    def direction(angle):
        angle = np.atleast_1d(angle)
        s, c, z = np.sin(angle), np.cos(angle), np.zeros_like(angle)
        wave = np.column_stack((c, z, s))
        polarization = np.column_stack((-s, z, c))
        return np.einsum("ni,nj->nij", polarization, wave), wave, polarization

    def minimize(relative):
        spacing = relative*epsilon/k

        def raw(theta):
            D, _, _ = direction(theta)
            values = []
            for j in (-2, -1, 1, 2):
                g = G+j*spacing*D
                values.append(nu*(g+g.swapaxes(1, 2))-covariance(g)[0])
            dv = (values[0]-8*values[1]+8*values[2]-values[3])/(12*spacing)
            return np.einsum("nij,nij->n", dv, D)

        grid = np.arange(240)*np.pi/240
        values = raw(grid)
        i = int(values.argmin())
        opt = minimize_scalar(lambda x: raw(x)[0],
            bounds=(grid[i]-np.pi/240, grid[i]+np.pi/240), method="bounded",
            options={"xatol": 1e-12})
        _, wave, polarization = direction(opt.x)
        return dict(spacing=spacing, angle=opt.x % np.pi, raw=opt.fun, over_mu=opt.fun/mu,
                    wave=wave[0].tolist(), polarization=polarization[0].tolist())

    minima = [minimize(relative) for relative in (1e-3, 5e-4, 2.5e-4)]
    final = minima[-1]["over_mu"]
    return dict(geometry=cell["geometry"], bounds=cell["bounds"], side=side, a=a, rho=radius,
                position=metric["position"][0].tolist(), outer=outer, normal=normal.tolist(),
                velocity=velocity[0].tolist(), gradient=G.tolist(), divergence=float(np.trace(G)),
                K=float(K), k=float(k), epsilon=float(epsilon), k_over_epsilon=float(k/epsilon),
                chi=chi, physical_auxiliary=float(aux), production=float(production[0]),
                mu=float(mu), wj=dict(root=float(root[0]), root_defect=float(root_defect[0]),
                    root_derivative=float(root_derivative[0]), strain_squared=inv_s,
                    rotation_trace_square=inv_o, mixed_invariant=inv_4),
                covariance=R[0].tolist(), covariance_over_k=float(np.linalg.eigvalsh(R[0])[0]/k),
                minima=minima, differentiation_spread=max(abs(v["over_mu"]-final) for v in minima))


def neighborhood(audit, cell, a, r, samples):
    xl, xr, rl, rr = cell["bounds"]
    foot = audit.hill.L
    near_foot = min(abs(xl-foot), abs(xr-foot)) < 1e-10
    def nearby(c):
        l, h, low, high = c["bounds"]
        r_near = low <= rr+1e-10 and high >= rl-1e-10
        if near_foot:
            return min(abs(l-foot), abs(h-foot)) < 1e-10 and r_near
        return l <= xr+1e-10 and h >= xl-1e-10 and r_near

    candidates = [c for c in audit.mesh if nearby(c)]
    near = []
    best = float("inf")
    witness = None
    for c in candidates:
        audit.clock()
        l, h, low, high = c["bounds"]
        if near_foot and low <= r <= high:
            side = "foot-right" if abs(l-foot)<1e-10 else "foot-left"
            near.append((c, foot, r, side))
        if l <= a <= h:
            for v in (rl, rr):
                if abs(v-low)<1e-10 or abs(v-high)<1e-10:
                    near.append((c, a, v, "radial-lower-owner" if abs(v-high)<1e-10 else "radial-upper-owner"))
        axis = np.linspace(l, h, samples)
        radial = np.linspace(low, high, samples)
        ax, rad = (array.ravel() for array in np.meshgrid(axis, radial))
        metric = audit.hill.metrics(c, ax, rad, outer=bool(c["epsilon"]))
        _, _, G, _, _, _, scalar, _, _, _ = audit.sample(c, ax, rad, metric)
        metric["K"] = audit.trace(c, "k", ax)
        val, index, *_ = response(audit.ref, G, scalar["k"], scalar["epsilon"], rad, metric, bool(c["epsilon"]))
        if val < best:
            best = val
            witness = (c, float(ax[index]), float(rad[index]), "local-grid-minimum")
    return near, witness, len(candidates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("native", "table", "out"):
        parser.add_argument("--"+arg, type=Path, required=True)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--samples", type=int, default=25)
    parser.add_argument("--seconds", type=float, default=180.)
    args = parser.parse_args()
    deadline = time.monotonic()+args.seconds
    ref = po.load_reference(args.table)
    data, meta = load_c1_state(args.native, True)
    audit = Audit(data, ref, 12, deadline)
    match = re.search(r"geometry=(\d+), x=([^,]+), r=([^,]+), .*?principal/mu=([^,]+)", meta["reason"])
    po.require(match is not None, "packet does not report a native constitutive witness")
    gid, a, r, principal = map(float, match.groups())
    cell, = [e for e in audit.mesh if e["geometry"] == gid]
    first = point(audit, cell, a, r, "first-native-rejection")
    sides, dense, cells = neighborhood(audit, cell, a, r, args.samples)
    _, denser, _ = neighborhood(audit, cell, a, r, 2*args.samples-1)
    comparisons = [point(audit, *sample) for sample in (*sides, dense, denser)]
    previous = None
    substitution = {}
    if args.previous:
        old, _ = load_c1_state(args.previous)
        po.require(old["geometry"] == data["geometry"] and old["elements"] == data["elements"]
                   and old["state_sha256"] == meta["input_history_identity"]["current_sha256"],
                   "previous checkpoint does not match rejected history")
        old_audit = Audit(old, ref, 12, deadline)
        previous = point(old_audit, cell, a, r, "accepted-at-first-rejection-location")
        x, radius = np.array([a]), np.array([r])
        metric = audit.hill.metrics(cell, x, radius, outer=first["outer"])
        for label, gradients, scalars in (("rejected_gradient_only", first, previous),
                                           ("rejected_scalars_only", previous, first)):
            gradient = np.asarray(gradients["gradient"])[np.ix_([0, 2], [0, 2])][None, :, :]
            metric["K"] = np.array([scalars["K"]])
            margin, _, uncertainty, _, _ = response(ref, gradient, np.array([scalars["k"]]),
                np.array([scalars["epsilon"]]), radius, metric, first["outer"])
            substitution[label] = dict(principal_over_mu=margin, differentiation_spread=uncertainty)
    report = dict(native=identity(args.native), table=identity(args.table),
                  evaluator=identity(__file__), role=meta["role"], native_first_principal=principal,
                  coefficient_sha256=data["state_sha256"], first=first, neighboring_samples=comparisons,
                  neighborhood_cells=cells, grid_order=[args.samples, 2*args.samples-1],
                  previous=previous, pointwise_substitution=substitution)
    with args.out.open("x") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")


if __name__ == "__main__":
    main()
