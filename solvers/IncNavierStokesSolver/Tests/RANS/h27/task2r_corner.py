#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Diagnose the mixed ground–outlet corner of the unforced Task 2R reference."""

import argparse
import cmath
import math
from pathlib import Path
import time

import numpy as np
from scipy.optimize import brentq

from physical_oracle import PhysicalOracle, load_native, load_reference, require
from resources import atomic, identity


def pencil(nu, mu):
    # (x,r) in the right-angle wedge, with x increasing into the domain.
    # For u=psi_r and w=-psi_x the frozen lower operator is
    # mu*(psi_xxxx+psi_rrrr)+(4*nu-2*mu)*psi_xxrr.
    # Use conjugate pairs of psi=(x+a*r)**(1+lambda). The first two rows
    # impose psi=psi_r=0 at r=0. The other rows impose shear and normal
    # traction zero at x=0 after a common nonzero monomial is removed.
    theta = math.asin(math.sqrt(nu/mu))

    def boundary_matrix(exponent):
        columns = []
        for angle in (theta, math.pi-theta):
            a = cmath.exp(1j*angle)
            traction = (a*a-1)*cmath.exp(1j*(exponent-1)*angle)
            column = np.array([1, a, traction, -a*traction])
            columns.extend((column.real, column.imag))
        return np.column_stack(columns)

    lambdas = np.linspace(.002, .998, 249)
    dets = [float(np.linalg.det(boundary_matrix(v))) for v in lambdas]
    roots = []
    for left, right, fa, fb in zip(lambdas[:-1], lambdas[1:], dets[:-1], dets[1:]):
        if fa*fb < 0:
            root = brentq(lambda v: np.linalg.det(boundary_matrix(v)), left, right, xtol=1e-14)
            sigma = np.linalg.svd(boundary_matrix(root), compute_uv=False)
            roots.append(dict(velocity_power=root, gradient_pressure_power=root-1,
                              boundary_pencil_relative_null_residual=float(sigma[-1]/sigma[0])))
    return dict(nu=nu, frozen_wall_shear_coefficient=mu, root_interval=[.002,.998], real_roots=roots)


def measure(path, ref, deadline):
    data, metadata = load_native(path)
    require(not data["manufactured"] and metadata["steady_mode"] and data["time"] == 0.,
            "corner requires a converged unforced zero-time reference")
    origin = PhysicalOracle(data, ref, 20, deadline)
    eps = 1e-12
    x_right = max(e["bounds"][1] for e in origin.elements)
    exterior = [e for e in origin.elements if abs(e["bounds"][1]-x_right) <= eps]
    ground = [e for e in origin.elements if abs(e["bounds"][2]) <= eps]
    vertex = max(ground, key=lambda e:e["bounds"][1])
    x_vertex = np.array([x_right])
    K = origin.trace(vertex, "k", x_vertex)
    mu = ref.nu + float(ref.aux(K, np.zeros(1))[0])

    def evaluate(owner, x, r):
        a = origin.sample(owner, np.array([x]), np.array([r]))
        grad = a["G"][0]
        traction_x = a["V"][0,0,0]-a["Pi"][0]
        return dict(x=x, r=r, U=float(a["U"][0]), W=float(a["W"][0]), Ux=float(grad[0,0]),
                    Ur=float(grad[0,2]), Wx=float(grad[2,0]), Wr=float(grad[2,2]),
                    Pi_over_tau=float(a["Pi"][0]/ref.tau), traction_x_over_tau=float(traction_x/ref.tau),
                    traction_r_over_tau=float(a["V"][0,0,2]/ref.tau), bounds=owner["bounds"])

    distance = [.032/2**n for n in range(13)] + [0.]
    along_outlet, along_ground = [], []
    for d in distance:
        outlet_owner = next(e for e in exterior if e["bounds"][2]-eps <= d <= e["bounds"][3]+eps)
        corner_owner = next(e for e in ground if e["bounds"][0]-eps <= x_right-d <= e["bounds"][1]+eps)
        along_outlet.append(evaluate(outlet_owner, x_right, d))
        along_ground.append(evaluate(corner_owner, x_right-d, 0.))
    return dict(state=identity(path), vertex_owner=vertex["bounds"], sampled_interface_K=float(K[0]),
                frozen_lower_pencil=pencil(ref.nu, mu), outflow=along_outlet, ground=along_ground)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", required=True, type=Path)
    p.add_argument("--native", type=Path, required=True, action="append")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--seconds", type=float, required=True)
    a = p.parse_args()
    deadline = time.monotonic()+a.seconds
    ref = load_reference(a.table)
    require(ref.exact_inner, "corner requires the exact-inner reference")
    atomic(a.out, dict(schema="WJH27Task2RCornerPencilDiagnosticV1", table=identity(a.table),
        kinetic_scope="Frozen local lower-wall velocity-pressure principal operator. A real root in (0,1) permits a divergent stress mode.",
        limitation="The mode amplitude in the coupled variable-coefficient solution and the effect of sampled interface forcing on it are unverified.",
        measurements={n.parent.name:measure(n, ref, deadline) for n in a.native}))


if __name__ == "__main__":
    main()
