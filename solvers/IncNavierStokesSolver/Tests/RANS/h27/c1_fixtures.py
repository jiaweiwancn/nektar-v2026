#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare small, immutable standard-main C1 cases."""

import argparse
import json
from pathlib import Path

from coupled_fixtures import generate, mesh
from resources import identity
from task2r_verify import edit


def bisect(axis):
    return sorted(axis + [(lo + hi) / 2 for lo, hi in zip(axis, axis[1:])])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("flat", "hill", "manufactured"))
    parser.add_argument("out", type=Path)
    parser.add_argument("--height", type=float)
    parser.add_argument("--initial-slowdown", type=float, default=0.)
    parser.add_argument("--level", type=int, choices=range(4), default=0)
    parser.add_argument("--dt", type=float, default=.02)
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--restart", type=Path)
    parser.add_argument("--quadrature", type=int, nargs=2)
    parser.add_argument("--recovery", type=int, default=1)
    parser.add_argument("--science", type=int, default=0)
    parser.add_argument("--step-seconds", type=float, default=90.)
    parser.add_argument("--steady-seconds", type=float, default=0.)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    parser.add_argument("--flat-map-control", action="store_true")
    parser.add_argument("--trace-derivative-check", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True)
    plans = Path(__file__).resolve().parents[5] / "docs/plans/wj-h27-linux"
    table = plans / "task2r-r1/tables-exact-inner-v1.json"
    manufactured = args.case == "manufactured"
    if manufactured:
        xs = [0., .5, 1.]
        rs = [.0875 + r for r in [0., .02435, .06, .1125, .15, .25, .375, .5]]
    else:
        admission = json.loads((plans / "task2-outlet-qual/cp2/outlet-qual-bf-admission-1.json").read_text())
        rs = admission["physical_z_edges"]
        xs = ([-10., -5., 0., 5., 15.] if args.case == "flat" else
              [-10., -5., -2.5, -1.875, -1.25, -.625, 0., .625, 1.25,
               1.875, 2.5, 4., 6., 9., 12., 15.])
    for _ in range(args.level):
        xs, rs = bisect(xs), bisect(rs)
    geometry = out / "geometry.xml"
    mesh(geometry, xs, rs, quadrature=12)
    session = out / "run.xml"
    qx, qr = args.quadrature or ((12, 16) if manufactured else (6, 8))
    generate(session, table, out / "state", mode="Manufactured" if manufactured else "Physical",
             template=geometry, dt=args.dt, steps=0 if args.steady_seconds else args.steps, science=args.science,
             recovery=args.recovery, quadrature_x=qx, quadrature_r=qr,
             solve="SteadyReference" if args.steady_seconds else "Advance", steady_seconds=args.steady_seconds,
             nonlinear=args.tolerance, check_current_actions=args.trace_derivative_check,
             backend="FastFull3DV1", outlet_policy="FullTractionV1" if manufactured else "ZeroTransverseV1",
             restart=args.restart.resolve() if args.restart else "")
    info = dict(WJReferencePolicy="ExactInnerFiniteNuV1", WJImplicitSolver="CoupledJFNKV1",
                WJNewtonAdvection="Implicit", WJStepPreconditioner="CurrentDataDriftV1")
    parameters = dict(WJStepNewtonSeconds=args.step_seconds, WJStepPCMaximumUses=3,
                      WJStepPCRelativeDrift=.02, WJStepFDStep=4e-6,
                      WJStepForcingMultiplier=4, WJStepForcingCeiling=.3)
    if args.case != "flat" and not args.flat_map_control:
        info["WJGeometry"] = "NormalCosineHillV1"
        parameters["WJHillHeight"] = args.height if args.height is not None else (.1 if manufactured else 1.)
        parameters["WJHillInitialSlowdown"] = args.initial_slowdown
    if manufactured:
        info["WJMMSVariant"] = "NormalTubeV1"
    edit(session, info, parameters)
    files = [geometry, session, out / "run.xml.opt", table]
    (out / "inputs.json").write_text(json.dumps([identity(file) for file in files], indent=2) + "\n")


if __name__ == "__main__":
    main()
