#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare and assess unforced Task 2R flat references on the retained geometry."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

import coupled_fixtures as fixtures
import physical_oracle as po
from physical_reference_refinement import axes
from resources import atomic, identity
from task2r_verify import edit, read, require
from task2r_reference import D0, H, POLICY


def bisect(edges, selection):
    return sorted(edges + [(a + b) / 2 for a, b in zip(edges[:-1], edges[1:]) if selection(a, b)])


def prepare(args):
    out, table = args.out, args.table
    require(po.load_reference(table).exact_inner, "exact-inner trajectory requires the qualified exact policy")
    xs, zs = axes(args.geometry)
    if args.refine_x or args.x_band:
        for a, b in args.x_band:
            require(a in xs and b in xs and a < b, "x refinement band must use original mesh edges")
        xs = bisect(xs, lambda a, b: args.refine_x or any(lo <= a < b <= hi for lo, hi in args.x_band))
    if args.terminal_levels:
        require(1 <= args.terminal_levels <= 8, "bounded geometric refinement requires one to eight terminal levels")
        right, width = xs[-1], xs[-1]-xs[-2]
        xs = sorted(xs+[right-width/2**level for level in range(1,args.terminal_levels+1)])
    if args.wall_levels:
        require(1 <= args.wall_levels <= 6, "bounded geometric refinement requires one to six wall levels")
        wall, width = zs[0], zs[1]-zs[0]
        zs = sorted(zs+[wall+width/2**level for level in range(1,args.wall_levels+1)])
    if args.refine_to_radius:
        zs = bisect(zs, lambda a, b: b <= D0 + args.refine_to_radius + 1e-13)
    require(abs(zs[0] - D0) <= 1e-14 and any(abs(z - D0 - H) <= 1e-14 for z in zs),
            "modeled wall/matching plane differs")
    out.mkdir()
    template = args.geometry
    if args.refine_x or args.x_band or args.refine_to_radius or args.terminal_levels or args.wall_levels:
        require(args.refinement_witness and args.seed, "nested refinement needs the measured base witness and same-model seed")
        template = out / "geometry.xml"
        fixtures.mesh(template, xs, zs, quadrature=12)
    seed, source = "", ""
    if args.seed:
        state, _ = po.load_native(args.seed) if args.seed.name.endswith(".result.json") else (
            json.loads(args.seed.read_text())["state"], {})
        require(state["reference_identity"] == po.EXACT_TABLE_SHA and state["step"] == state["time"] == 0.,
                "seed must be bound to the new reference with an unadvanced clock")
        seed, source = args.seed, state["source_identity"]
    session = out / "run.xml"
    fixtures.generate(session, table, out / "state", mode="Physical", dt=.05, steps=0,
                      solve="SteadyReference", steady_seconds=args.solve_seconds,
                      steady_pc_mass=args.preconditioner_mass, template=template,
                      quadrature_x=args.qx, quadrature_r=args.qr,
                      linear=args.linear, nonlinear=args.nonlinear, maximum_iterations=60,
                      backend="FastFull3DV1", initial_fields=seed, initial_source_identity=source)
    edit(session, {"WJReferencePolicy": POLICY, "WJSteadyVelocityPC": args.velocity_preconditioner},
         {"WJSteadyKrylovDimension": args.krylov_dimension})
    atomic(out / "inputs.json", dict(
        schema="WJH27Task2RUnforcedReferenceInputsV1", geometry=identity(template),
        session=identity(session), table=identity(table), cells=[len(xs)-1, len(zs)-1],
        modes=[4, 3], quadrature=[args.qx, args.qr], tolerances=[args.linear, args.nonlinear],
        observations=po.STATIONS, physical_limits=po.ENVELOPES, numerical_limits=po.NUMERICAL,
        radius_interface_code=H, radius_interface_m=.04*H,
        plane_elevation_m=.04*(D0+H), initial_guess=identity(seed) if seed else None,
        maximum_seconds=args.solve_seconds, preconditioner_mass=args.preconditioner_mass,
        velocity_preconditioner=args.velocity_preconditioner, krylov_dimension=args.krylov_dimension,
        refinement_witness=identity(args.refinement_witness) if args.refinement_witness else None,
        x_refinement_bands=args.x_band, terminal_levels=args.terminal_levels, wall_levels=args.wall_levels,
        radial_refinement_upper=args.refine_to_radius,
        time_code=0.))


def act(args):
    deadline = time.monotonic() + args.seconds
    ref = po.load_reference(args.table)
    def oracle(path, candidate=False):
        if candidate:
            packet = read(path)
            require(packet["schema"] == "WJH27SteadyCandidateV1" and not packet["qualified_reference"],
                    "expected a typed, unqualified native steady candidate")
            data = packet["state"]
        else:
            data, _ = po.load_native(path)
        return po.PhysicalOracle(data, ref, args.order, deadline)

    first = oracle(args.native, args.candidate)
    if args.other:
        require(not args.candidate, "comparison requires two completed registered runs")
        second = oracle(args.other)
        require(first.data["source_identity"] == second.data["source_identity"], "different native model sources")
        result = first.compare(second)
        require(abs(first.data["time"] - second.data["time"]) <= 1e-13, "unequal simulation times")
        if args.control in ("quadrature", "algebraic"):
            for name, norm in result["norms"].items():
                if name.endswith("/k") and "/inner/" in "/" + name:
                    norm["limits"] = [1e-4, 1e-4]
                    norm["status"] = "MET" if max(norm["RMS"], norm["MAX"]) <= 1e-4 else "NOT_MET"
            result["status"] = "MET" if all(row["status"] == "MET" for row in result["norms"].values()) else "NOT_MET"
        result["control"] = args.control
        result["paired_native"] = identity(args.other)
    else:
        current = first.assemble()
        bc = first.boundary_checks()
        interfaces = first.interfaces()
        profiles = first.profiles()
        prior = po.load_reference(args.historical)
        current_report = current["physical_envelope"]
        # Compare the same station and owner inventory against the historical
        # inner profile for diagnostics, without changing the acceptance targets.
        old = po.Norms()
        gauss, weights = leggauss(32)
        for x in po.STATIONS:
            for e in first.elements:
                xl, xr, rl, rr = e["bounds"]
                if x < xl - 1e-12 or x > xr + 1e-12 or rl >= H - 1e-12:
                    continue
                first.clock()
                r = np.r_[rl, (rl + rr) / 2 + gauss * (rr - rl) / 2, rr]
                weight = np.r_[0., weights * (rr - rl) / 2, 0.]
                frozen = prior.inlet(r)
                found = first.sample(e, np.full_like(r, x), r)
                for field, scale in (("U", 1.), ("k", ref.Kin)):
                    old.add(f"x={x:g}/inner", field, (found[field] - frozen[field]) / scale, weight)
        result = dict(core=current, boundary_checks=bc, interfaces=interfaces,
                      profiles={key: profiles[key] for key in ("status", "stations", "norms", "ground_stations")},
                      historical_inner_diagnostic=old.report({}),
                      equations="MET" if current["steady_residual_status"] == current["accounting_identity_status"] == bc["status"] == "MET" else "NOT_MET",
                      physical="MET" if all(row["status"] == "MET" for row in current_report.values())
                      and profiles["status"] == "MET" else "NOT_MET")
    result.update(native=identity(args.native), source=first.data["source_identity"], table=identity(args.table),
                  order=args.order, observations=po.STATIONS, typed_diagnostic_candidate=args.candidate)
    atomic(args.out, result)
    print({k: result[k] for k in ("status", "equations", "physical") if k in result}, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_subparsers(dest="mode", required=True)
    fixture = modes.add_parser("prepare")
    for name in ("table", "geometry", "out"):
        fixture.add_argument("--"+name, required=True, type=Path)
    fixture.add_argument("--solve-seconds", required=True, type=float)
    fixture.add_argument("--seed", type=Path)
    fixture.add_argument("--refine-x", action="store_true")
    fixture.add_argument("--terminal-levels", type=int, default=0)
    fixture.add_argument("--wall-levels", type=int, default=0)
    fixture.add_argument("--x-band", type=float, nargs=2, action="append", default=[])
    fixture.add_argument("--refine-to-radius", type=float, default=0.)
    fixture.add_argument("--refinement-witness", type=Path)
    fixture.add_argument("--qx", type=int, default=8)
    fixture.add_argument("--qr", type=int, default=12)
    fixture.add_argument("--linear", type=float, default=1e-12)
    fixture.add_argument("--nonlinear", type=float, default=1e-10)
    fixture.add_argument("--preconditioner-mass", type=float, default=0.)
    fixture.add_argument("--velocity-preconditioner", choices=("FullGradientV1", "PositiveFrozenV1"),
                         default="FullGradientV1")
    fixture.add_argument("--krylov-dimension", type=int, default=80)
    audit = modes.add_parser("audit")
    for name in ("table", "native", "out"):
        audit.add_argument("--"+name, required=True, type=Path)
    audit.add_argument("--historical", type=Path)
    audit.add_argument("--candidate", action="store_true")
    audit.add_argument("--other", type=Path)
    audit.add_argument("--order", type=int, default=24)
    audit.add_argument("--seconds", type=float, required=True)
    audit.add_argument("--control", choices=("mesh", "quadrature", "algebraic", "time"), default="mesh")
    args = parser.parse_args()
    if args.mode == "prepare":
        prepare(args)
    else:
        require(args.other is not None or args.historical is not None, "physical audit needs the historical comparator")
        act(args)


if __name__ == "__main__":
    main()
