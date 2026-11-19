#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare immutable zero-time outlet cases and exercise the registered history path."""

import argparse
from pathlib import Path
import subprocess
import xml.etree.ElementTree as ET

import numpy as np

import coupled_fixtures as fixtures
import physical_oracle as po
from resources import atomic, identity
from task2_outlet_baseline import mesh_axes, read, source_identity, verified
from task2r_verify import edit, require


FULL = "FullTractionV1"
TRANSVERSE = "ZeroTransverseV1"
# Directory, frozen mesh source, outlet tensor law.
CASES = {
    "A0-replay": ("outlet-a0-replay", "A0", FULL),
    "A1-replay": ("outlet-a1-replay", "A1", FULL),
    "B0": ("outlet-b0", "A0", TRANSVERSE),
    "B1": ("outlet-b1", "A1", TRANSVERSE),
    "A1c": ("outlet-a1c", "M1c", FULL),
    "B1c": ("outlet-b1c", "M1c", TRANSVERSE),
}


def registered_command(runtime, session):
    return ["mpirun", "-np", "1", str(runtime), "--use-opt-file", str(session)+".opt", str(session)]


def seal(root, out):
    source, files = source_identity(root / "source-outlet-v1")
    runtime = root / "runtime-outlet-v1/IncNavierStokesSolver"
    require(source.encode() in runtime.read_bytes(), "compiled source identity absent from executable")
    atomic(out, dict(schema="WJH27OutletRuntimeV1", source=source, application=identity(runtime),
                     source_inventory=files, build=identity(root/"operations/outlet-native-build/receipt.json"),
                     source_stage=identity(root/"operations/outlet-source-stage/receipt.json")))


def get_context(root):
    baseline = read(root/"evidence/outlet-baseline.json")
    seal = read(root/"evidence/outlet-runtime.json")
    return baseline, seal, verified(baseline["table"]), verified(seal["application"])


def accepted_directory(root, case):
    name = CASES[case][0]
    science = root/"science"
    matches = [p for p in [science/name, *science.glob(name+"-continue-*")]
               if (p/"state.result.json").is_file()]
    require(len(matches) == 1, f"expected one accepted {case} result, found {len(matches)}")
    return matches[0]


def create_session(session, geometry, table, policy, *, seconds=0, initial_fields="",
                   initial_source="", restart="", export_actions=False):
    fixtures.generate(session, table, session.parent/"state", mode="Physical", template=geometry,
                      dt=.05, steps=0, quadrature_x=6, quadrature_r=8, linear=1e-12,
                      nonlinear=1e-10, backend="FastFull3DV1",
                      solve="SteadyReference" if seconds else "Advance", steady_seconds=seconds,
                      initial_fields=initial_fields, initial_source_identity=initial_source,
                      restart=restart, outlet_policy=policy, check_outlet_actions=export_actions,
                      check_current_actions=export_actions)
    edit(session, {"WJReferencePolicy":"ExactInnerFiniteNuV1", "WJSteadyVelocityPC":"FullGradientV1"},
         {"WJSteadyKrylovDimension":64})


def prepare(root, case, seconds, continuation):
    baseline, seals, table, runtime = get_context(root)
    subdir, mesh, policy = CASES[case]
    parent = root/"science"/subdir
    out = parent.with_name(f"{subdir}-continue-{continuation}") if continuation else parent
    out.mkdir()
    if mesh == "M1c":
        if case == "B1c":
            control = read(root/"science/outlet-a1c/inputs.json")
            geometry = verified(control["geometry"])
        else:
            geometry = out/"geometry.xml"
            coords = baseline["M1c"]
            fixtures.mesh(geometry, coords["x_edges"], coords["physical_z_edges"], quadrature=12)
        xs, zs, _ = mesh_axes(geometry)
        require(xs == baseline["M1c"]["x_edges"] and
                zs == baseline["M1c"]["physical_z_edges"], "collar-only bisection changed mesh axes")
    else:
        geometry = verified(baseline["cases"][mesh]["geometry"])
        xs, zs, _ = mesh_axes(geometry)
    if continuation:
        prior = parent.with_name(f"{subdir}-continue-{continuation-1}") if continuation > 1 else parent
        seed = prior/"state.steady-candidate.json"
        packet = read(seed)
        require(packet["schema"] == "WJH27SteadyCandidateV1" and not packet["qualified_reference"],
                "expected an explicitly typed unqualified steady continuation")
        prior_inputs = read(prior/"inputs.json")
        require(prior_inputs["case"] == case and prior_inputs["application"] == seals["application"]
                and prior_inputs["geometry"]["sha256"] == identity(geometry)["sha256"],
                "continuation input identity differs")
        seed_state = packet["state"]
        require(seed_state["source_identity"] == seals["source"] and seed_state["outlet_policy"] == policy,
                "continuation is from a different solver or outlet policy")
    elif case == "B1c":
        seed = accepted_directory(root,"B1")/"state.result.json"
        seed_state, _ = po.load_native(seed)
        require(seed_state["source_identity"] == seals["source"] and seed_state["outlet_policy"] == policy,
                "candidate collar guess is not from the same B application")
    else:
        seed = verified(baseline["cases"][mesh if mesh != "M1c" else "A1"]["state"])
        seed_state, _ = po.load_native(seed)
        require(seed_state["source_identity"] == baseline["source_identity"], "archive source differs")
    session = out/"run.xml"
    create_session(session, geometry, table, policy, seconds=seconds, initial_fields=seed,
                   initial_source=seed_state["source_identity"])
    atomic(out/"inputs.json", dict(schema="WJH27OutletCaseInputsV1", case=case,
        policy=policy, source=seals["source"], application=seals["application"],
        geometry=identity(geometry), session=identity(session), options=identity(Path(str(session)+".opt")),
        table=identity(table), initial_guess=identity(seed), initial_coefficient_sha256=seed_state["state_sha256"],
        cells=[len(xs)-1,len(zs)-1], modes=[4,3], quadrature=[6,8], tolerances=[1e-12,1e-10],
        steady_time=0., internal_cap_seconds=seconds, command=registered_command(runtime,session)))


def check_history(root, out):
    baseline, seals, table, runtime = get_context(root)
    out.mkdir()
    geometry = out/"geometry.xml"
    # The independent probe reduces the streamwise DOFs. Use retained radial
    # resolution to preserve native admissibility of the exact profile lift.
    fixtures.mesh(geometry, [-10., 0., 15.], baseline["cases"]["A1"]["physical_z_edges"], quadrature=12)
    attempts = {}

    def run(name, policy, seed=None, error=None, wrong_mask=False):
        location = out/name
        location.mkdir()
        session = location/"run.xml"
        history = out/seed/"state.final.bin" if seed else ""
        create_session(session, geometry, table, policy, restart=history, export_actions=True)
        if wrong_mask:
            xml = ET.parse(session)
            region = xml.find("CONDITIONS/BOUNDARYCONDITIONS/REGION[@REF='1']")
            region.find("D[@VAR='v']").tag = "N"
            xml.write(session, encoding="unicode")
        command = registered_command(runtime, session)
        with (location/"process.log").open("x") as log:
            status = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False).returncode
        result = location/"state.result.json"
        log = location/"process.log"
        record = dict(command=command, session=identity(session), options=identity(Path(str(session)+".opt")),
                      log=identity(log), expected_error=error, exit=status,
                      input_history=identity(history) if history else None)
        attempts[name] = record
        if error:
            require(status != 0 and not result.exists() and error in log.read_text(), f"{name} failed its expected rejection")
            return None
        require(status == 0, f"{name} did not run through the application: {log}")
        state, _ = po.load_native(result)
        require(state["source_identity"] == seals["source"] and state["outlet_policy"] == policy
                and state["step"] == state["time"] == 0. and not state["gauge"], "probe case binding differs")
        record.update(signature=state["signature"], digest=state["state_sha256"],
                      outlet_identity=state["outlet_identity"], native=identity(result),
                      history=identity(location/"state.final.bin"))
        if seed:
            require(attempts[seed]["signature"] == record["signature"]
                    and attempts[seed]["history"]["sha256"] == record["history"]["sha256"],
                    "same-policy cold history differs")
        return state

    a = run("A", FULL)
    b = run("B", TRANSVERSE)
    for key in ("signature", "outlet_identity"):
        require(a[key] != b[key], "outlet policies unexpectedly share an identity")
    run("A-to-A", FULL, "A")
    run("B-to-B", TRANSVERSE, "B")
    for policy, seed, name in ((TRANSVERSE, "A", "A-to-B"), (FULL, "B", "B-to-A")):
        run(name, policy, seed, "coupled checkpoint scientific identity differs")
    run("B-wrong-mask", TRANSVERSE, error="native component boundary mask disagrees", wrong_mask=True)
    atomic(out/"history-checks.json",dict(status="PASS", geometry=identity(geometry),
           application=seals["application"], attempts=attempts))


def parity(root, out):
    baseline, seals, _, _ = get_context(root)
    rows = {}
    for name in ("A0", "A1"):
        archive = baseline["cases"][name]
        first, _ = po.load_native(verified(archive["state"]))
        directory = root/"science"/CASES[name+"-replay"][0]
        replay = directory/"state.result.json"
        second, meta = po.load_native(replay)
        transfer = second["initial_field_transfer"]
        require(second["source_identity"] == seals["source"] and second["outlet_policy"] == FULL
                and second["boundary_identity"] == first["boundary_identity"]
                and second["time"] == second["step"] == first["time"] == first["step"] == 0.
                and meta["steady_mode"] and meta["steady_solver"]["converged"]
                and not transfer["projected"] and not transfer["source_was_candidate"]
                and transfer["source_state_sha256"] == archive["coefficient_sha256"]
                and transfer["packet_sha256"] == archive["state"]["sha256"],
                "A replay inputs, typed seed, or metadata differ")
        require(first["elements"] == second["elements"] and
                first["free_velocity_ids"] == second["free_velocity_ids"],
                "A replay coefficients are on different signed maps")
        scales = dict(flow=1., k=.008809015145515192, epsilon=.0003473286082198882/.15)
        diffs = {field:float(np.max(abs(np.asarray(first["coefficients"][field])-
                                         np.asarray(second["coefficients"][field]))))/scale
                 for field,scale in scales.items()}
        require(max(diffs.values()) <= 1e-10 and max(second["residuals"]) <= 1e-10,
                "replayed A differs from archived coefficients or current equations")
        rows[name] = dict(original=archive["state"], reproduced=identity(replay),
                          original_source=first["source_identity"], reproduced_source=second["source_identity"],
                          original_signature=first["signature"], reproduced_signature=second["signature"],
                          original_coefficient_sha256=first["state_sha256"],
                          reproduced_coefficient_sha256=second["state_sha256"],
                          maximum_scaled_coefficient_differences=diffs, original_residuals=first["residuals"],
                          fresh_native_residuals=second["residuals"],
                          residual_change_MAX=max(abs(x-y) for x,y in zip(first["residuals"],second["residuals"])),
                          fresh_advection_sha256=second["advection_sha256"], typed_transfer=transfer)
    atomic(out,dict(schema="WJH27OutletArchiveParityV1",status="PASS",cases=rows,
                    new_application=seals["application"], old_application=baseline["application"],
                    reference=baseline["table"]))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    modes = parser.add_subparsers(dest="mode", required=True)
    for name in ("seal", "probe-history", "parity"):
        modes.add_parser(name).add_argument("out", type=Path)
    fixture = modes.add_parser("prepare")
    fixture.add_argument("case", choices=CASES)
    fixture.add_argument("seconds", type=float)
    fixture.add_argument("--continuation", type=int, default=0)
    args = parser.parse_args()
    if args.mode == "prepare":
        require(args.seconds > 0 and args.continuation >= 0, "invalid steady case cap")
        prepare(args.root, args.case, args.seconds, args.continuation)
    else:
        {"seal":seal, "probe-history":check_history, "parity":parity}[args.mode](args.root, args.out)


if __name__ == "__main__":
    main()
