#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare immutable zero-time outlet cases and exercise the registered history path."""

import argparse
import hashlib
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
    "B2": ("outlet-b2", "M1cc", TRANSVERSE),
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
                   initial_source="", restart="", export_actions=False,
                   quadrature=(6, 8), tolerances=(1e-12, 1e-10)):
    fixtures.generate(session, table, session.parent/"state", mode="Physical", template=geometry,
                      dt=.05, steps=0, quadrature_x=quadrature[0], quadrature_r=quadrature[1],
                      linear=tolerances[0], nonlinear=tolerances[1], backend="FastFull3DV1",
                      solve="SteadyReference" if seconds else "Advance", steady_seconds=seconds,
                      initial_fields=initial_fields, initial_source_identity=initial_source,
                      restart=restart, outlet_policy=policy, check_outlet_actions=export_actions,
                      check_current_actions=export_actions)
    edit(session, {"WJReferencePolicy":"ExactInnerFiniteNuV1", "WJSteadyVelocityPC":"FullGradientV1"},
         {"WJSteadyKrylovDimension":64})


def tapered_guess(state, output, original):
    # A direct change of the endpoint W coefficient can create a steep
    # derivative in the last graded column. Only the starting W polynomial
    # is tapered. The native recipe still rebuilds all lifts and histories.
    n = state["sizes"][0]
    before = np.asarray(state["coefficients"]["flow"])
    transverse = before[n:2*n].copy()
    axial_nodes = np.array([0.,1.,1./3.,2./3.])
    bubble = axial_nodes*(1.-axial_nodes)
    modes = np.array([1.-axial_nodes,axial_nodes,bubble,(4.*axial_nodes-2.)*bubble]).T
    inverse_bubbles = np.linalg.inv(modes[2:,2:].T)
    left, right = 14., 15.
    claims = {}
    for e in state["elements"]:
        xl,xr = e["bounds"][:2]
        if xr <= left:
            continue
        x = xl+(xr-xl)*axial_nodes
        # Quintic matches values and the first two derivatives at both ends.
        t = np.clip((right-x)/(right-left),0.,1.)
        taper = t**3*(10.-15.*t+6.*t*t)
        ids = np.asarray(e["u"]).reshape(4,4)
        signs = np.asarray(e["us"]).reshape(4,4)
        local = before[n+ids]*signs
        values = (local@modes.T)*taper
        ends = values[:,:2]
        bubbles = (values[:,2:]-ends@modes[2:,:2].T)@inverse_bubbles
        changed = np.column_stack((ends,bubbles))*signs
        for dof,value in zip(ids.flat,changed.flat):
            require(abs(claims.setdefault(dof,value)-value)<=1e-12,"tapered W has disjoint signed traces")
    for dof,value in claims.items():
        transverse[dof] = value
    flow = before.copy()
    flow[n:2*n] = transverse
    coefficients = dict(state["coefficients"],flow=flow.tolist())
    digest = hashlib.sha256()
    for field in ("flow","k","epsilon"):
        digest.update(np.asarray(coefficients[field],dtype="<f8").tobytes())
    packet = {key:state[key] for key in ("schema","step","time","basis","model","source_identity",
              "reference_identity","mms_variant","manufactured","gauge","geometry","sizes","elements")}
    packet.update(coefficients=coefficients,state_sha256=digest.hexdigest())
    atomic(output,dict(schema="WJH27SteadyCandidateV1",qualified_reference=False,
        purpose="Smooth W transfer into B; original fields and guards are retained in the governing solve",
        source_packet=original, taper_band_code=[left,right],
        tapered_W_dofs=len(claims), maximum_W_coefficient_change=float(np.max(abs(transverse-before[n:2*n]))),
        state=packet))
    return packet


def prepare(root, case, seconds, continuation, guess):
    baseline, seals, table, runtime = get_context(root)
    subdir, mesh, policy = CASES[case]
    parent = root/"science"/subdir
    out = parent.with_name(f"{subdir}-continue-{continuation}") if continuation else parent
    out.mkdir()
    plan = None
    if mesh == "M1cc":
        path = root/"evidence/outlet-second-collar-plan.json"
        plan = identity(path)
        coords = read(path)
        if continuation:
            geometry = verified(read(parent/"inputs.json")["geometry"])
        else:
            geometry = out/"geometry.xml"
            fixtures.mesh(geometry,coords["x_edges"],coords["physical_z_edges"],quadrature=12)
        xs,zs,_ = mesh_axes(geometry)
        require(xs == coords["x_edges"] and zs == coords["physical_z_edges"],
                "second collar bisection changed the frozen axes")
    elif mesh == "M1c":
        if case == "B1c" or continuation:
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
    if guess == "smooth-outlet-archive":
        require(case in ("B0","B1"),"archive W taper is restricted to the matched B diagnostics")
        original = baseline["cases"][mesh]["state"]
        source_packet, _ = po.load_native(verified(original))
        require(source_packet["source_identity"] == baseline["source_identity"],"archive source differs")
        seed = out/"tapered-initial.json"
        seed_state = tapered_guess(source_packet,seed,original)
    elif guess == "b1c-after-worker-stop":
        require(case=="B2" and continuation>0,"fresh worker-stop retry requires a bounded B2 retry")
        previous = parent if continuation==1 else parent.with_name(f"{subdir}-continue-{continuation-1}")
        failed = read(root/"operations"/(previous.name+"-steady")/"receipt.json")
        require(failed["result"]["status"]=="FAIL"
                and failed["result"]["reason"] in ("HOST_WORKER_LIMIT","HOST_PRESSURE_REMAINS")
                and not (previous/"state.steady-candidate.json").exists(),
                "a fresh B2 retry requires the recorded worker stop before a native candidate")
        seed = accepted_directory(root,"B1c")/"state.result.json"
        seed_state, meta = po.load_native(seed)
        require(seed_state["source_identity"]==seals["source"] and seed_state["outlet_policy"]==policy
                and meta["steady_mode"] and meta["steady_solver"]["converged"],
                "worker-stop retry seed is not an accepted same-source B1c")
    elif continuation:
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
    elif case in ("B1c","B2") or guess == "accepted-b0":
        require(case in ("B1","B1c","B2"),"unexpected nested B guess")
        previous = {"B1":"B0","B1c":"B1","B2":"B1c"}[case]
        seed = accepted_directory(root,previous)/"state.result.json"
        seed_state, meta = po.load_native(seed)
        require(seed_state["source_identity"] == seals["source"] and seed_state["outlet_policy"] == policy
                and meta["steady_mode"] and meta["steady_solver"]["converged"],
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
        steady_time=0., initial_guess_method=guess, internal_cap_seconds=seconds, refinement_plan=plan,
        command=registered_command(runtime,session)))


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
    fixture.add_argument("--guess", choices=("inherited","smooth-outlet-archive","accepted-b0",
                                              "b1c-after-worker-stop"), default="inherited")
    args = parser.parse_args()
    if args.mode == "prepare":
        require(args.seconds > 0 and args.continuation >= 0, "invalid steady case cap")
        prepare(args.root, args.case, args.seconds, args.continuation, args.guess)
    else:
        {"seal":seal, "probe-history":check_history, "parity":parity}[args.mode](args.root, args.out)


if __name__ == "__main__":
    main()
