#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind and admit physical advances on the qualified Task 2R Bf mesh."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time
import xml.etree.ElementTree as ET

from resources import atomic, identity, update
from task2r_run import (GIB, balance, check_storage, choose_rank_cpu, require,
                        storage)


AUTHORIZATION = "0f979f00b03a83e76b4e06b49d9c0f03cbd96afd06d845fa94431215498fd631"
HANDOFF = "b1b7f70393b147652d9c1e4d0c83e0d18510b5890ffb6d82bb68b16fbac629c8"
BASE_SHA = "4b74e7e230bc25fe1e31bc32151228af8dfbeba6"
BINDING = "r3-bf-binding.json"
AMENDMENT = "r3-time-amendment"
OUTLET = "ZeroTransverseV1"
STEP_PRECONDITIONERS = ("CurrentDataDriftV1", "CurrentDataAMDOrderingV1")
INFO = dict(WJNativeMode="Physical", WJReferencePolicy="ExactInnerFiniteNuV1",
            WJConstitutiveBackend="FastFull3DV1", WJStressSplit="PositiveShearV1",
            WJOutletPolicy=OUTLET, WJCheckOutletActions="False",
            WJCheckCurrentActions="False", WJImplicitSolver="CoupledJFNKV1",
            WJStepPreconditioner="CurrentDataDriftV1")
RECIPE = dict(WJQuadratureX=6, WJQuadratureR=8, WJLinearTolerance=1e-12,
              WJCoupledTolerance=1e-10, WJDiagnosticOutletX=15, Time=0,
              Kinvis=6.423173803526447e-5, WJStepPCRelativeDrift=.02,
              WJStepPCMaximumUses=8)


def read(path):
    return json.loads(Path(path).read_text())


def verified(record):
    path = Path(record["path"])
    require(identity(path) == record, f"bound file differs: {path}")
    return path


def publish_event(path, event):
    if not path.exists():
        atomic(path, event)
    require(read(path) == event, f"stored accounting receipt differs: {path}")


def admit_scope(root, book, authorization, began):
    mandate = identity(authorization)
    handoff = identity(root / "evidence/outlet-campaign-final-handoff.json")
    require(mandate["sha256"] == AUTHORIZATION and handoff["sha256"] == HANDOFF,
            "R3 authority or Bf handoff differs")
    ledger = read(book)
    before = balance(ledger)
    label = "r3-scope"
    receipt = root / "evidence/r3-scope.json"
    if label in ledger["entries"]:
        event = ledger["entries"][label]
        require(event["authorization"] == mandate and event["handoff"] == handoff,
                "existing R3 scope binding differs")
    else:
        settlement = root / "operations/outlet-campaign-final-publish/receipt.json"
        prior = read(settlement)
        require(not ledger["closed"] and ledger["ceiling"] == 108000
                and ledger["protected_seconds"] == 3600 and before["reserved_seconds"] == 0
                and ledger["entries"]["outlet-campaign-final-publish"] == prior["result"]
                and ledger["entries"]["outlet-worker-amendment"]["worker_limit"] == 80
                and before["charged_seconds"] >= prior["balance"]["charged_seconds"],
                "R3 requires the reconciled outlet settlement")
        require(not receipt.exists(), "unexpected R3 scope receipt")
        for name in ("science-r3", "tmp-r3"):
            (root / name).mkdir(exist_ok=True)
        info = storage(root, advancing=True)
        check_storage(info)
        require(info["r3_science_bytes"] == 0, "R3 storage is occupied before admission")

        def record(value):
            require(value == ledger, "ledger changed under R3 scope lock")
            value["entries"][label] = event = dict(
                status="PASS", phase="accounting", seconds=time.monotonic() - began,
                authorization=mandate, handoff=handoff, settlement=identity(settlement),
                opening_ledger=identity(book), before=before, storage=info,
                science_cap_bytes=8 * GIB, combined_science_cap_bytes=16 * GIB,
                free_floor_bytes=15 * GIB, build_cap_bytes=8 * GIB, build_workers=4,
                build_rss_bytes=16 * GIB, named_science_maximum_rss_bytes=20 * GIB,
                worker_limit=80, science_ranks=1, threads=1)
            event["after"] = balance(value)

        event = update(book, record)["entries"][label]
    publish_event(receipt, event)


def amend_time(root, book, forecast_path, began):
    ledger = read(book)
    scope = ledger["entries"]["r3-scope"]
    forecast_file = identity(forecast_path)
    receipt = root / "evidence" / f"{AMENDMENT}.json"
    if AMENDMENT in ledger["entries"]:
        event = ledger["entries"][AMENDMENT]
        require(event["forecast"] == forecast_file and ledger["ceiling"] == event["ceiling"],
                "R3 time top-up was already applied to a different forecast")
    else:
        forecast = read(forecast_path)
        before = balance(ledger)
        require(not ledger["closed"] and ledger["ceiling"] == 108000 and
                ledger["protected_seconds"] == 3600 and before["reserved_seconds"] == 0
                and forecast["schema"] == "WJH27R3TimeNeedV1" and not receipt.exists(),
                "R3 top-up requires a settled campaign with a measured forecast")
        for anchor in forecast["measured_prefixes"]:
            prior = read(verified(anchor))
            require(prior["request"]["phase"] == "r3-advance"
                    and prior["result"] == ledger["entries"][prior["request"]["label"]],
                    "time forecast anchor is not a recorded R3 prefix")
        require(forecast["measured_prefixes"] and forecast["milestone"].strip(),
                "time top-up requires an observed prefix and named milestone")
        need = forecast["remaining_cost_seconds"]
        extra = forecast["additional_seconds"]
        require(math.isfinite(need) and 0 < extra <= 21600
                and 0 < need - before["nonprotected_remaining"] <= extra,
                "R3 time request exceeds authorization or has no measured need")

        def amend(value):
            require(value == ledger, "ledger changed under R3 top-up lock")
            value["ceiling"] += extra
            value["entries"][AMENDMENT] = event = dict(
                status="PASS", phase="accounting", seconds=time.monotonic() - began,
                authorization=scope["authorization"], forecast=forecast_file,
                opening_ledger=identity(book), before=before, prior_ceiling=108000,
                added_seconds=extra, ceiling=value["ceiling"], protected_seconds=3600)
            event["after"] = balance(value)

        event = update(book, amend)["entries"][AMENDMENT]
    publish_event(receipt, event)


def identity_of_state(path):
    import numpy as np

    state = json.loads(Path(path).read_text(),
                       parse_int=lambda token: -0.0 if token == "-0" else int(token))
    state = state.get("state", state)
    digest = hashlib.sha256()
    for name in ("flow", "k", "epsilon"):
        digest.update(np.asarray(state["coefficients"][name], dtype="<f8").tobytes())
    require(digest.hexdigest() == state["state_sha256"], "native coefficient digest differs")
    return state


def bind(root):
    from task2_outlet_baseline import mesh_axes

    scope = root / "evidence/r3-scope.json"
    handoff = read(verified(read(scope)["handoff"]))
    selected = handoff["verify_inputs"]["Bf"]
    record = selected["binding"]
    for key in ("state", "mesh", "executable", "inputs", "operation", "admission"):
        verified(record[key])
    for key in ("audit", "table", "session", "options"):
        verified(selected[key])
    x, z, _ = mesh_axes(record["mesh"]["path"])
    require(record["cells"] == [len(x)-1, len(z)-1] == [111, 92]
            and x[-1] == 15 and record["policy"] == OUTLET and record["quadrature"] == [6, 8]
            and record["tolerances"] == [1e-12, 1e-10], "Bf mesh or numerical recipe differs")
    state = identity_of_state(record["state"]["path"])
    require(state["state_sha256"] == record["coefficient_sha256"]
            and state["source_identity"] == record["source"]
            and state["outlet_policy"] == OUTLET
            and state["signature"] == record["signature"]
            and state["step"] == state["time"] == 0
            and record["source"].encode() in Path(record["executable"]["path"]).read_bytes(),
            "Bf coefficient, source, clock or model binding differs")
    anchor = read(verified(handoff["parity"]))
    atomic(root / "evidence" / BINDING, dict(
        schema="WJH27R3BfBindingV1", settled_revision=BASE_SHA,
        scope=identity(scope), handoff=read(scope)["handoff"], recipe=RECIPE, solver=INFO,
        inherited=selected, x15_parity=handoff["parity"], parity_replay=anchor["replay"],
        base_identity={key: state[key] for key in ("model", "reference_identity", "basis",
                                                  "outlet_policy", "state_sha256")},
        affinity=choose_rank_cpu()))


def prepare(root, label, role, dt, steps, cap, step_cap, rss_gib, peak_mib,
            science_stride, recovery_stride, basis_path=None, certificate=None, recovery=None,
            phase="r3-advance", native_seconds=None, initial_step=0, represented_proof=None,
            maximum_pause=None, step_preconditioner=STEP_PRECONDITIONERS[0]):
    import coupled_fixtures as fixtures
    from task2r_verify import edit

    require(re.fullmatch(r"r3-[a-z0-9]+(?:-[a-z0-9]+)*", label), "invalid R3 label")
    bound = root / "evidence" / BINDING
    base = read(bound)
    selected = base["inherited"]
    old = selected["binding"]
    prepared = read(basis_path) if basis_path else None
    native = prepared or dict(application=old["executable"], source=old["source"], state=old["state"])
    runtime = verified(native["application"])
    geometry = verified(old["mesh"])
    table = verified(selected["table"])
    require(role in ("quiet", "startup", "perturbed", "restart", "same-mesh-parity"),
            "invalid physical start")
    seed = verified(native["state"]) if role in ("quiet", "perturbed", "same-mesh-parity") else ""
    source = native.get("seed_source", native["source"]) if seed else ""
    parity = role == "same-mesh-parity"
    initial_only = phase == "r3-initial"
    require(maximum_pause is None or 0 < maximum_pause < cap,
            "pause allowance must fit inside the charged operation cap")
    require((phase == "r3-parity") == parity and dt > 0
            and (steps == 0 if parity or initial_only else steps > 0)
            and (role == "perturbed") == bool(certificate) and (not initial_only or role == "perturbed")
            and (role == "restart") == bool(recovery) and (phase == "r3-restart") == bool(recovery)
            and (initial_step >= 2 if recovery else initial_step == 0)
            and (parity or not seed or source == native["source"]), "invalid R3 start identity")
    if certificate and not initial_only:
        require(represented_proof is not None, "advancement requires independently certified Bf coefficients")
        proof = read(represented_proof)
        require(proof["schema"] == "WJH27R3PhysicalBfPerturbationAuditV1" and proof["status"] == "MET"
                and proof["base"] == native["state"] and proof["certificate"] == identity(certificate),
                "represented perturbation proof differs from the advancing inputs")
    out = root / "science-r3" / label
    out.mkdir()
    session = out / "run.xml"
    fixtures.generate(session, table, out/"state", mode="Physical", template=geometry,
                      dt=dt, steps=steps, quadrature_x=6, quadrature_r=8, backend="FastFull3DV1",
                      recovery=recovery_stride, science=science_stride, initial_fields=seed,
                      initial_source_identity=source, restart=recovery or "", outlet_policy=OUTLET,
                      solve="SteadyReference" if parity else "Advance",
                      steady_seconds=native_seconds if parity else 0)
    require(step_preconditioner in STEP_PRECONDITIONERS,"unrecognized tested advancing preconditioner")
    info = dict(INFO, WJStepPreconditioner=step_preconditioner)
    parameters = dict(RECIPE, WJStepNewtonSeconds=step_cap)
    if parity:
        info["WJSteadyVelocityPC"] = "FullGradientV1"
        parameters["WJSteadyKrylovDimension"] = 64
    else:
        info["WJTerminalOutput"] = "True"
        parameters.update(WJTerminalBegin=dt*initial_step, WJTerminalEnd=dt*(initial_step+max(steps, 1)))
    if certificate:
        info.update(WJPerturbationCertificate=str(certificate),
                    WJPerturbationCertificateSHA256=identity(certificate)["sha256"])
    edit(session, info, parameters)
    opt = Path(str(session) + ".opt")
    command = ["mpirun", "-np", "1", str(runtime), "--use-opt-file", str(opt), str(session)]
    inputs = [identity(path) for path in
              (geometry, table, seed, certificate, recovery, basis_path, represented_proof) if path]
    plan = dict(schema="WJH27R3AdmissionV1", label=label, phase=phase, role=role,
                seconds=cap, peak_new_bytes=peak_mib * 1024 ** 2, rss_limit_bytes=rss_gib * GIB,
                maximum_paused_seconds=maximum_pause,
                step_preconditioner=step_preconditioner,
                cpu=base["affinity"]["cpu"], bf_binding=identity(bound), command=command,
                timestep=dt, initial_step=initial_step, requested_steps=steps, native_step_cap_seconds=step_cap,
                science_stride=science_stride, recovery_stride=recovery_stride,
                source=native["source"], application=identity(runtime), geometry=identity(geometry),
                table=identity(table), initial_fields=identity(seed) if seed else None,
                restart=identity(recovery) if recovery else None, inputs=inputs,
                certificate=identity(certificate) if certificate else None,
                represented_proof=identity(represented_proof) if represented_proof else None,
                resource_anchors=[old["operation"], base["x15_parity"]],
                fixture=identity(session), options=identity(opt))
    atomic(out / "admission.json", plan)


def admission(root, label, phase, command, cap):
    out = root / "science-r3" / label
    path = out / "admission.json"
    plan = read(path)
    base = read(verified(plan["bf_binding"]))
    require(plan["schema"] == "WJH27R3AdmissionV1" and plan["phase"] == phase
            and plan["label"] == label and plan["command"] == command and plan["seconds"] == cap
            and plan["cpu"] == base["affinity"]["cpu"]
            and 0 < plan["rss_limit_bytes"] <= 20 * GIB
            and 0 < plan["peak_new_bytes"] <= 8 * GIB, "R3 command or resource forecast differs")
    for record in plan["inputs"]:
        verified(record)
    if phase == "r3-audit":
        require(command[:2] == ["python3", "-B"] and Path(command[2]).is_file()
                and Path(command[2]).is_relative_to(root / "tools"),
                "R3 audit requires a recorded independent analyzer")
    else:
        session, opt = verified(plan["fixture"]), verified(plan["options"])
        runtime = verified(plan["application"])
        require(command == ["mpirun", "-np", "1", str(runtime), "--use-opt-file",
                            str(opt), str(session)] and session == out / "run.xml",
                "R3 requires the named standard-main one-rank command")
        require(plan["geometry"] == base["inherited"]["binding"]["mesh"],
                "R3 solve has departed from the Bf geometry")
        xml = ET.parse(session)
        info = {i.get("PROPERTY"): i.get("VALUE") for i in xml.findall("CONDITIONS/SOLVERINFO/I")}
        parameters = dict(p.text.replace(" ", "").split("=") for p in xml.findall("CONDITIONS/PARAMETERS/P"))
        if phase == "r3-stock":
            stock_switch = "False" if plan["stock_svv"] == "Off" else plan["stock_svv"]
            require(info["SolverType"] == "VCSWeakPressure" and info["Driver"] == "Standard"
                    and int(parameters["NumSteps"]) == plan["requested_steps"] > 0
                    and float(parameters["TimeStep"]) == plan["timestep"]
                    and float(parameters["Kinvis"]) == RECIPE["Kinvis"]
                    and info.get("SpectralVanishingViscosity") == stock_switch,
                    "stock admission requires the conventional weak-pressure solver")
        else:
            step_pc = plan.get("step_preconditioner",STEP_PRECONDITIONERS[0])
            solver_info = dict(INFO,WJStepPreconditioner=step_pc)
            require(info["SolverType"] == "WJNativeCoupledV2" and info["Driver"] == "Standard"
                    and step_pc in STEP_PRECONDITIONERS and all(info.get(k) == v for k, v in solver_info.items())
                    and all(float(parameters[k]) == v for k, v in RECIPE.items())
                    and info["WJReferenceTable"] == plan["table"]["path"]
                    and plan["table"] == base["inherited"]["table"]
                    and info.get("WJInitialFields") == (plan["initial_fields"] or {}).get("path")
                    and info.get("WJRestart") == (plan["restart"] or {}).get("path"),
                    "R3 changed the frozen coupled recipe or start")
            parity = phase == "r3-parity"
            zero = parity or phase == "r3-initial"
            certificate = plan.get("certificate")
            require(info["WJNativeSolve"] == ("SteadyReference" if parity else "Advance")
                    and int(parameters["NumSteps"]) == plan["requested_steps"]
                    and float(parameters["TimeStep"]) == plan["timestep"]
                    and (plan["requested_steps"] == 0 if zero else plan["requested_steps"] > 0)
                    and info.get("WJPerturbationCertificate") == (certificate or {}).get("path")
                    and info.get("WJPerturbationCertificateSHA256") == (certificate or {}).get("sha256")
                    and (phase != "r3-advance" or certificate is None or plan.get("represented_proof"))
                    and (phase == "r3-restart") == bool(plan["restart"]),
                    "R3 timeline, start history or zero-time scope differs")
    return dict(plan, admission=identity(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    sub = parser.add_subparsers(dest="action", required=True)
    sub.add_parser("bind")
    case = sub.add_parser("prepare")
    case.add_argument("label")
    case.add_argument("role")
    case.add_argument("--dt", type=float, required=True)
    case.add_argument("--steps", type=int, required=True)
    case.add_argument("--cap", type=int, required=True)
    case.add_argument("--step-cap", type=int, required=True)
    case.add_argument("--rss-gib", type=int, default=16)
    case.add_argument("--peak-mib", type=int, default=550)
    case.add_argument("--science-stride", type=int, default=1)
    case.add_argument("--recovery-stride", type=int, default=2)
    case.add_argument("--basis-path", type=Path)
    case.add_argument("--certificate", type=Path)
    case.add_argument("--recovery", type=Path)
    case.add_argument("--phase", default="r3-advance")
    case.add_argument("--native-seconds", type=float)
    case.add_argument("--initial-step", type=int, default=0)
    case.add_argument("--represented-proof", type=Path)
    case.add_argument("--maximum-pause", type=int)
    case.add_argument("--step-preconditioner", choices=STEP_PRECONDITIONERS, default=STEP_PRECONDITIONERS[0])
    args = vars(parser.parse_args())
    action = args.pop("action")
    if action == "bind":
        bind(args["root"])
    else:
        prepare(**args)


if __name__ == "__main__":
    main()
