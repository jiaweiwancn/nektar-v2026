#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Prepare, admit and independently check the flat outlet B qualification ladder."""

import argparse
import json
import math
from pathlib import Path
import time

import coupled_fixtures as fixtures
import numpy as np
import physical_oracle as po

from resources import atomic, identity
from task2_outlet_analysis import INTERIOR, audit_state, bind_elements, open_case, pair
from task2_outlet_baseline import mesh_axes, read, source_identity, verified
from task2_outlet_case import TRANSVERSE, create_session, get_context, registered_command
from task2r_run import GIB, QUAL_CASES, balance, memory_headroom, memory_snapshot, require, workers
from task2r_verify import edit


CP4 = "aabffc99f67a2dc941d02b356c32e1c496ad9cfb"
LAUNCH_HOLD = "742e63ec01096965983316838b06a31ced766f2f33c590bf64380045f9e558c5"
QUADRATURE = [6, 8]
TOLERANCES = [1e-12, 1e-10]
HIGHER_QUADRATURE = [8, 12]
TIGHTER_TOLERANCES = [3e-13, 3e-12]
SHIFT = .0875
BLEND = (.06, .1125)
PLAN = "outlet-qual-plan.json"
PUBLICATION_SECONDS = 300
GRACE = 7


def slug(case):
    return "outlet-qual-" + QUAL_CASES[case]


def planned_attempt(native, watch, iterations):
    return dict(native_seconds=native, watch_seconds=watch, newton_iterations=iterations,
                preparation_seconds=18)


def case_dir(root, case, continuation=0):
    stem = slug(case)
    return root / "science" / (f"{stem}-continue-{continuation}" if continuation else stem)


def accepted_path(root, case):
    parent = case_dir(root, case)
    matches = [p for p in [parent, *parent.parent.glob(parent.name + "-continue-*")]
               if (p / "state.result.json").is_file()]
    require(len(matches) == 1, f"{case} must have one accepted native output")
    return matches[0] / "state.result.json"


def split_band(edges, lo, hi):
    extra = [(a + b) / 2 for a, b in zip(edges, edges[1:])
             if lo - 1e-12 <= a - SHIFT and b - SHIFT <= hi + 1e-12]
    return sorted(edges + extra), [z - SHIFT for z in extra]


def pass_audit(root, case):
    path = root / "evidence" / f"outlet-{case}-audit.json"
    audit = read(path)
    require(all(status == "MET" for status in audit["status"].values()),
            f"{case} did not pass its independent gates")
    return audit, identity(path)


def pin(root, out, provenance):
    evidence = root / "evidence"
    baseline, runtime, table, application = get_context(root)
    hold_path = evidence / "task3-launch-hold.json"
    hold = read(hold_path)
    require(identity(hold_path)["sha256"] == LAUNCH_HOLD and hold["TASK3_READY"] == "NO",
            "CP4 launch hold differs")
    source, inventory = source_identity(root / "source-outlet-v1")
    require(runtime["source"] == hold["native_source"] == source
            and runtime["source_inventory"] == inventory
            and hold["application"] == identity(application)
            and hold["reference"] == identity(table), "CP4 runtime or exact-inner reference differs")
    bound = {}
    for case in ("B1", "B1c", "B2"):
        audit, digest = pass_audit(root, case)
        witness = audit["provenance"]
        path = verified(witness["state"])
        geometry = verified(witness["mesh"])
        state, meta = po.load_native(path)
        cells = bind_elements(geometry, state)
        require(state["step"] == state["time"] == 0.
                and meta["steady_solver"]["converged"] and state["outlet_policy"] == TRANSVERSE
                and state["source_identity"] == source and state["inner_k_transport"] == "ScalarCdV2"
                and state["reference_identity"] == hold["reference"]["sha256"]
                and state["state_sha256"] == witness["coefficient_sha256"] and cells == witness["cells"]
                and [state["quadratureX"], state["quadratureR"]] == QUADRATURE
                and [state["linear_tolerance"], state["nonlinear_tolerance"]] == TOLERANCES,
                f"{case} differs from its archived scientific identity")
        bound[case] = dict(state=witness["state"], mesh=witness["mesh"], audit=digest,
                           coefficients=state["state_sha256"], signature=state["signature"], cells=cells)
    require(bound["B2"]["state"] == hold["scientific"]["B"]["latest_state"]
            and bound["B2"]["mesh"] == hold["scientific"]["B"]["latest_mesh"],
            "B2 does not bind the exact CP4 launch hold")
    cp4_pair = verified(hold["scientific"]["B"]["second_collar_comparison"])
    archived = read(cp4_pair)["comparison"]
    require(archived["first"]["state"] == bound["B1c"]["state"]
            and archived["second"]["state"] == bound["B2"]["state"]
            and len(archived["rows"]) == 690 and len(archived["failures"]) == 122,
            "CP4 comparison states or row inventory differ")

    xs, zs, _ = mesh_axes(Path(bound["B2"]["mesh"]["path"]))
    lo, _, hi = baseline["M1c"]["original_collar_owner_r"]
    full, full_r = split_band(zs, lo, hi)
    lower, lower_r = split_band(zs, lo, .1)
    require((len(xs), len(zs), len(full_r), len(lower_r)) == (112, 70, 8, 4),
            "expected the frozen B2 and both B3 nested designs")
    b2_labels = ["outlet-b2-continue-3-steady", "outlet-b2-continue-4-steady"]
    baseline_ops = [root / "operations" / label / "receipt.json" for label in b2_labels]
    result = [read(p)["result"] for p in baseline_ops]
    base_peak = max(r["peak_tree_rss_bytes"] for r in result)
    base_wall = sum(r["seconds"] for r in result)
    b3 = {}
    for label, edges, inserted in (("full", full, full_r), ("lower", lower, lower_r)):
        ratio = (len(edges) - 1) / (len(zs) - 1)
        b3[label] = dict(x_edges=xs, physical_z_edges=edges, added_r=inserted,
                         cells=[len(xs)-1, len(edges)-1], reference_case="B2",
                         support_r=[lo, hi if label=="full" else .1],
                         new_comparison_sequence=label!="full",
                         forecast_peak_tree_rss_bytes=math.ceil(base_peak * ratio * 1.10),
                         scaled_native_wall_seconds=math.ceil(base_wall * ratio))
    memory = memory_snapshot()
    live = balance(read(root / "accounting.json"))
    attempts = {case_dir(root, "B3", i).name: planned_attempt(inner, wall, iterations)
                for i, (inner, wall, iterations) in enumerate(((240, 355, 1), (740, 855, 50), (460, 575, 50)))}
    b3_cost = sum(item["watch_seconds"] + item["preparation_seconds"] + 2 * GRACE
                  for item in attempts.values())
    stage_caps = dict(cp0_publication=62, archived_replay=147, cp1_admission=47,
                      B3_preparation_native_and_continuation=b3_cost, B3_independent_audit=207,
                      B3_comparison=132, cp1_publication=62, final_publication=PUBLICATION_SECONDS)
    first_stage = sum(stage_caps.values())
    later = dict(broader_radial_witness=2200, bounded_streamwise_witness=1750,
                 higher_quadrature=1600, tighter_algebraic=1200, evidence_analysis=240,
                 conditional_extent=2464-PUBLICATION_SECONDS)
    # The extension and additional witness branches get their own admissions once their inputs exist.
    chosen = "full" if b3["full"]["forecast_peak_tree_rss_bytes"] <= 12*GIB else "lower"
    require(first_stage <= live["nonprotected_remaining"]
            and memory_headroom(memory) >= 14*GIB
            and b3[chosen]["forecast_peak_tree_rss_bytes"] <= 12*GIB, "B3 lacks full-stage resources")
    atomic(out, dict(schema="WJH27FlatQualificationPlanV1", prior_commit=CP4,
        old_launch_hold=identity(hold_path), cp4_pair=identity(cp4_pair), references=read(provenance),
        source=source, application=runtime["application"], table=hold["reference"], cases=bound,
        quadrature=QUADRATURE, tolerances=TOLERANCES, meshes=b3, chosen_B3=chosen, attempts=attempts,
        numerical_limits=po.NUMERICAL, higher_quadrature=HIGHER_QUADRATURE,
        tighter_tolerances=TIGHTER_TOLERANCES, both_control_inner_k_limit=1e-4, blend_r=list(BLEND),
        range_to_inspect_r=[.06875, .1125], experimental_x=list(INTERIOR),
        complete_priority_forecast=dict(first_stage=stage_caps, subsequent_estimates=later,
            staged_first_seconds=first_stage, all_if_supported_seconds=first_stage+sum(later.values()),
            native_cost_anchors=[identity(p) for p in baseline_ops]),
        extent_frozen=identity(evidence / "outlet-location-eligibility.json"),
        amendment=identity(evidence / "outlet-qual-amendment.json"),
        memory_at_pin=memory, workers_at_pin=workers(), live_with_pin_reservation=live))
    print("pinned", chosen, b3[chosen]["cells"], "staged forecast", first_stage, flush=True)


def admit_b3(root, out):
    path = root / "evidence" / PLAN
    plan = read(path)
    mesh = plan["meshes"][plan["chosen_B3"]]
    estimate = mesh["forecast_peak_tree_rss_bytes"]
    limit = 12*GIB if estimate > 10*GIB else 10*GIB
    live = balance(read(root / "accounting.json"))
    memory = memory_snapshot()
    caps = plan["complete_priority_forecast"]["first_stage"]
    pending = sum(caps[key] for key in ("B3_preparation_native_and_continuation", "B3_independent_audit",
                                       "B3_comparison", "cp1_publication", "final_publication"))
    require(pending <= min(live["nonprotected_remaining"], live["tranche_remaining"])
            and memory_headroom(memory) >= limit + 2*GIB and estimate <= limit,
            "live B3 reservation does not fit")
    atomic(out, dict(schema="WJH27FlatQualificationAdmissionV1", status="ADMITTED", case="B3",
        admission_plan=identity(path), source=plan["source"], application=plan["application"],
        prior_state=plan["cases"]["B2"]["state"], prior_mesh=plan["cases"]["B2"]["mesh"],
        quadrature=QUADRATURE, tolerances=TOLERANCES, **mesh, attempts=plan["attempts"],
        required_seconds=pending, reservation_caps=caps, live=live, worker_census=workers(),
        memory=memory, rss_limit_bytes=limit))
    print("admitted B3",mesh["cells"],estimate,flush=True)


def forecast_b3_recovery(root, out, start):
    evidence = root/"evidence"
    suffix = {2:"",4:"-1",6:"-2"}[start]
    first = evidence/f"outlet-qual-b3-admission{suffix}.json"
    original = read(first)
    previous = root/"operations"/(case_dir(root,"B3",start-1).name+"-steady")/"receipt.json"
    stop = read(previous)["result"]
    native_timeout = start==6
    packet = case_dir(root,"B3",start-1 if native_timeout else 0)/"state.steady-candidate.json"
    steady = read(packet.with_name("state.steady.json"))
    allowed = (stop["reason"] is None and steady["reason"].startswith("DEADLINE")
               if native_timeout else retry_before_checkpoint(stop))
    require(allowed and steady["best_from_accepted_iterate"] and steady["accepted_newton_steps"] > 0
            and packet.is_file(), "recovery forecast needs the sealed B3 iterate and bounded interruption")
    schedule = [(470,590,1),(560,680,50)]
    if start==2:
        schedule.insert(0,(380,500,1))
    elif native_timeout:
        require(TOLERANCES[1] < steady["true_norm"] <= 2*TOLERANCES[1],
                "the short recovery requires an iterate within twice the unchanged acceptance norm")
        schedule = [(365,480,50)]
    attempts = {case_dir(root,"B3",i).name:planned_attempt(inner,wall,count)
                for i,(inner,wall,count) in enumerate(schedule,start=start)}
    caps = dict(admission=27, preparation_and_native=sum(item["preparation_seconds"]+
                item["watch_seconds"]+2*GRACE for item in attempts.values()),
                audit=207,comparison=132,blend_probe=31,checkpoint_evidence=62,
                final_publication=PUBLICATION_SECONDS)
    live = balance(read(root/"accounting.json"))
    memory = memory_snapshot()
    required = sum(caps.values())
    blockers = [label for label, fits in (
        ("wall_time", required<=live["nonprotected_remaining"]),
        ("live_job_headroom", memory_headroom(memory)>=original["rss_limit_bytes"]+2*GIB)
        ) if not fits]
    atomic(out,dict(original,schema="WJH27FlatQualificationForecastV1",
        status="FITS" if not blockers else "DEFERRED",blockers=blockers,
        previous_admission=identity(first),interruption=identity(previous),
        recovery_checkpoint=identity(packet),recovery_residual=steady["true_norm"],
        recovery_first_label=case_dir(root,"B3",start).name,
        attempts=attempts,required_seconds=required,reservation_caps=caps,
        memory=memory,live=live,audit_cap_seconds=200,compare_cap_seconds=125))
    print("B3 recovery budget",required,"prior stop seconds",stop["seconds"],blockers,flush=True)


def retry_before_checkpoint(result):
    return (result["reason"] in ("HOST_WORKER_LIMIT","HOST_PRESSURE_REMAINS")
            or (result["reason"] == "native rank affinity differs" and result["seconds"] < 15.))


def recovery_packet(root, case, continuation):
    for index in range(continuation-1,-1,-1):
        previous = case_dir(root,case,index)
        receipt = read(root/"operations"/(previous.name+"-steady")/"receipt.json")
        require(receipt["result"]["status"] == "FAIL", "retry expected a settled failed prior segment")
        if (previous/"state.steady-candidate.json").is_file():
            return previous
        require(retry_before_checkpoint(receipt["result"]),
                "a sealed iterate or a recorded launch or worker interruption is required for retry")
    return None


def prepare(root, case, continuation, admission_path):
    context = get_context(root)
    _, runtime, table, application = context
    specification = read(admission_path)
    require(specification["status"] == "ADMITTED" and specification["case"] == case
            and specification["source"] == runtime["source"], "qualification case lacks a bound admission")
    target = case_dir(root, case, continuation)
    attempt = specification["attempts"][target.name]
    xs, zs = specification["x_edges"], specification["physical_z_edges"]
    initial = verified(specification["prior_state"])
    parent_mesh = verified(specification["prior_mesh"])
    compiled = runtime["source"]

    previous = recovery_packet(root,case,continuation)
    if previous is not None:
        binding = read(previous / "inputs.json")
        candidate = previous / "state.steady-candidate.json"
        packet = read(candidate)
        prior = read(previous / "state.steady.json")
        require(packet["schema"] == "WJH27SteadyCandidateV1" and not packet["qualified_reference"]
                and prior["best_from_accepted_iterate"] and prior["accepted_newton_steps"] > 0
                and prior["reason"].startswith(("DEADLINE", "NEWTON_NOT_CONVERGED"))
                and binding["quadrature"] == specification["quadrature"]
                and binding["tolerances"] == specification["tolerances"]
                and binding["x_edges"] == xs and binding["physical_z_edges"] == zs,
                "native packet lacks a matching, accepted zero-time Newton checkpoint")
        initial = candidate
        state = packet["state"]
        parent_mesh = verified(binding["geometry"])
    else:
        state, _ = po.load_native(initial)
    if target.name == specification.get("recovery_first_label"):
        require(identity(initial) == specification["recovery_checkpoint"],
                "recovery seed differs from the frozen B3 Newton iterate")
    require(state["step"] == state["time"] == 0. and state["source_identity"] == compiled
            and state["outlet_policy"] == TRANSVERSE, "wrong scope for the typed qualification seed")
    target.mkdir()
    x0, z0, _ = mesh_axes(parent_mesh)
    geometry = parent_mesh
    if (x0, z0) != (xs, zs):
        geometry = target / "geometry.xml"
        fixtures.mesh(geometry, xs, zs, quadrature=12)
        x, z, _ = mesh_axes(geometry)
        require(x == xs and z == zs, "session refinement differs from the admitted axes")
    session = target / "run.xml"
    create_session(session, geometry, table, TRANSVERSE, seconds=attempt["native_seconds"],
                   initial_fields=initial, initial_source=compiled, quadrature=specification["quadrature"],
                   tolerances=specification["tolerances"])
    edit(session, {}, {"WJSteadyNewtonIterations": attempt["newton_iterations"]})
    atomic(target / "inputs.json", dict(schema="WJH27FlatQualificationInputsV1", case=case,
        policy=TRANSVERSE, source=compiled, admission=identity(admission_path), application=identity(application),
        geometry=identity(geometry), x_edges=xs, physical_z_edges=zs,
        session=identity(session), options=identity(Path(str(session)+".opt")), table=identity(table),
        initial_guess=identity(initial), initial_coefficient_sha256=state["state_sha256"], modes=[4,3],
        cells=[len(xs)-1, len(zs)-1], quadrature=specification["quadrature"],
        tolerances=specification["tolerances"], steady_time=0., command=registered_command(application,session)))


def load_case(root, name, context, ref, order, deadline):
    if name not in QUAL_CASES:
        return open_case(root, name, context, ref, order, deadline)
    _, runtime, table, application = context
    path = accepted_path(root, name)
    manifest_path = path.parent / "inputs.json"
    manifest = read(manifest_path)
    forecast = read(verified(manifest["admission"]))
    receipt_path = root / "operations" / (path.parent.name+"-steady") / "receipt.json"
    receipt = read(receipt_path)
    witness = identity(path)
    geometry = verified(manifest["geometry"])
    x, z, _ = mesh_axes(geometry)
    require(manifest["case"] == name == forecast["case"] and manifest["source"] == runtime["source"]
            and manifest["application"] == identity(application) and manifest["table"] == identity(table)
            and manifest["quadrature"] == forecast["quadrature"] and manifest["tolerances"] == forecast["tolerances"]
            and [x,z] == [forecast["x_edges"],forecast["physical_z_edges"]]
            and receipt["request"]["command"] == manifest["command"]
            and all(manifest[key] in receipt["request"]["input_files"] for key in ("application","session","options"))
            and receipt["result"]["status"] == "PASS" and witness in receipt["result"]["output_files"],
            "accepted qualification evidence does not bind its manifest and RSS-watchdog receipt")
    for key in ("session","options","initial_guess"):
        verified(manifest[key])
    data, meta = po.load_native(path)
    transfer = data["initial_field_transfer"]
    require(data["source_identity"] == manifest["source"]
            and data["reference_identity"] == ref.identity and data["outlet_policy"] == TRANSVERSE
            and data["basis"] == "Modified_A" and data["inner_k_transport"] == "ScalarCdV2"
            and data["constitutive_backend"] == "FastFull3DV1" and data["stress_split"] == "PositiveShearV1"
            and data["step"] == data["time"] == 0. and meta["steady_mode"] and meta["steady_solver"]["converged"]
            and max(data["residuals"]) <= manifest["tolerances"][1]
            and [data["quadratureX"],data["quadratureR"]] == manifest["quadrature"]
            and [data["linear_tolerance"],data["nonlinear_tolerance"]] == manifest["tolerances"]
            and transfer["packet_sha256"] == manifest["initial_guess"]["sha256"]
            and transfer["source_state_sha256"] == manifest["initial_coefficient_sha256"],
            "qualified native field settings or accepted seed provenance differ")
    cells = bind_elements(geometry, data)
    state = po.PhysicalOracle(data,ref,order,deadline)
    return state, dict(state=witness, mesh=manifest["geometry"], source=runtime["source"],
        executable=manifest["application"], inputs=identity(manifest_path), operation=identity(receipt_path),
        admission=manifest["admission"], coefficient_sha256=data["state_sha256"], signature=data["signature"],
        outlet_identity=state.outlet_identity, policy=TRANSVERSE, cells=cells,
        quadrature=manifest["quadrature"], tolerances=manifest["tolerances"])


def open_pair(root, base, fine, seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    deadline = time.monotonic()+seconds
    a, pa = load_case(root,base,context,ref,20,deadline)
    b, pb = load_case(root,fine,context,ref,20,deadline)
    return a, b, pair(a,b,(pa,pb))


def raw_owner(state, owner, x, r):
    raw = state.sample(owner,np.array([x]),np.array([r]))
    G = raw["G"][0]
    return dict(owner=owner["geometry"], bounds=owner["bounds"], k=raw["k"][0], K=raw["K"][0],
        E=raw["E"][0], chi=raw["chi"][0], E_type="transported" if raw["outer"] else "dependent",
        Ux=G[0,0], Ur=G[0,2], Wx=G[2,0], Wr=G[2,2], D=raw["div"][0],
        R=raw["R"][0].tolist(), Pi=raw["Pi"][0], Vxr=raw["V"][0,0,2],
        production=raw["P"][0])


def locate_failures(a, b, comparison):
    owners = ({e["geometry"]: e for e in state.elements} for state in (a,b))
    coarse, fine = owners
    data = {}
    for name in comparison["failures"]:
        row = comparison["rows"][name]
        limits = row["limits"]
        locations = []
        for rms, maximum, peak in zip(row["component_RMS"],row["component_MAX"],row["component_locations"]):
            if rms <= limits[0] and maximum <= limits[1]:
                continue
            x, r = peak["x"], peak["r"]
            e, f = coarse[peak["coarse_owner"]], fine[peak["fine_owner"]]
            first = raw_owner(a,e,x,r)
            second = raw_owner(b,f,x,r)
            locations.append(dict(RMS=rms, MAX=maximum, limit=limits, peak=peak,
                controls_row=peak["component"] == row["peak"]["component"],
                in_experimental_x=INTERIOR[0] <= x <= INTERIOR[1], first=first, second=second))
        data[name] = locations
    return data


def outlet_shear_split(a,b,peak):
    x,r = peak["x"],peak["r"]
    require(abs(x-15.)<1e-12 and abs(r-.08125)<1e-12,"CP4 signed shear peak moved")
    samples = []
    for state,owner in zip((a,b),(peak["coarse_owner"],peak["fine_owner"])):
        e, = [e for e in state.elements if e["geometry"] == owner]
        samples.append(raw_owner(state,e,x,r))
    p,q = samples
    molecular = a.ref.nu*((q["Ur"]+q["Wx"])-(p["Ur"]+p["Wx"]))
    reynolds = q["R"][0][2]-p["R"][0][2]
    delta = q["Vxr"]-p["Vxr"]
    require(abs(delta-molecular+reynolds) <= 1e-12*a.ref.tau
            and abs(delta/a.ref.tau-peak["signed_difference"]) <= 1e-10,
            "direct shear reconstruction or signed peak differs")
    return dict(x=x,r=r,first=p,second=q,delta_molecular_over_tau=molecular/a.ref.tau,
                delta_Rxr_over_tau=reynolds/a.ref.tau,delta_Vxr_over_tau=delta/a.ref.tau)


def replay(root, out, seconds):
    a,b,reproduced = open_pair(root,"B1c","B2",seconds)
    reference = root/"evidence/outlet-second-collar-result.json"
    old = read(reference)["comparison"]
    # Native owner bounds are tuples in memory and JSON arrays on disk.
    for part in ("first","second","rows","interior_rows","failures"):
        require(json.loads(json.dumps(reproduced[part])) == old[part],
                f"CP4 replay differs in {part}")
    require(len(reproduced["failures"]) == 122,"CP4 replay has a different failed-row inventory")
    data = locate_failures(a,b,reproduced)
    shear = outlet_shear_split(a,b,reproduced["rows"]["whole/traction"]["peak"])
    atomic(out,dict(schema="WJH27FlatQualificationCP4ReplayV1",status="PASS",old=identity(reference),
        comparison=reproduced,failed_components=data,shear_identity=shear,
        hidden_component_maxima=sum(not p["controls_row"] for rows in data.values() for p in rows)))
    print("exact 690-row replay; 122 failed rows", shear["delta_Vxr_over_tau"],flush=True)


def audit(root, case, out, seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    state,binding = load_case(root,case,context,ref,24,time.monotonic()+seconds)
    audit_state(state,binding,context[2],case,out)


def stricter_inner_k(result):
    for group in ("rows","interior_rows"):
        for name,row in result[group].items():
            if name.split("/")[-2:] != ["inner","k"]:
                continue
            row["limits"] = [1e-4,1e-4]
            if max(row["RMS"],row["MAX"]) > 1e-4:
                row["status"] = "NOT_MET"
                row["peak"] = row["component_locations"][0]
                row["RMS_controlling_component"] = row["components"][0]
                row["RMS_component_peak"] = row["peak"]
    result["failures"] = [name for name,row in result["rows"].items() if row["status"] != "MET"]


def compare(root, base, fine, out, seconds):
    a,b,result = open_pair(root,base,fine,seconds)
    if fine in ("Bq","Bt"):
        require(result["first"]["mesh"] == result["second"]["mesh"],
                "a final control must use the exact canonical geometry")
        stricter_inner_k(result)
    fixed_misses = [name for name,row in result["interior_rows"].items() if row["status"] != "MET"]
    atomic(out,dict(schema="WJH27FlatQualificationSignedComparisonV1",cases=[base,fine],
                    comparison=result,failed_components=locate_failures(a,b,result),
                    fixed_region_failures=fixed_misses))
    print(base,fine,len(result["rows"]),len(result["failures"]),len(fixed_misses),flush=True)


def blend_probe(root, base, fine, out, seconds):
    ctx = get_context(root)
    ref = po.load_reference(ctx[2])
    deadline = time.monotonic()+seconds
    a,pa = load_case(root,base,ctx,ref,20,deadline)
    b,pb = load_case(root,fine,ctx,ref,20,deadline)
    lo,mid,hi = ctx[0]["M1c"]["original_collar_owner_r"]
    radii = [.05,.06,.06875,lo,mid,hi,.1125,.11541061001735,ref.H]
    x_stations = [-10.,-2.5,0.,3.75,5.,8.,15.]
    data = {}
    for name,state in ((base,a),(fine,b)):
        rows = []
        for x in x_stations:
            columns = [e for e in state.elements if e["bounds"][0]-1e-12 <= x <= e["bounds"][1]+1e-12]
            for r in radii:
                for e in columns:
                    if e["bounds"][2]-1e-12 <= r <= e["bounds"][3]+1e-12:
                        rows.append(dict(x=x,r=r,**raw_owner(state,e,x,r)))
        data[name] = rows
    atomic(out,dict(schema="WJH27FlatBlendOwnerProbeV1",first=pa,second=pb,
                    x_stations=x_stations,r_stations=radii,raw_owners=data))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    action = parser.add_subparsers(dest="mode",required=True)
    setup = action.add_parser("pin")
    setup.add_argument("out",type=Path)
    setup.add_argument("context",type=Path)
    admission = action.add_parser("admit-b3")
    admission.add_argument("out",type=Path)
    recovery = action.add_parser("forecast-b3-recovery")
    recovery.add_argument("out",type=Path)
    recovery.add_argument("--start",type=int,choices=(2,4,6),default=2)
    inputs = action.add_parser("prepare")
    inputs.add_argument("case",choices=QUAL_CASES)
    inputs.add_argument("admission",type=Path)
    inputs.add_argument("--continuation",type=int,default=0)
    for mode in ("replay","audit","compare","blend-probe"):
        task = action.add_parser(mode)
        task.add_argument("out",type=Path)
        task.add_argument("--seconds",required=True,type=float)
        if mode != "replay":
            task.add_argument("--case",required=True)
        if mode in ("compare","blend-probe"):
            task.add_argument("--base",required=True)
    args = parser.parse_args()
    root = args.root
    if args.mode == "pin":
        pin(root,args.out,args.context)
    elif args.mode == "admit-b3":
        admit_b3(root,args.out)
    elif args.mode == "forecast-b3-recovery":
        forecast_b3_recovery(root,args.out,args.start)
    elif args.mode == "prepare":
        prepare(root,args.case,args.continuation,args.admission)
    elif args.mode == "replay":
        replay(root,args.out,args.seconds)
    elif args.mode == "audit":
        audit(root,args.case,args.out,args.seconds)
    elif args.mode == "blend-probe":
        blend_probe(root,args.base,args.case,args.out,args.seconds)
    else:
        compare(root,args.base,args.case,args.out,args.seconds)


if __name__ == "__main__":
    main()
