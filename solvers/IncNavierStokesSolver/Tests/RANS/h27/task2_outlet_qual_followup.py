#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze and admit one follow-up flat radial, streamwise, or algebraic control."""

import argparse
import math
from pathlib import Path

import physical_oracle as po
from resources import atomic, identity
from task2_outlet_baseline import mesh_axes, read, verified
from task2_outlet_qualification import (GRACE, HIGHER_QUADRATURE, PLAN, PUBLICATION_SECONDS,
    QUADRATURE, SHIFT, TIGHTER_TOLERANCES, TOLERANCES, case_dir, pass_audit, planned_attempt,
    retry_before_checkpoint)
from task2r_run import GIB, balance, memory_headroom, memory_snapshot, require, workers

CONTROLS = ("B4", "Bn", "Bn2", "Bf", "Bw", "Bw2", "Btip", "Bq", "Bt")


def audit_binding(root, case):
    audit, witness = pass_audit(root, case)
    binding = audit["provenance"]
    verified(binding["state"])
    return binding, witness


def split_bands(zs, bands):
    added = [(a+b)/2 for lo,hi in bands for a,b in zip(zs,zs[1:])
             if lo-1e-12 <= a-SHIFT and b-SHIFT <= hi+1e-12]
    require(added and len(added)==len(set(added)), "radial bands overlap or have no native intervals")
    return sorted(zs+added)


def adjacent_H(zs):
    j, = [i for i, z in enumerate(zs) if abs(z-SHIFT-.15)<1e-12]
    return [zs[j-1]-SHIFT,zs[j+1]-SHIFT]


def station_columns(xs, positions):
    return {(a,b) for a,b in zip(xs,xs[1:])
            if any(a-1e-12<=x<=b+1e-12 for x in positions)}


def horizontal_owners(xs, base, comparison, capacity):
    stations = [x for x in po.STATIONS if -2.5 <= x <= 8.]
    tip = set(zip(xs[-3:-1],xs[-2:]))
    full = tip | station_columns(xs,stations)
    if len(full) <= capacity:
        selected = full
    else:
        require(comparison is not None,"bounded x witness needs measured peak-owner evidence")
        paired = read(comparison)["comparison"]
        require(base in (paired["first"]["state"],paired["second"]["state"]),
                "bounded streamwise control is missing a direct comparison with its base")
        # Preserve paired station sides and cover the locations of the main
        # stress, pressure, shear and gradient peaks in the experimental strip.
        observed = (paired["interior_rows"]["interior/whole/"+f]
                    for f in ("R","Pi","traction","G"))
        peaks = [r["component_locations"][r["component_MAX"].index(r["MAX"])]["x"]
                 for r in observed]
        sentinels = [-2.5,0.,3.75,5.,6.25,7.25,8.,*peaks]
        selected = tip | station_columns(xs,sentinels)
        if len(selected) <= capacity:
            for x in sorted(stations, key=lambda x:
                    min(abs(x-v) for v in sentinels), reverse=True):
                expanded = selected | station_columns(xs,[x])
                if len(expanded) <= capacity:
                    selected = expanded
                    sentinels.append(x)
    coverage = {f"{x:g}":dict(challenged=len(station_columns(xs,[x]) & selected),
                            total=len(station_columns(xs,[x]))) for x in stations}
    added = [(a+b)/2 for a,b in sorted(selected)]
    return sorted(xs+added), sorted(selected), dict(full_columns=sorted(full), stations=coverage,
                full_coverage=selected==full, allowable_insertions=capacity)


def native_chain(root, case):
    binding,audit = audit_binding(root,case)
    cells,seconds,peak,anchors = native_cost_chain(root,binding)
    return cells, seconds, peak, anchors+[audit]


def native_cost_chain(root, binding):
    path = verified(binding["state"]).parent
    anchors, seconds, peak = [], 0., 0
    visited = set()
    while True:
        require(path not in visited,"cyclic accepted Newton provenance")
        visited.add(path)
        manifest = read(path/"inputs.json")
        op = root/"operations"/(path.name+"-steady")/"receipt.json"
        receipt = read(op)["result"]
        require(receipt["reason"] is None,"native cost anchor contains an external interruption")
        anchors.append(identity(op))
        seconds += receipt["seconds"]
        peak = max(peak,receipt["peak_tree_rss_bytes"])
        guess = Path(manifest["initial_guess"]["path"])
        if guess.name != "state.steady-candidate.json":
            break
        require(guess.parent.parent == root/"science","unexpected Newton seed location")
        path = guess.parent
    return binding["cells"], seconds, peak, anchors


def controlling_bands(comparison_path):
    report = read(comparison_path)
    require(report["cases"][-1] == "Bw" and report["comparison"]["failures"],
            "repeated broad-band test requires an observed miss on Bw")
    admission = read(verified(report["comparison"]["second"]["admission"]))
    points = [row["peak"] for rows in report["failed_components"].values() for row in rows]
    parents = {tuple(p["owner_bounds"][0][2:4]) for p in points
               if p["owner_bounds"][0][2:4] != p["owner_bounds"][1][2:4]}
    bands = sorted((lo,hi) for lo,hi in parents
                   if any(a-1e-12 <= lo < hi <= b+1e-12 for a,b in admission["radial_bands"]))
    require(bands,"no failed owner was in the newly split radial bands")
    return bands


def terminal_control(root, xs, binding, comparison_path):
    require(comparison_path is not None,"terminal control needs a measured Bn2 outlet peak")
    report = read(comparison_path)
    require(report["cases"]==["Bn","Bn2"] and report["comparison"]["second"]["state"]==binding["state"]
            and report["comparison"]["failures"] and not report["fixed_region_failures"]
            and all(abs(r["peak"]["x"]-15.)<1e-12
                    for rows in report["failed_components"].values() for r in rows),
            "terminal-only control requires localized outlet misses")
    old = [("B0","outlet-b0-continue-1-steady"),("B1","outlet-b1-steady")]
    densities,anchors = {},[]
    for case,label in old:
        audit,witness = audit_binding(root,case)
        path = root/"operations"/label/"receipt.json"
        result = read(path)["result"]
        require(result["reason"] is None and result["peak_tree_rss_bytes"]>0,
                "original streamwise resource anchor differs")
        densities[case] = result["peak_tree_rss_bytes"]/math.prod(audit["cells"])
        anchors.extend((witness,identity(path)))
    spread = max(densities.values())/min(densities.values())-1
    factor = max(1.04,1+spread)
    columns = list(zip(xs[-3:-1],xs[-2:]))
    nodes = sorted(xs+[(a+b)/2 for a,b in columns])
    return nodes,columns,factor,anchors,dict(scope="TWO_TERMINAL_COLUMNS",
        full_experimental_coverage=False, original_x_mesh_bytes_per_cell=densities,
        original_x_normalized_rss_spread=spread)


def freeze(root, case, base, out, previous_comparison):
    evidence = root/"evidence"
    plan_path = evidence/PLAN
    plan = read(plan_path)
    binding, accepted_audit = audit_binding(root,base)
    xs,zs,_ = mesh_axes(verified(binding["mesh"]))
    saved_x,saved_z = xs,zs
    patch = plan["meshes"]["full"]["support_r"]
    bands,columns = [],[]
    x_coverage = None
    q,tol = QUADRATURE,TOLERANCES
    b3,wall,peak,anchors = native_chain(root,"B3")
    b2,_,peak_b2,anchors_b2 = native_chain(root,"B2")
    densities = {"B3":peak/math.prod(b3),"B2":peak_b2/math.prod(b2)}
    margin = 1.10
    plan_anchors = anchors+anchors_b2
    if base not in ("B2","B3"):
        base_cells,base_wall,base_peak,base_anchors = native_chain(root,base)
        densities[base] = base_peak/math.prod(base_cells)
        wall = max(wall,base_wall*math.prod(b3)/math.prod(base_cells))
        plan_anchors += base_anchors
    rss_per_cell = max(densities.values())
    if case == "Bw":
        bands = [[.05,patch[0]],[patch[1],.11541061001735],adjacent_H(zs)]
        zs = split_bands(zs,bands)
    elif case == "Bw2":
        require(base=="Bw" and previous_comparison is not None,"repeat must start from the first broad witness")
        bands = controlling_bands(previous_comparison)
        zs = split_bands(zs,bands)
    elif case == "B4":
        require(base=="B3" and previous_comparison is not None,"further collar test needs a measured B3 miss")
        observed = read(previous_comparison)["comparison"]
        require(observed["second"]["state"] == binding["state"] and observed["failures"],
                "B4 requires the exact preceding collar comparison")
        bands = [patch]
        zs = split_bands(zs,bands)
    elif case == "Bn":
        require(base=="B4" and previous_comparison is not None,
                "neighbor control requires the measured lower-patch comparison")
        observed = read(previous_comparison)["comparison"]
        peak = observed["rows"]["x=15/whole/traction"]["peak"]
        lo,hi = .06875,patch[0]
        require(observed["second"]["state"]==binding["state"]
                and abs(peak["r"]-hi)<1e-12 and observed["failures"],
                "neighbor control requires the measured lower collar-edge shear peak")
        first_split = split_bands(zs,[[lo,hi]])
        mid = (lo+hi)/2
        zs = split_bands(first_split,[[mid,hi]])
        bands = [[lo,hi]]
    elif case == "Bn2":
        require(base=="Bn" and previous_comparison is not None,
                "repeat of the lower neighbor requires the accepted Bn comparison")
        observed = read(previous_comparison)
        require(observed["cases"]==["B4","Bn"]
                and observed["comparison"]["second"]["state"]==binding["state"]
                and observed["comparison"]["failures"], "repeat needs measured neighbor sensitivity")
        lo,hi = .06875,patch[0]
        peaks = [r["peak"]["r"] for rows in observed["failed_components"].values() for r in rows]
        require(any(lo-1e-12<=r<=hi+1e-12 for r in peaks),
                "no failed maximum lies on the repeated radial support")
        bands = [[lo,hi]]
        zs = split_bands(zs,bands)
        class4,_,peak4,anchors4 = native_chain(root,"B4")
        densities["B4"] = peak4/math.prod(class4)
        plan_anchors += anchors4
        rss_per_cell = max(densities.values())
        # An empirical margin for the same lower-neighbor support, using
        # twice the observed peak-per-cell spread across the measured ladder.
        margin = max(1.05,1+2*(rss_per_cell/min(densities.values())-1))
    elif case == "Bf":
        require(base=="Bn2" and previous_comparison is not None,
                "paired flank test needs its accepted Bn2 base and the terminal control")
        x_test = read(previous_comparison)
        lower_test_path = evidence/"outlet-qual-Bn2-Bn.json"
        lower = read(lower_test_path)["comparison"]
        peak = lower["rows"]["x=15/whole/traction"]["peak"]
        edge = peak["r"]
        require(x_test["cases"]==["Bn2","Btip"]
                and x_test["comparison"]["first"]["state"]==binding["state"]
                and lower["second"]["state"]==binding["state"] and abs(edge-.06875)<1e-12,
                "paired flank test does not bind the measured controlling edge")
        i, = [i for i,z in enumerate(zs) if abs(z-SHIFT-edge)<1e-12]
        bands = [[zs[i-1]-SHIFT,zs[i+1]-SHIFT]]
        zs = split_bands(zs,bands)
        require(len(zs)==len(saved_z)+2,"expected one division on each side of the flank")
        # Retain the terminal case as a cost measurement even when its
        # independent audit limits its use to a diagnostic comparison.
        tip_audit_path = evidence/"outlet-Btip-audit.json"
        tip_audit = read(tip_audit_path)
        tip_binding = x_test["comparison"]["second"]
        require(tip_audit["name"]=="Btip" and tip_binding==tip_audit["provenance"],
                "terminal resource record and independent audit must bind the same native state")
        xtip,_,xtip_peak,xtip_anchors = native_cost_chain(root,tip_binding)
        densities["Btip"] = xtip_peak/math.prod(xtip)
        plan_anchors += [identity(lower_test_path),identity(tip_audit_path),*xtip_anchors]
        rss_per_cell = max(densities.values())
        # One added element on either side of the face. Cover the complete
        # observed density spread from the audited radial and x controls.
        margin = max(1.04,rss_per_cell/min(densities.values()))
    elif case == "Btip":
        if base=="Bn2":
            xs,columns,margin,x_anchors,x_coverage = terminal_control(root,xs,binding,previous_comparison)
            plan_anchors += x_anchors
        else:
            max_nx = math.floor(12*GIB/(1.10*rss_per_cell*(len(zs)-1)))
            xs,columns,x_coverage = horizontal_owners(xs,binding["state"],previous_comparison,
                                                      max_nx-(len(xs)-1))
    elif case == "Bq":
        q = HIGHER_QUADRATURE
    else:
        tol = TIGHTER_TOLERANCES
    require(set(saved_x) <= set(xs) and set(saved_z) <= set(zs),
            "follow-up would remove existing owner edges")

    dims = [len(xs)-1,len(zs)-1]
    count = math.prod(dims)
    from_b3 = count / math.prod(b3)
    predicted_time = math.ceil(wall * from_b3)
    reserve_points = count*(math.prod(q)-math.prod(QUADRATURE))
    peak_estimate = math.ceil(margin*rss_per_cell*count + reserve_points*2048)
    if case == "B4" and peak_estimate > 12*GIB:
        # Preserve both sides of every lower-face peak with the affordable lower patch variant.
        radii = [p["peak"]["r"] for rows in read(previous_comparison)["failed_components"].values() for p in rows]
        require(radii and all(patch[0]-1e-12<=r<=.1+1e-12 for r in radii),
                "full third refinement exceeds memory and the lower patch misses a failed peak")
        bands = [[patch[0],.1]]
        zs = split_bands(saved_z,bands)
        dims = [len(xs)-1,len(zs)-1]
        count = math.prod(dims)
        from_b3 = count/math.prod(b3)
        predicted_time = math.ceil(wall*from_b3)
        peak_estimate = math.ceil(margin*rss_per_cell*count)
    if case == "Bf":
        observed = max(read(verified(item))["result"]["seconds"] for item in base_anchors+xtip_anchors
                       if item["path"].endswith("-steady/receipt.json"))
        proposed = [(240,1),(math.ceil(1.12*observed),50),(350,50)]
    elif case in ("Bq","Bt"):
        anchor = root/"operations"/("r3-quadrature-steady" if case=="Bq" else "r3-algebraic-steady")/"receipt.json"
        rec = read(anchor)
        xx,zz,_ = mesh_axes(Path(rec["request"]["command"][-1]))
        scaled = math.ceil(rec["result"]["seconds"]*count/((len(xx)-1)*(len(zz)-1)))
        plan_anchors.append(identity(anchor))
        if case == "Bq":
            proposed = [(max(240,scaled),1),(680,50),(350,50)]
        else:
            proposed = [(max(640,math.ceil(scaled*1.35)),50),(390,50)]
    else:
        proposed = [(240,1),(max(570, math.ceil(predicted_time*.95+50)),50),(390,50)]
    attempts = {case_dir(root,case,i).name:planned_attempt(native,native+115,steps)
                for i,(native,steps) in enumerate(proposed)}
    scale_analysis = max(1, count/math.prod(b2))
    audit_cap = math.ceil(190*scale_analysis)
    compare_cap = math.ceil(110*scale_analysis)
    spends = dict(admission=27, native_and_preparations=sum(item["preparation_seconds"]+
            item["watch_seconds"]+2*GRACE for item in attempts.values()),
            audit=audit_cap+GRACE,comparison=compare_cap+GRACE, checkpoint_evidence=70,
            final_publication=PUBLICATION_SECONDS)
    required = sum(spends.values())
    limit = 10*GIB if peak_estimate <= 10*GIB else 12*GIB
    live = balance(read(root/"accounting.json"))
    memory = memory_snapshot()
    blockers = [name for name, fits in (
        ("wall_time",required<=live["nonprotected_remaining"]),
        ("forecast_RSS",peak_estimate<=limit),
        ("live_job_headroom",memory_headroom(memory)>=limit+2*GIB)) if not fits]
    atomic(out,dict(schema="WJH27FlatQualificationForecastV1",
        status="FITS" if not blockers else "DEFERRED", blockers=blockers, case=case, base=base,
        prior_state=binding["state"], prior_mesh=binding["mesh"], prior_audit=accepted_audit,
        x_edges=xs, physical_z_edges=zs, cells=dims, radial_bands=bands, x_columns=columns,
        x_owner_coverage=x_coverage,
        quadrature=q, tolerances=tol, source=plan["source"], application=plan["application"],
        forecast_peak_tree_rss_bytes=peak_estimate, rss_limit_bytes=limit,
        memory_estimate_inputs=dict(refinement_margin=margin, measured_bytes_per_cell=densities,
            extra_native_points=reserve_points,
            bytes_per_extra_point=2048), anchors=plan_anchors, native_wall_anchor_scaled_seconds=predicted_time,
        attempts=attempts, audit_cap_seconds=audit_cap, compare_cap_seconds=compare_cap,
        reservation_caps=spends,required_seconds=required,memory=memory,live=live,
        prior_comparison=identity(previous_comparison) if previous_comparison else None,
        admission_plan=identity(plan_path)))
    print(case,base,dims,"seconds",required,"RSS GiB",peak_estimate/GIB,blockers,flush=True)


def recover_flank(root, admission, out):
    original = read(admission)
    first = case_dir(root,"Bf").name
    path = root/"operations"/f"{first}-steady"/"receipt.json"
    stop = read(path)
    require(original["case"]=="Bf" and original["status"]=="ADMITTED"
            and stop["request"]["label"]==f"{first}-steady"
            and stop["result"]["status"]=="FAIL" and retry_before_checkpoint(stop["result"])
            and not (case_dir(root,"Bf")/"state.steady-candidate.json").exists()
            and not (case_dir(root,"Bf")/"state.result.json").exists(),
            "flank readmission requires an interrupted launch with no typed state")
    expected = original["attempts"]
    require(list(expected)==[case_dir(root,"Bf",i).name for i in range(3)],
            "flank launch recovery may be staged once")
    attempts = {case_dir(root,"Bf",i).name:a for i,a in enumerate(expected.values(),start=1)}
    live = balance(read(root/"accounting.json"))
    memory = memory_snapshot()
    blockers = [name for name, fits in (
        ("wall_time",original["required_seconds"]<=live["nonprotected_remaining"]),
        ("live_job_headroom",memory_headroom(memory)>=original["rss_limit_bytes"]+2*GIB)) if not fits]
    atomic(out,dict(original,schema="WJH27FlatQualificationForecastV1",
        status="FITS" if not blockers else "DEFERRED",blockers=blockers, attempts=attempts,
        previous_admission=identity(admission),interruption=identity(path),memory=memory,live=live))
    print("Bf launch recovery",original["required_seconds"],stop["result"]["seconds"],blockers,flush=True)


def admit(root, forecast, out):
    expected = read(forecast)
    live = balance(read(root/"accounting.json"))
    memory = memory_snapshot()
    remaining = expected["required_seconds"]-expected["reservation_caps"]["admission"]
    require(expected["schema"]=="WJH27FlatQualificationForecastV1" and expected["status"]=="FITS"
            and remaining<=min(live["tranche_remaining"],live["nonprotected_remaining"])
            and memory_headroom(memory)>=expected["rss_limit_bytes"]+2*GIB,
            "live budget or job memory no longer supports the complete admitted case")
    atomic(out,dict(expected, schema="WJH27FlatQualificationAdmissionV1",
        status="ADMITTED", forecast=identity(forecast), live_at_admission=live,
        memory_at_admission=memory, worker_census=workers()))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    sub = parser.add_subparsers(dest="mode",required=True)
    design = sub.add_parser("freeze")
    design.add_argument("case",choices=CONTROLS)
    design.add_argument("base")
    design.add_argument("out",type=Path)
    design.add_argument("--comparison",type=Path)
    admission = sub.add_parser("admit")
    admission.add_argument("forecast",type=Path)
    admission.add_argument("out",type=Path)
    recovery = sub.add_parser("recover-flank")
    recovery.add_argument("admission",type=Path)
    recovery.add_argument("out",type=Path)
    args = parser.parse_args()
    if args.mode=="freeze":
        freeze(args.root,args.case,args.base,args.out,args.comparison)
    elif args.mode=="recover-flank":
        recover_flank(args.root,args.admission,args.out)
    else:
        admit(args.root,args.forecast,args.out)


if __name__ == "__main__":
    main()
