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
    QUADRATURE, SHIFT, TIGHTER_TOLERANCES, TOLERANCES, case_dir, pass_audit, planned_attempt)
from task2r_run import GIB, balance, memory_headroom, memory_snapshot, require, workers

CONTROLS = ("B4", "Bw", "Bw2", "Btip", "Bq", "Bt")


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


def horizontal_owners(xs):
    stations = [x for x in po.STATIONS if -2.5 <= x <= 8.]
    columns = [(a,b) for a,b in zip(xs,xs[1:])
               if a>=xs[-3] or any(a-1e-12<=x<=b+1e-12 for x in stations)]
    return sorted(xs+[(a+b)/2 for a,b in columns]), columns


def native_chain(root, case):
    binding,audit = audit_binding(root,case)
    path = Path(binding["state"]["path"]).parent
    anchors, seconds, peak = [], 0., 0
    while True:
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
    return binding["cells"], seconds, peak, anchors+[audit]


def controlling_bands(comparison_path):
    report = read(comparison_path)
    require(report["cases"][-1] == "Bw" and report["comparison"]["failures"],
            "repeated broad-band test requires an observed miss on Bw")
    admission = read(verified(report["comparison"]["second"]["admission"]))
    points = [row["peak"] for rows in report["failed_components"].values() for row in rows]
    bands = []
    for lo, hi in admission["radial_bands"]:
        if any(lo-1e-12 <= p["owner_bounds"][1][2]
               and p["owner_bounds"][1][3] <= hi+1e-12 for p in points):
            bands.append([lo,hi])
    require(bands,"no failed owner was in the newly split radial bands")
    return bands


def freeze(root, case, base, out, previous_comparison):
    evidence = root/"evidence"
    plan_path = evidence/PLAN
    plan = read(plan_path)
    binding, accepted_audit = audit_binding(root,base)
    xs,zs,_ = mesh_axes(verified(binding["mesh"]))
    saved_x,saved_z = xs,zs
    patch = plan["meshes"]["full"]["support_r"]
    bands,columns = [],[]
    q,tol = QUADRATURE,TOLERANCES
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
    elif case == "Btip":
        xs,columns = horizontal_owners(xs)
    elif case == "Bq":
        q = HIGHER_QUADRATURE
    else:
        tol = TIGHTER_TOLERANCES
    require(set(saved_x) <= set(xs) and set(saved_z) <= set(zs),
            "follow-up would remove existing owner edges")

    dims = [len(xs)-1,len(zs)-1]
    count = math.prod(dims)
    b3,wall,peak,anchors = native_chain(root,"B3")
    b2,_,peak_b2,anchors_b2 = native_chain(root,"B2")
    from_b3 = count / math.prod(b3)
    base_peak_estimate = max(peak * from_b3, peak_b2 * count / math.prod(b2))
    predicted_time = math.ceil(wall * from_b3)
    plan_anchors = anchors+anchors_b2
    reserve_points = count*(math.prod(q)-math.prod(QUADRATURE))
    peak_estimate = math.ceil(base_peak_estimate*1.10 + reserve_points*2048)
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
        peak_estimate = math.ceil(max(peak*from_b3,peak_b2*count/math.prod(b2))*1.10)
    if case in ("Bq","Bt"):
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
        quadrature=q, tolerances=tol, source=plan["source"], application=plan["application"],
        forecast_peak_tree_rss_bytes=peak_estimate, rss_limit_bytes=limit,
        memory_estimate_inputs=dict(refinement_margin=1.10, extra_native_points=reserve_points,
            bytes_per_extra_point=2048), anchors=plan_anchors, native_wall_anchor_scaled_seconds=predicted_time,
        attempts=attempts, audit_cap_seconds=audit_cap, compare_cap_seconds=compare_cap,
        reservation_caps=spends,required_seconds=required,memory=memory,live=live,
        prior_comparison=identity(previous_comparison) if previous_comparison else None,
        admission_plan=identity(plan_path)))
    print(case,base,dims,"seconds",required,"RSS GiB",peak_estimate/GIB,blockers,flush=True)


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
    args = parser.parse_args()
    if args.mode=="freeze":
        freeze(args.root,args.case,args.base,args.out,args.comparison)
    else:
        admit(args.root,args.forecast,args.out)


if __name__ == "__main__":
    main()
