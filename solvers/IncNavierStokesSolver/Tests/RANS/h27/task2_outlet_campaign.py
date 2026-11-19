#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze measured admissions for the steady flat outlet campaign."""

import argparse
import math
from pathlib import Path

import physical_oracle as po
from resources import atomic, identity, size
from task2_outlet_baseline import mesh_axes, read, verified
from task2_outlet_qualification import (GRACE, HIGHER_QUADRATURE, PUBLICATION_SECONDS,
    QUADRATURE, SHIFT, TIGHTER_TOLERANCES, TOLERANCES, case_dir, pass_audit, planned_attempt)
from task2_outlet_qual_followup import adjacent_H, admit, native_chain, split_bands, station_columns
from task2r_run import (CAMPAIGN_AMENDMENT, CAMPAIGN_RSS_GIB, GIB, balance, campaign_case,
    memory_headroom, memory_snapshot, require, storage)

PRIOR = "d28a8f9c46d600e70cee18172f6452faab19f08d"
SHEAR_TRACK = [(15., .06), (15., .06875), (-9.939670769620442, .059965642995925474)]
BLEND_SUPPORT = [.05, .11541061001735]


def radial_grading(edges, window, ratio):
    """Bisect the wider neighbor at faces where the frozen ratio is exceeded."""
    lo,hi = window
    require(0<=lo<hi and math.isfinite(ratio) and ratio>=2., "invalid dyadic grading prospect")
    while True:
        additions = set()
        for left,face,right in zip(edges,edges[1:],edges[2:]):
            if not lo-1e-12<=face-SHIFT<=hi+1e-12:
                continue
            left_width,right_width = face-left,right-face
            if left_width>ratio*right_width+1e-12:
                additions.add((left+face)/2)
            elif right_width>ratio*left_width+1e-12:
                additions.add((face+right)/2)
        if not additions:
            return edges
        require(not additions.intersection(edges), "radial grading cannot add a distinct edge")
        edges = sorted([*edges,*additions])


def axes_for(root, design, x, z):
    case = design["case"]
    require(campaign_case(case), "case is outside the outlet campaign")
    quadrature, tolerances = QUADRATURE, TOLERANCES
    kind = case[1]
    new_x, new_z, columns = x, z, []
    if case == "Bfl":
        pinned = read(root/"evidence/outlet-qual-next-experiment.json")
        require(design["base"]=="Bf", "lower-band repeat needs its registered Bf state")
        new_z = split_bands(z, [[.06,.06875]])
        require([new_x,new_z]==[pinned["x_edges"],pinned["physical_z_edges"]]
                and len(new_x)==112 and len(new_z)==95, "archived 111x94 axes differ")
    elif kind in ("r","c"):
        if design.get("radial_bands"):
            new_z = split_bands(z, design["radial_bands"])
        if design.get("radial_grading"):
            grading = design["radial_grading"]
            new_z = radial_grading(new_z,grading["faces_r"],grading["max_neighbor_ratio"])
        require(new_z!=z, "radial witness must refine its declared support")
    elif kind == "h":
        lo,hi = adjacent_H(z)
        support = design.get("radial_bands", [[lo,hi]])
        require(len(support)==1 and support[0][0]<=lo<.15<hi<=support[0][1],
                "H support must contain both immediate owners")
        new_z = split_bands(z, support)
    elif kind == "q":
        quadrature = HIGHER_QUADRATURE
    elif kind == "t":
        tolerances = TIGHTER_TOLERANCES
    if kind in ("x","c"):
        if design.get("x_scope") == "experimental-and-terminal":
            experimental = [s for s in po.STATIONS if -2.5 <= s <= 8.]
            columns = sorted(station_columns(x, experimental) | set(zip(x[-3:-1],x[-2:])))
            if design.get("include_inlet"):
                columns = sorted(set(columns) | {(x[0],x[1])})
        else:
            columns = [tuple(ends) for ends in design["x_columns"]]
        existing = set(zip(x,x[1:]))
        require(columns and len(columns)==len(set(columns)) and set(columns)<=existing,
                "streamwise witness requires exact distinct original owners")
        new_x = sorted(x+[(a+b)/2 for a,b in columns])
    require(len(new_x)>len(x) or len(new_z)>len(z) or kind in ("q","t"),
            "named solve has no registered spatial or equation change")
    return new_x, new_z, columns, quadrature, tolerances


def measure(root, cases, cells, quadrature):
    """Scale successful checkpoints on the same steady path and preserve each anchor."""
    count = math.prod(cells)
    startup, convergence, total = [], [], []
    densities, disk_rates, anchors, timings = {}, [], [], {}
    for case in dict.fromkeys(cases):
        dims, charged, peak, chain = native_chain(root,case)
        n = math.prod(dims)
        densities[case] = peak/n
        total.append(charged*count/n)
        receipts = []
        for link in chain:
            path = verified(link)
            if path.name != "receipt.json":
                continue
            receipt = read(path)
            op = receipt["result"]
            session = Path(receipt["request"]["command"][-1])
            disk_rates.append(size(session.parent)/n)
            receipts.append(op["seconds"])
            if op["status"]=="PASS":
                convergence.append(op["seconds"]*count/n)
            else:
                startup.append(op["seconds"]*count/n)
        timings[case] = dict(cells=dims,native_seconds=receipts,peak_tree_rss_bytes=peak)
        anchors.extend(chain)
    # Preserve the measured factor-storage spread and reserve native storage
    # for the extra (8x12) panel points on the same polynomial grid.
    margin = max(1.04,max(densities.values())/min(densities.values()))
    extra = count*(math.prod(quadrature)-math.prod(QUADRATURE))
    peak = math.ceil(max(densities.values())*count*margin+extra*2048)
    scale_q = math.sqrt(math.prod(quadrature)/math.prod(QUADRATURE))
    first = math.ceil(max(240, 0.82*max(startup,default=300))*scale_q)
    again = math.ceil(max(640, 1.13*max(convergence), .72*max(total))*scale_q)
    retry = max(380,math.ceil(first*1.25))
    disk = math.ceil(1.10*max(disk_rates)*count*3)
    return dict(peak=peak,first=first,continuation=again,retry=retry,disk=disk,
                anchors=anchors,timings=timings,densities=densities,margin=margin,
                extra_quadrature_points=extra,quadrature_cost_factor=scale_q)


def freeze(root, design_path, out):
    design = read(design_path)
    case, base = design["case"], design["base"]
    base_audit, audit_receipt = pass_audit(root,base)
    prior = base_audit["provenance"]
    prior_mesh, prior_state = verified(prior["mesh"]), verified(prior["state"])
    x,z,_ = mesh_axes(prior_mesh)
    new_x,new_z,columns,quadrature,tolerances = axes_for(root,design,x,z)
    cells = [len(new_x)-1,len(new_z)-1]
    for comparator in design.get("comparison_bases",[base]):
        pass_audit(root,comparator)
    data = measure(root, ["Bn2","Bf",base,*design.get("extra_cost_cases",[])], cells, quadrature)
    attempts = {case_dir(root,case,i).name:planned_attempt(seconds,seconds+115,n)
                for i,(seconds,n) in enumerate(((data["first"],1),(data["continuation"],50),
                                                 (data["retry"],50)))}
    volume_scale = math.prod(cells)/(111*92)
    audit_seconds = math.ceil(254*max(1,volume_scale))
    compare_seconds = math.ceil(147*max(1,volume_scale))
    comparison_bases = design.get("comparison_bases",[base])
    extra_analysis = 90 if case=="Bfl" else design.get("extra_analysis_seconds",0)
    caps = dict(admission=27, preparations_and_native=sum(a["preparation_seconds"]+
        a["watch_seconds"]+2*GRACE for a in attempts.values()), independent_audit=audit_seconds+GRACE,
        direct_comparisons=len(comparison_bases)*(compare_seconds+GRACE),
        diagnostic_evidence=extra_analysis, checkpoint_publication=57,
        final_evidence=PUBLICATION_SECONDS)
    required = sum(caps.values())
    live = balance(read(root/"accounting.json"))
    memory = memory_snapshot()
    live_storage = storage(root)
    rss_limit = max(10, math.ceil(data["peak"]/GIB))*GIB
    blockers = [label for label, ok in (
        ("wall_time", required<=live["nonprotected_remaining"]),
        ("named_peak_RSS", rss_limit<=CAMPAIGN_RSS_GIB*GIB),
        ("host_or_cgroup", memory_headroom(memory)>=rss_limit+2*GIB),
        ("science_storage", live_storage["science_bytes"]+data["disk"]<=8*GIB),
        ("build_storage", live_storage["source_build_runtime_bytes"]<=8*GIB),
        ("free_space",live_storage["free_bytes"]-data["disk"]>=15*GIB)) if not ok]
    amended = root/"evidence"/f"{CAMPAIGN_AMENDMENT}.json"
    require(read(amended)["status"]=="PASS", "missing campaign amendment")
    atomic(out,dict(schema="WJH27FlatQualificationForecastV1",
        case=case,base=base,source=prior["source"],application=prior["executable"],
        prior_state=identity(prior_state),prior_mesh=identity(prior_mesh),prior_audit=audit_receipt,
        design=identity(design_path),campaign_amendment=identity(amended),
        cells=cells,x_edges=new_x,physical_z_edges=new_z,x_columns=columns,
        quadrature=quadrature,tolerances=tolerances,comparison_bases=comparison_bases,
        forecast_peak_tree_rss_bytes=data["peak"],rss_limit_bytes=rss_limit,
        memory_estimate_inputs=dict(margin=data["margin"],bytes_per_cell=data["densities"],
            extra_quadrature_points=data["extra_quadrature_points"],extra_bytes_per_point=2048),
        measured_native_times=data["timings"],quadrature_cost_factor=data["quadrature_cost_factor"],
        anchors=data["anchors"],attempts=attempts,audit_cap_seconds=audit_seconds,
        compare_cap_seconds=compare_seconds,required_seconds=required,reservation_caps=caps,
        diagnostic_cap_seconds=max(0,extra_analysis-GRACE),live=live,memory=memory,
        storage=live_storage,forecast_additional_science_bytes=data["disk"],
        status="DEFERRED" if blockers else "FITS",blockers=blockers))
    print(case,base,cells,"seconds",required,"RSS GiB",data["peak"]/GIB,blockers,flush=True)


def pin(root, initial_forecast, out):
    ev = root/"evidence"
    hold_path = ev/"outlet-qual-task3-launch-hold.json"
    hold = read(hold_path)
    require(hold["status"]=="NOT_READY_DO_NOT_LAUNCH" and hold["TASK2_COMPLETE"]==hold["TASK3_READY"]=="NO",
            "prior task launch hold differs")
    verified(hold["application"])
    verified(hold["reference"])
    for name, item in hold["scientific_observations"].items():
        for key in ("state","mesh"):
            verified(item[key])
        require(read(ev/f"outlet-{name}-audit.json")["provenance"]["state"]==item["state"],
                "prior hold and audited state differ")
    first = read(initial_forecast)
    require(first["case"]=="Bfl" and first["status"]=="FITS" and first["cells"]==[111,94],
            "complete first experiment is not admitted")
    atomic(out,dict(schema="WJH27FlatQualificationCampaignV1", baseline_commit=PRIOR,
        amendment=identity(ev/f"{CAMPAIGN_AMENDMENT}.json"),prior_hold=identity(hold_path),
        first_forecast=identity(initial_forecast),application=hold["application"],
        reference=hold["reference"],native_source=hold["native_source"],
        resumed_observations=hold["scientific_observations"],
        prior_flank_comparison=identity(ev/"outlet-qual-Bf-Bn2.json"),
        comparison_rows=690,stations=po.STATIONS,experimental_x=[-2.5,8.],
        spatial_limits=po.NUMERICAL,current_limits=po.RESIDUAL_LIMITS,
        physical_limits=dict(po.ENVELOPES,wall_shear=[.05,.10]),
        inlet_limits=dict(U=1e-5,k=1e-4,G=1e-3),track_samples=SHEAR_TRACK,
        higher_quadrature=HIGHER_QUADRATURE,tighter_tolerances=TIGHTER_TOLERANCES,
        controls_inner_k_RMS_MAX=1e-4, H=.15,blend=[.06,.1125],blend_with_flanks=BLEND_SUPPORT,
        task_status="OPEN",TASK2_COMPLETE="NO",TASK3_READY="NO"))
    print("campaign pinned",first["cells"],first["required_seconds"],flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    modes = parser.add_subparsers(dest="mode",required=True)
    for mode in ("freeze","pin","admit"):
        task = modes.add_parser(mode)
        task.add_argument("input",type=Path)
        task.add_argument("output",type=Path)
    args=parser.parse_args()
    {"freeze":freeze,"pin":pin,"admit":admit}[args.mode](args.root,args.input,args.output)


if __name__ == "__main__":
    main()
