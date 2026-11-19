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
    QUADRATURE, SHIFT, TIGHTER_TOLERANCES, TOLERANCES, accepted_path, case_dir, pass_audit, planned_attempt)
from task2_outlet_qual_followup import (adjacent_H, admit, native_chain, native_cost_chain,
    split_bands, station_columns)
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


def measure(root, cases, cells, quadrature, diagnostic=None):
    """Scale successful checkpoints on the same steady path and preserve each anchor."""
    count = math.prod(cells)
    startup, convergence, total = [], [], []
    densities, disk_rates, anchors, timings = {}, [], [], {}
    for case in dict.fromkeys(cases):
        if diagnostic and case==diagnostic["case"]:
            dims,charged,peak,chain = native_cost_chain(root,diagnostic["binding"])
            chain += [diagnostic["audit"],diagnostic["probe"]]
        else:
            dims,charged,peak,chain = native_chain(root,case)
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
                extra_quadrature_points=extra,quadrature_cost_factor=scale_q,
                scaled_first_charge=max(startup,default=300))


def readmit_memory_stop(root, previous, base, x, z, q, t, data):
    """Reforecast the same mesh from an observed, charged RSS interruption."""
    stopped = case_dir(root,previous)
    receipt = root/"operations"/f"{stopped.name}-steady"/"receipt.json"
    event = read(receipt)["result"]
    inputs = read(stopped/"inputs.json")
    admission = read(verified(inputs["admission"]))
    require(event["status"]=="FAIL" and event["reason"]=="RSS_LIMIT"
            and event["peak_tree_rss_bytes"]>event["rss_limit_bytes"]
            and not (stopped/"state.result.json").exists()
            and not (stopped/"state.steady-candidate.json").exists()
            and base==admission["base"] and [x,z]==[inputs["x_edges"],inputs["physical_z_edges"]]
            and [q,t]==[inputs["quadrature"],inputs["tolerances"]],
            "RSS retry must preserve the stopped unqualified mesh and numerical recipe")
    operations = root/"operations"
    growth_anchors = ("outlet-qual-btip", "outlet-qual-bfl")
    increases = {}
    for label in growth_anchors:
        anchors = [operations/f"{label}{suffix}-steady"/"receipt.json"
                   for suffix in ("","-continue-1")]
        first,final = (read(a)["result"] for a in anchors)
        require(first["reason"] is None and final["status"]=="PASS"
                and first["peak_tree_rss_bytes"]>0, "measured factor growth anchor differs")
        increases[label] = final["peak_tree_rss_bytes"]/first["peak_tree_rss_bytes"]
        data["anchors"].extend(identity(a) for a in anchors)
    memory_growth = max(1.,*increases.values())
    measured_overrun = max(1.,event["seconds"]/data["scaled_first_charge"])
    data["peak"] = math.ceil(event["peak_tree_rss_bytes"]*memory_growth*data["margin"])
    data["first"] = max(data["first"],math.ceil(event["seconds"]*1.15))
    data["continuation"] = math.ceil(data["continuation"]*measured_overrun)
    data["retry"] = max(data["retry"],data["first"])
    data["anchors"].append(identity(receipt))
    return dict(stop=identity(receipt),observed_peak_tree_rss_bytes=event["peak_tree_rss_bytes"],
        observed_wall_seconds=event["seconds"],observed_factor_growth=increases,
        memory_growth_factor=memory_growth,runtime_overrun_factor=measured_overrun)


def scale_radius_after_stop(root, previous, base, x, z, q, t, data):
    """Use both the early and late memory from a wider stopped radial solve."""
    first = case_dir(root,previous)
    later = case_dir(root,previous,1)
    root_ops = root/"operations"
    start_path = root_ops/f"{first.name}-steady"/"receipt.json"
    late_path = root_ops/f"{later.name}-steady"/"receipt.json"
    start,late = (read(p)["result"] for p in (start_path,late_path))
    inputs = read(later/"inputs.json")
    original = read(verified(inputs["admission"]))
    require(start["reason"] is None and late["reason"]=="RSS_LIMIT"
            and late["peak_tree_rss_bytes"]>late["rss_limit_bytes"]==CAMPAIGN_RSS_GIB*GIB
            and original["base"]==base and x==inputs["x_edges"]
            and set(z)<set(inputs["physical_z_edges"])
            and [q,t]==[inputs["quadrature"],inputs["tolerances"]],
            "staging forecast requires a nested subset of the stopped radial grid")
    scale = (len(z)-1)/inputs["cells"][1]
    past_growth = original["readmission"]["memory_growth_factor"]
    # The late stop already includes much of the growth from the first
    # Newton factor. Charge the larger of the observed late peak and the
    # early peak carried through the measured prior growth.
    large_peak = max(late["peak_tree_rss_bytes"],start["peak_tree_rss_bytes"]*past_growth)
    data["peak"] = max(data["peak"],math.ceil(large_peak*scale*data["margin"]))
    data["first"] = max(data["first"],math.ceil(start["seconds"]*scale*1.04))
    data["continuation"] = max(data["continuation"],
        math.ceil(original["attempts"][later.name]["native_seconds"]*scale))
    data["retry"] = max(data["retry"],data["first"])
    data["anchors"].extend(identity(p) for p in (start_path,late_path))
    return dict(early=identity(start_path),late=identity(late_path),
        early_peak_tree_rss_bytes=start["peak_tree_rss_bytes"],
        late_peak_tree_rss_bytes=late["peak_tree_rss_bytes"],
        modelled_larger_peak_bytes=large_peak,scaled_cell_count_ratio=scale,
        measured_first_to_final_growth=past_growth)


def diagnostic_anchor(root, base, case, probe_name):
    audit_path = root/"evidence"/f"outlet-{base}-audit.json"
    audit = read(audit_path)
    prior = audit["provenance"]
    probe_path = root/"evidence"/probe_name
    probe = read(probe_path)
    require(campaign_case(case) and case[1] in ("q","t") and probe_path.parent==root/"evidence"
            and probe["case"]==base and probe["state"]==prior["state"]
            and [k for k,v in audit["status"].items() if v!="MET"]==["independent_current"]
            and audit["core"]["residuals"]["momentum"]>po.RESIDUAL_LIMITS["momentum"]
            and identity(accepted_path(root,base))==prior["state"],
            "equation diagnostic seed must be the bound native result of the audited momentum miss")
    native = read(verified(prior["operation"]))["result"]
    require(native["status"]=="PASS" and prior["state"] in native["output_files"],
            "equation diagnostic seed lacks an accepted native receipt")
    return dict(case=base,binding=prior,audit=identity(audit_path),probe=identity(probe_path))


def streamwise_schedule(root, case, anchor, x, z, quadrature, tolerances, data):
    require(campaign_case(case) and case[1]=="x" and campaign_case(anchor)
            and anchor[1]=="x" and case!=anchor,
            "segmented timing requires a previous measured streamwise case")
    audit_path = root/"evidence"/f"outlet-{anchor}-audit.json"
    binding = read(audit_path)["provenance"]
    mx,mz,_ = mesh_axes(verified(binding["mesh"]))
    require(identity(accepted_path(root,anchor))==binding["state"] and set(x)<=set(mx) and z==mz
            and binding["quadrature"]==quadrature and binding["tolerances"]==tolerances,
            "timing anchor must enclose the staged mesh at its exact radial and numerical settings")
    dims,_,peak,operations = native_cost_chain(root,binding)
    count = (len(x)-1)*(len(z)-1)
    scale = count/math.prod(dims)
    measured, attempts, reports = [], {}, []
    for i,link in enumerate(reversed(operations)):
        receipt = read(verified(link))
        session = Path(receipt["request"]["command"][-1])
        path = session.parent/"state.steady.json"
        report = read(path)
        steps = report["accepted_newton_steps"]
        require(steps in (1,2) and report["best_from_accepted_iterate"],
                "timing anchor must expose short useful Newton segments")
        wall = receipt["result"]["seconds"]
        cap = math.ceil(1.22*wall*max(1.,scale))
        # Allow the accepted packet and native inventories to finish writing
        # after the solver's cooperative Newton deadline.
        attempts[case_dir(root,case,i).name] = planned_attempt(cap,cap+244,steps)
        measured.append(dict(steps=steps,seconds=wall,peak_tree_rss_bytes=receipt["result"]["peak_tree_rss_bytes"]))
        reports.append(identity(path))
        data["disk"] = max(data["disk"],math.ceil(1.10*3*size(session.parent)*scale))
    data["peak"] = max(data["peak"],math.ceil(peak*scale*data["margin"]))
    data["anchors"].extend([identity(audit_path),*operations,*reports])
    return attempts,dict(anchor=anchor,measured=measured,mesh=binding["mesh"],
                         cell_factor=scale,watch_gap_seconds=244)


def freeze(root, design_path, out):
    design = read(design_path)
    case, base = design["case"], design["base"]
    diagnostic = None
    if design.get("equation_diagnostic"):
        diagnostic = diagnostic_anchor(root,base,case,design["equation_diagnostic"])
        prior,audit_receipt=diagnostic["binding"],diagnostic["audit"]
    else:
        base_audit,audit_receipt=pass_audit(root,base)
        prior=base_audit["provenance"]
    prior_mesh, prior_state = verified(prior["mesh"]), verified(prior["state"])
    x,z,_ = mesh_axes(prior_mesh)
    new_x,new_z,columns,quadrature,tolerances = axes_for(root,design,x,z)
    cells = [len(new_x)-1,len(new_z)-1]
    for comparator in design.get("comparison_bases",[base]):
        if not diagnostic or comparator!=base:
            pass_audit(root,comparator)
    data = measure(root, ["Bn2","Bf",base,*design.get("extra_cost_cases",[])],
                   cells,quadrature,diagnostic)
    retry = None
    if design.get("readmit_after_rss"):
        retry = readmit_memory_stop(root,design["readmit_after_rss"],base,new_x,new_z,
                                    quadrature,tolerances,data)
    elif design.get("scaled_from_rss_stop"):
        retry = scale_radius_after_stop(root,design["scaled_from_rss_stop"],base,new_x,new_z,
                                        quadrature,tolerances,data)
    attempts = {case_dir(root,case,i).name:planned_attempt(seconds,seconds+115,n)
                for i,(seconds,n) in enumerate(((data["first"],1),(data["continuation"],50),
                                                 (data["retry"],50)))}
    schedule = None
    if design.get("streamwise_timing_anchor"):
        attempts,schedule = streamwise_schedule(root,case,design["streamwise_timing_anchor"],
                                                new_x,new_z,quadrature,tolerances,data)
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
        readmission=retry,equation_diagnostic=diagnostic,segmented_timing=schedule,
        anchors=data["anchors"],attempts=attempts,audit_cap_seconds=audit_seconds,
        compare_cap_seconds=compare_seconds,required_seconds=required,reservation_caps=caps,
        diagnostic_cap_seconds=max(0,extra_analysis-GRACE),live=live,memory=memory,
        storage=live_storage,forecast_additional_science_bytes=data["disk"],
        status="DEFERRED" if blockers else "FITS",blockers=blockers))
    print(case,base,cells,"seconds",required,"RSS GiB",data["peak"]/GIB,blockers,flush=True)


def resume_time(root, admission_path, out):
    initial = read(admission_path)
    case = initial["case"]
    first,later = (case_dir(root,case,i) for i in (0,1))
    o = root/"operations"
    start_path = o/f"{first.name}-steady"/"receipt.json"
    late_path = o/f"{later.name}-steady"/"receipt.json"
    start,late = (read(p)["result"] for p in (start_path,late_path))
    earlier_state = first/"state.steady-candidate.json"
    earlier = read(first/"state.steady.json")
    require(case=="Bx1" and initial["status"]=="ADMITTED" and late["reason"]=="TIME_LIMIT"
            and start["reason"] is None and earlier["best_from_accepted_iterate"]
            and earlier["accepted_newton_steps"]==1 and earlier_state.is_file()
            and not (later/"state.steady-candidate.json").exists()
            and not (later/"state.result.json").exists()
            and read(later/"inputs.json")["initial_guess"]==identity(earlier_state)
            and read(later/"inputs.json")["admission"]==identity(admission_path)
            and late["peak_tree_rss_bytes"]<initial["rss_limit_bytes"],
            "streamwise readmission requires an intact earlier checkpoint and a bound wall stop")
    first_duration=math.ceil(min(2.9*start["seconds"],.85*late["seconds"]))
    second_duration=math.ceil(2.2*start["seconds"])
    backup_duration=math.ceil(1.3*start["seconds"])
    serial_padding=max(225,math.ceil(late["seconds"]/6))
    attempts={case_dir(root,case,i).name:planned_attempt(seconds,seconds+serial_padding,steps)
              for i,(seconds,steps) in enumerate(((first_duration,2),(second_duration,50),
                                                     (backup_duration,50)),start=2)}
    caps=dict(initial["reservation_caps"])
    caps["preparations_and_native"]=sum(a["preparation_seconds"]+a["watch_seconds"]+
                                         2*GRACE for a in attempts.values())
    required=sum(caps.values())
    live=balance(read(root/"accounting.json"))
    memory=memory_snapshot()
    disk=storage(root)
    blockers=[k for k,ok in (
        ("wall_time",required<=live["nonprotected_remaining"]),
        ("host_or_cgroup",memory_headroom(memory)>=initial["rss_limit_bytes"]+2*GIB),
        ("science_storage",disk["science_bytes"]+initial["forecast_additional_science_bytes"]<=8*GIB),
        ("free_space",disk["free_bytes"]-initial["forecast_additional_science_bytes"]>=15*GIB)) if not ok]
    atomic(out,dict(initial,schema="WJH27FlatQualificationForecastV1",
        status="DEFERRED" if blockers else "FITS",blockers=blockers,
        prior_admission=identity(admission_path),time_stops=[identity(late_path)],
        recovery_checkpoint=identity(earlier_state),recovery_first_label=case_dir(root,case,2).name,
        recovery_residual=earlier["true_norm"],early_timing_anchor=identity(start_path),
        extra_watch_padding_seconds=serial_padding,attempts=attempts,reservation_caps=caps,
        required_seconds=required,live=live,memory=memory,storage=disk))
    print(case,"bound readmission from first Newton",required,"seconds",blockers,flush=True)


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
    for mode in ("freeze","pin","admit","resume-time"):
        task = modes.add_parser(mode)
        task.add_argument("input",type=Path)
        task.add_argument("output",type=Path)
    args=parser.parse_args()
    {"freeze":freeze,"pin":pin,"admit":admit,"resume-time":resume_time}[args.mode](args.root,args.input,args.output)


if __name__ == "__main__":
    main()
