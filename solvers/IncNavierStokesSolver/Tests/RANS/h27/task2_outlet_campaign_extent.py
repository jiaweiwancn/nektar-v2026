#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze, audit, and compare a matched extension of the qualified flat recipe."""

import argparse
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

import physical_oracle as po
from resources import THREADS, atomic, identity, size
from task2_outlet_analysis import INTERIOR
from task2_outlet_baseline import mesh_axes, read, source_identity, verified
from task2_outlet_campaign import measure
from task2_outlet_case import get_context
from task2_outlet_qualification import GRACE, PUBLICATION_SECONDS, case_dir, load_case, pass_audit, planned_attempt
from task2r_run import (CAMPAIGN_AMENDMENT, CAMPAIGN_RSS_GIB, EXTENT_BINDING, EXTENT_RUNTIME,
    GIB, balance, extent_case, memory_headroom, memory_snapshot, require, storage)
from task2s_observations import overlapping

EXTENT_PREFIX = "outlet-campaign-extent"
TERMINAL_WIDTH = .625
FROZEN_BREAKS = [14.6875,15.,15.375,15.75,16.125,16.5,16.875]


def extend(x, exit_x):
    old = x[-1]
    start = old-TERMINAL_WIDTH
    require(start in x and math.isfinite(exit_x) and exit_x>old+TERMINAL_WIDTH,
            "extent needs the inherited terminal grading and a farther exit")
    shift = exit_x-old
    if old==15. and exit_x==17.5:
        inserted=FROZEN_BREAKS
    else:
        count=math.ceil(shift/.375)
        inserted=sorted({old,*(start+i*shift/count for i in range(1,count))})
    new = [a for a in x if a<=start]+inserted+[a+shift for a in x if a>start]
    require(new==sorted(set(new)) and old in new and new[-1]==exit_x,
            "extent axes do not retain the original section")
    return new,start


def forecast(root, selection, case, exit_x, out):
    require(extent_case(case) and case[1]=="e", "expected a new extent case")
    proof=read(selection)
    require(proof["canonical_qualified"] and proof["spatial_status"]=="MET",
            "extension requires the selected audited flat recipe")
    name=proof["selected_case"]
    state=pass_audit(root,name)[0]["provenance"]
    require(state==proof["canonical"], "selected state and independent audit differ")
    x,z,_=mesh_axes(verified(state["mesh"]))
    new_x,start=extend(x,exit_x)
    q,t=state["quadrature"],state["tolerances"]
    cells=[len(new_x)-1,len(z)-1]
    tracked=measure(root,["Bn2","Bf",name],cells,q)
    limit=max(10,math.ceil(tracked["peak"]/GIB))*GIB
    live=balance(read(root/"accounting.json"))
    memory=memory_snapshot()
    disk=storage(root)
    anchors={}
    for label in ("outlet-source-stage","outlet-configure","outlet-native-build",
                  "outlet-a1-replay-steady","outlet-qual-bf-audit","outlet-qual-bf-compare"):
        receipt=root/"operations"/label/"receipt.json"
        result=read(receipt)["result"]
        require(result["status"]=="PASS","the matched extent forecast requires a passing resource anchor")
        anchors[label]=(identity(receipt),result["seconds"])
    projected_build=math.ceil(1.5*anchors["outlet-native-build"][1])
    parity_attempt=planned_attempt(110,250,2)
    parity_audit,parity_compare=265,160
    parity_native=parity_attempt["preparation_seconds"]+parity_attempt["watch_seconds"]+2*GRACE
    stages=dict(setup=45,configure=32,build=projected_build+GRACE,binding=32,
        parity_phase_forecast=32,parity_admission=27,parity_native_and_preparation=parity_native,
        parity_audit_and_comparison=parity_audit+parity_compare+2*GRACE,
        parity_verification=32,extent_phase_forecast=32,extent_native_and_preparations=0,extent_admission=27,
        extent_audit=math.ceil(280*math.prod(cells)/(111*92)),
        extent_comparison=math.ceil(175*math.prod(cells)/(111*92)),
        extent_uncertainty=180,
        checkpoint_publication=100,final_evidence=PUBLICATION_SECONDS)
    native_seconds=max(tracked["continuation"],math.ceil(2.6*tracked["first"]))
    attempts={case_dir(root,case,i).name:planned_attempt(native_seconds,native_seconds+244,2)
              for i in range(2)}
    science_disk=math.ceil(tracked["disk"]*len(attempts)/3)+180_000_000
    stages["extent_native_and_preparations"]=sum(value["preparation_seconds"]+value["watch_seconds"]+
                                                 2*GRACE for value in attempts.values())
    stages["extent_audit"]+=GRACE
    stages["extent_comparison"]+=GRACE
    cost=sum(stages.values())
    build_disk=math.ceil(1.25*sum(size(root/name) for name in (
        "source-outlet-v1","build-outlet-v1","runtime-outlet-v1")))
    blockers=[key for key,ok in (
        ("wall_time",cost<=live["nonprotected_remaining"]),
        ("named_RSS",limit<=CAMPAIGN_RSS_GIB*GIB),
        ("host_or_cgroup",memory_headroom(memory)>=limit+2*GIB),
        ("science_storage",disk["science_bytes"]+science_disk<=8*GIB),
        ("build_storage",disk["source_build_runtime_bytes"]+build_disk<=8*GIB),
        ("free_space",disk["free_bytes"]-science_disk-build_disk>=15*GIB)) if not ok]
    limits=read(root/"evidence/outlet-location-eligibility.json")["intended_stationary_test"]["physical_RMS_MAX_limits"]
    atomic(out,dict(schema="WJH27OutletExtentCampaignV1",
        status="DEFERRED" if blockers else "FITS",blockers=blockers,
        selected_case=name,canonical=state,selection=identity(selection),
        canonical_comparisons=proof["canonical_comparisons"],case=case,cells=cells,
        x_edges=new_x,physical_z_edges=z,original_x_edges=x,former_outlet=x[-1],outlet_x=exit_x,
        retained_upstream_through=start,stations=sorted(set([*po.STATIONS,exit_x])),
        quadrature=q,tolerances=t,attempts=attempts,required_seconds=cost,reservation_caps=stages,
        forecast_peak_tree_rss_bytes=tracked["peak"],rss_limit_bytes=limit,
        footprint=dict(science=science_disk,build=build_disk),
        parity_schedule=dict(native=parity_attempt,audit=parity_audit,comparison=parity_compare),
        anchors=[*tracked["anchors"],*(link for link,_ in anchors.values())],
        live=live,memory=memory,storage=disk,limits=limits,experimental_x=INTERIOR))
    print(case,cells,cost,tracked["peak"]/GIB,blockers,flush=True)


def seal(root, out):
    source,files=source_identity(root/"source-outlet-campaign-extent")
    path=root/EXTENT_RUNTIME/"IncNavierStokesSolver"
    compilation=root/"operations/outlet-campaign-extent-build/receipt.json"
    require(read(compilation)["result"]["status"]=="PASS" and source.encode() in path.read_bytes(),
            "compiled extent source differs from the archived inventory")
    atomic(out,dict(source=source, application=identity(path),compiled_inventory=files,
                   compile_receipt=identity(compilation)))


def phase_forecast(root, design_path, case, out, parity_path):
    design=read(design_path)
    require(design["status"]=="FITS" and extent_case(case), "extent stage lacks complete admission")
    ev=root/"evidence"
    bind=ev/EXTENT_BINDING
    runtime=read(bind)
    canonical=design["canonical"]
    cap=design["reservation_caps"]
    is_parity=case[1]=="p"
    if is_parity:
        seed=canonical
        x,exit_x=design["original_x_edges"],design["former_outlet"]
        timings=design["parity_schedule"]
        attempts={case_dir(root,case).name:timings["native"]}
        estimate,limit=8*GIB,10*GIB
        audit,comparison=timings["audit"],timings["comparison"]
        science_disk=180_000_000
        reservations=dict(admission=cap["parity_admission"],
            preparations_and_native=cap["parity_native_and_preparation"],
            independent_audit=audit+GRACE,direct_comparisons=comparison+GRACE,
            parity_verification=cap["parity_verification"])
    else:
        require(case==design["case"] and parity_path is not None, "unregistered extended case")
        parity=read(parity_path)
        require(parity["status"]=="PASS" and parity["original"]==canonical["state"]
                and parity["application"]==runtime["application"], "x=15 parity proof differs")
        seed=parity["replay"]
        x,exit_x=design["x_edges"],design["outlet_x"]
        attempts=design["attempts"]
        estimate,limit=design["forecast_peak_tree_rss_bytes"],design["rss_limit_bytes"]
        audit,comparison=cap["extent_audit"]-GRACE,cap["extent_comparison"]-GRACE
        science_disk=design["footprint"]["science"]-180_000_000
        reservations=dict(admission=cap["extent_admission"],
            preparations_and_native=cap["extent_native_and_preparations"],
            independent_audit=cap["extent_audit"],direct_comparisons=cap["extent_comparison"],
            extent_uncertainty=cap["extent_uncertainty"])
    live=balance(read(root/"accounting.json"))
    memory=memory_snapshot()
    disk=storage(root)
    reservations.update(checkpoint_publication=cap["checkpoint_publication"],final_evidence=PUBLICATION_SECONDS)
    required=sum(reservations.values())
    require(required<=live["nonprotected_remaining"] and memory_headroom(memory)>=limit+2*GIB
            and disk["science_bytes"]+science_disk<=8*GIB and disk["free_bytes"]-science_disk>=15*GIB
            and disk["source_build_runtime_bytes"]<=8*GIB,
            "extent resources changed before the named admission")
    atomic(out,dict(schema="WJH27FlatQualificationForecastV1",status="FITS",blockers=[],
        case=case,base=design["selected_case"],comparison_bases=[design["selected_case"]],
        design=identity(design_path),extent_runtime=identity(bind),source=runtime["source"],
        application=runtime["application"],seed_source=seed["source"],
        prior_state=seed["state"],prior_mesh=seed["mesh"],outlet_x=exit_x,
        x_edges=x,physical_z_edges=design["physical_z_edges"],
        cells=[len(x)-1,len(design["physical_z_edges"])-1],
        forecast_peak_tree_rss_bytes=estimate,rss_limit_bytes=limit,anchors=design["anchors"],
        quadrature=design["quadrature"],tolerances=design["tolerances"],attempts=attempts,
        required_seconds=required,reservation_caps=reservations,
        diagnostic_cap_seconds=0,audit_cap_seconds=audit,compare_cap_seconds=comparison,
        forecast_additional_science_bytes=science_disk,storage=disk,
        campaign_amendment=identity(ev/f"{CAMPAIGN_AMENDMENT}.json"),live=live,memory=memory))


def verify_parity(root, case, design_path, out):
    design=read(design_path)
    original=design["canonical"]
    old,_=po.load_native(verified(original["state"]))
    audit,receipt=pass_audit(root,case)
    replay=audit["provenance"]
    new,meta=po.load_native(verified(replay["state"]))
    first=pass_audit(root,design["selected_case"])[0]
    report_path=root/"evidence"/f"outlet-qual-{case}-{design['selected_case']}.json"
    c=read(report_path)["comparison"]
    require(new["state_sha256"]==old["state_sha256"] and new["residuals"]==old["residuals"]
            and new["initial_field_transfer"]["packet_sha256"]==original["state"]["sha256"]
            and not new["initial_field_transfer"]["projected"] and meta["steady_solver"]["converged"]
            and replay["mesh"]==original["mesh"] and first["core"]["residuals"]==audit["core"]["residuals"]
            and c["second"]["state"]==replay["state"] and c["first"]["state"]==original["state"]
            and len(c["rows"])==690 and all(row["MAX"]==0. for row in c["rows"].values()),
            "native or evaluator x=15 parity differs")
    atomic(out,dict(status="PASS",original=original["state"],replay=replay,
        application=replay["executable"],audit=receipt,comparison=identity(report_path),
        coefficient_sha256=new["state_sha256"],native_residuals=new["residuals"],
        independent_residuals=audit["core"]["residuals"]))
    print(case,"x=15 native and independent parity PASS",flush=True)


def wall_effect(a,b):
    def ground(state):
        return {(e["bounds"][0],e["bounds"][1]):e for e,side in state.faces["u"] if side==2}
    left,right=ground(a),ground(b)
    xs,zs=sorted(left),sorted(right)
    t,w=leggauss(24)
    norms=po.Norms()
    peak=None
    for p,q in overlapping(xs,zs):
        lo,hi=max(p[0],q[0],INTERIOR[0]),min(p[1],q[1],INTERIOR[1])
        if hi-lo<=1e-12:
            continue
        x=np.r_[lo,(lo+hi)/2+(hi-lo)*t/2,hi]
        values=[state.sample(e,x,np.zeros(len(x)))["V"][:,0,2]
                for state,e in ((a,left[p]),(b,right[q]))]
        change=(values[1]-values[0])/a.ref.tau
        weights=np.r_[0.,w*(hi-lo)/2,0.]
        norms.add("interior","wall_shear",change,weights)
        i=int(np.argmax(np.abs(change)))
        if peak is None or abs(change[i])>abs(peak["signed"]):
            peak=dict(x=float(x[i]),signed=float(change[i]),x_owner_intervals=[p,q],
                      included_with_quadrature_weight=bool(weights[i]>0.))
    row=norms.report({})["interior/wall_shear"]
    require(abs(row["measure"]-(INTERIOR[1]-INTERIOR[0]))<1e-10, "wrong shared wall interval")
    return dict(row,peak=peak)


def window_departures(state):
    """Use the same physical strip on each tested extent."""
    norms = po.Norms()
    for owner in state.elements:
        left,right,bottom,top = owner["bounds"]
        lo,hi = max(left,INTERIOR[0]),min(right,INTERIOR[1])
        if hi-lo<=1e-12:
            continue
        x,r,weights = state.rule(owner,(lo,hi,bottom,top))
        for edge in (lo,hi):
            cuts = state.cuts(bottom,top,(state.ref.cutoff(
                state.trace(owner,"k",np.array([edge])))[0],state.ref.inlet_cutoff))
            x = np.r_[x,np.full(len(cuts),edge)]
            r = np.r_[r,cuts]
            weights = np.r_[weights,np.zeros(len(cuts))]
        values = state.sample(owner,x,r)
        state.inlet_norms(norms,"whole",values,r,weights)
        state.inlet_norms(norms,"outer" if values["outer"] else "inner",values,r,weights)
        if abs(bottom)>1e-12:
            continue
        nodes,weights = state.gauss
        xx = np.r_[lo,(lo+hi)/2+(hi-lo)*nodes/2,hi]
        shear = state.sample(owner,xx,np.zeros(len(xx)))["V"][:,0,2]
        wall_weights = np.r_[0.,(hi-lo)*weights/2,0.]
        norms.add("ground","shear_over_tau0",shear/state.ref.tau-1.,wall_weights)
        norms.add("ground","shear_minus_Tref_over_tau0",
                  (shear-state.ref.reference_shear)/state.ref.tau,wall_weights)
    limits = dict(po.ENVELOPES,shear_over_tau0=(.05,.10))
    rows = norms.report(limits)
    require(abs(rows["whole/U"]["measure"]-(INTERIOR[1]-INTERIOR[0])*state.ref.top)<1e-8
            and abs(rows["ground/shear_over_tau0"]["measure"]-(INTERIOR[1]-INTERIOR[0]))<1e-10,
            "independent inlet departures do not cover the fixed experimental strip")
    return rows


def assess(root, design_path, out, seconds):
    design=read(design_path)
    context=get_context(root)
    ref=po.load_reference(context[2])
    deadline=time.monotonic()+seconds
    name,case=design["selected_case"],design["case"]
    a,pa=load_case(root,name,context,ref,24,deadline)
    b,pb=load_case(root,case,context,ref,24,deadline)
    require(a.right==design["former_outlet"] and b.right==design["outlet_x"],
            "tested outlet locations differ")
    base_audit,current_audit=(pass_audit(root,n)[0] for n in (name,case))
    report_path = root/"evidence"/f"outlet-qual-{case}-{name}.json"
    report = read(report_path)
    comparison = report["comparison"]
    require(report["cases"]==[name,case] and comparison["first"]==pa and comparison["second"]==pb
            and len(comparison["rows"])==690 and comparison["interior_x"]==list(INTERIOR),
            "measured matched comparison does not bind both raw states and the fixed strip")
    rows=comparison["interior_rows"]
    anchors=[read(verified(item)) for item in design["canonical_comparisons"]]
    reference={}
    wall_reference=[]
    for witness in anchors:
        coarse,fine=witness["cases"]
        require(coarse==name and witness["comparison"]["first"]["state"]==pa["state"],
                "extent uncertainty requires direct comparisons with the selected state")
        f,p=load_case(root,fine,context,ref,24,deadline)
        require(p["state"]==witness["comparison"]["second"]["state"],
                "extent uncertainty witness differs")
        wall_reference.append(wall_effect(a,f))
    wall=wall_effect(a,b)
    for name_field in ("U","W","k","epsilon","R","Pi","traction","wall_shear"):
        key="interior/whole/"+name_field
        expected=design["limits"].get(name_field)
        effect=wall if name_field=="wall_shear" else rows[key]
        prior=wall_reference if name_field=="wall_shear" else [v["comparison"]["interior_rows"][key] for v in anchors]
        values={}
        for i,norm in enumerate(("RMS","MAX")):
            high=sum(v[norm] for v in prior)
            actual=effect[norm]
            values[norm]=dict(effect=actual,within_B_indicator=high,
                resolved_by_magnitude=actual>high,limit=expected[i] if expected else None,
                decision=("MET" if actual+high<=expected[i] else "NOT_MET" if actual-high>expected[i] else
                          "UNCERTAIN") if expected else "NO_PHYSICAL_THRESHOLD")
        reference[name_field]=values
    rows_passed=all(info["decision"]=="MET" for k in design["limits"] for info in reference[k].values())
    atomic(out,dict(schema="WJH27OutletLocationCampaignDecisionV1",
        status="ACCEPTABLE" if rows_passed else "OPEN",matched_comparison=comparison,
        matched_comparison_receipt=identity(report_path),
        fixed_region_limits=design["limits"],fixed_region_wall_shear=wall,
        achieved_within_B_uncertainty=reference,both_physical_departures={
            state.right:dict(fixed_window=window_departures(state),
                whole_domain=audit["core"]["physical_envelope"],station_norms=audit["profile_norms"],
                ground_stations=audit["ground_stations"],Tref=state.ref.reference_shear)
            for state,audit in ((a,base_audit),(b,current_audit))},
        new_outlet_and_former_station_profiles=b.profiles(stations=sorted({a.right,b.right})),
        reference_stations=design["stations"],source=identity(design_path),
        uncertainty_interpretation="Sum of directly measured within-canonical numerical witnesses. The interval brackets each physical acceptance threshold."))
    print(case,"ACCEPTABLE" if rows_passed else "OPEN",flush=True)


def handoff(root, design_path, assessment_path, out):
    design = read(design_path)
    selection = read(verified(design["selection"]))
    assessed = read(assessment_path)
    first,second = (assessed["matched_comparison"][part] for part in ("first","second"))
    parity_path = root/"evidence/outlet-campaign-x15-parity.json"
    parity = read(parity_path)
    require(selection["canonical_qualified"] and assessed["status"]=="ACCEPTABLE"
            and assessed["source"]==identity(design_path) and first==selection["canonical"]
            and parity["status"]=="PASS" and parity["original"]==first["state"]
            and parity["application"]==second["executable"]
            and all(row["decision"]=="MET" for field in design["limits"]
                    for row in assessed["achieved_within_B_uncertainty"][field].values()),
            "closure requires the passing qualified selection, exact parity and acceptable matched outlet")
    witnesses = {}
    for name,binding in ((design["selected_case"],first),(design["case"],second)):
        audit,receipt = pass_audit(root,name)
        require(audit["provenance"]==binding,"handoff audit and tested field differ")
        data,_ = po.load_native(verified(binding["state"]))
        require(data["time"]==data["step"]==0 and data["state_sha256"]==binding["coefficient_sha256"],
                "handoff requires the exact zero-time typed coefficients")
        manifest = read(verified(binding["inputs"]))
        witnesses[name] = dict(binding=binding,audit=receipt,table=manifest["table"],
            session=manifest["session"],options=manifest["options"],command=manifest["command"])
    reference = witnesses[design["selected_case"]]["table"]
    require(all(item["table"]==reference for item in witnesses.values()),"matched reference table differs")
    verified(reference)
    atomic(out,dict(schema="WJH27SteadyFlatOutletCampaignHandoffV1",campaign_status="CLOSED",
        recommended_case=design["selected_case"],recommended_outlet_x=design["former_outlet"],
        farther_tested_outlet_x=design["outlet_x"],shortest_tested_supported_extent=True,
        physical_time_advanced=0.,gate_matrix=dict(selection["gate_matrix"],outlet_location="MET"),
        qualified_selection=design["selection"],parity=identity(parity_path),
        extent_decision=identity(assessment_path),verify_inputs=witnesses,
        reference=reference,default_outlet_policy="FullTractionV1",
        linux_host="jay@100.88.70.60",experiment_root=str(root),
        linux_environment=str(root.parent.parent/"env.sh"),mpi_ranks=1,thread_settings=THREADS,
        worker_limit=80,protected_seconds=3600,cumulative_ceiling_seconds=108000,
        advancing_requirements=["physical timestep accuracy","cold restart","matched stock timing",
                                "complete advancing trajectory cost"],
        production_cost_ratio=1.25,production_wall_limit_seconds=86400,
        TASK2_COMPLETE="NO",TASK3_READY="NO"))
    print("Steady campaign CLOSED",design["selected_case"],design["former_outlet"],flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    action=parser.add_subparsers(dest="mode",required=True)
    for mode in ("forecast","seal","phase-forecast","verify-parity","assess","handoff"):
        p=action.add_parser(mode)
        p.add_argument("out",type=Path)
        if mode!="seal":
            p.add_argument("design",type=Path)
        if mode in ("forecast","phase-forecast","verify-parity"):
            p.add_argument("case")
        if mode=="forecast":
            p.add_argument("exit_x",type=float)
        if mode=="phase-forecast":
            p.add_argument("--parity",type=Path)
        if mode=="assess":
            p.add_argument("--seconds",type=float,required=True)
        if mode=="handoff":
            p.add_argument("assessment",type=Path)
    args=parser.parse_args()
    if args.mode=="forecast":
        forecast(args.root,args.design,args.case,args.exit_x,args.out)
    elif args.mode=="seal":
        seal(args.root,args.out)
    elif args.mode=="phase-forecast":
        phase_forecast(args.root,args.design,args.case,args.out,args.parity)
    elif args.mode=="verify-parity":
        verify_parity(args.root,args.case,args.design,args.out)
    elif args.mode=="assess":
        assess(args.root,args.design,args.out,args.seconds)
    else:
        handoff(args.root,args.design,args.assessment,args.out)
