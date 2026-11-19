#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently audit the outlet matrix and compare its signed owner observations."""

import argparse
import gzip
import json
import os
from pathlib import Path
import time

import numpy as np

import physical_oracle as po
from resources import atomic, identity
from task2_outlet_baseline import mesh_axes, read, verified
from task2_outlet_case import CASES, FULL, TRANSVERSE, accepted_directory, get_context
from task2r_verify import projection, require
from task2s_local_audit import continuity, weak_traces
from task2s_observations import compare, orientation_audit, check_archived


INTERIOR = (-2.5, 8.)
NEAR_RADII = (.0025,.003,.0035,.004,.0045,.005,.0055,.006)
FAR_RADII = (.0125,.025,.05,.075)
NEAR_ANGLES = (12,18,24,30,36,42,48)
FAR_ANGLES = (15,30,45,60,75)
ORIGINAL_PAIR = "task2rs-s2-nested-finish/nested-comparison.json"


def bind_elements(mesh, data):
    xs, zs, geometry = mesh_axes(mesh)
    vertices = {int(v.attrib["ID"]):tuple(map(float,v.text.split()))[:2]
                for v in geometry.find("VERTEX")}
    edges = {int(e.attrib["ID"]):tuple(map(int,e.text.split())) for e in geometry.find("EDGE")}
    expected = {}
    d0 = data["geometry"]["d0"]
    for e in geometry.find("ELEMENT"):
        verts = {v for edge in map(int,e.text.split()) for v in edges[edge]}
        x,z = zip(*(vertices[v] for v in verts))
        expected[int(e.attrib["ID"])] = (min(x),max(x),min(z)-d0,max(z)-d0)
    actual = {e["geometry"]:e["bounds"] for e in data["elements"]}
    require(actual.keys() == expected.keys() and all(
        max(abs(x-y) for x,y in zip(actual[key],value)) <= 1e-12 for key,value in expected.items()),
        "accepted native element IDs or owner coordinates differ from the bound geometry")
    return [len(xs)-1,len(zs)-1]


def open_case(root, name, context, ref, order, deadline):
    baseline, runtime, table, _ = context
    mesh, policy = name, FULL
    if name in ("A0", "A1"):
        archive = baseline["cases"][name]
        proof = read(root/"evidence/outlet-archive-parity.json")
        row = proof["cases"][name]
        require(proof["status"] == "PASS" and proof["new_application"] == runtime["application"]
                and proof["old_application"] == baseline["application"] and proof["reference"] == baseline["table"]
                and row["original"] == archive["state"] and row["original_signature"] == archive["signature"]
                and row["original_coefficient_sha256"] == archive["coefficient_sha256"]
                and row["original_source"] == baseline["source_identity"]
                and row["reproduced_source"] == runtime["source"], "archived A provenance differs")
        parity = verified(row["reproduced"])
        copied, _ = po.load_native(parity)
        require(copied["signature"] == row["reproduced_signature"] and copied["outlet_policy"] == FULL,
                "reproduced A binding differs")
        path, source = verified(archive["state"]), baseline["source_identity"]
        witnesses = dict(state=archive["state"], mesh=archive["geometry"], executable=baseline["application"],
                         parity=identity(root/"evidence/outlet-archive-parity.json"))
    else:
        path = accepted_directory(root,name)/"state.result.json"
        case = path.parent
        manifest = read(case/"inputs.json")
        _, mesh, policy = CASES[name]
        expected_geometry = (baseline["cases"][mesh]["geometry"] if mesh != "M1c" else
                             read((root/"science/outlet-a1c/inputs.json"))["geometry"])
        require(manifest["case"] == name and manifest["policy"] == policy and
                manifest["source"] == runtime["source"] and manifest["application"] == runtime["application"] and
                manifest["table"] == baseline["table"] and manifest["geometry"] == expected_geometry and
                manifest["modes"] == [4,3] and manifest["quadrature"] == [6,8] and
                manifest["tolerances"] == [1e-12,1e-10] and manifest["steady_time"] == 0.,
                "new outlet case binding differs")
        receipt = root/"operations"/(case.name+"-steady")/"receipt.json"
        evidence = read(receipt)
        require(evidence["request"]["command"] == manifest["command"]
                and evidence["result"]["status"] == "PASS" and
                all(manifest[key] in evidence["request"]["input_files"] for key in ("application","session","options")),
                "new case lacks the exact registered run receipt")
        result_file = identity(path)
        if "output_files" in evidence["result"]:
            require(result_file in evidence["result"]["output_files"],"accepted output digest differs from the receipt")
        for key in ("geometry","session","options","initial_guess","table"):
            verified(manifest[key])
        source = runtime["source"]
        witnesses = dict(state=result_file, mesh=manifest["geometry"], executable=runtime["application"],
                         inputs=identity(case/"inputs.json"), operation=identity(receipt))
    data, meta = po.load_native(path)
    require(data["source_identity"] == source and data.get("outlet_policy", FULL) == policy
            and data["step"] == data["time"] == 0. and meta["steady_mode"] and
            meta["steady_solver"]["converged"] and max(data["residuals"]) <= 1e-10 and
            data["basis"] == "Modified_A" and data["constitutive_backend"] == "FastFull3DV1" and
            data["stress_split"] == "PositiveShearV1" and
            [data["quadratureX"],data["quadratureR"]] == [6,8] and
            [data["linear_tolerance"],data["nonlinear_tolerance"]] == [1e-12,1e-10] and
            ref.identity == identity(table)["sha256"],
            "outlet comparison requires a converged, explicitly bound state")
    cells = bind_elements(verified(witnesses["mesh"]),data)
    if name not in ("A0","A1"):
        transfer = data["initial_field_transfer"]
        require(cells == manifest["cells"] and transfer["packet_sha256"] == manifest["initial_guess"]["sha256"]
                and transfer["source_state_sha256"] == manifest["initial_coefficient_sha256"],
                "accepted steady result differs from the signed mesh or typed guess inventory")
    state = po.PhysicalOracle(data,ref,order,deadline)
    witnesses.update(source=source, coefficient_sha256=data["state_sha256"], signature=data["signature"],
                     outlet_identity=state.outlet_identity, policy=policy, cells=cells)
    return state, witnesses


def raw_points(state, e, x, r):
    x,r = np.atleast_1d(x),np.atleast_1d(r)
    q = state.sample(e,x,r)
    return [dict(x=float(x[i]), r=float(r[i]), owner=e["geometry"], bounds=e["bounds"],
                 U=float(q["U"][i]), W=float(q["W"][i]), k=float(q["k"][i]),
                 E=float(q["E"][i]), E_type="transported" if q["outer"] else "dependent",
                 sampled_K=float(q["K"][i]), production=float(q["P"][i]), D=float(q["div"][i]),
                 Pi=float(q["Pi"][i]), R=q["R"][i].tolist(), gradient=q["G"][i].tolist(),
                 Vxr=float(q["V"][i,0,2])) for i in range(len(x))]


def off_vertex(state):
    groups = {}
    last = [e for e,side in state.faces["u"] if side==1]
    ground = [e for e,side in state.faces["u"] if side==2]
    near = [e for e in state.elements if e["bounds"][1]>=15.-max(FAR_RADII)
            and e["bounds"][2]<=max(FAR_RADII)]

    def at(records, x, r, owners, **metadata):
        matched = [e for e in owners if e["bounds"][0]-1e-12<=x<=e["bounds"][1]+1e-12
                   and e["bounds"][2]-1e-12<=r<=e["bounds"][3]+1e-12]
        require(matched,"fixed boundary observation has no owner")
        for e in matched:
            records.append(dict(**raw_points(state,e,x,r)[0],**metadata))

    for name,radii,angles in (("near_annulus",NEAR_RADII,NEAR_ANGLES),
                               ("broad_annulus",FAR_RADII,FAR_ANGLES)):
        records = []
        for radius in radii:
            for angle in angles:
                radians = np.deg2rad(angle)
                at(records,15.-radius*np.cos(radians),radius*np.sin(radians),
                   near, angular_degrees=angle, corner_distance_code=radius)
        groups[name] = records
    for name,owners in (("outlet",last),("ground",ground)):
        records = []
        for radius in NEAR_RADII+FAR_RADII:
            x,r = (15.,radius) if name=="outlet" else (15.-radius,0.)
            at(records,x,r,owners,corner_distance_code=radius)
        groups[name] = records
    return {name:dict(samples=records, owner_count=len({row["owner"] for row in records}),
                       D_sample_RMS=float(np.sqrt(np.mean([row["D"]**2 for row in records]))),
                       D_sample_MAX=max(abs(row["D"]) for row in records))
            for name,records in groups.items()}


def audit(root,name,out,seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    state,binding = open_case(root,name,context,ref,24,time.monotonic()+seconds)
    projection_file = out.with_name(out.stem+"-inlet.json")
    projection(context[2],Path(binding["state"]["path"]),projection_file)
    core = state.assemble()
    boundary = state.boundary_checks()
    interfaces = state.interfaces()
    profiles = state.profiles()
    weak,shift = continuity(state)
    traces = weak_traces(state,shift)
    if state.outlet_policy == TRANSVERSE:
        vertex = traces["raw_vertex"]
        require(abs(vertex["D"]) <= 1e-10 and max(abs(vertex["G"][i][i]) for i in (0,1)) <= 1e-10,
                "ground/outlet zero divergence trace identity differs")
    shared = interfaces["norms"]
    owner_flux_pass = max(shared[f"{zone}/{q}"]["MAX"] for zone in ("H","interior")
                          for q in ("mass_shared","k_advection_shared")) <= 1e-10
    guards = {key: state.data[key] for key in ("minimum_k","minimum_epsilon","minimum_covariance_over_k",
                                             "minimum_principal_over_mu")}
    limits = {
        "independent_current":core["steady_residual_status"], "weak_accounting":core["accounting_identity_status"],
        "boundary":boundary["status"], "shared_interface_flux":"MET" if owner_flux_pass else "NOT_MET",
        "physical":"MET" if profiles["status"] == "MET" and
                       all(v["status"]=="MET" for v in core["physical_envelope"].values()) else "NOT_MET",
        "inlet_projection":"MET" if read(projection_file)["status"]=="PASS" else "NOT_MET",
        "native_guards":"MET" if min(guards.values())>0 and guards["minimum_covariance_over_k"]>1e-10
                       and guards["minimum_principal_over_mu"]>1e-8 else "NOT_MET",
    }
    atomic(out,dict(schema="WJH27OutletPhysicalAuditV1",name=name, provenance=binding,
                    core=core,boundary=boundary,interfaces=interfaces,
                    profile_norms=profiles["norms"],ground_stations=profiles["ground_stations"],
                    retained_physical_stations=profiles["stations"],fig3_source=profiles["fig3_source"],
                    physical_7mm=profiles["physical_7mm"], continuity=weak, traces=traces,
                    off_vertex=off_vertex(state), inlet_projection=identity(projection_file),
                    native_guards=guards, status=limits))
    print(name,limits,flush=True)


def write_station_witness(root,name,out,seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    state,binding = open_case(root,name,context,ref,24,time.monotonic()+seconds)
    profiles = state.profiles(detailed=True)
    weak,_ = continuity(state)
    old_audit = root/"evidence"/f"outlet-{name}-audit.json"
    require(read(old_audit)["provenance"]["state"] == binding["state"],"station and current audit states differ")
    payload = dict(schema="WJH27OutletStationWitnessV1", case=name, provenance=binding,
                   station_count=len(profiles["stations"]), stations=profiles["stations"],
                   owner_rows=profiles["rows"], wall_shear=profiles["ground_stations"],
                   continuity=weak, current_audit=identity(old_audit))
    require(not out.exists(),"station witness already exists")
    temporary = out.with_name(out.name+".tmp")
    with gzip.open(temporary,"xt",compresslevel=4) as stream:
        json.dump(payload,stream,separators=(",",":"),allow_nan=False)
    os.replace(temporary,out)
    print(name,len(profiles["rows"]),weak["pressure_tested_Riesz_RMS"],flush=True)


def archive_witness(root):
    published = read(root/"tools/outlet-control-v1/summary.json")
    return verified(published["comparisons"][ORIGINAL_PAIR]["witness"])


def pair(first,second,bindings,replay=None):
    rows,interior,sampling = compare(first,second,interior=INTERIOR,reverse=True,detailed=True)
    result = dict(first=bindings[0],second=bindings[1],difference="second minus first",
                  rows=rows,interior_rows=interior,interior_x=list(INTERIOR), sampling=sampling,
                  failures=[k for k,v in rows.items() if v["status"] != "MET"],
                  four_corner_rows={k:rows[f"x=15/whole/{k}"] for k in ("R","Pi","traction","G")})
    if replay:
        result["archived_A_replay"] = check_archived(rows,read(replay))
        require(result["archived_A_replay"]["status"]=="PASS","A pair did not replay archived spatial rows")
    return result


def station_uncertainty(spatial_a,spatial_b,effect_0,effect_1,
                        differences=("A1_minus_A0","B1_minus_B0","B0_minus_A0","B1_minus_A1")):
    rows = dict(effect_1["interior_rows"])
    rows.update({k:v for k,v in effect_1["rows"].items() if k.startswith("x=")})
    results = {}
    for name,row in rows.items():
        kind = "interior_rows" if name.startswith("interior/") else "rows"
        a,b,e0,e1 = [item[kind][name] for item in (spatial_a,spatial_b,effect_0,effect_1)]
        data = {}
        for i,comp in enumerate(row["components"]):
            data[comp] = {}
            for norm in ("RMS","MAX"):
                field = "component_"+norm
                aa,bb,d0,d1 = [item[field][i] for item in (a,b,e0,e1)]
                mean = [e["component_signed_mean"][i] for e in (e0,e1)]
                signed = mean if norm=="RMS" else [e["component_locations"][i]["signed_difference"] for e in (e0,e1)]
                resolved = d0>aa+bb and d1>aa+bb and signed[0]*signed[1]>0.
                data[comp][norm] = dict(zip(differences,(aa,bb,d0,d1)),
                                       uncertainty_indicator=aa+bb, status="RESOLVED" if resolved else "UNRESOLVED",
                                       effect_directions=signed)
        results[name] = data
    return results


def compare_matrix(root,out,seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    deadline = time.monotonic()+seconds
    states = {n:open_case(root,n,context,ref,20,deadline) for n in ("A0","A1","B0","B1")}
    names = (("A0","A1"),("B0","B1"),("A0","B0"),("A1","B1"))
    pairs = {}
    for a,b in names:
        first,s1 = states[a]
        second,s2 = states[b]
        key = b+"-"+a
        pairs[key] = pair(first,second,(s1,s2),archive_witness(root) if (a,b)==names[0] else None)
        print(key,len(pairs[key]["failures"]),flush=True)
    atomic(out,dict(schema="WJH27OutletMatrixV1", comparisons=pairs,
        station_uncertainty=station_uncertainty(*(pairs[b+"-"+a] for a,b in names)),
        measurement_inventory=po.STATIONS, signed_basis={n:orientation_audit(v[0]) for n,v in states.items()}))


def collar(root,base,fine,out,seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    end = time.monotonic()+seconds
    a,pa = open_case(root,base,context,ref,20,end)
    b,pb = open_case(root,fine,context,ref,20,end)
    result = pair(a,b,(pa,pb))
    lo,mid,hi = context[0]["M1c"]["original_collar_owner_r"]
    radii = ((lo+mid)/2,mid,(hi+mid)/2)
    samples = {}
    for name,state in ((base,a),(fine,b)):
        relevant = [e for e in state.elements if e["bounds"][2]<hi-1e-12 and e["bounds"][3]>lo+1e-12]
        samples[name] = [raw_points(state,e,x,r)[0] for x in po.STATIONS for r in radii for e in relevant
                         if e["bounds"][0]-1e-12<=x<=e["bounds"][1]+1e-12
                         and e["bounds"][2]-1e-12<=r<=e["bounds"][3]+1e-12]
    atomic(out,dict(schema="WJH27OutletCollarControlV1",case_pair=[base,fine],comparison=result,
                    original_owner_radial_edges=[lo,mid,hi],samples=samples))
    print(fine+"-"+base,len(result["failures"]),flush=True)


def collar_effect(root,out,seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    deadline = time.monotonic()+seconds
    a,pa = open_case(root,"A1c",context,ref,20,deadline)
    b,pb = open_case(root,"B1c",context,ref,20,deadline)
    require(pa["mesh"]==pb["mesh"],"A and B collar effects require precisely the same mesh")
    paths = [root/"evidence"/name for name in ("outlet-a-collar.json","outlet-b-collar.json","outlet-matrix.json")]
    a_refine,b_refine = (read(path)["comparison"] for path in paths[:2])
    effect_m1 = read(paths[2])["comparisons"]["B1-A1"]
    require(a_refine["second"]["state"]==pa["state"] and b_refine["second"]["state"]==pb["state"]
            and a_refine["first"]["state"]==effect_m1["first"]["state"]
            and b_refine["first"]["state"]==effect_m1["second"]["state"],
            "original and radial indicators do not bind the same paired states")
    effect = pair(a,b,(pa,pb))
    atomic(out,dict(schema="WJH27OutletCollarBoundaryEffectV1",effect=effect,
          reference_comparisons=[identity(path) for path in paths],
          station_uncertainty=station_uncertainty(a_refine,b_refine,effect_m1,effect,
                ("A1c_minus_A1","B1c_minus_B1","B1_minus_A1","B1c_minus_A1c"))))
    print("published radial uncertainty on all 29 stations and the fixed interior",flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("mode",choices=("audit","matrix","collar","collar-effect","station-witness"))
    parser.add_argument("out",type=Path)
    parser.add_argument("--case")
    parser.add_argument("--seconds",type=float,required=True)
    args = parser.parse_args()
    if args.mode=="audit":
        audit(args.root,args.case,args.out,args.seconds)
    elif args.mode=="station-witness":
        write_station_witness(args.root,args.case,args.out,args.seconds)
    elif args.mode=="matrix":
        compare_matrix(args.root,args.out,args.seconds)
    elif args.mode=="collar-effect":
        collar_effect(args.root,args.out,args.seconds)
    else:
        require(args.case in ("A","B"),"collar comparison needs a policy")
        collar(args.root,args.case+"1",args.case+"1c",args.out,args.seconds)
