#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bind direct spatial witnesses and equation controls to one flat candidate."""

import argparse
from pathlib import Path

import physical_oracle as po
from resources import atomic, identity
from task2_outlet_baseline import mesh_axes, read, verified
from task2_outlet_qual_followup import station_columns, native_chain
from task2_outlet_campaign import BLEND_SUPPORT
from task2_outlet_qualification import HIGHER_QUADRATURE, SHIFT, TIGHTER_TOLERANCES, TRANSVERSE, pass_audit
from task2r_run import balance, require


def axes(binding):
    x,z,_ = mesh_axes(verified(binding["mesh"]))
    return x,z


def inserted_in(coarse, fine, a, b):
    return a in coarse and b in coarse and any(a<x<b for x in fine)


def pair_proof(root, path, candidate, canonical, kind):
    document = read(path)
    first,second = document["cases"]
    require(first==candidate and first!=second, "witness must compare directly with the candidate")
    audit,_ = pass_audit(root,second)
    result = document["comparison"]
    binding = audit["provenance"]
    require(result["first"]==canonical and result["second"]==binding
            and all(binding[key]==canonical[key] for key in ("source","executable","outlet_identity","policy")),
            "comparison does not bind the named candidates")
    require(len(result["rows"])==690 and result["interior_x"]==[-2.5,8.]
            and not document["fixed_region_failures"] and not result["failures"],
            "signed comparison is missing required passing observations")
    for rows in (result["rows"],result["interior_rows"]):
        for key,row in rows.items():
            region,field = key.rsplit("/",2)[-2:]
            inner_control = kind in ("q","t") and (region,field)==("inner","k")
            allowed = 1e-4 if inner_control else po.NUMERICAL[field]
            require(row["limits"]==[allowed,allowed] and row["status"]=="MET"
                    and max(row["RMS"],row["MAX"])<=allowed,
                    f"{second} {key} differs from its declared acceptance rule")
    x,z = axes(canonical)
    fx,fz = axes(binding)
    if kind in ("q","t"):
        require(binding["mesh"]==canonical["mesh"], "numerical controls require the exact canonical mesh")
        q = HIGHER_QUADRATURE if kind=="q" else canonical["quadrature"]
        t = TIGHTER_TOLERANCES if kind=="t" else canonical["tolerances"]
    else:
        require(set(x)<=set(fx) and set(z)<=set(fz), "witness must preserve the candidate's physical owners")
        q,t=canonical["quadrature"],canonical["tolerances"]
    require(binding["quadrature"]==q and binding["tolerances"]==t,
            "witness changes an unauthorized numerical control")
    return dict(witness=second,comparison=identity(path),binding=binding,x=fx,z=fz,
                audit=identity(root/"evidence"/f"outlet-{second}-audit.json"))


def select(root, config, out):
    raw=read(config)
    ev=root/"evidence"
    candidate=raw["candidate"]
    best,audit=pass_audit(root,candidate)
    canonical=best["provenance"]
    verified(canonical["state"])
    require(canonical["policy"]==TRANSVERSE,"selected recipe must use the flat B boundary policy")
    x,z=axes(canonical)
    cases={}
    for group,files in raw["direct_comparisons"].items():
        cases[group]=[pair_proof(root,ev/name,candidate,canonical,group) for name in files]
    require(len(cases.get("q",[]))==len(cases.get("t",[]))==1, "both final numerical controls are required")
    radial=[(a,b) for a,b in zip(z,z[1:])
            if BLEND_SUPPORT[0]-1e-12<=a-SHIFT and b-SHIFT<=BLEND_SUPPORT[1]+1e-12]
    r_cases=cases.get("radial",[])
    require(r_cases and all(v["x"]==x for v in r_cases), "radial comparison must hold the x grid fixed")
    covered_r=[(a,b) for a,b in radial if any(inserted_in(z,v["z"],a,b) for v in r_cases)]
    require(covered_r==radial, "a declared radial owner was never independently challenged")

    experimental=[v for v in po.STATIONS if -2.5<=v<=8.]
    cols=station_columns(x,experimental) | set(zip(x[-3:-1],x[-2:]))
    if raw["include_inlet_column"]:
        cols.add((x[0],x[1]))
    x_cases=cases.get("streamwise",[])
    require(x_cases and all(v["z"]==z for v in x_cases), "streamwise comparison must hold the radial grid fixed")
    covered_x=[(a,b) for a,b in zip(x,x[1:]) if any(inserted_in(x,v["x"],a,b) for v in x_cases)]
    require(cols<=set(covered_x), "a declared experimental or terminal column remains unchanged")

    direction_path=ev/raw["direction_probe"]
    direction=read(direction_path)
    witnesses={candidate:canonical}
    witnesses.update({v["witness"]:v["binding"] for v in (*r_cases,*x_cases)})
    require(direction["schema"]=="WJH27OutletDirectionalOwnerProbeV1"
            and direction["bindings"][candidate]==canonical
            and all(witnesses.get(name)==proof for name,proof in direction["bindings"].items())
            and (not direction["coupled_follow_up_indicated"] or cases.get("coupled")),
            "directional owner evidence or the required coupled comparison is absent")

    path=ev/raw["H_proof"]
    h=read(path)
    require(h["status"]=="MET" and h["levels"][0]==candidate and h["bindings"][0]["state"]==canonical["state"]
            and h["direct_RMS_MAX_limit"]==1e-3 and h["successive_ratio_limit"]==.8
            and h["successive_resolved_floor"]==1e-5 and h["shared_MAX_limit"]==1e-10,
            "frozen H convergence rule or candidate does not match")
    h_states=[v["binding"] for v in cases.get("H",[])]
    require(h["bindings"][0]==canonical and len(h["bindings"])==3
            and all(binding in h_states for binding in h["bindings"][1:]),
            "H convergence needs direct from-candidate comparisons against both fine levels")

    sizes,native_wall,peak,anchors=native_chain(root,candidate)
    references=[case["comparison"] for group in cases.values() for case in group]
    atomic(out,dict(schema="WJH27FlatOutletSelectedRecipeV1", canonical_qualified=True,
        spatial_status="MET",selected_case=candidate,canonical=canonical,candidate_audit=audit,
        physical_status=best["status"]["physical"],current_status=best["status"]["independent_current"],
        quadrature=canonical["quadrature"],tolerances=canonical["tolerances"],
        gate_matrix=dict(radial="MET",streamwise="MET",H="MET",quadrature="MET",algebraic="MET"),
        quadrature_and_algebraic={"q":cases["q"][0]["comparison"],"t":cases["t"][0]["comparison"]},
        radial_owner_coverage=[[a-SHIFT,b-SHIFT] for a,b in covered_r],
        streamwise_owner_coverage=covered_x,required_streamwise_owners=sorted(cols),
        experimental_stations=experimental,H_proof=identity(path),direction_probe=identity(direction_path),
        canonical_comparisons=references,direct_comparisons=cases,
        measured_candidate_cost=dict(cells=sizes,successful_native_chain_seconds=native_wall,
            peak_tree_rss_bytes=peak,anchors=anchors),
        verified_checkpoints=raw["verified_checkpoints"],live=balance(read(root/"accounting.json")),
        outlet_location="PENDING_MATCHED_TEST",simulation_time_advanced=0.,
        TASK2_COMPLETE="NO",TASK3_READY="NO"))
    print(candidate,"passed all steady qualification gates",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("config",type=Path)
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    select(args.root,args.config,args.out)
