#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compare complete one-sided fields on frozen core, patches and coverage."""

import argparse
import csv
import json
import math
from pathlib import Path
import time

import numpy as np

from c1_oracle import Audit, coefficients, load_c1_state
from c1_paired_grid import APPROACHES, PARENT_NAMES, core_rows
from c1_support_controls import coincides
from c1_support_witness import (PROFILE_COLUMNS, FIXED_AXIAL, choose, extent, sampled_action_order,
                                section_points, write_profiles)
import physical_oracle as po
from physical_transient_oracle import state_hash
from resources import identity


GRADIENT = ("Gxx","Gxz","Gzx","Gzz")
TRACTIONS = ("traction_axial_x","traction_axial_z","traction_radial_x","traction_radial_z",
             "k_flux_axial","k_flux_radial")
FIELDS = ("margin","raw_response","mu",*GRADIENT,"gradient_frobenius","divergence",
          "u","w","k","epsilon","K","Pi_H","u_t_H","u_n_H","dPi_da","dK_da","dUt_da",
          "matching_area","pressure","k_load","convection","total_load","applied_total_load",*TRACTIONS,
          *("jump_"+key for key in TRACTIONS))
SCALED = {"gradient_frobenius":"G", **{k:k for k in ("pressure","k_load","convection","total_load")}}


def approached(records,H,domain):
    views = {}
    for xy,rows in records.items():
        x,r = xy
        a_left,a_right,top = domain
        q = {(sx,sr):choose(rows,sx,sr) for sx,sr in APPROACHES
             if not ((sx<0 and coincides(x,a_left,H)) or (sx>0 and coincides(x,a_right,H)) or
                     (sr<0 and coincides(r,0.,H)) or (sr>0 and coincides(r,top,H)))}
        for (sx,sr),owner in q.items():
            row = dict(owner,gradient_frobenius=math.hypot(*(owner[k] for k in GRADIENT)))
            for field in TRACTIONS:
                km,kp = (((-1,sr),(1,sr)) if "axial" in field else ((sx,-1),(sx,1)))
                # A boundary has a traction but no interior jump.
                row["jump_"+field] = q[kp][field]-q[km][field] if km in q and kp in q else 0.
            views[x,r,sx,sr] = row
    return views


def spec(row,value):
    return dict(value=value,owner=row["owner"],
                bounds=[row[k] for k in ("a_left","a_right","r_low","r_high")],
                cartesian=[row["x"],row["z"]])


def compare(a,b,scales,noise):
    if a.keys()!=b.keys():
        raise ValueError("approaching rays differ at the shared physical coordinates")
    maxima = {}
    max_position_error = 0.
    flips = common_negative = 0
    for key, u in a.items():
        v = b[key]
        error = math.hypot(u["x"]-v["x"],u["z"]-v["z"])
        max_position_error = max(error,max_position_error)
        if error>5e-12:
            raise ValueError(f"Cartesian coordinates disagree at {key}")
        ua,ub = (row["angular_difference"] for row in (u,v))
        flips += (u["margin"]<0.) != (v["margin"]<0.)
        common_negative += max(u["margin"]+ua,v["margin"]+ub)<0.
        for field in FIELDS:
            if field=="gradient_frobenius":
                diff = math.hypot(*(u[k]-v[k] for k in GRADIENT))
                va,vb = ([row[k] for k in GRADIENT] for row in (u,v))
            else:
                va,vb = u[field],v[field]
                diff = abs(va-vb)
            uncertainty = noise.get(field,0.)
            if field=="margin":
                uncertainty += ua+ub
            if field=="raw_response":
                uncertainty += ua*u["mu"]+ub*v["mu"]
            upper = diff+uncertainty
            absolute = maxima.get(field)
            bound = maxima.get("upper_"+field) if field in ("margin","raw_response") else None
            replace_absolute = absolute is None or diff>absolute["difference"]
            replace_upper = field in ("margin","raw_response") and (bound is None or upper>bound["upper"])
            if replace_absolute or replace_upper:
                x,r,sx,sr = key
                item = dict(difference=diff,evaluation_uncertainty=uncertainty,upper=upper,
                            chart=[x,r],approach=[sx,sr],a=spec(u,va),b=spec(v,vb))
                if field in SCALED:
                    item["normalized_difference"] = diff/scales[SCALED[field]]
                if replace_absolute:
                    maxima[field] = item
                if replace_upper:
                    maxima["upper_"+field] = item
    return dict(approach_samples=len(a),maxima=maxima,max_cartesian_position_difference=max_position_error,
                sampled_sign_changes=int(flips),common_negative_upper_samples=int(common_negative))


def historical_noise(profile, approached_views):
    old = core_rows(profile)
    noise = {key:0. for key in FIELDS}
    for old_name,new_name in zip(PARENT_NAMES,("A0","B0")):
        for (a,r,sx,sr),v in approached_views[new_name].items():
            u = choose(old[old_name][a,r],sx,sr)
            for key in FIELDS:
                if key=="gradient_frobenius":
                    diff = math.hypot(*(u[k]-v[k] for k in GRADIENT))
                elif key in u:
                    diff = abs(u[key]-v[key])
                else:
                    continue
                noise[key] = max(noise[key],diff)
    return noise


def patches(reports,H,a_interval,h,roffset,domain):
    found = {}
    for label,report in reports.items():
        groups = report["regions"]
        for prefix in ("inner","outer"):
            eligible = [v for k,v in groups.items() if ("outer" in k)==(prefix=="outer")]
            smallest = min(eligible,key=lambda g:g["principal_over_mu"]-g["max_directional_difference"])
            w = smallest["witness"]
            a,r = w["a"],w["rho"]
            if a_interval[0]<=a<=a_interval[1] and abs(r-h)<=roffset+1e-15:
                continue
            key = next((q for q in found if coincides(q[0],a,H) and coincides(q[1],r,H)), (a,r))
            found.setdefault(key,[]).append(dict(field=label,region=prefix,
                                                 reported_lower=smallest["principal_over_mu"]-
                                                 smallest["max_directional_difference"]))
    min_a,max_a,top = domain
    return [dict(centre=list(q),half_a=.05,half_r=.01,witnesses=found[q],
                 coords=sorted({(float(a),float(r)) for a in np.clip(np.linspace(q[0]-.05,q[0]+.05,9),min_a,max_a)
                                                     for r in np.clip(np.linspace(q[1]-.01,q[1]+.01,9),0.,top)}))
            for q in found]


def coverage_coordinates(audits,envelope):
    section,fixed,radial = section_points(audits,True)
    ref = next(iter(audits.values())).ref
    lo,hi = envelope
    edges = sorted({float(a) for audit in audits.values() for bounds in audit.columns
                    for a in bounds if lo<=a<=hi})
    coordinates = set(map(tuple,fixed))
    for a in edges:
        if not any(coincides(a,v,ref.H) for v in section.values()):
            section[f"axial face {a:.17g}"] = a
        coordinates.update((a+off,ref.h+r) for off in FIXED_AXIAL if lo<=a+off<=hi for r in radial)
    return np.array(sorted(coordinates)),section


def save(path,payload):
    with path.open("x") as stream:
        json.dump(payload,stream,indent=2,allow_nan=False)
        stream.write("\n")


def transfer_fields(audit_file,root):
    receipt=identity(audit_file)
    root_identity=identity(root)
    report=json.loads(audit_file.read_text(),parse_int=lambda token:-0.0 if token=="-0" else int(token))
    if (report["schema"]!="WJH27NativeOriginalTransferAuditV1" or report["physical_history"] is not False
            or report["eta"]!=0. or report["root_packet_sha256"]!=root_identity["sha256"]):
        raise ValueError("nested transfer analysis requires the bound native original-equation probe")
    data,meta=load_c1_state(root,True)
    if (report["reference_identity"]!=data["reference_identity"] or
            report["stages"]["solved"]["state_sha256"]!=data["state_sha256"]):
        raise ValueError("the native transfer probe used a different root")
    stages={}
    for name,stage in report["stages"].items():
        digest=stage["state_sha256"]
        if state_hash(coefficients(stage["coefficients"]))!=digest:
            raise ValueError(f"coefficient identity differs at {name}")
        stages[name]=(
            dict(data,coefficients=stage["coefficients"],state_sha256=digest,source_identity=report["source_identity"]),
            dict(packet=root_identity,coefficient_sha256=digest,source_identity=report["source_identity"],
                 native_original_residuals=meta["original_residuals"] if name=="solved" else None,
                 audit=receipt,stage=name,weak_continuity=stage["weak_continuity"]))
    return stages, report.get("root_sampled_residual_rows")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frozen",type=Path,required=True)
    p.add_argument("--table",type=Path,required=True)
    p.add_argument("--field",nargs=2,action="append",required=True,metavar=("LABEL","PACKET"))
    p.add_argument("--transfer-audit",type=Path,nargs=2,metavar=("NATIVE_AUDIT","ROOT_INPUT"))
    p.add_argument("--compare",nargs=2,action="append",required=True,metavar=("A","B"))
    p.add_argument("--margins",nargs=2,action="append",default=[],metavar=("LABEL","REPORT"))
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--seconds",type=float,default=550.)
    p.add_argument("--coverage",action="store_true")
    p.add_argument("--sampled-action",action="store_true")
    p.add_argument("--archived-noise",type=Path)
    p.add_argument("--noise-receipt",type=Path)
    p.add_argument("--noise-pair",nargs=2,action="append",default=[],metavar=("A","B"))
    args = p.parse_args()
    deadline=time.monotonic()+args.seconds
    core=json.loads(args.frozen.read_text())
    if len(core["chart_coordinates"])!=2175 or [tuple(v) for v in core["approaches"]]!=list(APPROACHES):
        raise ValueError("fixed comparison does not match the declared core")
    if args.archived_noise and identity(args.archived_noise)["sha256"]!=core["profile"]["sha256"]:
        raise ValueError("archival profile changed after data freeze")
    ref=po.load_reference(args.table)
    audits,states={},{}
    common_geometry=None
    for label,packet_string in args.field:
        if label in audits:
            raise ValueError("duplicate study field")
        packet=Path(packet_string)
        data,meta=load_c1_state(packet,True)
        geom=data["geometry"]
        if (meta["regularization"]["eta"]!=0. or geom["hill_height"]!=.0024 or
                not coincides(ref.h,geom["h"],ref.H)):
            raise ValueError("paired profiles require the pinned original height, reference and support")
        geometry={k:geom[k] for k in ("H","h","d0","left","right","top","hill_height","hill_half_length","normal_map_signature")}
        if common_geometry is not None and geometry!=common_geometry:
            raise ValueError("the physical chart differs across the compared fields")
        common_geometry=geometry
        audits[label]=Audit(data,ref,28,deadline)
        states[label]=dict(packet=identity(packet),coefficient_sha256=data["state_sha256"],
                           source_identity=data["source_identity"],native_original_residuals=meta["original_residuals"])
    domain=tuple(common_geometry[k] for k in ("left","right","top"))
    native_sampled=None
    if args.transfer_audit:
        stages,native_sampled=transfer_fields(*args.transfer_audit)
        for name,(data,state) in stages.items():
            if name in audits or any(data["geometry"][k]!=v for k,v in common_geometry.items()):
                raise ValueError("transfer stage label or physical chart differs across fields")
            audits[name]=Audit(data,ref,28,deadline)
            states[name]=state
    reports={}
    report_ids={}
    for label,file in args.margins:
        path=Path(file)
        report=json.loads(path.read_text())
        if report["native"]["sha256"]!=states[label]["packet"]["sha256"]:
            raise ValueError("independent margins refer to different coefficients")
        reports[label]=report
        report_ids[label]=identity(path)
    study_patches=patches(reports,ref.H,core["critical_interval"],ref.h,core["radial_absolute_offset"],domain)
    datasets=[("core",np.array(core["chart_coordinates"]),{})]
    datasets.extend((f"limiting-patch-{i+1}",np.array(q["coords"]),{}) for i,q in enumerate(study_patches))
    if args.coverage:
        a,section=coverage_coordinates(audits,core["envelope"])
        datasets.append(("coverage",a,section))
    noise={}
    if args.noise_receipt:
        previous=json.loads(args.noise_receipt.read_text())
        if previous["frozen"]["sha256"]!=identity(args.frozen)["sha256"]:
            raise ValueError("replay noise was measured on a different core")
        noise=previous["noise"]
    args.out.mkdir()
    results={}
    for label,points,section in datasets:
        csv_path=args.out/(label+".csv")
        views,extents={},{}
        with csv_path.open("x",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=PROFILE_COLUMNS)
            writer.writeheader()
            for name,audit in audits.items():
                rec,rel=write_profiles(writer,audit,name,points,section,include_outer=True)
                views[name]=approached(rec,ref.H,domain)
                extents[name]=extent(rec,rel)
        if label=="core" and args.archived_noise:
            noise=historical_noise(args.archived_noise,views)
        if label=="core":
            for a,b in args.noise_pair:
                if states[a]["coefficient_sha256"]!=states[b]["coefficient_sha256"]:
                    raise ValueError("noise-pair inputs are not identical coefficient fields")
                pairs=compare(views[a],views[b],core["scales"],{})["maxima"]
                for field in FIELDS:
                    noise[field]=max(noise.get(field,0.),pairs[field]["difference"])
        results[label]=dict(profile=identity(csv_path),fixed_coordinates=len(points),
                            extents=extents,comparisons={
            f"{a}/{b}":compare(views[a],views[b],core["scales"],noise) for a,b in args.compare})
    save(args.out/"paired-comparison.json",dict(
        schema="WJH27PairedAxialComparisonV1",frozen=identity(args.frozen),table=identity(args.table),
        analysis=identity(__file__),inputs=states,scales=core["scales"],noise=noise,margin_reports=report_ids,
        reported_margin_witnesses={name:{k:{n:v[n] for n in ("principal_over_mu","max_directional_difference","witness","samples")}
             for k,v in report["regions"].items()} for name,report in reports.items()},
        patches=study_patches,
        sampled_action_orders={k:sampled_action_order(a) for k,a in audits.items()} if args.sampled_action else {},
        transfer_root_signed_sampled_action=sampled_action_order(audits["solved"],native_sampled)
            if native_sampled is not None else None,
        gradient_layout="Cartesian 3x3 with nonzero entries allowed only in the x-z block",
        noise_pairs=args.noise_pair,
        noise_definition="Empirical archive-to-replay discrepancy for R4 columns and current independent repeated evaluation with identical coefficients.",
        datasets=results))


if __name__=="__main__":
    main()
