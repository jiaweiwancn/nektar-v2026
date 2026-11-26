#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Localize solved field changes and split the signed Q3 derivative contributions."""

import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from c1_oracle import Audit, load_c1_state
from c1_paired_profiles import approached, compare, save
from c1_support_controls import contains
from c1_support_witness import PROFILE_COLUMNS, extent, write_profiles
import physical_oracle as po
from resources import identity


def directional_modes(audit, cell, a, r):
    chart = (np.array([a]), np.array([r]))
    metric = audit.hill.metrics(cell, *chart)
    _, da, dr = po.basis(cell, "u", *chart)
    ids = np.asarray(cell["u"])
    velocity = np.stack([audit.s["flow"][ids], audit.s["flow"][audit.n+ids]], axis=1)
    g_a = da[0,:,None,None]*velocity[:,:,None]*metric["da"][0]
    g_r = dr[0,:,None,None]*velocity[:,:,None]*metric["dr"][0]
    local = g_a + g_r
    _, _, gradient, *_ = audit.sample(cell, *chart, metric)
    closure = float(np.abs(np.sum(local,axis=0)-gradient[0]).max())
    if closure > 2e-12:
        raise ValueError("independent signed velocity modes did not reconstruct Cartesian gradient")
    return dict(owner=cell["geometry"], bounds=cell["bounds"], closure=closure,
                cartesian=metric["position"][0].tolist(),
                local_velocity_coefficients=(velocity*np.asarray(cell["us"])[:,None]).tolist(),
                axial_metric=metric["da"][0].tolist(),radial_metric=metric["dr"][0].tolist(),
                xi_eta=[2*(v-lo)/(hi-lo)-1 for v,(lo,hi) in zip((a,r),
                        (cell["bounds"][:2],cell["bounds"][2:]))],
                a_gradient=g_a.sum(axis=0).tolist(), r_gradient=g_r.sum(axis=0).tolist(),
                row_groups=local.reshape(4,4,2,2).sum(axis=1).tolist(),
                axial_groups=local.reshape(4,4,2,2).sum(axis=0).tolist(),
                radial_derivative_modes=g_r.reshape(4,4,2,2).sum(axis=1).tolist())


def sections(h, H, foot):
    return dict(foot=(foot,h+.0015,.05,.015),
                response=(foot-.005,h,.025,.007),
                upstream=(-6.25,H,.09,.012),
                off_core=(.17578125,h-.015,.04,.008))


def patch_geometry(station, audits, height):
    a0,r0,da,dr = station
    ax = set(a0 + np.linspace(-da,da,11))
    ax.add(a0)
    rs = set(r0 + np.linspace(-dr,dr,41))
    rs.add(r0)
    if abs(height-r0)<=dr:
        rs.add(height)
    relative = set()
    for audit in audits.values():
        for cell in audit.mesh:
            if not contains(cell,a0,r0,audit.ref.H):
                continue
            for i,(value,lo,hi) in enumerate(((a0,*cell["bounds"][:2]),(r0,*cell["bounds"][2:]))):
                for fraction in (0.,.001,.01,.05,.2):
                    for edge in (lo,hi):
                        sample = edge + fraction*(value-edge)
                        limit = da if i==0 else dr
                        if abs(sample-value) <= limit:
                            relative.add((sample,r0) if i==0 else (a0,sample))
    base={(a,r) for a in ax for r in rs}
    return np.array(sorted(base|relative)), sorted(relative)


def plot(output, cases, cuts, lines):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig,axes=plt.subplots(3,4,figsize=(13,8),constrained_layout=True)
    fields=(("margin","Principal response / μ"),("divergence","Strong divergence"),
            ("Gxz","Cartesian Gxz"))
    for col,(label,(_,r0,_,_)) in enumerate(cuts.items()):
        for e,(key,title) in enumerate(fields):
            axis=axes[e,col]
            for case,color in zip(cases,("C0","C1","C2","C3","C4","C5")):
                for i,rows in enumerate(lines[label][case].values()):
                    rows.sort(key=lambda row:row["r"])
                    axis.plot([1e3*(row["r"]-r0)*.04 for row in rows],
                              [row[key] for row in rows],color=color,alpha=.85,
                              label=case if i==0 else None)
            axis.axvline(0,color=".70",lw=.6)
            axis.grid(alpha=.2)
            if e==0:
                axis.set_title(label.replace("_"," "))
                axis.legend(fontsize=7,frameon=False,loc="best")
            if col==0:
                axis.set_ylabel(title)
            if e==2:
                axis.set_xlabel("Radial offset (mm)")
    fig.savefig(output/"local-variation.png",dpi=155)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table",type=Path,required=True)
    p.add_argument("--frozen",type=Path,required=True)
    p.add_argument("--field",action="append",nargs=2,required=True)
    p.add_argument("--compare",action="append",nargs=2,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--seconds",type=float,default=150.)
    args=p.parse_args()
    ref=po.load_reference(args.table)
    core=json.loads(args.frozen.read_text())
    audits,inputs={},{}
    deadline=time.monotonic()+args.seconds
    for name,path in args.field:
        source=Path(path)
        data,meta=load_c1_state(source,True)
        if (name in audits or data["geometry"]["hill_height"]!=.0024 or meta["regularization"]["eta"]!=0.
                or max(meta["original_residuals"])>1e-10 or meta.get("accepted_checkpoint")):
            raise ValueError("localization requires distinct solved original fixed-height fields")
        audits[name]=Audit(data,ref,28,deadline)
        inputs[name]=dict(packet=identity(source),coefficient=data["state_sha256"])
    first=next(iter(audits.values()))
    domain=(first.hill.left,first.hill.right,first.top)
    cuts=sections(ref.h,ref.H,first.hill.L)
    results,lines={},{}
    args.out.mkdir()
    for label,station in cuts.items():
        points,relative=patch_geometry(station,audits,ref.h)
        a0,r0,*_=station
        path=args.out/(label+".csv")
        rows,extents,lines[label]={},{},{}
        with path.open("x",newline="") as stream:
            writer=csv.DictWriter(stream,fieldnames=PROFILE_COLUMNS)
            writer.writeheader()
            for name,audit in audits.items():
                physical,_=write_profiles(writer,audit,name,points,{},include_outer=True)
                rows[name]=approached(physical,ref.H,domain)
                extents[name]=extent(physical,[])
                lines[label][name]={}
                for (a,r),owners in physical.items():
                    if a==a0:
                        for row in owners:
                            if row["axial_side"] in (0,-1):
                                lines[label][name].setdefault(row["owner"],[]).append(row)
        diffs,modal={},{}
        for left,right in args.compare:
            key=f"{left}/{right}"
            diffs[key]=compare(rows[left],rows[right],core["scales"],{})
            modal[key]={}
            for field in ("margin","gradient_frobenius","divergence","total_load"):
                w=diffs[key]["maxima"][field]
                a,r=w["chart"]
                owners=[]
                for case,side in ((left,"a"),(right,"b")):
                    cell,=(e for e in audits[case].mesh if e["geometry"]==w[side]["owner"])
                    owners.append(directional_modes(audits[case],cell,a,r))
                delta={axis:(np.array(owners[1][axis])-np.array(owners[0][axis])).tolist()
                       for axis in ("a_gradient","r_gradient","row_groups","axial_groups","radial_derivative_modes")}
                modal[key][field]=dict(owners=owners,right_minus_left=delta)
        results[label]=dict(station=station,fixed_and_relative_points=len(points),
                            cell_relative_points=relative,profile=identity(path),
                            extents=extents,comparisons=diffs,mode_changes=modal)
    plot(args.out, list(audits), cuts, lines)
    save(args.out/"local-variation.json",dict(schema="WJH27LocalVariationV1",table=identity(args.table),
         frozen=identity(args.frozen),inputs=inputs,stations=results,
         figure=identity(args.out/"local-variation.png"),
         mode_order="Modified_A, axial index fastest; signed modes in Cartesian physical derivatives"))


if __name__=="__main__":
    main()
