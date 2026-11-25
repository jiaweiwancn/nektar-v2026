#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze R4 physical samples and choose the common nested axial envelope."""

import argparse
import csv
import json
import math
from pathlib import Path

from c1_support_witness import choose
from resources import identity
from task2_outlet_baseline import mesh_axes


CORE_INTERVAL = (.875, 2.5)
ABSOLUTE_RADIAL_OFFSET = .0015
APPROACHES = ((-1,-1),(-1,1),(1,-1),(1,1))
SCALARS = ("margin","Gxz","Gxx","Gzx","Gzz","divergence","pressure","k_load","convection")
PARENT_NAMES = ("Mref","Mshift-ref")
ARCHIVED = (*PARENT_NAMES,"Mxr","Mshift")


def core_rows(profile):
    views = {name:{} for name in ARCHIVED}
    with profile.open(newline="") as stream:
        for row in csv.DictReader(stream):
            name = row["mesh"]
            if name not in views or row["fixed"]!="True":
                continue
            a,r = float(row["a"]),float(row["r"])
            if not (CORE_INTERVAL[0]<=a<=CORE_INTERVAL[1]
                    and abs(float(row["r_minus_h"]))<=ABSOLUTE_RADIAL_OFFSET+1e-15):
                continue
            values = {field:float(row[field]) for field in SCALARS}
            values["G"] = math.hypot(*(values[k] for k in ("Gxx","Gxz","Gzx","Gzz")))
            values["total_load"] = sum(values[k] for k in ("pressure","k_load","convection"))
            values.update(axial_side=int(row["axial_side"]),radial_side=int(row["radial_side"]))
            views[name].setdefault((a,r),[]).append(values)
    return views


def buffered_interval(axis):
    identified = [i for i,(a,b) in enumerate(zip(axis,axis[1:]))
                  if max(a,CORE_INTERVAL[0]) < min(b,CORE_INTERVAL[1])]
    if not identified or identified[0]==0 or identified[-1]+2>=len(axis):
        raise ValueError("critical strip needs one existing axial neighbor on each side")
    return axis[identified[0]-1],axis[identified[-1]+2]


def subdivisions(axis, envelope):
    if not axis[0]<=envelope[0]<envelope[1]<=axis[-1]:
        raise ValueError("refinement envelope outside the chart")
    inserted = [v for v in envelope if v not in axis]
    edges = sorted([*axis,*inserted])
    midpoints = [a+(b-a)/2 for a,b in zip(edges,edges[1:]) if envelope[0]<=a<b<=envelope[1]]
    if not midpoints or any(v in edges for v in midpoints):
        raise ValueError("axial subdivisions are empty or rounded onto existing faces")
    return inserted,midpoints


def freeze(profile, parent_paths, out):
    views = core_rows(profile)
    points = sorted(views[PARENT_NAMES[0]])
    if len(points)!=2175 or any(rows.keys()!=views[PARENT_NAMES[0]].keys() for rows in views.values()):
        raise ValueError("did not reconstruct the archived 2175 common physical coordinates")

    def samples(name):
        return [choose(views[name][key],sx,sr) for key in points for sx,sr in APPROACHES]

    values = {name:samples(name) for name in ARCHIVED}
    comparisons = {}
    for a,b in (PARENT_NAMES,("Mxr","Mref"),("Mshift","Mshift-ref")):
        comparisons[f"{a}/{b}"] = {
            field:max(abs(x[field]-y[field]) for x,y in zip(values[a],values[b]))
            for field in (*SCALARS,"total_load")}
        comparisons[f"{a}/{b}"]["gradient_frobenius"] = max(
            math.hypot(*(x[k]-y[k] for k in ("Gxx","Gxz","Gzx","Gzz")))
            for x,y in zip(values[a],values[b]))
    comparison = comparisons["Mref/Mshift-ref"]
    expected = dict(margin=.03868827439139451,Gxz=.13669091414769152,divergence=.06626143662171902)
    if any(abs(comparison[k]-v)>5e-15 for k,v in expected.items()):
        raise ValueError("archived comparison differs on the frozen coordinates")
    scales = {key:max(abs(row[key]) for name in PARENT_NAMES for row in values[name])
              for key in ("G","total_load","pressure","k_load","convection")}
    parent_buffers = [buffered_interval(mesh_axes(path)[0]) for path in parent_paths]
    envelope = [min(s[0] for s in parent_buffers),max(s[1] for s in parent_buffers)]
    with out.open("x") as stream:
        json.dump(dict(schema="WJH27FrozenPairedAxialV1",profile=identity(profile),
            chart_coordinates=points,approaches=APPROACHES,
            critical_interval=CORE_INTERVAL,radial_absolute_offset=ABSOLUTE_RADIAL_OFFSET,
            parent_names=PARENT_NAMES,parent_sessions=[identity(p) for p in parent_paths],
            parent_buffers=parent_buffers,envelope=envelope,scales=scales,
            archived_comparisons=comparisons),stream,allow_nan=False,separators=(",",":"))
        stream.write("\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("profile",type=Path)
    p.add_argument("--parents",type=Path,nargs=2,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    freeze(a.profile,a.parents,a.out)
