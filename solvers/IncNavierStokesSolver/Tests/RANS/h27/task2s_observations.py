#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replay H27 spatial observations and record original owner extrema."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import time

import numpy as np

import physical_oracle as oracle
from resources import atomic, identity


COMPONENTS = dict(U=("U",), W=("W",), k=("k",), R=("Rxx","Rxy","Rxr","Ryy","Ryr","Rrr"),
                  Pi=("Pi",), traction=("Vxr",), G=("Ux","Ur","Wx","Wr"), epsilon=("epsilon",))


def load_state(path, ref, order, deadline):
    data, meta = oracle.load_native(path)
    oracle.require(meta["steady_mode"] and data["time"] == data["step"] == 0,
                   "expected a registered stationary zero-time state")
    return oracle.PhysicalOracle(data, ref, order, deadline)


def polynomial_basis(element, field, x, r):
    """Evaluate Modified_A in [0,1] and differentiate in physical code units."""
    def modes(y, delta):
        p = y*(1-y)
        dp = 1-2*y
        values = [1-y, y, p]
        slopes = [-np.ones_like(y), np.ones_like(y), dp]
        if field != "p":
            values.append((4*y-2)*p)
            slopes.append((4*y-2)*dp+4*p)
        return np.column_stack(values), np.column_stack(slopes)/delta

    xl, xr, rl, rr = element["bounds"]
    bx, dx = modes((x-xl)/(xr-xl), xr-xl)
    br, dr = modes((r-rl)/(rr-rl), rr-rl)
    signs = np.array(element[oracle.SIGNS[field]])
    return [(left[:, :, None]*right[:, None, :]).reshape(len(x), len(signs))*signs
            for left,right in ((br,bx),(br,dx),(dr,bx))]


def orientation_audit(state):
    results = {}
    for field, key in oracle.SIGNS.items():
        usable = [e for e in state.elements if e[key]]
        negative = [e for e in usable if -1. in e[key]]
        # Cover both exterior and interior and every owner with a reversed map.
        stride = max(1,len(usable)//24)
        witnesses = usable[::stride]+negative
        maximum = 0.
        for e in witnesses:
            xl, xr, rl, rr = e["bounds"]
            x = xl+(xr-xl)*np.array([.17,.29,.63,.88])
            r = rl+(rr-rl)*np.array([.25,.68,.19,.74])
            explicit = polynomial_basis(e,field,x,r)
            observed = oracle.basis(e,field,x,r)
            maximum = max(maximum, *(float(np.max(abs(a-b)))/(1+np.max(abs(a)))
                                     for a,b in zip(explicit,observed)))
        test = dict(usable[len(usable)//2])
        sign = list(test[key])
        width = 3 if field == "p" else 4
        # Flip an interior radial mode. The explicit and production paths must
        # both change the derivative column by the signed map.
        column = width*(width-1)
        sign[column] *= -1
        test[key] = sign
        xl,xr,rl,rr = test["bounds"]
        x = np.array([xl+.31*(xr-xl)])
        r = np.array([rl+.24*(rr-rl)])
        rest = polynomial_basis(test,field,x,r)
        got = oracle.basis(test,field,x,r)
        flipped = float(np.max(abs(rest[2][0,column]-got[2][0,column]))/(1+abs(rest[2][0,column])))
        original = oracle.basis(usable[len(usable)//2],field,x,r)
        oracle.require(abs(rest[2][0,column]+original[2][0,column])
                       <= 2e-12*(1+abs(rest[2][0,column])), "reversed radial mode did not reverse")
        results[field] = dict(elements=len(usable), elements_with_negative_sign=len(negative),
                              negative_entries=sum(e[key].count(-1.) for e in negative),
                              relative_basis_derivative_disagreement=maximum,
                              flipped_radial_relative_disagreement=flipped)
    oracle.require(max(max(v["relative_basis_derivative_disagreement"],v["flipped_radial_relative_disagreement"])
                       for v in results.values()) <= 5e-12, "independent signed basis disagreement")
    return results


def cell_index(elements):
    owners = {}
    xa, ra = set(), set()
    x_keys = oracle.coordinate_aliases(e["bounds"][i] for e in elements for i in (0,1))
    r_keys = oracle.coordinate_aliases(e["bounds"][i] for e in elements for i in (2,3))
    for e in elements:
        a,b,c,d = e["bounds"]
        xi,ri = (x_keys[a],x_keys[b]),(r_keys[c],r_keys[d])
        xa.add(xi)
        ra.add(ri)
        oracle.require((xi,ri) not in owners,"duplicate rectangular owner")
        owners[xi,ri] = e
    oracle.require(len(owners) == len(xa)*len(ra),"expected the registered tensor grid")
    return sorted(xa), sorted(ra), owners


def overlapping(a, b):
    return [(left,right) for left in a for right in b if min(left[1],right[1])-max(left[0],right[0]) > 1e-12]


def common_cells(first, second):
    ax, ar, aa = cell_index(first.elements)
    bx, br, bb = cell_index(second.elements)
    radial = overlapping(ar,br)
    for ia,ib in overlapping(ax,bx):
        for ja,jb in radial:
            a, b = aa[ia,ja], bb[ib,jb]
            al, ah, ab, at = a["bounds"]
            bl, bh, blr, bhr = b["bounds"]
            yield a,b,(max(al,bl),min(ah,bh),max(ab,blr),min(at,bhr))


class Rows:
    def __init__(self):
        self.records = {}

    def add(self, groups, left, right, x, r, weights, a, b, ref):
        av, bv = oracle.PhysicalOracle.observable(a), oracle.PhysicalOracle.observable(b)
        for field in oracle.NUMERICAL:
            if field == "epsilon" and not a["outer"]:
                continue
            scale = (ref.Kin if field == "k" else ref.tau if field in ("R","Pi","traction")
                     else ref.A/r if field == "epsilon" else 1.)
            delta = ((av[field]-bv[field])/scale).reshape(len(x),-1)
            pos = np.argmax(np.abs(delta),axis=0)
            squares = np.einsum("i,ij,ij->j",weights,delta,delta)
            maxima = abs(delta[pos,np.arange(len(pos))])
            for group in groups:
                row = self.records.setdefault(group+"/"+field, dict(measure=0.,
                    squares=np.zeros(len(pos)), maxima=np.full(len(pos),-1.), locations=[None]*len(pos)))
                row["measure"] += float(np.sum(weights))
                row["squares"] += squares
                for i in np.flatnonzero(maxima > row["maxima"]):
                    j = pos[i]
                    row["maxima"][i] = maxima[i]
                    row["locations"][i] = dict(component=COMPONENTS[field][i], x=float(x[j]), r=float(r[j]),
                        coarse_owner=left["geometry"], fine_owner=right["geometry"],
                        owner_bounds=[left["bounds"],right["bounds"]], signed_difference=float(delta[j,i]),
                        included_with_quadrature_weight=bool(weights[j]>0))

    def finish(self):
        records = {}
        for name, row in sorted(self.records.items()):
            field = name.rsplit("/",1)[1]
            limit = oracle.NUMERICAL[field]
            rms = np.sqrt(row["squares"]/row["measure"])
            maximum = row["maxima"]
            status = "MET" if max(np.max(rms),np.max(maximum)) <= limit else "NOT_MET"
            records[name] = dict(RMS=float(np.max(rms)), MAX=float(np.max(maximum)),
                component_RMS=rms.tolist(), component_MAX=maximum.tolist(), limits=[limit,limit],
                status=status, measure=row["measure"])
            if status == "NOT_MET":
                worst = int(np.argmax(maximum))
                integral = int(np.argmax(rms))
                records[name].update(peak=row["locations"][worst],
                                     RMS_controlling_component=COMPONENTS[field][integral],
                                     RMS_component_peak=row["locations"][integral])
        return records


def cuts(first, second, e, f, x, lower, upper):
    radii = [model.ref.cutoff(model.trace(owner,"k",np.array([x])))[0]
             for model,owner in ((first,e),(second,f))]
    return first.cuts(lower,upper,radii)


def station_points(a, w, breaks, multiplicity):
    mids = (np.asarray(breaks[:-1])+breaks[1:])/2
    half = np.diff(breaks)/2
    rad = (mids[:,None]+half[:,None]*a).ravel()
    weights = (half[:,None]*w/multiplicity).ravel()
    return np.r_[rad,breaks], np.r_[weights,np.zeros(len(breaks))]


def compare(first,second,quiet=()):
    rows, count = Rows(), 0
    local = Rows()
    edges = sorted({e["bounds"][i]
                    for model in (first,second) for e in model.elements for i in (0,1)})
    multiplicity = {x:2 if -10. < x < 15. and any(abs(x-edge)<1e-12 for edge in edges) else 1
                    for x in oracle.STATIONS}
    for e,f,bounds in common_cells(first,second):
        first.clock()
        count += 1
        left,right,bottom,top = bounds
        region = "outer" if e["epsilon"] else "inner"
        oracle.require(bool(e["epsilon"]) == bool(f["epsilon"]), "common cell crossed H")
        x,r,weights = first.rule(e,bounds,lambda x:second.ref.cutoff(second.trace(f,"k",x)),
                                 second.cutoff_x_splits(f,bounds))
        corners = [cuts(first,second,e,f,xe,bottom,top) for xe in (left,right)]
        extra_radius = np.concatenate(corners)
        x = np.r_[x,np.repeat([left,right],[len(v) for v in corners])]
        r = np.r_[r,extra_radius]
        weights = np.r_[weights,np.zeros(len(extra_radius))]
        a,b = first.sample(e,x,r),second.sample(f,x,r)
        rows.add(("whole",region),e,f,x,r,weights,a,b,first.ref)
        for lo,hi in quiet:
            if lo-1e-12 <= left < right <= hi+1e-12:
                local.add((f"{lo:g}..{hi:g}/whole",f"{lo:g}..{hi:g}/{region}"),
                          e,f,x,r,weights,a,b,first.ref)
        for x in oracle.STATIONS:
            if left-1e-12 <= x <= right+1e-12:
                r,w = station_points(*first.gauss,cuts(first,second,e,f,x,bottom,top),multiplicity[x])
                xx = np.full(len(r),x)
                rows.add((f"x={x:g}/whole",f"x={x:g}/{region}"),e,f,xx,r,w,
                         first.sample(e,xx,r),second.sample(f,xx,r),first.ref)
    result = rows.finish()
    expected = {"whole":first.ref.top,"inner":first.ref.H,"outer":first.ref.top-first.ref.H}
    weight_error = 0.
    for key,row in result.items():
        group, field = key.split("/")[-2:]
        if field == "epsilon":
            group = "outer"
        measure = expected[group] * (1 if key.startswith("x=") else 25.)
        weight_error = max(weight_error,abs(row["measure"]-measure)/measure)
    oracle.require(len(result) == 690 and weight_error < 1e-10, "original row inventory or station measure changed")
    local_rows = local.finish()
    oracle.require(len(local_rows) == len(quiet)*23, "redistributed x-strip observation missing")
    return result,local_rows,dict(common_owner_subcells=count, relative_measure_error=weight_error,
                                 station_multiplicity={str(x):m for x,m in multiplicity.items()})


def check_archived(rows, archived):
    old = archived["norms"]
    oracle.require(old.keys() == rows.keys(), "archived observations differ")
    metrics = ("RMS","MAX","measure","component_RMS","component_MAX")
    difference = max(float(np.max(abs(np.asarray(value[key])-np.asarray(old[name][key]))))
                     for name,value in rows.items() for key in metrics)
    mismatched = [name for name,value in rows.items() if value["status"] != old[name]["status"]]
    return dict(status="PASS" if difference <= 5e-10 and not mismatched else "NOT_MET",
                maximum_absolute_replay_difference=difference, mismatched_rows=mismatched,
                old_failed_rows=[name for name,value in old.items() if value["status"]!="MET"])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("table","first","second","out"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--archived",type=Path)
    p.add_argument("--quiet-plan",type=Path,help="Add every reallocated upstream strip as a separate observation.")
    p.add_argument("--seconds",type=float,required=True)
    args=p.parse_args()
    end=time.monotonic()+args.seconds
    ref=oracle.load_reference(args.table)
    a=load_state(args.first,ref,20,end)
    b=load_state(args.second,ref,20,end)
    oracle.require(a.data["source_identity"] == b.data["source_identity"], "comparison crosses native sources")
    signed_basis=dict(first=orientation_audit(a),second=orientation_audit(b))
    quiet = [(row["x"][0],row["x"][2]) for row in json.loads(args.quiet_plan.read_text())["removed"]] if args.quiet_plan else []
    rows,extra,sampling=compare(a,b,quiet)
    failures=[key for key,value in rows.items() if value["status"]!="MET"]
    info=dict(schema="WJH27Task2RSOwnerComparisonV1", first=identity(args.first),
              second=identity(args.second), table=identity(args.table), source=a.data["source_identity"],
              station_inventory=oracle.STATIONS, signed_basis=signed_basis,
              rows=rows, sampling=sampling, failures=failures, comparison_status="NOT_MET" if failures else "MET")
    if args.quiet_plan:
        info["redistributed_strip_controls"]=dict(plan=identity(args.quiet_plan),rows=extra,
            failures=[name for name,item in extra.items() if item["status"]!="MET"])
    if args.archived:
        info["archive_replay"]=dict(artifact=identity(args.archived),
                                   **check_archived(rows,json.loads(args.archived.read_text())))
    atomic(args.out,info)
    print(len(rows),"rows,",len(failures),"misses",info.get("archive_replay",{}).get("status","NEW_PAIR"),flush=True)


if __name__=="__main__":
    main()
