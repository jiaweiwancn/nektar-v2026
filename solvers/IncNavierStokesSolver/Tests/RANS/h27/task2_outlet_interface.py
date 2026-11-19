#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Measure one-sided H flux convergence on the frozen common streamwise grid."""

import argparse
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss

import physical_oracle as po
from resources import atomic
from task2_outlet_baseline import mesh_axes, verified
from task2_outlet_case import get_context
from task2_outlet_qual_followup import adjacent_H
from task2_outlet_qualification import SHIFT, load_case, pass_audit, split_band
from task2r_run import require

SIDES = ("lower", "upper")
FIELDS = ("k_diffusion_over_A", "x_traction_over_tau0", "r_traction_over_tau0")
LIMIT = 1e-3
RESOLUTION_FLOOR = 1e-5
CONTRACTION = .8


def owners(state):
    faces = {side:[] for side in SIDES}
    for e in state.elements:
        bounds = e["bounds"]
        for side,edge in (("lower",3),("upper",2)):
            if abs(bounds[edge]-state.ref.H)<1e-12:
                faces[side].append(e)
    for side in SIDES:
        faces[side].sort(key=lambda e:e["bounds"][0])
    require(faces["lower"] and len(faces["lower"])==len(faces["upper"]),
            "H does not have one physical owner per side and x interval")
    return faces


def traces(state, face, x):
    data = state.sample(face,x,np.full(len(x),state.ref.H))
    ref = state.ref
    return np.vstack((-data["Dk"]*data["gradk"][:,1]/ref.A,
        data["V"][:,0,2]/ref.tau, (data["V"][:,2,2]-data["Pi"])/ref.tau,
        data["W"],data["W"]*data["k"]/ref.Kin))


def measure(values, weights, x, physical_owners):
    i = int(np.argmax(np.abs(values)))
    return dict(RMS=float(np.sqrt(np.dot(weights,values**2)/weights.sum())),
        MAX=float(abs(values[i])),signed_MAX=float(values[i]),x=float(x[i]),
        owners=physical_owners[i],included_with_quadrature_weight=bool(weights[i]>0))


def freeze_check(axes, H, support):
    x,z = axes[0]
    require(all(x==grid[0] for grid in axes), "H comparison requires unchanged streamwise axes")
    widths = []
    for _,r in axes:
        h, = [i for i,v in enumerate(r) if abs(v-SHIFT-H)<1e-12]
        widths.append([r[h]-r[h-1],r[h+1]-r[h]])
    require(all(abs(b-a/2)<1e-12 for old,new in zip(widths,widths[1:]) for a,b in zip(old,new)),
            "H owners must halve on both sides in each witness")
    require(support==adjacent_H(z), "frozen H support differs from coarse owner intervals")
    expected = z
    for _,actual in axes[1:]:
        expected,_ = split_band(expected,*support)
        require(actual==expected, "H witness changed a second radial support")
    return widths


def compare(root, levels, support, out, seconds):
    context = get_context(root)
    ref = po.load_reference(context[2])
    deadline = time.monotonic()+seconds
    bindings, states, axes, audits = [], [], [], []
    for level in levels:
        audit, proof = pass_audit(root,level)
        state, manifest = load_case(root,level,context,ref,20,deadline)
        require(manifest["state"]==audit["provenance"]["state"], "H audit describes a different state")
        x,z,_ = mesh_axes(verified(manifest["mesh"]))
        bindings.append(manifest)
        states.append(state)
        axes.append((x,z))
        audits.append(proof)
    widths = freeze_check(axes,ref.H,support)
    x_edges = axes[0][0]
    faces = [owners(s) for s in states]
    t,w = leggauss(20)
    x_stack, w_stack, flux = [], [], []
    owner_stack = []
    for i,(lo,hi) in enumerate(zip(x_edges,x_edges[1:])):
        probes = sorted({lo,hi,*(x for x in po.STATIONS if lo<x<hi)})
        x = np.r_[(lo+hi)/2+(hi-lo)*t/2, probes]
        weight = np.r_[(hi-lo)*w/2,np.zeros(len(probes))]
        x_stack.append(x)
        w_stack.append(weight)
        records = [[face[side][i] for side in SIDES] for face in faces]
        require(all(max(abs(a-b) for a,b in zip(e["bounds"][:2],(lo,hi)))<1e-12
                    for row in records for e in row), "H face owner on wrong x segment")
        flux.append(np.asarray([[traces(s,e,x) for e in row] for s,row in zip(states,records)]))
        indices = [{side:dict(owner=e["geometry"],bounds=e["bounds"])
                    for side,e in zip(SIDES,row)} for row in records]
        owner_stack.extend([indices]*len(x))
    values = np.concatenate(flux,axis=-1)
    x,weight = np.concatenate(x_stack),np.concatenate(w_stack)
    side_results, shared, jumps = {}, {}, {}
    for i,side in enumerate(SIDES):
        side_results[side] = {}
        for j,field in enumerate(FIELDS):
            first,second,third = values[:,i,j,:]
            rows = {tag:measure(delta,weight,x,owner_stack) for tag,delta in (
                ("level_1_minus_0",second-first),("level_2_minus_1",third-second),
                ("level_2_minus_0",third-first))}
            ratios = {}
            for norm in ("RMS","MAX"):
                a,b,c = (rows[k][norm] for k in rows)
                ratios[norm] = dict(ratio=b/a if a else None,
                    contraction=b<=CONTRACTION*a or max(a,b)<=RESOLUTION_FLOOR,
                    direct_gate=c<=LIMIT)
            side_results[side][field] = dict(rows=rows,gates=ratios)
    for i,name in enumerate(levels):
        delta = values[i,1,:,:]-values[i,0,:,:]
        jumps[name] = {k:measure(delta[j],weight,x,owner_stack) for j,k in enumerate(FIELDS)}
        shared[name] = {k:measure(delta[j+3],weight,x,owner_stack)
                        for j,k in enumerate(("mass","k_advection_over_Kin"))}
    passed = all(v["direct_gate"] and v["contraction"] for items in side_results.values()
                 for field in items.values() for v in field["gates"].values())
    passed &= all(v["MAX"]<=LIMIT for v in jumps[levels[-1]].values())
    passed &= all(f["MAX"]<=1e-10 for row in shared.values() for f in row.values())
    atomic(out,dict(schema="WJH27OneSidedHOwnerFluxConvergenceV1",
        status="MET" if passed else "NOT_MET",levels=levels,support=support,H=ref.H,
        bindings=bindings,audits=audits,successive_widths=widths,
        common_normal=[0,1],paired_side_flux=side_results,shared=shared,one_sided_jumps=jumps,
        direct_RMS_MAX_limit=LIMIT,successive_ratio_limit=CONTRACTION,
        successive_resolved_floor=RESOLUTION_FLOOR,shared_MAX_limit=1e-10,
        sample_rule="20-point streamwise Gauss with both endpoints and retained stations",
        x_columns=len(x_edges)-1,weighted_measure=float(weight.sum())))
    print("H owner flux gates",levels,"MET" if passed else "NOT_MET",flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root",type=Path)
    parser.add_argument("out",type=Path)
    parser.add_argument("levels",nargs=3)
    parser.add_argument("--support",required=True,nargs=2,type=float)
    parser.add_argument("--seconds",required=True,type=float)
    args=parser.parse_args()
    compare(args.root,args.levels,args.support,args.out,args.seconds)
