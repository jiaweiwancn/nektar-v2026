#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Design a joint outlet refinement using verified quiet upstream x strips."""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from numpy.polynomial import Polynomial

from coupled_fixtures import mesh
import physical_oracle as po
from physical_reference_refinement import axes
from resources import atomic, identity
from task2s_observations import cell_index, load_state, polynomial_basis


def trace(cardinal, at, derivative=False):
    return np.column_stack([(p.deriv() if derivative else p)(at) for p in cardinal])


def cardinal(nodes):
    polynomials = [Polynomial.fromroots(np.delete(nodes,i)) for i in range(len(nodes))]
    return [p/p(node) for p,node in zip(polynomials,nodes)]


def coefficients(model, owner, x, r):
    b,dx,dr = polynomial_basis(owner,"u",x,r)
    ids = np.array(owner["u"])
    u,w = model.s["flow"][ids],model.s["flow"][model.n+ids]
    k = polynomial_basis(owner,"k",x,r)[0] @ model.s["k"][owner["k"]]
    p = polynomial_basis(owner,"p",x,r)[0] @ model.s["flow"][2*model.n+np.array(owner["p"])]
    outer = bool(owner["epsilon"])
    e = polynomial_basis(owner,"epsilon",x,r)[0] @ model.s["epsilon"][owner["epsilon"]] if outer else np.zeros(len(x))
    K = model.trace(owner,"k",x)
    return np.column_stack((b@u,b@w,k,e,p,K,dx@u,dr@u,dx@w,dr@w))


def stresses(ref, a, r, outer):
    k, E, K = a[:,2],a[:,3],a[:,5]
    G = np.zeros((len(r),3,3))
    G[:,(0,0,2,2),(0,2,0,2)] = a[:,6:]
    strain = G+G.transpose(0,2,1)
    R = k[:,None,None]*(2/3)*np.eye(3)
    if not outer:
        E = K**1.5/ref.cl/np.maximum(r,ref.cutoff(K))
        R[:,0,2] = R[:,2,0] = -ref.aux(K,r)*strain[:,0,2]
    chi = np.clip((r/ref.H-.4)/.35,0.,1.)
    chi = chi**3*(10-15*chi+6*chi*chi)
    active = chi>0
    if np.any(active):
        full = po.full_wj(k[active],E[active],G[active])[0]
        mix = chi[active,None,None]
        R[active] = (1-mix)*R[active]+mix*full
    return R,ref.nu*strain[:,0,2]-R[:,0,2]


def merge_errors(model, left, right, radial, owners):
    lo, seam, hi = left[0],left[1],right[1]
    width = hi-lo
    q3 = np.array([0.,(1-1/np.sqrt(5))/2,(1+1/np.sqrt(5))/2,1.])
    q2 = np.array([0.,.5,1.])
    nodes = np.linspace(0.,1.,9)
    x_two = np.r_[lo+(seam-lo)*nodes,seam+(hi-seam)*nodes]
    d = {key:0. for key in po.NUMERICAL}
    p3, p2 = cardinal(q3), cardinal(q2)
    az = np.array([0.,.166666666666667,.5,.833333333333333,1.])
    f3,f2 = trace(p3,(x_two-lo)/width),trace(p2,(x_two-lo)/width)
    grad = trace(p3,(x_two-lo)/width,derivative=True)/width
    for rspan in radial:
        model.clock()
        e,f = owners[left,rspan],owners[right,rspan]
        rb,rt = e["bounds"][2:]
        r = rb+(rt-rb)*az
        def values(points, left_only=None):
            chunks = []
            for i,x in enumerate(points):
                side = x <= seam if left_only is None else i < left_only
                chunks.append(coefficients(model,e if side else f,np.full(len(r),x),r))
            return np.asarray(chunks)
        old = values(x_two,left_only=len(nodes))
        sample3, sample2 = values(lo+width*q3), values(lo+width*q2)
        merged = np.zeros_like(old)
        for i in (0,1,2,3,5,7,9):
            merged[:,:,i] = f3@sample3[:,:,i]
        merged[:,:,4] = f2@sample2[:,:,4]
        merged[:,:,6] = grad@sample3[:,:,0]
        merged[:,:,8] = grad@sample3[:,:,1]
        dif = abs(old-merged)
        d["U"] = max(d["U"],float(dif[:,:,0].max()))
        d["W"] = max(d["W"],float(dif[:,:,1].max()))
        d["k"] = max(d["k"],float(dif[:,:,2].max())/model.ref.Kin)
        d["Pi"] = max(d["Pi"],float(dif[:,:,4].max())/model.ref.tau)
        d["G"] = max(d["G"],float(dif[:,:,6:].max()))
        radii = np.tile(r,len(x_two))
        outer = bool(e["epsilon"])
        if outer:
            d["epsilon"] = max(d["epsilon"],float((dif[:,:,3]*r[None,:]/model.ref.A).max()))
        orig_R, orig_T = stresses(model.ref,old.reshape(-1,10),radii,outer)
        merge_R, merge_T = stresses(model.ref,merged.reshape(-1,10),radii,outer)
        d["R"] = max(d["R"],float(abs(orig_R-merge_R).max())/model.ref.tau)
        d["traction"] = max(d["traction"],float(abs(orig_T-merge_T).max())/model.ref.tau)
    return d


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("table","geometry","native","quiet-witness","out"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--max-cells",type=int,required=True)
    p.add_argument("--x-levels",type=int,required=True)
    p.add_argument("--wall-levels",type=int,required=True)
    p.add_argument("--collar",action="store_true")
    p.add_argument("--relocate-tip",action="store_true",
                   help="Shorten the two corner cells by moving the adjacent grid lines.")
    p.add_argument("--seconds",type=float,required=True)
    a=p.parse_args()
    end=time.monotonic()+a.seconds
    ref=po.load_reference(a.table)
    b=load_state(a.native,ref,20,end)
    xs,zs=axes(a.geometry)
    bx,br,owners=cell_index(b.elements)
    dx,dr = xs[-1]-xs[-2],zs[1]-zs[0]
    if a.relocate_tip:
        po.require(a.x_levels>0 and a.wall_levels>0 and not a.collar,
                   "the isolated corner control moves only the two adjacent grid lines")
        x,z = list(xs),list(zs)
        x[-2] = xs[-1]-dx/2**a.x_levels
        z[1] = zs[0]+dr/2**a.wall_levels
    else:
        x = sorted(xs+[xs[-1]-dx/2**j for j in range(1,a.x_levels+1)])
        z = sorted(zs+[zs[0]+dr/2**j for j in range(1,a.wall_levels+1)])
    if a.collar:
        i = min(range(1,len(z)-1),key=lambda j:abs(z[j]-ref.d0-.1))
        po.require(abs(z[i]-ref.d0-.1)<1e-12,"collar witness is not a native edge")
        z = sorted(z+[(z[i-1]+z[i])/2,(z[i]+z[i+1])/2])
    allowed=a.max_cells//(len(z)-1)
    need=max(0,len(x)-1-allowed)
    witness=json.loads(a.quiet_witness.read_text())
    quiet=[v["x"] for v in witness["inner_k_bands"] if v["x"][0]>=-8.75 and v["x"][1]<=4.75
           and v["RMS_over_Kin"]<=1.1e-4]
    options = []
    for l,r in zip(bx[:-1],bx[1:]):
        if need and any(lo<=l[0] and r[1]<=hi for lo,hi in quiet):
            errs=merge_errors(b,l,r,br,owners)
            score=max(errs[key]/lim for key,lim in po.NUMERICAL.items())
            options.append(dict(x=[l[0],l[1],r[1]], projected_MAX=errs, fraction_of_limit=score))
    removed=[]
    used=set()
    for item in sorted(options,key=lambda row:row["fraction_of_limit"]):
        lo,mid,hi=item["x"]
        if len(removed)<need and item["fraction_of_limit"]<=.20 and not ({lo,mid,hi}&used):
            removed.append(item)
            used.update((lo,mid,hi))
    po.require(len(removed)==need,f"found {len(removed)} safe independent quiet-strip merges for {need} required")
    banned={round(row["x"][1],12) for row in removed}
    x=[v for v in x if round(v,12) not in banned]
    po.require((len(x)-1)*(len(z)-1)<=a.max_cells,"resource design bound exceeded")
    a.out.mkdir()
    template=a.out/"geometry.xml"
    mesh(template,x,z,quadrature=12)
    atomic(a.out/"forecast.json",dict(schema="WJH27Task2RSJointMeshV1",geometry=identity(template),
        source_mesh=identity(a.geometry),seed=identity(a.native),quiet_observation=identity(a.quiet_witness),
        table=identity(a.table),cells=[len(x)-1,len(z)-1],terminal_dx=x[-1]-x[-2],wall_dr=z[1]-z[0],
        x=x,r=[v-ref.d0 for v in z],max_cells=a.max_cells,removed=removed,
        screened_merges=options,collar_h_split=a.collar,checked_fraction=.20,
        relocated_tip=a.relocate_tip,
        previous_terminal_dx_dr=[dx,dr],previous_adjacent_dx_dr=[xs[-2]-xs[-3],zs[2]-zs[1]],
        next_adjacent_dx_dr=[x[-2]-x[-3],z[2]-z[1]],
        meaning="Resampling verifies representation on the quiet source fields. Fresh coupled equations and all owner observations determine the actual effect."))


if __name__=="__main__":
    main()
