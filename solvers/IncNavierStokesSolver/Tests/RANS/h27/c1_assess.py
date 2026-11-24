#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reduce C1 metric, rotation, spatial and common-time controls."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.chebyshev import chebder, chebval
from numpy.polynomial.legendre import leggauss

from c1_geometry import Hill
from c1_mms import Manufactured
import physical_oracle as po
from resources import identity


FIELDS = ("velocity", "pressure", "k", "epsilon")


def geometry_check(path, reference):
    data=json.loads(path.read_text())
    hill=Hill(dict(left=-10,right=15,H=.15,d0=.0875,tunnel_top=20.,
                   hill_half_length=2.5,hill_height=data["height"]))
    defects={key: 0. for key in ("position","inverse","normals","jacobian",
                                  "gradient","divergence","matching_measure")}
    for row in data["points"]:
        a,r=np.array([row["a"]]),np.array([row["rho"]])
        m=hill.metrics(row,a,r)
        g=(.23+.07*r)*m["da"]+(.07*a-.24*r)*m["dr"]
        div=g[0,0]+(.054*a*m["da"][:,1]+.11*m["dr"][:,1])[0]
        expected=dict(position=m["position"][0], inverse=np.vstack([m["da"],m["dr"]]),
                      normals=m["n"][0], jacobian=m["det"][0],gradient=g[0],divergence=div)
        inverse=np.asarray(row["inverse"])[[0,2]][:,[0,2]]
        actual=dict(position=np.asarray(row["position"])[[0,2]],inverse=inverse,
                    normals=np.asarray(row["normal"])[[0,2]],jacobian=row["det"],
                    gradient=row["grad"],divergence=row["divergence"])
        for key in expected:
            defects[key]=max(defects[key],float(np.max(np.abs(expected[key]-actual[key]))))
    nodes,weights=leggauss(32)
    for face in data["faces"]:
        lo,hi=face["bounds"][:2]
        a=lo+(hi-lo)*(nodes+1)/2
        length,_=hill.face(face,face["side"],a,np.full_like(a,hill.H))
        expected=float(weights@length)*(hi-lo)/2
        defects["matching_measure"]=max(defects["matching_measure"],abs(expected-face["measure"]))
    rotation={key: 0. for key in ("R","V","traction","production","covariance","principal","principal_raw")}
    base={row["radius"]:row for row in data["rotation"] if row["angle"]==0}
    for row in data["rotation"]:
        theta=row["angle"]
        c,s=math.cos(theta),math.sin(theta)
        Q=np.array([[c,0,-s],[0,1,0],[s,0,c]])
        for key in rotation:
            value=base[row["radius"]][key]
            if key in ("R","V"):
                expected=Q@np.asarray(value)@Q.T
            elif key=="traction":
                expected=Q@value
            else:
                expected=value
            rotation[key]=max(rotation[key],float(np.max(np.abs(expected-row[key]))))
    derivative=chebder(reference.coefficients)*2/(reference.hi-reference.lo)
    cutoff={key: 0. for key in ("c","dc","EwK_below","EwK_above","traceK")}
    for row in data["cutoff_derivatives"]:
        K=row["K"]
        x=(2*K-reference.lo-reference.hi)/(reference.hi-reference.lo)
        c=reference.H*math.exp(chebval(x,reference.coefficients))
        log_slope=chebval(x,derivative)
        A=K**1.5/reference.cl
        expected=dict(c=c,dc=c*log_slope,EwK_below=A/c*(1.5/K-log_slope),
                      EwK_above=1.5*A/(.04*K),traceK=1.5*A/(reference.H*K))
        for key in cutoff:
            cutoff[key]=max(cutoff[key],abs(expected[key]-row[key])/abs(expected[key]))
    limits=dict(metrics=2e-11,rotation=5e-12,cutoff_derivative=1e-8)
    return dict(native=identity(path),height=hill.height,certificate=hill.certificate(),
                defects=defects,rotation_defects=rotation,cutoff_relative_defects=cutoff,
                limits=limits,passed=(max(defects.values())<=limits["metrics"] and
                    max(rotation.values())<=limits["rotation"] and max(cutoff.values())<=limits["cutoff_derivative"]))


def spatial(paths):
    rows=[]
    for p in paths:
        d=json.loads(p.read_text())
        rows.append(dict(file=identity(p),error=d["error"],
                         raw_traction=d["interface"]["raw_traction_jump_rms"],
                         raw_k_flux=d["interface"]["raw_k_diffusion_jump_rms"],
                         weak_divergence_work=d["weak_divergence_work"]))
    rates=[]
    for a,b in zip(rows,rows[1:]):
        rates.append({key: math.log2(a["error"][key]/b["error"][key]) for key in FIELDS})
    return dict(rows=rows,refinement_rates=rates,
                nominal_smooth_L2_rates=dict(velocity=4,pressure=3,k=4,epsilon=4),
                gradient_rate=3, chart_geometry="analytic at every quadrature point")


def field_norms(data, ref, other=None, order=20):
    geometry=data["geometry"]
    hill=Hill(geometry)
    exact=Manufactured(hill,ref)
    n,n_p,*_=data["sizes"]
    state={key:np.asarray(value) for key,value in data["coefficients"].items()}
    if other is not None:
        po.require(other["geometry"]==geometry and other["elements"]==data["elements"],
                   "temporal differences require identical native spaces")
        state={key:value-np.asarray(other["coefficients"][key]) for key,value in state.items()}
    nodes,weights=leggauss(order)
    err=dict.fromkeys(FIELDS,0.)
    mean=total=outer=0.
    for e in data["elements"]:
        xl,xr,rl,rr=e["bounds"]
        a,r=np.meshgrid(xl+(xr-xl)*(nodes+1)/2,rl+(rr-rl)*(nodes+1)/2)
        a,r=a.ravel(),r.ravel()
        metric=hill.metrics(e,a,r)
        w=np.outer(weights,weights).ravel()*(xr-xl)*(rr-rl)/4*metric["det"]
        fields=exact.local(*exact.to_tube(e,a,r),data["time"])
        xy=exact.rotate(fields,"u")
        bu=po.basis(e,"u",a,r)[0]; ids=np.asarray(e["u"])
        u=bu@state["flow"][ids]
        v=bu@state["flow"][n+ids]
        p=po.basis(e,"p",a,r)[0]@state["flow"][2*n+np.asarray(e["p"])]
        k=po.basis(e,"k",a,r)[0]@state["k"][e["k"]]
        if other is None:
            u-=xy[:,0];v-=xy[:,1];p-=fields["Pi"];k-=fields["k"]
        err["velocity"]+=float(w@(u*u+v*v))
        err["pressure"]+=float(w@(p*p))
        err["k"]+=float(w@((k/ref.Kin)**2))
        mean+=float(w@p)
        total+=float(w.sum())
        if e["epsilon"]:
            eps=po.basis(e,"epsilon",a,r)[0]@state["epsilon"][e["epsilon"]]
            if other is None:
                eps-=fields["epsilon"]
            err["epsilon"]+=float(w@((eps/fields["epsilon"])**2))
            outer+=float(w.sum())
    err["pressure"]=max(0.,err["pressure"]-mean*mean/total)
    return {key:math.sqrt(value/(outer if key=="epsilon" else total))
            for key,value in err.items()}


def temporal(case_dirs,ref):
    rows=[]
    snapshots=[]
    for p in case_dirs:
        data,meta=po.load_native(p/"state.result.json")
        dt=data["dt"]
        frames,raw=[],[]
        for t in (.16,.32):
            f=p/f"state.science-{round(t/dt)}.json"
            frame,_=po.load_native(f)
            po.require(abs(frame["time"]-t)<1e-13,"different common positive clock")
            frames.append(dict(time=t,step=frame["step"],file=identity(f),
                               error=field_norms(frame,ref),native_residuals=frame["residuals"]))
            raw.append(frame)
        rows.append(dict(dt=dt,frames=frames,timing=meta))
        snapshots.append(raw)
    comparisons=[]
    for a,b in zip(rows,rows[1:]):
        comparisons.append({t["time"]: {key:math.log(t["error"][key]/u["error"][key],a["dt"]/b["dt"])
                            for key in FIELDS} for t,u in zip(a["frames"],b["frames"])})
    differences=[{a["time"]:field_norms(a,ref,b) for a,b in zip(older,newer)}
                 for older,newer in zip(snapshots,snapshots[1:])]
    self_rates=[{t:{key:math.log2(a[t][key]/b[t][key]) for key in FIELDS} for t in a}
                for a,b in zip(differences,differences[1:])]
    return dict(rows=rows,rates=comparisons,same_space_differences=differences,
                same_space_rates=self_rates)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table",type=Path,required=True)
    parser.add_argument("--geometry",type=Path,nargs="*")
    parser.add_argument("--spatial",type=Path,nargs="*")
    parser.add_argument("--temporal",type=Path,nargs="*")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    ref=po.load_reference(args.table)
    result=dict(table=identity(args.table))
    if args.geometry:
        result["geometry"]=[geometry_check(p,ref) for p in args.geometry]
    if args.spatial:
        result["spatial"]=spatial(args.spatial)
    if args.temporal:
        result["temporal"]=temporal(args.temporal,ref)
    with args.out.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False)
        stream.write("\n")


if __name__=="__main__":
    main()
