#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Screen prescribed-amplitude sign and support alternatives at Bf owners."""

import argparse
from itertools import combinations
from pathlib import Path
import re
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import BSpline
from scipy.optimize import minimize_scalar

import perturbation_certificate as pc
import physical_oracle as po
from resources import atomic, identity
from task2r_advancing import read, require
from task2s_local_audit import angular_direction, principal_tensor


def extended_layout(old, radial):
    # First native edge giving at least 2.0 square code units, with the
    # original centered axial treatment and quartile radial construction.
    result = dict(old, shape="BfExtendedQuartileV1")
    x = old["x_knots"]
    lower = old["r_knots"][0]
    upper = next(r for r in radial if (r-lower)*(x[-1]-x[0]) >= 2.)
    inside = [r for r in radial if lower < r < upper]
    targets = (.25,.5,.75)
    knots = min(combinations(inside,3),key=lambda group:
        sum(((z-lower)/(upper-lower)-t)**2 for z,t in zip(group,targets)))
    H = old["H"]
    peak = min((r for r in inside if r > H),key=lambda r:abs(r-(H+upper)/2))
    result.update(r_knots=[lower,*knots,upper],epsilon_peak_radius=peak)
    return result


def inventory(base, layout, rejected):
    x0,x1 = layout["x_knots"][::4]
    top = layout["r_knots"][-1]
    point_lists = {}
    qx = (leggauss(6)[0]+1)/2
    qr = (leggauss(8)[0]+1)/2
    nodes = np.linspace(0.,1.,5)
    for owner in base.elements:
        xl,xr,rl,rr = owner["bounds"]
        if xr < x0 or xl > x1 or rl > top:
            continue
        gauss_x,gauss_r = np.meshgrid(xl+(xr-xl)*qx,rl+(rr-rl)*qr)
        grid_x,grid_r = np.meshgrid(xl+(xr-xl)*nodes,rl+(rr-rl)*nodes)
        x = np.r_[gauss_x.ravel(),grid_x.ravel()]
        r = np.r_[gauss_r.ravel(),grid_r.ravel()]
        if owner["geometry"] == rejected["owner"]:
            x,r = np.r_[x,rejected["x"]],np.r_[r,rejected["r"]]
        values = base.sample(owner,x,r)
        lift = (rr-r)/(rr-rl) if owner["epsilon"] and rl == base.ref.H else np.zeros(len(r))
        record = dict(x=x,r=r,owner=np.full(len(r),owner["geometry"]),lift=lift,
                      outer=np.full(len(r),bool(owner["epsilon"])))
        record.update({key:values[key] for key in ("k","E","K","G","chi")})
        for key,value in record.items():
            point_lists.setdefault(key,[]).append(value)
    return {key:np.concatenate(blocks) for key,blocks in point_lists.items()}


def hermite(r, low, peak, high):
    t = np.clip(np.where(r <= peak,(r-low)/(peak-low),(high-r)/(high-peak)),0.,1.)
    return t*t*(3-2*t)


def velocity_shape(layout, points, normalization):
    spline = [BSpline.basis_element(layout[key],extrapolate=False) for key in ("x_knots","r_knots")]
    x,r = points["x"],points["r"]
    a,b = ([np.nan_to_num(s(v,nu=d)) for d in (0,1,2)] for s,v in zip(spline,(x,r)))
    G = np.zeros_like(points["G"])
    G[:,0,0] = a[1]*b[1]
    G[:,0,2] = a[0]*b[2]
    G[:,2,0] = -a[2]*b[0]
    G[:,2,2] = -G[:,0,0]
    return 1.5*a[0],G*normalization


def candidate(ref, layout, points, coefficients, signs, rejected):
    shape_x, delta_g = coefficients
    velocity_sign,k_sign,epsilon_sign = signs
    r = points["r"]
    lower,upper = layout["r_knots"][::4]
    H = ref.H
    G = points["G"]+velocity_sign*delta_g
    k = points["k"]+k_sign*.01*ref.Kin*shape_x*hermite(r,lower,H,upper)
    K = points["K"]+k_sign*.01*ref.Kin*shape_x
    E = K**1.5/ref.cl/np.maximum(r,ref.cutoff(K))
    outer = points["outer"]
    e = points["E"][outer]
    gamma = (K[outer]**1.5-points["K"][outer]**1.5)/(ref.cl*H)
    bump = hermite(r[outer],H,layout["epsilon_peak_radius"],upper)
    E[outer] = e+gamma*points["lift"][outer]+epsilon_sign*.01*ref.A/r[outer]*shape_x[outer]*bump
    chi = points["chi"]
    root,tensor = np.full(len(r),1.8),2*k[:,None,None]*np.eye(3)/3
    inner = ~outer
    shear = ref.aux(K[inner],r[inner])*(G[inner,0,2]+G[inner,2,0])
    tensor[inner,0,2] = tensor[inner,2,0] = -shear
    active = chi > 0
    wj,root[active],root_error,_ = po.full_wj(k[active],E[active],G[active])
    tensor[active] = (1-chi[active,None,None])*tensor[active]+chi[active,None,None]*wj
    cov = np.linalg.eigvalsh(tensor)[:,0]/k
    args = (ref,k,E,K,r,G,chi,root)
    matrices,mu = principal_tensor(*args,step=1e-5)
    angles = np.arange(128)*np.pi/128
    directions = angular_direction(angles)
    symbol = np.einsum("ti,qij,tj->qt",directions,matrices,directions,optimize=True)
    worst_angle = np.argmin(symbol,axis=1)
    values = symbol[np.arange(len(r)),worst_angle]
    selected = np.unique(np.r_[np.argsort(values)[:20],
        np.flatnonzero((points["owner"] == rejected["owner"]) &
                       (points["x"] == rejected["x"]) & (r == rejected["r"]))])
    refined,_ = principal_tensor(ref,*(v[selected] for v in args[1:]),step=5e-6)
    findings = []
    for j,i in enumerate(selected):
        centre = angles[worst_angle[i]]
        def action(angle,matrix):
            d = angular_direction(angle)
            return float(d@matrix@d)
        sample = minimize_scalar(lambda t:action(t,refined[j]),method="bounded",
            bounds=(centre-np.pi/128,centre+np.pi/128),options={"xatol":1e-11})
        findings.append(dict(principal_over_mu=float(sample.fun),x=float(points["x"][i]),
            r=float(r[i]),owner=int(points["owner"][i]),half_step_change=abs(sample.fun-action(sample.x,matrices[i])),
            k=float(k[i]),E=float(E[i]),Gxz=float(G[i,0,2]),sampled_K=float(K[i])))
    worst = min(findings,key=lambda row:row["principal_over_mu"])
    anchor = next(row for row in findings if row["owner"] == rejected["owner"]
                  and row["x"] == rejected["x"] and row["r"] == rejected["r"])
    return dict(signs=list(signs),layout=layout,points=len(r),worst=worst,rejection_owner=anchor,
        minimum_covariance_over_k=float(np.min(cov)),maximum_root_error=float(np.max(root_error)),
        minimum_auxiliary_principal=float(np.min(ref.nu/mu[~active])),
        promising=bool(worst["principal_over_mu"] > max(1e-5,10*worst["half_step_change"]) and
                       np.min(cov) > 1e-10))


def screen(base_path, table, certificate, rejection_log, out, seconds):
    began = time.monotonic()
    deadline = began+seconds
    proof = read(certificate)
    match = re.search(r"admissibility failed: geometry=(\d+), x=([^,]+), r=([^,]+)",rejection_log.read_text())
    require(match,"native failure location is missing")
    rejected = dict(owner=int(match[1]),x=float(match[2]),r=float(match[3]))
    data,_ = po.load_native(base_path)
    ref = po.load_reference(table)
    base = po.PhysicalOracle(data,ref,16,deadline)
    require(proof["source_identity"] == data["source_identity"]
            and proof["base_state_sha256"] == data["state_sha256"],"shape pre-screen base differs")
    old = proof["normalization"]["layout"]
    results = []
    for layout in (old,extended_layout(old,proof["r_edges"])):
        points = inventory(base,layout,rejected)
        peak = pc.maximum(pc.shape_pieces(layout,"velocity"),deadline)
        cen = np.sqrt(.5*(peak["upper"]+peak["lower"]))
        coefficients = velocity_shape(layout,points,.001/cen)
        for signs in ((1,1,1),(1,-1,1),(-1,1,1),(-1,-1,1)):
            base.clock()
            results.append(candidate(ref,layout,points,coefficients,signs,rejected))
    atomic(out,dict(schema="WJH27R3BfPerturbationPreScreenV1",base=identity(base_path),
        certificate=identity(certificate),rejected=identity(rejection_log),table=identity(table),
        elapsed_seconds=time.monotonic()-began,results=results,
        scope="Pointwise analytical candidate screen on the actual Bf field. The outer epsilon bump is the "
              "unprojected target and the dependent change is lifted on the first outer owner. A promising "
              "screen still needs the native projected state, exact signed norm and full guard checks."))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("base_path","table","certificate","rejection_log","out"):
        p.add_argument("--"+key.replace("_","-"),type=Path,required=True)
    p.add_argument("--seconds",type=float,required=True)
    screen(**vars(p.parse_args()))
