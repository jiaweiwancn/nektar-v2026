#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Audit current continuity, natural trace actions and planar margins in Task 2R-S."""

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar
from scipy.sparse.linalg import splu

import physical_oracle as po
from resources import atomic, identity
from task2s_observations import load_state, polynomial_basis


PLANE = (0, 2)
TEST_GRADIENTS = np.zeros((3,3,3))
TEST_GRADIENTS[0,0,0], TEST_GRADIENTS[0,2,2] = 1., -1.
TEST_GRADIENTS[1,0,2], TEST_GRADIENTS[2,2,0] = 1., 1.


def quadrature(order):
    x,w = leggauss(order)
    return (x+1)/2,w/2


def continuity(state):
    points, weights = quadrature(6)
    volume_map = np.outer(weights,weights).ravel()
    mass, weak = po.SparseMass(state.np),np.zeros(state.np)
    pressure_shift = np.zeros((2,state.n))
    volume, absolute, squared, maximum = np.zeros(4)
    location = None
    for element in state.elements:
        state.clock()
        xl,xr,rl,rr = element["bounds"]
        x = np.repeat(xl+(xr-xl)*points,len(points))
        r = np.tile(rl+(rr-rl)*points,len(points))
        q = volume_map*(xr-xl)*(rr-rl)
        p = polynomial_basis(element,"p",x,r)[0]
        b, dx, dr = polynomial_basis(element,"u",x,r)
        ids = np.asarray(element["u"])
        divergence = dx@state.s["flow"][ids]+dr@state.s["flow"][state.n+ids]
        pid = np.asarray(element["p"])
        mass.add(pid,p.T@(q[:,None]*p))
        weak[pid] += p.T@(q*divergence)
        pressure_shift[0,ids] -= q@dx
        pressure_shift[1,ids] -= q@dr
        volume += float(q.sum())
        absolute += float(q@abs(divergence))
        squared += float(q@(divergence*divergence))
        peak = int(np.argmax(abs(divergence)))
        if abs(divergence[peak]) > maximum:
            maximum = float(abs(divergence[peak]))
            location = dict(owner=element["geometry"],x=float(x[peak]),r=float(r[peak]))
    y = splu(mass.finish()).solve(weak)
    return dict(pressure_tested_dual=math.sqrt(max(0.,float(weak@y)))/volume,
                pressure_test_action_MAX=float(np.max(abs(weak))),
                strong_L1_over_area=absolute/volume,strong_RMS=math.sqrt(squared/volume),
                strong_interior_Gauss_MAX=maximum,strong_interior_peak=location,area=volume),pressure_shift


def weak_traces(state, shift):
    from_faces = np.zeros_like(shift)
    prescribed, reproduced = np.zeros_like(shift),np.zeros_like(shift)
    boundary_norms = po.Norms()
    wall_and_outlet_owners = []
    corner = None
    for element, side in state.faces["u"]:
        n = np.array(po.NORMALS[side])
        x,r,q = state.face(element,side,endpoints=True)
        b = polynomial_basis(element,"u",x,r)[0]
        ids = np.asarray(element["u"])
        from_faces[:,ids] -= n[:,None]*(q@b)[None,:]
        if side not in (1,2):
            continue
        a = state.sample(element,x,r)
        V = a["V"][:,PLANE,:][:,:,PLANE]
        actual = V@n-a["Pi"][:,None]*n
        if side == 1:
            reference = state.ref.outlet(r)
            for comp in (0,1):
                prescribed[comp,ids] -= (q*reference[:,comp])@b
                reproduced[comp,ids] -= (q*actual[:,comp])@b
            boundary_norms.add("outlet","raw_minus_prescribed_traction",(actual-reference)/state.ref.tau,q)
        boundary_norms.add(po.SIDES[side],"D",a["div"],q)
        if side == 2:
            boundary_norms.add("ground","tangential_derivatives",a["G"][:,PLANE,0],q)
            boundary_norms.add("ground","Wr",a["G"][:,2,2],q)
        i = int(np.argmax(abs(a["div"])))
        wall_and_outlet_owners.append(dict(side=po.SIDES[side], owner=element["geometry"],
            bounds=element["bounds"], D_RMS=math.sqrt(float(q@(a["div"]**2))/q.sum()),
            D_at_MAX=float(a["div"][i]),x=float(x[i]),r=float(r[i])))
        if side == 1 and element["bounds"][2] == 0.:
            j = int(np.argmin(r))
            oracle_values = []
            for comp in (0,1):
                keep = np.isin(ids,state.ufree[comp])
                oracle_values.append(float(np.max(abs(b[j,keep]))))
            corner = dict(owner=element["geometry"],free_test_trace_MAX=oracle_values,
                x=float(x[j]),r=float(r[j]),Pi_over_tau=float(a["Pi"][j]/state.ref.tau),
                raw_traction_over_tau=(actual[j]/state.ref.tau).tolist(),
                prescribed_traction_over_tau=(reference[j]/state.ref.tau).tolist(),
                D=float(a["div"][j]),G=a["G"][j][:,PLANE][PLANE,:].tolist())
    err = [float(np.max(abs(shift[c,ids]-from_faces[c,ids]))) for c,ids in enumerate(state.ufree)]
    max_action = [float(np.max(abs(shift[c,ids]))) for c,ids in enumerate(state.ufree)]
    natural = [float(np.max(abs(prescribed[c,ids]))) for c,ids in enumerate(state.ufree)]
    defect = [float(np.max(abs(reproduced[c,ids]-prescribed[c,ids]))) for c,ids in enumerate(state.ufree)]
    po.require(corner and max(corner["free_test_trace_MAX"]) < 1e-12, "shared no-slip corner has a nonzero admissible trace")
    po.require(max(err) < 1e-11 and max_action[0] > 1e-6, "open normal trace fails pressure level test")
    return dict(plane_normals={"outlet":[1,0],"ground":[0,-1]},
        restricted_pressure_constant_action_MAX=max_action,
        restricted_pressure_volume_minus_face_MAX=err,
        pressure_constant_AWF_derivative="Zero derivative of represented constant at H",
        restricted_prescribed_outlet_residual_action_MAX=natural,
        restricted_raw_traction_minus_prescribed_action_MAX=defect,
        raw_vertex=corner,boundary_norms=boundary_norms.report({}),
        wall_and_outlet_owners=wall_and_outlet_owners)


def angular_direction(angle):
    c,s = np.cos(angle),np.sin(angle)
    return np.stack((-s*c,-s*s,c*c),axis=-1)


def principal_tensor(ref,k,E,K,r,G,chi,root,step):
    nu = ref.nu
    kE = k/E
    rotation = kE[:,None,None]*(G-G.transpose(0,2,1))/2
    iw = np.einsum("qij,qji->q",rotation,rotation)
    aux = np.where(r < ref.H,ref.aux(K,np.minimum(r,ref.H)),0.)
    mu = nu+(1-chi)*aux+chi*.6*root*k*kE/(root*root-2*iw)
    active = chi > 0
    matrix = np.zeros((len(r),3,3))
    h = step*np.maximum(1.,np.max(abs(G[active]),axis=(1,2)))
    G_active = G[active]
    for j,direction in enumerate(TEST_GRADIENTS):
        symmetric = direction+direction.T
        delta = nu*np.broadcast_to(symmetric,G.shape).copy()
        shear = (1-chi)*aux*symmetric[0,2]
        delta[:,0,2] += shear
        delta[:,2,0] += shear
        if np.any(active):
            perturb = h[:,None,None]*direction
            plus = po.full_wj(k[active],E[active],G_active+perturb)[0]
            minus = po.full_wj(k[active],E[active],G_active-perturb)[0]
            delta[active] -= chi[active,None,None]*(plus-minus)/(2*h[:,None,None])
        matrix[:,:,j] = np.einsum("qmn,imn->qi",delta,TEST_GRADIENTS)/mu[:,None]
    return matrix,mu


def point_inventory(state):
    fields = ("k","E","K","G","chi","root","covariance","div")
    pieces = {key:[] for key in fields+("x","r","owner")}
    standard = np.array([0.,.5,1.])
    dense = np.linspace(0.,1.,7)
    for owner in state.elements:
        xl,xr,rl,rr = owner["bounds"]
        fine = xr >= 14.80 and rl < .125
        points = dense if fine else standard
        x = np.repeat(xl+(xr-xl)*points,len(points))
        r = np.tile(rl+(rr-rl)*points,len(points))
        obs = state.sample(owner,x,r)
        pieces["x"].append(x)
        pieces["r"].append(r)
        pieces["owner"].append(np.full(len(x),owner["geometry"],dtype=int))
        for key in fields:
            pieces[key].append(obs[key])
    return {key:np.concatenate(blocks) for key,blocks in pieces.items()}


def admissibility(state):
    points = point_inventory(state)
    k,E,K,r,G,chi,root = (points[k] for k in ("k","E","K","r","G","chi","root"))
    args = (state.ref,k,E,K,r,G,chi,root)
    m,mu = principal_tensor(*args,step=1e-5)
    count=64
    theta = np.arange(count)*np.pi/count
    v = angular_direction(theta)
    p = np.einsum("ti,qij,tj->qt",v,m,v,optimize=True)
    indices = np.argmin(p,axis=1)
    minima = p[np.arange(len(r)),indices]
    # Exact across all directions in the auxiliary-only region.
    minima[chi == 0] = state.ref.nu/mu[chi == 0]
    x = points["x"]
    groups = dict(whole=np.ones(len(r),dtype=bool),lower_aux=chi==0,
        corner_box=(x>=14.96)&(r<=.04),outlet_collar=(x>=14.8)&(.085<=r)&(r<=.115))
    cov, div = points["covariance"],points["div"]
    output = {}
    for key,mask in groups.items():
        ids = np.flatnonzero(mask)
        if not len(ids):
            continue
        # Refine the strongest sampled spatial locations in angle and repeat the
        # constitutive tangent with half the perturbation at those locations.
        chosen = ids[np.argsort(minima[ids])[:24]]
        mm,_ = principal_tensor(state.ref,k[chosen],E[chosen],K[chosen],r[chosen],
                                G[chosen],chi[chosen],root[chosen],step=5e-6)
        angular = []
        for j,i in enumerate(chosen):
            centre = theta[indices[i]] if chi[i]>0 else np.pi/4
            def form(a, mat):
                d = angular_direction(a)
                return float(d@mat@d)
            solved = minimize_scalar(lambda a:form(a,mm[j]),bounds=(centre-np.pi/count,centre+np.pi/count),
                                     method="bounded",options={"xatol":1e-11})
            angular.append((solved.fun, i, solved.x, abs(solved.fun-form(solved.x,m[i]))))
        smallest,i,angle,fd = min(angular)
        j = ids[np.argmin(cov[ids])]
        d = ids[np.argmax(abs(div[ids]))]
        def locate(i,value):
            return dict(value=float(value), x=float(x[i]), r=float(r[i]), owner=int(points["owner"][i]),
                        D=float(div[i]), k=float(k[i]), E=float(E[i]), K=float(K[i]))
        output[key] = dict(points=len(ids),minimum_planar_over_mu=locate(i,smallest),
                           direction=angle,half_step_change=fd,minimum_covariance_over_k=locate(j,cov[j]),
                           largest_sampled_D=locate(d,div[d]))
    lower = chi==0
    support = np.flatnonzero(lower)
    margin = 2*k[lower]/3-state.ref.aux(K[lower],r[lower])*abs(G[lower,0,2]+G[lower,2,0])
    n = int(np.argmin(margin/k[lower]))
    i = support[n]
    return dict(pointwise_sample_only=True,points=len(r),regions=output,
        native_reports={k:state.data[k] for k in ("minimum_covariance_over_k","minimum_principal_over_mu")},
        auxiliary_shear_covariance_reserve_over_k_min=dict(value=float(margin[n]/k[i]),
            x=float(x[i]),r=float(r[i]),owner=int(points["owner"][i])))


def compare_collar(first,second,report):
    if report is None:
        return None
    observations = json.loads(report.read_text())["rows"]
    key = "x=15/inner/R"
    peak = observations[key].get("peak")
    if peak is None:
        return dict(row=key,status=observations[key]["status"])
    x,r = (np.array([peak[k]]) for k in ("x","r"))
    owners = [next(e for e in a.elements if e["geometry"]==peak[k])
              for a,k in ((first,"coarse_owner"),(second,"fine_owner"))]
    a,b = [s.sample(e,x,r) for s,e in zip((first,second),owners)]
    def summary(v):
        return dict(k=float(v["k"][0]),epsilon=float(v["E"][0]),K=float(v["K"][0]),
            Ryy_over_tau=float(v["R"][0,1,1]/first.ref.tau),D=float(v["div"][0]),
            raw_gradient=v["G"][0][:,PLANE][PLANE,:].tolist(),chi=float(v["chi"][0]),
            epsilon_ratio_k_over_E=float(v["k"][0]/v["E"][0]),
            pressure_over_tau=float(v["Pi"][0]/first.ref.tau))
    def yy(k,E,G,chi):
        return float(((1-chi)*2*k/3+chi*po.full_wj(k,E,G)[0][:,1,1])[0]/first.ref.tau)
    def remove_div(v):
        return v["G"]-v["div"][:,None,None]*np.diag([.5,0.,.5])
    path = [yy(a["k"],a["E"],a["G"],a["chi"]),yy(b["k"],a["E"],a["G"],a["chi"]),
            yy(b["k"],b["E"],a["G"],a["chi"]),yy(b["k"],b["E"],b["G"],b["chi"])]
    filtered = yy(b["k"],b["E"],remove_div(b),b["chi"])-yy(a["k"],a["E"],remove_div(a),a["chi"])
    return dict(peak=peak,first=summary(a),second=summary(b),diagnostic_only=True,
        raw_fine_minus_coarse=float(path[-1]-path[0]),change_from_k_then_E_then_G=np.diff(path).tolist(),
        delta_using_trace_free_planar_gradients=filtered,
        interpretation="Telescoping finite substitutions depend on the order. All governing observations use the original gradients.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table",type=Path,required=True)
    p.add_argument("--native",type=Path,action="append",required=True)
    p.add_argument("--comparison",type=Path)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--seconds",type=float,required=True)
    a=p.parse_args()
    deadline=time.monotonic()+a.seconds
    ref=po.load_reference(a.table)
    fields = {}
    for path in a.native:
        state=load_state(path,ref,24,deadline)
        balance,shift=continuity(state)
        fields[path.parent.name]=dict(state=identity(path),continuity=balance,
            weak_natural_boundary=weak_traces(state,shift),samples=admissibility(state))
    collar = None
    if a.comparison:
        first,second = (load_state(path,ref,24,deadline) for path in a.native[-2:])
        collar=compare_collar(first,second,a.comparison)
    atomic(a.out,dict(schema="WJH27Task2RSLocalAuditV1",table=identity(a.table),
                      cases=fields,collar=collar))


if __name__=="__main__":
    main()
