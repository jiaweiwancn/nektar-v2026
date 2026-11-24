#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independently integrate the WJ normal-hill equations and budgets."""

import argparse
from collections import defaultdict
import json
import math
from pathlib import Path
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import brentq
from scipy.sparse.linalg import splu

from c1_geometry import Hill
from c1_mms import Manufactured
import physical_oracle as po
from physical_transient_oracle import scheme, state_hash
from resources import identity


def dot(a, b):
    return np.einsum("qi,qi->q", a, b)


def coefficients(data):
    return {key: np.array(value) for key, value in data.items()}


def rate(data, state):
    if data["step"] == 0:
        return {k: np.zeros_like(x) for k, x in state.items()}
    integrator = scheme(data)
    po.require(integrator["advection_mode"] == "Implicit", "current advection required")
    history = data["input_history"]
    old = coefficients(history["current"])
    older = coefficients(history["previous"])
    po.require(state_hash(old) == data["previous_sha256"], "history identity changed")
    dt = integrator["dt"]
    return {k: ((x-old[k]) if integrator["startup"] else
                1.5*(x-old[k])-.5*(old[k]-older[k]))/dt for k,x in state.items()}


class Audit:
    # These are independent polynomial and panel definitions from the flat
    # oracle. Their weights are in the chart and get physical metrics below.
    trace = po.PhysicalOracle.trace
    cuts = po.PhysicalOracle.cuts
    cutoff_x_splits = po.PhysicalOracle.cutoff_x_splits
    rule = po.PhysicalOracle.rule
    face = po.PhysicalOracle.face

    def __init__(self, data, ref, order, deadline):
        self.data, self.ref, self.order, self.deadline = data, ref, order, deadline
        self.s = coefficients(data["coefficients"])
        po.require(state_hash(self.s) == data["state_sha256"], "accepted coefficient hash differs")
        self.n, self.np, self.nk, self.ne = map(int,data["sizes"])
        self.mesh = data["elements"]
        self.geometry = data["geometry"]
        po.require(self.geometry["mapping"] == "NormalCosineHillV1" and
                   data["reference_identity"] == ref.identity, "wrong geometry or reference identity")
        self.hill = Hill(self.geometry)
        self.certificate = self.hill.certificate()
        self.exact = Manufactured(self.hill,ref) if data["manufactured"] else None
        self.top = self.hill.top-self.hill.d
        self.gauss = leggauss(order)
        self.columns = defaultdict(list)
        for e in self.mesh:
            self.columns[tuple(e["bounds"][:2])].append(e)
        self.owner, self.interfaces = {}, []
        for column in self.columns.values():
            inner, = [e for e in column if abs(e["bounds"][3]-ref.H)<1e-12]
            outer, = [e for e in column if abs(e["bounds"][2]-ref.H)<1e-12]
            po.require(not inner["epsilon"] and bool(outer["epsilon"]), "wrong epsilon ownership")
            self.interfaces.append((inner,outer))
            for e in column:
                self.owner[e["geometry"]] = outer
        self.free = {"u": np.asarray(data["free_velocity_ids"][0]),
                     "w": np.asarray(data["free_velocity_ids"][1]),
                     "p": np.arange(self.np),
                     "k": np.asarray(data["k_free_ids"]),
                     "epsilon": np.asarray(data["epsilon_free_ids"])}
        self.size = dict(u=self.n,w=self.n,p=self.np,k=self.nk,epsilon=self.ne)
        self.fixed = {k: np.setdiff1d(np.arange(self.size[k]), ids)
                      for k, ids in self.free.items()}
        self.one = {k: np.full(size,np.nan) for k,size in self.size.items()}
        for e in self.mesh:
            for key in ("u","p","k","epsilon"):
                if e[key]:
                    stride = 3 if key=="p" else 4
                    for local,(i,s) in enumerate(zip(e[key],e[po.SIGNS[key]])):
                        self.one[key][i] = s if local%stride<2 and local//stride<2 else 0.
        self.one["w"] = self.one["u"]
        po.require(all(np.isfinite(x).all() for x in self.one.values()), "incomplete signed map")
        self.time = data["time"]
        self.rates = rate(data,self.s)
        self.temporal = None
        if data["step"]:
            h = data["input_history"]
            self.temporal = (coefficients(h["current"]),coefficients(h["previous"]),
                             scheme(data)["mass_coefficients"])

    def clock(self):
        if time.monotonic()>self.deadline:
            raise TimeoutError("normal-hill audit reached its wall cap")

    def zeros(self):
        return {key: np.zeros(size) for key,size in self.size.items()}

    @staticmethod
    def basis(e, key, a, r, metric):
        b, da, dr = po.basis(e,key,a,r)
        g = da[:,:,None]*metric["da"][:,None,:]+dr[:,:,None]*metric["dr"][:,None,:]
        return b,g

    def velocity(self, b, ids, state):
        return np.column_stack((b@state["flow"][ids],b@state["flow"][self.n+ids]))

    def sample(self,e,a,r,metric):
        ref=self.ref
        basis={key: self.basis(e,key,a,r,metric) for key in ("u","p","k","epsilon") if e[key]}
        b, grad=basis["u"]
        ids=np.asarray(e["u"])
        v=self.velocity(b,ids,self.s)
        G=np.stack((np.einsum("qij,i->qj",grad,self.s["flow"][ids]),
                    np.einsum("qij,i->qj",grad,self.s["flow"][self.n+ids])),axis=1)
        k=basis["k"][0]@self.s["k"][e["k"]]
        K=self.trace(e,"k",a)
        cutoff=ref.cutoff(K)
        outer=bool(e["epsilon"])
        E=basis["epsilon"][0]@self.s["epsilon"][e["epsilon"]] if outer else (
            K**1.5/ref.cl/np.maximum(r,cutoff))
        M=G+G.transpose(0,2,1)
        normal=metric["n"]
        tangent=metric["t"]
        aux=np.zeros_like(E) if outer else ref.aux(K,r)
        bA=-aux*np.einsum("qi,qij,qj->q",tangent,M,normal)
        R=np.zeros((len(a),3,3))
        for i in range(3):
            R[:,i,i]=2*k/3
        wall=bA[:,None,None]*(tangent[:,:,None]*normal[:,None,:] +
                               normal[:,:,None]*tangent[:,None,:])
        for index,i in enumerate((0,2)):
            for column,j in enumerate((0,2)):
                R[:,i,j]+=wall[:,index,column]
        q=(r/ref.H-.4)/.35
        chi=np.where(q<=0,0.,np.where(q>=1,1.,q**3*(10-15*q+6*q*q)))
        mixed=chi>0
        if np.any(mixed):
            full_G=np.zeros_like(R)
            for index,i in enumerate((0,2)):
                for column,j in enumerate((0,2)):
                    full_G[:,i,j]=G[:,index,column]
            wj,*_=po.full_wj(k[mixed],E[mixed],full_G[mixed])
            R[mixed]=(1-chi[mixed,None,None])*R[mixed]+chi[mixed,None,None]*wj
        planar=R[:,(0,2),:][:,:,(0,2)]
        V=ref.nu*M-planar+2*k[:,None,None]*np.eye(2)/3
        P=-np.einsum("qij,qij->q",planar,G)
        Pi=basis["p"][0]@self.s["flow"][2*self.n+np.asarray(e["p"])]
        D=dict(k=ref.nu+ref.cd*k*k/(ref.sk*E),epsilon=ref.nu+ref.cd*k*k/(ref.se*E))
        scalars=dict(k=k,epsilon=E)
        gradient={key: np.einsum("qij,i->qj",basis[key][1],self.s[key][e[key]])
                  for key in D if key in basis}
        return basis,v,G,Pi,V,P,scalars,D,gradient,np.min(np.linalg.eigvalsh(R)[:,0]/k)

    def sampled(self,e,a,r,metric):
        owner=self.owner[e["geometry"]]
        section=np.full_like(a,self.ref.H)
        h=self.hill.metrics(owner,a,section,outer=False)
        b,da,_=po.basis(owner,"u",a,section)
        ids=np.asarray(owner["u"])
        U=self.velocity(b,ids,self.s)
        ua=self.velocity(da,ids,self.s)
        ut=dot(h["t"],U)
        dudl=dot(h["t"],ua)/h["matching_area"]-h["kh"]*dot(h["n"],U)
        pressure=self.trace(e,"p",a,1)/h["matching_area"]
        scalar=-2*self.trace(e,"k",a,1)/(3*h["matching_area"])
        return {(name): (q*(r<self.ref.h))[:,None]*metric["t"] for name,q in
                (("pressure",pressure),("k",scalar),("convection",2*ut*dudl))}

    def integrate(self):
        ref, data=self.ref,self.data
        lumped={key: po.SparseMass(self.size[key]) for key in ("u","p","k","epsilon")}
        volume=defaultdict(float)
        regions={name: self.zeros() for name in ("inner","outer")}
        bulk={name: defaultdict(float) for name in regions}
        physical, raw, errors=defaultdict(float),defaultdict(float),defaultdict(float)
        minima=dict(k=math.inf,epsilon=math.inf,covariance_over_k=math.inf)
        pressure_mean=0.
        self.metrics_min = math.inf
        for e in self.mesh:
            self.clock()
            cutoff=(lambda x:ref.cutoff(self.exact.K(x,self.time))) if self.exact else None
            exact_crossings=[]
            if cutoff:
                xl,xr,rl,rr=e["bounds"]
                left,right=cutoff(np.array([xl,xr]))
                for face in (rl,rr):
                    if min(left,right)<face<max(left,right):
                        exact_crossings.append(brentq(lambda a:cutoff(np.array([a]))[0]-face,xl,xr))
            a,r,w=self.rule(e,extra_cutoff=cutoff,extra_x_splits=exact_crossings)
            metric=self.hill.metrics(e,a,r)
            weight=w*metric["det"]
            self.metrics_min=min(self.metrics_min,float(np.min(metric["det"])))
            B,U,G,Pi,V,P,q,D,grad,cov=self.sample(e,a,r,metric)
            outer=bool(e["epsilon"])
            region="outer" if outer else "inner"
            action, integral=regions[region],bulk[region]
            volume[region]+=float(weight.sum())
            minima["k"]=min(minima["k"],float(q["k"].min()))
            minima["epsilon"]=min(minima["epsilon"],float(q["epsilon"].min()))
            minima["covariance_over_k"]=min(minima["covariance_over_k"],float(cov))
            all_sampled=self.sampled(e,a,r,metric) if not outer else {}
            awf=sum(all_sampled.values(),np.zeros_like(U))
            if self.exact:
                f,source=self.exact.sources(e,a,r,self.time,steady=not data["step"])
                pt=self.exact.to_tube(e,a,r)
                exact=self.exact.local(*pt,self.time)
                delta=U-self.exact.rotate(exact,"u")
                errors["velocity"]+=float(weight@dot(delta,delta))
                pressure_mean+=float(weight@(Pi-exact["Pi"]))
                errors["pressure"]+=float(weight@((Pi-exact["Pi"])**2))
                errors["k"]+=float(weight@(((q["k"]-exact["k"])/ref.Kin)**2))
                if outer:
                    errors["epsilon"]+=float(weight@(((q["epsilon"]-exact["epsilon"])/exact["epsilon"])**2))
            else:
                f=np.zeros_like(U)
                source={}
            ids=np.asarray(e["u"])
            Dt=self.velocity(B["u"][0],ids,self.rates)
            adv=np.einsum("qij,qj->qi",G,U)
            div=np.trace(G,axis1=1,axis2=2)
            stress=V-Pi[:,None,None]*np.eye(2)
            b,g=B["u"]
            for i,key in enumerate(("u","w")):
                action[key][ids]+=b.T@(weight*(Dt[:,i]+adv[:,i]-awf[:,i]-f[:,i]))
                action[key][ids]+=np.einsum("qij,qj->i",g,weight[:,None]*stress[:,i])
            action["p"][e["p"]]+=B["p"][0].T@(weight*div)
            molecular=.5*ref.nu*np.einsum("qij,qij->q",G+G.transpose(0,2,1),G+G.transpose(0,2,1))
            values=dict(molecular=molecular,production=P,destruction=q["epsilon"],
                        velocity_BDF_work=dot(U,Dt),manufactured_work=dot(U,f),
                        weak_divergence_work=(dot(U,U)/2+Pi-2*q["k"]/3)*div,
                        mass_divergence=div,strong_divergence_square=div**2)
            if self.temporal:
                old,older,fs=self.temporal
                uo,upp=self.velocity(b,ids,old),self.velocity(b,ids,older)
                energy_now,energy_old,energy_older=[dot(u,u)/2 for u in (U,uo,upp)]
                values["BDF_kinetic_energy_rate"]=(
                    fs[0]*(energy_now-energy_old)+fs[2]*(energy_older-energy_old))/data["dt"]
            values.update({name+"_sampled_work": dot(U,force) for name,force in all_sampled.items()})
            for key, val in values.items():
                integral[key]+=float(weight@val)
            for key, reaction in (("k",P-q["epsilon"]),
                    ("epsilon",(1.44*P-1.82*q["epsilon"])*q["epsilon"]/q["k"])):
                if not e[key]:
                    continue
                b,g=B[key]; ids=np.asarray(e[key])
                storage=b@self.rates[key][ids]
                rhs=reaction+source.get(key,0.)
                action[key][ids]+=b.T@(weight*(storage-rhs)) + np.einsum(
                    "qij,qj->i",g,weight[:,None]*(D[key][:,None]*grad[key]-q[key][:,None]*U))
                integral[key+"_storage"]+=float(weight@storage)
                integral[key+"_source"]+=float(weight@rhs)
                integral[key+"_manufactured_source"]+=float(weight@np.broadcast_to(source.get(key,0.),storage.shape))
            for key,m in lumped.items():
                if e[key]:
                    b=B[key][0]
                    scale=ref.Kin if key=="k" else ref.A/r if key=="epsilon" else 1.
                    m.add(e[key],b.T@((weight*scale**2)[:,None]*b))
            self.faces(e,action,integral,physical,raw)
        total={key: regions["inner"][key]+regions["outer"][key] for key in self.size}
        total_area=volume["inner"]+volume["outer"]
        rr=lumped["u"].finish()
        masses={k:m.finish() for k,m in lumped.items() if k!="u"}
        masses.update(u=rr,w=rr)
        residual={}
        for key,ids in self.free.items():
            A=masses[key][ids][:,ids].tocsc()
            residual[key]=float(total[key][ids] @ splu(A).solve(total[key][ids]))
        limits=po.RESIDUAL_LIMITS
        norms=dict(momentum=math.sqrt(max(0.,residual["u"]+residual["w"])/total_area),
                   continuity=math.sqrt(max(0.,residual["p"])/total_area),
                   k=math.sqrt(max(0.,residual["k"])/total_area),
                   epsilon=math.sqrt(max(0.,residual["epsilon"])/volume["outer"]))
        stage=defaultdict(float,{key: bulk["inner"][key]+bulk["outer"][key]
                                  for key in set(bulk["inner"])|set(bulk["outer"])})
        work_free, work_reaction=0.,0.
        for i,key in enumerate(("u","w")):
            U=self.s["flow"][i*self.n:(i+1)*self.n]
            work_free+=float(U[self.free[key]]@total[key][self.free[key]])
            work_reaction+=float(U[self.fixed[key]]@total[key][self.fixed[key]])
        sampled_work=sum(stage[key+"_sampled_work"] for key in ("pressure","k","convection"))
        power=(stage["velocity_BDF_work"]+physical["kinetic_flux"]-stage["weak_divergence_work"]+
               stage["molecular"]+stage["production"]-sampled_work-stage["manufactured_work"]-
               physical["prescribed_work"]-work_reaction)
        eps_trace=self.trace_error()
        intr=self.interface_checks(regions,total)
        inner=bulk["inner"]
        for key in ("k","epsilon"):
            fixed=self.fixed[key]
            stage[key+"_essential_diffusion_reaction"]=-float(self.one[key][fixed]@total[key][fixed])
        inner["k_matching_diffusion_reaction"]=intr["weak_inner_k_diffusion_reaction"]
        inner["k_inlet_diffusion_reaction"]=-float(
            self.one["k"][self.fixed["k"]]@regions["inner"]["k"][self.fixed["k"]])
        inner["k_weak_balance"]=(inner["k_storage"]+inner["k_advective_boundary"]+
            inner["k_prescribed_diffusion"]+inner["k_matching_diffusion_reaction"]+
            inner["k_inlet_diffusion_reaction"]-inner["k_source"])
        if self.exact:
            errors["pressure"]=max(0.,errors["pressure"]-pressure_mean**2/total_area)
        return dict(residuals=norms, limits=limits,
                    residual_pass=all(norms[key]<=limit for key,limit in limits.items()),
                    integration_minimum_jacobian=self.metrics_min, chart_certificate=self.certificate,
                    minima=minima, outer_volume=volume["outer"], volume=total_area,
                    weak_divergence_work=stage["weak_divergence_work"],
                    mass_integration_identity_error=stage["mass_divergence"]-physical["mass_flux"],
                    kinetic_time_consistency=stage["velocity_BDF_work"]-stage["BDF_kinetic_energy_rate"],
                    power=dict(balance=power, free_test_work=work_free, algebraic_error=power-work_free,
                               weak_boundary_reaction_work=work_reaction,
                               pressures_sampling_work=stage["pressure_sampled_work"],
                               scalar_sampling_work=stage["k_sampled_work"],
                               convection_sampling_work=stage["convection_sampled_work"], **physical),
                    bulk=dict(stage), raw_boundary=dict(raw), layer=dict(inner), interface=intr,
                    epsilon_trace=eps_trace, **({"error": {k: math.sqrt(v/(
                        volume["outer"] if k=="epsilon" else total_area)) for k,v in errors.items()}}
                        if self.exact else {}))

    def faces(self,e,action,integral,physical,raw):
        xl,xr,rl,rr=e["bounds"]
        ext=(abs(xl-self.hill.left)<1e-12,abs(xr-self.hill.right)<1e-12,
             abs(rl)<1e-12,abs(rr-self.top)<1e-12)
        matching=(False,False,abs(rl-self.ref.H)<1e-12,abs(rr-self.ref.H)<1e-12)
        for side in range(4):
            if not (ext[side] or matching[side]):
                continue
            a,r,w=self.face(e,side)
            length,n=self.hill.face(e,side,a,r)
            weight=w*length
            B,U,_,Pi,V,_,q,D,grad,_=self.sample(e,a,r,self.hill.metrics(e,a,r))
            un=dot(n,U)
            if ext[side]:
                traction=np.einsum("qij,qj->qi",V-Pi[:,None,None]*np.eye(2),n)
                T=np.zeros_like(U)
                if not self.exact and side==1:
                    T=self.ref.outlet(r)
                    if self.data["outlet_policy"]=="ZeroTransverseV1":
                        T[:,1]=0.
                    raw["outlet_minimum_normal_velocity"]=min(
                        raw.get("outlet_minimum_normal_velocity",math.inf),float(un.min()))
                elif not self.exact and side==3:
                    T[:,0]=self.ref.top_tangent
                for i,key in enumerate(("u","w")):
                    action[key][e["u"]]-=B["u"][0].T@(weight*T[:,i])
                physical["kinetic_flux"]+=float(weight@(dot(U,U)*un/2))
                physical["prescribed_work"]+=float(weight@dot(U,T))
                physical["raw_traction_work"]+=float(weight@dot(U,traction))
                physical["mass_flux"]+=float(weight@un)
            for key in grad:
                flux=un*q[key]
                if self.exact and ext[side]:
                    neumann=self.exact.normal_flux(e,a,r,self.time,n,key)
                else:
                    neumann=(self.ref.top_flux if not self.exact and key=="epsilon" and ext[side] and side==3 else 0.)
                action[key][e[key]]+=B[key][0].T@(weight*(flux+neumann))
                integral[key+"_advective_boundary"]+=float(weight@flux)
                if matching[side]:
                    integral[key+"_matching_advection"]+=float(weight@flux)
                integral[key+"_prescribed_diffusion"]+=float(weight@np.broadcast_to(neumann,weight.shape))

    def trace_error(self):
        exact_error=normal_error=0.
        for _,e in self.interfaces:
            a,r,w=self.face(e,2)
            length,_=self.hill.face(e,2,a,r)
            B=po.basis(e,"epsilon",a,r)[0]
            E=B@self.s["epsilon"][e["epsilon"]]
            exact=self.trace(e,"k",a)**1.5/(self.ref.cl*self.ref.H)
            error=E-exact
            exact_error=max(exact_error,float(np.max(np.abs(error))))
            # The two edge interior modes are the free projection tests
            # with interpolated endpoints held fixed on this trace.
            for j in po.EDGES[2][2:]:
                phi=B[:,j]
                normal_error=max(normal_error,abs(float((w*length)@(phi*error)))/
                                 (float((w*length)@np.abs(phi))*self.ref.A/self.ref.H))
        return dict(pointwise_approximation=exact_error, interior_projection_defect=normal_error)

    def along_layer_k_flux(self,column,a):
        values=0.
        nodes,weights=self.gauss
        for e in column:
            if e["epsilon"]:
                continue
            cut=self.ref.cutoff(self.trace(e,"k",np.array([a])))[0]
            cuts=self.cuts(*e["bounds"][2:],(cut,))
            for lo,hi in zip(cuts[:-1],cuts[1:]):
                r=lo+(hi-lo)*(1+nodes)/2
                ax=np.full_like(r,a)
                metric=self.hill.metrics(e,ax,r)
                b,g=self.basis(e,"k",ax,r,metric)
                local=self.s["k"][e["k"]]
                K=self.trace(e,"k",ax)
                k=b@local
                E=K**1.5/self.ref.cl/np.maximum(r,cut)
                D=self.ref.nu+self.ref.cd*k*k/(self.ref.sk*E)
                u=self.velocity(po.basis(e,"u",ax,r)[0],np.asarray(e["u"]),self.s)
                grad=np.einsum("qij,i->qj",g,local)
                values+=(hi-lo)/2*float(weights@dot(k[:,None]*u-D[:,None]*grad,metric["t"]))
        return values

    def interface_checks(self,regions,total):
        area=jump_traction=jump_k=face_error=trace_jump=0.
        raw_mass=destruction_error=layer_mass_identity=0.
        ids={k:set() for k in ("u","w","k")}
        fluxes=defaultdict(float)
        lower_K,upper_K=math.inf,-math.inf
        crossings=[]
        layer_k_square,layer_k_measure=0.,0.
        for lo,hi in self.interfaces:
            a,r,w=self.face(lo,3)
            Li,n=self.hill.face(lo,3,a,r)
            Lo,no=self.hill.face(hi,2,a,r)
            weight=w*Li
            area+=float(weight.sum())
            face_error=max(face_error,float(np.max(np.abs(Li-Lo))),float(np.max(np.abs(n+no))))
            samples=[self.sample(e,a,r,self.hill.metrics(e,a,r)) for e in (lo,hi)]
            _,u0,_,p0,v0,_,q0,d0,g0,_=samples[0]
            _,u1,_,p1,v1,_,q1,d1,g1,_=samples[1]
            traction=np.einsum("qij,qj->qi",v1-v0-(p1-p0)[:,None,None]*np.eye(2),n)
            jk=dot(d0["k"][:,None]*g0["k"]-d1["k"][:,None]*g1["k"],n)
            trace_jump=max(trace_jump,float(np.max(np.abs(u1-u0))),
                           float(np.max(np.abs(q1["k"]-q0["k"]))),float(np.max(np.abs(p1-p0))))
            jump_traction+=float(weight@dot(traction,traction))
            jump_k+=float(weight@(jk*jk))
            for side,(_,U,_,Pi,V,_,q,D,g,_) in zip(("inner","outer"),samples):
                un=dot(U,n)
                stress=np.einsum("qij,qj->qi",V-Pi[:,None,None]*np.eye(2),n)
                momentum=un[:,None]*U-stress
                fluxes[side+"_mass"]+=float(weight@un)
                fluxes[side+"_axial_momentum"]+=float(weight@momentum[:,0])
                fluxes[side+"_vertical_momentum"]+=float(weight@momentum[:,1])
                fluxes[side+"_k"]+=float(weight@(un*q["k"]-dot(D["k"][:,None]*g["k"],n)))
                if side=="outer":
                    fluxes["outer_epsilon_dirichlet_owner_flux"]+=float(weight@dot(
                        D["epsilon"][:,None]*g["epsilon"],n))
            for key in ids:
                field="u" if key=="w" else key
                ids[key].update(hi[field][i] for i in po.EDGES[2])
            JH=self.hill.metrics(lo,a,r)
            through=JH["matching_area"]*dot(u0,n)
            K=self.trace(hi,"k",a)
            cutoff=self.ref.cutoff(K)
            destroyed=np.zeros(len(a))
            along=np.zeros(len(a))
            divint=np.zeros(len(a))
            storage_k,source_k,flux_k=np.zeros((3,len(a)))
            column=self.columns[tuple(lo["bounds"][:2])]
            xl,xr=lo["bounds"][:2]
            xfit=np.linspace(xl,xr,4)
            polynomial=np.polynomial.Chebyshev.fit(xfit,self.trace(hi,"k",xfit),3,domain=(xl,xr))
            extrema=[xl,xr]+[float(z.real) for z in polynomial.deriv().roots()
                             if abs(z.imag)<1e-12 and xl<z.real<xr]
            bounds=polynomial(np.asarray(extrema))
            lower_K=min(lower_K,float(bounds.min()))
            upper_K=max(upper_K,float(bounds.max()))
            for e in column:
                crossings.extend(self.cutoff_x_splits(e))
            for i,(aa,c) in enumerate(zip(a,cutoff)):
                for e in column:
                    if e["epsilon"]:
                        continue
                    cuts=self.cuts(*e["bounds"][2:],(c,))
                    for l,h in zip(cuts[:-1],cuts[1:]):
                        nodes,gauss=self.gauss
                        rr=l+(h-l)*(1+nodes)/2
                        ax=np.full_like(rr,aa)
                        metric=self.hill.metrics(e,ax,rr)
                        basis,U,G,_,_,production,q,_,_,_=self.sample(e,ax,rr,metric)
                        wr=gauss*(h-l)/2
                        destroyed[i]+=float(wr@(metric["det"]*q["epsilon"]))
                        distance_curve=metric["curvature"]*metric["effective_area"]
                        ut_a=metric["det"]*np.einsum("qi,qij,qj->q",metric["t"],G,metric["t"])
                        ut_a-=distance_curve*dot(U,metric["n"])
                        along[i]+=float(wr@ut_a)
                        divint[i]+=float(wr@(metric["det"]*np.trace(G,axis1=1,axis2=2)))
                        storage_k[i]+=float(wr@(metric["det"]*(basis["k"][0]@self.rates["k"][e["k"]])))
                        forcing=self.exact.sources(e,ax,rr,self.time,steady=not self.data["step"])[1]["k"] if self.exact else 0.
                        source_k[i]+=float(wr@(metric["det"]*(production-q["epsilon"]+forcing)))
                step=min(1e-4*(xr-xl),.2*(aa-xl),.2*(xr-aa))
                fa=[self.along_layer_k_flux(column,aa+step*s) for s in (-2,-1,1,2)]
                flux_k[i]=(fa[0]-8*fa[1]+8*fa[2]-fa[3])/(12*step)
            target=JH["effective_area"]*K**1.5/self.ref.cl*(
                1+np.log(self.ref.H/cutoff)+JH["curvature"]*(self.ref.H-cutoff/2))
            destruction_error=max(destruction_error,float(np.max(np.abs(destroyed-target))))
            raw_mass=max(raw_mass,float(np.max(np.abs(through+along))))
            layer_mass_identity=max(layer_mass_identity,float(np.max(np.abs(through+along-divint))))
            exchange=JH["matching_area"]*(dot(u0,n)*q0["k"]-dot(d0["k"][:,None]*g0["k"],n))
            ds=JH["effective_area"]
            residual=(storage_k+flux_k+exchange-source_k)/ds
            layer_k_square+=float((w*ds)@(residual*residual))
            layer_k_measure+=float((w*ds).sum())
        free={key: np.intersect1d(list(i),self.free[key]) for key,i in ids.items()}
        weak={key: float(np.max(np.abs(total[key][i]),initial=0.)) for key,i in free.items()}
        flux_reaction=-float(self.one["k"][free["k"]]@regions["inner"]["k"][free["k"]])
        return dict(shared_trace_max=trace_jump, physical_face_measure_max_difference=face_error,
                    measure=area, raw_traction_jump_rms=math.sqrt(jump_traction/area),
                    raw_k_diffusion_jump_rms=math.sqrt(jump_k/area), weak_free_interface_rows_max=weak,
                    weak_inner_k_diffusion_reaction=flux_reaction, raw_fluxes=dict(fluxes),
                    K_over_Kin=[lower_K/self.ref.Kin,upper_K/self.ref.Kin],
                    cutoff_split_abscissae=sorted(set(crossings)),
                    destruction_formula_max_error=destruction_error,
                    raw_integrated_normal_mass_defect=raw_mass,
                    raw_column_k_budget_rms=math.sqrt(layer_k_square/layer_k_measure),
                    layer_mass_coordinate_identity_max_error=layer_mass_identity)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("native","table","out"):
        parser.add_argument("--"+key,type=Path,required=True)
    parser.add_argument("--order",type=int,default=20)
    parser.add_argument("--seconds",type=float,default=180.)
    args=parser.parse_args()
    ref=po.load_reference(args.table)
    data,_=po.load_native(args.native)
    began=time.monotonic()
    evidence=Audit(data,ref,args.order,began+args.seconds).integrate()
    output=dict(native=identity(args.native),table=identity(args.table),
                audit_source=identity(__file__),source_identity=data["source_identity"],
                order=args.order,time=data["time"],step=data["step"],
                elapsed=time.monotonic()-began,**evidence)
    with args.out.open("x") as stream:
        json.dump(output,stream,indent=2,allow_nan=False)
        stream.write("\n")
    print(json.dumps(dict(residuals=evidence["residuals"],passed=evidence["residual_pass"],
                          elapsed=output["elapsed"])))


if __name__=="__main__":
    main()
