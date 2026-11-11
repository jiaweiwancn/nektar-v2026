#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Measured rough-case completion and independent finite-nu flat calibration.

The data plane is the hill-absent virtual crest. Spanwise variance, epsilon,
free-stream completion, and fetch identification are explicitly assumed.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from scipy.interpolate import PchipInterpolator, CubicHermiteSpline
from scipy.integrate import quad
from scipy.optimize import minimize_scalar

QE=81/(2*math.sqrt(143)); CD=572/6561
KAPPA=math.sqrt(1.3*(1.82-1.44)*math.sqrt(CD))
NU=.051*.005/3.97; US=.051; D0=.0875; TOP=20.; DELTA=6.25
DATA=Path(__file__).with_name('cao-tamura-fig3')/'rough-fig3.csv'


def wall(K,H,h):
    a=KAPPA*math.sqrt(K/QE); A=(K/QE)**1.5/KAPPA
    x=h*math.sqrt(K)/(70*NU)
    b=10.7 if h==0 else -10.7*math.expm1(max(.5-.4*x**.7,1-.79*x**(-.28))*math.log(x))
    zv=NU*b/math.sqrt(K); start=min(H,max(0.,zv)); d=NU+a*max(start-zv,0.)
    def J(r):
        r=np.asarray(r); s=np.minimum(r,start)
        return s/NU+np.log1p(a*np.maximum(0.,r-start)/d)/a
    j0=float(J(H)); j2=start/NU**2+(H-start)/(d*(d+a*(H-start)))
    budget=(K/QE)**2*(j0-NU*j2); lam=budget/A
    rc=H*math.exp(1-lam)
    def Q(r):
        s=min(r,start); dd=NU+a*max(start-zv,0)
        j2r=s/NU**2+max(0.,r-start)/(dd*(dd+a*max(0.,r-start)))
        E=A*(r/rc if r<rc else 1+math.log(r/rc))
        return E-(K/QE)**2*(float(J(r))-NU*j2r)
    def kref(r):
        cuts=sorted(set([float(r),H]+[v for v in (rc,h,start) if r<v<H]))
        return K-sum(quad(lambda z:Q(z)/(NU+a*max(z-zv,0)),l,u,
            epsabs=2e-13,epsrel=2e-12)[0] for l,u in zip(cuts[:-1],cuts[1:]))
    return dict(a=a,A=A,zv=zv,J=J,J0=j0,Lambda=lam,rc=rc,Q=Q,kref=kref,
                molecular_fraction_H=NU/(NU+a*max(H-zv,0)))


def smooth_tail(x,y,target,end=12.5):
    """Keep digitized knots, C1 Hermite tail to a constant at declared height."""
    p=PchipInterpolator(x,y,extrapolate=False)
    tail=CubicHermiteSpline([x[-1],end],[y[-1],target],[float(p.derivative()(x[-1])),0.])
    def value(z,derivative=0):
        a=np.asarray(z,dtype=float); inside=np.clip(a,x[0],x[-1]); outer=np.clip(a,x[-1],end)
        v=p(inside,derivative); t=tail(outer,derivative)
        result=np.where(a<x[0],y[0] if derivative==0 else 0.,
            np.where(a<=x[-1],v,np.where(a<end,t,target if derivative==0 else 0.)))
        return float(result) if result.ndim==0 else result
    return value


class Completion:
    def __init__(self,span=.5,length_fraction=.09,H=None,h=None):
        if span not in (0.,.5,1.) or length_fraction not in (.05,.09,.15):raise ValueError('declared missing-input family')
        self.span=span;self.length_fraction=length_fraction
        rows=list(csv.DictReader(DATA.open()));self.rows=rows
        def panel(c):
            r=[v for v in rows if v['panel']==c]
            return np.array([float(v['y_over_H']) for v in r]),np.array([float(v['value']) for v in r])
        x,U=panel('a');self.Udata=(x,U)
        self.mean=smooth_tail(x,U,1.)
        x,su=panel('c');self.su=smooth_tail(x,su*US,.003)
        x,sw=panel('d');self.sw=smooth_tail(x,sw*US,.003)
        candidates=[v for v in (.375,.5,.75) if np.count_nonzero(self.Udata[0]-D0<=v)>=3]
        self.H=min(candidates) if H is None else H
        if h is None:
            inside=self.Udata[0]-D0<=self.H
            r=self.Udata[0][inside]-D0;U=self.Udata[1][inside];K=self.k(D0+self.H)
            def objective(h):return float(np.sum((US**2*wall(K,self.H,h)['J'](r)-U)**2))
            opt=minimize_scalar(objective,bounds=(1e-7,self.H*(1-1e-9)),method='bounded',options={'xatol':1e-14})
            if not opt.success:raise RuntimeError('rough-flat calibration failed')
            self.h=float(opt.x)
        else:self.h=h
        self.K=self.k(D0+self.H);self.w=wall(self.K,self.H,self.h);self.UH=self.mean(D0+self.H)

    def k(self,z):
        u=self.su(z);w=self.sw(z);span=(1-self.span)*w*w+self.span*u*u
        return .5*(u*u+w*w+span)

    def epsilon(self,z):
        r=np.asarray(z)-D0;L=self.length_fraction*DELTA
        mixing=KAPPA*r*(1+KAPPA*self.H/L)/(1+KAPPA*r/L)
        return CD**.75*self.k(z)**1.5/mixing

    def U(self,z,derivative=0):
        r=np.asarray(z)-D0
        if derivative==0:inner=self.UH*self.w['J'](np.minimum(np.maximum(r,0),self.H))/self.w['J0']
        else:inner=self.UH/self.w['J0']/(NU+self.w['a']*np.maximum(r-self.w['zv'],0))
        return np.where(r<self.H,inner,self.mean(z,derivative))

    def summary(self):
        x,U=self.Udata;inside=x-D0<=self.H
        fit=US**2*self.w['J'](x[inside]-D0)
        ref=[self.w['kref'](float(r)) for r in np.linspace(0,self.H,33)]
        return dict(schema='WJH2FiniteDepthCompletionV1',span_weight=self.span,
            asymptotic_length_fraction_delta=self.length_fraction,H_w=self.H,h_s=self.h,
            data_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),K_H=self.K,
            K_H_over_utau2=self.K/US**2,epsilon_H=float(self.epsilon(D0+self.H)),
            epsilon_H_over_neutral_tail=float(self.epsilon(D0+self.H)/(US**3/KAPPA/self.H)),
            U_H=self.UH,implied_initial_traction=self.UH/self.w['J0'],reported_traction=US**2,
            fit=dict(y=x[inside].tolist(),measured=U[inside].tolist(),model=fit.tolist(),
                residual=(fit-U[inside]).tolist(),rms=float(np.sqrt(np.mean((fit-U[inside])**2)))),
            overlap_checks=[dict(y=float(z),measured=float(u),AWF=float(US**2*self.w['J'](min(z-D0,self.H))))
                for z,u in zip(x,U) if z-D0<=self.H],
            cutoff=dict(Lambda=self.w['Lambda'],rc=self.w['rc'],zv=self.w['zv'],
                reference_k_min=min(ref),reference_k_ground=ref[0],
                molecular_fraction_H=self.w['molecular_fraction_H'],Q_endpoints=[self.w['Q'](0),self.w['Q'](self.H)]),
            assumptions=dict(data_plane_x=0.,inlet_x=0.,outlet_x=8.,upstream_fetch_approximation=True,
                free_stream_TI_per_component=.003,completion_height=12.5,
                epsilon_model='Blackadar finite mixing length, normalized at H to A(K)/H',
                prescribed_wall_function_initial_velocity='scaled to measured U(H), not stress-equilibrium restoration'))


def constitutive(k,e,G,Nguess=None):
    """Independent note algebra; complex continuation differentiates physical root."""
    B=k/e*(G-np.trace(G)*np.eye(3)/3); S=(B+B.T)/2;O=(B-B.T)/2;I=np.eye(3)
    s=np.trace(S@S);w=np.trace(O@O);i4=np.trace(S@O@O)
    N=max(z.real for z in np.roots([1,-1.8,-2.7*s-2*w,3.6*w]) if abs(z.imag)<1e-9) if Nguess is None else Nguess
    for _ in range(12):
        f=(N-1.8)*(N*N-2*w)-2.7*s*N
        N-=f/(3*N*N-3.6*N-2*w-2.7*s)
    Q=5/6*(N*N-2*w)*(2*N*N-w)
    b1=-N*(2*N*N-7*w)/Q;b3=-12*i4/N/Q;b4=-2*(N*N-2*w)/Q;b6=-6*N/Q;b9=6/Q
    a=b1*S+b3*(O@O-w*I/3)+b4*(S@O-O@S)+b6*(S@O@O+O@O@S-2*i4*I/3)+b9*(O@S@O@O-O@O@S@O)
    return k*(a+2*I/3),3*N/(5*(N*N-2*w))*k*k/e,N


def audit_inlet(c,probe=None,count=401):
    angles=np.linspace(0,math.pi,721);xi=np.c_[np.cos(angles),np.zeros(len(angles)),np.sin(angles)]
    ap=np.c_[-np.sin(angles),np.zeros(len(angles)),np.cos(angles)];at=np.tile([0,1,0],(len(angles),1))
    rank=[np.einsum('ni,nj->nij',a,xi).reshape(-1,9) for a in (ap,at)]
    minima=dict(covariance=(float('inf'),None),principal_in_plane=(float('inf'),None),principal_transverse=(float('inf'),None))
    rows=[]
    for z in sorted(set(np.linspace(D0+c.H,TOP,count).tolist()+[float(v['y_over_H']) for v in c.rows if float(v['y_over_H'])>=D0+c.H])):
        k=float(c.k(z));e=float(c.epsilon(z));G=np.zeros((3,3));G[0,2]=c.mean(z,1)
        R,nut,N=constitutive(k,e,G);eig=np.linalg.eigvalsh(R)
        L=np.empty((9,9))
        for j,v in enumerate(np.eye(9)):
            L[:,j]=-constitutive(k,e,G.astype(complex)+1e-28j*v.reshape(3,3),N)[0].imag.ravel()/1e-28/nut
        values=[np.einsum('ni,ij,nj->n',h,L,h) for h in rank]
        record=dict(z=z,k=k,epsilon=e,G=G.tolist(),q=k/e*G[0,2],root=float(N),
            covariance=eig.tolist(),principal=[float(v.min()) for v in values],
            angle=[float(angles[v.argmin()]) for v in values],nu_eff=float(nut))
        for name,val in zip(minima,(float(eig.min()),float(values[0].min()),float(values[1].min()))):
            if val<minima[name][0]:minima[name]=(val,record)
        if probe is not None:
            v=probe.closure(k,e,G)
            if np.max(abs(v['R']-R))>2e-12:raise AssertionError('compiled inlet covariance differs')
            for i,h in enumerate(rank):
                j=values[i].argmin();action=probe.closure(k,e,G,(0,0,h[j].reshape(3,3)))
                actual=-np.sum(np.array(action['R'])*h[j].reshape(3,3))/nut
                if abs(actual-values[i][j])>2e-10:raise AssertionError('compiled principal replay differs')
        rows.append(record)
    return dict(status='PASS' if all(v[0]>=0 for v in minima.values()) else 'OBSTRUCTED',
        scope='continuous prescribed outer inlet; finite 721-angle xz/transverse rank-one sampling; not represented-state/global proof',
        minima={k:dict(value=v[0],witness=v[1]) for k,v in minima.items()},rows=rows)


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--probe',type=Path)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    probe=None
    if a.probe:
        from validate_wj1 import Probe
        probe=Probe(str(a.probe))
    summaries=[]
    for span,length in ((.5,.09),(0.,.09),(1.,.09),(.5,.05),(.5,.15)):
        c=Completion(span,length);summary=c.summary();audit=audit_inlet(c,probe)
        key=f'span-{span:g}-length-{length:g}'
        (a.output/(key+'.json')).write_text(json.dumps(dict(completion=summary,audit=audit),indent=2,allow_nan=False)+'\n')
        summaries.append(dict(key=key,status=audit['status'],completion=summary,minima=audit['minima']))
    if probe:probe.close()
    (a.output/'summary.json').write_text(json.dumps(summaries,indent=2,allow_nan=False)+'\n')
    print(json.dumps([{k:v for k,v in s.items() if k!='minima'} for s in summaries],indent=2))
if __name__=='__main__':main()
