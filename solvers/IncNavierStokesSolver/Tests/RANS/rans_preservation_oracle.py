#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent analytical reference and weak balances in exported native bases.

Never calls a production closure or subtracts a native residual. Native inputs
are only coordinates, integration weights and test-function values/gradients.
"""
import math
from functools import lru_cache
from rans_qualification import dependencies,compare_samples
dependencies('numpy','scipy')
import numpy as np
from scipy.integrate import quad
from scipy.special import spence
from numpy.polynomial.legendre import leggauss
from validate_a93_a95_e94e95 import constants,oracle
from a93_a95_e94e95_note_wall_oracle import Kernel

class Reference:
    def __init__(self,model,policy,wall):
        self.model,self.policy,self.wall=model,policy,wall
        self.q,self.kappa,self.sigma=constants(model,policy)
        self.us=.4;self.nu=1.5e-5;self.H=2.;self.Z=4.;self.K=self.q*self.us**2
        self.A=self.us**3/self.kappa;self.t0=self.us**2
        self.kernel=Kernel(self.K,self.nu,self.H,.1,.2295 if wall=='EP' else self.kappa/math.sqrt(self.q),
                           2.55 if wall=='EP' else self.kappa*self.q**1.5)
        k=self.kernel;self.Lambda=k.Pint(self.H,self.t0)/self.A
        if self.Lambda<1:raise ValueError('independent reference cutoff inadmissible')
        self.zc=self.H*math.exp(1-self.Lambda)
        self.UH=self.t0*k.J0(self.H)
        G=np.zeros((3,3));G[0,2]=self.us/(self.kappa*self.H)
        self.pressure_jump=float(oracle(model,self.K,self.A/self.H,G)['Q'][2,2])
        self.pi=(self.pressure_jump*(self.Z-self.H)/self.Z,-self.pressure_jump*self.H/self.Z)
    def Q(self,z):
        sink=self.A*(z/self.zc if z<=self.zc else 1+math.log(z/self.zc))
        return sink-self.kernel.Pint(z,self.t0)
    def primitive(self,z):
        k=self.kernel;d=k.nu-k.a*k.zv;x=k.a*z/d;c=k.a*self.zc/d;beta=k.nu/d
        def low(v):
            l=math.log1p(v)
            return (v-l)/c-l*l/2+beta*(l+1/(1+v)-1)
        def high(v):
            l=math.log1p(v)
            return (1+math.log(v/c)+beta)*l+float(spence(1+v))-l*l/2+beta/(1+v)
        return low(x) if x<=c else low(c)+high(x)-high(c)
    @lru_cache(maxsize=None)
    def inner_k(self,z):
        if self.wall=='EP':return self.K
        if self.kernel.zv<=0:
            return self.K+self.A/self.kernel.a*(self.primitive(z)-self.primitive(self.H))
        points=sorted({z,self.H,*[p for p in self.kernel.breaks([self.zc]) if z<p<self.H]})
        return self.K-sum(quad(lambda s:self.Q(s)/self.kernel.D(s),a,b,epsabs=2e-13,epsrel=2e-13)[0] for a,b in zip(points[:-1],points[1:]))
    @lru_cache(maxsize=None)
    def at(self,z,region):
        G=np.zeros((3,3));gradient=np.zeros((5,3));values=np.zeros(5)
        if region==0:
            U=self.t0*self.kernel.J0(z);g=self.t0/self.kernel.D(z)
            k=self.inner_k(z);kp=0. if self.wall=='EP' else self.Q(z)/self.kernel.D(z)
            e=self.A/max(z,self.zc);nut=self.kernel.nt(z);D=self.kernel.D(z)
            P=nut*g*g;destroy=P if self.wall=='EP' else e
            G[0,2]=g;R=2*k*np.eye(3)/3-nut*(G+G.T)
            stress=D*(G+G.T)-self.pi[0]*np.eye(3)
            ksource=destroy-P;esource=0.;De=0.;ep=0.
        else:
            U=self.UH+self.us/self.kappa*math.log(z/self.H);g=self.us/(self.kappa*z)
            k=self.K;kp=0.;e=self.A/z;ep=-self.A/z**2;G[0,2]=g
            c=oracle(self.model,k,e,G);R=c['R'];nut=c['nut'];D=nut
            stress=-R+(2*k/3-self.pi[1])*np.eye(3)
            ksource=-c['k_source'];esource=-c['eps_source'];De=nut/self.sigma
        values[:]=[U,0.,k,e if region else 0.,self.pi[region]]
        gradient[0,2]=g;gradient[2,2]=kp;gradient[3,2]=ep if region else 0.
        return dict(values=values,gradient=gradient,G=G,R=R,stress=stress,D=D,De=De,
                    ksource=ksource,esource=esource,dependent_epsilon=e)
    def knots(self):
        base=[0.,.0025,.005,.01,.02,.04,.1,.25,.5,1.,2.,2.25,2.5,3.,3.5,4.]
        return sorted(set(base+[v for v in (self.zc,self.kernel.zv,self.kernel.ze,.1) if 0<v<self.H]))
    def self_check(self):
        # Independently integrate the defining k ODE instead of its dilogarithm.
        errors=[]
        if self.wall=='Accepted':
            for z in (0.,self.zc/2,self.zc,.7,self.H):
                points=sorted(set([z,self.H]+[x for x in self.kernel.breaks([self.zc]) if z<x<self.H]))
                answer=self.K-sum(quad(lambda s:self.Q(s)/self.kernel.D(s),a,b,epsabs=2e-13,epsrel=2e-13)[0] for a,b in zip(points[:-1],points[1:]))
                errors.append(abs(answer-self.inner_k(z)))
        return dict(k_closed_vs_independent_quadrature=max(errors,default=0.),Q0=self.Q(0.),QH=self.Q(self.H),
            k_min=min(self.inner_k(float(z)) for z in np.linspace(0,self.H,101)),Lambda=self.Lambda,cutoff=self.zc,
            pressure_jump=self.pressure_jump,outer_molecular_diffusion=0.,inner_molecular_diffusion=self.nu)

def physical_grid(knots,order=10,xorder=8):
    x,w=leggauss(order);xx,wx=leggauss(xorder);points=[];weights=[]
    for a,b in zip(knots[:-1],knots[1:]):
        for xi,wi in zip(xx,wx):
            for z,v in zip((a+b)/2+(b-a)*x/2,w*(b-a)/2):
                points.append([float(xi),0.,float(z)]);weights.append(float(wi*v))
    return points,weights

def representation(ref,packet,weights):
    grid=packet['grid'];regions=packet['grid_regions']
    target=dict(time=packet['time'],grid=grid,grid_regions=regions,grid_sides=packet.get('grid_sides',[0]*len(grid)),grid_weights=weights,grid_values=[ref.at(float(p[2]),r)['values'].tolist() for p,r in zip(grid,regions,strict=True)])
    if "grid_weights" in packet and not np.array_equal(np.asarray(packet["grid_weights"]),np.asarray(weights)):raise ValueError("reference comparison weights differ")
    numerical_packet=dict(packet,grid_weights=weights)
    result=compare_samples(numerical_packet,target,weights)
    numerical=np.asarray(packet['grid_gradients']);exact=np.asarray([ref.at(float(p[2]),r)['gradient'] for p,r in zip(grid,regions,strict=True)])
    w=np.asarray(weights);result['gradient_weighted_RMS']={}
    for j,name in enumerate(('u','v','k','epsilon_outer','Pi')):
        mask=np.asarray(regions)==1 if j==3 else np.ones(len(w),bool)
        result['gradient_weighted_RMS'][name]=float(np.sqrt(np.dot(w[mask],np.sum((numerical[mask,j,:]-exact[mask,j,:])**2,axis=1))/w[mask].sum()))
    return result

def weak_balance(ref,packet):
    """Complete weak equations, including exterior fluxes and one shared row.

    All rows (even constrained boundary modes) are tested with analytical natural
    fluxes. Interface faces enter opposite region balances and cancel in the
    assembled common row. No independent inner epsilon equation is introduced.
    """
    n=packet['sizes'][0];ne=packet['sizes'][3];np0=packet['sizes'][4];np1=packet['sizes'][5]
    accum={k:np.zeros(n) for k in ('momentum_x','momentum_z','inner_k','outer_k','inner_source','inner_diffusion')}
    scale={k:np.zeros_like(v) for k,v in accum.items()}
    accum['outer_epsilon']=np.zeros(ne);scale['outer_epsilon']=np.zeros(ne)
    accum['continuity']=np.zeros(np0+np1);scale['continuity']=np.zeros(np0+np1)
    def add(key,row,v):
        ids=np.asarray(row['ids'],int);accum[key][ids]+=v;scale[key][ids]+=np.abs(v)
    for p in packet['volume']:
        region=p['region'];z=float(p['coordinate'][2]);a=ref.at(z,region);w=p['weight'];row=p['velocity_k']
        phi=np.asarray(row['phi']);grad=np.asarray(row['derivative']);U=a['values'][0];k=a['values'][2]
        transport=a['stress'].copy();transport[0,0]-=U*U
        for component,key in ((0,'momentum_x'),(2,'momentum_z')):add(key,row,w*(transport[component,:]@grad))
        dif=w*a['D']*a['gradient'][2,2]*grad[2];src=w*a['ksource']*phi;adv=-w*U*k*grad[0]
        add('inner_k' if region==0 else 'outer_k',row,dif+src+adv)
        if region==0:
            add('inner_source',row,src);add('inner_diffusion',row,dif)
        else:
            e=p['epsilon'];pe=np.asarray(e['phi']);ge=np.asarray(e['derivative'])
            add('outer_epsilon',e,w*(a['De']*a['gradient'][3,2]*ge[2]+a['esource']*pe-U*a['values'][3]*ge[0]))
        # Continuity is evaluated in its actual divergence form with exact G.
        pressure=p['pressure'];pr=dict(pressure);pr['ids']=[i+(np0 if region else 0) for i in pressure['ids']]
        add('continuity',pr,w*np.trace(a['G'])*np.asarray(pressure['phi']))
    interface_momentum=[];interface_flux=[]
    for p in packet['faces']:
        r=p['region'];a=ref.at(float(p['coordinate'][2]),r);normal=np.asarray(p['area']);w=p['weight'];row=p['velocity_k'];phi=np.asarray(row['phi']);U=a['values'][0]
        transport=a['stress'].copy();transport[0,0]-=U*U
        flux=transport@normal
        for component,key in ((0,'momentum_x'),(2,'momentum_z')):add(key,row,-w*flux[component]*phi)
        kflux=a['D']*a['gradient'][2,2]*normal[2]-U*a['values'][2]*normal[0]
        add('inner_k' if r==0 else 'outer_k',row,-w*kflux*phi)
        if r:
            e=p['epsilon'];eflux=a['De']*a['gradient'][3,2]*normal[2]-U*a['values'][3]*normal[0]
            add('outer_epsilon',e,-w*eflux*np.asarray(e['phi']))
        if p['interface']:
            interface_momentum.append((r,flux.tolist()));interface_flux.append((r,kflux))
    accum['global_k']=accum['inner_k']+accum['outer_k'];scale['global_k']=scale['inner_k']+scale['outer_k']
    # A global row supported on both regions is counted once here.
    support=[set(),set()]
    for p in packet['volume']:support[p['region']].update(p['velocity_k']['ids'])
    shared=sorted(support[0]&support[1]);accum['shared_interface_k']=accum['global_k'][shared];scale['shared_interface_k']=scale['global_k'][shared]
    rows={}
    for key,v in accum.items():
        s=scale[key];maximum=float(np.max(np.abs(v),initial=0.));den=max(float(np.max(s,initial=0.)),1e-12)
        rows[key]=dict(max_absolute=maximum,L2_rows=float(np.linalg.norm(v)),scale=den,scaled=maximum/den,rows=len(v))
    gauge=sum(p['weight']*ref.pi[p['region']] for p in packet['pressure_mass'])
    ai=ref.at(ref.H,0);ao=ref.at(ref.H,1)
    traction=(ai['stress']-ao['stress'])@np.array([0.,0.,1.])
    traces=dict(full_traction_jump=traction.tolist(),k_trace=ai['values'][2]-ao['values'][2],
                k_flux_jump=ai['D']*ai['gradient'][2,2]-ao['D']*ao['gradient'][2,2],
                dependent_epsilon=ai['dependent_epsilon']-ao['values'][3],common_velocity=ai['values'][0]-ao['values'][0],global_pressure_integral=gauge)
    required=('momentum_x','momentum_z','inner_k','outer_k','outer_epsilon','global_k','shared_interface_k','continuity')
    checks={k:rows[k]['scaled']<1e-9 or rows[k]['max_absolute']<1e-12 for k in required}
    checks['interface_and_gauge']=max(abs(float(x)) for v in traces.values() for x in (v if isinstance(v,list) else [v]))<1e-11
    checks['accepted_nonconstant_moments_present']=ref.wall!='Accepted' or rows['inner_source']['max_absolute']>1e-6
    checks['EP_plateau_pointwise_source_zero']=ref.wall!='EP' or rows['inner_source']['max_absolute']==0.
    check=ref.self_check();checks['independent_closed_reference']=check['k_closed_vs_independent_quadrature']<1e-11
    return dict(status='pass' if all(checks.values()) else 'fail',checks=checks,rows=rows,traces=traces,
                independent_reference_check=check,shared_rows=shared,
                definition='analytical fields/gradients at native quadrature, native tests; complete boundary fluxes; no production residual used')
