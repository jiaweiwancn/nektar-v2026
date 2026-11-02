#!/usr/bin/env python3
"""Source-pinned Linux M5 qualification: independent geometry and weak equations.

Native assembly maps/signs describe the FE space only. Independent Jacobi bases,
NumPy normal geometry, dimensional indexed CLS and dimensional AWF primitives
construct the oracle; no exported stress/production/volume weight is its input.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
THREADS={k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import scipy
from validate_cls_epawf_m1 import indexed
from validate_cls_epawf_m2 import Profiles
from validate_cls_epawf_m4 import Space, Runner as FlatRunner
from validate_cls_epawf_m5_geometry import geometry
from cls_epawf_m5_cases import mapped_session, plane_reference, configure_map
from cls_epawf_m5_boundary_oracle import boundary_projection
from cls_epawf_m3_cases import write, set_parameter, set_info

checks=0
errors={}

def require(condition,message):
    global checks
    checks+=1
    if not condition:raise RuntimeError(message)

def near(actual,expected,name,scale=1.,tolerance=1e-10):
    a,b=np.asarray(actual),np.asarray(expected)
    require(a.shape==b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)),'invalid '+name)
    value=float(np.max(np.abs(a-b),initial=0.))/max(scale,float(np.max(np.abs(b),initial=0.)))
    errors[name]=max(value,errors.get(name,0.))
    require(value<=tolerance,f'{name}: {value} > {tolerance}')
    return value

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

class MappedSpace(Space):
    def __init__(self,record,dimension,mapping):
        super().__init__(record,dimension);self.mapping=mapping
    def rows(self,e,points):
        ids,p,g=super().rows(e,points)
        inverse=geometry(self.mapping,np.atleast_2d(points))['inverse']
        return ids,p,np.einsum('qbi,qba->qai',g,inverse)
    def quadrature(self,e,order=20):
        x,w=super().quadrature(e,order)
        return x,w*geometry(self.mapping,x)['determinant']
    @property
    def volume(self):return sum(np.sum(self.quadrature(e)[1]) for e in self.elements)


def spaces(raw):return [MappedSpace(record,raw['dimension'],raw['map_parameters']) for record in raw['fields']]

def physical_rows(field,state,offset,points,geometry_id=None,side=0):
    points=np.atleast_2d(points)
    element=field.by_id[geometry_id] if geometry_id is not None else field.element(points[0],side)
    ids,p,g=field.rows(element,points);values=np.asarray(state)[offset+ids]
    return ids,p,g,p@values,np.einsum('qai,i->qa',g,values)

def bodies(points,dimension,inner):
    x,y,r=np.asarray(points).T
    if inner:
        values=[.01+.02*x]+([.015-.01*y+.02*r] if dimension==3 else [])+[-.02+.03*r]
    else:
        values=[.04*x+.01*r]+([.02*y-.01*r] if dimension==3 else [])+[-.03+.02*r]
    out=np.zeros((len(x),3));out[:,[0,2] if dimension==2 else [0,1,2]]=np.array(values).T
    return out


def forcing(raw,state,point,fields=None,include_body=True):
    d=raw['dimension'];fields=fields or spaces(raw);offset=raw['offsets']
    p=np.array(point,dtype=float);p[2]=raw['parameters'][1]
    geom=geometry(raw['map_parameters'],p);n=geom['normal'];P=np.eye(3)-np.outer(n,n)
    U=np.zeros(3);G=np.zeros((3,3))
    for b,axis in enumerate((0,2) if d==2 else (0,1,2)):
        _,_,_,value,gradient=physical_rows(fields[0],state,offset[b],[p],side=1)
        U[axis]=value[0];G[axis]=gradient[0]
    _,_,_,_,gp=physical_rows(fields[d+3],state,offset[d+3],[p],side=1)
    _,_,_,_,gk=physical_rows(fields[d],state,offset[d],[p],side=1)
    ut=P@U;w=U@n;gradut=P@G-np.outer(n,U@geom['Br'])-w*geom['Br']
    value=gradut@ut+ut*np.trace(gradut)+gp[0]-2*gk[0]/3
    if include_body:value-=bodies([p],d,False)[0]
    return P@value


def residual_oracle(raw,normal_order=24,tangential_order=None,outer_order=None,verify=True,weight_function=None,return_loads=False):
    """Complete all-field conservative weak operator, including unmasked rows."""
    f=spaces(raw);d=raw['dimension'];off=raw['offsets'];mapping=raw['map_parameters']
    nu,H,rough,rK,rm,K0,qstar,kappa=raw['parameters']
    alpha=raw.get('alpha',5.3);chi=raw.get('chi',0.);state=np.asarray(raw['lifted_state']);history=np.asarray(raw['history'])
    axes=(0,2) if d==2 else (0,1,2);result=np.zeros(len(state));parts=np.zeros((2,len(state)))
    roughness=np.zeros(5);energy=np.zeros(9)
    tangential_order=tangential_order or raw.get('quadrature_orders',[12,16])[0]
    outer_order=outer_order or raw.get('quadrature_orders',[12,16])[0]
    def add(region,block,ids,p,g,w,source,flux):
        np.add.at(parts[region],off[block]+ids,p.T@(w*source)+np.einsum('qai,qa,q->i',g,flux,w))
    def integrate(eid,points,w,inner,profile=None,force=None):
        region=0 if inner else 1;metrics=geometry(mapping,points);r=metrics['distance']
        velocities=[physical_rows(f[0],state,off[b],points,eid) for b in range(d)]
        ids,p,g=velocities[0][:3];U=np.zeros((len(points),3));G=np.zeros((len(points),3,3))
        for b,axis in enumerate(axes):U[:,axis]=velocities[b][3];G[:,axis]=velocities[b][4]
        ki,kp,kg,K,gK=physical_rows(f[d],state,off[d],points,eid)
        pb=d+2 if inner else d+3
        pi,pp,pg,pressure,_=physical_rows(f[pb],state,off[pb],points,eid)
        x,y,rho=points.T;S=G+G.transpose(0,2,1)
        body=bodies(points,d,inner)
        if inner:
            nt=float(profile.a)*np.maximum(r-float(profile.zv),0)
            reynolds=2*K[:,None,None]/3*np.eye(3)-nt[:,None,None]*S
            production=-np.sum(reynolds*G,axis=(1,2));stress=(nu+nt)[:,None,None]*S-pressure[:,None,None]*np.eye(3)
            ea=float(profile.k)**1.5/(2.55*np.maximum(r,float(profile.ze)))
            p0=np.asarray(profile.fields(r)[1],dtype=float)
            local_chi=float(profile.chi)
            sk=production-(1-local_chi)*ea-local_chi*p0+.03+.01*x+.04*y+.02*rho
            dk=nu+nt
            below=r<=rough
            body+=below[:,None]*force
            roughness[:3]+=np.sum(w[below,None]*force,axis=0)
            roughness[3]+=np.sum(w[below]*np.einsum('qi,i->q',U[below],force))
            roughness[4]=max(roughness[4],abs(force@metrics['normal'][0]))
        else:
            ei,ep,eg,E,gE=physical_rows(f[d+1],state,off[d+1],points,eid)
            closure=[indexed(a,k,e,'physical') for a,k,e in zip(G,K,E)]
            reynolds=np.asarray([c['stress'].real for c in closure]);nt=np.asarray([c['nut'].real for c in closure])
            production=-np.sum(reynolds*G,axis=(1,2));stress=-reynolds+(2*K/3-pressure)[:,None,None]*np.eye(3)
            dk,de=nt,nt/1.3;sk=production-E-.02+.01*x+.02*y+.03*rho
            ratio=K**1.5/(2.5*E*r)
            se=1.44*E/K*production-1.92*E*E/K+.83*E*E/K*np.maximum(ratio-1,0)*ratio**2+.01*(1+x+y+rho)
            add(region,d+1,ei,ep,eg,w,-se,de[:,None]*gE-U*E[:,None])
        for b,axis in enumerate(axes):
            add(region,b,ids,p,g,w,-body[:,axis],stress[:,axis]-U[:,axis,None]*U)
        add(region,d,ki,kp,kg,w,-sk,dk[:,None]*gK-U*K[:,None])
        add(region,pb,pi,pp,pg,w,np.trace(G,axis1=1,axis2=2),np.zeros_like(U))
        energy[2]+=np.sum(w*(np.einsum('qab,qab->q',G,stress)-sk))
        energy[3]-=np.sum(w*np.einsum('qab,qa,qb->q',G,U,U))
        energy[4]-=np.sum(w*np.sum(U*body,axis=1))
    # Independent normal panels use factor three, and split at every actual
    # native element/ray breakpoint. Chart geometry is not an integration input.
    zt,wt=np.polynomial.legendre.leggauss(tangential_order)
    zn,wn=np.polynomial.legendre.leggauss(normal_order)
    columns=sorted({tuple(e['lower'][:2]+e['upper'][:2]) for e in f[d+2].elements})
    for x0,y0,x1,y1 in columns:
        for index in itertools.product(range(len(zt)),repeat=d-1):
            x=x0+(zt[index[0]]+1)*(x1-x0)/2;y=0.;wh=wt[index[0]]*(x1-x0)/2
            if d==3:y=y0+(zt[index[1]]+1)*(y1-y0)/2;wh*=wt[index[1]]*(y1-y0)/2
            chart=np.array([x,y,0.]);metric=geometry(mapping,chart);scale=metric['scale'];height=metric['height']
            sampled=f[d].evaluate(state,off[d],[x,y,rK/scale])[0]
            local_chi=chi if weight_function is None else weight_function(chart,state,f)
            profile=Profiles(dict(geometry=[nu,height,rough,rK,rm],metric=[metric['kappa1'],metric['kappa2']],k=sampled,q=qstar,kappa=kappa,chi=local_chi))
            FA=forcing(raw,state,[x,y,H],f)
            breaks={0.,float(height),rough,float(profile.rc),float(profile.ze)}
            if profile.zv>0:breaks.add(float(profile.zv))
            start=min(float(profile.rc),nu/float(profile.a),float(profile.ze))
            while start<height:breaks.add(start);start*=3
            for e in f[d+2].elements:
                if e['lower'][0]<=x<=e['upper'][0] and (d==2 or e['lower'][1]<=y<=e['upper'][1]):
                    breaks.update((e['lower'][2]*scale,e['upper'][2]*scale))
            breaks=sorted(v for v in breaks if 0<=v<=height)
            for lo,hi in zip(breaks[:-1],breaks[1:]):
                r=lo+(zn+1)*(hi-lo)/2
                points=np.column_stack([np.full_like(r,x),np.full_like(r,y),r/scale])
                eid=f[d+2].element([x,y,(lo+hi)/(2*scale)])['geometry']
                w=wh*wn*(hi-lo)/2*metric['baseArea']*np.asarray(profile.jacobian(r),dtype=float)
                integrate(eid,points,w,True,profile,FA)
    for e in f[d+3].elements:
        points,w=f[d+3].quadrature(e,outer_order)
        integrate(e['geometry'],points,w,False)
    # Physical consistent mass, separate regional pressure means, one gauge.
    for b in range(d+2):
        value=alpha*f[b].mass()@(state[off[b]:off[b+1]]-history[off[b]:off[b+1]])
        result[off[b]:off[b+1]]+=value
        if b<d:energy[1]+=state[off[b]:off[b+1]]@value
        if b==d:energy[1]+=f[b].constant()@value
        if b<=d:
            current=state[off[b]:off[b+1]];test=current*.5 if b<d else f[b].constant()
            energy[8]+=test@f[b].mass()@current
    volume=f[0].volume
    for region,b in enumerate((d+2,d+3)):
        load=f[b].mass()@f[b].constant()/volume
        parts[region,off[b]:off[b+1]]+=state[off[-1]]*load
        parts[region,off[-1]]+=state[off[b]:off[b+1]]@load
    result+=parts.sum(axis=0)
    for b in range(d):energy[0]+=state[off[b]:off[b+1]]@result[off[b]:off[b+1]]
    energy[0]+=f[d].constant()@result[off[d]:off[d+1]]
    energy[6]=sum(energy[1:5]);energy[7]=energy[0]-energy[6]
    block_errors=[]
    for b in range(d+4):
        actual=np.asarray(raw['raw_residual'])[off[b]:off[b+1]];expected=result[off[b]:off[b+1]]
        block_errors.append(float(np.max(np.abs(actual-expected),initial=0.)))
        if verify:near(actual,expected,'mapped_residual_block_'+str(b),tolerance=2e-10)
    if verify:
        near(raw['raw_residual'][-1],result[-1],'mapped_pressure_gauge')
        near(raw['roughness_budget'],roughness,'roughness_force_and_work',tolerance=2e-10)
        near(raw['energy_budget'],energy,'resolved_plus_k_energy_budget',tolerance=2e-10)
        require(abs(energy[7])<1e-10*max(1.,abs(energy[0])),'independent energy identity')
        for sample in raw['forcing_samples']:
            near(sample['force'],forcing(raw,state,sample['chart'],f),'rotating_frame_forcing',tolerance=1e-10)
    receipt=dict(dimension=d,maximum_block_absolute_errors=block_errors,roughness_budget=roughness.tolist(),energy_budget=energy.tolist(),
                independent_quadrature_orders=[tangential_order,normal_order,outer_order],
                comparison='same declared finite outer Gauss rule, independently reconstructed; higher-order inner quadrature')
    if return_loads:receipt['residual']=result.tolist()
    return receipt


def dependent_epsilon(raw,state,chart,fields=None):
    fields=fields or spaces(raw);d=raw['dimension'];off=raw['offsets'];p=np.asarray(chart,dtype=float)
    m=geometry(raw['map_parameters'],p);r=m['distance'];nu,H,rough,rK,rm,*rest=raw['parameters']
    sample=np.array([p[0],p[1],rK/m['scale']]);sm=geometry(raw['map_parameters'],sample)
    _,_,_,K,gK=physical_rows(fields[d],state,off[d],[sample]);K=float(K[0])
    # Differentiate the fixed-physical-r sampler before following the H trace.
    chart_gradient=gK[0]@sm['jacobian'];hg=m['height_gradient']
    sampled_gradient=np.array([chart_gradient[0]-chart_gradient[2]*sample[2]*hg[0]/m['height'],chart_gradient[1]-chart_gradient[2]*sample[2]*hg[1]/m['height'],0.])
    ze=5.1*nu/math.sqrt(K)
    if r>ze:
        value=K**1.5/(2.55*r)
        r_gradient=np.array([p[2]*hg[0]/H,p[2]*hg[1]/H,m['scale']])
        gradient=1.5*math.sqrt(K)/(2.55*r)*sampled_gradient-K**1.5/(2.55*r*r)*r_gradient
    else:
        value=K*K/(2.55*5.1*nu);gradient=2*K/(2.55*5.1*nu)*sampled_gradient
    return value,gradient@m['inverse'],K


def interface_oracle(raw):
    d=raw['dimension'];off=raw['offsets'];state=raw['lifted_state'];f=spaces(raw);mapping=raw['map_parameters']
    nu,H,rough,rK,rm,K0,qstar,kappa=raw['parameters'];axes=(0,2) if d==2 else (0,1,2)
    def evaluate(q):
        metric=geometry(mapping,q);area=metric['matching_area'];norm=float(np.linalg.norm(area));N=area/norm
        U=np.zeros(3);Ui=np.zeros(3);Gi=np.zeros((3,3));Go=np.zeros((3,3))
        for b,axis in enumerate(axes):
            _,_,_,value,gradient=physical_rows(f[0],state,off[b],[q],side=-1);Ui[axis]=value[0];Gi[axis]=gradient[0]
            _,_,_,value,gradient=physical_rows(f[0],state,off[b],[q],side=1);U[axis]=value[0];Go[axis]=gradient[0]
        _,_,_,k,gki=physical_rows(f[d],state,off[d],[q],side=-1)
        _,_,_,_,gko=physical_rows(f[d],state,off[d],[q],side=1)
        _,_,_,e,ge=physical_rows(f[d+1],state,off[d+1],[q],side=1)
        _,_,_,pi,_=physical_rows(f[d+2],state,off[d+2],[q],side=-1)
        _,_,_,po,_=physical_rows(f[d+3],state,off[d+3],[q],side=1)
        k,e,pi,po=map(lambda a:float(a[0]),(k,e,pi,po))
        Ew,gEw,K=dependent_epsilon(raw,state,q,f)
        profile=Profiles(dict(geometry=[nu,metric['height'],rough,rK,rm],metric=[metric['kappa1'],metric['kappa2']],k=K,q=qstar,kappa=kappa,chi=0))
        nt=float(profile.a)*max(metric['height']-float(profile.zv),0.)
        closure=indexed(Go,k,e,'physical');nto=float(closure['nut'].real)
        Ti=(nu+nt)*(Gi+Gi.T)-pi*np.eye(3)
        To=-closure['stress'].real+(2*k/3-po)*np.eye(3)
        tractionI=Ti@area;tractionO=To@area;mass=U@area
        ki=(U*k-(nu+nt)*gki[0])@area;ko=(U*k-nto*gko[0])@area
        eo=(U*e-nto/1.3*ge[0])@area;ei=(U*Ew)@area
        base=metric['base_tangents'];gradH=base@np.linalg.solve(base.T@base,metric['height_gradient'])
        w=U@metric['normal'];ut=U-w*metric['normal'];v=np.linalg.solve(np.eye(3)+metric['height']*metric['B'],ut)
        m=mass/metric['baseArea'];formula=metric['J']*(w-v@gradH);naive=metric['J']*w
        near(m,formula,'variable_height_relative_matching_mass')
        record=dict(area=area,N_H=N,normal=metric['normal'],U=U,G_inner=Gi,G_outer=Go,traction_inner_area=tractionI,traction_outer_area=tractionO,
                    shared_advective_momentum=U*mass,k=k,sampled_k=K,epsilon=e,dependent_epsilon=Ew,physical_r=metric['distance'],J=metric['J'],base_area=metric['baseArea'],
                    m_H=m,m_H_formula=formula,m_H_naive=naive,shared_advective_mass=mass,shared_advective_k=k*mass,k_flux_inner=ki,k_flux_outer=ko,
                    epsilon_flux_outer=eo,epsilon_inner_advective_flux=ei,epsilon_flux_reaction=eo-ei)
        integral=np.array([norm,mass,k*mass,*(U*mass),*tractionI,*tractionO,tractionI@N,tractionO@N,ki,ko,eo,ei,eo-ei,metric['baseArea']*naive])
        return record,integral
    for sample in raw['interface']['samples']:
        result,_=evaluate(sample['chart'])
        for name,value in result.items():near(sample[name],value,'independent_interface_'+name,tolerance=1e-10)
    z,w=np.polynomial.legendre.leggauss(raw['interface']['quadrature_order']);total=np.zeros(20)
    for elem in f[d+2].elements:
        if abs(elem['upper'][2]-H)>1e-12:continue
        for index in itertools.product(range(len(z)),repeat=d-1):
            q=np.array([elem['lower'][0]+(z[index[0]]+1)*(elem['upper'][0]-elem['lower'][0])/2,0.,H])
            weight=w[index[0]]*(elem['upper'][0]-elem['lower'][0])/2
            if d==3:
                q[1]=elem['lower'][1]+(z[index[1]]+1)*(elem['upper'][1]-elem['lower'][1])/2
                weight*=w[index[1]]*(elem['upper'][1]-elem['lower'][1])/2
            total+=weight*evaluate(q)[1]
    near(raw['interface']['integrals'],total,'independent_integrated_matching_budgets',tolerance=1e-10)
    require(abs(total[18])>1e-6,'epsilon constraint reaction was silently omitted')
    if np.any(np.abs(np.asarray(mapping)[3:6])>0):require(abs(total[1]-total[19])>1e-5,'variable-H naive mass omission control inactive')
    diagnostic=[]
    for sample in raw['inner_epsilon_diagnostic']:
        q=np.asarray(sample['chart']);value,gradient,K=dependent_epsilon(raw,state,q,f)
        old=dependent_epsilon(raw,raw['history'],q,f)[0]
        U=np.zeros(3);G=np.zeros((3,3))
        for b,axis in enumerate(axes):
            _,_,_,v,g=physical_rows(f[0],state,off[b],[q]);U[axis]=v[0];G[axis]=g[0]
        reaction=raw['alpha']*(value-old)+value*np.trace(G)+U@gradient
        near(sample['epsilon'],value,'inner_algebraic_epsilon')
        near(sample['history_epsilon'],old,'inner_history_epsilon')
        near(sample['gradient'],gradient,'inner_epsilon_analytic_gradient_vs_fd',tolerance=2e-7)
        near(sample['reaction'],reaction,'inner_epsilon_diagnostic_reaction',tolerance=2e-7)
        diagnostic.append(reaction)
    require(max(map(abs,diagnostic))>1e-5,'diagnostic epsilon residual was forced to zero')
    return dict(integrals=total.tolist(),epsilon_diagnostic_reactions=diagnostic,normal_traction=[float(total[12]),float(total[13])],generic_state_is_not_traction_balanced=True)


def sum_load(field,function):
    result=np.zeros(field.size)
    for e in field.elements:
        points,w=field.quadrature(e);ids,p,g=field.rows(e,points)
        np.add.at(result,ids,p.T@(w*np.array([function(q) for q in points])))
    return result


def native_history_difference(a,b,name,tolerance=1e-8):
    """Compare every native history slot and the full coupled physical state."""
    require(a['dimension']==b['dimension'] and a['order']==b['order'],'different compared history shapes')
    near(a['map_parameters'],b['map_parameters'],name+'_geometry',tolerance=0.)
    near(a['native_history_times'],b['native_history_times'],name+'_history_times')
    near([a['time'],a['steps'],a['dt']],[b['time'],b['steps'],b['dt']],name+'_time')
    d,us,H=a['dimension'],a['ustar'],a['parameters'][1]
    K=a['parameters'][5];E=us**3/(a['parameters'][7]*H)
    scales=[us]*d+[K,E,us*us,us*us];maximum=0.
    for block,(field,scale) in enumerate(zip(spaces(a),scales)):
        start=a['offsets'][block];stop=start+field.size;mass=field.mass();volume=field.volume
        delta=np.asarray(a['state'][start:stop])-b['state'][start:stop]
        value=math.sqrt(max(0.,float(delta@mass@delta)))/(scale*math.sqrt(volume));maximum=max(maximum,value)
        require(value<=tolerance,name+' full coupled field '+str(block)+': '+str(value))
        if block<d+2:
            require(len(a['native_history'])==len(b['native_history'])==2*a['order'],'incomplete native history')
            for slot,(first,second) in enumerate(zip(a['native_history'],b['native_history'])):
                delta=np.asarray(first[start:stop])-second[start:stop]
                stored_scale=scale if slot<a['order'] else scale*a['dt']/(H/us)
                value=math.sqrt(max(0.,float(delta@mass@delta)))/(stored_scale*math.sqrt(volume));maximum=max(maximum,value)
                require(value<=tolerance,name+' native history '+str(slot)+': '+str(value))
    near(a['state'][-1],b['state'][-1],name+'_gauge',tolerance=tolerance)
    near(a['fixed_reference'],b['fixed_reference'],name+'_fixed_reference',tolerance=1e-10)
    near(a['fixed_defect'],b['fixed_defect'],name+'_fixed_defect',tolerance=1e-10)
    errors[name]=maximum
    return maximum


def stage_budgets(raw):
    d=raw['dimension'];off=raw['offsets'];loads=np.asarray(raw['final_stage_loads']);regions=np.asarray(raw['regional_stage_loads'])
    require(regions.shape==(2,len(loads)) and len(loads)==len(raw['state']),'missing actual curved stage budgets')
    near(regions.sum(axis=0),loads,'regional_conservative_stage_sum',tolerance=1e-10)
    constraints={int(row) for row,value in raw['final_stage_constraints']}
    f=spaces(raw);interface={}
    for b in list(range(d))+[d]:
        ids=f[b].boundary(lambda axis,side,e: axis==2 and abs((e['upper'][2] if side else e['lower'][2])-raw['parameters'][1])<1e-11)
        rows=sorted(off[b]+i for i in ids if off[b]+i not in constraints)
        require(rows,'missing free shared trace rows')
        near(regions[0,rows],-regions[1,rows],'shared_stage_trace_'+str(b),tolerance=1e-9)
        interface[str(b)]=dict(free_rows=len(rows),inner_reaction_max=float(np.max(np.abs(regions[0,rows]),initial=0)),sum_max=float(np.max(np.abs(loads[rows]),initial=0)))
    require(abs(raw['state'][-1])<1e-9,'pressure multiplier conceals incompatible net flux')
    near(raw['anchor_energy_work'],0.,'curved_no_anchor_work')
    require(raw['eta']==0 and raw['chi']==0,'curved stage activated equilibrium correction')
    if raw['energy_budget']:near(raw['energy_budget'][7],0.,'native_energy_budget_closure',tolerance=1e-10)
    return interface


def mms_error(raw,swirl=False,stationary=False):
    fields=spaces(raw);d=raw['dimension'];off=raw['offsets'];s=np.asarray(raw['state']);t=raw['time'];mp=raw['map_parameters']
    K0=raw['parameters'][5];K=K0 if stationary else K0*(1+math.sin(t)/50);A=K**1.5/2.55
    H0=raw['parameters'][1];norms=np.zeros(d+4);volumes=np.zeros(d+4);normal=0.;strongdiv=0.;normal_max=0.;divergence_max=0.;pressure_mean=0.
    for b in (d+2,d+3):
        for e in fields[b].elements:
            p,w=fields[b].quadrature(e,20);_,_,_,value,_=physical_rows(fields[b],s,off[b],p,e['geometry'])
            target=0*p[:,0] if int(mp[0])==4 else .002*geometry(mp,p)['point'][:,2]
            pressure_mean+=np.sum(w*(value-target))
    pressure_mean/=fields[0].volume
    for b in range(d+4):
        for e in fields[b].elements:
            p,w=fields[b].quadrature(e,20);_,_,_,value,gradient=physical_rows(fields[b],s,off[b],p,e['geometry']);g=geometry(mp,p)
            if b==0:
                U=np.zeros((len(p),3));G=np.zeros((len(p),3,3))
                for component,axis in enumerate((0,2) if d==2 else (0,1,2)):
                    _,_,_,v,dv=physical_rows(fields[0],s,off[component],p,e['geometry']);U[:,axis]=v;G[:,axis]=dv
                un=np.einsum('qa,qa->q',U,g['normal']);divergence=np.trace(G,axis1=1,axis2=2)
                normal+=np.sum(w*un**2);strongdiv+=np.sum(w*divergence**2)
                normal_max=max(normal_max,float(np.max(np.abs(un))));divergence_max=max(divergence_max,float(np.max(np.abs(divergence))))
            if b<d:
                target=np.zeros(len(p))
                if swirl and e['upper'][2]<=H0+1e-12:
                    r=g['distance'];ss=r/H0;speed=.3*(1+.1*math.sin(t))*ss**2*(1-ss)**2;angle=p[:,0]/mp[10]
                    target=speed*(np.cos(angle) if b==0 else -np.sin(angle))
            elif b==d:target=np.full(len(p),K)
            elif b==d+1:target=A/g['distance']
            else:target=(0*p[:,0] if int(mp[0])==4 else .002*g['point'][:,2])+pressure_mean
            norms[b]+=np.sum(w*(value-target)**2);volumes[b]+=np.sum(w)
    norms=np.sqrt(norms/volumes)
    return dict(field_rms=norms.tolist(),velocity_rms=float(np.linalg.norm(norms[:d])),k_rms=float(norms[d]),epsilon_rms=float(norms[d+1]),pressure_rms=float(np.linalg.norm(norms[d+2:])),pressure_gauge_offset=float(pressure_mean),normal_velocity_rms=math.sqrt(normal/volumes[0]),normal_velocity_sample_max=normal_max,strong_divergence_rms=math.sqrt(strongdiv/volumes[0]),strong_divergence_sample_max=divergence_max)

class Runner(FlatRunner):
    def curved(self,name,**kwargs):return self.evolution(name,generator=mapped_session,**kwargs)
    def wall(self,name,family,dimension=2,degree=3,split=1,ranks=1,ownership_shift=0,quadrature=12,invariant=False):
        directory=self.output/name;directory.mkdir()
        xml=mapped_session(self.template,directory/'case.xml',family=family,dimension=dimension,degree=degree,split=split,ownership_shift=ownership_shift,quadrature=quadrature)
        if invariant:
            tree=ET.parse(xml);set_parameter(tree.getroot(),'RANSWallProbeInvariant',1);write(tree.getroot(),xml)
        report=directory/'probe.json';self.command(name,'RANSCurvedWallProbe',[xml,report],ranks=ranks)
        return json.loads(report.read_text())


def qualify(run,sections):
    directory=run.output/'geometry';directory.mkdir()
    report=directory/'geometry.json'
    run.command('geometry','RANSCurvedGeometryProbe',['--output',report])
    geometry_check=run.output/'geometry-check.json'
    run.command('geometry-oracle',sys.executable,[str(Path(__file__).with_name('validate_cls_epawf_m5_geometry.py')),'--geometry-report',report,'--output',geometry_check])
    sections['geometry']=json.loads(geometry_check.read_text())
    sections['weak_operators']=[]
    for family,dimension in [('plane',2),('circle',2),('smooth',2),('variable',2),('surface3d',3),('feet',2)]:
        raw=run.wall('operator-'+family,family,dimension)
        independent=residual_oracle(raw)
        independent['interface']=interface_oracle(raw)
        independent['boundary_projection']=boundary_projection(raw,spaces(raw),geometry,near,require)
        if family=='plane':
            independent['higher_outer_order_diagnostic']=residual_oracle(raw,outer_order=16,verify=False)
            gap=[];field=spaces(raw)[0];state=raw['lifted_state'];axes=(0,2) if dimension==2 else (0,1,2)
            for e in spaces(raw)[dimension+3].elements:
                points,w=field.quadrature(field.by_id[e['geometry']],20);G=np.zeros((len(points),3,3))
                for b,axis in enumerate(axes):G[:,axis]=physical_rows(field,state,raw['offsets'][b],points,e['geometry'])[4]
                gap.extend(np.einsum('qab,qba->q',G,G))
            independent['invariant_max_switch']=dict(trace_G_squared_range=[float(min(gap)),float(max(gap))],sign_change=bool(min(gap)<0<max(gap)))
        for mass in raw['physical_masses']:
            near(mass['values'],spaces(raw)[mass['block']].mass(),'independent_physical_mass',tolerance=1e-10)
        path=run.output/('operator-'+family)/'probe.json'
        audit=run.output/('operator-'+family)/'geometry-check.json'
        run.command('wall-audit-'+family,sys.executable,[str(Path(__file__).with_name('validate_cls_epawf_m5_geometry.py')),'--geometry-report',report,'--wall-report',path,'--output',audit])
        independent.update(family=family,geometry_audit=json.loads(audit.read_text()))
        sections['weak_operators'].append(independent)
    # The same global chart polynomials lie exactly in both native spaces.
    # Change the actual lateral patch subdivision and compare identical weak
    # tests, without confusing spatial refinement with coefficient equality.
    partition_levels=[]
    for order in (12,24,48):
        first=run.wall('partition-one-q'+str(order),'circle',invariant=True,quadrature=order)
        second=run.wall('partition-two-q'+str(order),'circle',split=2,invariant=True,quadrature=order)
        first_fields,second_fields=spaces(first),spaces(second);moments=[];maximum=0.
        for block,(a,b) in enumerate(zip(first_fields,second_fields)):
            for test,name in [(lambda p:1.,'one'),(lambda p:p[0],'q1'),(lambda p:p[2],'rho')]:
                ca=np.linalg.solve(a.mass(),sum_load(a,test));cb=np.linalg.solve(b.mass(),sum_load(b,test))
                ra=np.asarray(first['raw_residual'])[first['offsets'][block]:first['offsets'][block+1]]
                rb=np.asarray(second['raw_residual'])[second['offsets'][block]:second['offsets'][block+1]]
                value=float(ca@ra);other=float(cb@rb);error=abs(value-other)/max(1.,abs(other));maximum=max(maximum,error)
                moments.append(dict(block=block,test=name,one_patch=value,two_patches=other,error=error))
                if order==48:near(value,other,'lateral_patch_partition_'+str(block)+'_'+name,tolerance=1e-10)
            for point in ([.17,0.,.37],[.61,0.,1.7],[.83,0.,2.9]):
                if not any(np.all(np.asarray(point)>=e['lower']) and np.all(np.asarray(point)<=e['upper']) for e in a.elements):continue
                va,ga=a.evaluate(first['lifted_state'],first['offsets'][block],point);vb,gb=b.evaluate(second['lifted_state'],second['offsets'][block],point)
                near([va,*ga],[vb,*gb],'partition_identical_physical_field',tolerance=1e-10)
        partition_levels.append(dict(outer_order=order,maximum_weak_moment_difference=maximum,moments=moments))
    require(partition_levels[-1]['maximum_weak_moment_difference']<1e-10,'lateral partition loads did not converge')
    sections['patch_partition']=dict(identical_global_polynomials=True,actual_lateral_patch_counts=[1,2],quadrature_refinement=partition_levels)
    # Rotation covariance of the entire projected target and actual IMEX map.
    # The plane is an isolated, constant-H neutral verification case (eta=1).
    sections['planar_reference']=[]
    for dimension,order in [(2,1),(2,2),(3,2)]:
        options=dict(degree=3,dimension=dimension,steps=2,order=order,dt=.02,eta=1,chi=1,length=1.,span=1.,inner_knots=[0.,.1,2.],outer_knots=[2.,4.],tangential_order=8)
        if dimension==3:options['direction']=(.8,.6)
        flat=run.evolution(f'flat-reference-d{dimension}-o{order}',**options)
        plane=run.evolution(f'plane-reference-d{dimension}-o{order}',generator=plane_reference,**options)
        foff,poff=flat['offsets'],plane['offsets'];axes=(0,2) if dimension==2 else (0,1,2)
        g=geometry(plane['map_parameters'],np.zeros(3));n=g['normal'];T=np.eye(3)-np.outer(n,n)
        # One tangent rotates x into the inclined plane; the extruded y axis is
        # unchanged. Compare coefficients before/after steps, including pi jump.
        t1=np.array([1.,0.,.25]);t1/=np.linalg.norm(t1)
        for key in ('initial','state'):
            for b,axis in enumerate(axes):
                target=t1[axis]*np.asarray(flat[key][foff[0]:foff[1]])
                if dimension==3 and axis==1:target=np.asarray(flat[key][foff[1]:foff[2]])
                near(plane[key][poff[b]:poff[b+1]],target,'planar_rotated_'+key+'_'+str(b),tolerance=1e-9)
            near(plane[key][poff[dimension]:],flat[key][foff[dimension]:],'planar_scalars_pressure_'+key,tolerance=1e-9)
        require(max(a[2] for a in plane['step_audits'])<1e-9,'rotated planar neutral target drifted')
        sections['planar_reference'].append(dict(dimension=dimension,order=order,maximum_drift=max(a[2] for a in plane['step_audits']),eta=1,chi=1))
    sections['manufactured']=[]
    curved_cases={}
    for family,dimension,swirl in [('circle',2,True),('variable',2,False),('surface3d',3,False),('feet',2,False)]:
        degrees=[3,4] if dimension==3 else ([3,5,6] if family=='circle' else [3,4,5])
        records=[]
        for degree in degrees:
            name=f'mms-{family}-p{degree}'
            raw=run.curved(name,family=family,dimension=dimension,degree=degree,steps=1,order=2,dt=.005,swirl=swirl,stationary=not swirl,quadrature=8 if dimension==3 else 12)
            balance=stage_budgets(raw);error=mms_error(raw,swirl=swirl,stationary=not swirl)
            records.append(dict(degree=degree,error=error,interface=balance,boundary_surface_rms_correction=raw['boundary_flux_correction']))
            curved_cases[family]=raw
        # Fixed physical H,rK,rm and domain throughout each p sequence. For the
        # C1-foot case pressure/velocity are identically zero; epsilon tests the
        # genuine piecewise geometry, including its curvature jump.
        require(records[-1]['error']['epsilon_rms']<records[0]['error']['epsilon_rms']*.6,'outer epsilon p refinement did not converge: '+family)
        require(records[-1]['error']['velocity_rms']<1e-3,'manufactured velocity error too large: '+family)
        if family!='feet':
            require(records[-1]['error']['pressure_rms']<max(records[0]['error']['pressure_rms']*.8,1e-7),'normal-momentum pressure p refinement failed: '+family)
            require(records[-1]['error']['strong_divergence_rms']<max(records[0]['error']['strong_divergence_rms']*.8,1e-8),'strong-divergence p diagnostic failed to improve: '+family)
        sections['manufactured'].append(dict(family=family,dimension=dimension,swirl=swirl,sequence=records))
    quad=run.curved('variable-mms-quadrature',family='variable',degree=5,steps=1,order=2,dt=.005,stationary=True,quadrature=16)
    stage_budgets(quad)
    qerror=native_history_difference(curved_cases['variable'],quad,'smooth_mms_quadrature_12_to_16',tolerance=1e-8)
    sections['smooth_quadrature_refinement']=dict(outer_orders=[12,16],maximum_scaled_physical_state_and_history_change=qerror)
    # Actual stage solver, ownership and history comparisons at one identical
    # discretization. These compare full coupled state and native IMEX history.
    options=dict(family='variable',dimension=2,degree=4,steps=3,order=2,dt=.02,quadrature=12)
    baseline=run.curved('curved-evolving',**options);stage_budgets(baseline)
    condensed=run.curved('curved-condensed',**dict(options,condensed=True));stage_budgets(condensed)
    native_history_difference(baseline,condensed,'curved_full_vs_condensed',tolerance=1e-8)
    sections['ownership']=[]
    for ranks in (2,4):
        raw=run.curved(f'curved-mpi{ranks}',ranks=ranks,ownership_shift=1,**options);stage_budgets(raw)
        native_history_difference(baseline,raw,'curved_mpi'+str(ranks),tolerance=1e-8)
        require(set(int(owner) for _,owner in raw['cell_ownership'])==set(range(ranks)),'empty or repeated MPI ownership')
        sections['ownership'].append(dict(ranks=ranks,shift=1,full_history_equal=True))
    checkpoint=run.output/'curved-save'/'history'
    saved=run.curved('curved-save',checkpoint=checkpoint,**dict(options,steps=1));stage_budgets(saved)
    restored=run.curved('curved-restart',restart=saved['checkpoint'],ranks=4,ownership_shift=2,**dict(options,steps=2));stage_budgets(restored)
    native_history_difference(baseline,restored,'curved_restart_changed_ownership',tolerance=1e-8)
    sections['restart']=dict(changed_ranks=4,changed_ownership_shift=2,complete_native_history_equal=True)
    # h study changes only element subdivision, retaining physical model data.
    coarse=curved_cases['variable']
    href=run.curved('variable-href',family='variable',degree=3,split=2,radial_split=2,steps=1,dt=.005,stationary=True)
    stage_budgets(href);eh=mms_error(href,stationary=True)
    e0=sections['manufactured'][1]['sequence'][0]['error']
    require(eh['epsilon_rms']<e0['epsilon_rms']*.6,'fixed-height h refinement did not converge')
    sections['h_refinement']=dict(coarse=e0,refined=eh,physical_heights_unchanged=True)
    # Explicit hard rejections: cannot claim a smooth tube from vertex signs,
    # hide a C1 foot in an element, or activate curved equilibrium anchoring.
    sections['rejections']=[]
    for name,kind in [('curved-anchor','eta'),('curved-chi','chi'),('folded-between-nodes','fold'),('unsplit-c1-foot','foot')]:
        path=run.output/name;path.mkdir()
        xml=mapped_session(run.template,path/'case.xml',family='feet' if kind=='foot' else 'variable',steps=0)
        root=ET.parse(xml).getroot()
        if kind=='eta':set_parameter(root,'RANSWallEta',1)
        if kind=='chi':set_parameter(root,'RANSWallChi',1)
        if kind=='fold':set_parameter(root,'RANSWallMapAmplitude',1.);set_parameter(root,'RANSWallMapLengthX',.1)
        if kind=='foot':
            from cls_epawf_m4_cases import tensor_mesh
            tensor_mesh(root,[0.,1.,2.,2.1,4.],3,2,3.,1.,xknots=[-1.5,1.5])
        write(root,xml);run.command(name,'IncNavierStokesSolver',[xml],expect_failure=True)
        sections['rejections'].append(name)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--build',type=Path,required=True);parser.add_argument('--output',type=Path)
    args=parser.parse_args();output=args.output or args.build/('cls-epawf-m5-'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S')+'-'+str(os.getpid()));run=Runner(args.build,output);root=Path(__file__).resolve().parents[4]
    revision=subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip()
    report=dict(schema=1,milestone='M5',source_commit=revision,started_utc=datetime.now(timezone.utc).isoformat(),
                python=sys.version,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,threads=THREADS,
                frozen_contract_sha256=digest(root/'docs/plans/cls-epawf-revised-contract-v1.json'),sections={},commands=[],verification_status='failed',milestone_closure=False,
                scope='analytic normal-coordinate maps on native CG charts; replicated dense reference backend; eta=chi=0 on curved cases')
    try:
        frozen={
            'docs/plans/cls-epawf-revised-contract-v1.json':'4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a',
            'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex':'d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f',
            'library/SolverUtils/Turbulence/RANSModels.cpp':'86a7154ba3cb27e4d3319aa42d156a228fd1b180f35f9b217159822ec9d66b05',
            'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp':'5baf29819fed6b81ad1b042bc09258a3e0779d15de13076b688bf32f80ddcf7e'}
        report['frozen_sources']={name:digest(root/name) for name in frozen}
        for name,expected in frozen.items():require(report['frozen_sources'][name]==expected,'frozen source changed: '+name)
        qualify(run,report['sections'])
        report.update(verification_status='passed',milestone_closure=True)
    except Exception as exc:
        report['error']=str(exc)
        raise
    finally:
        report.update(checks=checks,maximum_errors=errors,commands=run.commands,finished_utc=datetime.now(timezone.utc).isoformat())
        (run.output/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['verification_status'],checks=checks,source_commit=revision)),flush=True)

if __name__=='__main__':main()
