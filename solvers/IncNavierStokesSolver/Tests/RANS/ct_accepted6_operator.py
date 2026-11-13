#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT6 all-row conservative load reconstruction, independent original laws.

Uses exported native geometry/basis and prescribed data, never subtracts a
reference residual or modifies an evolving state. Fixed body panels and moving
wall-source panels remain separate.
"""
import argparse,json,math,os
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
from ct_accepted_feasibility import Column,parameters,require
from ct_accepted_open import row_value
from validate_a93_a95_e94e95 import oracle


def geometry(q,r,amplitude=.005):
    h=amplitude*math.exp(-q*q);hp=-2*q*h;hpp=(4*q*q-2)*h
    a=math.sqrt(1+hp*hp);gamma=-hpp/a**3;J=1+r*gamma
    t=np.array([1/a,0,hp/a]);n=np.array([-hp/a,0,1/a])
    return a,J,gamma,t,n


def reconstruct(packet,model,amplitude=.005):
    o=packet['observation'];require('CT6' in o['identity'],'foreign_CT6_packet')
    state=np.asarray(o['state'],float);off=np.asarray(o['offsets'],int)
    total=np.zeros(len(state));absolute=np.zeros(len(state));columns={};ownership=0;geom_error=0.
    p=parameters(model);nu=1.5e-5;H=.5;hr=.05
    def add(block,row,scalar,flux,weight):
        ids=np.asarray(row['ids'],int)+off[block]
        parts=np.concatenate([np.asarray(row['phi'])[None,:]*scalar,np.asarray(row['derivative'])*np.asarray(flux)[:,None]])
        total[ids]+=weight*np.sum(parts,axis=0);absolute[ids]+=abs(weight)*np.sum(abs(parts),axis=0)
    for point in o['volume']:
        v=point['velocity_k'];pr=point['pressure'];region=int(point['region']);data=point['CT6_prescribed_data']
        require(point['generating_geometry']==v['geometry']==pr['geometry'],'volume_owner_mismatch')
        if region:require(point['epsilon']['geometry']==v['geometry'],'epsilon_owner_mismatch')
        ownership+=1
        chart=point['coordinate'];a,J,gamma,t,n=geometry(chart[0],chart[2],amplitude)
        geom_error=max(geom_error,abs(a*J-point['volume_metric']),abs(J-point['normal_J']))
        U=np.zeros(3);G=np.zeros((3,3))
        for b,axis in enumerate((0,2)):U[axis],G[axis]=row_value(v,state,off[b])
        k,gk=row_value(v,state,off[2]);pi,_=row_value(pr,state,off[4+region]);weight=point['weight']
        if region:
            e,ge=row_value(point['epsilon'],state,off[3]);law=oracle(model,k,e,G)
            T=-law['R']+np.eye(3)*(2*k/3-pi);Dk=law['nut'];De=Dk/p['sigma_epsilon']
            sk=law['k_source']+data['k_source'];se=law['eps_source']+data['epsilon_source'];rough=np.zeros(3)
        else:
            hrow=data['H_velocity_k_row'];hpressure=data['H_pressure_row']
            K,_=row_value(hrow,state,off[2]);key=(chart[0],K)
            if key not in columns:columns[key]=Column(model,K,H,hr,nu,curvature=-gamma,order=96)
            c=columns[key];require(c.cutoff is not None,'independent_current_cutoff_inadmissible')
            nt=float(c.nt(chart[2]));Dk=nu+nt;T=Dk*(G+G.T)-np.eye(3)*pi
            production=nt*float(np.sum((G+G.T)*G))-2*k*np.trace(G)/3
            sk=production-c.A/max(chart[2],c.cutoff)+data['k_source'];rough=np.zeros(3)
            if chart[2]<=hr:
                _,Jh,gamma,t,n=geometry(chart[0],H,amplitude);P=np.eye(3)-np.outer(n,n);Br=gamma/Jh*np.outer(t,t)
                geom_error=max(geom_error,float(np.max(abs(np.asarray(data['H_normal'])-n))),float(np.max(abs(np.asarray(data['H_Br'])-Br))))
                Uh=np.zeros(3);Gh=np.zeros((3,3))
                for b,axis in enumerate((0,2)):Uh[axis],Gh[axis]=row_value(hrow,state,off[b])
                _,gp=row_value(hpressure,state,off[5]);_,gkh=row_value(hrow,state,off[2])
                w=float(Uh@n);ut=P@Uh;gradient=P@Gh-np.outer(n,Br.T@Uh)-w*Br
                rough=P@(gp-2*gkh/3+ut*np.trace(gradient)+gradient@ut)
        for b,axis in enumerate((0,2)):add(b,v,-rough[axis],T[axis]-U[axis]*U,weight)
        add(4+region,pr,np.trace(G),[0,0,0],weight)
        add(2,v,-sk,Dk*gk-U*k,weight)
        if region:add(3,point['epsilon'],-se,De*ge-U*e,weight)
    require(len(o['CT6_body_volume'])>0,'fixed_body_panels_missing')
    for point in o['CT6_body_volume']:
        data=point['CT6_prescribed_data'];body=np.array([data['body'][0],0,data['body'][1]])
        if not point['region'] and point['coordinate'][2]<=hr:
            _,_,_,_,n=geometry(point['coordinate'][0],H,amplitude)
            out=np.array([data['outer_body_at_H'][0],0,data['outer_body_at_H'][1]])
            body-=out-n*float(out@n)
        for b,axis in enumerate((0,2)):add(b,point['velocity_k'],-body[axis],[0,0,0],point['weight'])
    # The physical outlet has both advective and natural flux contributions.
    for face in packet['open_boundary']['faces']:
        area=np.asarray(face['area']);norm=np.linalg.norm(area);weight=face['weight']
        U=np.array([face['velocity'][0],0,face['velocity'][1]]);un=float(U@area)
        for b in (0,1):add(b,face['velocity_k_row'],un*face['velocity'][b]-norm*face['traction'][b],[0,0,0],weight)
        add(2,face['velocity_k_row'],un*face['k']-norm*face['k_diffusive_datum'],[0,0,0],weight)
        if face['region']:add(3,face['epsilon_row'],un*face['epsilon']-norm*face['epsilon_diffusive_datum'],[0,0,0],weight)
    native=np.asarray(packet['resolution']['native']['spatial_unmasked'])
    require(native.shape==total.shape,'native_load_size_mismatch')
    # Same physical scaling on every row; accumulated absolute primitives are
    # retained to distinguish cancellation from a small physical equation.
    equation_scales=np.empty(len(total));scales=[.2,.2,p['q']*.16,.4**3/(p['kappa']*.5),.2,.2,1.]
    for b in range(len(off)-1):equation_scales[off[b]:off[b+1]]=scales[b]
    equation_scales[-1]=1.
    denominator=np.maximum(equation_scales*2.4,absolute)
    scaled=abs(native-total)/denominator
    constrained=total.copy()
    for row,value in packet['resolution']['constraints']:constrained[int(row)]=state[int(row)]-value
    cn=np.asarray(packet['resolution']['native']['stage_constrained'])
    constrained_defect=float(np.max(abs(constrained-cn)/denominator))
    primitive_sum=np.asarray(packet['resolution']['terms']['spatial_sum'])
    reconstruction=float(np.max(abs(native-primitive_sum)/denominator))
    maximum=float(np.max(scaled));status='pass' if max(maximum,constrained_defect,reconstruction,geom_error)<=2e-10 else 'fail'
    return dict(schema='CT6IndependentAllRowOperatorV1',status=status,model=model,
        maximum_scaled_defect=maximum,constrained_scaled_defect=constrained_defect,
        native_term_sum_scaled_defect=reconstruction,geometry_primitive_defect=geom_error,
        witness=dict(row=int(np.argmax(scaled)),native=float(native[np.argmax(scaled)]),independent=float(total[np.argmax(scaled)]),scale=float(denominator[np.argmax(scaled)])),
        actual_stored_slots=len(state),consumed_volume_points=ownership,independent_rays=len(columns),
        norm='all-row maximum absolute difference / max(physical equation scale*2.4, sum absolute primitive increments)',
        scope='current represented state, conservative unmasked operator; prescribed-row replacement separately checked; not analytical equilibrium cancellation or spatial/temporal convergence')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--model',choices=['A93','A95','E94E95'],required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--amplitude',type=float,default=.005);a=p.parse_args()
    result=reconstruct(json.loads(a.packet.read_text()),a.model,a.amplitude);write_json(a.output,result)
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
