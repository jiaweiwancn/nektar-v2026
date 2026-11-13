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
from scipy.integrate import quad


def geometry(q,r,amplitude=.005):
    h=amplitude*math.exp(-q*q);hp=-2*q*h;hpp=(4*q*q-2)*h
    a=math.sqrt(1+hp*hp);gamma=-hpp/a**3;J=1+r*gamma
    t=np.array([1/a,0,hp/a]);n=np.array([-hp/a,0,1/a])
    return a,J,gamma,t,n


def geometry_qualification(packet,amplitude=.005):
    o=packet['observation'];sizes=o['sizes'];H=.5
    integrals={};faces={};volumes={};worst=0.;constant=0.;cofactor=0.
    def check_row(row,coordinate):
        nonlocal worst,constant
        a,J,gamma,t,n=geometry(coordinate[0],coordinate[2],amplitude)
        chart=np.asarray(row['native_chart_derivative']);d=np.asarray(row['derivative'])
        expected=t[:,None]*chart[0]/(a*J)+n[:,None]*chart[2]
        worst=max(worst,float(np.max(abs(expected-d)))/max(1.,float(np.max(abs(expected)))),float(np.max(abs(np.asarray(row['native_chart_phi'])-row['phi']))))
        c=np.asarray(row['constant']);constant=max(constant,abs(float(c@row['phi'])-1),float(np.max(abs(d@c)))/max(1.,float(np.max(abs(d)))))
    def volume(block,row,coordinate,weight):
        check_row(row,coordinate);key=(block,row['geometry']);ids=np.asarray(row['ids'],int)
        if key not in integrals:integrals[key]=np.zeros((sizes[block],3))
        integrals[key][ids]+=weight*np.asarray(row['derivative']).T
        volumes[block]=volumes.get(block,0.)+weight
    for point in o['volume']:
        volume(0,point['velocity_k'],point['coordinate'],point['weight'])
        a,J,*_=geometry(point['coordinate'][0],point['coordinate'][2],amplitude)
        worst=max(worst,abs(a*J-point['volume_metric']),abs(J-point['normal_J']))
        if point['region']:volume(3,point['epsilon'],point['coordinate'],point['weight'])
    for point in o['pressure_mass']:volume(4+point['region'],point['pressure'],point['coordinate'],point['weight'])
    for f in o['CT6_cell_faces']:
        check_row(f['row'],f['coordinate']);a,J,gamma,t,n=geometry(f['coordinate'][0],f['coordinate'][2],amplitude)
        expected=f['side']*(t if f['axis']==0 else a*J*n)
        cofactor=max(cofactor,float(np.max(abs(expected-np.asarray(f['area'])))))
        block=f['block'];key=(block,f['row']['geometry']);ids=np.asarray(f['row']['ids'],int)
        if key not in faces:faces[key]=np.zeros((sizes[block],3))
        faces[key][ids]+=f['weight']*np.asarray(f['row']['phi'])[:,None]*np.asarray(f['area'])
    require(set(integrals)==set(faces),'missing_native_element_face_inventory')
    stokes=max(float(np.max(abs(integrals[key]-faces[key])))/max(1.,float(np.max(abs(integrals[key]))),float(np.max(abs(faces[key])))) for key in integrals)
    bounds={0:(0.,1.),3:(H,1.),4:(0.,H),5:(H,1.)};volume_defects={}
    for block,(lo,hi) in bounds.items():
        expected=quad(lambda x:geometry(x,0,amplitude)[0]*((hi-lo)+geometry(x,0,amplitude)[2]*(hi*hi-lo*lo)/2),-1.2,1.2,epsabs=1e-13,epsrel=1e-13)[0]
        volume_defects[str(block)]=abs(volumes[block]-expected)/expected
    maximum=max(worst,constant,cofactor,stokes,*volume_defects.values())
    return dict(schema='CT6NativeMappedGeometryV1',status='pass' if maximum<=2e-11 else 'fail',
        basis_transform=worst,partition_and_constant_gradient=constant,cofactor=cofactor,
        all_element_all_mode_Stokes=stokes,physical_volume_defects=volume_defects,
        tested_element_field_pairs=len(integrals),scope='native chart rows independently transformed to physical Cartesian gradients; all cell faces, oriented cofactors and all-mode element Stokes identities')


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
    physical_volume=math.fsum(point['weight'] for point in o['volume'])
    require(math.isfinite(physical_volume) and physical_volume>0,'invalid_physical_volume')
    denominator=np.maximum(equation_scales*physical_volume,absolute)
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
        actual_stored_slots=len(state),physical_volume=physical_volume,consumed_volume_points=ownership,independent_rays=len(columns),
        norm='all-row maximum absolute difference / max(physical equation scale*actual integrated physical volume, sum absolute primitive increments)',
        scope='current represented state, conservative unmasked operator; prescribed-row replacement separately checked; not analytical equilibrium cancellation or spatial/temporal convergence')


def continuous_weak(packet,model):
    """Continuous fields/derivatives against the actual native test functions.

    Prescribed physical face fluxes are independently evaluated on *all*
    boundary test rows. This includes dependent outer-epsilon interface rows;
    it supplies no extra interface equation to the solver.
    """
    o=packet['observation'];off=np.asarray(o['offsets'],int);size=len(o['state'])
    r=np.zeros(size);absolute=np.zeros(size);p=parameters(model)
    def add(block,row,value,flux,weight):
        ids=np.asarray(row['ids'],int)+off[block]
        primitives=np.concatenate([np.asarray(row['phi'])[None,:]*value,
            np.asarray(row['derivative'])*np.asarray(flux)[:,None]])
        r[ids]+=weight*np.sum(primitives,axis=0);absolute[ids]+=abs(weight)*np.sum(abs(primitives),axis=0)
    for point in o['volume']:
        d=point['CT6_prescribed_data'];region=point['region'];v=point['velocity_k'];w=point['weight']
        value=np.asarray(d['continuous_values']);gradient=np.asarray(d['continuous_gradients']);T=np.asarray(d['continuous_T'])
        U=np.array([value[0],0,value[1]]);k,e=value[2:4]
        dt=np.asarray(d['continuous_dtU']);body=np.asarray(d['continuous_required_body'])
        for b,axis in enumerate((0,2)):add(b,v,dt[axis]-body[axis],T[axis]-U[axis]*U,w)
        add(4+region,point['pressure'],gradient[0,0]+gradient[1,2],[0,0,0],w)
        sk=d['continuous_production']-e+d['k_source']
        add(2,v,d['continuous_dtk']-sk,d['continuous_Dk']*gradient[2]-U*k,w)
        if region:
            G=np.zeros((3,3));G[0]=gradient[0];G[2]=gradient[1]
            law=oracle(model,k,e,G)
            se=law['eps_source']+d['epsilon_source']
            add(3,point['epsilon'],d['continuous_dtepsilon']-se,d['continuous_Depsilon']*gradient[3]-U*e,w)
    interface_k=np.zeros(size);interface_momentum=np.zeros(size)
    for face in o['faces']:
        value=np.asarray(face['CT6_continuous']['values']);g=np.asarray(face['CT6_continuous']['gradients']);T=np.asarray(face['CT6_continuous']['T'])
        U=np.array([value[0],0,value[1]]);area=np.asarray(face['area']);un=float(U@area);w=face['weight'];d=face['CT6_continuous']
        for b,axis in enumerate((0,2)):
            flux=un*U[axis]-float(T[axis]@area)
            add(b,face['velocity_k'],flux,[0,0,0],w)
            if face['interface']:
                ids=np.asarray(face['velocity_k']['ids'],int)+off[b]
                interface_momentum[ids]+=w*np.asarray(face['velocity_k']['phi'])*flux
        kflux=un*value[2]-d['Dk']*float(g[2]@area)
        add(2,face['velocity_k'],kflux,[0,0,0],w)
        if face['interface']:
            ids=np.asarray(face['velocity_k']['ids'],int)+off[2]
            interface_k[ids]+=w*np.asarray(face['velocity_k']['phi'])*kflux
        if face['region']:add(3,face['epsilon'],un*value[3]-d['Depsilon']*float(g[3]@area),[0,0,0],w)
    scale=np.empty(size);physical_volume=math.fsum(point['weight'] for point in o['volume'])
    require(math.isfinite(physical_volume) and physical_volume>0,'invalid_physical_volume')
    for b,s in enumerate([.2,.2,p['q']*.16,.4**3/(p['kappa']*.5),.2,.2]):scale[off[b]:off[b+1]]=s*physical_volume
    scale[-1]=1.;denominator=np.maximum(scale,absolute)
    defects=abs(r)/denominator
    maximum=float(np.max(defects));interface=max(float(np.max(abs(interface_k)/scale)),float(np.max(abs(interface_momentum)/scale)))
    return dict(schema='CT6ContinuousNativeWeakV1',model=model,status='pass' if max(maximum,interface)<=2e-10 else 'fail',
        maximum_scaled_all_row_defect=maximum,shared_interface_cancellation=interface,
        witness=dict(row=int(np.argmax(defects)),signed_load=float(r[np.argmax(defects)]),scale=float(denominator[np.argmax(defects)])),
        block_maximum=[float(max(defects[off[b]:off[b+1]],default=0)) for b in range(6)],
        scope='analytical continuous fields and independently differentiated forcing on native test functions; all physical face fluxes included; projected/evolved residual is a separate object')


def frozen_systems(packet,model,amplitude=.005):
    """Original frozen Oseen/reaction matrices and actual free-mass Schur rank."""
    o=packet['observation'];state=np.asarray(o['state']);off=np.asarray(o['offsets'],int)
    n=len(state);A=np.zeros((n,n));p=parameters(model);columns={};mass=packet['open_frozen_systems']['masses']
    for b in range(4):A[off[b]:off[b+1],off[b]:off[b+1]]=np.asarray(mass[b])
    for point in o['volume']:
        v=point['velocity_k'];pr=point['pressure'];region=point['region'];weight=point['weight']
        data=point['CT6_prescribed_data'];ids=np.asarray(v['ids'],int);phi=np.asarray(v['phi']);d=np.asarray(v['derivative'])
        U=np.zeros(3);G=np.zeros((3,3))
        for b,axis in enumerate((0,2)):U[axis],G[axis]=row_value(v,state,off[b])
        k,_=row_value(v,state,off[2]);adv=d.T@U
        if region:
            e,_=row_value(point['epsilon'],state,off[3]);law=oracle(model,k,e,G)
            viscosity=law['nut'];rates={2:(e-min(0.,law['production']))/k,3:(1.9 if model=='E94E95' else 1.92)*e/k}
            diffusion={2:viscosity,3:viscosity/p['sigma_epsilon']}
        else:
            K,_=row_value(data['H_velocity_k_row'],state,off[2]);key=(point['coordinate'][0],K)
            if key not in columns:
                *_,gamma,tangent,normal=geometry(point['coordinate'][0],0,amplitude)
                columns[key]=Column(model,K,.5,.05,1.5e-5,curvature=-gamma,order=96)
            viscosity=float(columns[key].D(point['coordinate'][2]));rates={2:0.};diffusion={2:viscosity}
        pp=np.asarray(pr['ids'],int)+off[4+region];pphi=np.asarray(pr['phi'])
        for b,axis in enumerate((0,2)):
            ii=ids+off[b]
            for other,jaxis in enumerate((0,2)):
                block=viscosity*np.outer(d[jaxis],d[axis])
                if b==other:block+=viscosity*d.T@d-np.outer(adv,phi)
                A[np.ix_(ii,ids+off[other])]+=weight*block
            A[np.ix_(ii,pp)]-=weight*np.outer(d[axis],pphi)
            A[np.ix_(pp,ii)]+=weight*np.outer(pphi,d[axis])
            if not region and point['coordinate'][2]<=.05:
                hp=data['H_pressure_row'];_,_,_,_,normal=geometry(point['coordinate'][0],.5,amplitude)
                pressure_gradient=np.asarray(hp['derivative']);P=np.eye(3)-np.outer(normal,normal)
                A[np.ix_(ii,np.asarray(hp['ids'],int)+off[5])]-=weight*np.outer(phi,(P@pressure_gradient)[axis])
        for b in ((2,3) if region else (2,)):
            row=point['epsilon'] if b==3 else v;ix=np.asarray(row['ids'],int)+off[b]
            ph=np.asarray(row['phi']);der=np.asarray(row['derivative'])
            A[np.ix_(ix,ix)]+=weight*(rates[b]*np.outer(ph,ph)-np.outer(der.T@U,ph)+diffusion[b]*der.T@der)
    for face in packet['open_boundary']['faces']:
        U=np.array([face['velocity'][0],0,face['velocity'][1]]);un=float(U@np.asarray(face['area']))
        for b in ((0,1,2,3) if face['region'] else (0,1,2)):
            row=face['epsilon_row'] if b==3 else face['velocity_k_row'];ii=np.asarray(row['ids'],int)+off[b];ph=np.asarray(row['phi'])
            A[np.ix_(ii,ii)]+=face['weight']*un*np.outer(ph,ph)
    blocks={};flow=None;constraints={int(i):v for i,v in packet['resolution']['constraints']}
    for block in packet['open_frozen_systems']['blocks']:
        b=block['block'];ids=np.asarray(block['rows'],int);actual=np.asarray(block['matrix']);expected=A[np.ix_(ids,ids)].copy()
        if b==3:
            for j,row in enumerate(ids):
                if row in constraints:expected[j,:]=0;expected[j,j]=1
        normal=max(1.,float(np.max(abs(actual))));error=float(np.max(abs(actual-expected)))/normal
        residual=np.asarray(block['residual']);raw=np.asarray(packet['open_frozen_systems']['raw_residual'])[ids]
        for j,row in enumerate(ids):
            if b==3 and row in constraints:raw[j]=state[row]-constraints[row]
        blocks[str(b)]=dict(matrix_scaled_defect=error,residual_scaled_defect=float(np.max(abs(residual-raw)))/max(1.,float(np.max(abs(raw)))))
        if b==0:flow=(ids,actual)
    ids,S=flow;nv=int(packet['mixed_space']['velocity_size']);B=S[:nv,nv:];C=S[nv:,:nv]
    M=np.zeros((nv,nv))
    for i,gi in enumerate(ids[:nv]):
        b=0 if gi<off[1] else 1
        for j,gj in enumerate(ids[:nv]):
            if off[b]<=gj<off[b+1]:M[i,j]=mass[b][gi-off[b]][gj-off[b]]
    schur=C@np.linalg.solve(M,B);row_scale=np.maximum(np.linalg.norm(schur,axis=1),np.finfo(float).tiny)
    normalized=schur/row_scale[:,None];col_scale=np.maximum(np.linalg.norm(normalized,axis=0),np.finfo(float).tiny);normalized/=col_scale[None,:]
    singular=np.linalg.svd(normalized,compute_uv=False);relative=float(singular[-1]/singular[0])
    constant=np.asarray(packet['mixed_space']['pressure_constant_with_zero_gauge'])[nv:]
    constant_action=float(np.linalg.norm(B@constant)/np.linalg.norm(constant))
    maximum=max(max(row.values()) for row in blocks.values())
    return dict(schema='CT6IndependentFrozenSystemsV1',status='pass' if maximum<=2e-10 and relative>1e-11 and constant_action>1e-8 else 'fail',
        blocks=blocks,pressure_unknowns=len(constant),free_velocity_unknowns=nv,
        free_mass_symmetry=float(np.max(abs(M-M.T))),minimum_free_mass_eigenvalue=float(np.linalg.eigvalsh(M)[0]),
        normalized_consistent_mass_pressure_minimum_singular_ratio=relative,
        constant_pressure_action=constant_action,additional_pressure_null_modes=int(np.sum(singular<=singular[0]*1e-11)),
        scope='actual constrained matrix, projected roughness-volume pressure columns in both Cartesian components, prescribed outlet advection and frozen scalar reaction; consistent free-mass pressure action, no mean multiplier; finite fixture only')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--packet',type=Path,required=True)
    p.add_argument('--model',choices=['A93','A95','E94E95'],required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--amplitude',type=float,default=.005);a=p.parse_args()
    packet=json.loads(a.packet.read_text());result=reconstruct(packet,a.model,a.amplitude)
    result['continuous_weak']=continuous_weak(packet,a.model)
    result['geometry']=geometry_qualification(packet,a.amplitude)
    result['status']='pass' if result['status']==result['continuous_weak']['status']==result['geometry']['status']=='pass' else 'fail'
    write_json(a.output,result)
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
