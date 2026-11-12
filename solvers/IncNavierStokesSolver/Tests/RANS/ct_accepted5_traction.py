# SPDX-License-Identifier: MIT
"""Read-only actual/stage-frozen tractions on native owned face rows."""
import math,os
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted_open import row_value,require
from ct_accepted4_reference import parameters
from ct_accepted3_metrics import point

def endpoint_trace(initial,state,c):
    o=initial['observation'];traces=[];p=c['degree'];off=o['offsets'];H=c['H']
    for lo,hi in zip(c['x'][:-1],c['x'][1:]):
        xs=[lo+(hi-lo)*j/p for j in range(p+1)];ks=[]
        for x in xs:
            i=next(i for i,(q,side) in enumerate(zip(o['grid'],o['grid_sides'],strict=True)) if q==[x,0.,H] and side==1)
            ks.append(row_value(initial['test_basis'][i]['velocity_k'],state,off[2])[0])
        traces.append((lo,hi,np.polynomial.polynomial.polyfit((np.array(xs)-lo)/(hi-lo),ks,p)))
    def at(x):
        lo,hi,coef=next(t for t in traces if t[0]<=x<=t[1]);v=float(np.polynomial.polynomial.polyval((x-lo)/(hi-lo),coef));require(v>0,'invalid_actual_endpoint_trace');return v
    return at

def diffusivity(K,z,c):
    v=c['roughness_input']*math.sqrt(K)/(70*c['nu']);b=10.7 if v==0 else 10.7*(1-v**max(.5-.4*v**.7,1-.79*v**(-.28)))
    return c['nu']+parameters(c['model'])['alpha']*math.sqrt(K)*max(0.,z-c['nu']*b/math.sqrt(K))

def profiles(initial,state,c,frozen_endpoint_state=None):
    o=initial['observation'];off=o['offsets'];s=np.array(state);K=endpoint_trace(initial,s,c)
    frozen=endpoint_trace(initial,np.array(frozen_endpoint_state),c) if frozen_endpoint_state is not None else K
    ground=[];interface=[];undefined=[]
    for f in o['faces']:
        if not f['interface'] and f['coordinate'][2]!=0:continue
        x,_,z=f['coordinate'];G=np.zeros((3,3))
        for b,axis in enumerate((0,2)):G[axis]=row_value(f['velocity_k'],s,off[b])[1]
        kH=K(x);D=diffusivity(kH,z,c);Df=diffusivity(frozen(x),z,c)
        row=dict(coordinate=f['coordinate'],weight=f['weight'],geometry=f['velocity_k']['geometry'],region=f['region'],side=-1 if not f['region'] else 1,
            endpoint_K=kH,D=D,S=float(G[0,2]+G[2,0]),raw_gradient=G.tolist(),signed_divergence=float(np.trace(G)),
            actual_traction=D*float(G[0,2]+G[2,0]),counterfactual_frozen_D_traction=Df*float(G[0,2]+G[2,0]))
        if not f['interface']:ground.append(row);continue
        pi=row_value(f['pressure'],s,off[4+f['region']])[0]
        if f['region']:
            law=point(f,s,off,c['model'])
            if 'total_stress' not in law:row.update(status='undefined_outer_stress',invariants=law);undefined.append(row);interface.append(row);continue
            T=np.array(law['total_stress'])
        else:T=D*(G+G.T)-np.eye(3)*pi
        row.update(status='defined',Pi=pi,total_traction=T[:,2].tolist());interface.append(row)
    w=np.array([q['weight'] for q in ground]);t=np.array([q['actual_traction'] for q in ground]);mean=float(w@t/w.sum());ref=c['ustar']**2
    pairs={}
    for q in interface:pairs.setdefault(round(q['coordinate'][0],14),[]).append(q)
    jumps=[]
    for x,qs in sorted(pairs.items()):
        require(len(qs)==2 and {q['side'] for q in qs}=={-1,1},'missing_one_sided_interface_pair')
        if any(q['status']!='defined' for q in qs):continue
        inner=next(q for q in qs if q['side']==-1);outer=next(q for q in qs if q['side']==1)
        jumps.append(dict(coordinate=inner['coordinate'],weight=inner['weight'],jump=(np.array(inner['total_traction'])-outer['total_traction']).tolist(),geometries=[inner['geometry'],outer['geometry']]))
    jr=math.sqrt(sum(q['weight']*sum(v*v for v in q['jump']) for q in jumps)/sum(q['weight'] for q in jumps)) if jumps and not undefined else None
    return dict(state_sha256=exact_vector_hash(state),frozen_endpoint_state_sha256=exact_vector_hash(frozen_endpoint_state) if frozen_endpoint_state is not None else None,
        ground=ground,interface=interface,interface_jump=jumps,undefined_interface=undefined,mean_wall_shear=mean,
        normalized_wall_traction_error=abs(mean-ref)/ref,normalized_interface_traction_RMS=jr/ref if jr is not None else None,
        wall_spatial_standard_deviation=float(np.sqrt(w@((t-mean)**2)/w.sum())),
        minimum_wall=min(ground,key=lambda q:q['actual_traction']),maximum_wall=max(ground,key=lambda q:q['actual_traction']),
        scope='actual current endpoint K and FE gradients; frozen-D evaluation explicitly counterfactual; no constitutive evaluation needed for ground traction')

def delta(a,b,c):
    require(len(a['ground'])==len(b['ground']),'ground_inventory_length_mismatch');w=[];ds=[];dk=[];cross=[];actual=[]
    witness=[]
    for x,y in zip(a['ground'],b['ground'],strict=True):
        require(x['coordinate']==y['coordinate'] and x['weight']==y['weight'],'ground_inventory_identity_mismatch')
        D,S=x['D'],x['S'];dD=y['D']-D;dS=y['S']-S;v=(D*dS,S*dD,dD*dS);change=y['actual_traction']-x['actual_traction']
        require(abs(sum(v)-change)<2e-10*c['ustar']**2,'traction_stage_decomposition_failed')
        w.append(x['weight']);ds.append(v[0]);dk.append(v[1]);cross.append(v[2]);actual.append(change)
        witness.append(dict(coordinate=x['coordinate'],geometry=x['geometry'],D0=D,S0=S,DeltaD=dD,DeltaS=dS,gradient=v[0],coefficient=v[1],cross=v[2],actual=change,identity_defect=sum(v)-change))
    w=np.array(w);mean=lambda v:float(w@np.array(v)/w.sum())
    return dict(gradient=mean(ds),coefficient=mean(dk),cross=mean(cross),actual=mean(actual),
        identity_defect=mean(ds)+mean(dk)+mean(cross)-mean(actual),witness=max(witness,key=lambda q:abs(q['actual'])),
        first_state_sha256=a['state_sha256'],second_state_sha256=b['state_sha256'])

def difference(a,b,c):
    require(len(a['ground'])==len(b['ground']) and len(a['interface'])==len(b['interface']),'traction_profile_coverage')
    groups={}
    for name in ('ground','interface'):
        num=den=0.;maximum=0.
        for p,q in zip(a[name],b[name],strict=True):
            require(p['coordinate']==q['coordinate'] and p['weight']==q['weight'] and p['side']==q['side'],'traction_sample_identity')
            d=np.array([p['actual_traction']-q['actual_traction']]) if name=='ground' else np.array(p['total_traction'])-q['total_traction']
            num+=p['weight']*float(d@d);den+=p['weight'];maximum=max(maximum,float(np.linalg.norm(d)))
        groups[name]=dict(normalized_weighted_RMS=math.sqrt(num/den)/c['ustar']**2,normalized_maximum=maximum/c['ustar']**2)
    return dict(signed_mean_difference=(a['mean_wall_shear']-b['mean_wall_shear'])/c['ustar']**2,profiles=groups)
