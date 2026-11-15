# SPDX-License-Identifier: MIT
"""Physical-time sensitivity, distinct from H18's exact trajectory replay."""
import copy
import math
import xml.etree.ElementTree as ET
import numpy as np
from study_wj_h18 import protected
from wj_h4_points import aligned
from wj_h7_analysis import context
from wj_h3_analysis import algebra

GATE=.002068763923454
FLOOR=4e-11
LIMITS=dict(U=1e-4,k=1e-3,R=1e-3,pi=1e-3,G=1e-2,epsilon=1e-2)
MEASURES=dict(full=497.8125,inner=3.75,outer=494.0625)

def spatial_xml(x):
    x=copy.deepcopy(x)
    nodes=[p for p in x.findall('./CONDITIONS/PARAMETERS/P')if p.text.split('=')[0].strip()=='TimeStep']
    assert len(nodes)==1
    x.find('./CONDITIONS/PARAMETERS').remove(nodes[0])
    return protected(x)

def physical_pair(coarse,fine,coarse_step):
    assert coarse_step in (1,2)and coarse['dt']==.005 and fine['dt']==.0025
    a=context(coarse,coarse_step);b=context(fine,2*coarse_step)
    assert a['step']==coarse_step and b['step']==2*coarse_step
    aligned(coarse_step,.005,a['time']);aligned(2*coarse_step,.0025,b['time'])
    assert abs(a['time']-b['time'])<=128*np.finfo(float).eps*max(1.,abs(a['time']),abs(b['time']))
    assert a['alpha']==200 and b['alpha']==400,'each completed callback has its own alpha'
    return a,b

def fixed_action(k,e,G,theta):
    x=np.cos(theta);z=np.sin(theta)
    xi=np.stack((x,np.zeros(len(x)),z),axis=1)
    a=np.stack((-z,np.zeros(len(x)),x),axis=1);H=a[:,:,None]*xi[:,None,:]
    R,nt,N=algebra(k,e,G)
    tangent=algebra(k,e,G.astype(complex)+1e-27j*H)[0].imag/1e-27
    d=-np.einsum('nij,nij->n',H,tangent)
    return dict(dimensional=d,nu_eff=nt,normalized=d/nt,xi=xi,a=a)

def classification(value,limit,floor=FLOOR):
    assert np.isfinite(value)and limit>0
    return dict(value=float(value),limit=float(limit),utilization=float(value/limit),floor=floor,
        resolved=bool(value>floor),status='PASS'if value<=limit else'EXCEEDED')

def principal_screen(coarse,fine,outer):
    ids=np.flatnonzero(outer);budget=.1*(fine['margin'][ids]-GATE)
    records={};arrays={};unresolved=bool(np.any(budget<=FLOOR));failed=False
    fields=lambda d:(d['k'][ids],d['epsilon'][ids],d['G'][ids])
    differences={'minimized':coarse['margin'][ids]-fine['margin'][ids]}
    for name,angle in (('fine_direction',fine['angle'][ids]),('coarse_direction',coarse['angle'][ids])):
        a=fixed_action(*fields(coarse),angle);b=fixed_action(*fields(fine),angle)
        differences[name]=a['normalized']-b['normalized']
        for prefix,v in (('coarse',a),('fine',b)):
            for key,value in v.items():arrays[f'{name}_{prefix}_{key}']=value
    arrays.update(outer_ids=ids,budget=budget)
    for name,delta in differences.items():
        utilization=np.divide(abs(delta),budget,out=np.full(len(ids),np.inf),where=budget>FLOOR)
        j=int(np.argmax(utilization));arrays[name+'_signed_difference']=delta
        failed|=bool(np.any((budget>FLOOR)&(abs(delta)>budget)))
        records[name]=dict(max_utilization=float(utilization[j])if np.isfinite(utilization[j])else None,row=int(ids[j]),signed_difference=float(delta[j]),budget=float(budget[j]),
            coarse_minimized=float(coarse['margin'][ids[j]]),fine_minimized=float(fine['margin'][ids[j]]))
    return dict(status='UNRESOLVED'if unresolved else'EXCEEDED'if failed else'PASS',gate=GATE,floor=FLOOR,checks=records),arrays

def scales(rows):
    c=rows.config
    return dict(U=1.,k=c.K,R=c.u_tau**2,pi=c.u_tau**2,
        G=np.maximum(np.linalg.norm(rows.exact['G'],axis=(1,2)),c.u_tau/19.9125),epsilon=rows.exact['epsilon'])

def metrics(rows,a,b):
    out={};n=len(rows.rows);scale=scales(rows)
    for region,mask in (('full',np.ones(n,bool)),('inner',~rows.outer),('outer',rows.outer)):
        out[region]={}
        for name in LIMITS:
            if name=='epsilon'and region=='inner':continue
            active=mask&(rows.outer if name=='epsilon'else True)
            delta=np.asarray(a[name])-b[name];mag=np.linalg.norm(delta.reshape(n,-1),axis=1)
            s=scale[name]if np.isscalar(scale[name])else scale[name][active]
            normalized=mag[active]/s;weights=rows.weights[active];measure=float(np.sum(weights))
            expected=MEASURES['outer'if name=='epsilon'else region]
            assert abs(measure-expected)<=1e-10*expected,'common physical region measure differs'
            rms=math.sqrt(float(weights@normalized**2)/measure);j=int(np.argmax(normalized));ids=np.flatnonzero(active)
            out[region][name]=dict(rms=classification(rms,LIMITS[name]),maximum=classification(normalized[j],LIMITS[name]),
                row=int(ids[j]),physical_measure=measure,unnormalized_scaled_L2=rms*math.sqrt(measure),L2_floor=FLOOR*math.sqrt(measure))
    return out

def compare(rows,coarse,fine,initial):
    temporal=metrics(rows,coarse,fine);principal,arrays=principal_screen(coarse,fine,rows.outer)
    exceeded=any(v[m]['status']=='EXCEEDED'for r in temporal.values()for v in r.values()for m in ('rms','maximum'))
    status='UNRESOLVED'if principal['status']=='UNRESOLVED'else'EXCEEDED'if exceeded or principal['status']=='EXCEEDED'else'PASS'
    result=dict(schema='WJH19ProvisionalSensitivityV1',status=status,temporal=temporal,principal=principal,
        reference_error=dict(coarse=metrics(rows,coarse,rows.exact),fine=metrics(rows,fine,rows.exact)),
        drift=dict(coarse=metrics(rows,coarse,initial),fine=metrics(rows,fine,initial)),
        representation=metrics(rows,initial,rows.exact),scales=dict(U=1.,k=rows.config.K,R=rows.config.u_tau**2,pi=rows.config.u_tau**2,
            G='max(||G_eq||_F,.051/19.9125)',epsilon='positive epsilon_eq(r), outer only'),
        coarse_interface=rows.interface_history(coarse),fine_interface=rows.interface_history(fine),
        scope='provisional numerical review; no application error bound, production timestep or temporal order')
    # Signed component means and finite trace/profile diagnostics supplement the
    # vector/tensor screens. They are not extra physical acceptance thresholds.
    signed={}
    for region,mask in (('full',np.ones(len(rows.rows),bool)),('inner',~rows.outer),('outer',rows.outer)):
        signed[region]={}
        for name in ('U','G','k','R','pi','production','K','rc','destruction'):
            delta=coarse[name]-fine[name];w=rows.weights[mask];v=delta[mask]
            signed[region][name]=dict(mean=(np.tensordot(w,v,axes=1)/np.sum(w)).tolist(),
                minimum=np.min(v,axis=0).tolist(),maximum=np.max(v,axis=0).tolist())
    result['signed_components']=signed
    profiles=[];boundary=[];npnt=rows.point_count
    xs=sorted({float(e[k][0])for e in rows.report['fields'][0]['mesh']for k in ('lower','upper')})
    for x in xs:
        ids=[i for i,(p,o)in enumerate(rows.rows[:npnt])if p[0]==x]
        if not ids:continue
        profiles.append(dict(x=x,owned_rows=ids,max_delta_U=float(np.max(np.linalg.norm(coarse['U'][ids]-fine['U'][ids],axis=1))),
            max_delta_G=float(np.max(np.linalg.norm(coarse['G'][ids]-fine['G'][ids],axis=(1,2)))),
            evidence='all signed physical-height/owner profile components in the hash-bound NPZ; no trace averaging'))
    xmin,xmax=min(xs),max(xs)
    for i,(p,o)in enumerate(rows.rows[:npnt]):
        roles=[]
        if p[0]==xmin:roles.append(('inlet',0,-1.))
        if p[0]==xmax:roles.append(('outlet',0,1.))
        if o['chart'][2]==0:roles.append(('ground',2,-1.))
        if abs(p[2]-20)<1e-12:roles.append(('top',2,1.))
        for role,axis,sign in roles:
            boundary.append(dict(row=i,role=role,physical=p,owner=o,
                delta_mass_flux=float(sign*(coarse['U'][i,axis]-fine['U'][i,axis])),
                delta_total_traction=(sign*(coarse['T'][i,:,axis]-fine['T'][i,:,axis])).tolist(),
                delta_represented_k_diffusive_flux=float(sign*(coarse['nu_d'][i]*coarse['grad_k'][i,axis]-fine['nu_d'][i]*fine['grad_k'][i,axis])),
                delta_represented_epsilon_diffusive_flux=float(sign*(coarse['nu_d'][i]*coarse['grad_epsilon'][i,axis]-fine['nu_d'][i]*fine['grad_epsilon'][i,axis])/1.3)if rows.outer[i]else None))
    result.update(station_profiles=profiles,exterior_trace_differences=boundary,
        exterior_scope='represented one-sided flux/traction traces; prescribed data fixed. These are not independently integrated weak boundary loads.')
    return result,arrays
