#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Physical-point, all-owner H4 comparisons on the flat composite tunnel.

Coordinates define comparison identity. Element bounds describe provenance,
never whether two different meshes may be compared. No governing state changes.
"""
from fractions import Fraction
import hashlib
import itertools
import json
import math
import numpy as np
from validate_cls_epawf_m3 import basis
from wj_h2_profiles import CD,NU,wall,constitutive
from wj_h3_analysis import principal

EPS=np.finfo(float).eps
OFFSET=1e-6
def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def vector_hash(v):return hashlib.sha256(np.asarray(v,dtype='<f8').tobytes()).hexdigest()


def mapping(report):
    signature=report['geometry_signature']
    if ':CompositeTunnelV1:' not in signature:raise ValueError('H4 requires the declared flat composite map')
    d0,top=map(float,signature.rsplit(':',2)[-2:]);H=float(report['parameters'][1])
    if 'map_parameters' in report and report['map_parameters'][11]!=0:raise ValueError('H4 is flat only')
    if not 0<H<top-d0:raise ValueError('invalid composite height')
    return H,d0,top,(top-d0-H)/(top-H)


def chart(report,physical):
    x,y,z=map(float,physical);H,d0,top,J=mapping(report)
    if not all(map(math.isfinite,(x,y,z))) or y!=0:raise ValueError('finite planar physical point required')
    rho=z-d0 if z<=d0+H else H+(z-d0-H)/J
    return np.array([x,0.,rho])


def physical(report,q):
    H,d0,top,J=mapping(report);x,_,r=map(float,q)
    return [x,0.,d0+r if r<=H else d0+H+(r-H)*J]


def owner_point(report,element,p):
    q=chart(report,p);lo=np.asarray(element['lower']);hi=np.asarray(element['upper'])
    tolerance=64*EPS*np.maximum(1,np.maximum(abs(lo),abs(hi)))
    for axis in (0,2):
        if q[axis]<lo[axis]-tolerance[axis] or q[axis]>hi[axis]+tolerance[axis]:
            raise ValueError('physical point outside declared owner; extrapolation forbidden')
    snapped=q.copy()
    for axis in (0,2):
        if abs(q[axis]-lo[axis])<=tolerance[axis]:snapped[axis]=lo[axis]
        elif abs(q[axis]-hi[axis])<=tolerance[axis]:snapped[axis]=hi[axis]
    side=lambda axis:'lower' if snapped[axis]==lo[axis] else 'upper' if snapped[axis]==hi[axis] else 'interior'
    H=mapping(report)[0];region='inner' if hi[2]<=H else 'outer'
    return dict(geometry=element['geometry'],element=element['geometry'],region=region,
        x_side=side(0),normal_side=side(2),interface_role=('matching_'+region) if snapped[2]==H else 'volume_or_exterior',
        chart=snapped.tolist(),requested_chart=q.tolist(),roundoff_snap=(snapped-q).tolist(),
        lower=lo.tolist(),upper=hi.tolist(),physical_lower=physical(report,lo),physical_upper=physical(report,hi))


def owners(report,p):
    result=[]
    for element in report['fields'][0]['mesh']:
        try:result.append(owner_point(report,element,p))
        except ValueError:pass
    if not result:raise ValueError('physical point has no native owner')
    return result


def field(report,state,block,owner):
    element=next((v for v in report['fields'][block]['mesh'] if v['geometry']==owner['geometry']),None)
    if element is None:raise ValueError('field absent in this region (no inner epsilon)')
    # Revalidate the owner even when called with externally supplied provenance.
    p=physical(report,owner['chart']);o=owner_point(report,element,p);x,_,rho=o['chart'];lo,hi=element['lower'],element['upper']
    px,dx=basis(element['degree'],[2*(x-lo[0])/(hi[0]-lo[0])-1]);pr,dr=basis(element['degree'],[2*(rho-lo[2])/(hi[2]-lo[2])-1])
    c=np.asarray(state)[np.asarray(element['ids'],int)+int(report['offsets'][block])]*element['sign']
    phi=np.outer(pr[0],px[0]).ravel();gx=np.outer(pr[0],dx[0]).ravel()*2/(hi[0]-lo[0])
    J=1. if o['region']=='inner' else mapping(report)[3]
    gz=np.outer(dr[0],px[0]).ravel()*2/(hi[2]-lo[2])/J
    return float(phi@c),np.array([gx@c,0.,gz@c])


def matching(report,owner):
    H=mapping(report)[0]
    e=next(v for v in report['fields'][0]['mesh'] if v['lower'][2]==H and
        v['lower'][0]==owner['lower'][0] and v['upper'][0]==owner['upper'][0])
    return owner_point(report,e,physical(report,[owner['chart'][0],0,H]))


def evaluate(report,state,p,owner):
    # Full3D tensors use physical Cartesian x,y,z axes despite planar means.
    U=np.zeros(3);G=np.zeros((3,3))
    for block,axis in ((0,0),(1,2)):U[axis],G[axis]=field(report,state,block,owner)
    k,gk=field(report,state,2,owner);pi,gpi=field(report,state,4 if owner['region']=='inner' else 5,owner)
    H=mapping(report)[0];K=field(report,state,2,matching(report,owner))[0]
    v=dict(physical=list(map(float,p)),owner=owner,U=U.tolist(),G=G.tolist(),k=k,pi=pi,mechanical_pressure=pi-2*k/3,
        grad_k=gk.tolist(),sampled_K=K,raw_divergence=float(np.trace(G)))
    if not k>0 or not K>0:raise ValueError('nonpositive scalar at expanded point')
    if owner['region']=='inner':
        w=wall(K,H,report['parameters'][2]);nt=w['a']*max(owner['chart'][2]-w['zv'],0.)
        R=2*k/3*np.eye(3)-nt*(G+G.T);T=(NU+nt)*(G+G.T)-pi*np.eye(3)
        v.update(epsilon=None,epsilon_scope='not_evaluated:no_inner_epsilon_field',tau=None,principal=None,N=None,F_N=None,root_residual=None,
            nu_eff=nt,nu_d=NU+nt,cutoff=w['rc'],cutoff_margin=H-w['rc'],destruction=w['A']/max(owner['chart'][2],w['rc']))
    else:
        e,ge=field(report,state,3,owner)
        if not e>0:raise ValueError('nonpositive outer epsilon')
        R,nt,N=constitutive(k,e,G);T=-R+(2*k/3-pi)*np.eye(3)
        B=k/e*(G-np.trace(G)*np.eye(3)/3);S=(B+B.T)/2;O=(B-B.T)/2;s=float(np.trace(S@S));w=float(np.trace(O@O))
        F=(N-1.8)*(N*N-2*w)-2.7*s*N;FN=1.8*(N*N-2*w)/N+2*N*(N-1.8)
        v.update(epsilon=e,epsilon_scope='outer_physical_field',tau=k/e,principal=principal(k,e,G),N=float(N),F_N=float(FN),
            root_residual=float(F/(1+N**3+abs(2.7*s*N)+abs(w*N))),nu_eff=float(nt),nu_d=CD*k*k/e,grad_epsilon=ge.tolist())
    v.update(R=R.tolist(),T=T.tolist(),covariance=np.linalg.eigvalsh(R).tolist(),Praw=float(-np.sum(R*G)),
        Pdev=float(-np.sum(R*G)+2*k*np.trace(G)/3))
    return v


def inventory(reports,previous=None):
    H,d0,top,J=mapping(reports[0]);points={}
    def add(p,tag,near=False):
        p=list(map(float,p));key=','.join(v.hex() for v in p)
        if not (0<=p[0]<=8 and d0<=p[2]<=top):return
        item=points.setdefault(key,dict(physical=p,tags=[]))
        if tag not in item['tags']:item['tags'].append(tag)
        if near:
            for axis in (0,2):
                for sign in (-1,1):
                    q=p.copy();q[axis]+=sign*OFFSET;add(q,f'physical_offset:{axis}:{sign}:{key}')
    if previous:
        for point in previous['points']:
            for tag in point['tags']:add(point['physical'],tag)
    for x in (0.,8/3,3.,4.,16/3,6.,7.5,7.9,7.99,8.):
        add([x,0,d0+H],'declared_H',True)
        for z in (.6720111978440669,.8818157002137355,.9080394908593612,1.1947077409162705):add([x,0,z],'declared_limiting_neighborhood',True)
    for r in reports:
        packets=r['wj_audit_packets']+[e['wj_packet'] for e in r.get('physical_range_events',[]) if e.get('wj_packet')]
        for packet in packets:
            for name,w in packet['principal_witnesses'].items():add(w['position'],f'native_witness:{packet["state_sha256"]}:{name}',True)
        for row in r['open_boundary']['interface']:add(row['physical'],'existing_interface')
        for row in r['open_boundary']['boundary_rows']:
            if row['role']=='outlet':
                add(row['physical'],'existing_outlet')
                q=row['physical'].copy();q[0]-=OFFSET;add(q,'outlet_physical_inward')
    result=dict(schema='WJH4PhysicalInventoryV1',points=[dict(id=i,**v) for i,v in enumerate(points.values())],
        offset=OFFSET,scope='finite physical union; every legitimate owner retained')
    result['sha256']=digest(result);return result


def scan(report,state,inv):
    rows=[];failures=[];minimum=None
    for point in inv['points']:
        for o in owners(report,point['physical']):
            try:
                v=evaluate(report,state,point['physical'],o);v['point_id']=point['id'];rows.append(v)
                if o['region']=='outer':
                    m=min(v['principal']['in_plane'],v['principal']['transverse'])
                    if minimum is None or m<minimum['margin']:minimum=dict(margin=m,row=v)
                    if m<=2e-11:failures.append(dict(category='principal_negative' if m< -2e-11 else 'principal_sign_unresolved',row=v))
                if min(v['covariance'])< -2e-13:failures.append(dict(category='covariance_negative',row=v))
            except (ValueError,FloatingPointError,np.linalg.LinAlgError) as e:failures.append(dict(category='invalid_state',point=point,owner=o,reason=str(e)))
    return dict(status='PASS' if not failures else 'FAIL' if any(v['category']!='principal_sign_unresolved' for v in failures) else 'UNRESOLVED',
        inventory_sha256=inv['sha256'],state_sha256=vector_hash(state),minimum=minimum,failures=failures,rows=rows)


def aligned(step,dt,recorded):
    target=Fraction(str(dt))*int(step);value=float(target);tol=128*EPS*max(1,abs(value))
    if abs(recorded-value)>tol:raise ValueError('recorded time inconsistent with integer steps')
    return target


def extend_scan(report,state,old,old_inventory,new_inventory):
    """Evaluate only new physical points; retain every old owner/result."""
    n=len(old_inventory['points'])
    if [v['physical']for v in old_inventory['points']]!=[v['physical']for v in new_inventory['points'][:n]]:
        raise ValueError('inventory is not a physical-coordinate prefix extension')
    if old['state_sha256']!=vector_hash(state):raise ValueError('stale state for incremental scan')
    addition=scan(report,state,dict(new_inventory,points=new_inventory['points'][n:]))
    minima=[v for v in (old['minimum'],addition['minimum']) if v is not None]
    failures=old['failures']+addition['failures']
    result=dict(old,rows=old['rows']+addition['rows'],failures=failures,
        minimum=min(minima,key=lambda v:v['margin']) if minima else None,inventory_sha256=new_inventory['sha256'],
        status='PASS' if not failures else 'FAIL' if any(v['category']!='principal_sign_unresolved' for v in failures) else 'UNRESOLVED')
    result.pop('compiled_in_plane_replay',None)
    return result


def pair_rows(a,b,step_a,dt_a,time_a,step_b,dt_b,time_b):
    ta=aligned(step_a,dt_a,time_a);tb=aligned(step_b,dt_b,time_b)
    if ta!=tb:raise ValueError('different intended times; no interpolation')
    if a['inventory_sha256']!=b['inventory_sha256']:raise ValueError('different comparison inventories')
    def group(scan):
        result={}
        for row in scan['rows']:result.setdefault((row['point_id'],row['owner']['region']),[]).append(row)
        return result
    aa,bb=group(a),group(b)
    if aa.keys()!=bb.keys():raise ValueError('missing physical point or region')
    result=[]
    for key in sorted(aa):
        # Explicit one-to-many owner relations; never select a lowest ID or
        # require matching element bounds on independently refined meshes.
        for x,y in itertools.product(aa[key],bb[key]):
            # If BOTH meshes have a face here, compare the same physical side.
            # Opposite-side pairs would confuse a gradient jump with a time/grid
            # difference. Face-versus-interior remains explicitly one-to-many.
            if any(x['owner'][axis]!='interior' and y['owner'][axis]!='interior' and
                   x['owner'][axis]!=y['owner'][axis] for axis in ('x_side','normal_side')):continue
            result.append(dict(point_id=key[0],region=key[1],physical=x['physical'],owner_a=x['owner'],owner_b=y['owner'],
                intended_time=float(ta),recorded_time_a=time_a,recorded_time_b=time_b,
                margin_a=None if x['principal'] is None else x['principal']['in_plane'],
                margin_b=None if y['principal'] is None else y['principal']['in_plane'],
                delta_G=(np.asarray(y['G'])-x['G']).tolist(),delta_tau=None if x['tau'] is None else y['tau']-x['tau'],
                delta_k=y['k']-x['k'],delta_epsilon=None if x['epsilon'] is None else y['epsilon']-x['epsilon'],delta_pi=y['pi']-x['pi']))
    return result
