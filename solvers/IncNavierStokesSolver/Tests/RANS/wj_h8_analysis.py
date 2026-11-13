#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H8 local owners, passive pressure mass algebra and three-level metrics."""
import math
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from wj_h4_points import owners,physical,vector_hash,aligned,digest
from wj_h7_analysis import Rows,principal_batch,context,EXACT_MARGIN
from wj_h5_neutral import CD,KAPPA
from validate_cls_epawf_m3 import basis

FACE=.25607517888672315

def add_local(report,inventory,c):
    """Frozen deterministic physical rule, no choice based on its results."""
    inv={tuple(float(v).hex()for v in p['physical']):dict(physical=p['physical'],tags=list(p['tags']))for p in inventory['points']}
    local=set()
    def add(p,tag):
        if not(0<=p[0]<=8 and c.d0<=p[2]<=c.top):return
        key=tuple(float(v).hex()for v in p);inv.setdefault(key,dict(physical=list(map(float,p)),tags=[]))['tags'].append(tag);local.add(key)
    intervals=[e for e in report['fields'][3]['mesh']if abs(physical(report,e['lower'])[2]-FACE)<1e-13 or abs(physical(report,e['upper'])[2]-FACE)<1e-13]
    if len(intervals)!=2:raise ValueError('two declared adjacent outer elements required')
    heights=[c.d0+c.H,FACE,FACE-1e-6,FACE+1e-6]
    for e in intervals:heights.append(physical(report,(np.array(e['lower'])+e['upper'])/2)[2])
    for x in (0.,4.,7.9,7.99,8.-1e-6,8.):
        for z in heights:add([x,0,z],'H8 fixed face/H/interior/outlet')
    g=np.array([-1.,-math.sqrt(3/7),0.,math.sqrt(3/7),1.])
    for e in intervals:
        lo=np.array(e['lower']);hi=np.array(e['upper'])
        for x in (lo[0]+hi[0])/2+(hi[0]-lo[0])/2*g:
            for r in (lo[2]+hi[2])/2+(hi[2]-lo[2])/2*g:add(physical(report,[x,0,r]),'H8 adjacent-element5x5Lobatto')
    result=dict(points=list(inv.values()),local_keys=[list(k)for k in sorted(local)],rule='H7 complete union + frozen H8 physical rule/all legitimate owners',adjacent_geometry=[e['geometry']for e in intervals])
    result['sha256']=digest(result);return result


def extend_witnesses(report,inv,c):
    table={tuple(float(v).hex()for v in p['physical']):p for p in inv['points']}
    for packet in report['wj_audit_packets']+[e['wj_packet']for e in report.get('physical_range_events',[])if e.get('wj_packet')]:
        for w in packet['principal_witnesses'].values():
            p=w['position'];ps=[p]
            for axis in (0,2):
                for sign in (-1,1):q=p.copy();q[axis]+=sign*1e-6;ps.append(q)
            for q in ps:
                if 0<=q[0]<=8 and c.d0<=q[2]<=c.top:table.setdefault(tuple(float(v).hex()for v in q),dict(physical=q,tags=['H8 new native witness/neighborhood']))
    out=dict(inv,points=list(table.values()));out.pop('sha256',None);out['sha256']=digest(out);return out


class PressureProjection:
    """Actual continuous regional Q1 and broken-Q1 physical L2 projections.

    Solves only diagnostic mass equations, never a governing pressure equation.
    Input rows already contain the physical quadrature weights/Jacobians.
    """
    def __init__(self,rows):
        self.rows=rows;self.regions={};self.factorizations=0
        for block,mask in ((4,~rows.outer),(5,rows.outer)):
            active=np.flatnonzero(mask & (rows.weights>0));w=rows.weights[active];mesh={e['geometry']:e for e in rows.report['fields'][block]['mesh']};n=rows.report['fields'][block]['size']
            B=np.zeros((len(active),n));geometry=[]
            for j,i in enumerate(active):
                p,o=rows.rows[i];e=mesh[o['geometry']];lo,hi=e['lower'],e['upper'];x,_,r=o['chart']
                phi=np.outer(basis(1,[2*(r-lo[2])/(hi[2]-lo[2])-1])[0][0],basis(1,[2*(x-lo[0])/(hi[0]-lo[0])-1])[0][0]).ravel()*e['sign']
                B[j,e['ids']]=phi;geometry.append(e['geometry'])
            M=B.T@(w[:,None]*B);factor=cho_factor(M);self.factorizations+=1;elements=[]
            for g in sorted(set(geometry)):
                ii=np.flatnonzero(np.array(geometry)==g);ids=np.asarray(mesh[g]['ids']);local=B[np.ix_(ii,ids)];mf=cho_factor(local.T@(w[ii,None]*local));self.factorizations+=1
                elements.append((g,ii,ids,local,mf))
            self.regions[str(block)]=(active,w,B,factor,elements)
    def measure(self,data):
        out={};solves=0
        for region,(idx,w,B,factor,elements)in self.regions.items():
            div=data['divergence'][idx];rhs=B.T@(w*div);coeff=cho_solve(factor,rhs);solves+=1;norm=float(w@(div*div));proj=float(rhs@coeff);comp=float(w@((div-B@coeff)**2));broken=0.;broken_comp=0.;detail=[]
            for g,ii,ids,local,mf in elements:
                b=local.T@(w[ii]*div[ii]);v=cho_solve(mf,b);solves+=1;bp=float(b@v);dn=float(w[ii]@(div[ii]**2));broken+=bp;broken_comp+=float(w[ii]@((div[ii]-local@v)**2))
                detail.append(dict(geometry=g,divergence_L2_squared=dn,broken_Q1_projection_squared=bp,raw_complement_squared=dn-bp,mode_coefficients=v.tolist()))
            out[region]=dict(divergence_L2=math.sqrt(max(0,norm)),continuous_Q1_projection_L2=math.sqrt(max(0,proj)),continuous_Q1_complement_L2=math.sqrt(comp),
                broken_Q1_projection_L2=math.sqrt(max(0,broken)),broken_Q1_complement_L2=math.sqrt(broken_comp),raw_continuous_complement_squared=norm-proj,
                pressure_test_moments=rhs.tolist(),maximum_pressure_test_moment=float(max(abs(rhs))),integrated_net_flux=float(sum(rhs)),elements=detail,
                measure=float(sum(w)),identity_roundoff=abs(norm-proj-comp),broken_identity_roundoff=abs(norm-broken-broken_comp),gauge_contribution=0.,scope='outlet fixes level; no gauge unknown; continuous-Q1 complement includes continuity restrictions beyond broken polynomial modes')
        return dict(regions=out,mass_solves=solves,governing_solves=0,scope='passive physical pressure-space mass algebra, state/history unchanged')


def local_rows(cache,data,indices):
    oo=[i for i in indices if cache.outer[i]];details={}
    if oo:
        pv=principal_batch(data['k'][oo],data['epsilon'][oo],data['G'][oo],details=True)
        for j,i in enumerate(oo):details[i]=(pv[3][j],pv[4][j],pv[5][j],pv[6]['transverse_xi'][j])
    out=[]
    for i in indices:
        v=cache.witness(data,i);v['epsilon_scope']='outer_physical'if cache.outer[i]else 'absent:no_inner_epsilon'
        if cache.outer[i]:
            m,t,angle,tx=details[i];xi=np.array([math.cos(angle),0.,math.sin(angle)]);a=np.array([-xi[2],0.,xi[0]])
            G=data['G'][i];N=data['root'][i];B=data['k'][i]/data['epsilon'][i]*(G-np.trace(G)*np.eye(3)/3);S=(B+B.T)/2;O=(B-B.T)/2;s=float(np.trace(S@S));w=float(np.trace(O@O))
            v.update(margin=float(m),transverse=float(t),angle=float(angle),xi=xi.tolist(),polarization=a.tolist(),a_dot_xi=float(a@xi),transverse_xi=tx.tolist(),F_N=float(1.8*(N*N-2*w)/N+2*N*(N-1.8)),
                scaled_root_residual=float(((N-1.8)*(N*N-2*w)-2.7*s*N)/(1+N**3+abs(2.7*s*N)+abs(w*N))),Pdev=float(v['production']+2*v['k']*v['divergence']/3))
        else:
            for k in ('epsilon','grad_epsilon','margin','transverse','angle','root'):v[k]=None
        out.append(v)
    return out


def ratio(d1,d2,floor):
    if not all(math.isfinite(x)for x in (d1,d2,floor))or min(d1,d2,floor)<0:raise ValueError('finite nonnegative comparison required')
    if d1<=floor:return dict(D1=d1,D2=d2,floor=floor,rho=None,status='UNRESOLVED',reason='denominator at comparison floor')
    r=d2/d1;return dict(D1=d1,D2=d2,floor=floor,rho=r,status='CONTRACTS'if r<1 else 'DOES_NOT_CONTRACT')


def strict_context(report,step):
    c=context(report,step)
    if c['role']!='accepted' or c['step']!=step:raise ValueError('strict view requires genuine accepted context')
    return c
