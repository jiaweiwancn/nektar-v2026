#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Batched actual-phase weak loads with recovered native projected tests.

No native Project/Lift, pressure/momentum assembly or governing solve. Inner
quadrature is rebuilt at each state's endpoint-K cutoff, not borrowed from A.
"""
import xml.etree.ElementTree as ET
import numpy as np
from validate_cls_epawf_m3 import basis
from wj_h4_points import mapping
from wj_h6_projection import expression
from wj_h7_analysis import wall_batch
from wj_h3_analysis import algebra
from wj_h5_neutral import CD


def shapes(el,x,r,J):
    lo,hi=el['lower'],el['upper'];p=el['degree']
    px,dx=basis(p,2*(x-lo[0])/(hi[0]-lo[0])-1);pr,dr=basis(p,2*(r-lo[2])/(hi[2]-lo[2])-1)
    make=lambda a,b:np.einsum('ni,nj->nij',a,b).reshape(len(x),-1)*el['sign']
    return make(pr,px),make(pr,dx)*2/(hi[0]-lo[0]),make(dr,px)*2/(hi[2]-lo[2])/J


def reaction(role,T,n,given):
    """Current-state reactions; outlet and top tangential data alone are fixed."""
    actual=np.einsum('nij,j->ni',T,n)
    if role=='outlet':return np.broadcast_to(given,actual.shape).copy()
    if role=='top':return given+(actual@n-given@n)[:,None]*n
    return actual


def pressure_terms(pi,grad,phi,grad_pi_H,rough):
    volume=-pi[None,:,None]*grad
    forcing=np.zeros_like(volume);forcing[:,:,0]=-phi*grad_pi_H[None,:]*rough[None,:]
    return volume,forcing


class ProjectedTests:
    def __init__(self,report,packet):
        self.report=report;self.H,self.d0,self.top,self.J=mapping(report);self.coeff={};self.max_value_error=0.;self.max_gradient_error=0.;self.fits=0
        groups={}
        for row in packet['samples']:groups.setdefault(row['geometry'],[]).append(row)
        for geom,rows in groups.items():
            x=np.array([v['chart'][0]for v in rows]);r=np.array([v['chart'][2]for v in rows]);inner=rows[0]['inner'];J=1 if inner else self.J
            for typ,b in [('v',0),('p',4 if inner else 5)]+([]if inner else [('e',3)]):
                el=next(e for e in report['fields'][b]['mesh']if e['geometry']==geom);B,Bx,Bz=shapes(el,x,r,J)
                y=np.array([[v['modes'][m][typ]for m in range(3)]for v in rows]);co=np.linalg.lstsq(B,y,rcond=None)[0];self.fits+=1
                self.max_value_error=max(self.max_value_error,float(np.max(abs(B@co-y))))
                if typ!='p':
                    g=np.array([[v['modes'][m]['grad_'+typ]for m in range(3)]for v in rows])
                    self.max_gradient_error=max(self.max_gradient_error,float(np.max(abs(Bx@co-g[:,:,0]))),float(np.max(abs(Bz@co-g[:,:,2]))))
                self.coeff[(geom,typ)]=(el,co)
        assert self.max_value_error<2e-11 and self.max_gradient_error<2e-10
        self.boundary_error=self.validate_boundary(packet)
    def values(self,geom,typ,x,r):
        el,c=self.coeff[(geom,typ)];J=1 if el['upper'][2]<=self.H else self.J;B,Bx,Bz=shapes(el,x,r,J)
        grad=np.zeros((3,len(x),3));grad[:,:,0]=(Bx@c).T;grad[:,:,2]=(Bz@c).T
        return (B@c).T,grad
    def validate_boundary(self,packet):
        error=0.;groups={}
        for row in packet['boundary_mode_samples']:groups.setdefault(row['geometry'],[]).append(row)
        for g,rows in groups.items():
            x=np.array([v['chart'][0]for v in rows]);r=np.array([v['chart'][2]for v in rows]);v,_=self.values(g,'v',x,r)
            for b in (0,1,2):
                expected=np.array([[z['modes'][m][b]for z in rows]for m in range(3)]);error=max(error,float(np.max(abs(v-expected))))
            if(g,'e')in self.coeff:
                v,_=self.values(g,'e',x,r);expected=np.array([[z['modes'][m][3]for z in rows]for m in range(3)]);error=max(error,float(np.max(abs(v-expected))))
        assert error<2e-11;return error
    def export(self):
        return dict(fits=self.fits,value_max_error=self.max_value_error,gradient_max_error=self.max_gradient_error,boundary_max_error=self.boundary_error,
            coefficients=[dict(geometry=g,type=t,coefficients=c.tolist())for(g,t),(e,c)in self.coeff.items()],scope='reconstructed from actual exported native polynomial values/gradients/traces, not analytic hats')


class WeakAudit:
    def __init__(self,report,xml,tests):
        self.r=report;self.tests=tests;self.H,self.d0,self.top,self.J=mapping(report);self.nu=report['parameters'][0];self.h=report['parameters'][2]
        self.el=[{e['geometry']:e for e in f['mesh']}for f in report['fields']];self.fixed={};self.rules={};self.evaluations=0
        for f in ET.parse(xml).findall('./CONDITIONS/FUNCTION'):
            if f.get('NAME')in ('OutletTractionInner','OutletTractionOuter','TopTractionOuter','TopDiffusiveFluxOuter'):
                for e in f.findall('E'):self.fixed[(f.get('NAME'),e.get('VAR'))]=expression(e.get('VALUE'))
        self.chart_top=max(e['upper'][2]for e in report['fields'][0]['mesh']);self.match=next(e['geometry']for e in report['fields'][3]['mesh']if e['lower'][2]==self.H)
    def field(self,state,b,g,x,r):
        el=self.el[b][g];J=1 if el['upper'][2]<=self.H else self.J;B,Bx,Bz=shapes(el,x,r,J);c=np.asarray(state)[np.asarray(el['ids'])+self.r['offsets'][b]]
        grad=np.zeros((len(x),3));grad[:,0]=Bx@c;grad[:,2]=Bz@c;return B@c,grad
    def fields(self,s,g,x,r):
        n=len(x);inner=self.el[0][g]['upper'][2]<=self.H;U=np.zeros((n,3));G=np.zeros((n,3,3))
        U[:,0],G[:,0]=self.field(s,0,g,x,r);U[:,2],G[:,2]=self.field(s,1,g,x,r);k,gk=self.field(s,2,g,x,r);pi,gpi=self.field(s,4 if inner else 5,g,x,r)
        h=np.full(n,self.H);K,gkh=self.field(s,2,self.match,x,h);uh,gh=self.field(s,0,self.match,x,h);_,gph=self.field(s,5,self.match,x,h)
        a,A,zv,lam,rc=wall_batch(K,self.H,self.h,self.nu);M=G+G.swapaxes(-1,-2)
        if inner:
            nt=a*np.maximum(r-zv,0);R=2*k[:,None,None]*np.eye(3)/3-nt[:,None,None]*M;dk=self.nu+nt;de=np.zeros(n);e=A/np.maximum(r,rc);ge=np.zeros((n,3));T=dk[:,None,None]*M-pi[:,None,None]*np.eye(3)
        else:
            e,ge=self.field(s,3,g,x,r);R,nt,N=algebra(k,e,G);dk=CD*k*k/e;de=dk/1.3;T=-R+(2*k/3-pi)[:,None,None]*np.eye(3)
        P=-np.sum(R*G,axis=(1,2))
        return dict(U=U,G=G,k=k,gk=gk,epsilon=e,ge=ge,pi=pi,R=R,T=T,P=P,dk=dk,de=de,nt=nt,K=K,rc=rc,Lambda=lam,
            grad_pi_H=gph[:,0],force_H=gph[:,0]-2*gkh[:,0]/3+2*uh*gh[:,0])
    def quadrature(self,s,q):
        z,w=np.polynomial.legendre.leggauss(q)
        for g,el in self.el[0].items():
            lo,hi=el['lower'],el['upper'];xx=(lo[0]+hi[0])/2+(hi[0]-lo[0])/2*z;wx=w*(hi[0]-lo[0])/2;inner=hi[2]<=self.H
            cuts=sorted(set([lo[2],hi[2]]+([self.h]if lo[2]<self.h<hi[2]else [])))
            K,_=self.field(s,2,self.match,xx,np.full(q,self.H));rc=wall_batch(K,self.H,self.h,self.nu)[4]if inner else None
            for a,b in zip(cuts[:-1],cuts[1:]):
                mid=np.clip(rc,a,b)if inner else None
                segments=[(np.full(q,a),mid),(mid,np.full(q,b))]if inner else [(np.full(q,a),np.full(q,b))]
                for aa,bb in segments:
                    good=bb>aa
                    if not np.any(good):continue
                    x=np.broadcast_to(xx[good,None],(sum(good),q)).ravel();r=((aa[good]+bb[good])[:,None]/2+(bb[good]-aa[good])[:,None]/2*z).ravel()
                    weight=(wx[good,None]*(bb[good]-aa[good])[:,None]/2*w).ravel()*(1 if inner else self.J)
                    yield g,x,r,weight,inner
    def exterior(self,q):
        z,w=np.polynomial.legendre.leggauss(q)
        for g,el in self.el[0].items():
            lo,hi=el['lower'],el['upper'];inner=hi[2]<=self.H;J=1 if inner else self.J
            for role,x,n in (('inlet',lo[0],[-1.,0,0]),('outlet',hi[0],[1.,0,0])):
                rr=(lo[2]+hi[2])/2+(hi[2]-lo[2])/2*z
                yield role,g,np.full(q,x),rr,w*(hi[2]-lo[2])/2*J,np.array(n),inner
            for role,r,n in [('ground',lo[2],[0.,0,-1.]) ]if lo[2]==0 else [('top',hi[2],[0.,0,1.])]if hi[2]==self.chart_top else []:
                xx=(lo[0]+hi[0])/2+(hi[0]-lo[0])/2*z
                yield role,g,xx,np.full(q,r),w*(hi[0]-lo[0])/2,np.array(n),inner
    def given(self,fam,var,x,r):
        f=self.fixed.get((fam,var));physical=self.d0+np.where(r<=self.H,r,self.H+(r-self.H)*self.J)
        return np.zeros_like(x)if f is None else f(x,physical)
    def evaluate(self,s,history,alpha,q):
        self.evaluations+=1;parts={};split={};roles={};cutoffs=[]
        def add(table,name,region,block,values,w):table.setdefault(name,np.zeros((2,3,6)))[region,:,block]+=values@w
        for g,x,r,w,inner in self.quadrature(s,q):
            side=int(not inner);d=self.fields(s,g,x,r);phi,grad=self.tests.values(g,'v',x,r);pp,_=self.tests.values(g,'p',x,r);u=d['U'];k=d['k'];e=d['epsilon'];div=np.trace(d['G'],axis1=-2,axis2=-1)
            cutoffs.append([float(min(d['rc'])),float(max(d['rc']))])
            hist=[self.field(history,b,g,x,r)[0]for b in ((0,1,2)if inner else (0,1,2,3))]
            adv=np.einsum('mni,ni->mn',grad,u)
            for b,axis in ((0,0),(1,2)):
                add(parts,'time_mass',side,b,alpha*phi*(u[:,axis]-hist[b]),w);add(parts,'conservative_advection',side,b,-adv*u[:,axis],w)
                stress=np.einsum('mni,ni->mn',grad,d['T'][:,axis]);add(parts,'inner_momentum_stress'if inner else 'outer_momentum_stress',side,b,stress,w)
                pressure=-grad[:,:,axis]*d['pi'];add(split,'reduced_pressure_volume',side,b,pressure,w);add(split,'nonpressure_stress_volume',side,b,stress-pressure,w)
                if inner:
                    add(parts,'roughness_force',side,b,-phi*d['force_H']*(r<=self.h)if axis==0 else np.zeros_like(phi),w)
                    add(split,'roughness_pressure',side,b,-phi*d['grad_pi_H']*(r<=self.h)if axis==0 else np.zeros_like(phi),w)
            add(parts,'time_mass',side,2,alpha*phi*(k-hist[2]),w);add(parts,'conservative_advection',side,2,-adv*k,w);add(parts,'k_diffusion',side,2,d['dk']*np.einsum('mni,ni->mn',grad,d['gk']),w)
            add(parts,'inner_k_production'if inner else 'outer_k_production',side,2,-phi*d['P'],w);add(parts,'inner_k_destruction'if inner else 'outer_k_destruction',side,2,phi*e,w)
            add(parts,'continuity',side,4 if inner else 5,pp*div,w)
            if not inner:
                ep,eg=self.tests.values(g,'e',x,r);ea=np.einsum('mni,ni->mn',eg,u)
                for name,val in [('time_mass',alpha*ep*(e-hist[3])),('conservative_advection',-ea*e),('epsilon_diffusion',d['de']*np.einsum('mni,ni->mn',eg,d['ge'])),('epsilon_production',-ep*1.44*d['P']*e/k),('epsilon_destruction',ep*1.82*e*e/k)]:add(parts,name,side,3,val,w)
        for role,g,x,r,w,n,inner in self.exterior(q):
            side=int(not inner);d=self.fields(s,g,x,r);given=np.zeros((len(x),3))
            if role=='outlet':
                fam='OutletTractionInner'if inner else 'OutletTractionOuter';given[:,0]=self.given(fam,'u',x,r);given[:,2]=self.given(fam,'v',x,r)
            elif role=='top':given[:,0]=self.given('TopTractionOuter','u',x,r)
            applied=reaction(role,d['T'],n,given);un=d['U']@n;phi,_=self.tests.values(g,'v',x,r)
            for b,axis in ((0,0),(1,2)):
                flux=d['U'][:,axis]*un-applied[:,axis];add(parts,'exterior_boundary',side,b,phi*flux,w);add(roles,role,side,b,phi*flux,w)
                if role!='outlet':add(split,'pressure_reaction_'+role,side,b,phi*d['pi']*n[axis],w)
            for b,name,gradname,diffname in [(2,'k','gk','dk')]+([]if inner else [(3,'epsilon','ge','de')]):
                ph=phi if b==2 else self.tests.values(g,'e',x,r)[0]
                j=-d[diffname]*(d[gradname]@n)if role=='inlet'else self.given('TopDiffusiveFluxOuter','k'if b==2 else 'epsilon',x,r)if role=='top'else np.zeros(len(x))
                val=ph*(d[name]*un+j);add(parts,'exterior_boundary',side,b,val,w);add(roles,role,side,b,val,w)
        total=sum(parts.values());return dict(order=q,components={k:v.sum(axis=0).tolist()for k,v in parts.items()},regional=total.tolist(),regional_components={k:v.tolist()for k,v in parts.items()},
            splits={k:v.tolist()for k,v in split.items()},exterior_roles={k:v.tolist()for k,v in roles.items()},cutoff_range=[min(a for a,b in cutoffs),max(b for a,b in cutoffs)],scope='same history/alpha/time, actual projected tests, all reactions recomputed; uncancelled weak loads, not frozen free-row residual')
    def pressure_action(self,direction,q):
        # No constitutive or scalar evaluation of this pressure-only field.
        pieces={}
        g,wg=np.polynomial.legendre.leggauss(q)
        def add(name,side,b,v,w):pieces.setdefault(name,np.zeros((2,3,6)))[side,:,b]+=v@w
        for geom,el in self.el[0].items():
            lo,hi=el['lower'],el['upper'];inner=hi[2]<=self.H;side=int(not inner);cuts=sorted(set([lo[2],hi[2]]+([self.h]if lo[2]<self.h<hi[2]else [])))
            for a,b in zip(cuts[:-1],cuts[1:]):
                xx,rr=np.meshgrid((lo[0]+hi[0])/2+(hi[0]-lo[0])/2*g,(a+b)/2+(b-a)/2*g);x=xx.ravel();r=rr.ravel();w=(np.outer(wg,wg)*(hi[0]-lo[0])*(b-a)/4*(1 if inner else self.J)).ravel()
                pi,_=self.field(direction,4 if inner else 5,geom,x,r);phi,grad=self.tests.values(geom,'v',x,r);gp=self.field(direction,5,self.match,x,np.full(len(x),self.H))[1][:,0]
                vol,rough=pressure_terms(pi,grad,phi,gp,(r<=self.h)if inner else np.zeros(len(x),bool))
                for bl,axis in ((0,0),(1,2)):add('volume',side,bl,vol[:,:,axis],w);add('roughness',side,bl,rough[:,:,axis],w)
        for role,geom,x,r,w,n,inner in self.exterior(q):
            if role=='outlet':continue
            pi,_=self.field(direction,4 if inner else 5,geom,x,r);phi,_=self.tests.values(geom,'v',x,r)
            for b,axis in ((0,0),(1,2)):add('reaction_'+role,int(not inner),b,phi*pi*n[axis],w)
        return dict(pieces={k:v.tolist()for k,v in pieces.items()},regional=sum(pieces.values()).tolist(),scope='unrestricted projected-test weak pressure functional; NOT reduced free-row Gp')


def compare(a,b):
    names=set(a['components'])|set(b['components']);conditions={};errors={}
    for k in names:
        x=np.array(a['components'].get(k,np.zeros((3,6))));y=np.array(b['components'].get(k,np.zeros((3,6))));err=abs(x-y);conditions[k]=bool(np.all(np.isfinite(err))and np.all(err<=2e-9+2e-8*abs(y)));errors[k]=float(np.max(err))
    x=np.asarray(a['regional']);y=np.asarray(b['regional']);conditions['regional_sum']=bool(np.all(abs(x-y)<=4e-9+2e-8*abs(y)))
    return dict(status='PASS'if all(conditions.values())else 'UNRESOLVED',conditions=conditions,component_max_errors=errors,regional_max_error=float(np.max(abs(x-y))))
