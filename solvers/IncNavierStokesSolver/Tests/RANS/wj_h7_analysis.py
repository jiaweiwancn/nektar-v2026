#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H7 immutable row cache and independent per-state full-gradient diagnostics."""
import math
import numpy as np
from validate_cls_epawf_m3 import basis
from wj_h4_points import owners,mapping,physical,matching,vector_hash,aligned
from wj_h3_analysis import algebra,state_context
from wj_h5_neutral import QE,KAPPA,CD,NE

EXACT_MARGIN=.004137527846908
FIELDS=('U','G','k','epsilon','R','pi','mechanical_pressure')

def context(r,step=None,rejected=False):
    if rejected:
        if not r.get('physical_range_events'):raise ValueError('no retained rejected-stage context')
        s,h,a,t,p=state_context(r,rejected=True);role='rejected'
    else:
        step=r['steps']if step is None else step
        first=r['steps']+1-len(r['accepted_states'])
        if not first<=step<=r['steps']:raise ValueError('accepted step outside genuine producer')
        s,h,a,t,p=state_context(r,step);role='accepted'
    if vector_hash(s)!=p['state_sha256']:raise ValueError('state/audit binding differs')
    if a and (h is None or vector_hash(h)!=p['stage_history_sha256']):raise ValueError('callback binding differs')
    if a and a!=1/r['dt']:raise ValueError('callback alpha/dt differs')
    n=round(t/r['dt']);aligned(n,r['dt'],t)
    return dict(state=s,history=h,alpha=a,time=t,step=n,role=role,packet=p)


def principal_batch(k,e,G,*,details=False):
    """Same complex implicit-root continuation and planar quartic as H3.

    Only independent points are batched; no direction or I4 term is omitted.
    Geometry caches do not contain closure values or state-dependent admission.
    """
    R,nt,N=algebra(k,e,G);L=np.empty((len(k),9,9))
    for j,h in enumerate(np.eye(9)):
        L[:,:,j]=-algebra(k,e,G.astype(complex)+1e-27j*h.reshape(3,3))[0].imag.reshape(-1,9)/1e-27/nt[:,None]
    A=np.zeros((3,3));A[2,0]=1;B=np.diag([-1.,0,1]);C=np.zeros((3,3));C[0,2]=-1
    V=np.array([A.ravel(),B.ravel(),C.ravel()]);Q=np.einsum('ai,nij,bj->nab',V,(L+L.transpose(0,2,1))/2,V)
    v=np.zeros((2,9));v[0,3]=1;v[1,5]=1
    transverse_matrix=np.einsum('ai,nij,bj->nab',v,(L+L.transpose(0,2,1))/2,v)
    if details:
        et,vt=np.linalg.eigh(transverse_matrix);transverse=et[:,0]
    else:transverse=np.linalg.eigvalsh(transverse_matrix)[:,0]
    result=np.empty(len(k));angles=np.empty(len(k))
    for i,q in enumerate(Q):
        a,b,c,d,f=q[0,0],2*q[0,1],q[1,1]+2*q[0,2],2*q[1,2],q[2,2]
        co=np.array([-d,4*f-2*c,3*d-3*b,2*c-4*a,b]);roots=np.roots(np.trim_zeros(co,'f'))if np.any(co)else []
        ts=np.array([0.,math.pi/2]+[math.atan(float(z.real))%math.pi for z in roots if abs(z.imag)<1e-8])
        x=np.cos(ts);z=np.sin(ts);h=np.array([x*x,x*z,z*z]).T
        forms=np.einsum('ni,ij,nj->n',h,q,h);j=np.argmin(forms);result[i]=forms[j];angles[i]=ts[j]
    answer=(R,nt,N,result,transverse,angles)
    if details:return answer+(dict(transverse_xi=np.stack((vt[:,0,0],np.zeros(len(k)),vt[:,1,0]),axis=1)),)
    return answer


def wall_batch(K,H,h,nu):
    if np.any(~np.isfinite(K))or np.any(K<=0):raise ValueError('nonpositive sampled K')
    a=KAPPA*np.sqrt(K/QE);A=(K/QE)**1.5/KAPPA;x=h*np.sqrt(K)/(70*nu)
    b=-10.7*np.expm1(np.maximum(.5-.4*x**.7,1-.79*x**(-.28))*np.log(x));zv=nu*b/np.sqrt(K)
    start=np.minimum(H,np.maximum(0,zv));d=nu+a*np.maximum(start-zv,0)
    j0=start/nu+np.log1p(a*(H-start)/d)/a;j2=start/nu**2+(H-start)/(d*(d+a*(H-start)))
    lam=(K/QE)**2*(j0-nu*j2)/A;rc=H*np.exp(1-lam)
    if np.any(lam<1)or np.any(rc<=0)or np.any(rc>H):raise ValueError('counterfactual cutoff inadmissible')
    return a,A,zv,lam,rc


class Rows:
    """Sparse physical observation rows, fixed by mesh/map and explicit owners."""
    def __init__(self,r,points,c,*,owner_lookup=None,geometry_cache=None,include_volume=True):
        self.report=r;self.config=c;rows=[];weights=[]
        for p in points:
            for o in (owner_lookup(p)if owner_lookup is not None else owners(r,p)):rows.append((list(p),o));weights.append(0.)
        self.point_count=len(rows);g,w=np.polynomial.legendre.leggauss(8)
        for el in (r['fields'][0]['mesh']if include_volume else ()):
            lo,hi=el['lower'],el['upper'];za=physical(r,lo)[2];zb=physical(r,hi)[2]
            cuts=sorted(set([za,zb]+[c.d0+v for v in (c.w['rc'],c.h,max(0.,c.w['zv']))if za<c.d0+v<zb]))
            for a,b in zip(cuts[:-1],cuts[1:]):
                for i,x in enumerate((lo[0]+hi[0])/2+(hi[0]-lo[0])/2*g):
                    for j,z in enumerate((a+b)/2+(b-a)/2*g):
                        p=[float(x),0,float(z)];o=next(o for o in owners(r,p)if o['geometry']==el['geometry'])
                        rows.append((p,o));weights.append(float(w[i]*w[j]*(hi[0]-lo[0])*(b-a)/4))
        point_basis=basis if geometry_cache is None else geometry_cache.basis
        point_matching=(lambda o:matching(r,o))if geometry_cache is None else geometry_cache.matching
        self.rows=rows;self.weights=np.array(weights);self.outer=np.array([o['region']=='outer'for p,o in rows]);self.rphys=np.array([p[2]-c.d0 for p,o in rows]);self.maps={}
        for block in range(6):
            active=np.flatnonzero(self.outer if block in (3,5)else ~self.outer if block==4 else np.ones(len(rows),bool));compiled=[]
            elements={e['geometry']:e for e in r['fields'][block]['mesh']}
            for i in active:
                p,o=rows[i];el=elements[o['geometry']];lo,hi=el['lower'],el['upper'];x,_,rho=o['chart']
                px,dx=point_basis(el['degree'],[2*(x-lo[0])/(hi[0]-lo[0])-1]);pr,dr=point_basis(el['degree'],[2*(rho-lo[2])/(hi[2]-lo[2])-1]);J=1 if not self.outer[i]else mapping(r)[3]
                compiled.append((np.asarray(el['ids'])+r['offsets'][block],np.outer(pr[0],px[0]).ravel()*el['sign'],
                    np.outer(pr[0],dx[0]).ravel()*el['sign']*2/(hi[0]-lo[0]),np.outer(dr[0],px[0]).ravel()*el['sign']*2/(hi[2]-lo[2])/J))
            nc=(r['fields'][block]['mesh'][0]['degree']+1)**2
            self.maps[block]=(active,tuple(np.array(v)for v in zip(*compiled))if compiled else
                (np.empty((0,nc),dtype=int),np.empty((0,nc)),np.empty((0,nc)),np.empty((0,nc))))
        # Endpoint K uses its own outer interface value, never local k.
        self.kmap=[]
        for p,o in rows:
            m=point_matching(o);el=next(e for e in r['fields'][2]['mesh']if e['geometry']==m['geometry']);lo,hi=el['lower'],el['upper']
            px=point_basis(3,[2*(m['chart'][0]-lo[0])/(hi[0]-lo[0])-1])[0][0];pr=point_basis(3,[-1])[0][0]
            self.kmap.append((np.asarray(el['ids'])+r['offsets'][2],np.outer(pr,px).ravel()*el['sign']))
        self.kids=np.array([v[0]for v in self.kmap]);self.kphi=np.array([v[1]for v in self.kmap]);self.exact=self.reference()

    def field(self,state,b):
        active,(ids,phi,gx,gz)=self.maps[b];cs=np.asarray(state)[ids];v=np.zeros(len(self.rows));G=np.zeros((len(self.rows),3))
        v[active]=np.sum(phi*cs,axis=1);G[active,0]=np.sum(gx*cs,axis=1);G[active,2]=np.sum(gz*cs,axis=1);return v,G

    def evaluate(self,state,principal=True):
        n=len(self.rows);outer=self.outer;inner=~outer;c=self.config;U=np.zeros((n,3));G=np.zeros((n,3,3))
        U[:,0],G[:,0]=self.field(state,0);U[:,2],G[:,2]=self.field(state,1);k,gk=self.field(state,2);eps,ge=self.field(state,3)
        if np.any(k<=0)or np.any(eps[outer]<=0):raise ValueError('expanded nonpositive scalar')
        pi=self.field(state,4)[0]+self.field(state,5)[0];K=np.sum(np.asarray(state)[self.kids]*self.kphi,axis=1)
        a,A,zv,lam,rc=wall_batch(K,c.H,c.h,c.nu);nt=a*np.maximum(self.rphys-zv,0)
        R=2*k[:,None,None]/3*np.eye(3)-nt[:,None,None]*(G+G.transpose(0,2,1));T=(c.nu+nt)[:,None,None]*(G+G.transpose(0,2,1))-pi[:,None,None]*np.eye(3)
        root=np.zeros(n);m=np.full(n,np.nan);trans=m.copy();angles=m.copy();nu_eff=nt.copy();nu_d=c.nu+nt
        if principal:
            ro,no,N,mo,to,ao=principal_batch(k[outer],eps[outer],G[outer]);m[outer]=mo;trans[outer]=to;angles[outer]=ao
        else:ro,no,N=algebra(k[outer],eps[outer],G[outer])
        R[outer]=ro;nu_eff[outer]=no;nu_d[outer]=CD*k[outer]**2/eps[outer];root[outer]=N
        T[outer]=-ro+(2*k[outer]/3-pi[outer])[:,None,None]*np.eye(3)
        return dict(U=U,G=G,k=k,epsilon=eps,pi=pi,mechanical_pressure=pi-2*k/3,R=R,T=T,K=K,rc=rc,Lambda=lam,nu_eff=nu_eff,nu_d=nu_d,
            root=root,margin=m,transverse=trans,angle=angles,covariance=np.linalg.eigvalsh(R),divergence=np.trace(G,axis1=1,axis2=2),
            production=-np.sum(R*G,axis=(1,2)),turbulent_production=nt*np.sum((G+G.transpose(0,2,1))*G,axis=(1,2)),reference_production=(K/QE)**2*nt/(c.nu+nt)**2,destruction=A/np.maximum(self.rphys,rc),
            grad_k=gk,grad_epsilon=ge)

    def reference(self):
        c=self.config;n=len(self.rows);rr=self.rphys;outer=self.outer;U=np.zeros((n,3));G=np.zeros((n,3,3));k=np.full(n,c.K);eps=np.zeros(n);R=np.tile(c.Rc,(n,1,1));pi=np.where(outer,0.,-c.K/NE)
        for z in np.unique(rr):
            # Owner/region selects a one-sided reference at H. Subtracting d0
            # can round an outer H coordinate just below H; that must never
            # substitute the finite-nu inner derivative or covariance.
            ii=(rr==z)&~outer
            if np.any(ii):
                U[ii,0]=c.U_inner(z);G[ii,0,2]=c.u_tau**2/c.D(z);k[ii]=c.k_inner(z);R[ii]=c.covariance_inner(z)
            ii=(rr==z)&outer
            if np.any(ii):U[ii,0]=c.U_outer(z);G[ii,0,2]=c.u_tau/(KAPPA*z)
        eps[outer]=c.A/rr[outer]
        return dict(U=U,G=G,k=k,epsilon=eps,pi=pi,mechanical_pressure=pi-2*k/3,R=R)

    def norms(self,a,b):
        out={};c=self.config
        for scope,mask in (('full',np.ones(len(self.rows),bool)),('inner',~self.outer),('outer',self.outer)):
            out[scope]={}
            for name in FIELDS:
                if name=='epsilon'and scope=='inner':continue
                active=mask & (self.outer if name=='epsilon'else True);delta=np.asarray(a[name])-b[name];d2=np.sum(delta.reshape(len(delta),-1)**2,axis=1)
                out[scope][name]=dict(L2=float(np.sqrt(self.weights[active]@d2[active])),sampled_max=float(np.sqrt(np.max(d2[active]))))
        return out

    def witness(self,data,i):
        p,o=self.rows[i];return dict(physical=p,owner=o,**{k:np.asarray(v[i]).tolist()for k,v in data.items()})

    def interface_history(self,data):
        c=self.config;groups={};wall=[]
        for i,(p,o)in enumerate(self.rows[:self.point_count]):
            if o['chart'][2]==0:wall.append(i)
            if o['interface_role'].startswith('matching_'):groups.setdefault(float(p[0]).hex(),{})[o['region']]=i
        rows=[]
        for group in groups.values():
            if set(group)!=set(('inner','outer')):continue
            i,j=group['inner'],group['outer'];K=data['K'][j];a,A,zv,lam,rc=wall_batch(np.array([K]),c.H,c.h,c.nu)
            rows.append(dict(x=self.rows[j][0][0],K=float(K),rc=float(rc[0]),Lambda=float(lam[0]),
                shared_U=float(np.max(abs(data['U'][i]-data['U'][j]))),shared_k=float(data['k'][i]-data['k'][j]),
                epsilon=float(data['epsilon'][j]),epsilon_target=float(A[0]/c.H),
                traction_jump=(data['T'][i,:,2]-data['T'][j,:,2]).tolist(),
                k_flux_jump=float(data['nu_d'][i]*data['grad_k'][i,2]-data['nu_d'][j]*data['grad_k'][j,2]),
                pressure_jump=float(data['pi'][i]-data['pi'][j])))
        ii=~self.outer;w=self.weights
        return dict(rows=rows,wall_shear_range=[float(np.min(data['T'][wall,0,2])),float(np.max(data['T'][wall,0,2]))]if wall else None,
            wall_shear_scope='sampled_ground_rows'if wall else 'not_evaluated:no_ground_rows_in_requested_inventory',
            inner_integrals=dict(actual_turbulent_production=float(w[ii]@data['turbulent_production'][ii]),raw_stress_work=float(w[ii]@data['production'][ii]),
                counterfactual_production=float(w[ii]@data['reference_production'][ii]),cutoff_destruction=float(w[ii]@data['destruction'][ii])),
            note='actual evolving production separate from current-K counterfactual cutoff budget; no stationary restoration')

    def summary(self,data,initial):
        c=self.config;outer=self.outer;point=np.arange(len(self.rows))<self.point_count;active=np.flatnonzero(outer & point);i=active[np.argmin(data['margin'][active])]
        minm=float(data['margin'][i]);mint=float(np.min(data['transverse'][active]));cov=float(np.min(data['covariance']));minscalar=min(float(np.min(data['k'])),float(np.min(data['epsilon'][outer])))
        status=classify(minm,mint,cov,minscalar)
        errors=dict(E_G=float(np.max(abs(data['G'][active,0,2]/self.exact['G'][active,0,2]-1))),E_epsilon=float(np.max(abs(data['epsilon'][active]/self.exact['epsilon'][active]-1))),
            E_U=float(np.max(abs(data['U'][active]-self.exact['U'][active]))),E_k=float(np.max(abs(data['k'][active]/c.K-1))))
        return dict(status=status,min_in_plane=minm,min_transverse=mint,min_covariance=cov,min_scalar=minscalar,
            witness=self.witness(data,i),errors=errors,interface_budget=self.interface_history(data),reference_error=self.norms(data,self.exact),initial_drift=self.norms(data,initial),
            raw_divergence_max=float(np.max(abs(data['divergence']))),row_count=self.point_count,quadrature_rows=len(self.rows)-self.point_count)


def classify(m,t,c,s):
    if not all(math.isfinite(v)for v in (m,t,c,s)):return 'UNRESOLVED'
    if c<0 or s<=0 or min(m,t)<-2e-11:return 'OBSTRUCTED'
    if min(m,t)<=2e-11:return 'UNRESOLVED'
    return 'QUALITY_LOST'if m<.5*EXACT_MARGIN else 'PASS'
