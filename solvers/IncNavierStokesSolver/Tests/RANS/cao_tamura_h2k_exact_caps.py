# SPDX-License-Identifier: MIT
"""Pure-slack native candidate construction with exact physical quadratics.

Rayleigh/cap separation remains complete at proposed accepted states. Scalar
min/max constraints are equivalent to every declared scalar row, grouped only
to avoid a large dense SLSQP Jacobian. Numerical status cannot certify a seed.
"""
import time
import numpy as np
from scipy import linalg,optimize
from cao_tamura_h2h_scalar import BudgetExhausted
from cao_tamura_h2i_convex import summary,phase as linear_phase


def exact_caps_phase(P,delta,band,budget,save,initial=None,lp_solver=None):
    n=P.dim;nv=P.nv;ne=P.ne;started=time.monotonic();x=np.zeros(n) if initial is None else np.array(initial);cuts=[];rhs=[];keys=set();active=set();best=None;best_summary=None;best_slack=float('inf');rounds=[];termination='inconclusive';terminal=x.copy();last_full_slack=None
    G=.5*(P.mass[1:,1:]+P.mass[1:,1:].T);h=.5*(P.mass[0,1:]+P.mass[1:,0]);centre=-linalg.solve(G,h,assume_a='pos');radius=np.sqrt(.12**2-P.mass[0,0]-h@centre)
    T=np.zeros((n+1,n+1));T[:nv,:nv]=radius*linalg.solve_triangular(np.linalg.cholesky(G).T,np.eye(nv),lower=False);scale_e=np.array([max(abs(a),abs(b)) for a,b in P.bounds[nv:]]);T[nv:n,nv:n]=np.diag(scale_e);T[n,n]=1.;shift=np.r_[centre,np.zeros(ne+1)]
    scalar=P.E/P.ec[:,None];ratio=P.e/P.ec;groups=[np.flatnonzero((~P.inner)&(P.geometry==g)) for g in np.unique(P.geometry[~P.inner])];cache={}
    def scalars(y):
        e=y[nv:n]*scale_e;key=e.tobytes()
        if cache.get('key')!=key:
            q=ratio+scalar@e;idslo=np.array([ix[np.argmin(q[ix])] for ix in groups]);idshi=np.array([ix[np.argmax(q[ix])] for ix in groups]);values=np.r_[q[idslo]-band[0],band[1]-q[idshi]];jac=np.zeros((len(values),n+1));jac[:len(groups),nv:n]=scalar[idslo]*scale_e;jac[len(groups):,nv:n]=-scalar[idshi]*scale_e;cache.update(key=key,values=values,jac=jac)
        return cache['values']
    quads=[(P.mass,.12**2)]+list(zip(P.div,P.div_limits));symmetric=[(.5*(Q+Q.T),limit) for Q,limit in quads]
    def physical(y):
        a=shift[:nv]+T[:nv,:nv]@y[:nv];u=np.r_[1.,a];return np.array([1-u@Q@u/limit for Q,limit in symmetric])
    def physical_jac(y):
        a=shift[:nv]+T[:nv,:nv]@y[:nv];u=np.r_[1.,a];jac=np.zeros((4,n+1))
        for j,(Q,limit) in enumerate(symmetric):jac[j,:nv]=-2*(Q@u)[1:]@T[:nv,:nv]/limit
        return jac
    def add(row,b,key):
        if key in keys:return
        keys.add(key);s=max(1,float(np.max(abs(row))),abs(b));cuts.append(row/s);rhs.append(b/s)
    def separate(s):
        eligible=np.flatnonzero(P.cov);selected=[]
        for inner in (True,False):
            ids=eligible[P.inner[eligible]==inner];selected.extend(ids[np.argsort(s['lambda_transformed'][ids])[:48]])
        for g in np.unique(P.geometry):
            ids=np.flatnonzero(P.cov&(P.geometry==g))
            if len(ids):selected.append(ids[np.argmin(s['lambda_transformed'][ids])])
        active.update(int(i) for i in selected)
        dv=s['velocity_increment'];ids=np.argsort(abs(dv).ravel())[-96:]
        for ij in ids:
            i,j=divmod(int(ij),2);sign=1 if dv[i,j]>=0 else -1;add(*P.cap_cut(i,j,sign),('cap',i,j,sign))
    def keep(s,v):
        nonlocal best,best_summary,best_slack
        valid=s['hard_violation']<=1e-9 and min(s['epsilon_ratio'])>=band[0]-1e-10 and max(s['epsilon_ratio'])<=band[1]+1e-10
        if valid and s['slack']<best_slack:
            best=v.copy();best_summary=summary(s);best_slack=s['slack'];save(dict(event='best_hard_diagnostic',coefficients=v.tolist(),summary=best_summary))
        return valid
    try:
        s=P.scan(x,delta,budget);keep(s,x);tmax=s['slack']+1.;z=np.r_[x,s['slack']+1e-6];separate(s)
        bounds=[(-1-1e-8,1+1e-8)]*nv+[(-1.,1.)]*ne+[(0.,tmax)]
        for iteration in range(48):
            budget.check();A=np.array(cuts);b=np.array(rhs);Ay=A@T;by=b-A@shift;y=linalg.solve(T,z-shift)
            ids=np.array(sorted(active));k=P.k[ids];e=P.e[ids];inner=P.inner[ids];c=2*k/3-delta;F=np.where(inner,P.nt[ids]/P.K0,.09*k*k/(P.K0*P.ec[ids]));iso=np.where(inner,c/P.K0,c*e/(P.K0*P.ec[ids]));g=P.g[ids];VB=P.VB[ids];dxx=np.zeros((len(ids),n+1));dzz=dxx.copy();dxz=dxx.copy();di=dxx.copy();dxx[:,:nv]=-2*F[:,None]*VB[:,:,0];dzz[:,:nv]=-2*F[:,None]*VB[:,:,3];dxz[:,:nv]=-F[:,None]*(VB[:,:,1]+VB[:,:,2]);di[:,nv:n]=np.where(inner,0,c/(P.K0*P.ec[ids]))[:,None]*P.E[ids];dxx+=di;dzz+=di;bx=iso-2*F*g[:,0];bz=iso-2*F*g[:,3];bxz=-F*(g[:,1]+g[:,2]);xx0=bx+dxx@shift;zz0=bz+dzz@shift;xz0=bxz+dxz@shift;ii0=iso+di@shift;Jx=dxx@T;Jz=dzz@T;Jxz=dxz@T;Ji=di@T;soc_cache={}
            def soc(y):
                key=y.tobytes()
                if soc_cache.get('key')!=key:
                    xx=xx0+Jx@y;zz=zz0+Jz@y;xz=xz0+Jxz@y;ii=ii0+Ji@y;norm=np.hypot(xx-zz,2*xz);den=np.maximum(norm,1e-30);grad=.5*(Jx+Jz)-.5*((xx-zz)[:,None]*(Jx-Jz)+4*xz[:,None]*Jxz)/den[:,None];jac=np.vstack((grad,Ji));jac[:,-1]+=1;values=np.r_[.5*(xx+zz-norm)+y[-1],ii+y[-1]];soc_cache.update(key=key,values=values,jac=jac)
                return soc_cache['values']
            def fun(y):budget.check();return y[-1]
            def jac(y):v=np.zeros(n+1);v[-1]=1.;return v
            constraints=[optimize.LinearConstraint(Ay,-np.inf,by),dict(type='ineq',fun=soc,jac=lambda y:(soc(y),soc_cache['jac'])[1]),dict(type='ineq',fun=physical,jac=physical_jac),dict(type='ineq',fun=scalars,jac=lambda y:(scalars(y),cache['jac'])[1])]
            budget.round();tic=time.monotonic();r=optimize.minimize(fun,y,jac=jac,bounds=bounds,constraints=constraints,method='SLSQP',options=dict(ftol=1e-12,maxiter=120,disp=False));z=shift+T@r.x;x=z[:n];terminal=x.copy()
            rec=dict(method='SLSQP_exact_PSD_and_physical_quadratics',status=int(r.status),message=r.message,seconds=time.monotonic()-tic,iterations=int(r.nit),function_evaluations=int(r.nfev),slack_objective=float(z[-1]),active_PSD_points=len(active),linear_cap_rows=len(cuts));rounds.append(rec);save(dict(event='exact_cap_relaxation',round=budget.rounds,receipt=rec,coefficients=x.tolist()))
            if np.any(P.e+P.E@x[nv:]<=0):termination='numerical_scalar_domain_failure';break
            s=P.scan(x,delta,budget);valid=keep(s,x);old_active=len(active);old_cuts=len(cuts);separate(s);save(dict(event='exact_cap_complete_scan',summary=summary(s),new_active_points=len(active)-old_active,new_cap_rows=len(cuts)-old_cuts))
            if valid and min(s['inner_min'],s['outer_min'])>=delta-1e-12:termination='seen_feasible';break
            z[-1]=max(0,s['slack'])+1e-7
            stagnant=last_full_slack is not None and abs(s['slack']-last_full_slack)<1e-9 and old_active==len(active) and old_cuts==len(cuts)
            last_full_slack=s['slack']
            if iteration>=2 and stagnant and s['hard_violation']<=1e-8:break
        if termination!='seen_feasible':
            # A local numerical status is never an obstruction. At an apparent
            # positive-slack optimum, ask the independent necessary-cut LP for
            # a checked lower bound using actual original coefficients.
            oldmax=budget.maxrounds;budget.maxrounds=min(oldmax,budget.rounds+4)
            try:r=linear_phase(P,delta,band,budget,save,initial=best if best is not None else x,lp_solver=lp_solver)
            finally:budget.maxrounds=oldmax
            r['exact_cap_relaxations']=rounds;r['coordinate_transform']=dict(mass_unit_ball=True,velocity_rank=nv,scalar_rank=ne,physical_projection=False);return r
    except BudgetExhausted as error:termination='budget_exhausted:'+str(error)
    result=dict(termination=termination,delta=delta,epsilon_band=band,rounds=rounds,best_coefficients=best.tolist() if best is not None else None,best_summary=best_summary,terminal_coefficients=terminal.tolist(),budget=budget.report(),cut_count=len(cuts),proof=None,coordinate_transform=dict(mass_unit_ball=True,velocity_rank=nv,scalar_rank=ne,physical_projection=False),seconds=time.monotonic()-started)
    save(dict(event='phase_terminal',result=result));return result
