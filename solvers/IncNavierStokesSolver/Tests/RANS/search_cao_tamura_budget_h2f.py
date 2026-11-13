#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded outer phase I with independently hard affine inner covariance cuts."""
import argparse,os,time,resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from scipy.optimize import minimize
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,require
from search_cao_tamura_budget_h2d import tensors

class WorkLimit(Exception):pass

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('context','preflight','start','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--divergence-limit',type=float);p.add_argument('--stage',choices=['A','B'],required=True);p.add_argument('--prior-full-scans',type=int,default=0);a=p.parse_args();a.output.mkdir(exist_ok=False);setup=time.monotonic();c=NativeContext(a.context);pre=load(a.preflight);K0=pre['heights'][-1]['K0'];target=1e-4;tol=1e-10
    from cao_tamura_h2g_constitutive import tensors_for_spec
    def tensors(g,k,e,inner,nt,physical,us,kap):return tensors_for_spec(g,k,e,inner,nt,physical,us,kap,c.model)
    T,Ts,conditioning=c.transforms();nv=T.shape[1];nk=c.counts['k'];ne=c.counts['epsilon'];ns=Ts.shape[1];dim=nv+(ns if a.stage=='B' else 0);G=T.T@c.array('mass_gram')@T
    transform=np.zeros((len(T)+1,nv+1));transform[0,0]=1;transform[1:,1:]=T;D=transform.T@c.array('divergence')@transform;strain=transform.T@c.array('strain')@transform;roof=transform.T@c.array('roof')@transform;Spen=Ts.T@(c.array('scalar_gram')+c.array('scalar_smoothness'))@Ts
    basis=np.array(c.array('basis'));Z=np.vstack([T.T@basis[:len(T)],Ts.T@basis[len(T):]]);base=np.array(c.array('base')[0]);raw=c.block('velocity_basis');B=np.einsum('nbc,bj->ncj',raw[:,:,:4],T);V=np.einsum('nbc,bj->ncj',raw[:,:,4:],T);g0=c.gradient();k0=c.block('k');e0=c.block('epsilon');nt=c.block('inner_nu_t');physical=c.block('physical');inner=c.block('inner').astype(bool);cov=c.block('covariance').astype(bool);dep=c.block('dependent_epsilon').astype(bool);innerids=np.flatnonzero(inner&cov);outerids=np.flatnonzero(~inner&cov)
    rawS=np.column_stack([c.block('k_basis'),c.block('epsilon_basis')]);SB=rawS@Ts;BK=np.zeros_like(SB);BE=np.zeros_like(SB);BK[:,:nk]=SB[:,:nk];BE[:,nk:]=SB[:,nk:];require(np.max(abs(BK[inner]))==0 and np.max(abs(BE[inner|dep]))==0,'frozen scalar functional changed')
    start=load(a.start);x=np.zeros(dim);oldx=start['conditioned_velocity_coefficients'];x[:len(oldx)]=oldx
    if a.stage=='B' and start.get('conditioned_scalar_coefficients'):x[nv:]=start['conditioned_scalar_coefficients']
    free=np.array(c.meta['rows']['velocity']+(c.meta['rows']['k']+c.meta['rows']['epsilon'] if a.stage=='B' else []),int);Dmax=np.sqrt(D[0,0])*(1+1e-6)+1e-9
    if a.divergence_limit is not None:
        require(abs(a.divergence_limit-Dmax)<1e-10,'shared divergence limit disagrees with native base');Dmax=a.divergence_limit
    full=a.prior_full_scans;inner_scans=active=0;history=[];cuts=[];cutkeys=set();responsecaps=set();scalarcaps=set();maximum_cuts=2000;maxiter=20 if a.stage=='A' else 40;seconds=300 if a.stage=='A' else 900;fullcap=180 if a.stage=='A' else 320;activecap=3000 if a.stage=='A' else 6000;began=None;deadline=None
    def check_work():
        if began is not None and (time.monotonic()>deadline or full+inner_scans>=fullcap-2 or active>=activecap):raise WorkLimit()
    def values(v,ids):
        dg=np.einsum('ncj,j->nc',B[ids],v[:nv]);g=g0[ids].copy()
        for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):g[:,i,j]+=dg[:,col]
        k=np.array(k0[ids]);e=np.array(e0[ids])
        if a.stage=='B':k+=BK[ids]@v[nv:];e+=BE[ids]@v[nv:]
        return g,k,e
    def stress(v,ids):g,k,e=values(v,ids);return tensors(g,k,e,inner[ids],nt[ids],physical[ids],pre['u_tau'],pre['kappa_CLS'])
    def scan(v):
        nonlocal full
        full+=1;Ri=stress(v,innerids);Ro=stress(v,outerids);vi,ei=np.linalg.eigh(Ri);vo,eo=np.linalg.eigh(Ro);return dict(inner_values=vi/K0,inner_vectors=ei,outer_values=vo/K0,outer_vectors=eo)
    def inner_scan(v):
        nonlocal inner_scans
        inner_scans+=1;ev,vec=np.linalg.eigh(stress(v,innerids));return ev/K0,vec
    def quad(A,v):z=np.r_[1.,v[:nv]];return float(z@A@z)
    def distance(v):return float(v[:nv]@G@v[:nv])
    def penalty(v):return .001*distance(v)+.001*quad(strain,v)+.01*quad(roof,v)+(.01*v[nv:]@Spen@v[nv:] if a.stage=='B' else 0)
    def pgrad(v):return np.r_[.002*G@v[:nv]+.002*(strain@np.r_[1,v[:nv]])[1:]+.02*(roof@np.r_[1,v[:nv]])[1:],.02*Spen@v[nv:] if a.stage=='B' else []]
    def bounds(v):
        response=np.einsum('ncj,j->nc',V,v[:nv]);kr=np.ones(c.N);er=np.ones(c.N);outside=~inner
        if a.stage=='B':kr[outside]+=BK[outside]@v[nv:]/k0[outside];er[outside]+=BE[outside]@v[nv:]/e0[outside]
        return response,kr,er
    def within_bounds(v,metrics=None):
        response,kr,er=bounds(v) if metrics is None else metrics
        return distance(v)<=.12**2*(1+1e-8) and quad(D,v)<=Dmax**2+1e-12 and np.max(abs(response))<=.5+1e-10 and min(kr.min(),er.min())>=.5-1e-10 and max(kr.max(),er.max())<=2+1e-10
    def add_inner(ev,evec,threshold=None):
        selected=np.argsort(ev[:,0])[:24] if threshold is None else np.flatnonzero(ev[:,0]<threshold)
        if len(selected)>24:selected=selected[np.argsort(ev[selected,0])[:24]]
        added=0
        for local in selected:
            q=int(innerids[local]);R0=stress(np.zeros(dim),[q])[0]
            for eigen in np.flatnonzero(ev[local]-ev[local,0]<=1e-7):
                v=evec[local,:,eigen];vv=np.outer(v,v);key=(q,vv.tobytes())
                if key in cutkeys:continue
                cutkeys.add(key);row=-2*nt[q]/K0*(vv[0,0]*B[q,0]+vv[0,2]*(B[q,1]+B[q,2])+vv[2,2]*B[q,3]);constant=float(np.sum(vv*R0))/K0
                cuts.append((q,v.tolist(),row,constant));added+=1
        return added
    def minimum(s):return float(s['outer_values'][:,0].min())
    def merit(v,s):return max(0,target-minimum(s))+penalty(v)
    def save(v,s,label):
        require(float(s['inner_values'][:,0].min())>=target-tol,'attempt to save inner-infeasible state');require(within_bounds(v),'attempt to save physical bound violation')
        state=base.copy();state[free]+=(v@Z[:dim])[free];fixed=np.ones(len(state),bool);fixed[free]=False;require(state[fixed].tobytes()==base[fixed].tobytes(),'fixed coefficient bits changed')
        response,kr,er=bounds(v);j=int(np.argmin(s['outer_values'][:,0]));q=int(outerids[j]);g,k,e=values(v,[q]);out=dict(model_specification=c.model,version='H2fHardInnerCandidateV1',stage=a.stage,state=state.tolist(),state_sha256=state_hash(state),base_sha256=state_hash(base),conditioned_velocity_coefficients=v[:nv].tolist(),conditioned_scalar_coefficients=v[nv:].tolist(),min_lambda_over_K0=minimum(s),inner_min_lambda_over_K0=float(s['inner_values'][:,0].min()),divergence_limit=Dmax,strong_divergence_L2=float(np.sqrt(max(0,quad(D,v)))),velocity_mass_distance=float(np.sqrt(max(0,distance(v)))),maximum_sampled_component_increment=float(np.max(abs(response))),K_envelope=[float(kr.min()),float(kr.max())],epsilon_envelope=[float(er.min()),float(er.max())],context_sha256=sha(a.context),selection=label,witness=dict(index=q,chart=c.block('chart')[q].tolist(),physical=physical[q].tolist(),gradient=g[0].tolist(),k=float(k[0]),epsilon=float(e[0]),inner=False))
        atomic(a.output/'best-candidate.json',out)
    # Setup and verification are separate from the bounded optimizer clock.
    s=scan(x);require(within_bounds(x) and s['inner_values'][:,0].min()>=target-tol,'starting state not hard-inner-feasible')
    reconstructed=base.copy();reconstructed[free]+=(x@Z[:dim])[free];initial_reconstruction=float(np.max(abs(reconstructed-np.array(start['state']))));require(initial_reconstruction<1e-10,'starting physical field changed beyond coordinate roundoff')
    t=time.monotonic();scan(x);scan_cost=time.monotonic()-t;add_inner(s['inner_values'],s['inner_vectors']);best=x.copy();bs=s;save(best,bs,'initial');began=time.monotonic();deadline=began+seconds;radius=.01;stagnant=0;reason='iteration_budget';accepted_margins=[float(s['inner_values'][:,0].min())];subproblems=0;retained=set()
    report=dict(model_specification=c.model,version='H2fHardInnerPhaseIV1',stage=a.stage,conditioning=conditioning,initial_coefficient_reconstruction_max=initial_reconstruction,context_setup_seconds=began-setup,seconds_per_full_scan=scan_cost,inner_target=target,hard_tolerance_unchanged=True,prior_full_scan_attempts=a.prior_full_scans,history=history,status='running',K0=K0)
    try:
        for iteration in range(maxiter):
            check_work();predictor=iteration==0;s=scan(x);before=merit(x,s);before_outer=minimum(s);before_inner=float(s['inner_values'][:,0].min());candidates=np.argsort(s['outer_values'][:,0])[:64];retained.update(map(int,candidates));local=np.array(sorted(retained),int)
            if len(local)>192:local=local[np.argsort(s['outer_values'][local,0])[:192]]
            ids=outerids[local];vcrit=s['outer_vectors'][local,:,0];av=s['outer_values'][local,0];g,k,e=values(x,ids);J=np.zeros((len(ids),dim))
            # Fixed-eigenvector Rayleigh derivatives respect the actual max
            # invariant branch at each perturbed constitutive evaluation.
            for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):
                if c.model['name']=='StandardKEpsilon':
                    deriv=-2*.09*k*k/e*vcrit[:,i]*vcrit[:,j]/K0;active+=1
                else:
                    h=2e-6*np.maximum(1,abs(g[:,i,j]));gp=g.copy();gm=g.copy();gp[:,i,j]+=h;gm[:,i,j]-=h
                    Rp=tensors(gp,k,e,inner[ids],nt[ids],physical[ids],pre['u_tau'],pre['kappa_CLS']);Rm=tensors(gm,k,e,inner[ids],nt[ids],physical[ids],pre['u_tau'],pre['kappa_CLS']);deriv=np.einsum('ni,nij,nj->n',vcrit,Rp-Rm,vcrit)/(2*h*K0);active+=2
                J[:,:nv]+=deriv[:,None]*B[ids,col]
            if a.stage=='B':
                for scalar,fieldsB in [('k',BK),('epsilon',BE)]:
                    h=2e-6*(k if scalar=='k' else e);kp=k.copy();km=k.copy();ep=e.copy();em=e.copy()
                    if scalar=='k':kp+=h;km-=h
                    else:ep+=h;em-=h
                    Rp=tensors(g,kp,ep,inner[ids],nt[ids],physical[ids],pre['u_tau'],pre['kappa_CLS']);Rm=tensors(g,km,em,inner[ids],nt[ids],physical[ids],pre['u_tau'],pre['kappa_CLS']);J[:,nv:]+=np.einsum('ni,nij,nj->n',vcrit,Rp-Rm,vcrit)[:,None]/(2*h[:,None]*K0)*fieldsB[ids];active+=2
            response,kr,er=bounds(x);flat=response.reshape(-1);near=np.flatnonzero(abs(flat)>.4);responsecaps.update(map(int,near[np.argsort(-abs(flat[near]))[:128]]))
            if a.stage=='B':
                for block,r in enumerate((kr,er)):
                    near=np.flatnonzero((r<.65)|(r>1.8));scalarcaps.update((block,int(i)) for i in near[:128])
            newcuts=0;capadds=0;solution=None
            for separation in range(4):
                check_work();Ain=np.vstack([v[2] for v in cuts]);cin=np.array([v[3] for v in cuts]);ci=np.array(sorted(responsecaps),int);C=np.zeros((len(ci),dim));C[:,:nv]=V.reshape(-1,nv)[ci];cb=flat[ci];si=sorted(scalarcaps);As=np.zeros((len(si),dim));sb=np.zeros(len(si))
                for j,(block,i) in enumerate(si):As[j,nv:]=BK[i]/k0[i] if block==0 else BE[i]/e0[i];sb[j]=(kr if block==0 else er)[i]
                def state(y):return x+y[:-1]
                def padded(A):return np.column_stack([A,np.zeros(len(A))])
                def objective(y):return (0 if predictor else y[-1])+penalty(state(y))
                def jac(y):return np.r_[pgrad(state(y)),0 if predictor else 1.]
                cons=[dict(type='ineq',fun=lambda y:cin+Ain@state(y)[:nv]-target,jac=lambda y:np.column_stack([Ain,np.zeros((len(Ain),dim-nv+1))])),dict(type='ineq',fun=lambda y:1-distance(state(y))/.12**2,jac=lambda y:np.r_[-2*G@state(y)[:nv],np.zeros(dim-nv+1)]/.12**2),dict(type='ineq',fun=lambda y:(Dmax**2-quad(D,state(y)))/Dmax**2,jac=lambda y:np.r_[-2*(D@np.r_[1,state(y)[:nv]])[1:],np.zeros(dim-nv+1)]/Dmax**2),dict(type='ineq',fun=lambda y:1-(y[:nv]@G@y[:nv]+y[nv:-1]@y[nv:-1])/radius**2,jac=lambda y:np.r_[-2*G@y[:nv],-2*y[nv:-1],0.]/radius**2)]
                if not predictor:cons.append(dict(type='ineq',fun=lambda y:av+J@y[:-1]+y[-1]-target,jac=lambda y:np.column_stack([J,np.ones(len(J))])))
                if len(ci):cons += [dict(type='ineq',fun=lambda y:.5-cb-C@y[:-1],jac=lambda y:-padded(C)),dict(type='ineq',fun=lambda y:.5+cb+C@y[:-1],jac=lambda y:padded(C))]
                if len(si):cons +=[dict(type='ineq',fun=lambda y:2-sb-As@y[:-1],jac=lambda y:-padded(As)),dict(type='ineq',fun=lambda y:sb+As@y[:-1]-.5,jac=lambda y:padded(As))]
                scale=np.r_[np.full(dim,radius),1.];scaled=[dict(type='ineq',fun=lambda w,c=c:c['fun'](w*scale),jac=lambda w,c=c:c['jac'](w*scale)*scale) for c in cons];subproblems+=1
                solution=minimize(lambda w:objective(w*scale),np.r_[np.zeros(dim),max(0,target-float(av.min()))],jac=lambda w:jac(w*scale)*scale,method='SLSQP',bounds=[(None,None)]*dim+[(0,None)],constraints=scaled,callback=lambda w:check_work(),options=dict(maxiter=70,ftol=1e-10));solution.x*=scale;proposed=x+solution.x[:-1];metrics=bounds(proposed);added=False;rm,km,em=metrics
                bad=np.flatnonzero(abs(rm.reshape(-1))>.5+1e-10)
                for j in bad[np.argsort(-abs(rm.reshape(-1)[bad]))[:64]]:
                    if int(j) not in responsecaps:responsecaps.add(int(j));added=True;capadds+=1
                if a.stage=='B':
                    for block,r in enumerate((km,em)):
                        for j in np.flatnonzero((r<.5-1e-10)|(r>2+1e-10))[:64]:
                            if (block,int(j)) not in scalarcaps:scalarcaps.add((block,int(j)));added=True;capadds+=1
                iv,ie=inner_scan(proposed)
                if iv[:,0].min()<target-tol:newcuts+=add_inner(iv,ie,target);added=True
                if len(cuts)>maximum_cuts:reason='inner_cut_budget';raise WorkLimit()
                if not added:break
            accepted=False;trial_s=s;fraction_used=0
            for fraction in (1.,.5,.25,.125,.0625):
                check_work();trial=x+fraction*solution.x[:-1]
                if not within_bounds(trial):continue
                iv,ie=inner_scan(trial)
                if iv[:,0].min()<target-tol:newcuts+=add_inner(iv,ie,target);continue
                ts=scan(trial)
                if merit(trial,ts)<before-1e-8 and minimum(ts)>=before_outer-1e-5:x=trial;trial_s=ts;accepted=True;fraction_used=fraction;accepted_margins.append(float(iv[:,0].min()));break
            after=merit(x,trial_s);oldradius=radius;improvement=minimum(trial_s)-before_outer
            if accepted:radius=min(.04,radius*1.25);stagnant=0 if improvement>=1e-5 else stagnant+1
            else:radius*=.5;stagnant+=1
            if minimum(trial_s)>minimum(bs):best=x.copy();bs=trial_s;save(best,bs,'improvement')
            history.append(dict(iteration=iteration,predictor=predictor,accepted=accepted,fraction=fraction_used,true_objective_before=before,true_objective_after=after,outer_before=before_outer,outer_after=minimum(trial_s),inner_before=before_inner,inner_after=float(trial_s['inner_values'][:,0].min()),inner_cuts=len(cuts),new_inner_cuts=newcuts,active_response_caps=len(responsecaps),new_cap_rows=capadds,active_scalar_caps=len(scalarcaps),radius_before=oldradius,radius_after=radius,local_solver_success=bool(solution.success),local_solver_message=str(solution.message),mass_distance=float(np.sqrt(max(0,distance(x)))),strong_divergence_L2=float(np.sqrt(max(0,quad(D,x)))),full_tensor_scans=full,complete_inner_scans=inner_scans,active_derivative_batches=active,seconds=time.monotonic()-began))
            report.update(history=history,best_outer_minimum=minimum(bs),minimum_accepted_inner_margin=min(accepted_margins));atomic(a.output/'progress.json',report)
            if minimum(bs)>=target:reason='complete_interior_target';break
            if iteration+1>=maxiter:reason='iteration_budget';break
            if not predictor and (stagnant>=6 or radius<1e-6):reason='relative_progress_threshold' if accepted else 'stagnation';break
    except WorkLimit:
        if reason!='inner_cut_budget':reason='work_budget'
    # Reserved full scans retain the last actual accepted state before exit.
    terminal=scan(x);require(terminal['inner_values'][:,0].min()>=target-tol,'terminal inner protection failed')
    if minimum(terminal)>minimum(bs):best=x.copy();bs=terminal
    bs=scan(best);save(best,bs,'terminal_best');report.update(status='sampled_candidate' if minimum(bs)>=target else 'bounded_search_rejected',termination=reason,best_outer_minimum=minimum(bs),best_inner_minimum=float(bs['inner_values'][:,0].min()),minimum_accepted_inner_margin=min(accepted_margins),accepted_states=len(accepted_margins),full_tensor_scans=full,complete_inner_scans=inner_scans,conservative_full_scan_work=full+inner_scans,active_derivative_batches=active,subproblems=subproblems,search_seconds=time.monotonic()-began,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,classification=('standard fixed-scalar affine-tensor convex finite-inventory constraints; bounded solver stop is not an infeasibility certificate' if c.model['name']=='StandardKEpsilon' and a.stage=='A' else 'bounded local nonconvex search; finite inventory; no global infeasibility certificate'));atomic(a.output/'progress.json',report)
if __name__=='__main__':main()
