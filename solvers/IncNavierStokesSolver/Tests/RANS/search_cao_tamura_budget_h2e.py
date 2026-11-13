#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Mass-conditioned, divergence-constrained native H2e phase-I search.

No governing solve or state acceptance occurs here. Every saved candidate needs
an independent native production audit. Local SQP models are not certificates.
"""
import argparse,hashlib,json,math,os,time,resource
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from scipy.optimize import minimize
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,require
from search_cao_tamura_budget_h2d import tensors
from validate_cls_epawf_m1 import indexed

def condition(G):
    e,V=np.linalg.eigh((G+G.T)/2);keep=e>e[-1]*1e-10
    require(e[0]>-1e-10*e[-1],'indefinite physical Gram')
    T=V[:,keep]/np.sqrt(e[keep]);return T,dict(eigenvalues=e.tolist(),rank=int(keep.sum()),condition_before=float(e[-1]/e[0]),orthogonality_error=float(np.max(abs(T.T@G@T-np.eye(keep.sum())))),relative_rank_threshold=1e-10)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--context',type=Path,required=True);p.add_argument('--preflight',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--stage',choices=['A','B'],required=True);p.add_argument('--start',type=Path);p.add_argument('--max-major',type=int)
    a=p.parse_args();a.output.mkdir(exist_ok=False);started=time.monotonic();doc=load(a.context);f=doc['context'];pre=load(a.preflight);K0=pre['heights'][-1]['K0'];root=Path(doc['array_directory'])
    def array(name):
        d=f['files'][name];require(sha(root/d['file'])==d['sha256'],'context hash mismatch '+name)
        return np.memmap(root/d['file'],dtype=d['dtype'],mode='r',shape=tuple(d['shape']))
    base=np.array(array('base')[0]);Z0=np.array(array('basis'));G=np.array(array('mass_gram'));T,vc=condition(G);nv=T.shape[1];sg=np.array(array('scalar_gram'));Tk,kc=condition(sg[:8,:8]);Te,ec=condition(sg[8:,8:]);ns=Tk.shape[1]+Te.shape[1]
    Tv=np.zeros((33,nv+1));Tv[0,0]=1;Tv[1:,1:]=T
    Q=Tv.T@array('divergence')@Tv;Qs=Tv.T@array('strain')@Tv;Qr=Tv.T@array('roof')@Tv
    Ts=np.zeros((16,ns));Ts[:8,:Tk.shape[1]]=Tk;Ts[8:,Tk.shape[1]:]=Te;Ps=Ts.T@(sg+array('scalar_smoothness'))@Ts
    Z=np.vstack([T.T@Z0[:32],Ts.T@Z0[32:]])
    obs=array('observations');N=len(obs);physical=obs[:,3:6];g0=np.zeros((N,3,3));g0[:,0,0]=obs[:,6];g0[:,0,2]=obs[:,7];g0[:,2,0]=obs[:,8];g0[:,2,2]=obs[:,9]
    raw=obs[:,18:210].reshape(N,32,6)
    B=np.einsum('nbc,bj->ncj',raw[:,:,:4],T);V=np.einsum('nbc,bj->ncj',raw[:,:,4:],T)
    SB=obs[:,210:]@Ts;k0=obs[:,12];e0=obs[:,13];inner=obs[:,15].astype(bool);cov=obs[:,16].astype(bool);dependent=obs[:,17].astype(bool);nt=obs[:,14]
    BK=np.zeros((N,ns));BE=np.zeros((N,ns));BK[:,:Tk.shape[1]]=SB[:,:Tk.shape[1]];BE[:,Tk.shape[1]:]=SB[:,Tk.shape[1]:]
    # Scalar K rows and epsilon rows share this storage only after native exact
    # row evaluation. All inner and dependent-trace scalar directions are zero.
    require(np.max(abs(BK[inner]))==0 and np.max(abs(BE[inner|dependent]))==0,'scalar support/trace changed')
    Dmax=math.sqrt(Q[0,0])*(1+1e-6)+1e-9;size=nv+(ns if a.stage=='B' else 0);x=np.zeros(size);target=1e-4
    if a.start:
        previous=load(a.start);x[:nv]=previous['conditioned_velocity_coefficients']
    free=np.array(f['rows']['velocity']+(f['rows']['k']+f['rows']['epsilon'] if a.stage=='B' else []),int)
    # Conditioning is a change of coordinates, not a change of physical span.
    old=np.sin(np.arange(32)+.3)*1e-4;y=np.linalg.lstsq(T,old,rcond=1e-12)[0];reconstruction=np.max(abs(old@Z0[:32]-y@Z[:nv]))
    require(reconstruction<1e-11,'old physical increment not recovered')
    checks=dict(velocity=vc,k=kc,epsilon=ec,old_increment_coefficient_error=float(reconstruction),homogeneous_continuity_max=f['homogeneous_continuity_max'],frozen_direction_max=f['frozen_direction_max'],base_state_sha256=state_hash(base),context_sha256=sha(a.context))
    full=active_batches=0
    def fields(v,idx):
        dg=np.einsum('ncj,j->nc',B[idx],v[:nv]);g=g0[idx].copy();g[:,0,0]+=dg[:,0];g[:,0,2]+=dg[:,1];g[:,2,0]+=dg[:,2];g[:,2,2]+=dg[:,3]
        k=np.array(k0[idx]);e=np.array(e0[idx]);
        if a.stage=='B':k+=BK[idx]@v[nv:];e+=BE[idx]@v[nv:]
        return g,k,e
    def evaluate(v,idx=None):
        nonlocal full,active_batches
        if idx is None:idx=np.flatnonzero(cov);full+=1
        else:active_batches+=1
        g,k,e=fields(v,idx);R=tensors(g,k,e,inner[idx],nt[idx],physical[idx],pre['u_tau'],pre['kappa_CLS']);ev=np.linalg.eigvalsh(R)[:,0]/K0
        return ev
    # Fresh same-state batch/indexed comparison on both invariant branches.
    batch=tensors(g0,k0,e0,inner,nt,physical,pre['u_tau'],pre['kappa_CLS']);err=0
    for j in np.flatnonzero(~inner)[::max(1,int((~inner).sum())//30)]:
        t=indexed(g0[j],k0[j],e0[j],'physical');S=g0[j]+g0[j].T;sr=np.zeros((3,3));sr[0,2]=sr[2,0]=pre['u_tau']/(pre['kappa_CLS']*physical[j,2]);r=2*k0[j]*np.eye(3)/3-float(t['nut'].real)*S+t['quadratic'].real+.4*t['cubic'].real-pre['kappa_CLS']*pre['u_tau']*physical[j,2]*(S-sr);err=max(err,float(np.max(abs(r-batch[j]))))
    require(err<1e-11,'batch/indexed stress mismatch');del batch
    def quad(Q,v):z=np.r_[1.,v[:nv]];return float(z@Q@z)
    def penalty(v):return .001*(v[:nv]@v[:nv])+.001*quad(Qs,v)+.01*quad(Qr,v)+(.01*(v[nv:]@Ps@v[nv:]) if a.stage=='B' else 0)
    def pgrad(v):r=np.r_[2*.001*v[:nv]+2*.001*(Qs@np.r_[1.,v[:nv]])[1:]+2*.01*(Qr@np.r_[1.,v[:nv]])[1:],np.zeros(size-nv)];r[nv:]=.02*Ps@v[nv:] if a.stage=='B' else r[nv:];return r
    def bounds_metrics(v):
        response=np.einsum('ncj,j->nc',V,v[:nv]);kr=np.ones(N);er=np.ones(N);valid=~inner
        if a.stage=='B':kr[valid]+=BK[valid]@v[nv:]/k0[valid];er[valid]+=BE[valid]@v[nv:]/e0[valid]
        return response,kr,er
    def valid(v,metrics=None):
        response,kr,er=bounds_metrics(v) if metrics is None else metrics
        return np.linalg.norm(v[:nv])<=.12*(1+1e-8) and quad(Q,v)<=Dmax**2+1e-12 and np.max(abs(response))<=.5+1e-10 and min(kr.min(),er.min())>=.5-1e-10 and max(kr.max(),er.max())<=2+1e-10
    require(valid(x),'initial controlled state invalid');history=[];began=time.monotonic();primer=time.monotonic()
    for _ in range(3):evaluate(x)
    scan_cost=(time.monotonic()-primer)/3;covids=np.flatnonzero(cov);val=evaluate(x);best=x.copy();bestval=float(val.min());bestphi=max(0,target-bestval)+penalty(x)
    out=dict(version='H2eConditionedPhaseIV1',stage=a.stage,status='running',conditioning=checks,K0=K0,base_divergence_L2=math.sqrt(Q[0,0]),divergence_limit=Dmax,native_unique_observations=N,covariance_observations=len(covids),native_membership_shape=f['files']['membership']['shape'],seconds_per_full_tensor_scan=scan_cost,batch_indexed_error=err,context_load_condition_seconds=time.monotonic()-started,history=history)
    radius=.01;stagnant=0;maxiter=min(a.max_major or (20 if a.stage=='A' else 40),20 if a.stage=='A' else 40);seconds=300 if a.stage=='A' else 900;full_limit=180 if a.stage=='A' else 320;active_limit=3000 if a.stage=='A' else 6000;reason='iteration_budget';caps=set();scalarcaps=set();retained=set()
    def save(v,vals):
        state=base.copy();delta=v@Z[:size];state[free]+=delta[free]
        # Check every frozen bit before serializing, including signed zeros.
        fixed=np.ones(len(base),bool);fixed[free]=False;require(state[fixed].tobytes()==base[fixed].tobytes(),'frozen state bits changed')
        response,kr,er=bounds_metrics(v);j=int(np.argmin(vals));q=int(covids[j]);g,k,e=fields(v,[q])
        c=dict(version='OuterScalarInitializerV1',stage=a.stage,state=state.tolist(),state_sha256=state_hash(state),base_sha256=state_hash(base),conditioned_velocity_coefficients=v[:nv].tolist(),conditioned_scalar_coefficients=v[nv:].tolist(),min_lambda_over_K0=float(vals.min()),divergence_limit=Dmax,strong_divergence_L2=math.sqrt(max(0,quad(Q,v))),velocity_mass_distance=float(np.linalg.norm(v[:nv])),maximum_sampled_component_increment=float(np.max(abs(response))),K_envelope=[float(kr.min()),float(kr.max())],epsilon_envelope=[float(er.min()),float(er.max())],time_scale_ratio=[float((kr/er).min()),float((kr/er).max())],length_scale_ratio=[float((kr**1.5/er).min()),float((kr**1.5/er).max())],witness=dict(index=q,chart=obs[q,:3].tolist(),physical=physical[q].tolist(),gradient=g[0].tolist(),k=float(k[0]),epsilon=float(e[0]),inner=bool(inner[q])),context_sha256=sha(a.context),reference_immutable=True)
        atomic(a.output/'best-candidate.json',c)
    save(best,val)
    # A smooth constrained predictor uses the same conditioned coordinates and
    # physical divergence quadratic as every subsequent covariance trial.
    pred=minimize(lambda v:penalty(v),x,jac=pgrad,method='SLSQP',bounds=[(-.02,.02)]*nv+([(-.1,.1)]*ns if a.stage=='B' else []),constraints=[dict(type='ineq',fun=lambda v:(Dmax**2-quad(Q,v))/Dmax**2,jac=lambda v:np.r_[-2*(Q@np.r_[1,v[:nv]])[1:],np.zeros(size-nv)]/Dmax**2),dict(type='ineq',fun=lambda v:1-v[:nv]@v[:nv]/.12**2,jac=lambda v:np.r_[-2*v[:nv],np.zeros(size-nv)]/.12**2)],options=dict(maxiter=60,ftol=1e-11))
    for fraction in (1.,.5,.25,.125):
        trial=x+fraction*(pred.x-x)
        if valid(trial):
            pv=evaluate(trial)
            if pv.min()>val.min() and max(0,target-float(pv.min()))+penalty(trial)<bestphi:x=trial;val=pv;best=x.copy();bestval=float(val.min());bestphi=max(0,target-bestval)+penalty(x);save(best,val);break
    for iteration in range(maxiter):
        if time.monotonic()-began>=seconds or full>=full_limit or active_batches>=active_limit:reason='work_budget';break
        val=evaluate(x);minimum=float(val.min());before=max(0,target-minimum)+penalty(x)
        if minimum>=target:reason='sampled_interior_target';break
        selected=covids[np.argsort(val)[:64]];retained.update(map(int,selected));ids=np.array(sorted(retained),int)
        if len(ids)>192:
            vv=evaluate(x,ids);ids=ids[np.argsort(vv)[:192]]
        av=evaluate(x,ids);J=np.empty((len(ids),size));h=2e-6
        for b in range(size):
            vp=x.copy();vm=x.copy();vp[b]+=h;vm[b]-=h;J[:,b]=(evaluate(vp,ids)-evaluate(vm,ids))/(2*h)
        response,kr,er=bounds_metrics(x);flat=response.reshape(-1);near=np.flatnonzero(abs(flat)>.4);caps.update(map(int,near[np.argsort(-abs(flat[near]))[:128]]))
        if a.stage=='B':
            for block,r in enumerate((kr,er)):
                near=np.flatnonzero((r<.65)|(r>1.8));scalarcaps.update((block,int(i)) for i in near[:128])
        solution=None;newcap=0
        for separation in range(3):
            ci=np.array(sorted(caps),int);cr=V.reshape(-1,nv)[ci];cbase=flat[ci];si=sorted(scalarcaps)
            A=np.zeros((len(si),size));sbase=np.zeros(len(si))
            for j,(block,i) in enumerate(si):A[j,nv:]=(BK[i]/k0[i] if block==0 else BE[i]/e0[i]);sbase[j]=(kr if block==0 else er)[i]
            def append_zero(A):return np.column_stack([A,np.zeros(len(A))])
            def total(y):return x+y[:-1]
            def obj(y):return y[-1]+penalty(total(y))
            def jac(y):return np.r_[pgrad(total(y)),1.]
            cons=[dict(type='ineq',fun=lambda y:av+J@y[:-1]+y[-1]-target,jac=lambda y:np.column_stack([J,np.ones(len(J))])),dict(type='ineq',fun=lambda y:1-total(y)[:nv]@total(y)[:nv]/.12**2,jac=lambda y:np.r_[-2*total(y)[:nv],np.zeros(size-nv+1)]/.12**2),dict(type='ineq',fun=lambda y:(Dmax**2-quad(Q,total(y)))/Dmax**2,jac=lambda y:np.r_[-2*(Q@np.r_[1,total(y)[:nv]])[1:],np.zeros(size-nv+1)]/Dmax**2),dict(type='ineq',fun=lambda y:1-y[:-1]@y[:-1]/radius**2,jac=lambda y:np.r_[-2*y[:-1],0.]/radius**2)]
            if len(ci):
                C=np.column_stack([cr,np.zeros((len(cr),size-nv))]);cons += [dict(type='ineq',fun=lambda y:.5-cbase-C@y[:-1],jac=lambda y:-append_zero(C)),dict(type='ineq',fun=lambda y:.5+cbase+C@y[:-1],jac=lambda y:append_zero(C))]
            if len(si):cons +=[dict(type='ineq',fun=lambda y:2-sbase-A@y[:-1],jac=lambda y:-append_zero(A)),dict(type='ineq',fun=lambda y:sbase+A@y[:-1]-.5,jac=lambda y:append_zero(A))]
            # Solve in unit trust coordinates. The physical conditioned
            # amplitudes are small; presenting them unscaled to SLSQP made
            # its stopping tests insensitive to trust/volume constraints.
            scale=np.r_[np.full(size,radius),1.]
            scaled=[dict(type='ineq',fun=lambda w,c=c:c['fun'](w*scale),jac=lambda w,c=c:c['jac'](w*scale)*scale) for c in cons]
            solution=minimize(lambda w:obj(w*scale),np.r_[np.zeros(size),max(0,target-float(av.min()))],jac=lambda w:jac(w*scale)*scale,method='SLSQP',bounds=[(None,None)]*size+[(0,None)],constraints=scaled,options=dict(maxiter=70,ftol=1e-10))
            solution.x=solution.x*scale
            proposed=x+solution.x[:-1];rm,km,em=bounds_metrics(proposed);bad=np.flatnonzero(abs(rm.reshape(-1))>.5+1e-10);added=False
            for i in bad[np.argsort(-abs(rm.reshape(-1)[bad]))[:64]]:
                if int(i) not in caps:caps.add(int(i));added=True;newcap+=1
            if a.stage=='B':
                for block,r in enumerate((km,em)):
                    bads=np.flatnonzero((r<.5-1e-10)|(r>2+1e-10))
                    for i in bads[:64]:
                        if (block,int(i)) not in scalarcaps:scalarcaps.add((block,int(i)));added=True;newcap+=1
            if not added:break
        accepted=False;fraction_used=0;tv=val
        for fraction in (1.,.5,.25,.125,.0625):
            trial=x+fraction*solution.x[:-1]
            if not valid(trial):continue
            trialval=evaluate(trial);after=max(0,target-float(trialval.min()))+penalty(trial)
            if after<before-1e-8 and trialval.min()>=minimum-1e-5:x=trial;tv=trialval;accepted=True;fraction_used=fraction;break
        after=max(0,target-float(tv.min()))+penalty(x);oldradius=radius
        improvement=float(tv.min())-minimum
        if accepted:radius=min(.04,radius*1.25);stagnant=0 if improvement>=1e-5 else stagnant+1
        else:radius*=.5;stagnant+=1
        # Save immediately, including an improvement on the last major iteration.
        if float(tv.min())>bestval:best=x.copy();bestval=float(tv.min());bestphi=after;save(best,tv)
        history.append(dict(iteration=iteration,true_objective_before=before,true_objective_after=after,minimum_before=minimum,minimum_after=float(tv.min()),accepted=accepted,fraction=fraction_used,radius_before=oldradius,radius_after=radius,active_covariance=len(ids),active_response_caps=len(caps),active_scalar_caps=len(scalarcaps),new_cap_rows=newcap,local_solver_success=bool(solution.success),local_solver_message=str(solution.message),proposed_divergence_L2=math.sqrt(max(0,quad(Q,proposed))),proposed_velocity_L2=float(np.linalg.norm(proposed[:nv])),proposed_step_L2=float(np.linalg.norm(solution.x[:-1])),strong_divergence_L2=math.sqrt(max(0,quad(Q,x))),velocity_mass_distance=float(np.linalg.norm(x[:nv])),full_scans=full,active_derivative_batches=active_batches,seconds=time.monotonic()-began))
        out.update(history=history,best_minimum_lambda_over_K0=bestval);atomic(a.output/'progress.json',out)
        if stagnant>=6 or radius<1e-6:reason='stagnation';break
    # Always evaluate the terminal accepted state, even at a budget exit.
    terminal=evaluate(x)
    if terminal.min()>bestval:best=x.copy();bestval=float(terminal.min());save(best,terminal)
    best_scan=evaluate(best);save(best,best_scan)
    out.update(status='sampled_candidate' if bestval>=target else 'bounded_search_rejected',termination=reason,best_minimum_lambda_over_K0=bestval,terminal_minimum_lambda_over_K0=float(terminal.min()),full_scans=full,active_derivative_batches=active_batches,search_seconds=time.monotonic()-began,total_seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,classification='bounded nonconvex finite-basis search; independent native certification required; no global infeasibility certificate')
    atomic(a.output/'progress.json',out)
if __name__=='__main__':main()
