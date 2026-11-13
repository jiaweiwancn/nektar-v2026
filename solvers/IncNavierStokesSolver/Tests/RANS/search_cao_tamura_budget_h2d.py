#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded nonconvex velocity-only phase I in the native continuity nullspace."""
import argparse,json,os,time,math
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from scipy.optimize import minimize
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import state_hash
from validate_cls_epawf_m1 import indexed

from cao_tamura_h2d_state import load

def tensors(g,k,e,inner,nt,physical,us,kap):
    """Batch independent dimensional indexed CLS law, physical max branch."""
    S=g+g.transpose(0,2,1);O=g-g.transpose(0,2,1);eye=np.eye(3)
    safe=np.where(inner,1.,e);ss=np.einsum('nij,nij->n',S,S);oo=np.einsum('nij,nij->n',O,O)
    q=k/safe*np.sqrt(np.maximum(ss,oo)/2)
    with np.errstate(over='ignore',under='ignore'):
        cm=.3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(np.minimum(.75*q,700))))
    nu=cm*k*k/safe
    t1=S@S-ss[:,None,None]*eye/3;t2=O@S-S@O;t3=-(O@O)-oo[:,None,None]*eye/3
    t4=S@S@O-O@S@S
    Q=(k*nu/safe)[:,None,None]*(-.1*t1+.1*t2+.26*t3)
    C=(nu*k*k/safe**2)[:,None,None]*((-10*cm*cm)[:,None,None]*t4+(5*cm*cm*(oo-ss))[:,None,None]*S)
    z=np.where(inner,1.,physical[:,2]);sr=np.zeros_like(g);sr[:,0,2]=sr[:,2,0]=us/(kap*z);gamma=kap*us*z
    R=2*k[:,None,None]*eye/3-nu[:,None,None]*S+Q+.4*C-gamma[:,None,None]*(S-sr)
    R[inner]=2*k[inner,None,None]*eye/3-nt[inner,None,None]*S[inner]
    return (R+R.transpose(0,2,1))/2

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--setup',type=Path,required=True);p.add_argument('--preflight',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    doc=load(a.setup);f=doc['feasibility'];pre=json.loads(a.preflight.read_text());K0=pre['heights'][-1]['K0'];obs=f['observations'];base=np.array(f['base']);Z=np.array(f['basis']);gram=np.array(f['physical_mass_Gram'])
    g0=np.array([o['gradient'] for o in obs]);B=np.array([o['gradient_basis'] for o in obs]);V=np.array([o['velocity_basis'] for o in obs]);inner=np.array([o['inner'] for o in obs]);k=np.array([o['k'] for o in obs]);e=np.array([o['epsilon'] for o in obs]);nt=np.array([o['inner_nu_t'] for o in obs]);physical=np.array([o['physical'] for o in obs]);us=pre['u_tau'];kap=pre['kappa_CLS'];N=len(obs)
    # Validate batching against the separately indexed single-state formula.
    batch=tensors(g0,k,e,inner,nt,physical,us,kap);error=0
    for j in np.flatnonzero(~inner)[::max(1,int(np.sum(~inner))//30)]:
        o=indexed(g0[j],k[j],e[j],'physical');S=g0[j]+g0[j].T;sr=np.zeros((3,3));sr[0,2]=sr[2,0]=us/(kap*physical[j,2]);r=2*k[j]*np.eye(3)/3-float(o['nut'].real)*S+o['quadratic'].real+.4*o['cubic'].real-kap*us*physical[j,2]*(S-sr);error=max(error,float(np.max(abs(r-batch[j]))))
    assert error<1e-11,error
    trials=0;began=time.monotonic();history=[];best=np.zeros(32);best_value=-np.inf
    def evaluate(x,indices=None):
        nonlocal trials
        trials+=1;idx=np.arange(N) if indices is None else np.asarray(indices)
        g=g0[idx]+np.einsum('nbij,b->nij',B[idx],x)
        R=tensors(g,k[idx],e[idx],inner[idx],nt[idx],physical[idx],us,kap)
        values,vec=np.linalg.eigh(R)
        return values[:,0]/K0,vec[:,:,0],g
    t=time.monotonic()
    for _ in range(3):evaluate(np.zeros(32))
    cost=(time.monotonic()-t)/3
    out=dict(status='started',phase='velocity_only',native_observations=N,duplicate_membership_preserved_in_native_audit=True,provider='AnalyticalNeutralInflowV1',batch_indexed_error=error,seconds_per_full_tensor_scan=cost,basis_sizes=[16,32],K0=K0,budget=dict(iterations=40,trials=5000,seconds_per_basis=900,full_audits=6),history=history)
    atomic(a.output/'progress.json',out)
    limit=.15*.8
    for size in (16,32):
        phase_started=time.monotonic();phase_trials=trials;x=best.copy();radius=.01;stagnant=0
        # Linear smooth strain-extension predictor, with distance regularization.
        outer=np.flatnonzero(~inner);take=outer[::max(1,len(outer)//1500)]
        S=g0[take]+g0[take].transpose(0,2,1);SB=B[take,:size]+B[take,:size].transpose(0,1,3,2);sr=np.zeros_like(S);sr[:,0,2]=sr[:,2,0]=us/(kap*physical[take,2]);w=np.sqrt(kap*us*physical[take,2]/len(take))
        C=(SB*w[:,None,None,None]).transpose(0,2,3,1).reshape(-1,size);d=((sr-S)*w[:,None,None]).reshape(-1)
        predictor=np.zeros(32);predictor[:size]=np.linalg.lstsq(np.vstack([C,.1*np.eye(size)]),np.r_[d,np.zeros(size)],rcond=1e-10)[0]
        for _ in range(16):
            if math.sqrt(max(0,predictor@gram@predictor))<=limit and np.max(abs(np.einsum('nbj,b->nj',V,predictor)))<=.5:break
            predictor*=.5
        v,_,_=evaluate(x);pv,_,_=evaluate(predictor)
        if pv.min()>v.min():x=predictor
        retained=set();reason='iteration_budget'
        for iteration in range(40):
            if time.monotonic()-phase_started>900 or trials-phase_trials>5000:reason='work_budget';break
            val,eigen,g=evaluate(x);minimum=float(val.min());worst=int(val.argmin())
            if minimum>best_value:
                best_value=minimum;best=x.copy();state=base.copy();free=np.array(f['free_velocity_rows'],int);state[free]+=(best@Z)[free]
                atomic(a.output/'best-candidate.json',dict(state=state.tolist(),state_sha256=state_hash(state),base_sha256=doc['state_sha256'],coefficients=best.tolist(),min_lambda_over_K0=minimum,witness=obs[worst],initializer='NativeVelocityPhaseIV1',setup_sha256=__import__('hashlib').sha256(a.setup.read_bytes()).hexdigest()))
            if minimum>=1e-4:reason='sampled_interior_target';break
            active=np.argsort(val)[:48];retained.update(map(int,active));active=np.array(sorted(retained));
            # Bound active algebra without forgetting constraints in full scans.
            if len(active)>192:active=active[np.argsort(val[active])[:192]]
            J=np.empty((len(active),size));h=2e-6
            for b in range(size):
                plus=x.copy();minus=x.copy();plus[b]+=h;minus[b]-=h
                vp,_,_=evaluate(plus,active);vm,_,_=evaluate(minus,active);J[:,b]=(vp-vm)/(2*h)
            def penalty(y):
                z=x.copy();z[:size]+=y[:-1];res=C@z[:size]-d
                return y[-1]+1e-4*(z@gram@z)+1e-4*(res@res)
            def jac(y):
                z=x.copy();z[:size]+=y[:-1];return np.r_[2e-4*(gram[:size]@z+C.T@(C@z[:size]-d)),1.]
            constraints=[dict(type='ineq',fun=lambda y:val[active]+J@y[:-1]+y[-1]-1e-4,jac=lambda y:np.column_stack([J,np.ones(len(active))])),dict(type='ineq',fun=lambda y:limit**2-(x+np.r_[y[:-1],np.zeros(32-size)])@gram@(x+np.r_[y[:-1],np.zeros(32-size)]))]
            y0=np.r_[np.zeros(size),max(0.,1e-4-minimum)]
            solution=minimize(penalty,y0,jac=jac,bounds=[(-radius,radius)]*size+[(0,None)],constraints=constraints,method='SLSQP',options=dict(maxiter=60,ftol=1e-9))
            step=np.r_[solution.x[:-1],np.zeros(32-size)];accepted=False
            for fraction in (1.,.5,.25,.125,.0625):
                trial=x+fraction*step
                if trial@gram@trial>limit**2*(1+1e-6) or np.max(abs(np.einsum('nbj,b->nj',V,trial)))>.5:continue
                tv,_,_=evaluate(trial)
                if tv.min()>minimum+1e-7:x=trial;accepted=True;break
            row=dict(basis=size,iteration=iteration,minimum_lambda_over_K0=minimum,witness_index=worst,inner_min=float(val[inner].min()),outer_min=float(val[~inner].min()),radius=radius,accepted=accepted,linear_subproblem_success=bool(solution.success),mass_distance=float(math.sqrt(max(0,x@gram@x))),trials=trials,seconds=time.monotonic()-began)
            history.append(row);out.update(history=history,best_minimum_lambda_over_K0=best_value,trials=trials);atomic(a.output/'progress.json',out)
            if accepted:radius=min(.025,radius*1.2);stagnant=0 if float(tv.min())-minimum>=1e-5 else stagnant+1
            else:radius*=.5;stagnant+=1
            if stagnant>=8 or radius<1e-6:reason='stagnation';break
        out.setdefault('termination',[]).append(dict(basis=size,reason=reason,seconds=time.monotonic()-phase_started,trial_batches=trials-phase_trials));atomic(a.output/'progress.json',out)
        if best_value>=1e-4:break
    out.update(status='sampled_candidate' if best_value>=1e-4 else 'bounded_search_rejected',seconds=time.monotonic()-began,best_minimum_lambda_over_K0=best_value,trials=trials,classification='nonconvex finite-space search; native complete certificate required; no global infeasibility implication');atomic(a.output/'progress.json',out)
if __name__=='__main__':main()
