# SPDX-License-Identifier: MIT
"""Fixed-native-gradient standard covariance; signed scalar coordinates only.

The complete sampling sets are constraints, not an everywhere-positivity proof.
The phase-I dual bound certifies only this finite affine scalar space.
"""
import time,warnings
from pathlib import Path
import numpy as np
from scipy import sparse,linalg
from scipy.optimize import linprog,OptimizeWarning
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import sha,require

BLOCKS=dict(chart=(0,3),physical=(3,3),gradient=(6,9),k=(15,1),epsilon=(16,1),cold_k=(17,1),cold_epsilon=(18,1),inner=(19,1),covariance=(20,1),dependent=(21,1),weight=(22,1),inner_nu=(23,1),k_basis=(24,8),epsilon_basis=(32,8),k_error=(40,1),epsilon_error=(41,1),basis_error=(42,16),row_norm=(58,2),row_amplification=(60,2))
class Context:
 def __init__(self,path):
  self.path=Path(path);self.doc=load(path);self.meta=self.doc['context'];require(self.meta['version']=='FixedBScalarContextV1','unsupported scalar context')
  require(self.meta['model_specification']['name']=='StandardKEpsilon' and not self.meta['model_specification']['physical_regularization_active'],'wrong scalar model')
  receipt=self.doc.get('receipt',{});native=Path(receipt.get('output',''))
  require(native.is_file() and sha(native)==receipt.get('output_sha256'),'H2h exact native context receipt required')
  original=load(native);require(original['context']==self.meta,'H2h context differs from native export')
  for key in ('patch_signature','geometry_identity','inventory_identity','reference_sha256','defect_sha256'):require(self.doc[key]==original[key],'H2h native identity '+key)
  require(self.meta['model_specification']['coefficients']['Cmu']==.09,'H2h fixed standard Cmu differs')
  require(sha(Path(self.doc['basis_file']))==self.meta['basis_sha256']==self.doc['basis_sha256'],'H2h saved scalar basis differs')
  self.arrays={}
  for name,info in self.meta['files'].items():
   file=Path(self.doc['array_directory'])/info['file'];require(sha(file)==info['sha256'],'scalar context checksum '+name);self.arrays[name]=np.memmap(file,dtype=info['dtype'],mode='r',shape=tuple(info['shape']))
  require(self.arrays['points'].shape[1]==62,'scalar schema width')
 def block(self,name):
  start,width=BLOCKS[name];v=self.arrays['points'][:,start:start+width];return v[:,0] if width==1 else v
 def problem(self):
  outer=self.block('inner')==0;gradient=self.block('gradient')[outer].reshape(-1,3,3)
  return Problem(self.block('k')[outer],self.block('epsilon')[outer],self.block('cold_k')[outer],self.block('cold_epsilon')[outer],self.block('k_basis')[outer],self.block('epsilon_basis')[outer],np.linalg.eigvalsh(gradient+gradient.transpose(0,2,1))[:,-1],self.block('covariance')[outer]>0,indices=np.flatnonzero(outer))

class BudgetExhausted(RuntimeError):pass
class Budget:
 def __init__(self,rounds=24,scans=80,seconds=300):self.maxrounds=rounds;self.maxscans=scans;self.seconds=seconds;self.start=time.monotonic();self.rounds=0;self.scans=0
 def check(self):
  if time.monotonic()-self.start>=self.seconds:raise BudgetExhausted('optimizer wall-time budget')
 def scan(self):
  self.check()
  if self.scans>=self.maxscans:raise BudgetExhausted('complete separation-scan budget')
  self.scans+=1
 def round(self):
  self.check()
  if self.rounds>=self.maxrounds:raise BudgetExhausted('LP/constraint-generation round budget')
  self.rounds+=1
 def report(self):return dict(rounds=self.rounds,complete_scans=self.scans,seconds=time.monotonic()-self.start,limits=dict(rounds=self.maxrounds,scans=self.maxscans,seconds=self.seconds))

def dual_check(c,A,b,bounds,result):
 """Conservative finite-box dual bound, allowing nonzero stationarity error.

 y<=0 is dual feasible for A x<=b; lower/upper marginals are signed likewise.
 Residual*box subtraction yields a valid lower bound without asserting that a
 rounded floating-point dual vector is exactly stationary.
 """
 if result.x is None or not hasattr(result,'ineqlin'):return dict(verified=False,reason='no primal/dual solution')
 x=np.array(result.x);y=np.minimum(0,np.array(result.ineqlin.marginals));lo=np.maximum(0,np.array(result.lower.marginals));hi=np.minimum(0,np.array(result.upper.marginals))
 lower=np.array([v[0] for v in bounds]);upper=np.array([v[1] for v in bounds]);require(np.isfinite(lower).all() and np.isfinite(upper).all(),'dual check requires justified finite bounds')
 residual=c-A.T@y-lo-hi
 dual=float(b@y+lower@lo+upper@hi)
 stationarity_allowance=float(np.abs(residual)@np.maximum(abs(lower),abs(upper)))
 # Conservative accumulation allowance; data-rounding allowance is separate.
 arithmetic_allowance=128*np.finfo(float).eps*(float(abs(b)@abs(y)+abs(lower)@abs(lo)+abs(upper)@abs(hi))+stationarity_allowance+1)
 bound=dual-stationarity_allowance-arithmetic_allowance;primal=float(c@x);violation=float(max(0,np.max(A@x-b),np.max(lower-x),np.max(x-upper)))
 return dict(verified=bool(violation<=2e-9 and np.max(abs(residual))<=1e-8 and primal-bound<=1e-7),primal=primal,dual=dual,lower_bound=bound,gap=primal-bound,primal_max_violation=violation,stationarity_inf=float(np.max(abs(residual))),stationarity_allowance=stationarity_allowance,arithmetic_allowance=arithmetic_allowance,dual_signs_checked=True,dual_nonzero_rows=int(np.count_nonzero(y)),active_dual_rows=np.flatnonzero(abs(y)>1e-9).tolist(),active_dual_multipliers=y[abs(y)>1e-9].tolist())

def lp(c,A,b,bounds,budget):
 budget.round();start=time.monotonic()
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',category=OptimizeWarning,message='(Unrecognized options detected.*|Unknown solver options: threads)')
  result=linprog(c,A_ub=A,b_ub=b,bounds=bounds,method='highs-ds',options=dict(primal_feasibility_tolerance=1e-10,dual_feasibility_tolerance=1e-10,threads=1,time_limit=max(.01,budget.seconds-(time.monotonic()-budget.start))))
 receipt=dict(status=int(result.status),message=result.message,seconds=time.monotonic()-start,iterations=int(result.nit),rows=A.shape[0],columns=A.shape[1])
 if result.status==0:receipt['certificate']=dual_check(c,A,b,bounds,result)
 return result,receipt

class Problem:
 def __init__(self,k,e,kc,ec,Bk,Be,s,covariance=None,indices=None):
  self.k=np.array(k);self.e=np.array(e);self.kc=np.array(kc);self.ec=np.array(ec);self.n=len(k);self.K=np.zeros((self.n,16));self.E=self.K.copy();self.K[:,:8]=Bk;self.E[:,8:]=Be;self.s=np.array(s);self.b=.09*self.s;self.cov=np.ones(self.n,bool) if covariance is None else np.array(covariance,bool);self.indices=np.arange(self.n) if indices is None else indices
  require(np.min(self.kc)>0 and np.min(self.ec)>0,'invalid cold scalars');require(np.isfinite(self.K).all() and np.isfinite(self.E).all(),'invalid scalar response')
 def fields(self,a):return self.k+self.K@a,self.e+self.E@a
 def scan(self,a,delta,budget):
  budget.scan();k,e=self.fields(a);lam=(2/3)*k-self.b*k*k/e;g=delta-lam
  return dict(k=k,e=e,lambda_min=lam,g=g,minimum=float(np.min(lam[self.cov])),max_violation=float(max(0,np.max(g[self.cov]))),k_ratio=[float(np.min(k/self.kc)),float(np.max(k/self.kc))],epsilon_ratio=[float(np.min(e/self.ec)),float(np.max(e/self.ec))],worst=int(np.flatnonzero(self.cov)[np.argmin(lam[self.cov])]))
 def hard(self,buffer=0):
  Ck=self.K/self.kc[:,None];Ce=self.E/self.ec[:,None];A=np.vstack((Ck,-Ck,Ce,-Ce));b=np.concatenate((2-buffer-self.k/self.kc,self.k/self.kc-.5-buffer,2-buffer-self.e/self.ec,self.e/self.ec-.5-buffer))
  scale=np.maximum(1,np.maximum(np.max(abs(A),axis=1),abs(b)));return sparse.csr_matrix(A/scale[:,None]),b/scale,scale
 def coefficient_box(self):
  # Sixteen independent physical scalar rows imply this finite coefficient box.
  C=np.vstack((self.K/self.kc[:,None],self.E/self.ec[:,None]));lo=np.r_[.5-self.k/self.kc,.5-self.e/self.ec];hi=np.r_[2-self.k/self.kc,2-self.e/self.ec]
  _,R,pivot=linalg.qr(C.T,mode='economic',pivoting=True);rows=pivot[:16];Q=C[rows];require(np.linalg.matrix_rank(Q)==16,'saved sixteen scalar directions lost rank');N=np.linalg.inv(Q);eta=float(np.linalg.norm(np.eye(16)-N@Q,np.inf));require(eta<.01,'unstable implied coefficient box')
  bound=float(np.max(abs(N)@np.maximum(abs(lo[rows]),abs(hi[rows])))/(1-eta));bound=np.nextafter(bound*(1+1e-10)+1e-10,np.inf)
  return bound,dict(method='independent hard response rows; inverse residual bound',rank=16,condition=float(np.linalg.cond(Q)),inverse_residual_inf=eta,selected_rows=rows.tolist(),bound=bound)
 def covariance_rows(self,ratios,indices,delta,K0,positive=True):
  ids=np.array(indices);b=self.b[ids]
  if positive:
   r=np.array(ratios);ak=2*b*r-2/3;ae=-b*r*r;den=np.full(len(ids),K0)
  else:ak=b;ae=np.full(len(ids),-2/3);den=self.ec[ids]
  A=(ak[:,None]*self.K[ids]+ae[:,None]*self.E[ids])/den[:,None];rhs=-(ak*self.k[ids]+ae*self.e[ids]+(delta if positive else 0))/den
  # Scale complete rows, including the phase-I slack coefficient.
  scale=np.maximum(1,np.maximum(np.max(abs(A),axis=1),abs(rhs)))
  return sparse.csr_matrix(A/scale[:,None]),rhs/scale,-1/scale,scale

def solve(problem,K0,budget=None,save=None,rounding_allowance=0):
 budget=budget or Budget();P=problem;box,box_info=P.coefficient_box();hard,hb,hs=P.hard();n=16;rounds=[];cuts=[];out=dict(version='H2hScalarFeasibilityV1',coefficient_box=box_info,phase='zero_margin',rounds=rounds,classification='inconclusive')
 # The bound on optimal slack follows from the physical envelopes, even when
 # the scalar-feasible witness has not yet been found.
 tau_limit=float(max(1,np.max(np.maximum(P.b,0)*2*P.kc/P.ec)+2));bounds=[(-box,box)]*n+[(0,tau_limit)];c=np.r_[np.zeros(n),1.]
 def run(C,cb,slack,name):
  A=sparse.vstack((sparse.hstack((hard,sparse.csr_matrix((hard.shape[0],1)))),sparse.hstack((C,sparse.csr_matrix(slack[:,None])))),format='csr');rhs=np.r_[hb,cb]
  result,receipt=lp(c,A,rhs,bounds,budget);receipt['phase']=name;receipt['covariance_rows']=len(cb);receipt['hard_rows']=len(hb);rounds.append(receipt);out['budget']=budget.report()
  if save:save(out)
  return result,receipt,A,rhs
 try:
  ids=np.flatnonzero(P.cov & (P.s>0));C,cb,slack,scale=P.covariance_rows(None,ids,0,K0,False);result,rec,A,rhs=run(C,cb,slack,'zero_margin_covariance_only')
  if result.status!=0:out.update(termination='LP_numerical_or_hard_feasibility_unresolved',solver_status=int(result.status));return out,None
  a=result.x[:n];terminal=a.copy();scan=P.scan(a,0,budget);out['zero_margin']=dict(minimum_covariance=scan['minimum'],slack=float(result.x[-1]),k_ratio=scan['k_ratio'],epsilon_ratio=scan['epsilon_ratio'])
  certificate=rec['certificate'];certificate['native_rounding_allowance']=rounding_allowance;lower=certificate['lower_bound']-rounding_allowance
  if certificate['verified'] and lower>1e-8:
   out.update(classification='verified_limited_zero_margin_obstruction',termination='positive_verified_phase_I_lower_bound',certified_lower_bound=lower,terminal_coefficients=terminal.tolist());return out,terminal
  if result.x[-1]>1e-8:out.update(termination='phase_I_bound_not_certified',terminal_coefficients=terminal.tolist());return out,terminal
  # Protect nonpositive-strain rows explicitly: their viscosity term raises R.
  out['nonpositive_s_lower_margin']=float(np.min((2/3)*.5*P.kc[P.cov & (P.s<=0)])) if np.any(P.cov & (P.s<=0)) else None
  delta=1e-4*K0;target=delta+1e-8*K0;out['phase']='positive_margin';cut_ids=[];cut_ratios=[];best=None
  for iteration in range(budget.maxrounds-budget.rounds):
   scanned=P.scan(a,target,budget);envelope=max(.5-scanned['k_ratio'][0],scanned['k_ratio'][1]-2,.5-scanned['epsilon_ratio'][0],scanned['epsilon_ratio'][1]-2)
   if scanned['minimum']>=target and envelope<=1e-12:
    out.update(classification='sampled_positive_margin_candidate',termination='complete_native_and_GL20_separation_passed',terminal_coefficients=a.tolist(),minimum_covariance=scanned['minimum'],k_ratio=scanned['k_ratio'],epsilon_ratio=scanned['epsilon_ratio']);best=a.copy();return out,best
   violations=ids[np.argsort(scanned['g'][ids])[-min(256,len(ids)):]];violations=violations[scanned['g'][violations]>-1e-10*K0]
   if not len(violations):out.update(termination='unresolved_numerical_envelope',terminal_coefficients=a.tolist());return out,a
   cut_ids.extend(violations.tolist());cut_ratios.extend((scanned['k'][violations]/scanned['e'][violations]).tolist());C,cb,slack,scale=P.covariance_rows(cut_ratios,cut_ids,target,K0,True)
   # Positive phase-I uses original scalar bounds. A certified obstruction is
   # valid only for the physical target, so the buffered system is never used
   # to claim physical infeasibility.
   result,rec,A,rhs=run(C,cb,slack,'positive_margin_supporting_cuts')
   if result.status!=0:out.update(termination='positive_LP_numerical_failure',terminal_coefficients=a.tolist());return out,a
   a=result.x[:n];terminal=a.copy()
   if result.x[-1]>1e-7:
    C0,cb0,slack0,_=P.covariance_rows(cut_ratios,cut_ids,delta,K0,True);r0,rec0,A0,rhs0=run(C0,cb0,slack0,'unbuffered_positive_margin_diagnostic');terminal=r0.x[:n] if r0.status==0 else terminal
    if r0.status==0 and rec0['certificate']['verified'] and rec0['certificate']['lower_bound']-rounding_allowance>1e-8:
     out.update(classification='verified_limited_positive_margin_obstruction',termination='positive_verified_unbuffered_cut_relaxation_bound',certified_lower_bound=rec0['certificate']['lower_bound']-rounding_allowance,terminal_coefficients=terminal.tolist());return out,terminal
  out.update(termination='LP_round_budget',terminal_coefficients=terminal.tolist());return out,terminal
 except BudgetExhausted as e:
  out.update(termination=str(e));return out,locals().get('terminal')
 finally:out['budget']=budget.report()
