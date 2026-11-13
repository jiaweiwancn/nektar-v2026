#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Replay the saved phase-I dual and localize its native scalar couplings.

No new LP or field tuning. Counts both complete diagnostic scans against H2h.
"""
import argparse,time,resource
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_h2h_scalar import Context
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,state_hash,sha

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('context','result','candidate','output'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();started=time.monotonic();previous=load(a.output) if a.output.exists() else {};ctx=Context(a.context);P=ctx.problem();result=load(a.result);candidate=load(a.candidate);x=np.array(candidate['scalar_increments']);K0=result['K0'];k,e=P.fields(x);kb,eb=P.fields(np.zeros(16));m=ctx.arrays['members'];families=ctx.meta['families'];ids=np.flatnonzero(P.cov&(P.s>0));N=P.n
 def point(i):
  record=int(P.indices[i]);member=m[record];return dict(record=record,set='native' if member[1]==0 else 'GL20',family=families[int(member[3])] if member[3]>=0 else 'physical_volume_GL20',geometry=int(member[4]),element=int(member[5]),radial_side=int(member[6]),tangential_side=int(member[7]),block=int(member[8]),dependent_epsilon=bool(member[9]),chart=ctx.block('chart')[record].tolist(),physical=ctx.block('physical')[record].tolist(),B_k=float(kb[i]),B_epsilon=float(eb[i]),cold_k=float(P.kc[i]),cold_epsilon=float(P.ec[i]),k=float(k[i]),epsilon=float(e[i]),k_over_cold=float(k[i]/P.kc[i]),epsilon_over_cold=float(e[i]/P.ec[i]),s_max=float(P.s[i]),minimum_covariance=float((2/3)*k[i]-P.b[i]*k[i]*k[i]/e[i]))
 H,hb,_=P.hard();C,cb,slack,_=P.covariance_rows(None,ids,0,K0,False);cert=result['rounds'][0]['certificate'];rows=cert['active_dual_rows'];y=np.array(cert['active_dual_multipliers']);A=[];rhs=[];details=[]
 for row,multiplier in zip(rows,y):
  if row<4*N:
   block=row//N;i=row%N;kind=['k_upper','k_lower','epsilon_upper','epsilon_lower'][block];v=np.r_[H.getrow(row).toarray()[0],0.];b=hb[row]
  else:
   j=row-4*N;i=int(ids[j]);kind='zero_covariance';v=np.r_[C.getrow(j).toarray()[0],slack[j]];b=cb[j]
  A.append(v);rhs.append(b);details.append(dict(kind=kind,LP_row=row,dual_multiplier=float(multiplier),scaled_row=v.tolist(),scaled_rhs=float(b),point=point(i)))
 A=np.array(A);rhs=np.array(rhs);stationarity=A.T@y-np.r_[np.zeros(16),1.];bound=float(rhs@y);box=result['coefficient_box']['bound'];tauupper=float(max(1,np.max(np.maximum(P.b,0)*2*P.kc/P.ec)+2));allowance=float(abs(stationarity)@np.r_[np.full(16,box),tauupper]);require(np.max(abs(stationarity))<1e-8 and bound-allowance>1e-8,'independent active dual replay failed')
 def summary(k,e):
  lam=(2/3)*k-P.b*k*k/e;out={}
  for setid,label in [(0,'native'),(1,'GL20')]:
   mask=(m[P.indices,1]==setid);cov=mask&P.cov;j=int(np.flatnonzero(cov)[np.argmin(lam[cov])]);out[label]=dict(minimum_covariance=float(lam[j]),normalized_covariance=float(lam[j]/K0),witness=point(j),k_ratio=[float(np.min(k[mask]/P.kc[mask])),float(np.max(k[mask]/P.kc[mask]))],epsilon_ratio=[float(np.min(e[mask]/P.ec[mask])),float(np.max(e[mask]/P.ec[mask]))],negative_covariance_memberships=int(np.count_nonzero(lam[cov]<0)))
  return out
 baseline=summary(kb,eb);terminal=summary(k,e)
 for label,fields in ((baseline,(kb,eb)),(terminal,(k,e))):
  for setid,entry in label.items():
   i=int(np.flatnonzero(P.indices==entry['witness']['record'])[0]);kk,ee=fields;entry['witness'].update(k=float(kk[i]),epsilon=float(ee[i]),k_over_cold=float(kk[i]/P.kc[i]),epsilon_over_cold=float(ee[i]/P.ec[i]),minimum_covariance=entry['minimum_covariance'])
 # Tests whether even arbitrary pointwise scalar choices inside the original
 # envelope could close a row. This is a diagnostic, not a new search.
 optimistic_zero=(2/3)*2*P.ec-P.b*.5*P.kc;worst=int(ids[np.argmin(optimistic_zero[ids]/P.ec[ids])]);frozen=np.max(abs(P.K),axis=1)+np.max(abs(P.E),axis=1)==0;frozenbad=np.flatnonzero(frozen&P.cov&((2/3)*eb-P.b*kb<0));fixed=[]
 for i in frozenbad[np.argsort(((2/3)*eb-P.b*kb)[frozenbad])[:12]]:fixed.append(point(int(i)))
 # Independent tensor contractions at the unchangeable matching witness,
 # retaining every distinct one-sided native row and actual/dependent context.
 paired=[]
 if fixed:
  target=np.asarray(fixed[0]['chart']);all_chart=ctx.block('chart');matching=np.flatnonzero(np.max(abs(all_chart-target),axis=1)<1e-14);seen=set()
  for j in matching:
   membership=m[j];key=(int(membership[4]),int(membership[5]),int(membership[6]),int(membership[7]),bool(membership[9]))
   if key in seen:continue
   seen.add(key);G=ctx.block('gradient')[j].reshape(3,3);S=G+G.T;kk=float(ctx.block('k')[j]);ee=float(ctx.block('epsilon')[j]);inner=bool(ctx.block('inner')[j]);nu=float(ctx.block('inner_nu')[j]) if inner else .09*kk*kk/ee;R=(2/3)*kk*np.eye(3)-nu*S;ev,V=np.linalg.eigh(R)
   paired.append(dict(geometry=key[0],element=key[1],radial_side=key[2],tangential_side=key[3],dependent_epsilon=key[4],inner=inner,chart=all_chart[j].tolist(),physical=ctx.block('physical')[j].tolist(),k=kk,epsilon=ee,eddy_viscosity=nu,gradient=G.tolist(),divergence=float(np.trace(G)),isotropic_tensor=((2/3)*kk*np.eye(3)).tolist(),linear_tensor=(-nu*S).tolist(),tensor=R.tolist(),eigenvalues=ev.tolist(),critical_eigenvector=V[:,0].tolist(),eigenvector_residual=float(np.max(abs(R@V[:,0]-ev[0]*V[:,0])))))
 # Unique native evaluations retain element, sides and exact binary64 location.
 unique={}
 for d in details:
  q=d['point'];key=(q['set'],q['geometry'],q['element'],q['radial_side'],q['tangential_side'],q['block'],q['dependent_epsilon'],np.asarray(q['chart'],dtype='<f8').tobytes());unique[key]=q
 report=dict(version='H2hActiveConstraintReplayV1',candidate_sha256=candidate['state_sha256'],B_sha256=ctx.meta['state_sha256'],cold_sha256=ctx.meta['base_sha256'],context_sha256=sha(a.context),active_rows=details,independent_dual=dict(nonzero_rows=len(rows),row_combination=bound,stationarity_inf=float(np.max(abs(stationarity))),finite_box_residual_allowance=allowance,lower_bound_before_data_rounding=bound-allowance,coefficient_bounds_inactive=bool(np.max(abs(x))<box)),B=baseline,terminal=terminal,optimistic_pointwise_zero_margin=dict(minimum_normalized=float(optimistic_zero[worst]/P.ec[worst]),witness=point(worst),explanation='k=.5*cold, epsilon=2*cold independently at each positive-s point; ignores trace/basis coupling'),unchangeable_negative_rows=dict(count=len(frozenbad),unique_contexts=len({(int(m[P.indices[i],4]),int(m[P.indices[i],5]),int(m[P.indices[i],6]),int(m[P.indices[i],7]),bool(m[P.indices[i],9]),np.asarray(ctx.block('chart')[P.indices[i]],dtype='<f8').tobytes()) for i in frozenbad}),examples=fixed),independent_matching_tensor_replay=paired,optimizer_diagnostic_work=dict(additional_LP_rounds=0,additional_complete_scans=2,seconds=time.monotonic()-started,previous_diagnostic_scans=previous.get('optimizer_diagnostic_work',{}).get('aggregate_complete_scans',result['budget']['complete_scans'])-result['budget']['complete_scans'],aggregate_LP_rounds=result['budget']['rounds'],aggregate_complete_scans=previous.get('optimizer_diagnostic_work',{}).get('aggregate_complete_scans',result['budget']['complete_scans'])+2,aggregate_optimizer_seconds=previous.get('optimizer_diagnostic_work',{}).get('aggregate_optimizer_seconds',result['budget']['seconds'])+time.monotonic()-started),physical_attempts=0,accepted_physical_steps=0)
 require(report['optimizer_diagnostic_work']['aggregate_complete_scans']<=80 and report['optimizer_diagnostic_work']['aggregate_optimizer_seconds']<=300,'diagnostic work exceeded budget');atomic(a.output,report)
if __name__=='__main__':main()
