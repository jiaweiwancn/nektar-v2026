#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused convex math, signed coefficients, dual bounds and fresh gate tests."""
import argparse,time
from pathlib import Path
import numpy as np
import cao_tamura_h2h_scalar as engine
from types import SimpleNamespace
from cao_tamura_h2h_scalar import Problem,Budget,solve,BudgetExhausted,dual_check
from run_cao_tamura_budget_h2h import require_fresh
from continue_cls_epawf_pimple import atomic

def fixture(s=10.,blocked=False):
 n=9 if blocked else 8;K=np.zeros((n,8));K[:8]=np.eye(8);E=K.copy();ss=np.full(n,s);return Problem(np.ones(n),np.ones(n),np.ones(n),np.ones(n),K,E,ss)
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();tests={};start=time.monotonic()
 P=fixture();out,x=solve(P,1.,Budget(seconds=15));assert out['classification']=='sampled_positive_margin_candidate',out;assert np.any(x<0),'signed solution expected';tests['feasible_signed_coefficients']=out
 out,x=solve(fixture(blocked=True),1.,Budget(seconds=15));assert out['classification']=='verified_limited_zero_margin_obstruction',out;assert abs(out['certified_lower_bound']-(.9-2/3))<1e-7;tests['known_infeasible_verified_lower_bound']=out
 out,x=solve(P,1.,Budget(rounds=0,seconds=15));assert out['classification']=='inconclusive' and 'budget' in out['termination'];tests['budget_is_not_infeasibility']=True
 out,x=solve(fixture(s=(2/3)/.09,blocked=True),1.,Budget(seconds=15));assert out['classification']=='verified_limited_positive_margin_obstruction',out;tests['zero_feasible_positive_margin_obstructed']=out
 out,x=solve(fixture(s=-1),1.,Budget(seconds=15));assert out['classification']=='sampled_positive_margin_candidate',out;tests['nonpositive_s_handled_without_positive_s_cuts']=True
 original_lp=engine.lp
 def failed_lp(*args):return SimpleNamespace(status=4,x=None),dict(status=4,message='injected numerical failure')
 engine.lp=failed_lp
 try:
  out,x=solve(P,1.,Budget(seconds=15));assert out['classification']=='inconclusive' and 'numerical' in out['termination'];tests['numerical_failure_is_not_infeasibility']=True
 finally:engine.lp=original_lp
 # Independent tensor diagonalization on non-equilibrium tuples; S=G+G^T.
 rng=np.random.default_rng(721);G=rng.normal(size=(400,3,3));S=G+G.transpose(0,2,1);k=np.exp(rng.normal(size=400));e=np.exp(rng.normal(size=400));R=(2/3)*k[:,None,None]*np.eye(3)-(.09*k*k/e)[:,None,None]*S
 analytic=(2/3)*k-.09*k*k/e*np.linalg.eigvalsh(S)[:,-1];error=float(np.max(abs(analytic-np.linalg.eigvalsh(R)[:,0])));assert error<1e-12;tests['independent_standard_tensor_formula_max_error']=error
 r=k/e;b=.09*np.linalg.eigvalsh(S)[:,-1];positive=b>0;g=b*k*k/e-(2/3)*k+1e-4;kp=k*.93;ep=e*1.04;tangent=(2*b*r-2/3)*kp-b*r*r*ep+1e-4;exact=b*kp*kp/ep-(2/3)*kp+1e-4;assert np.all(tangent[positive]<=exact[positive]+1e-12);tests['supporting_cut_is_convex_lower_bound']=True
 model={'name':'StandardKEpsilon','K0':1.};base=dict(quadrature_order=24,state_sha256='candidate',cold_base_sha256='cold',model_specification=model,scalar_positive_finite=True,scalar_envelopes_passed=True,integrals_complete=True,minimum_covariance=2e-4);require_fresh(base,'candidate','cold',model)
 for key,value in [('quadrature_order',20),('state_sha256','other'),('cold_base_sha256','other'),('model_specification',{'name':'CLSHighRe'}),('scalar_envelopes_passed',False),('scalar_positive_finite',False),('integrals_complete',False),('minimum_covariance',0.)]:
  trial=dict(base);trial[key]=value
  try:require_fresh(trial,'candidate','cold',model);raise AssertionError('fresh gate accepted '+key)
  except ValueError:pass
 try:require_fresh(None,'candidate','cold',model);raise AssertionError('missing fresh audit')
 except ValueError:pass
 tests['fresh_GL24_is_required_for_issuance_and_launch']=True
 atomic(a.output,dict(version='H2hFocusedScalarChecksV1',passed=True,tests=tests,seconds=time.monotonic()-start))
if __name__=='__main__':main()
