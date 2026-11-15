#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Joint LMI equivalence, necessary cuts, physical caps and fail-fast driver."""
import argparse,tempfile,time,json
from pathlib import Path
import numpy as np
from cao_tamura_h2i_convex import JointProblem,phase,convex_segment
from cao_tamura_h2h_scalar import Budget
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import state_hash
import run_cao_tamura_budget_h2f_startup as startup
from run_cao_tamura_budget_h2i import fresh_gates,COLD_HASH

NAMES=['coarse-h0','coarse8','fine-h0','fine16','coarse4','coarse4-restart','post-assignment-reject']
class FakeRun(startup.StartupRun):
 def __init__(self,root,fail=None,kind='rejection'):
  self.root=root;root.mkdir();self.calls=[];self.fail=fail;self.kind=kind;self.cert={'state_sha256':state_hash([1.])};self.out={'results':{}};self.cps={}
 def save(self):pass
 def checkpoint_steps(self,path):return self.cps[str(path)]
 def validate_checkpoint(self,d,expected,dt):
  assert self.cps[d['checkpoint']]==expected
  if d.get('bad_checkpoint'):raise ValueError('injected checkpoint failure')
 def native(self,name,dt,steps=0,checkpoint=None,reject=False):
  self.calls.append(name);count=(self.checkpoint_steps(checkpoint) if checkpoint else 0)+(0 if reject else steps);cp=self.root/(name+'-cp');cp.mkdir();(cp/'State.bin').write_bytes(b'unchanged fake parent');self.cps[str(cp)]=count
  attempts=[dict.fromkeys(('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls'),1) for _ in range(steps)]
  for v in attempts:v.update(scheduled_complete=True,solved_epsilon_retained_bitwise=True)
  d=dict(failure='deliberate post-assignment explicit rejection' if reject else None,state=[1.],fixed_reference=[1.],fixed_defect=[0.],physical_samples={},steps=count,checkpoint=str(cp),initialization_lift_calls_iterations=[0,0],pimple_attempts=attempts,implicit_callback_count=steps)
  report=self.root/(name+'.json');report.write_text('{}');r=dict(code=1 if reject else 0,report=str(report),report_parse_error=None)
  if name==self.fail:
   if self.kind=='rejection':d['failure']='injected numerical rejection'
   elif self.kind=='status':r['code']=2
   elif self.kind=='missing_report':report.unlink()
   elif self.kind=='invalid_report':r['report_parse_error']='injected parse failure'
   elif self.kind=='steps':d['steps']=count+1
   elif self.kind=='checkpoint':d['bad_checkpoint']=True
   elif self.kind=='counts':d['pimple_attempts']=[] if steps else [{}]
   elif self.kind=='epsilon':d['pimple_attempts'][0]['solved_epsilon_retained_bitwise']=False
  self.out['results'][name]=r;return d

def fixture(inner_bad=False,repair=False,between=False):
 n=9;nv=2;g=np.zeros((n,4));k=np.full(n,.02);e=np.full(n,.01);inner=np.zeros(n,bool);inner[-1]=True;nt=np.zeros(n);nt[-1]=1e-4;V=np.zeros((n,nv,6));E=np.zeros((n,8));E[:8]=np.eye(8)*.01
 if inner_bad:g[-1,0]=100
 if between:g[-1,0]=(2*.02/3-1.5e-4*.02)/(2*1e-4)
 if repair:g[0,0]=100;V[0,0,0]=1000
 M=np.diag([0.,1.,1.]);D=np.zeros((3,3))
 return JointProblem(g,np.zeros((n,2)),np.zeros((n,2)),k,e,k,e,nt,inner,np.ones(n,bool),V,E,M,[D,D,D],.02)

def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();tests={};start=time.monotonic()
 P=fixture();r=phase(P,2e-4*P.K0,[.501,1.999],Budget(40,120,20),lambda x:None);assert r['termination']=='seen_feasible';tests['known_feasible_zero_increment']=True
 P=fixture(inner_bad=True);r=phase(P,1e-4*P.K0,[.5,2],Budget(40,120,20),lambda x:None);assert r['termination']=='verified_relaxation_obstruction',r;tests['known_fixed_inner_obstruction']=r['proof']['receipt']['certificate']
 P=fixture(between=True);b=phase(P,2e-4*P.K0,[.501,1.999],Budget(40,120,20),lambda x:None);u=phase(P,1e-4*P.K0,[.5,2],Budget(40,120,20),lambda x:None);assert b['termination']=='verified_relaxation_obstruction' and u['termination']=='seen_feasible';tests['buffer_obstruction_does_not_certify_actual_infeasibility']=True
 P=fixture(repair=True);r=phase(P,2e-4*P.K0,[.501,1.999],Budget(40,120,20),lambda x:None);assert r['termination']=='seen_feasible',r;assert np.min(r['best_coefficients'])<0;tests['signed_joint_repair_with_true_mass_cap']=r['best_summary']
 P=fixture(inner_bad=True);r=phase(P,1e-4*P.K0,[.5,2],Budget(0,120,20),lambda x:None);assert r['termination'].startswith('budget_exhausted');tests['budget_is_inconclusive']=True
 # Nonzero trace and nonsymmetric G: compare the exact transformed tensor.
 rng=np.random.default_rng(809);P=fixture();P.g=rng.normal(size=P.g.shape);x=rng.normal(size=P.dim)*.001;delta=1e-4*P.K0;s=P.scan(x,delta,Budget(seconds=10));G=np.zeros((P.n,3,3));G[:,0,0]=s['gradient'][:,0];G[:,0,2]=s['gradient'][:,1];G[:,2,0]=s['gradient'][:,2];G[:,2,2]=s['gradient'][:,3];S=G+G.transpose(0,2,1);nu=np.where(P.inner,P.nt,.09*P.k**2/s['epsilon']);R=2*P.k[:,None,None]/3*np.eye(3)-nu[:,None,None]*S;ev=np.linalg.eigvalsh(R)[:,0];assert np.max(abs(ev-s['lambda_covariance']))<1e-14
 for i in range(P.n):
  v=rng.normal(size=3);v/=np.linalg.norm(v);row,rhs=P.rayleigh_cut(i,v,delta);M=(R[i]-delta*np.eye(3))/P.K0 if P.inner[i] else (R[i]-delta*np.eye(3))*s['epsilon'][i]/(P.K0*P.ec[i]);assert abs(rhs-row@np.r_[x,0]-v@M@v)<1e-12
 tests['full_S_no_deviatoric_projection_and_affine_Rayleigh']=True
 Q=np.array([[0.,0.,0.],[0.,2.,.3],[0.,.3,1.]]);row,rhs=P.quadratic_cut(Q,.02,np.r_[.02,-.01,np.zeros(8)]);assert row@np.zeros(P.dim+1)<=rhs+1e-15;tests['quadratic_tangent_is_necessary_not_sufficient']=True
 candidate=dict(state_sha256='candidate',divergence_limit=.06030809942199109);model=dict(K0=.02)
 fresh=dict(quadrature_order=28,state_sha256='candidate',cold_base_sha256=COLD_HASH,model_specification=model,inner_minimum_covariance=3e-6,outer_minimum_covariance=3e-6,scalar_positive_finite=True,scalar_envelopes_passed=True,integrals_complete=True,velocity_change_L2=.01,maximum_sampled_velocity_increment=.1,strong_divergence_L2=.01)
 assert all(fresh_gates(fresh,candidate,model).values())
 for key,value in dict(quadrature_order=24,state_sha256='different',cold_base_sha256='different',model_specification=dict(K0=.03),inner_minimum_covariance=1e-6,outer_minimum_covariance=-.001,scalar_envelopes_passed=False,integrals_complete=False,velocity_change_L2=.13,maximum_sampled_velocity_increment=.6,strong_divergence_L2=.07).items():
  assert not all(fresh_gates(dict(fresh,**{key:value}),candidate,model).values()),key
 tests['fresh_GL28_identity_margin_envelope_and_physical_caps_required']=True
 original=startup.same_history;startup.same_history=lambda a,b:a['state']==b['state'];driver=[]
 try:
  with tempfile.TemporaryDirectory() as temp:
   root=Path(temp)
   for kind in ('rejection','status','missing_report','invalid_report','steps','checkpoint','counts','epsilon'):
    for stop in NAMES[:-1]:
     if kind=='epsilon' and stop.endswith('h0'):continue
     f=FakeRun(root/(kind+'-'+stop),stop,kind)
     try:f.pair();raise AssertionError('driver continued after '+kind+' '+stop)
     except (ValueError,AssertionError):
      assert f.calls==NAMES[:NAMES.index(stop)+1],(kind,stop,f.calls)
     driver.append(dict(kind=kind,stop=stop,calls=f.calls))
   f=FakeRun(root/'success');f.pair();assert f.calls==NAMES and f.out['status']=='completed_bounded_pair'
   f=FakeRun(root/'rollback-unexpected','post-assignment-reject');
   try:f.pair();raise AssertionError('unexpected rollback rejection accepted')
   except ValueError:assert f.calls==NAMES
 finally:startup.same_history=original
 tests['driver_failure_prefixes']=driver;tests['deliberate_rollback_distinct_and_successful_order']=True
 atomic(a.output,dict(version='H2iFocusedChecksV1',passed=True,tests=tests,seconds=time.monotonic()-start,scope='math and synthetic driver control; no native governing solves'))
if __name__=='__main__':main()
