#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded fixed-k joint initializer, seen/fresh audit ledger and exact V2 pilot."""
import argparse,os,time,resource,hashlib,json
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_initializer_context import NativeContext
from cao_tamura_budget_h2g import case
from cao_tamura_h2i_convex import Context,JointProblem,phase,summary
from cao_tamura_h2h_scalar import Budget,BudgetExhausted
from run_cao_tamura_budget_h2h import runtime_identity,same_runtime
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import require,sha,state_hash,THREADS
from run_cao_tamura_budget_h2e import certificate_v2
from run_cao_tamura_budget_h2f_startup import StartupRun
from cao_tamura_startup_receipt import divergence_gate,complete_receipt,validate_receipt
from datetime import datetime
from zoneinfo import ZoneInfo
B_HASH='c83da50b65d2e9980325d365779a9d53394fbf10e6431d1b87f79728366ae13d'
COLD_HASH='038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c'
def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def fresh_gates(v,c,model):
 K0=model['K0'];return dict(identity=v.get('quadrature_order')==28 and v.get('state_sha256')==c['state_sha256'] and v.get('cold_base_sha256')==COLD_HASH and v.get('model_specification')==model,covariance=min(v.get('inner_minimum_covariance',-float('inf')),v.get('outer_minimum_covariance',-float('inf')))>=1e-4*K0-1e-12,scalar=bool(v.get('scalar_positive_finite') and v.get('scalar_envelopes_passed') and v.get('integrals_complete')),velocity=v.get('velocity_change_L2',float('inf'))<=.12*(1+1e-10) and v.get('maximum_sampled_velocity_increment',float('inf'))<=.5+1e-10,divergence=v.get('strong_divergence_L2',float('inf'))<=c['divergence_limit']+1e-8)
def prepare(a):
 start=time.monotonic();old=NativeContext(a.parent/'standard-context/context.json');B=load(a.parent/'standard-B/best-candidate.json');require(state_hash(B['state'])==B_HASH and old.meta['base_sha256']==COLD_HASH,'joint B/cold identity');T,Ts,condition=old.transforms();require(T.shape==(41,41) and Ts.shape==(16,16),'qualified joint ranks differ');raw=np.array(old.array('basis'));Z=np.vstack((T.T@raw[:41],Ts.T@raw[41:]));a.output.mkdir(exist_ok=False);Z.astype('<f8').tofile(a.output/'basis.f64');np.asarray(T,dtype='<f8').tofile(a.output/'T.f64');np.asarray(Ts,dtype='<f8').tofile(a.output/'Ts.f64');arrays=a.output/'arrays';arrays.mkdir();xml=case(a.output/'case',a.h1,load(a.preflight),flat=False)
 packet=dict(operation='h2i_context',expected_model='StandardKEpsilon',inner_stress_model='Linear',xml_sha256=sha(xml),state=B['state'],state_sha256=B_HASH,base_sha256=COLD_HASH,basis_file=str(a.output/'basis.f64'),basis_sha256=sha(a.output/'basis.f64'),B_scalar_coefficients=B['conditioned_scalar_coefficients'],context_directory=str(arrays))
 _,r=run_probe(a.probe,xml,packet,a.output/'native');d=load(Path(r['output']));d.update(array_directory=str(arrays),basis_file=str(a.output/'basis.f64'),basis_sha256=sha(a.output/'basis.f64'),conditioning=condition,transform_sha256=dict(velocity=sha(a.output/'T.f64'),scalar=sha(a.output/'Ts.f64')),parent_context=str(old.path),parent_context_sha256=sha(old.path),parent_B_file=str(a.parent/'standard-B/best-candidate.json'),receipt=r,runtime=runtime_identity(a.probe),setup_seconds=time.monotonic()-start,created_Asia_Shanghai=now());atomic(a.output/'context.json',d)

def form_anchor(ctx,B,cold,theta,label):
 s=B.copy();rows=ctx.meta['rows']
 for key in ('k','epsilon'):
  ids=np.array(rows[key],int);t=0 if key=='k' and label=='B' else theta;s[ids]=cold[ids] if t==0 else cold[ids]+t*(B[ids]-cold[ids])
 return s

def make_candidate(ctx,anchor,B,Z,coeff,label,theta):
 s=anchor.copy();rows=ctx.meta['rows'];x=np.array(coeff);nv=ctx.meta['direction_counts']['velocity']
 if np.any(x[:nv]):ids=np.array(rows['velocity'],int);s[ids]=B[ids]+x[:nv]@Z[:nv,ids]
 if np.any(x[nv:]):ids=np.array(rows['epsilon'],int);s[ids]=anchor[ids]+x[nv:]@Z[nv+8:,ids]
 koffset=ctx.meta['native_coefficient_offsets'][2];ksize=ctx.meta['native_coefficient_sizes'][2];k=s[koffset:koffset+ksize];require(k.tobytes()==anchor[koffset:koffset+ksize].tobytes(),'fixed-k reconstruction changed bits')
 allowed=set(rows['velocity']+rows['k']+rows['epsilon']);fixed=np.array([i for i in range(len(s)) if i not in allowed]);require(s[fixed].tobytes()==B[fixed].tobytes(),'joint reconstruction changed protected bits')
 return dict(version='H2iFixedKCandidateV1',state=s.tolist(),state_sha256=state_hash(s),base_sha256=COLD_HASH,B_sha256=B_HASH,anchor=label,theta=theta,anchor_k=k.tolist(),anchor_k_sha256=state_hash(k),basis_sha256=ctx.doc['basis_sha256'],context_sha256=sha(ctx.path),joint_coefficients=x.tolist(),divergence_limit=.06030809942199109,pressure_and_prescribed_bits=True)

def search(a):
 setup=time.monotonic();ctx=Context(a.context);Bdoc=load(a.parent/'standard-B/best-candidate.json');B=np.array(Bdoc['state']);cold=np.array(NativeContext(a.parent/'standard-context/context.json').array('base')[0]);nv=ctx.meta['direction_counts']['velocity'];Z=np.memmap(ctx.doc['basis_file'],dtype='<f8',mode='r',shape=(nv+16,len(B)));require(sha(Path(ctx.doc['basis_file']))==ctx.doc['basis_sha256'],'joint basis changed');theta,retention=ctx.retention();retention=ctx.doc.get('saved_retention',retention);theta=retention['theta'];a.output.mkdir(exist_ok=False);budget=Budget(40,120,600);budget.scans=1;out=dict(version='H2iJointSearchV1',retention=retention,phases=[],classification='inconclusive',created_Asia_Shanghai=now(),setup_seconds=time.monotonic()-setup);allanchors={label:form_anchor(ctx,B,cold,theta,label) for label in ('A','B')};koff=ctx.meta['native_coefficient_offsets'][2];ksize=ctx.meta['native_coefficient_sizes'][2]
 anchors={label:dict(state_sha256=state_hash(s),k_sha256=state_hash(s[koff:koff+ksize]),theta=theta if label=='A' else 0,epsilon_theta=theta,scalar_coefficients='same archived scalar space, coefficient-level cold/B retraction') for label,s in allanchors.items()};atomic(a.output/'anchors.json',dict(retention=retention,anchors=anchors,B_hash=B_HASH,cold_hash=COLD_HASH))
 # Confirm old scalar-span membership without reconstructing B velocity.
 old_scalar=np.array(Bdoc['conditioned_scalar_coefficients']);res=(B-cold)-old_scalar@Z[nv:];ids=np.array(ctx.meta['rows']['k']+ctx.meta['rows']['epsilon']);out['anchor_span_roundoff_max']=float(np.max(abs(res[ids])));require(out['anchor_span_roundoff_max']<1e-12,'anchor outside archived scalar span')
 selected=None;event_count=0
 def save():out['budget']=budget.report();atomic(a.output/'result.json',out)
 try:
  for label in ('A','B'):
   if label=='B':require(out.get('A_actual_obstruction') is True,'anchor B requires actual-limit A obstruction')
   P=ctx.problem(label,theta);out.setdefault('bounds',{})[label]=P.bound_report
   for buffered in (True,False):
    tag=label+('-buffered' if buffered else '-actual');folder=a.output/tag;folder.mkdir();delta=(2 if buffered else 1)*1e-4*P.K0;band=[.501,1.999] if buffered else [.5,2]
    def persist(event):
     nonlocal event_count
     event_count+=1;atomic(folder/('event-%03d.json'%event_count),event)
     if event.get('event')=='best_hard_diagnostic':atomic(folder/'best-candidate.json',make_candidate(ctx,allanchors[label],B,Z,event['coefficients'],label,theta))
    r=phase(P,delta,band,budget,persist);out['phases'].append(dict(anchor=label,buffered=buffered,result=r));save()
    if r['best_coefficients'] is not None:
     candidate=make_candidate(ctx,allanchors[label],B,Z,r['best_coefficients'],label,theta);candidate['optimizer_summary']=r['best_summary'];candidate['phase']=tag;atomic(a.output/'candidate.json',candidate);selected=candidate
    if r['termination']=='seen_feasible':out.update(classification='seen_feasible_candidate',selected_anchor=label,selected_phase=tag);break
    if r['termination']!='verified_relaxation_obstruction':out.update(classification='inconclusive',termination=r['termination']);break
    if not buffered:out[label+'_actual_obstruction']=True
   if out['classification']=='seen_feasible_candidate' or not out.get(label+'_actual_obstruction'):break
   if label=='A':
    # No blind anchor replay: only proof rows known to be identical can carry
    # the obstruction. Otherwise the second predeclared anchor is solved.
    proof=out['phases'][-1]['result']['proof'];outer=[m for m in proof['cut_metadata'] if m.get('kind')=='PSD' and not m.get('inner')]
    if not outer:
     out.update(classification='verified_both_anchors_obstruction',B_obstruction_replayed_from_A=True,replay_scope='all nonzero dual PSD cuts are inner; k, epsilon origins and physical caps identical');break
   else:out['classification']='verified_both_anchors_obstruction'
 except BudgetExhausted as e:out.update(classification='inconclusive',termination=str(e))
 finally:
  out.update(budget=budget.report(),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,native_candidate_checks=0,physical_attempts=0,selected_state_sha256=selected['state_sha256'] if selected else None);save()

def repeatable(a):
 c=load(a.candidate);ctx=Context(a.context);pre=load(a.preflight);a.output.mkdir(exist_ok=False);ledger=a.output.parent/'candidate-checks.json';entries=load(ledger) if ledger.exists() else [];require(len(entries)<3,'repeatable candidate budget');entries.append(dict(state_sha256=c['state_sha256'],started=now(),kind='seen_native_only',orders=[20,24]));atomic(ledger,entries);xml=case(a.output/'case',a.h1,pre,flat=False);packet=dict(operation='h2i_candidate',inner_stress_model='Linear',expected_model='StandardKEpsilon',xml_sha256=sha(xml),state=c['state'],state_sha256=c['state_sha256'],base_sha256=COLD_HASH,coupled=True,anchor_k=c['anchor_k'],anchor_k_sha256=c['anchor_k_sha256']);start=time.monotonic();_,rr=run_probe(a.probe,xml,packet,a.output/'native');d=load(Path(rr['output']));delta=1e-4*pre['models'][1]['K0'];au=d['audit'];seen=d['seen_sampling'];gates=dict(native=au['accepted_awf_hard_admissible']==[1] and au['admissible']==[1],native_margin=min(au['accepted_awf_covariance'])>=delta-1e-12,budget=au['accepted_awf_scalar_domain'][8]<=1e-10,continuity=max(au['coupled_norms'][4:])<1e-9,fixed_k=d['anchor_k_sha256']==c['anchor_k_sha256'])
 gates['native_GL16']=divergence_gate(d,c['divergence_limit'])['passed']
 for key,v in seen.items():gates.update({key+'_margin':min(v['inner_minimum_covariance'],v['outer_minimum_covariance'])>=delta-1e-12,key+'_scalar':v['scalar_positive_finite'] and v['scalar_envelopes_passed'],key+'_velocity':v['velocity_change_L2']<=.12*(1+1e-10) and v['maximum_sampled_velocity_increment']<=.5+1e-10,key+'_divergence':v['strong_divergence_L2']<=c['divergence_limit']+1e-8})
 r=dict(version='H2iRepeatableSeenV1',passed=all(gates.values()),gates=gates,state_sha256=c['state_sha256'],base_sha256=COLD_HASH,anchor_k_sha256=c['anchor_k_sha256'],basis_sha256=c['basis_sha256'],context_sha256=sha(a.context),model_specification=load(Path(rr['output']+'.model.json')),output=rr['output'],output_sha256=sha(Path(rr['output'])),runtime=runtime_identity(a.probe),seconds=time.monotonic()-start,seen_orders=[20,24],fresh_audits=0);atomic(a.output/'receipt.json',r);entries[-1].update(receipt=str(a.output/'receipt.json'),passed=r['passed']);atomic(ledger,entries)

def fresh(a):
 c=load(a.candidate);seen=load(a.seen);require(seen['passed'] and all(seen['gates'].values()) and seen['state_sha256']==c['state_sha256'] and seen['anchor_k_sha256']==c['anchor_k_sha256'],'fresh GL28 requires frozen seen-qualified candidate');ctx=Context(a.context);require(c['context_sha256']==seen['context_sha256']==sha(a.context) and c['basis_sha256']==seen['basis_sha256']==ctx.doc['basis_sha256'] and c['base_sha256']==seen['base_sha256']==COLD_HASH,'fresh context/base/basis binding differs');require(same_runtime(seen['runtime'],runtime_identity(a.probe)),'fresh runtime differs from seen validation');require(sha(Path(seen['output']))==seen['output_sha256'],'seen native packet changed');ledger=a.output.parent/'fresh-GL28.json';require(not ledger.exists(),'fresh GL28 can run only once, never as seen candidate check');atomic(ledger,dict(candidate_sha256=c['state_sha256'],anchor_k_sha256=c['anchor_k_sha256'],seen_receipt_sha256=sha(a.seen),status='reserved_once',created=now()));a.output.mkdir(exist_ok=False);xml=case(a.output/'case',a.h1,load(a.preflight),flat=False);packet=dict(operation='h2i_fresh',inner_stress_model='Linear',expected_model='StandardKEpsilon',xml_sha256=sha(xml),state=c['state'],state_sha256=c['state_sha256'],base_sha256=COLD_HASH,coupled=True,anchor_k=c['anchor_k'],anchor_k_sha256=c['anchor_k_sha256'],fresh_order=28);_,rr=run_probe(a.probe,xml,packet,a.output/'native');d=load(Path(rr['output']));v=d['fresh_sampling'];K0=seen['model_specification']['K0'];gates=fresh_gates(v,c,seen['model_specification']);gates['native_binding']=d['state_sha256']==c['state_sha256'] and d['anchor_k_sha256']==c['anchor_k_sha256'] and all(d[k]==ctx.doc[k] for k in ('patch_signature','geometry_identity','inventory_identity','reference_sha256','defect_sha256'))
 r=dict(version='H2iFreshGL28V1',passed=all(gates.values()),gates=gates,state_sha256=c['state_sha256'],base_sha256=COLD_HASH,anchor_k_sha256=c['anchor_k_sha256'],basis_sha256=c['basis_sha256'],seen_receipt_file=str(a.seen),seen_receipt_sha256=sha(a.seen),output=rr['output'],output_sha256=sha(Path(rr['output'])),fresh_sampling=v,model_specification=seen['model_specification'],probe_runtime=runtime_identity(a.probe),solver_runtime=runtime_identity(a.solver));atomic(a.output/'receipt.json',r);atomic(ledger,dict(status='completed_once',passed=r['passed'],receipt=str(a.output/'receipt.json'),candidate_sha256=c['state_sha256']))
 if not r['passed']:return
 native=load(Path(seen['output']));src=Path(__file__).resolve().parents[4];names=['library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp','library/UnitTests/SolverUtils/RANSWallJointFeasibilityContext.h','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionJointInitializer.cpp'];compat=complete_receipt(Path(seen['output']),seen['model_specification'],c['divergence_limit'],runtime_identity(a.probe),runtime_identity(a.solver),v,dict(file=str(a.output/'receipt.json'),sha256=sha(a.output/'receipt.json')));validate_receipt(compat);combined=a.output/'issuance-receipt.json';atomic(combined,compat);prov=dict(K0=K0,divergence_limit=c['divergence_limit'],basis_sha256=c['basis_sha256'],context_sha256=c['context_sha256'],native_receipt_file=str(combined),source_revision='H2i exact hashes',source_hashes={n:sha(src/n) for n in names});cert=certificate_v2(native,compat,native['fixed_reference'],runtime_identity(a.solver),runtime_identity(a.probe),prov);cert['model_specification']=seen['model_specification'];cert['H2i_fresh']=dict(receipt=str(a.output/'receipt.json'),sha256=sha(a.output/'receipt.json'),anchor_k_sha256=c['anchor_k_sha256'],required_order=28);atomic(a.output/'certificate.json',cert)

def pilot(a):
 c=load(a.certificate);binding=c.get('H2i_fresh');require(binding and binding['required_order']==28,'H2i GL28 binding required');r=load(Path(binding['receipt']));require(sha(Path(binding['receipt']))==binding['sha256'] and r['passed'] and all(r['gates'].values()),'H2i fresh receipt failed/changed');require(r['state_sha256']==c['state_sha256'] and r['anchor_k_sha256']==binding['anchor_k_sha256'] and r['base_sha256']==c['immutable_base_sha256'] and r['basis_sha256']==c['basis_sha256'],'H2i candidate/anchor/base/basis identity changed');require(all(fresh_gates(r['fresh_sampling'],dict(state_sha256=c['state_sha256'],divergence_limit=.06030809942199109),c['model_specification']).values()),'H2i actual fresh sampling gate failed');require(same_runtime(r['solver_runtime'],runtime_identity(a.solver)),'H2i solver closure differs');StartupRun(a.solver,a.certificate,a.preflight,a.h1,a.output,builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=flat)).pair()

def flat(a):
 a.output.mkdir(exist_ok=False);out={}
 for model in ('cls','standard'):
  xml=case(a.output/model,a.h1,load(a.preflight),flat=True,model=model);_,r=run_probe(a.probe,xml,dict(operation='h2i_flat',inner_stress_model='Linear',xml_sha256=sha(xml),expected_model='StandardKEpsilon' if model=='standard' else 'CLSHighRe'),a.output/(model+'-probe'));out[model]=dict(receipt=r,data=load(Path(r['output']))['context'])
 atomic(a.output/'results.json',out)

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('prepare','search','repeatable','fresh','pilot','flat'))
 for n in ('parent','probe','solver','h1','preflight','context','candidate','seen','certificate'):p.add_argument('--'+n,type=Path)
 p.add_argument('--output',type=Path,required=True);a=p.parse_args();os.environ.update(THREADS);os.environ['OMP_THREAD_LIMIT']='1';globals()[a.phase](a)
if __name__=='__main__':main()
