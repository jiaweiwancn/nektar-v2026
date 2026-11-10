#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only scalar extrema on fixed exported native quadrature/face rows.

Mandatory native inventory minima remain authoritative for admissibility.
These extra maxima are finite-grid observations, not global positivity proofs.
One-sweep intermediate scalar fields equal the recorded old/new scalar blocks.
"""
import argparse,gzip,json,os,sys
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from scipy.sparse import csr_matrix
from qualify_rans_preservation import probe,replay_request,load
from rans_qualification import THREADS,write_json,exact_vector_hash
from rans_request_manifest import RequestManifest
from rans_preservation_oracle import Reference,weak_balance
from rans_pilot_evidence import epsilon_budget_gate

def matrix(rows,n):
 data=[];ids=[];offset=[0]
 for r in rows:data.extend(r['phi']);ids.extend(r['ids']);offset.append(len(ids))
 return csr_matrix((data,ids,offset),shape=(len(rows),n))

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--through-time',type=float);p.add_argument('--output-label',default='scalar-observations');a=p.parse_args();r=a.root.resolve();folder=r/'campaigns'/a.label;req=load(folder/'request-manifest.json')['request'];admission=RequestManifest(folder,req,True)
 os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
 out=folder/a.output_label;out.mkdir(exist_ok=False);seed=load(folder/'initial-packet.json')
 # Use an unchanged session; no solver action or checkpoint mutation.
 import shutil
 shutil.copy2(folder/('T-'+format(req['horizons'][0],'.4f'))/'case.xml',out/'case.xml')
 packet=probe(r/'runtime/RANSPreservationProbe',out,replay_request(seed,[],export_basis=True,state=seed['reference']))['observation'];off=packet['offsets'];n=packet['sizes'];points=packet['volume']+packet['faces'];groups={}
 for region,name in ((0,'inner_k'),(1,'outer_k')):
  q=[x for x in points if x['region']==region];groups[name]=(matrix([x['velocity_k'] for x in q],n[2]),q,2)
 q=[x for x in points if x['region']==1];groups['epsilon_outer']=(matrix([x['epsilon'] for x in q],n[3]),q,3)
 q=[x for x in packet['faces'] if x['interface'] and x['region']==0];groups['sampled_K']=(matrix([x['velocity_k'] for x in q],n[2]),q,2)
 records=[];minima={};maxima={};failure=[];epsilon_weak_budgets=[]
 em=groups['epsilon_outer'][0].toarray();constant,_,rank,_=np.linalg.lstsq(em,np.ones(em.shape[0]),rcond=None);assert rank==n[3] and np.max(np.abs(em@constant-1))<1e-12
 assert np.max(np.abs(constant-np.rint(constant)))<1e-10;constant=np.rint(constant)
 def observe(state,t,stage):
  row=dict(time=t,stage=stage,fields={},state_sha256=exact_vector_hash(state))
  for name,(m,q,b) in groups.items():
   values=m@np.array(state[off[b]:off[b]+n[b]]);imin=int(np.argmin(values));imax=int(np.argmax(values));v={}
   for kind,i,allwin in [('minimum',imin,minima),('maximum',imax,maxima)]:
    item=dict(value=float(values[i]),time=t,stage=stage,coordinate=q[i]['coordinate'],region=q[i]['region'],geometry=q[i]['velocity_k' if b==2 else 'epsilon']['geometry'],side=(-1 if q[i].get('interface') and q[i]['region']==0 else 1 if q[i].get('interface') else 0));v[kind]=item
    if name not in allwin or ((item['value']<allwin[name]['value']) if kind=='minimum' else (item['value']>allwin[name]['value'])):allwin[name]=item
   row['fields'][name]=v
   if not np.all(np.isfinite(values)) or np.min(values)<0:failure.append(dict(time=t,stage=stage,field=name,minimum=float(np.min(values))))
  records.append(row)
 reports=[];report_sources=[]
 for f in sorted(folder.glob('T-*/case.evolution.json')):
  if not admission.completed(f.parent):continue
  report=load(f)
  if a.through_time is not None and report['time']>a.through_time+1e-13:continue
  reports.append(report);report_sources.append(dict(path=str(f),sha256=__import__('hashlib').sha256(f.read_bytes()).hexdigest()));states=report['accepted_states'];times=report['audit_times'];assert len(states)==len(times)
  for k,(t,state) in enumerate(zip(times,states,strict=True)):
   observe(state,t,'accepted')
   if k:
    previous=states[k-1];observe(previous,t,'input/flow scalars')
    hybrid=list(previous);hybrid[off[3]:off[3]+n[3]]=state[off[3]:off[3]+n[3]];observe(hybrid,t,'epsilon scalars');observe(state,t,'k scalars')
  if report.get('final_stage_loads'):
   raw=np.array(report['final_stage_loads'][off[3]:off[3]+n[3]]);fixed={int(v[0])-off[3] for v in report['final_stage_constraints'] if off[3]<=int(v[0])<off[3]+n[3]}
   total=float(constant@raw);reaction=sum(constant[i]*raw[i] for i in fixed);epsilon_weak_budgets.append(dict(time=report['time'],alpha=report['last_stage_alpha'],total_load=total,boundary_reaction=float(reaction),free_constant_test_defect=float(total-reaction),meaning='actual native time-discrete load and prescribed-row reaction; includes declared splitting defect'))
  # Check reconstructed input identity against the actual stage receipt.
  for old,attempt in zip(states[:-1],report['pimple_attempts'][:len(states)-1],strict=True):
   assert exact_vector_hash(old)==attempt['incoming_sha256'],'one-sweep stage reconstruction identity mismatch'
  for event in report.get('physical_range_events',[]):
   if event.get('state'):observe(event['state'],event['time'],'rejected_candidate')
 weak=weak_balance(Reference(req['model'],req['policy'],req['wall']),packet)
 budget=epsilon_budget_gate(epsilon_weak_budgets,Reference(req['model'],req['policy'],req['wall']))
 result=dict(status='pass' if not failure and weak['status']=='pass' and budget['status']=='pass' else 'fail',epsilon_budget_gate=budget,analytical_weak_balance=weak,scope='fixed reference-quadrature and native face rows; extrema over all retained accepted states and reconstructed one-sweep scalar stages',minimum=minima,maximum=maxima,observations=records,failures=failure,
  epsilon_weak_budgets=epsilon_weak_budgets,source_packet_sha256=__import__('hashlib').sha256((out/'packet.json').read_bytes()).hexdigest(),roundoff_floor=1e-12,mandatory_native_domain_audits_unchanged=True)
 from rans_duration_contract import scalar_coverage
 result['configuration']={k:req[k] for k in ('model','policy','wall','descriptor')};result['trajectory_request_sha256']=admission.payload['sha256'];result['sealed_native_reports']=report_sources
 result['coverage']=scalar_coverage(result,reports,result['configuration'],a.through_time if a.through_time is not None else max(x['time'] for x in reports))
 result['epsilon_constant_test_interpretation']='prescribed-row reaction identity; all nonzero constant coefficients constrained on this one-column patch; not an independent physical conservation theorem'
 write_json(out/'results.json',result)
 with (out/'packet.json').open('rb') as src,gzip.open(out/'native-basis.json.gz','wb',compresslevel=3) as dst:__import__('shutil').copyfileobj(src,dst)
 print(result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
