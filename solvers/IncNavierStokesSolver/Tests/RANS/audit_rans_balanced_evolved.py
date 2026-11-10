#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native model-aware term replay of sealed evolved states, without advancement."""
import argparse,gzip,hashlib,json,os,shutil
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import THREADS,write_json,exact_vector_hash
from rans_preservation_oracle import Reference,weak_balance
from qualify_rans_preservation import load,replay_request
from rans_request_manifest import RequestManifest,runtime_identity
from rans_spatial_resolution import analyze,common_test,test_parent_geometries,test_effective_coordinates
from rans_pilot_evidence import epsilon_budget_gate
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from diagnose_rans_balanced_spatial import run_probe

def epsilon_budget(packet,native,ref):
 o=packet['observation'];n=o['sizes'][3];off=o['offsets'][3];b=np.zeros(n)
 for q in o['volume']:
  if q['region']==1:b[np.array(q['epsilon']['ids'],int)]+=q['weight']*np.array(q['epsilon']['phi'])
 c=cho_solve(cho_factor(np.array(packet['resolution']['mass'][3]),lower=True),b)
 assert np.max(np.abs(c-np.rint(c)))<1e-9;c=np.rint(c)
 load=np.array(native['final_stage_loads'])[off:off+n];fixed={int(i)-off for i,v in native['final_stage_constraints'] if off<=int(i)<off+n}
 total=float(c@load);reaction=float(sum(c[i]*load[i] for i in fixed));free_nonzero=[i for i in range(n) if i not in fixed and c[i]!=0]
 gate=epsilon_budget_gate([dict(time=native['time'],alpha=native['last_stage_alpha'],free_constant_test_defect=total-reaction)],ref)
 gate.update(total_load=total,prescribed_row_reaction=reaction,nonzero_free_constant_coefficients=len(free_nonzero),
  interpretation='reaction identity only: constant coefficients are all prescribed' if not free_nonzero else 'free constant-test load/reaction balance',
  dependent_epsilon_boundary_rate='not inferred; this is the actual frozen-trace time-discrete load')
 return gate

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--analysis-runtime',type=Path);p.add_argument('--campaign-root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--output-label',required=True);p.add_argument('--times',type=float,nargs='+',default=[.0016,.016]);p.add_argument('--common-directory',type=Path);p.add_argument('--constitutive-samples',action='store_true');a=p.parse_args();r=a.root.resolve();runtime=a.analysis_runtime.resolve() if a.analysis_runtime else r/'runtime';folder=a.campaign_root/'campaigns'/a.label;manifest=load(folder/'request-manifest.json');req=manifest['request'];admission=RequestManifest(folder,req,True);ref=Reference(req['model'],req['policy'],req['wall']);out=r/'evolved-audits'/a.output_label;out.mkdir(parents=True,exist_ok=False)
 os.environ.update(THREADS,LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'))
 source=load(folder/'initial-packet.json');common=a.common_directory or r/'diagnosis'/('initial-'+req['wall'])/'common'
 with gzip.open(common/'packet.json.gz','rt') as f:cp=json.load(f)
 write_json(out/'immutable-request.json',dict(schema='BalancedEvolvedAuditV1',trajectory_request_sha256=manifest['sha256'],times=a.times,alpha=0,common_packet_sha256=hashlib.sha256((common/'packet.json.gz').read_bytes()).hexdigest(),analysis_runtime=runtime_identity(runtime/'RANSPreservationProbe')))
 rows=[]
 for t in a.times:
  endpoint=folder/('T-'+format(t,'.4f'));assert admission.completed(endpoint);native=load(endpoint/'case.evolution.json');assert abs(native['time']-t)<1e-13 and not native['failure'];destination=out/('T-'+format(t,'.4f'));destination.mkdir();shutil.copy2(endpoint/'case.xml',destination/'case.xml')
  request=replay_request(source,[],export_basis=True,state=native['state'],time=native['time']);request['fixed_boundary']=native['initialization_receipt']['live_fixed_boundary'];request['export_resolution']=True
  if a.constitutive_samples:request['grid']=req['evaluation']['grid'];request['export_constitutive']=True
  packet,seconds=run_probe(runtime,destination,request);assert exact_vector_hash(packet['observation']['state'])==exact_vector_hash(native['state'])
  cross=out/('T-'+format(t,'.4f')+'-common-cross');cross.mkdir();shutil.copy2(common/'case.xml',cross/'case.xml');request=replay_request(cp,test_effective_coordinates(packet));request['export_test_basis']=True;request['grid_test_geometries']=test_parent_geometries(cp,packet);xp,_=run_probe(runtime,cross,request)
  row=dict(status='pass',time=native['time'],nominal_target=t,trajectory_source=str(endpoint),trajectory_request_sha256=manifest['sha256'],native_report_sha256=hashlib.sha256((endpoint/'case.evolution.json').read_bytes()).hexdigest(),history=load(endpoint/'metrics.json')['history_identity'],
   state_sha256=exact_vector_hash(native['state']),live_boundary_sha256=native['initialization_receipt']['explicit_boundary_sha256'],
   spatial=analyze(packet,ref),common_test=common_test(cp,packet,xp,ref),analytical_weak_balance=weak_balance(ref,packet['observation']),
   native_alpha_zero=packet['observation']['audit_alpha_zero'],time_discrete_stage_residual=native['stage_audits'][-1:],splitting_maxima=load(endpoint/'metrics.json')['source_and_constraint_maxima'],read_only_export_seconds=seconds)
  row['epsilon_constant_test_budget']=epsilon_budget(packet,native,ref)
  row['configuration']={k:req[k] for k in ('model','policy','wall','descriptor')}
  row['status']='pass' if row['analytical_weak_balance']['status']=='pass' and row['epsilon_constant_test_budget']['status']=='pass' else 'fail'
  write_json(destination/'analysis.json',row);rows.append(row);write_json(out/'results.json',dict(status='running',checks=rows));print(t,'pass',flush=True)
 status='pass' if all(x['status']=='pass' for x in rows) else 'fail'
 write_json(out/'results.json',dict(status=status,checks=rows,scope='sealed actual evolved fields; native alpha-zero term replay, fixed common tests, separate native stage residuals'));return 0 if status=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
