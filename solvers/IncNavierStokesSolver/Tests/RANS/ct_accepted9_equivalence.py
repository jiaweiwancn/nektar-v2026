#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh sparse CT8 physical pair as a backend equivalence experiment only.

Its physical uncertainty outcome must reproduce the sealed dense failure.
Passing this job qualifies algebraic equivalence, never physical continuation.
"""
import argparse,copy,os,resource
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_controls import enable
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair,backend_controls,backend_trajectory
from ct_accepted_pilot import replay,constitutive_comparison,artifact_hash
from ct_accepted3_metrics import domain
from ct_accepted_open import run,probe,require
from ct_accepted5_traction import profiles,difference
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','ct8','controls','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    controls=sealed(a.controls);require(backend_controls(controls)['status']=='pass','small backend equivalence primitives not admitted');require(n.read(a.controls/'request-manifest.json')['request']['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver'),'foreign actual backend runtime');pre=sealed(a.ct8/'preflight-v6');old=sealed(a.ct8/'strip-T008-v1');c=pre['case'];count=pre['global_count'];weights=pre['weights'];n.reference=lru_cache(None)(n.reference)
    require(old['qualification']['reasons']==['temporal:Pi'],'historical classification changed')
    m=RequestManifest(a.output,dict(schema='CT9BackendTrajectoryEquivalenceV1',ct8_pair=sha(a.ct8/'strip-T008-v1/results.json'),controls=sha(a.controls/'results.json'),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source=sha(Path(__file__)),dt=[.002,.001],target=.008,backend='AcceptedSerialSparseV1',maximum_scaled_equivalence=1e-8,scope='new histories; old dense physical failure expected and independently recomputed'))
    result=dict(schema='CT9PhysicalPairV1',required_dt=[.002,.001],target=.008,case=c,members={},equivalence={},temporal={},traction_temporal={});packets=[]
    for dt in (.002,.001):
        name=str(dt);dense=n.read(Path(old['members'][name]['directory'])/'case.evolution.json');require(sha(Path(old['members'][name]['directory'])/'case.evolution.json')==old['members'][name]['native_report_sha256'],'changed dense native report')
        zero=a.output/('zero-'+name);n.make(zero,a.runtime,c,dt,0);enable(zero);e=run([a.runtime/'IncNavierStokesSolver','case.xml'],zero,'solver.log');z=n.read(zero/'case.evolution.json');require(e['exit_code']==0 and not z['failure'],'sparse zero initialization failed')
        ini,pe=probe(a.runtime,zero,n.request(c),'initial');require(pe['exit_code']==0,'sparse initial probe failed');n.binding(ini,c)
        require(z['state']==ini['observation']['state']==dense['accepted_states'][0],'different actual initialized coefficients across backends')
        dense_ini=n.read(Path(pre['initial'][name]['directory'])/'initial.json')
        for key in ('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256','open_traction_load_sha256','open_diffusive_load_sha256','open_advective_load_sha256'):
            require(z['initialization_receipt'][key]==pre['initial'][name]['initialization_receipt'][key],'effective physical boundary/reference mismatch:'+key)
        folder=a.output/name;n.make(folder,a.runtime,c,dt,round(.008/dt),Path(z['checkpoint']),zero/'case.xml');os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE');r=n.read(folder/'case.evolution.json')
        require(r['accepted_states'][0]==z['state'],'sparse consumed history mismatch')
        require('AcceptedSerialSparseV1' in r['signature'] and all(x.get('linear_backend')=='AcceptedSerialSparseV1' for x in r['pimple_attempts']),'native trajectory did not consume the requested backend')
        row=dict(dt=r['dt'],time=r['time'],directory=str(folder),failure=r['failure'],execution=e,signature=r['signature'],initialization_receipt=r['initialization_receipt'],actual_initial_sha256=exact_vector_hash(z['state']),native_report_sha256=sha(folder/'case.evolution.json'),history_sha256=artifact_hash(r['checkpoint']))
        row['accepted_metrics']=[dict(time=t,**n.physical_metrics(ini,s,c,weights,count)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['traction_history']=[dict(time=t,profile=profiles(ini,s,c)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['actual_initial_traction_profile']=profiles(ini,z['state'],c)
        row['domain']=domain(ini,r,c);row['budgets']=n.budgets(a.runtime,folder,ini,r,c) if r['steps'] else dict(gate=dict(status='missing'))
        grid,sides,_,_=n.inventory(c);req=replay(ini,r['state'],r['time'],grid,export_constitutive=True);req['grid_sides']=sides;q,pe=probe(a.runtime,folder,req,'endpoint');row['constitutive']=constitutive_comparison(q,c) if pe['exit_code']==0 else dict(status='fail',reason=q.get('failure'))
        row['corrections']=sum(x.get('bound_calls',0) for x in r['pimple_attempts']);row['extra_solves']=sum(x.get('additional_governing_equation_solves',0) for x in r['pimple_attempts']);row['native_audited_extrema']=audited_extrema(r)
        result['members'][name]=row;packets.append((ini,r));common=[]
        left={round(t,12):s for t,s in zip(dense['audit_times'],dense['accepted_states'],strict=True)};right={round(t,12):s for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)}
        require(left.keys()==right.keys(),'unequal dense/sparse accepted time coverage')
        for t in left:
            v=[n.observations(p,s,count) for p,s in ((dense_ini,left[t]),(ini,right[t]))];d=n.norms(*v,weights,ini['observation']['grid_regions'][:count],c);tr=difference(profiles(dense_ini,left[t],c),profiles(ini,right[t],c),c)
            common.append(dict(time=t,fields=d,traction=tr))
        largest=max((d['scaled'] for q in common for d in q['fields'].values()),default=float('inf'));traction_largest=max((q['traction']['profiles'][f]['normalized_weighted_RMS'] for q in common for f in ('ground','interface')),default=float('inf'))
        result['equivalence'][name]=dict(status='pass' if max(largest,traction_largest)<=1e-8 else 'fail',maximum_scaled_field_difference=largest,maximum_scaled_traction_difference=traction_largest,common_times=common,old_runtime=n.read(a.ct8/'controls-v2/request-manifest.json')['request']['runtime'],dense_report=old['members'][name]['native_report_sha256'],sparse_report=row['native_report_sha256'])
        write_json(a.output/'progress.json',result);print('CT9 trajectory backend comparison',name,result['equivalence'][name]['status'],flush=True)
    result['same_actual_initial']=packets[0][1]['accepted_states'][0]==packets[1][1]['accepted_states'][0]
    require(all(abs(r['time']-.008)<1e-13 for _,r in packets),'sparse pair incomplete')
    delta=n.norms(*(n.observations(p,r['state'],count) for p,r in packets),weights,packets[0][0]['observation']['grid_regions'][:count],c)
    for f,v in delta.items():
        drift=max(x['accepted_metrics'][-1]['drift'][f]['scaled'] for x in result['members'].values());result['temporal'][f]=dict(difference=v,measured_drift=drift,limit=max(1e-7,.1*drift))
    rows=list(result['members'].values());tr=difference(*(x['traction_history'][-1]['profile'] for x in rows),c);drifts=[difference(x['actual_initial_traction_profile'],x['traction_history'][-1]['profile'],c) for x in rows]
    for f in ('ground','interface','signed_mean'):
        d=abs(tr['signed_mean_difference']) if f=='signed_mean' else tr['profiles'][f]['normalized_weighted_RMS'];v=max(abs(x['signed_mean_difference']) if f=='signed_mean' else x['profiles'][f]['normalized_weighted_RMS'] for x in drifts);result['traction_temporal'][f]=dict(spread=d,measured_drift=v,limit=max(1e-7,.1*v))
    result['physical_qualification']=pair(result);result['backend_qualification']=backend_trajectory(result);result['backend_status']=result['backend_qualification']['status'];result['physical_continuation_admitted']=False
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT9 backend',result['backend_status'],'physical',result['physical_qualification'],flush=True);return 0 if result['backend_status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
