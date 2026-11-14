#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exactly one fresh dense .0005 CT9 history and real-identity comparison."""
import argparse,copy,json,os
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted8_reference import POLICY
from ct_accepted9_gate import pair
from ct_accepted_feasibility import write_json,require
from ct_accepted_pilot import replay,constitutive_comparison,artifact_hash
from ct_accepted3_metrics import domain
from ct_accepted_open import run,probe
from ct_accepted5_traction import profiles,difference
from rans_request_manifest import RequestManifest,runtime_identity,sha

def read(p):return json.loads(Path(p).read_text())
def sealed(d):
    q=read(d/'completed-artifacts.json');require(all(sha(d/k)==v for k,v in q['files'].items()),'changed sealed artifact:'+str(d));return read(d/'results.json')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('ct8','runtime','contract','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    pre=sealed(a.ct8/'preflight-v6');old=sealed(a.ct8/'strip-T008-v1');sealed(a.ct8/'controls-v2')
    require(old['qualification']['reasons']==['temporal:Pi'],'historical CT8 classification unexpectedly changed')
    expect=read(a.ct8/'controls-v2/request-manifest.json')['request']['runtime'];require(runtime_identity(a.runtime/'IncNavierStokesSolver')==expect,'changed governing dense runtime')
    c=pre['case'];m=RequestManifest(a.output,dict(schema='CT9PressureFinerHistoryV1',dt=.0005,steps=16,target=.008,case=c,contract=sha(a.contract),ct8_pair=sha(a.ct8/'strip-T008-v1/results.json'),preflight=sha(a.ct8/'preflight-v6/results.json'),runtime=expect,source=sha(Path(__file__)),comparison='new .001/.0005 pair; no old labels or histories edited'))
    ini=read(Path(pre['initial']['0.001']['directory'])/'initial.json');old_native=read(Path(old['members']['0.001']['directory'])/'case.evolution.json')
    require(sha(Path(old['members']['0.001']['directory'])/'case.evolution.json')==old['members']['0.001']['native_report_sha256'],'changed reused report')
    zero=a.output/'zero';n.make(zero,a.runtime,c,.0005,0);e=run([a.runtime/'IncNavierStokesSolver','case.xml'],zero,'solver.log');z=read(zero/'case.evolution.json')
    require(e['exit_code']==0 and not z['failure'] and z['time']==0 and z['state']==ini['observation']['state']==old_native['accepted_states'][0],'fresh finer actual initialization mismatch')
    require(exact_vector_hash(z['state'])==exact_vector_hash(ini['observation']['state'])==exact_vector_hash(old_native['accepted_states'][0]),'fresh finer initialized bytes differ')
    for key in ('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256','open_traction_load_sha256','open_diffusive_load_sha256','open_advective_load_sha256'):
        require(z['initialization_receipt'][key]==pre['initial']['0.001']['initialization_receipt'][key],'fresh finer effective boundary/reference differs:'+key)
    folder=a.output/'trajectory';n.make(folder,a.runtime,c,.0005,16,Path(z['checkpoint']),zero/'case.xml')
    os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE')
    r=read(folder/'case.evolution.json');require(r['accepted_states'][0]==z['state'],'changed consumed history')
    # Cache only immutable continuous comparison values. No governing evaluation
    # or state-dependent validity is cached or altered.
    n.reference=lru_cache(None)(n.reference)
    row=dict(dt=r['dt'],time=r['time'],directory=str(folder),failure=r['failure'],execution=e,signature=r['signature'],state_sha256=exact_vector_hash(r['state']),history_sha256=artifact_hash(r['checkpoint']),native_report_sha256=sha(folder/'case.evolution.json'),initialization_receipt=r['initialization_receipt'],actual_initial_sha256=exact_vector_hash(z['state']),zero_report_sha256=sha(zero/'case.evolution.json'),zero_history_sha256=artifact_hash(z['checkpoint']))
    row['accepted_metrics']=[dict(time=t,**n.physical_metrics(ini,s,c,pre['weights'],pre['global_count'])) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)]
    row['traction_history']=[dict(time=t,profile=profiles(ini,s,c)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)]
    row['actual_initial_traction_profile']=profiles(ini,z['state'],c);row['first_accuracy_crossing']=next((v for v in row['accepted_metrics'] if v['status']!='pass'),None)
    row['domain']=domain(ini,r,c);row['budgets']=n.budgets(a.runtime,folder,ini,r,c) if r['steps'] else dict(gate=dict(status='missing'))
    grid,sides,_,_=n.inventory(c);req=replay(ini,r['state'],r['time'],grid,export_constitutive=True);req['grid_sides']=sides
    q,pe=probe(a.runtime,folder,req,'endpoint');row['constitutive']=constitutive_comparison(q,c) if pe['exit_code']==0 else dict(status='fail',reason=q.get('failure'))
    row['corrections']=sum(v.get('bound_calls',0) for v in r['pimple_attempts']);row['extra_solves']=sum(v.get('additional_governing_equation_solves',0) for v in r['pimple_attempts']);row['native_audited_extrema']=audited_extrema(r)
    oldrow=copy.deepcopy(old['members']['0.001']);oldrow.update(dt=old_native['dt'],actual_initial_sha256=exact_vector_hash(old_native['accepted_states'][0]),reuse_receipt_sha256=sha(a.ct8/'strip-T008-v1/results.json'))
    result=dict(schema='CT9PhysicalPairV1',required_dt=[.001,.0005],target=.008,case=c,members={'0.001':oldrow,'0.0005':row},same_actual_initial=r['accepted_states'][0]==old_native['accepted_states'][0],temporal={},traction_temporal={},common_time_fields=[],historical_CT8='failed temporal:Pi unchanged')
    if abs(r['time']-.008)<1e-13 and abs(old_native['time']-.008)<1e-13:
        end=[n.observations(ini,v,pre['global_count']) for v in (old_native['state'],r['state'])];delta=n.norms(*end,pre['weights'],ini['observation']['grid_regions'][:pre['global_count']],c)
        for f,v in delta.items():
            drift=max(x['accepted_metrics'][-1]['drift'][f]['scaled'] for x in result['members'].values());result['temporal'][f]=dict(difference=v,measured_drift=drift,limit=max(1e-7,.1*drift))
        rows=list(result['members'].values());spread=difference(*(x['traction_history'][-1]['profile'] for x in rows),c);drifts=[difference(x['actual_initial_traction_profile'],x['traction_history'][-1]['profile'],c) for x in rows]
        for f in ('ground','interface','signed_mean'):
            d=abs(spread['signed_mean_difference']) if f=='signed_mean' else spread['profiles'][f]['normalized_weighted_RMS'];v=max(abs(x['signed_mean_difference']) if f=='signed_mean' else x['profiles'][f]['normalized_weighted_RMS'] for x in drifts);result['traction_temporal'][f]=dict(spread=d,measured_drift=v,limit=max(1e-7,.1*v))
    left={round(t,12):s for t,s in zip(old_native['audit_times'],old_native['accepted_states'],strict=True)};right={round(t,12):s for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)}
    for t in sorted(left.keys()&right.keys()):
        v=[n.observations(ini,s,pre['global_count']) for s in (left[t],right[t])];result['common_time_fields'].append(dict(time=t,difference=n.norms(*v,pre['weights'],ini['observation']['grid_regions'][:pre['global_count']],c),state_sha256=[exact_vector_hash(s) for s in (left[t],right[t])]))
    result['qualification']=pair(result);write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT9 fresh dense pair',result['qualification'],flush=True)
    return 0 if result['qualification']['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
