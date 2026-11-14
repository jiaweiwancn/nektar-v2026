#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT10 one frozen radial contrast, fresh sparse histories, primitive decisions."""
import argparse,copy,os,resource,shutil
from functools import lru_cache
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash,audited_extrema
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_controls import enable
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair
from ct_accepted3_metrics import domain
from ct_accepted_pilot import replay,constitutive_comparison,artifact_hash
from ct_accepted_open import run,probe,require
from ct_accepted5_traction import profiles,difference
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

RUNTIME_SHA='5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7'
BOUNDARY_KEYS=('raw_initial_sha256','effective_initial_sha256','fixed_reference_sha256','explicit_boundary_sha256','implicit_boundary_sha256','open_traction_load_sha256','open_diffusive_load_sha256','open_advective_load_sha256')

def common_request(c,initial=None,state=None,time=0.):
    req=n.request(c)
    if initial is not None:
        grid,sides,_,_=n.inventory(c);req=replay(initial,state,time,grid,grid_sides=sides,export_basis=True,export_test_basis=True,export_constitutive=True,export_neutral_composite_reference=True)
    return req

def preflight(a):
    old=sealed(a.ct8/'preflight-v6');diag=sealed(a.diagnosis);selection=diag['selection'];base=copy.deepcopy(old['case']);fine=copy.deepcopy(base);fine['knots']=selection['knots'];common=copy.deepcopy(fine)
    require(selection['status']=='selected_before_evolution' and set(base['knots'])<=set(fine['knots']) and fine['knots']!=base['knots'],'invalid single selected contrast')
    m=RequestManifest(a.output,dict(schema='CT10StaticContrastV1',contract=sha(a.contract),diagnosis=sha(a.diagnosis/'results.json'),selection=sha(a.diagnosis/'selection.json'),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__)),cases=dict(baseline=base,candidate=fine),timesteps=[.001,.0005]))
    grid,sides,weights,count=n.inventory(common);out=dict(schema='CT10StaticContrastV1',cases=dict(baseline=base,candidate=fine),common=common,weights=weights,global_count=count,initial={},selection=selection,status='incomplete')
    # Baseline coefficients/history stay sealed; only export common physical rows.
    folder=a.output/'baseline-common';folder.mkdir();oldzero=a.ct9/'equivalence-v2/zero-0.001';original=n.read(oldzero/'initial.json');shutil.copy2(oldzero/'case.xml',folder/'case.xml')
    q,e=probe(a.runtime,folder,common_request(common,original,original['observation']['state']),'initial');require(e['exit_code']==0,'baseline common-row export failed');n.binding(q,base)
    require(q['observation']['state']==original['observation']['state'],'baseline replay changed coefficients')
    out['initial']['baseline']=dict(directory=str(folder),packet_sha256=sha(folder/'initial.json'),actual_initial_sha256=exact_vector_hash(q['observation']['state']),projection_metrics=n.physical_metrics(q,q['observation']['state'],base,weights,count),source_packet_sha256=sha(oldzero/'initial.json'))
    hashes=[]
    for dt in (.001,.0005):
        folder=a.output/('candidate-zero-'+str(dt));n.make(folder,a.runtime,fine,dt,0);enable(folder)
        inventory,ce=probe(a.runtime,folder,dict(operation='native_layout_inventory'),'count');require(ce['exit_code']==0 and inventory['stored_slots']<=70000,'selected actual native count unsupported')
        e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');native=n.read(folder/'case.evolution.json');require(e['exit_code']==0 and not native['failure'] and native['time']==0.,'candidate static domain/startup failed')
        q,pe=probe(a.runtime,folder,common_request(common),'initial');require(pe['exit_code']==0,'candidate initial export failed');n.binding(q,fine);require(q['observation']['state']==native['state'],'candidate solver/probe initialization mismatch')
        mat=n.physical_metrics(q,native['state'],fine,weights,count);dom=domain(q,native,fine);ref=n.oracle_agreement(q,fine);h=exact_vector_hash(native['state']);hashes.append(h)
        out['initial']['candidate-'+str(dt)]=dict(directory=str(folder),native_count=inventory,execution=e,probe_execution=pe,actual_initial_sha256=h,history_sha256=artifact_hash(native['checkpoint']),native_report_sha256=sha(folder/'case.evolution.json'),packet_sha256=sha(folder/'initial.json'),initialization_receipt=native['initialization_receipt'],projection_metrics=mat,domain=dom,reference_agreement=ref,raw_spatial_norms=native['raw_reference_norms'],status='pass' if mat['status']==dom['status']==ref['status']=='pass' else 'fail')
        write_json(a.output/'progress.json',out)
    out['same_candidate_actual_initial']=len(set(hashes))==1
    out['status']='pass' if out['same_candidate_actual_initial'] and all(v.get('status','pass')=='pass' for v in out['initial'].values()) else 'fail'
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT10 static',out['status'],'candidate slots',inventory['stored_slots'],flush=True);return out


def collect(rt,folder,ini,r,c,weights,count,zero):
    row=dict(dt=r['dt'],time=r['time'],directory=str(folder),failure=r['failure'],signature=r['signature'],initialization_receipt=r['initialization_receipt'],actual_initial_sha256=exact_vector_hash(zero['state']),native_report_sha256=sha(folder/'case.evolution.json'),history_sha256=artifact_hash(r['checkpoint']))
    row['accepted_metrics']=[dict(time=t,**n.physical_metrics(ini,s,c,weights,count)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)]
    row['traction_history']=[dict(time=t,profile=profiles(ini,s,c)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)];row['actual_initial_traction_profile']=profiles(ini,zero['state'],c)
    row['first_accuracy_crossing']=next((x for x in row['accepted_metrics'] if x['status']!='pass'),None);row['domain']=domain(ini,r,c);row['budgets']=n.budgets(rt,folder,ini,r,c) if r['steps'] else dict(gate=dict(status='missing'))
    grid,sides,_,_=n.inventory(c);q,e=probe(rt,folder,replay(ini,r['state'],r['time'],grid,grid_sides=sides,export_constitutive=True),'endpoint');row['constitutive']=constitutive_comparison(q,c) if e['exit_code']==0 else dict(status='fail',reason=q.get('failure'))
    row['corrections']=sum(v.get('bound_calls',0) for v in r['pimple_attempts']);row['extra_solves']=sum(v.get('additional_governing_equation_solves',0) for v in r['pimple_attempts']);row['native_audited_extrema']=audited_extrema(r)
    return row


def comparison(result,initials,natives,c,weights,count):
    target=result['target'];keys=list(result['members']);result['same_actual_initial']=all(x['actual_initial_sha256']==result['members'][keys[0]]['actual_initial_sha256'] for x in result['members'].values());result.update(temporal={},traction_temporal={},common_times=[])
    if all(abs(natives[k]['time']-target)<1e-13 for k in keys):
        end=[n.observations(initials[k],natives[k]['state'],count) for k in keys];d=n.norms(*end,weights,initials[keys[0]]['observation']['grid_regions'][:count],c)
        for f,v in d.items():
            drift=max(r['accepted_metrics'][-1]['drift'][f]['scaled'] for r in result['members'].values());result['temporal'][f]=dict(difference=v,measured_drift=drift,limit=max(1e-7,.1*drift),fraction_of_absolute_tolerance=v['scaled']/(1e-4 if f=='w' else 1e-3))
        rows=list(result['members'].values());spread=difference(*(r['traction_history'][-1]['profile'] for r in rows),c);drifts=[difference(r['actual_initial_traction_profile'],r['traction_history'][-1]['profile'],c) for r in rows]
        for f in ('ground','interface','signed_mean'):
            x=abs(spread['signed_mean_difference']) if f=='signed_mean' else spread['profiles'][f]['normalized_weighted_RMS'];d=max(abs(r['signed_mean_difference']) if f=='signed_mean' else r['profiles'][f]['normalized_weighted_RMS'] for r in drifts);result['traction_temporal'][f]=dict(spread=x,measured_drift=d,limit=max(1e-7,.1*d),fraction_of_absolute_tolerance=x/1e-3)
    left={round(t,12):s for t,s in zip(natives[keys[0]]['audit_times'],natives[keys[0]]['accepted_states'],strict=True)};right={round(t,12):s for t,s in zip(natives[keys[1]]['audit_times'],natives[keys[1]]['accepted_states'],strict=True)}
    for t in sorted(left.keys()&right.keys()):
        d=n.norms(*(n.observations(initials[k],s,count) for k,s in zip(keys,(left[t],right[t]),strict=True)),weights,initials[keys[0]]['observation']['grid_regions'][:count],c);result['common_times'].append(dict(time=t,fields=d))
    result['qualification']=pair(result);return result


def advance(a):
    pre=sealed(a.preflight);require(pre['status']=='pass','single contrast static gate failed');c=pre['cases'][a.case];weights=pre['weights'];count=pre['global_count'];base=sealed(a.ct9/'equivalence-v2')
    parent=sealed(a.parent) if a.parent else None
    if a.target==.02:require(a.case=='candidate' and parent and parent['qualification']['status']=='pass' and parent['target']==.008,'candidate continuation not admitted')
    else:require(not parent and a.target==.008,'invalid primary history request')
    request=dict(schema='CT10CohortV1',case=c,kind=a.case,required_dt=[.001,.0005],target=a.target,preflight=sha(a.preflight/'results.json'),parent=sha(a.parent/'results.json') if parent else None,runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),source=sha(Path(__file__)),baseline_reuse=sha(a.ct9/'equivalence-v2/results.json') if a.case=='baseline' else None,contract=sha(a.contract))
    reused={}
    if a.reuse:
        original=n.read(a.reuse/'request-manifest.json')['request'];require(all(original[k]==v for k,v in request.items() if k!='source'),'partial-member reuse request mismatch')
        reused=n.read(a.reuse/'progress.json')['members'];request['explicit_completed_member_reuse']=dict(directory=str(a.reuse),request_sha256=sha(a.reuse/'request-manifest.json'),progress_sha256=sha(a.reuse/'progress.json'),reason='controller-only endpoint tolerance correction; no native artifact or history edited')
    m=RequestManifest(a.output,request)
    out=dict(schema='CT9PhysicalPairV1',ct10_schema='CT10CohortV1',case=c,required_dt=[.001,.0005],target=a.target,members={});initials={};natives={}
    for dt in (.001,.0005):
        name=str(dt);initial_folder=Path(pre['initial']['baseline' if a.case=='baseline' else 'candidate-'+name]['directory']);ini=n.read(initial_folder/'initial.json');initials[name]=ini
        if name in reused:
            row=copy.deepcopy(reused[name]);folder=Path(row['directory']);r=n.read(folder/'case.evolution.json')
            require(sha(folder/'case.evolution.json')==row['native_report_sha256'] and artifact_hash(r['checkpoint'])==row['history_sha256'] and r['dt']==dt and abs(r['time']-a.target)<1e-13,'changed/incomplete explicitly reused native member')
            require(row['actual_initial_sha256']==exact_vector_hash(ini['observation']['state']),'reused member initial mismatch')
            for label in ('native','q48f16','q64f16','q64f32'):require(sha(folder/label/'budget.json')==row['budgets']['budgets'][label]['packet_sha256'],'changed reused physical budget')
            row['explicit_reuse_request_sha256']=sha(a.reuse/'request-manifest.json')
        elif a.case=='baseline' and dt==.001:
            oldrow=copy.deepcopy(base['members'][name]);folder=Path(oldrow['directory']);r=n.read(folder/'case.evolution.json');require(sha(folder/'case.evolution.json')==oldrow['native_report_sha256'],'changed reused sparse baseline report');row=oldrow
            row.update(accepted_metrics=[dict(time=t,**n.physical_metrics(ini,s,c,weights,count)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)],traction_history=[dict(time=t,profile=profiles(ini,s,c)) for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)],actual_initial_traction_profile=profiles(ini,r['accepted_states'][0],c),reuse_receipt=sha(a.ct9/'equivalence-v2/results.json'))
        else:
            if a.case=='baseline':
                zero_folder=a.output/'zero';n.make(zero_folder,a.runtime,c,dt,0);enable(zero_folder);e=run([a.runtime/'IncNavierStokesSolver','case.xml'],zero_folder,'solver.log');zero=n.read(zero_folder/'case.evolution.json');require(e['exit_code']==0 and not zero['failure'],'baseline finer initialization failed')
                previous=n.read(a.ct9/'equivalence-v2/0.001/case.evolution.json')
                for k in BOUNDARY_KEYS:require(zero['initialization_receipt'][k]==previous['initialization_receipt'][k],'baseline finer changed physical initialization:'+k)
            else:zero_folder=initial_folder;zero=n.read(zero_folder/'case.evolution.json')
            require(zero['state']==ini['observation']['state'],'actual zero-step coefficient mismatch')
            source=Path(parent['members'][name]['directory']) if parent else zero_folder;old=n.read(source/'case.evolution.json');require(old['dt']==dt and abs(old['time']-(.008 if parent else 0.))<1e-13,'changed timestep/history lineage')
            folder=a.output/name;n.make(folder,a.runtime,c,dt,round((a.target-old['time'])/dt),Path(old['checkpoint']),source/'case.xml');os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);r=n.read(folder/'case.evolution.json');require(r['accepted_states'][0]==old['state'],'native did not consume unchanged restart state')
            row=collect(a.runtime,folder,ini,r,c,weights,count,zero);row['execution']=e
        out['members'][name]=row;natives[name]=r;write_json(a.output/'progress.json',out);print('CT10',a.case,name,r['time'],r['failure'],flush=True)
    comparison(out,initials,natives,c,weights,count)
    if a.case=='baseline':
        dense=sealed(a.ct9/'pressure-finer-v1');member=dense.get('members',{}).get('0.0005')
        if member is None:member=dense
        # CT9 finer-history results are a complete pair with the retained dense member.
        file=Path(member['directory'])/'case.evolution.json';r=n.read(file);require(sha(file)==member['native_report_sha256'],'changed sealed dense finer member')
        sparse=natives['0.0005'];require(r['accepted_states'][0]==sparse['accepted_states'][0] and r['audit_times']==sparse['audit_times'],'backend transfer time/initial mismatch');eq=[]
        for t,x,y in zip(r['audit_times'],r['accepted_states'],sparse['accepted_states'],strict=True):eq.append(dict(time=t,fields=n.norms(n.observations(initials['0.0005'],x,count),n.observations(initials['0.0005'],y,count),weights,initials['0.0005']['observation']['grid_regions'][:count],c),traction=difference(profiles(initials['0.0005'],x,c),profiles(initials['0.0005'],y,c),c)))
        mx=max(v['scaled'] for q in eq for v in q['fields'].values());mt=max(q['traction']['profiles'][f]['normalized_weighted_RMS'] for q in eq for f in ('ground','interface'));out['dense_fine_transfer']=dict(status='pass' if max(mx,mt)<=1e-8 else 'fail',maximum_scaled_field=mx,maximum_scaled_traction=mt,common_times=eq,dense_report=sha(file),sparse_report=out['members']['0.0005']['native_report_sha256'])
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT10',a.case,out['qualification'],flush=True);return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('operation',choices=['preflight','advance'])
    for key in ('ct8','ct9','runtime','contract','output'):p.add_argument('--'+key,type=Path,required=True)
    for key in ('diagnosis','preflight','parent','reuse'):p.add_argument('--'+key,type=Path)
    p.add_argument('--case',choices=['baseline','candidate']);p.add_argument('--target',type=float,default=.008);a=p.parse_args()
    resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');require(sha(a.runtime/'IncNavierStokesSolver')==RUNTIME_SHA,'foreign governing runtime');n.reference=lru_cache(None)(n.reference)
    q=preflight(a) if a.operation=='preflight' else advance(a);return 0 if q.get('qualification',q).get('status')=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
