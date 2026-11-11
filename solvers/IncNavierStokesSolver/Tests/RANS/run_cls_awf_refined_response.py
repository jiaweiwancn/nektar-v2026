#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded refined-P1 late inner-k response; seeds never reset epsilon/reference.

The native ConstrainedInnerMassV1 seed operation is setup algebra. Production
blocks use ordinary compatible restart and the unchanged one-sweep step map.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import platform
import shutil
import time

from continue_cls_epawf_pimple import CampaignLock, atomic, process_identity
from diagnose_cls_accepted_covariance import checkpoint
from run_cls_awf_late_ivp import native, packet, prepare
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash


def same_history(a,b):
    return ([state_hash(x) for x in a['native_history']]==[state_hash(x) for x in b['native_history']]
        and state_hash(a['native_history_times'])==state_hash(b['native_history_times']))


def ready(folder,steps,dt,**kw):
    return prepare(folder,steps,dt,refinement='cutoff-v1',**kw)


def validate(args,parent,lock):
    root=args.output/'validation';root.mkdir();q=checkpoint(parent);dt=q['dt']
    ordinary,ro=native(args.solver,ready(root/'ordinary',1,dt,restart=parent),args.timeout,lock)
    imported,ri=native(args.solver,ready(root/'zero-import',1,dt,parent=parent,inner_k=0.),args.timeout,lock)
    left,right=packet(ordinary),packet(imported)
    checks={name:state_hash(left[name])==state_hash(right[name]) for name in
        ('committed','post_flow','pre_k','raw','corrected','history','fixed_reference')}
    checks.update(zero_seed_bits=state_hash(imported['accepted_states'][0])==state_hash(q['coupled']),
        zero_reference_bits=state_hash(imported['fixed_reference'])==state_hash(q['reference']),
        zero_defect_bits=state_hash(imported['fixed_defect'])==state_hash(q['defect']),
        zero_final_bits=state_hash(imported['state'])==state_hash(ordinary['state']),
        zero_histories_times_bits=same_history(imported,ordinary),
        external_coefficients=left['fixed_external_coefficients']==right['fixed_external_coefficients'],
        zero_lift=imported['initialization_lift_calls_iterations']==[0,0],
        zero_projection_noop=imported['diagnostic_ivp']['inner_k_seed']['mass_solves']==0 and
            imported['diagnostic_ivp']['inner_k_seed']['changed_k_coefficients']==0)
    seed,rs=native(args.solver,ready(root/'perturbed-seed',0,dt,parent=parent,inner_k=.001),args.timeout,lock)
    identity=seed['diagnostic_ivp']['inner_k_seed'];meta=seed.get('accepted_awf_final') or seed['accepted_awf_initial']
    sizes=list(map(int,meta['native_degrees_of_freedom']));off=sizes[0]+sizes[1];n=sizes[2]
    changed=[i for i,(a,b) in enumerate(zip(seed['state'],q['coupled'])) if a!=b]
    checks.update(only_k_changed=bool(changed) and all(off<=i<off+n for i in changed),
        actual_increment_bits=state_hash([a-b for a,b in zip(seed['state'][off:off+n],q['coupled'][off:off+n])])==identity['actual_increment_sha256'],
        boundary_outer_matching_unchanged=all(identity[k]==0 for k in ('fixed_coefficient_increment_max','outer_value_increment_max')) and
            identity.get('matching_native_value_change_max',0)==0 and identity['matching_trace_increment_max']<=identity.get('matching_trace_evaluation_roundoff',0),
        reference_defect_exact=state_hash(seed['fixed_reference'])==state_hash(q['reference']) and state_hash(seed['fixed_defect'])==state_hash(q['defect']),
        seed_lift_zero=seed['initialization_lift_calls_iterations']==[0,0],
        initial_bound_noop=seed['initial_k_bound_check']['activations']==0,
        own_initial_history=state_hash(seed['native_history'][0])==state_hash(seed['state'][:len(seed['native_history'][0])]))
    whole,rw=native(args.solver,ready(root/'whole2',2,dt,parent=parent,inner_k=.001),args.timeout,lock)
    first,rf=native(args.solver,ready(root/'split1',1,dt,restart=Path(seed['checkpoint'])),args.timeout,lock)
    second,rr=native(args.solver,ready(root/'split2',1,dt,restart=Path(first['checkpoint'])),args.timeout,lock)
    checks.update(restart_final_bits=state_hash(whole['state'])==state_hash(second['state']),restart_history_bits=same_history(whole,second))
    rejected,re=native(args.solver,ready(root/'forced-rejection',1,dt,restart=Path(seed['checkpoint']),hard_reject=1),args.timeout,lock,False)
    checks.update(rollback_state=state_hash(rejected['state'])==state_hash(seed['state']),rollback_history=same_history(rejected,seed),
        rollback_time=rejected['time']==seed['time'] and rejected['steps']==0,
        rejected_candidate_retained=bool(rejected['rejected_state']) and state_hash(rejected['rejected_state'])!=state_hash(seed['state']))
    bad={}
    for label,xml in [
        ('ordinary_changed_dt',ready(root/'bad-dt',0,dt/2,restart=parent)),
        ('ordinary_changed_mesh',prepare(root/'bad-mesh',0,dt,restart=parent,refinement='none')),
        ('import_bad_parent_hash',prepare(root/'bad-parent-hash',0,dt,parent=parent,refinement='cutoff-v1',inner_k=0.,bad_hash=True))]:
        _,rc=native(args.solver,xml,args.timeout,lock,False);bad[label]=rc
    require(all(checks.values()),'refined response contract failed: '+str(checks))
    result=dict(status='passed',checks=checks,seed_identity=identity,incompatible_requests=bad,
        parent_sha256=sha(parent/'State.bin'),receipts=dict(ordinary=ro,zero_import=ri,seed=rs,whole=rw,first=rf,second=rr,rejection=re))
    atomic(root/'summary.json',result);return result


def level(label,amplitude,args,parent,lock):
    root=args.output/label;root.mkdir();q=checkpoint(parent);dt=q['dt'];started=time.monotonic()
    out=dict(status='initializing',amplitude_K0=amplitude,accepted_steps=0,elapsed_time=0.,physical_t0=q['time'],dt=dt,
        segments=[],observations={},parent_sha256=sha(parent/'State.bin'))
    atomic(root/'progress.json',out)
    report,rc=native(args.solver,ready(root/'initial',0,dt,parent=parent,inner_k=amplitude),args.timeout,lock)
    cp=Path(report['checkpoint']);out.update(status='running',seed_checkpoint=str(cp),seed_checkpoint_sha256=sha(cp/'State.bin'),
        seed_identity=report['diagnostic_ivp'],seed_receipt=rc)
    atomic(root/'progress.json',out)
    try:
        for elapsed in args.times:
            target=round(elapsed/dt);require(abs(target*dt-elapsed)<1e-13,'nonintegral response endpoint')
            while out['accepted_steps']<target:
                require(time.monotonic()-started<args.max_hours*3600,'declared response wall-cost limit reached')
                start=out['accepted_steps'];end=min(start+32,target);folder=root/f'steps-{start+1:06d}-{end:06d}'
                report,rc=native(args.solver,ready(folder,end-start,dt,restart=cp),args.timeout,lock)
                require(state_hash(report['accepted_states'][0])==state_hash(checkpoint(cp)['coupled']),'restart seed changed')
                cp=Path(report['checkpoint']);require(report['steps']==end and abs(report['time']-(q['time']+end*dt))<1e-10,'response clock/count changed')
                out.update(accepted_steps=end,elapsed_time=end*dt,physical_time=report['time'],checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'))
                attempts=report['pimple_attempts'];out['segments'].append(dict(start_step=start,end_step=end,receipt=rc,
                    bound_activations=sum(t['bound_activations'] for t in attempts),
                    correction_integral=sum(t['k_bound']['integral_total'] for t in attempts),
                    governing_counts={k:sum(t[k] for t in attempts) for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls')},
                    additional_governing_solves=sum(t['additional_governing_equation_solves'] for t in attempts)))
                if end==target:out['observations'][str(elapsed)]=dict(checkpoint=str(cp),xml=str(folder/'case.xml'),
                    report=str(folder/'case.evolution.json.gz'),elapsed_time=end*dt,physical_time=report['time'])
                # Native checkpoints and first/last packet per block survive.
                # Deletion occurs only AFTER complete native success, with hashes.
                with (folder/'case.evolution.json').open('rb') as src,gzip.open(folder/'case.evolution.json.gz','wb') as dst:shutil.copyfileobj(src,dst)
                (folder/'case.evolution.json').unlink()
                with (root/'discarded-packets.jsonl').open('a') as tomb:
                    for t in attempts[1:-1]:
                        p=Path(t['k_snapshot']);v=load(p);tomb.write(json.dumps(dict(path=str(p),sha256=sha(p),matrix_sha256=v['matrix_sha256']))+'\n')
                        Path(v['matrix_file']).unlink();p.unlink()
                atomic(root/'progress.json',out)
        out.update(status='targets_completed',seconds=time.monotonic()-started);atomic(root/'progress.json',out)
    except Exception as error:
        out.update(status='interrupted_or_rejected',failure=str(error),seconds=time.monotonic()-started)
        atomic(root/'progress.json',out);raise
    return label,out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('solver','parent','policy','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--times',type=float,nargs='+',default=[.05,.1,.25,.5]);p.add_argument('--timeout',type=float,default=3600)
    p.add_argument('--max-hours',type=float,default=8);p.add_argument('--validate-only',action='store_true')
    a=p.parse_args();require(platform.system()=='Linux','authoritative native Linux execution required')
    require(a.times==sorted(set(a.times)) and 0<max(a.times)<=1.25,'bounded chronological response targets required')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output)
    parent=a.output/'parent';shutil.copytree(a.parent,parent);require(sha(parent/'State.bin')==sha(a.parent/'State.bin'),'parent copy changed')
    shutil.copy2(a.policy,a.output/'declared-policy.json');os.environ.update(THREADS)
    out=dict(status='validating',started_UTC=datetime.now(timezone.utc).isoformat(),driver_identity=process_identity(os.getpid()),
        parent_source=str(a.parent),parent_sha256=sha(parent/'State.bin'),solver_sha256=sha(a.solver),script_sha256=sha(Path(__file__)),
        policy_sha256=sha(a.policy),threads=THREADS,workers=2,cases={})
    atomic(a.output/'summary.json',out);out['validation']=validate(a,parent,lock);atomic(a.output/'summary.json',out)
    if a.validate_only:out['status']='validation_completed'
    else:
        out['status']='running';atomic(a.output/'summary.json',out)
        with ThreadPoolExecutor(max_workers=2) as pool:
            jobs={name:pool.submit(level,name,amp,a,parent,lock) for name,amp in [('control',0.),('perturbed',.001)]}
            for name,job in jobs.items():
                try:_,out['cases'][name]=job.result()
                except Exception as error:out['cases'][name]=dict(status='failed',reason=str(error))
                atomic(a.output/'summary.json',out)
        out['status']='targets_completed' if all(v['status']=='targets_completed' for v in out['cases'].values()) else 'terminal_failure'
    out['completed_UTC']=datetime.now(timezone.utc).isoformat();atomic(a.output/'summary.json',out)


if __name__=='__main__':main()
