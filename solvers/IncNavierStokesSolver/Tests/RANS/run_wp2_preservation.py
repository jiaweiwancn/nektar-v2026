#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded fail-fast WP2 time trials, preservation and recovery. No hill runs."""
import argparse,math,time
from pathlib import Path
from wp2_native import Run
from wp2_qualification import verify_report,quarter_budget_differences,physical_time_tolerance
from run_cls_epawf_picard import require,state_hash,sha
from diagnose_cls_accepted_covariance import checkpoint
from run_cao_tamura_budget_wp2 import charge
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic

US=.07426793792419734
K0={'CLS':.017598891430879275,'standard':.01838575534504143}
TREQ=1.1942387693600194


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','h1','preflight','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True)
    p.add_argument('--mode',choices=['trial','recovery-trial','preserve','recover','transaction','seed'],required=True)
    p.add_argument('--level',type=int,default=4);p.add_argument('--dt',type=float,required=True)
    p.add_argument('--block-steps',type=int,choices=[8,16,32])
    p.add_argument('--handoff-reason')
    p.add_argument('--endpoint',type=float,default=1.2);p.add_argument('--expected-seed')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--parent-progress',type=Path)
    p.add_argument('--pressure-space',choices=['baseline-v1','inner-p2-v1'],default='baseline-v1')
    a=p.parse_args();require(a.dt>0,'positive timestep')
    block_steps=a.block_steps or (16 if a.pressure_space=='inner-p2-v1' else 32)
    recovery=a.mode in ['recover','recovery-trial']
    x=Run(a.root,a.solver,a.h1,a.preflight,a.output,a.model,a.level,16,'RadialEnergyV1',a.resume,a.pressure_space)
    previous=[]
    require(not(a.resume and a.parent_progress),'choose same-runtime resume or an explicit successor, not both')
    if a.parent_progress:
        require(a.mode in ['preserve','recover'],'successor phase must be bounded preservation/recovery')
        prefix=load(a.parent_progress)
        require(prefix['status']=='completed_bounded_trajectory' and prefix['dt']==a.dt and prefix['model']==a.model and
                prefix['case_definition']==x.result['case_definition'],'successor parent configuration differs')
        previous=prefix['observations'];x.result['prefix_observations']=previous
        x.result['predecessor']=dict(progress=str(a.parent_progress),progress_sha256=sha(a.parent_progress),
            solver=prefix['solver'],solver_sha256=prefix['solver_sha256'],reason=a.handoff_reason or 'explicit exact-checkpoint successor after diagnostic-only runtime revalidation')
    if a.resume:
        require(a.mode in ['preserve','recover'] and x.result['dt']==a.dt,'incompatible continuation phase/dt')
        previous=x.result['observations']
        require(a.endpoint>x.result['actual_time'],'continuation endpoint must advance')
        old=a.output/('phase-before-'+str(x.result['accepted_steps'])+'.json')
        require(not old.exists(),'phase record already exists');atomic(old,x.result)
        x.result.setdefault('phase_history',[]).append(dict(path=str(old),sha256=sha(old),resumed_epoch=time.time()))
    x.result.update(mode=a.mode,dt=a.dt,block_steps=block_steps,Treq=TREQ,physical_endpoint=a.endpoint,qualification='running; no pass inferred from launch')
    x.result['status']='running';x.save();reports=[];parent=None;accepted=0
    if previous:
        accepted=previous[-1]['accepted'];parent=Path(previous[-1]['checkpoint'])
        require(sha(parent/'State.bin')==previous[-1]['checkpoint_sha256'],'continuation checkpoint changed')
        reports=[load(Path(previous[0]['report']))]
    try:
        if a.mode=='transaction':
            x.fixture(a.dt);return
        if a.mode=='seed':
            d=x.native('h0',0,a.dt,perturb=True);checks=verify_report(d,US,K0[a.model],allow_recovery=True)
            require(d['diagnostic_ivp']['maximum_actual_relative_k_departure']>=.02,'actual recovery seed is unresolved')
            x.result.update(status='completed_bounded_trajectory',accepted_steps=0,actual_time=0.,
                observations=[dict(time=0.,accepted=0,checkpoint=d['checkpoint'],checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin'),
                    report=x.result['results']['h0']['report'],physical_budgets_passed=False,final_physical_checks=checks[-1])])
            return
        if a.mode in ['trial','recovery-trial']:
            charge(a.root,'time trial '+a.model+' dt='+str(a.dt),timestep_trials=1)
            targets=[a.endpoint]
        else:
            require(a.endpoint>=TREQ,'physical preservation horizon is too short')
            targets=sorted(set([.3,.6,.9,a.endpoint-.3,a.endpoint]))
        if recovery:
            require(a.endpoint<=4*TREQ,'recovery cap')
        for target in targets:
            total=round(target/a.dt);require(abs(total*a.dt-target)<1e-13,'requested observation is not a native time')
            while accepted<total:
                count=min(block_steps,total-accepted);label='step'+str(accepted+count)
                d=x.native(label,count,a.dt,parent,perturb=(recovery and parent is None))
                checks=verify_report(d,US,K0[a.model],allow_recovery=recovery)
                if previous and 'prefix_essential_absolute_error' not in x.result:
                    table=d['wp2_per_state_monitoring'][0]['physical_checks'];ids=[int(i) for i in table['physical_essential_rows']];values=table['physical_essential_values']
                    error=0.
                    for observation in previous:
                        old=load(Path(observation['report']))
                        error=max(error,max(abs(state[i]-v) for state in old['accepted_states'] for i,v in zip(ids,values)))
                    require(error==0,'predecessor changed immutable physical prescribed coefficients')
                    x.result['prefix_essential_absolute_error']=error
                if parent is None:
                    if a.expected_seed:require(d['wp2_per_state_monitoring'][0]['physical_checks']['state_sha256']==a.expected_seed,'fresh integrator changed qualified seed')
                    if recovery:
                        q=d['diagnostic_ivp'];require(q['mode']=='WP2OuterKBumpV1' and q['maximum_actual_relative_k_departure']>=.02,'recovery setup did not exceed twice the field budget')
                        require(state_hash(d['fixed_reference'])==q['immutable_base_sha256'] and state_hash(d['state'])!=q['immutable_base_sha256'],'perturbation replaced the immutable reference')
                if reports:
                    for key in ['fixed_reference','fixed_defect']:require(state_hash(d[key])==state_hash(reports[0][key]),'restart changed '+key)
                accepted=d['steps'];parent=Path(d['checkpoint']);reports.append(d)
                require(abs(d['time']-accepted*a.dt)<=physical_time_tolerance(accepted,a.dt,TREQ),'physical time drift exceeds declared allowance')
                x.result.setdefault('observations',[]).append(dict(time=d['time'],accepted=accepted,checkpoint=str(parent),checkpoint_sha256=sha(parent/'State.bin'),
                    report=x.result['results'][label]['report'],physical_budgets_passed=all(v['passed'] for v in checks),
                    final_physical_checks=checks[-1],seed_sha256=d['wp2_per_state_monitoring'][0]['physical_checks']['state_sha256']))
                x.save()
        x.result.update(status='completed_bounded_trajectory',accepted_steps=accepted,actual_time=reports[-1]['time'],
                        transit_units=reports[-1]['time']/TREQ,reference_defect_preserved=True,
                        final_physical_budgets_passed=x.result['observations'][-1]['final_physical_checks']['passed'])
    except Exception as e:x.result.update(status='terminal_stop',failure=str(e));raise
    finally:x.save()


if __name__=='__main__':main()
