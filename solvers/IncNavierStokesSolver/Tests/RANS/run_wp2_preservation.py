#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded fail-fast WP2 time trials, preservation and recovery. No hill runs."""
import argparse,math
from pathlib import Path
from wp2_native import Run
from wp2_qualification import verify_report,quarter_budget_differences
from run_cls_epawf_picard import require,state_hash,sha
from diagnose_cls_accepted_covariance import checkpoint
from run_cao_tamura_budget_wp2 import charge

US=.07426793792419734
K0={'CLS':.017598891430879275,'standard':.01838575534504143}
TREQ=1.1942387693600194


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','h1','preflight','output']:p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True)
    p.add_argument('--mode',choices=['trial','preserve','recover','transaction'],required=True)
    p.add_argument('--level',type=int,default=4);p.add_argument('--dt',type=float,required=True)
    p.add_argument('--endpoint',type=float,default=1.2);p.add_argument('--expected-seed')
    a=p.parse_args();require(a.dt>0,'positive timestep')
    x=Run(a.root,a.solver,a.h1,a.preflight,a.output,a.model,a.level,16,'RadialEnergyV1')
    x.result.update(mode=a.mode,dt=a.dt,Treq=TREQ,physical_endpoint=a.endpoint,qualification='running; no pass inferred from launch')
    x.save();reports=[];parent=None;accepted=0
    try:
        if a.mode=='transaction':
            x.fixture();return
        if a.mode=='trial':
            charge(a.root,'time trial '+a.model+' dt='+str(a.dt),timestep_trials=1)
            targets=[a.endpoint]
        else:
            require(a.endpoint>=TREQ,'physical preservation horizon is too short')
            targets=[.3,.6,.9,a.endpoint]
        if a.mode=='recover':
            require(a.endpoint<=4*TREQ,'recovery cap')
        for target in targets:
            total=round(target/a.dt);require(abs(total*a.dt-target)<1e-13,'requested observation is not a native time')
            while accepted<total:
                count=min(32,total-accepted);label='step'+str(accepted+count)
                d=x.native(label,count,a.dt,parent,perturb=(a.mode=='recover' and parent is None))
                checks=verify_report(d,US,K0[a.model],allow_recovery=a.mode=='recover')
                if parent is None:
                    if a.expected_seed and a.mode!='recover':require(d['wp2_per_state_monitoring'][0]['physical_checks']['state_sha256']==a.expected_seed,'fresh integrator changed qualified seed')
                    if a.mode=='recover':
                        q=d['diagnostic_ivp'];require(q['mode']=='WP2OuterKBumpV1' and q['maximum_actual_relative_k_departure']>=.02,'recovery setup did not exceed twice the field budget')
                        require(state_hash(d['fixed_reference'])==q['immutable_base_sha256'] and state_hash(d['state'])!=q['immutable_base_sha256'],'perturbation replaced the immutable reference')
                if reports:
                    for key in ['fixed_reference','fixed_defect']:require(state_hash(d[key])==state_hash(reports[0][key]),'restart changed '+key)
                accepted=d['steps'];parent=Path(d['checkpoint']);reports.append(d)
                require(abs(d['time']-accepted*a.dt)<=64*math.ulp(max(1.,TREQ)),'physical time drift exceeds declared allowance')
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
