#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused Linux single-sweep experiment: counts, refinement, restart, rollback.

Retains every failed fixed-dt result. Does not substitute numerical completion
for constitutive qualification or for convergence of the original M8 stage.
"""
import argparse
import json
from pathlib import Path
import platform
import subprocess

from run_cls_epawf_picard import load, require, write_json, state_hash, sha
from run_cls_epawf_pimple import ROOT, identity
from compare_cls_epawf_pimple import describe, differences, same_history


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--resolution-probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--comparison-dt',type=float,default=.0025)
    p.add_argument('--comparison-steps',type=int,default=2)
    a=p.parse_args()
    require(a.comparison_dt>0 and a.comparison_steps>=2 and a.comparison_steps%2==0,'even comparison step count required')
    require(platform.system()=='Linux','Linux is authoritative')
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    runner=Path(__file__).with_name('run_cls_epawf_pimple.py')
    summary=dict(schema=1,source=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        executable=identity(a.solver.resolve()),cases={},status='running')
    def save():write_json(a.output/'summary.json',summary)
    def case(name,dt,steps,*flags):
        command=['python3',str(runner),'--solver',str(a.solver.resolve()),
            '--resolution-probe',str(a.resolution_probe.resolve()),'--output',str(a.output/name),
            '--dt',str(dt),'--steps',str(steps),*map(str,flags)]
        with (a.output/(name+'.driver.log')).open('w') as out:
            subprocess.run(command,stdout=out,stderr=subprocess.STDOUT,check=True)
        m=load(a.output/name/'manifest.json');r=load(a.output/name/'case.evolution.json')
        summary['cases'][name]=dict(command=command,status=m['status'],seconds=m['execution']['seconds'],
            accepted_steps=r['steps'],time=r['time'],failure=r['failure'],xml_sha256=sha(a.output/name/'case.xml'))
        save();return r
    save()
    initial=case('startup',.0025,0,'--diagnostic')
    strict=case('strict-dt0025',.0025,3)
    summary['strict_qualification_first_exit']=[{k:v for k,v in e.items() if k not in ('state','history','audit')}
        for e in strict.get('physical_range_events',[])]
    # Independent physical-step runs: smaller dt is never a relaxed field update.
    stable=[];dt=a.comparison_dt
    for attempt in range(5):
        stable=[]
        for level in range(3):
            r=case(f'refine-{attempt}-{level}',dt/2**level,a.comparison_steps*2**level,'--diagnostic')
            if r['failure']:
                break
            stable.append(r)
        if len(stable)==3:break
        dt/=2
    require(len(stable)==3,'three meaningful stable common-time comparisons unavailable; retained failures in summary')
    require(all(r['initialization_lift_calls_iterations']==initial['initialization_lift_calls_iterations'] for r in stable),
        'hidden nonlinear trace lifting during physical trajectory')
    summary['refinement']=[describe(r) for r in stable]
    deltas=[differences(x,y) for x,y in zip(stable,stable[1:])]
    summary['successive_physical_L2_differences']=deltas
    import math
    summary['observed_difference_orders']={k:math.log2(deltas[0][k]/deltas[1][k])
        for k in deltas[0] if k!='pressure_common_gauge_shift' and deltas[0][k]>0 and deltas[1][k]>0}
    half=a.comparison_steps//2
    verification=case('verification-only',dt,half,'--diagnostic','--verify')
    unverified=case('unverified-one',dt,half,'--diagnostic','--checkpoint')
    require(verification['state']==unverified['state'],'verification changed the retained solution')
    resumed=case('restart-one',dt,a.comparison_steps-half,'--diagnostic','--restart',unverified['checkpoint'])
    require(state_hash(resumed['accepted_states'][0])==state_hash(unverified['state']),'restart initial state differs from checkpoint')
    summary['restart_bitwise_equal']=same_history(stable[0],resumed)
    rejected=case('rejected-second',dt,half+1,'--diagnostic','--reject-step',half+1)
    require('deliberate PIMPLE rejection' in rejected['failure'],'did not exercise deliberate rejection')
    summary['rejection_bitwise_equal']=same_history(unverified,rejected)
    require(rejected['implicit_callback_count']==half+1 and rejected['steps']==half,'rejection/callback count')
    require(initial['initial']==stable[0]['initial'] and
        all(r['native_history'][0]==r['state'][:int(r['offsets'][4])] for r in stable),
        'native history did not evolve with differential solution')
    require(state_hash(stable[0]['state'])!=state_hash(initial['state']),'unforced initial residual was cancelled')
    for r in stable:
        previous=r['accepted_states'][-2]
        require(r['last_stage_history'][:int(r['offsets'][4])]==previous[:int(r['offsets'][4])],
            'unforced native callback history differs from previous accepted solution')
    summary.update(status='completed_numerical_experiment',single_sweep_numerical_trajectory=True,
        fully_converged_M8_stage=False)
    save();print(json.dumps({k:summary[k] for k in ('status','observed_difference_orders','restart_bitwise_equal','rejection_bitwise_equal')}))


if __name__=='__main__':main()
