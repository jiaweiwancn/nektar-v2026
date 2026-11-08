#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable primary campaign supervisor; never retries numerical rejection."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from continue_cls_epawf_pimple import CampaignLock,atomic,process_identity,load_progress
from run_cls_epawf_picard import load,require
from run_cls_epawf_pimple import add_bound_work_arguments,bound_work_arguments,bound_work_command


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--wall-policy',choices=('Upstream','UniformEPFlat','CLSCalibratedBudgetAWF'),default='Upstream')
    p.add_argument('--height-case',choices=list('ABCDEFGH'))
    p.add_argument('--upstream-case',choices=('U05-graded','U05-one','U20-graded','U20-one'))
    p.add_argument('--accepted-awf-case',choices=('Accepted-AWF-2G','Accepted-AWF-05G','Accepted-AWF-05O','Accepted-AWF-2O'))
    p.add_argument('--accepted-awf-companion',choices=('faithful','regularized','saturated','regularized-saturated'),default='faithful')
    p.add_argument('--accepted-awf-companion-reason',default='')
    p.add_argument('--accepted-awf-perturbation',type=float,default=0.)
    p.add_argument('--inner-refinement','--accepted-awf-inner-refinement',dest='accepted_awf_inner_refinement',choices=('none','cutoff-v1'),default='none')
    p.add_argument('--inventory',choices=('ConstitutiveV2','ConstitutiveV3','AcceptedBudgetV1'),default='ConstitutiveV2')
    p.add_argument('--bound-solver',choices=('WorkingSet','RetainedPool'),default='WorkingSet')
    p.add_argument('--targets',type=int,nargs='+',default=[32,160,800,1600,3200,8000])
    p.add_argument('--block-steps',type=int,default=16);p.add_argument('--timeout',type=float,default=7200)
    p.add_argument('--checkpoint',type=Path);p.add_argument('--prefix-report',type=Path)
    p.add_argument('--prefix-campaign',type=Path)
    add_bound_work_arguments(p)
    a=p.parse_args();a.output=a.output.resolve()
    require(sum(bool(x) for x in (a.height_case,a.upstream_case,a.accepted_awf_case))<=1,'choose one declared case family')
    if a.accepted_awf_case:
        from run_cls_accepted_awf import WORK
        require(a.wall_policy=='CLSCalibratedBudgetAWF' and a.inventory=='AcceptedBudgetV1' and
            a.bound_solver=='RetainedPool' and bound_work_arguments(a)==WORK,'accepted-AWF family/method controls mismatch')
        require(a.accepted_awf_companion=='faithful' or a.accepted_awf_companion_reason.strip(),'companion needs its preserved faithful blocker')
        require(a.accepted_awf_inner_refinement=='none' or a.accepted_awf_case.endswith('G'),'focused cutoff refinement requires a graded accepted-AWF case')
    else:require(a.wall_policy!='CLSCalibratedBudgetAWF' and a.inventory!='AcceptedBudgetV1' and a.accepted_awf_inner_refinement=='none','new formulation requires its declared case family')
    a.output.mkdir(parents=True,exist_ok=True)
    ownership=a.output/'supervisor';ownership.mkdir(exist_ok=True);lock=CampaignLock(ownership)
    root=a.output/'trajectory';driver=Path(__file__).with_name('continue_cls_epawf_pimple.py')
    base=[sys.executable,'-B',str(driver),'--solver',str(a.solver.resolve()),'--output',str(root),
        '--inventory',a.inventory,'--bound-solver',a.bound_solver,'--wall-policy',a.wall_policy,'--block-steps',str(a.block_steps),'--timeout',str(a.timeout),
        '--targets',*map(str,a.targets),*bound_work_command(bound_work_arguments(a))]
    if a.height_case:base += ['--height-case',a.height_case]
    if a.upstream_case:base += ['--upstream-case',a.upstream_case]
    if a.accepted_awf_case:
        base += ['--accepted-awf-case',a.accepted_awf_case,'--accepted-awf-companion',a.accepted_awf_companion,
                 '--accepted-awf-companion-reason',a.accepted_awf_companion_reason,
                 '--accepted-awf-perturbation',str(a.accepted_awf_perturbation),
                 '--inner-refinement',a.accepted_awf_inner_refinement]
    retries=0
    while True:
        try:progress=load_progress(root/'progress.json') if root.exists() else {}
        except (FileNotFoundError,json.JSONDecodeError):
            atomic(a.output/'status.json',dict(status='waiting_for_progress_publication',process_identity=process_identity(os.getpid())))
            time.sleep(5);continue
        if progress and (a.accepted_awf_case or progress.get('accepted_awf_case')):
            require(progress.get('accepted_awf_case')==a.accepted_awf_case and
                progress.get('accepted_awf_companion','faithful')==a.accepted_awf_companion and
                progress.get('accepted_awf_companion_reason','')==a.accepted_awf_companion_reason and
                progress.get('accepted_awf_perturbation',0.)==a.accepted_awf_perturbation and
                progress.get('accepted_awf_inner_refinement','none')==a.accepted_awf_inner_refinement and
                progress.get('wall_policy')==a.wall_policy and progress.get('inventory')==a.inventory and
                progress.get('bound_solver')==a.bound_solver and
                progress.get('projection_work_controls')==bound_work_arguments(a),
                'accepted-AWF supervisor configuration differs from the existing campaign')
        if progress.get('status') in ('numerically_rejected','targets_completed'):
            atomic(a.output/'status.json',dict(status=progress['status'],progress=str(root/'progress.json'),
                accepted_steps=progress['accepted_steps'],accepted_time=progress['accepted_time'],
                checkpoint=progress['checkpoint'],next_target=progress.get('next_target'),
                failure=progress.get('failure'),supervisor_identity=process_identity(os.getpid())))
            return
        existing=process_identity(progress['pid']) if progress.get('pid') else None
        if existing and str(existing['start'])==str(progress.get('process_start')):
            atomic(a.output/'status.json',dict(status=progress['status'],accepted_steps=progress['accepted_steps'],
                accepted_time=progress['accepted_time'],checkpoint=progress['checkpoint'],next_target=progress.get('next_target'),
                native_identity=progress.get('native_identity'),driver_identity=existing,
                supervisor_identity=process_identity(os.getpid()),attached_to_existing_driver=True))
            time.sleep(5);continue
        if progress.get('native_identity') and process_identity(progress['native_identity']['pid'])==progress['native_identity']:
            # Never terminate or replace a surviving orphan. Its inherited
            # campaign lock is an independent barrier to concurrent mutation.
            atomic(a.output/'status.json',dict(status='waiting_for_recorded_native',native_identity=progress['native_identity'],
                accepted_steps=progress['accepted_steps'],accepted_time=progress['accepted_time']))
            time.sleep(5);continue
        if progress:initial=['--resume']
        elif a.checkpoint:
            require(a.prefix_report is not None,'checkpoint continuation requires a native prefix report')
            initial=['--checkpoint',str(a.checkpoint.resolve()),'--prefix-report',str(a.prefix_report.resolve())]
            if a.prefix_campaign:initial+=['--prefix-campaign',str(a.prefix_campaign.resolve())]
        else:
            require(a.prefix_report is None,'prefix report without checkpoint')
            initial=['--from-h0']
        command=base+initial
        atomic(a.output/'status.json',dict(status='launching' if not progress else 'recovering_infrastructure',command=command,retries=retries))
        with (a.output/'driver.log').open('a') as stream:
            child=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT)
            while child.poll() is None:
                try:state=load_progress(root/'progress.json')
                except (FileNotFoundError,json.JSONDecodeError):time.sleep(1);continue
                if state:
                    atomic(a.output/'status.json',dict(status=state['status'],accepted_steps=state['accepted_steps'],
                        accepted_time=state['accepted_time'],checkpoint=state['checkpoint'],next_target=state.get('next_target'),
                        native_identity=state.get('native_identity'),driver_identity=process_identity(child.pid),
                        supervisor_identity=process_identity(os.getpid())))
                time.sleep(5)
        progress=load_progress(root/'progress.json')
        if progress.get('status') in ('numerically_rejected','targets_completed'):continue
        retries+=1
        if retries>3:
            atomic(a.output/'status.json',dict(status='infrastructure_interrupted',progress=progress,
                reason='three infrastructure recoveries exhausted; no numerical retry'))
            return
        require(progress,'driver exited before durable progress; inspect setup instead of retrying from invented history')
        time.sleep(5)


if __name__=='__main__':main()
