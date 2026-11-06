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


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--targets',type=int,nargs='+',default=[32,160,800,1600,3200,8000])
    p.add_argument('--block-steps',type=int,default=16);p.add_argument('--timeout',type=float,default=7200)
    a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=True)
    ownership=a.output/'supervisor';ownership.mkdir(exist_ok=True);lock=CampaignLock(ownership)
    root=a.output/'trajectory';driver=Path(__file__).with_name('continue_cls_epawf_pimple.py')
    base=[sys.executable,'-B',str(driver),'--solver',str(a.solver.resolve()),'--output',str(root),
        '--inventory','ConstitutiveV2','--block-steps',str(a.block_steps),'--timeout',str(a.timeout),
        '--targets',*map(str,a.targets)]
    retries=0
    while True:
        try:progress=load_progress(root/'progress.json') if root.exists() else {}
        except (FileNotFoundError,json.JSONDecodeError):
            atomic(a.output/'status.json',dict(status='waiting_for_progress_publication',process_identity=process_identity(os.getpid())))
            time.sleep(5);continue
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
        command=base+(['--resume'] if progress else ['--from-h0'])
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
