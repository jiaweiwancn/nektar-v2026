#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One optional read-only diagnostic worker beside the supervised trajectory.

Only completed native blocks are analysed. A failure is recorded separately;
this monitor never edits numerical state, restarts the solver, or changes dt.
"""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import os
from continue_cls_epawf_pimple import CampaignLock,atomic,process_identity,load_progress
from run_cls_epawf_picard import load,THREADS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--resume',action='store_true');a=p.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=a.resume)
    lock=CampaignLock(a.output);root=a.campaign.resolve()/'trajectory';done=set();receipts=[]
    if a.resume:
        previous=load_progress(a.output/'status.json');done=set(map(int,previous.get('completed_diagnostics',[])));receipts=previous.get('receipts',[])
        for folder in a.output.glob('step-??????'):
            if (folder/'summary.json').exists():done.add(int(folder.name.split('-')[-1]))
    while True:
        try:progress=load_progress(root/'progress.json')
        except (FileNotFoundError,json.JSONDecodeError):
            time.sleep(5);continue
        accepted=int(progress['accepted_steps'])
        candidates=[v for v in range(32,accepted+1,32) if v not in done]
        terminal=progress['status'] in ('numerically_rejected','targets_completed')
        if terminal and accepted and accepted not in done:candidates.append(accepted)
        if not candidates:
            atomic(a.output/'status.json',dict(status='finished' if terminal else 'waiting_for_fixed_checkpoint',
                accepted_steps=accepted,accepted_time=progress['accepted_time'],completed_diagnostics=sorted(done),receipts=receipts,
                process_identity=process_identity(os.getpid())))
            if terminal:return
            time.sleep(10);continue
        n=min(candidates)
        text=(root/'steps.jsonl').read_text();lines=text.splitlines()
        if text and not text.endswith('\n'):lines=lines[:-1]
        ledger=[json.loads(v) for v in lines]
        row=next(v for v in ledger if v['global_step']==n);segment=Path(row['segment'])
        if not (segment/'receipt.json').exists():time.sleep(5);continue
        folder=a.output/f'step-{n:06d}'
        command=[sys.executable,'-B',str(Path(__file__).with_name('analyse_cls_epawf_constitutive_v2.py')),
            '--campaign',str(a.campaign.resolve()),'--probe',str(a.probe.resolve()),'--output',str(folder),
            '--through-step',str(n),'--steps',*map(str,sorted({v for v in (n-32,n) if v>0}))]
        atomic(a.output/'status.json',dict(status='analysing_immutable_prefix',through_step=n,command=command,process_identity=process_identity(os.getpid())))
        with (a.output/f'step-{n:06d}.log').open('w') as log:
            try:code=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,**THREADS),timeout=7200).returncode
            except subprocess.TimeoutExpired:code=-999
        receipts.append(dict(step=n,exit_code=code,command=command,
            status='analysis_completed' if code==0 else 'diagnostic_failed; trajectory status is separate'))
        done.add(n)


if __name__=='__main__':main()
