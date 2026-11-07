#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable read-only milestone analysis; never launches/retries a solver."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from continue_cls_epawf_pimple import atomic,load_progress
from run_cls_epawf_picard import THREADS,load,require


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--uniform',type=Path,required=True)
    p.add_argument('--step',type=int,default=80);p.add_argument('--timeout',type=float,default=7200)
    a=p.parse_args();r=a.root.resolve();source=Path(__file__).parent
    status=r/f'milestone{a.step}.status.json';start=time.monotonic();os.environ.update(THREADS)
    os.environ['LD_LIBRARY_PATH']=str(r/'probe-budget-runtime/lib')
    def save(**values):atomic(status,dict(pid=os.getpid(),updated_unix=time.time(),requested_step=a.step,**values))
    try:
        while True:
            d=load_progress(r/'repaired-campaign/trajectory/progress.json')
            if d['accepted_steps']>=a.step:
                records=[json.loads(x) for x in (r/'repaired-campaign/trajectory/steps.jsonl').read_text().splitlines()]
                row=next(v for v in records if v['global_step']==a.step)
                if (Path(row['segment'])/'receipt.json').exists():break
            if d['status']=='numerically_rejected' and d['accepted_steps']<a.step:
                save(status='numerically_rejected_before_milestone',accepted_steps=d['accepted_steps'],failure=d['failure']);return
            require(time.monotonic()-start<a.timeout,'milestone observer timeout; trajectory status remains separate')
            save(status='waiting_for_accepted_native_milestone',accepted_steps=d['accepted_steps'],accepted_time=d['accepted_time'])
            time.sleep(10)
        save(status='analysing_immutable_native_milestone',accepted_steps=d['accepted_steps'],accepted_time=d['accepted_time'])
        def run(name,args):
            with (r/f'milestone{a.step}-{name}.log').open('w') as log:
                subprocess.run([sys.executable,'-B',str(source/name),*map(str,args)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
        run('analyse_cls_epawf_constitutive_v2.py',['--campaign',r/'repaired-campaign','--probe',r/'probe-budget-runtime/RANSWallResolutionProbe',
            '--output',r/f'readback{a.step}','--through-step',a.step,'--steps',64,73,a.step])
        run('analyse_cls_epawf_projection_repair.py',['--root',r,'--uniform',a.uniform,'--output',r/'repair-verification.json'])
        run('plot_cls_epawf_projection_repair.py',['--root',r,'--uniform',a.uniform,'--readback',r/f'readback{a.step}','--output',r/'figures'])
        run('capture_cls_epawf_projection_repair.py',['--root',r,'--uniform',a.uniform,'--output',r/f'evidence-after{a.step}'])
        save(status='milestone_accepted_and_analysed',native= str(r/f'evidence-after{a.step}/native.json'))
    except BaseException as error:
        save(status='diagnostic_interrupted; trajectory status is separate',failure=repr(error));raise


if __name__=='__main__':main()
