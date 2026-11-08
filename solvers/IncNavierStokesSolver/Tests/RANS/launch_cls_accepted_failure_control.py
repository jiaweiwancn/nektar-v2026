#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One durable, fixed-horizon faithful control. Numerical rejection is terminal."""
import argparse
from datetime import datetime,timezone
import os
from pathlib import Path
import platform
import subprocess
import sys
from run_cls_epawf_picard import THREADS,require,sha,write_json,load
from continue_cls_epawf_pimple import CampaignLock,atomic,process_identity


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--kind',choices=('half-dt','pressure-p2-v1'),required=True)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--solver-sha256',required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--reason',required=True);a=p.parse_args()
    require(platform.system()=='Linux','authoritative fixed-horizon controls run on Linux');os.environ.update(THREADS)
    a.output=a.output.resolve();require(a.output.is_relative_to(Path('/media/jay/新加卷/nektar-v2026-runs')),'isolated heavy run volume required')
    a.output.mkdir(exist_ok=False);lock=CampaignLock(a.output)
    require(sha(a.solver)==a.solver_sha256,'frozen production runtime differs')
    dt=.0003125 if a.kind=='half-dt' else .000625;steps=round(.12/dt)
    command=[sys.executable,'-B',str(Path(__file__).with_name('run_cls_accepted_awf.py')),
        '--solver',str(a.solver.resolve()),'--output',str(a.output/'native'),
        '--cases','Accepted-AWF-2G','--companion','faithful','--inner-refinement','cutoff-v1',
        '--dt',str(dt),'--steps',str(steps),'--checkpoint-every','1','--timeout','21600']
    if a.kind=='pressure-p2-v1':command+=['--outer-refinement','pressure-p2-v1']
    data=dict(schema=1,status='running',started_UTC=datetime.now(timezone.utc).isoformat(),kind=a.kind,reason=a.reason,
        dt=dt,requested_time=.12,requested_steps=steps,command=command,solver_sha256=a.solver_sha256,
        rank=1,threads=THREADS,policy='fresh native h0 and histories; fixed dt; no retry or model switch',
        controller=process_identity(os.getpid()),script_sha256=sha(Path(__file__)))
    env=os.environ.copy();env['LD_LIBRARY_PATH']=str(a.solver.parent/'lib')
    with (a.output/'controller.log').open('w') as stream:
        child=subprocess.Popen(command,stdout=stream,stderr=subprocess.STDOUT,env=env,pass_fds=(lock.stream.fileno(),))
        data.update(child=process_identity(child.pid));atomic(a.output/'status.json',data)
        code=child.wait(timeout=22000)
    resultfile=a.output/'native/Accepted-AWF-2G/result.json'
    data.update(completed_UTC=datetime.now(timezone.utc).isoformat(),exit_code=code)
    if resultfile.exists():
        result=load(resultfile);data.update(status=result['status'],accepted_steps=result.get('steps',0),
            accepted_time=result.get('time',0),failure=result.get('failure'),result=str(resultfile),result_sha256=sha(resultfile))
    else:data.update(status='infrastructure_interrupted',failure='native wrapper did not publish result; inspect process/checkpoints before any recovery')
    atomic(a.output/'status.json',data);print(data['status'])


if __name__=='__main__':main()
