#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read the protected terminal again after actual owned recovery deletion."""
import argparse,os,subprocess
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from wj_h17_observer import native_load


def run(root):
    d=root/'resumption-r1/science/C-sparse-restart-v1'
    rotation=load(d/'rotation.json');assert rotation['status']=='PASS'and rotation['removed']
    q=load(d/'qualification.json');out=d/'strict-after-rotation';out.mkdir()
    terminal=q['terminal'];commit=load(d/f'monitor/step-{terminal}.qualification-commit.json')
    for k in ('checkpoint','callback','report'):verify(commit[k])
    request=load(d/'strict/request.json')
    assert request['checkpoint_sha256']==commit['checkpoint']['sha256']and request['native_report_sha256']==commit['report']['sha256']
    atomic(out/'request.json',request)
    runtime=load(verify(load(d/'execution.json')['runtime']))
    probe=verify(next(e['wrapper']for e in runtime['executables']if Path(e['wrapper']['path']).name=='RANSWallResolutionProbe'))
    with(out/'native.log').open('x')as f:
        subprocess.run([str(probe),str(d/'case.xml'),str(out/'request.json'),str(out/'native.json')],stdout=f,stderr=subprocess.STDOUT,check=True)
    r=native_load(out/'native.json');assert r['initialization_work']==[0,0]and r['initial_velocity_calls']==r['governing_solves']==0
    assert r['selected_absolute_step']==terminal and r['producer_failure']==''
    assert r['weak_checks']['alpha_mass']==200 and r['scalar_domain']['admissible']
    atomic(out/'qualification.json',dict(status='PASS_POST_ROTATION_TERMINAL_READ',terminal=terminal,
        time=native_load(verify(commit['report']))['time'],request=identity(out/'request.json'),native=identity(out/'native.json'),
        rotation=identity(d/'rotation.json'),protected_commit=identity(d/f'monitor/step-{terminal}.qualification-commit.json'),
        genuine_alpha=r['weak_checks']['alpha_mass'],physical_mass_dual_norms=r['physical_mass_dual_norms'],
        weak_checks=r['weak_checks'],scalar_domain=r['scalar_domain'],
        scope='fresh strict QualifiedV1 read after deletion; genuine producing callback/history, zero initialization and zero governing solves'))
    print('protected terminal strict read after rotation PASS',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
