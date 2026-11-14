#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded read-only small fixtures; never launches a governing solver."""
import argparse
from pathlib import Path
import subprocess
from run_wj_h4 import load
from wj_h15_observer import identity,verify,atomic,assess
from wj_h17_observer import load_compatible,native_load,equal
from study_wj_h17 import completed_request


def preflight(root):
    cache,m=load_compatible(root/'science/small-cache-compatibility.json')
    inv=load(root/'inputs/small/common-support.json');rows,initial,stats=cache.view(inv['points'])
    pin=load(root/'inputs/small/step-2.json');report=native_load(verify(pin['report']))
    result,_=assess(rows,report['state'],initial,2);old=load(verify(pin['result']))
    equal(result['all_regional_Q1_moments'],old['all_regional_Q1_moments'],'small complete Q1 signed moments')
    equal(result['summary'],old['summary'],'small same support exact summary')
    atomic(root/'science/observer-small-preflight.json',dict(status='PASS',result=result,inventory=identity(root/'inputs/small/common-support.json'),cache_stats=stats,governing_work=0))


def component(root):
    folder=root/'science/component-checks';folder.mkdir()
    completion=root/'science/small/A/monitor/step-1.completion.json';c=load(completion)
    anchor=load(root/'inputs/small/step-2.json')
    request=dict(operation='checkpoint_policy_checks',inner_stress_model='Linear',import_role='PendingNativeV1',native_completion=identity(completion),
        native_report=c['report']['path'],native_report_sha256=c['report']['sha256'],native_checkpoint=str(verify(c['checkpoint']).parent),checkpoint_sha256=c['checkpoint']['sha256'],qualified_anchor=completed_request(anchor),startup_qualified_anchor=dict(native_report=c['report']['path'],native_report_sha256=c['report']['sha256'],native_checkpoint=str(Path(c['checkpoint']['path']).parent),checkpoint_sha256=c['checkpoint']['sha256'],completed_sidecar=c['callback'],qualified_cursor=identity(root/'science/small/B/monitor/startup.verification-commit.json')))
    atomic(folder/'request.json',request)
    cmd=[str(root/'runtime/RANSWallResolutionProbe'),str(root/'science/small/A/case.xml'),str(folder/'request.json'),str(folder/'native.json')]
    with(folder/'native.log').open('x')as log:code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    assert code==0,'native guard/pending/qualified-import fixture failed'
    r=load(folder/'native.json');assert r['status']=='PASS'and r['pending_integrity_rejections']==11 and r['full_audit_calls']==0
    assert r['initialization_work']==[0,0]and r['governing_solves']==0
    atomic(root/'science/component-checks.json',dict(status='PASS',native=identity(folder/'native.json'),fixture='actual shared StopOnBound and scalar domain;11 pending integrity refusals; qualified importer remains strict',no_target_audit_batch=True))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('preflight','component'),required=True);a=p.parse_args();globals()[a.mode](a.root)
