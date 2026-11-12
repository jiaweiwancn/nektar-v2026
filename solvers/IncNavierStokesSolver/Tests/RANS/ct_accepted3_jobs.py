#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Sequential fixed-case milestones. Invoke inside ONE owned CT3 command."""
import argparse, json, os, subprocess, sys
from pathlib import Path
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
from ct_accepted3_run import CASES
from rans_request_manifest import RequestManifest, sha, runtime_identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for x in ('root','ct1','ct2','output'):p.add_argument('--'+x,type=Path,required=True)
    p.add_argument('--case',choices=CASES,required=True)
    p.add_argument('--initial',type=Path);p.add_argument('--parent',type=Path)
    p.add_argument('--first-target',type=float,choices=[.008,.02,.2],default=.008)
    a=p.parse_args();runtime=a.ct2/'runtime';s=Path(__file__).parent
    contract=a.root/'source/docs/plans/cao-tamura-accepted-3-evidence/contract.json'
    base=[sys.executable,'-B',str(s/'ct_accepted3_run.py'),'--case',a.case,'--runtime',str(runtime),'--ct1',str(a.ct1),'--contract',str(contract)]
    manifest=RequestManifest(a.output,dict(schema='CT3BoundedFixedCaseJobsV1',case=a.case,source_sha256=sha(Path(__file__)),
        runtime=runtime_identity(runtime/'IncNavierStokesSolver'),first_target=a.first_target,
        existing_initial_sha256=sha(a.initial/'results.json') if a.initial else None,
        parent_sha256=sha(a.parent/'results.json') if a.parent else None))
    result=dict(case=a.case,status='incomplete',jobs=[]);initial=a.initial or a.output/'initial';parent=a.parent
    commands=[]
    if a.initial is None:commands.append(('initial',base+['--mode','initial','--output',str(initial)]))
    for target in (.008,.02,.2):
        if target<a.first_target:continue
        destination=a.output/('T'+str(target).replace('.','p'))
        cmd=base+['--mode','advance','--initial',str(initial),'--target',str(target),'--output',str(destination)]
        if parent:cmd+=['--parent',str(parent)]
        commands.append((str(target),cmd));parent=destination
    for label,cmd in commands:
        with (a.output/(label+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
        destination=Path(cmd[cmd.index('--output')+1]);file=destination/'results.json'
        row=dict(label=label,command=cmd,exit_code=code,result_file=str(file),result_sha256=sha(file) if file.exists() else None)
        result['jobs'].append(row);write_json(a.output/'progress.json',result)
        if code:
            result['status']='stopped_at_failed_or_unresolved_gate';result['first_stop']=row;break
    else:result['status']='completed_required_horizons'
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json'])
    print(a.case,result['status']);return 0 if result['status']=='completed_required_horizons' else 1
if __name__=='__main__':raise SystemExit(main())
