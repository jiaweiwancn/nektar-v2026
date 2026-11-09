#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Cost-bounded independent IMEX1 cases from their actual native h0.

Each dt owns its fresh initialization and GLM history. The middle trajectory is
read from the archived matching campaign, never reinterpreted at another dt.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import gzip
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic,header
from run_cls_accepted_awf import case_xml,check_report
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,write


def run(folder,solver,xml,timeout):
    before=time.monotonic();env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(solver.parent/'lib'),TMPDIR=str(folder))
    with (folder/'solver.log').open('w') as log:
        child=subprocess.run(['/usr/bin/time','-f','%e %U %S %M','-o','usage.txt',str(solver),str(xml)],
            cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    path=folder/'case.evolution.json';r=load(path) if path.exists() else None
    receipt=dict(exit_code=child.returncode,seconds=time.monotonic()-before,report=str(path),
        report_sha256=sha(path) if path.exists() else None,xml_sha256=sha(xml),log_sha256=sha(folder/'solver.log'))
    atomic(folder/'receipt.json',receipt);return r,receipt


def prepare(case,folder,n,dt,restart=None):
    xml=case_xml(case,folder,n,dt,restart=restart,companion='regularized',
        companion_reason='fixed selected model temporal sensitivity; faithful bulk covariance exit archived',
        inner_refinement='cutoff-v1',checkpoint_every=1)
    tree=ET.parse(xml).getroot();set_info(tree,'RANSWallEvaluationReuse','ExactStepV1');write(tree,xml)
    return xml


def campaign(label,dt,args):
    root=args.output/label;root.mkdir();out=dict(status='initializing',dt=dt,model='Regularized/Linear',
        audit_revision='NativeElementSideV2',evaluation_reuse='ExactStepV1',segments=[],observations={},accepted_steps=0,accepted_time=0.)
    atomic(root/'progress.json',out)
    xml=prepare('Accepted-AWF-2G',root/'initial',0,dt);d,rc=run(xml.parent,args.solver,xml,args.timeout)
    require(d is not None and not d['failure'] and rc['exit_code']==0,'native initialization failed')
    check_report(d,False,load(xml.parent/'definition.json'))
    require(d['initial_k_bound_check']['activations']==0,'startup changed native h0')
    h0=state_hash(d['state']);require(h0==args.h0_hash,'physical native h0 differs from matching middle trajectory')
    out.update(status='running',actual_native_h0_sha256=h0,initialization_history_sha256=[state_hash(v) for v in d['native_history']],
        initialization_history_times_sha256=state_hash(d['native_history_times']),checkpoint=d['checkpoint'],
        checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin'),initial_report=rc)
    cp=Path(d['checkpoint']);atomic(root/'progress.json',out)
    for target_time in args.times:
        target=round(target_time/dt);require(abs(target*dt-target_time)<1e-14,'nonintegral physical observation target')
        while out['accepted_steps']<target:
            start=out['accepted_steps'];end=min(start+args.block_steps,target)
            folder=root/f'steps-{start+1:06d}-{end:06d}';xml=prepare('Accepted-AWF-2G',folder,end-start,dt,cp)
            d,rc=run(folder,args.solver,xml,args.timeout)
            if d is None:
                out.update(status='infrastructure_interrupted',failure='native exit without report',last_receipt=rc);atomic(root/'progress.json',out);return label,out
            attempts=check_report(d,True,load(folder/'definition.json'))
            require(state_hash(d['accepted_states'][0])==state_hash(header(cp,True)['coupled_state']),'restart changed selected state')
            cp=Path(d['checkpoint']);out.update(accepted_steps=int(d['steps']),accepted_time=d['time'],checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'))
            out['segments'].append(dict(start_step=start,end_step=d['steps'],time=d['time'],checkpoint=str(cp),receipt=rc,
                correction_integral=sum(a.get('k_bound',{}).get('integral_total',0.) for a in attempts if a['numerically_accepted']),
                bound_activations=sum(a.get('bound_activations',0) for a in attempts),cache=d.get('evaluation_reuse')))
            if rc['exit_code'] or d['failure']:
                out.update(status='numerically_rejected',failure=d['failure']);atomic(root/'progress.json',out);return label,out
            require(d['steps']==end,'completed native block missed target')
            if end==target:out['observations'][str(target_time)]=dict(step=end,time=d['time'],checkpoint=str(cp),
                checkpoint_sha256=sha(cp/'State.bin'),xml=str(xml),report=str(folder/'case.evolution.json.gz'))
            with (folder/'case.evolution.json').open('rb') as src,gzip.open(folder/'case.evolution.json.gz','wb') as dst:shutil.copyfileobj(src,dst)
            (folder/'case.evolution.json').unlink();atomic(root/'progress.json',out)
    out.update(status='targets_completed',completed_UTC=datetime.now(timezone.utc).isoformat());atomic(root/'progress.json',out)
    return label,out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--h0-hash',required=True);p.add_argument('--policy',type=Path,required=True)
    p.add_argument('--times',type=float,nargs='+',default=[.01,.02,.04,.05]);p.add_argument('--block-steps',type=int,default=32)
    p.add_argument('--levels',nargs='+',choices=('coarse','fine','extra-fine'),default=['coarse','fine'])
    p.add_argument('--purpose',default='predeclared primary finite-horizon triplet')
    p.add_argument('--timeout',type=float,default=7200);a=p.parse_args();require(platform.system()=='Linux','Linux numerical environment required')
    require(1<=a.block_steps<=32,'bounded block size required');a.output.mkdir(parents=True,exist_ok=False)
    shutil.copy2(a.policy,a.output/'declared-policy.json');out=dict(status='running',started_UTC=datetime.now(timezone.utc).isoformat(),
        solver_sha256=sha(a.solver),policy_sha256=sha(a.policy),threads=THREADS,independent_jobs=min(2,len(a.levels)),
        requested_times=a.times,purpose=a.purpose,cases={})
    atomic(a.output/'summary.json',out)
    with ThreadPoolExecutor(max_workers=2) as pool:
        levels={'coarse':.00125,'fine':.0003125,'extra-fine':.00015625}
        futures=[pool.submit(campaign,n,levels[n],a) for n in a.levels]
        for future in futures:
            name,result=future.result();out['cases'][name]=result;atomic(a.output/'summary.json',out)
    out.update(status='recorded',completed_UTC=datetime.now(timezone.utc).isoformat());atomic(a.output/'summary.json',out)

if __name__=='__main__':main()
