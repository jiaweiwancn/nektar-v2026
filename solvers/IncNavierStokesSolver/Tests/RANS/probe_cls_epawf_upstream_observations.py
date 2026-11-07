#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable actual h0/accepted-state readbacks; query failures are not run failures."""
import argparse
import os
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from run_cls_epawf_picard import load,run,require,write_json,sha,state_hash,THREADS
from continue_cls_epawf_pimple import CampaignLock,atomic
from run_cls_epawf_pimple import identity
from run_cls_epawf_height_comparison import UPSTREAM_CASES

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screening',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cases',nargs='+',choices=list(UPSTREAM_CASES),default=list(UPSTREAM_CASES))
    p.add_argument('--workers',type=int,default=1);a=p.parse_args();require(1<=a.workers<=2,'bounded diagnostic allocation')
    a.output.mkdir(parents=True,exist_ok=True);lock=CampaignLock(a.output);os.environ.update(THREADS)
    summary=dict(probe=identity(a.probe.resolve()),production_governing_solves=0,workers=a.workers,cases={})
    jobs=[]
    for c in a.cases:
        r=a.screening/c
        if not (r/'result.json').exists():summary['cases'][c]=dict(status='waiting_for_completed_native_report');continue
        native=load(r/'result.json');summary['cases'][c]=dict(native_status=native['status'],native_failure=native.get('failure'),observations={})
        if not native.get('initialization'):continue
        initial=load(r/'initial/case.evolution.json');e=load(r/'trajectory/case.evolution.json')
        require(state_hash(e['accepted_states'][0])==state_hash(initial['state']),'readback h0 identity differs')
        for n in (0,8,16,32):
            if n<=e['steps']:jobs.append((c,n,initial['state'] if n==0 else e['accepted_states'][n],n*e['dt'],r/'initial/case.xml'))
    def observe(job):
        c,n,state,t,xml=job;root=a.output/c/f'step-{n:06d}';root.mkdir(parents=True,exist_ok=True);out={}
        for operation in ('height_comparison','source_balance','branch_coverage'):
            folder=root/operation;manifest=folder/'receipt.json'
            if manifest.exists():
                old=load(manifest);require(old['state_sha256']==state_hash(state) and old['probe_sha256']==summary['probe']['sha256'],
                    'readback input/runtime changed; choose a distinct output directory')
                out[operation]=old;continue
            require(not folder.exists(),'unfinished readback requires explicit reconciliation: '+str(folder));folder.mkdir()
            packet=dict(operation=operation,state=state,history=state,time=t,alpha=0.,require_history_h0=False,
                inner_stress_model='SaturatedR2',stations=[-3.5,-2.625,-1.,1.6],
                heights=[0.,.0001,.001,.005,.01,.05,.1,.25,.5,1.,2.])
            write_json(folder/'input.json',packet)
            rc=run([str(a.probe.resolve()),str(xml.resolve()),'input.json','result.json'],folder,'probe.log',1800)
            receipt=dict(**rc,state_sha256=state_hash(state),probe_sha256=summary['probe']['sha256'],time=t,alpha=0.,xml_sha256=sha(xml),
                path=str(folder/'result.json'),sha256=sha(folder/'result.json') if (folder/'result.json').exists() else None,
                diagnostic_only=True,retained_epsilon_trace=True,production_governing_solves=0)
            # Preserve the native outcome even if this optional query fails.
            if rc['exit_code']:receipt['query_failure']=(folder/'probe.log').read_text()[-2500:]
            write_json(manifest,receipt);out[operation]=receipt
        return c,n,out
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(observe,j) for j in jobs]):
            c,n,out=f.result();summary['cases'][c]['observations'][str(n)]=out;atomic(a.output/'summary.json',summary)
    atomic(a.output/'summary.json',summary);print(a.output/'summary.json')

if __name__=='__main__':main()
