#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Alternating corrected-reference/exact-reuse immutable native replays.

Checkpoint and attempt-packet writes remain enabled. The replay driver times
native execution, checks, checkpoint/packet processing and report compression.
It never changes production files or the declared physical/temporal method.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import atomic,header
from profile_cls_awf_campaign import fingerprint,timing_summary
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash,set_info,set_parameter,write

STAGES=('committed','post_flow','pre_k','raw','corrected','history','fixed_reference')

def stable(obj):
    if isinstance(obj,dict):
        return {k:stable(v) for k,v in obj.items() if 'seconds' not in k and k not in
            ('k_snapshot','k_snapshot_sha256','matrix_file','matrix_sha256','profiling','elapsed_wall_time','evaluation_reuse')}
    if isinstance(obj,list):return [stable(v) for v in obj]
    return obj


def replay(fixture,folder,solver,reuse,timeout,steps=None,start=None):
    folder.mkdir(parents=True,exist_ok=False);capture=load(fixture/'capture.json');tree=ET.parse(fixture/'production.xml').getroot()
    initial=start or fixture/'start';n=steps or capture['segment']['end_step']-capture['segment']['start_step']
    set_info(tree,'RANSWallRestart',str(initial));set_info(tree,'RANSWallCheckpoint',str(folder/'checkpoint'))
    set_info(tree,'RANSWallPimpleSnapshotDirectory',str(folder/'pimple-states'));set_info(tree,'RANSWallEvaluationReuse',reuse)
    set_parameter(tree,'NumSteps',n);set_parameter(tree,'IO_CheckSteps',1);write(tree,folder/'case.xml')
    env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(solver.parent/'lib'),NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),TMPDIR=str(folder))
    command=['/usr/bin/time','-f','%e %U %S %M','-o','usage.txt',str(solver),'case.xml']
    before=time.monotonic()
    with (folder/'solver.log').open('w') as log:
        child=subprocess.run(command,cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
    native_wall=time.monotonic()-before;require(child.returncode==0,'native replay failed: '+str(folder))
    report=load(folder/'case.evolution.json');require(not report['failure'],'native replay rejected')
    require(report['audit_revision']=='NativeElementSideV2','wrong audit revision')
    cp=Path(report['checkpoint']);end=fingerprint(cp)
    require(end['steps']==header(initial)['steps']+n,'wrong replay end step')
    packets=[]
    for a in report['pimple_attempts']:
        require(all(a[k]==1 for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls')),'governing schedule changed')
        require(a['additional_governing_equation_solves']==0 and a['scheduled_complete'] and not a['failure'],'additional/failed work')
        packet=Path(a['k_snapshot']);q=load(packet)
        require(sha(packet)==a['k_snapshot_sha256'] and sha(Path(q['matrix_file']))==q['matrix_sha256'],'attempt packet durability/hash differs')
        packets.append(dict(time=a['time'],stages={k:state_hash(q[k]) for k in STAGES},
            external_coefficients=q['fixed_external_coefficients'],matrix_sha256=q['matrix_sha256'],attempt=stable(a)))
    checks={k:report[k] for k in ('stage_audits','state_audits','physical_range_events','physical_samples') if k in report}
    result=dict(status='passed',mode=reuse,solver=str(solver),solver_sha256=sha(solver),xml_sha256=sha(folder/'case.xml'),
        checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),start_checkpoint_sha256=sha(initial/'State.bin'),final=end,
        packets=packets,checks=checks,cache=report.get('evaluation_reuse'),timing=timing_summary(folder/'timing.jsonl'),
        native_wall_seconds=native_wall,command=command)
    # Preserve the same expensive report compression and durable dense packet
    # writes as the production driver. No output optimization is mixed in.
    with (folder/'case.evolution.json').open('rb') as src,gzip.open(folder/'case.evolution.json.gz','wb') as dst:shutil.copyfileobj(src,dst)
    (folder/'case.evolution.json').unlink()
    usage=list(map(float,(folder/'usage.txt').read_text().split()));result.update(native_process_wall_seconds=usage[0],
        native_process_cpu_seconds=usage[1]+usage[2],peak_rss_KiB=usage[3],driver_wall_seconds=time.monotonic()-before)
    atomic(folder/'receipt.json',result);return result


def compare(a,b):
    for name in ('final','checkpoint_sha256','packets','checks'):
        require(a[name]==b[name],f'exact numerical equivalence failed: {name}')
    return dict(bitwise_stages=True,bitwise_native_histories_and_time_labels=True,bitwise_checkpoint_payload=True,
        identical_hard_audits_and_correction=True,driver_speedup=a['driver_wall_seconds']/b['driver_wall_seconds'],
        native_cpu_speedup=a['native_process_cpu_seconds']/b['native_process_cpu_seconds'],
        baseline_driver_seconds=a['driver_wall_seconds'],optimized_driver_seconds=b['driver_wall_seconds'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fixtures',type=Path,required=True)
    p.add_argument('--reference-solver',type=Path,required=True);p.add_argument('--optimized-solver',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--timeout',type=float,default=7200)
    p.add_argument('--cases',nargs='+',default=['BudgetH2','UpstreamOne']);p.add_argument('--pairs',type=int,default=3)
    p.add_argument('--reference-results',type=Path,help='reuse immutable corrected-reference receipts for exactly these fixtures')
    p.add_argument('--long-blocks',action='store_true');a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    out=dict(started_UTC=datetime.now(timezone.utc).isoformat(),status='running',threads=THREADS,parallel_jobs=2,cases={})
    atomic(a.output/'summary.json',out)
    def case(name):
        base=a.output/name;base.mkdir();fixture=a.fixtures/name;data=dict(pairs=[])
        # The optimized executable's forced-Off path must first match the
        # separately frozen corrected uncached implementation.
        if a.reference_results:
            receipt=a.reference_results/name/'reference/receipt.json';ref=load(receipt)
            require(ref['solver_sha256']==sha(a.reference_solver) and ref['start_checkpoint_sha256']==sha(fixture/'start/State.bin') and
                ref['final']==fingerprint(fixture/'expected'),'reused reference identity differs')
            data['corrected_reference_receipt']=dict(path=str(receipt),sha256=sha(receipt))
        else:ref=replay(fixture,base/'reference',a.reference_solver,'Off',a.timeout)
        for i in range(a.pairs):
            runs={}
            for mode in (['Off','ExactStepV1'] if i%2==0 else ['ExactStepV1','Off']):
                runs[mode]=replay(fixture,base/(f'pair-{i+1}-'+mode),a.optimized_solver,mode,a.timeout)
            if i==0:data['new_off_vs_corrected_reference']=compare(ref,runs['Off'])
            data['pairs'].append(compare(runs['Off'],runs['ExactStepV1']))
            atomic(base/'summary.json',data)
        if a.long_blocks:
            long=replay(fixture,base/'long-32',a.optimized_solver,'ExactStepV1',a.timeout,steps=32)
            prior=fixture/'start';pieces=[]
            for i in range(4):
                piece=replay(fixture,base/f'blocks8-{i+1}',a.optimized_solver,'ExactStepV1',a.timeout,steps=8,start=prior)
                pieces.append(piece);prior=Path(piece['checkpoint'])
            require(long['final']==pieces[-1]['final'] and long['checkpoint_sha256']==pieces[-1]['checkpoint_sha256'],'block size changed native endpoint')
            require(long['packets']==[p for q in pieces for p in q['packets']],'block size changed a governing stage/correction')
            data['block_comparison']=dict(steps=32,bitwise_endpoint_and_all_stages=True,
                block32_seconds=long['driver_wall_seconds'],four_blocks8_seconds=sum(q['driver_wall_seconds'] for q in pieces),
                block32_peak_rss_KiB=long['peak_rss_KiB'],blocks8_peak_rss_KiB=max(q['peak_rss_KiB'] for q in pieces),
                long_cache=long['cache'])
        atomic(base/'summary.json',data);return name,data
    with ThreadPoolExecutor(max_workers=2) as pool:
        for name,data in pool.map(case,a.cases):out['cases'][name]=data;atomic(a.output/'summary.json',out)
    out.update(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat());atomic(a.output/'summary.json',out)
    print(a.output/'summary.json')

if __name__=='__main__':main()
