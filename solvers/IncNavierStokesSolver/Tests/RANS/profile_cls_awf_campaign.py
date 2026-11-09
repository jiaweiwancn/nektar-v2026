#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Time immutable completed campaign blocks without touching production paths.

Capture each block's actual XML and start/end native checkpoints. Replay with
an opt-in timed solver, then require bitwise state and native-history equality
with the archived production endpoint. All native hard checks stay enabled.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import gzip
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET

from continue_cls_epawf_pimple import atomic, freeze, header
from run_cls_epawf_picard import THREADS, load, require, sha, state_hash, set_info, set_parameter, write
from run_cls_epawf_pimple import check_report, identity


def fingerprint(checkpoint):
    h=header(checkpoint,True)
    return dict(steps=h['steps'],time=h['time'],dt=h['dt'],order=h['order'],
                coupled_state_sha256=state_hash(h['coupled_state']),
                native_history_sha256=[state_hash(v) for v in h['history']],
                native_history_times_sha256=state_hash(h['history_times']))


def timing_summary(path):
    records=[json.loads(line) for line in path.read_text().splitlines()]
    by_id={r['id']:r for r in records}
    require(len(by_id)==len(records),'duplicate timing ID')
    children=defaultdict(list)
    for r in records:
        require(r['parent_id']==0 or r['parent_id'] in by_id,'missing timing parent')
        require(r['wall_seconds']>=0 and r['process_cpu_seconds']>=0,'negative timing')
        children[r['parent_id']].append(r)
    for r in records:
        for total,exclusive in [('wall_seconds','exclusive_wall_seconds'),
                                ('process_cpu_seconds','exclusive_process_cpu_seconds')]:
            child_sum=sum(c[total] for c in children[r['id']])
            require(abs(r[total]-child_sum-r[exclusive])<.0001,'timing tree does not partition')
    steps=[r for r in records if r['name']=='native.step']
    require(steps and len({r['step'] for r in steps})==len(steps),'invalid physical step timing')
    groups=defaultdict(lambda:dict(calls=0,wall=0.,cpu=0.,exclusive_wall=0.,exclusive_cpu=0.))
    for r in records:
        chain=[];p=r
        while p['parent_id']:
            p=by_id[p['parent_id']];chain.append(p['name'])
        section='step' if r['name']=='native.step' or 'native.step' in chain else 'block_overhead'
        g=groups[section+'/'+r['name']];g['calls']+=1
        for field,key in [('wall_seconds','wall'),('process_cpu_seconds','cpu'),
                          ('exclusive_wall_seconds','exclusive_wall'),('exclusive_process_cpu_seconds','exclusive_cpu')]:
            g[key]+=r[field]
    for g in groups.values():
        for key in ('wall','cpu','exclusive_wall','exclusive_cpu'):g[key+'_per_step']=g[key]/len(steps)
    return dict(physical_steps=len(steps),tree_partition_passed=True,groups=dict(groups))


def capture(campaign,output):
    p=load(campaign/'progress.json');segments=p['segments'];require(len(segments)>1,'no complete replay block')
    end=segments[-1];previous=segments[-2]
    require(previous['end_step']==end['start_step'] and end['exit_code']==0,'noncontiguous/completed block')
    output.mkdir(parents=True,exist_ok=False)
    for name,source in [('start',Path(previous['checkpoint'])),('expected',Path(end['checkpoint']))]:
        (output/name).mkdir();shutil.copy2(source/'State.bin',output/name/'State.bin')
    directory=Path(end['checkpoint']).parent
    shutil.copy2(directory/'case.xml',output/'production.xml')
    initial=header(output/'start');final=header(output/'expected')
    require(initial['steps']==end['start_step'] and final['steps']==end['end_step'],'checkpoint steps differ')
    recent=[]
    for segment in segments[-4:]:
        metrics=Path(segment['checkpoint']).parent/'accepted-metrics.json.gz'
        with gzip.open(metrics,'rt') as stream:report=json.load(stream)
        recent.append(dict(segment=segment,metrics_path=str(metrics),metrics_sha256=sha(metrics),
            sweeps=[dict(seconds=r['seconds'],bound=r.get('k_bound',{}).get('seconds_total'),
                         bound_activations=r.get('bound_activations'),bound_iterations=r.get('bound_active_set_iterations'))
                    for r in (record['attempt'] for record in report)]))
    result=dict(captured_UTC=datetime.now(timezone.utc).isoformat(),campaign=str(campaign),
                case_definition=p.get('case_definition'),
                production_runtime=p['executable'],production_executable_sha256=p['executable_sha256'],
                segment=end,recent_completed_blocks=recent,start=fingerprint(output/'start'),
                expected=fingerprint(output/'expected'),
                start_checkpoint_sha256=sha(output/'start/State.bin'),expected_checkpoint_sha256=sha(output/'expected/State.bin'),
                production_xml_sha256=sha(output/'production.xml'))
    atomic(output/'capture.json',result);return result


def replay(directory,solver,timeout,timed=True):
    capture=load(directory/'capture.json');run=directory/('timed' if timed else 'untimed');run.mkdir(exist_ok=False)
    root=ET.parse(directory/'production.xml').getroot()
    # Only restart/output paths and the completed block length change. No
    # physical selector, inventory, work control, quadrature or source changes.
    set_info(root,'RANSWallRestart',str(directory/'start'))
    set_info(root,'RANSWallCheckpoint',str(run/'checkpoint'))
    set_info(root,'RANSWallPimpleSnapshotDirectory',str(run/'pimple-states'))
    set_parameter(root,'NumSteps',capture['segment']['end_step']-capture['segment']['start_step'])
    write(root,run/'case.xml')
    trace=run/'timing.jsonl';env=dict(os.environ,**THREADS,NEKTAR_RANS_TIMING_FILE=str(trace))
    if not timed:env.pop('NEKTAR_RANS_TIMING_FILE',None)
    # RUSAGE_CHILDREN is process-wide. Use it only for the single-child
    # untimed control; timed concurrent replays use native per-process clocks.
    usage=resource.getrusage(resource.RUSAGE_CHILDREN) if not timed else None
    started=time.monotonic()
    with (run/'solver.log').open('w') as stream:
        result=subprocess.run([str(solver),'case.xml'],cwd=run,env=env,stdout=stream,stderr=subprocess.STDOUT,timeout=timeout)
    native_seconds=time.monotonic()-started
    receipt=dict(exit_code=result.returncode,native_wall_seconds=native_seconds,
                 log_sha256=sha(run/'solver.log'),xml_sha256=sha(run/'case.xml'))
    if usage:
        after=resource.getrusage(resource.RUSAGE_CHILDREN)
        receipt['child_cpu_seconds']=after.ru_utime+after.ru_stime-usage.ru_utime-usage.ru_stime
    result_path=directory/('result.json' if timed else 'result-untimed.json')
    atomic(result_path,receipt)
    require(result.returncode==0,'native timing replay failed: '+str(directory))
    # Reading/parsing the full report is outside the native process timer.
    started=time.monotonic();report=load(run/'case.evolution.json');receipt['python_report_load_seconds']=time.monotonic()-started
    started=time.monotonic()
    if report.get('wall_policy')=='CLSCalibratedBudgetAWF':
        from run_cls_accepted_awf import check_report as accepted_check_report
        accepted_check_report(report,True,capture['case_definition'])
    else:check_report(report,True)
    receipt['python_report_checks_seconds']=time.monotonic()-started
    receipt['steps']=report['steps'];receipt['time']=report['time']
    got=dict(steps=report['steps'],time=report['time'],dt=report['dt'],order=report['order'],
             coupled_state_sha256=state_hash(report['state']),native_history_sha256=[state_hash(v) for v in report['native_history']],
             native_history_times_sha256=state_hash(report['native_history_times']))
    require(got==capture['expected'],'timed state/history differ from production block')
    receipt['bitwise_production_state_and_all_native_histories_equal']=True
    receipt['final']=got
    if timed:receipt['timing']=timing_summary(trace);receipt['trace_sha256']=sha(trace)
    else:require(not trace.exists(),'untimed control unexpectedly wrote a trace')
    receipt['attempts']=[{k:r.get(k) for k in ('time','seconds','outer_sweeps','momentum_predictor_solves',
        'pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls','bound_activations',
        'bound_active_set_iterations','additional_governing_equation_solves')} for r in report['pimple_attempts']]
    receipt['bound_work']=[{k:r['k_bound'].get(k) for k in ('activations','active_set_iterations',
        'separation_passes','inventory_rebuilds','mass_factor_applications','kkt_solves','seconds_total','seconds_algebra',
        'seconds_factorization','seconds_separation','seconds_inventory')} for r in report['pimple_attempts']]
    started=time.monotonic()
    with (run/'case.evolution.json').open('rb') as source,gzip.open(run/'case.evolution.json.gz','wb') as target:shutil.copyfileobj(source,target)
    receipt['python_report_gzip_seconds']=time.monotonic()-started
    (run/'case.evolution.json').unlink();atomic(result_path,receipt)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--campaign',action='append',default=[],help='name=/absolute/campaign/path')
    p.add_argument('--solver',type=Path)
    p.add_argument('--untimed-control',help='Replay an already captured case with the same frozen binary, timing disabled')
    p.add_argument('--parallel',type=int,default=1,choices=(1,2))
    p.add_argument('--timeout',type=float,default=1800)
    a=p.parse_args();require(platform.system()=='Linux','Linux authoritative timings required')
    a.output.mkdir(parents=True,exist_ok=True)
    for spec in a.campaign:
        name,path=spec.split('=',1);require(Path(name).name==name,'unsafe case name')
        capture(Path(path),a.output/name)
    if a.untimed_control:
        require(not a.solver and not a.campaign,'untimed control uses the existing frozen runtime/capture')
        require(Path(a.untimed_control).name==a.untimed_control,'unsafe control case name')
        solver=a.output/'runtime/IncNavierStokesSolver';pinned=load(a.output/'runtime/identity.json')['frozen']
        os.environ['LD_LIBRARY_PATH']=str(a.output/'runtime/lib')
        require(sha(solver)==pinned['sha256'],'frozen control executable changed')
        require(identity(solver)['libraries']==pinned['libraries'],'frozen control library closure changed')
        replay(a.output/a.untimed_control,solver,a.timeout,False)
        return
    if a.solver:
        solver,pinned=freeze(a.solver.resolve(),a.output/'runtime')
        cases=sorted(d for d in a.output.iterdir() if (d/'capture.json').exists())
        with ThreadPoolExecutor(max_workers=a.parallel) as pool:
            results=list(pool.map(lambda d:replay(d,solver,a.timeout),cases))
        atomic(a.output/'summary.json',dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),
            source_revision=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip(),
            threads=THREADS,rank=1,parallel_replays=a.parallel,runtime=pinned,
            cases={d.name:r for d,r in zip(cases,results)}))


if __name__=='__main__':main()
