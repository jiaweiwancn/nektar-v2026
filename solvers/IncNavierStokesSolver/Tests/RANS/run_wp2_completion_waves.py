#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable, bounded WP2 confirmation/recovery waves after the half-step legs.

Every ordinary phase must complete before its successor starts. This controller
never restarts a failed trajectory and never changes an executable or checkpoint.
"""
import argparse,json,os,subprocess,time
from pathlib import Path
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import require,sha,THREADS
from run_cao_tamura_budget_wp2 import charge
from wp2_all_state import Rows


def completed(path):
    d=load(path)
    require(d['status']!='terminal_stop','predecessor failed: '+str(path))
    if d['status']!='completed_bounded_trajectory':return False
    require(d['actual_time']>=1.1942387693600194 and d['final_physical_budgets_passed'],
            'predecessor has not passed its bounded physical horizon')
    return True


def main():
    p=argparse.ArgumentParser()
    for n in ['root','solver','probe','h1','preflight']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();r=a.root;record=r/'completion-waves.json'
    require(not record.exists(),'completion controller already has an owner')
    state=dict(version='WP2CompletionWavesV1',status='waiting_for_half_step',
               solver=str(a.solver),solver_sha256=sha(a.solver),probe=str(a.probe),probe_sha256=sha(a.probe),events=[])
    atomic(record,state);env=dict(os.environ,**THREADS)
    def event(name,**values):
        state['events'].append(dict(name=name,epoch=time.time(),**values));atomic(record,state)
    try:
        while not all(completed(r/f'preservation-fine-{m}/progress.json') for m in ['CLS','standard']):
            require(time.time()-load(r/'origin.json')['started_epoch']<43200,'WP2 aggregate wall envelope');time.sleep(15)
        event('both_half_step_legs_completed')
        # Geometry-only rows are constructed with no other WP2 native case
        # running. They do not initialize/evaluate a candidate or advance time.
        common=load(r/'preservation-coarse-standard/step60/case.evolution.json')['wp2_common_samples']['physical_z_partition']
        catalogs={}
        for name,xml in [('selected',r/'preservation-coarse-standard/step60/case.xml'),
                         ('confirmation',r/'confirmation-trial-standard/step2/case.xml')]:
            charge(r,'immutable retrospective geometry rows '+name,geometry_context_builds=1)
            packet=dict(operation='wp2_geometry_rows',inner_stress_model='Linear',xml_sha256=sha(xml),expected_model='StandardKEpsilon',z_partition=common)
            out=r/('geometry-rows-'+name)
            _,receipt=run_probe(a.probe,xml,packet,out)
            catalogs[name]=Rows(out/'result.json');event('geometry_rows_'+name,receipt=receipt)
        checks={}
        for m in ['CLS','standard']:
            checks[m]=dict(selected=catalogs['selected'].check_native(load(r/f'preservation-coarse-{m}/step60/case.evolution.json')),
                confirmation=catalogs['confirmation'].check_native(load(r/f'confirmation-trial-{m}/step2/case.evolution.json')))
        atomic(r/'geometry-row-replay.json',checks);event('both_model_native_tuple_replays_passed')
        for phase,mode,level,dt,parent in [('confirmation','preserve',5,.02,'confirmation-trial'),
                                         ('recovery','recover',4,.01,'recovery-seed')]:
            jobs=[];state['status']='running_'+phase;atomic(record,state)
            for m in ['CLS','standard']:
                output=r/f'preservation-{phase}-{m}';require(not output.exists(),'phase output already exists')
                command=['/usr/bin/python3',str(Path(__file__).with_name('run_wp2_preservation.py')),
                    '--root',str(r),'--solver',str(a.solver),'--h1',str(a.h1),'--preflight',str(a.preflight),
                    '--output',str(output),'--model',m,'--mode',mode,'--level',str(level),'--dt',str(dt),'--endpoint','1.2',
                    '--parent-progress',str(r/f'{parent}-{m}/progress.json')]
                atomic(r/f'preservation-{phase}-{m}-command.json',command)
                with (r/f'preservation-{phase}-{m}.log').open('x') as log:
                    job=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
                jobs.append((m,job,output));event('launched_'+phase,model=m,pid=job.pid,command=command)
            while any(job.poll() is None for _,job,_ in jobs):time.sleep(15)
            for m,job,output in jobs:
                require(job.returncode==0,'phase failed '+phase+' '+m)
                require(completed(output/'progress.json'),'phase did not complete '+phase+' '+m)
                event('completed_'+phase,model=m,progress=str(output/'progress.json'),progress_sha256=sha(output/'progress.json'))
        state['status']='completed_bounded_waves'
    except Exception as e:
        state.update(status='terminal_stop',failure=str(e));raise
    finally:atomic(record,state)


if __name__=='__main__':main()
