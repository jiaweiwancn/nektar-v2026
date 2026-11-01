#!/usr/bin/env python3
"""Execute the declared CV3 design with explicit bounded Linux CPU workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import fcntl
import json
import math
import os
from pathlib import Path
from queue import Queue,Empty
import sys
import threading

sys.dont_write_bytecode=True
from cao_tamura_experiment import ROOT,THREADS,config,environment,sha,write_json
from cao_tamura_calibration import VERSION,design,develop,select

BASE=dict(geometry='floor',diffusion='Molecular',dimension=2,degree=3,pressure_degree=2,
    layers=6,outer_blend_cells=3,x_columns=[-1.,0.,1.2],dt=.05,steps=1,checkpoint=1,
    velocity_scaling='incident_ReH',k_profile='MeasuredRms',epsilon_profile='ShearLength',
    inlet_wall_collar=True,inlet_reference_q=10/3,inlet_reference_kappa=math.sqrt(.48*1.3/(10/3)),
    advection='Implicit',startup_flow_times=0.)
NUMERICS=dict(GlobalSysSoln='DirectMultiLevelStaticCond',RANSNonlinearSolver='NewtonKrylov',RANSGradDiv=.01)
TERMINAL={'developed_provisional_2D','failed_or_timed_out','failed_input_or_audit','not_developed_within_budget'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--build',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cpus',required=True)
    p.add_argument('--controller-cpu',type=int,required=True);p.add_argument('--timeout',type=int,default=7200)
    p.add_argument('--max-extensions',type=int,default=20);p.add_argument('--resume',action='store_true')
    a=p.parse_args();cpus=list(map(int,a.cpus.split(',')))
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    if (not cpus or len(set(cpus))!=len(cpus) or a.controller_cpu in cpus or len(cpus)+1>16):
        raise ValueError('At most 15 distinct case CPUs and one separate controller CPU')
    os.environ.update(THREADS);os.sched_setaffinity(0,{a.controller_cpu});build=environment(a.build)
    output=a.output.resolve()
    if output.is_relative_to(ROOT):raise ValueError('Store runs outside the source checkout')
    output.mkdir(parents=True,exist_ok=a.resume)
    lock=(output/'campaign.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    recipe=dict(version=VERSION,base=BASE,numerics=NUMERICS,timeout=a.timeout,max_extensions=a.max_extensions,
        binary_sha256=sha(build/'IncNavierStokesSolver'),geometry_sha256=sha(build/'geometry-lib/libSpatialDomains.so'),
        source_sha256={n:sha(Path(__file__).parent/n) for n in ('cao_tamura_campaign.py','cao_tamura_calibration.py','cao_tamura_experiment.py','cao_tamura_cv3_budgets.py','terrain_continuity.py')},
        cpus=cpus,controller_cpu=a.controller_cpu,threads=THREADS,maximum_workers=len(cpus)+1,
        priority='k amplitude high, roughness low, epsilon high, then the remaining original design; order does not enter selection')
    path=output/'recipe.json'
    if path.exists() and json.loads(path.read_text())!=recipe:raise ValueError('Use a new run directory for a changed recipe or binary')
    write_json(path,recipe)
    tasks=[]
    for surface in ('S','R'):
        for c in design(surface):
            for model in ('standard_original','cls_original'):
                tasks.append((surface,c,model))
    priority={'k_scale-high':0,'ks-low':1,'epsilon_scale-high':2,'centre':3}
    tasks.sort(key=lambda t:(priority.get(t[1]['id'],4),t[1]['id'],t[0],t[2]))
    queue=Queue()
    for surface,c,model in tasks:
        folder=output/surface/c['id']/model
        saved=folder/'development.json'
        if saved.exists():
            status=json.loads(saved.read_text())
            if status['status'] not in TERMINAL:raise ValueError('Interrupted case requires explicit history recovery: '+str(folder))
            continue
        queue.put((surface,c,model))
    print_lock=threading.Lock()
    def worker(cpu):
        os.sched_setaffinity(0,{cpu})
        while True:
            try:surface,c,model=queue.get_nowait()
            except Empty:return
            folder=output/surface/c['id']/model
            value=config({**BASE,**c['parameters'],'surface':surface,'model':model})
            with print_lock:print('START',surface,c['id'],model,'CPU',cpu,flush=True)
            try:
                state=develop(value,folder,build,1,str(cpu),a.timeout,a.max_extensions,NUMERICS)
            except Exception as error:
                folder.mkdir(parents=True,exist_ok=True);path=folder/'development.json'
                state=json.loads(path.read_text()) if path.exists() else dict(version=VERSION,config=value,stages=[])
                if state.get('status') not in TERMINAL:state['status']='failed_input_or_audit'
                state.update(complete=False,error=repr(error));write_json(path,state)
            with print_lock:print('END',surface,c['id'],model,state['status'],'stages',len(state.get('stages',[])),flush=True)
            queue.task_done()
    with ThreadPoolExecutor(max_workers=len(cpus)) as pool:
        futures=[pool.submit(worker,cpu) for cpu in cpus]
        for future in futures:future.result()
    write_json(output/'selection.json',select(output))
    print('All declared design attempts concluded; selection eligibility remains separately gated.',flush=True)


if __name__=='__main__':main()
