#!/usr/bin/env python3
"""Explicit, unfitted CV3-R1 composite-floor feasibility operations.

One invocation owns one case allocation. No parameter selection or calibration
freeze is performed here. Use new output directories for every recipe change.
"""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
from cao_tamura_experiment import (ROOT,THREADS,allocation,config,environment,execute,
    generate,prepare,records,run,sha,write_csv,write_json)
from cao_tamura_calibration import develop,probes
from cao_tamura_routes import coverage
from generate_m3 import element
from generate_m4 import info,parameter

BASE=dict(geometry='floor',wall_route='EPTransition',diffusion='Molecular',dimension=2,
    degree=4,pressure_degree=2,quadrature_points=8,layers=12,outer_blend_cells=3,
    x_columns=[-1.,-.6,-.32,-.28,-.26,-.24,-.22,-.2,-.16,-.1,0.,.4,1.2],
    inlet_delta=.21,k_scale=1.,epsilon_scale=1.,ramp_end=-.2,ramp_width_h=2.,
    dt=.01,steps=6,checkpoint=6,velocity_scaling='incident_ReH',
    k_profile='MeasuredRms',epsilon_profile='ShearLength',inlet_wall_collar=True,
    inlet_reference_q=10/3,inlet_reference_kappa=math.sqrt(.48*1.3/(10/3)),
    advection='Implicit',startup_flow_times=0.)
NUMERICS=dict(GlobalSysSoln='DirectMultiLevelStaticCond',RANSNonlinearSolver='NewtonKrylov',
    RANSGradDiv=.01,IterativeSolverTolerance=1e-8,RANSMatrixFreeTolerance=1e-11,
    PressureInverseTolerance=1e-13)


def numerical_settings(root,numerics):
    for key,val in numerics.items():
        if key=='PressureInverseTolerance':
            group=root.find('CONDITIONS/GLOBALSYSSOLNINFO')
            if group is None:group=element(root.find('CONDITIONS'),'GLOBALSYSSOLNINFO')
            field=element(group,'V',VAR='p')
            element(field,'I',PROPERTY='IterativeSolverTolerance',VALUE=str(val))
        else:(info if key in ('RANSNonlinearSolver','GlobalSysSoln') else parameter)(root,key,val)
    parameter(root,'RANSStageAudit',1)


def route_probes(surface,path):
    rows=probes(surface,path)
    for x in (-.6,-.28,-.26,-.24,-.22,-.2,-.16,-.1,0.):
        for z in (.011,.022,.04,.1,.2,.3):
            rows.append(dict(id=f'route:{x:g}:{z:g}',x=x,span=0.,z=z))
    write_csv(path,rows,['id','x','span','z'])
    receipt=json.loads(Path(str(path)+'.json').read_text())
    receipt.update(scope='CV3-R1-unfitted-feasibility-no-reserved-data',csv_sha256=sha(path))
    write_json(str(path)+'.json',receipt)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('operation',choices=('pilot','develop'))
    p.add_argument('--surface',choices=('S','R'),required=True)
    p.add_argument('--model',choices=('standard_fixed','cls_fixed'),required=True)
    p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpus',required=True);p.add_argument('--ranks',type=int,default=1)
    p.add_argument('--timeout',type=int,default=7200);p.add_argument('--extensions',type=int,default=20)
    p.add_argument('--checkpoint-interval',type=int,default=5)
    p.add_argument('--override',type=Path,help='Preregistered JSON configuration and numerical overrides')
    a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    allocation(a.ranks,a.cpus);os.environ.update(THREADS)
    os.sched_setaffinity(0,set(map(int,a.cpus.split(','))))
    override=json.loads(a.override.read_text()) if a.override else {}
    c=config({**BASE,**override.get('config',{}),'surface':a.surface,'model':a.model,
              'ks':override.get('ks',0. if a.surface=='S' else .004)})
    numerics={**NUMERICS,**override.get('numerics',{})};build=environment(a.build)
    out=a.output.resolve()
    if out.is_relative_to(ROOT) or out.exists():raise ValueError('Use a fresh run directory outside the source tree')
    receipt=dict(version='cao-tamura-cv3-r1-v1',operation=a.operation,config=c,numerics=numerics,
        parent_revision=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source_sha256={n:sha(Path(__file__).parent/n) for n in
            ('cao_tamura_r1.py','cao_tamura_experiment.py','cao_tamura_calibration.py','cao_tamura_routes.py',
             'cao_tamura_cv3_budgets.py','terrain_continuity.py')},
        binaries={n:sha(build/n) for n in ('IncNavierStokesSolver','RANSWallFieldAudit','RANSObservationAudit')},
        geometry_sha256=sha(build/'geometry-lib/libSpatialDomains.so'),ranks=a.ranks,cpus=a.cpus,threads=THREADS,
        override_sha256=sha(a.override) if a.override else None,status='RUNNING',calibration_complete=False)
    receipt['checkpoint_interval']=a.checkpoint_interval
    out.parent.mkdir(parents=True,exist_ok=True)
    write_json(str(out)+'.start.json',receipt)
    try:
        if a.operation=='pilot':
            m=generate(c,out);write_json(out/'operation.json',receipt)
            tree=ET.parse(out/'case.xml');numerical_settings(tree.getroot(),numerics)
            tree.write(out/'case.xml',encoding='utf-8',xml_declaration=True)
            m.update(case_sha256=sha(out/'case.xml'),numerics=numerics);write_json(out/'manifest.json',m)
            receipt['preparation']=prepare(out/'case.xml',build)
            receipt['run']=run(out/'case.xml',build,a.ranks,a.cpus,a.timeout)
            receipt['coverage']=coverage(out/'case.xml')
            receipt['status']='PILOT_COMPLETE_NOT_DEVELOPED'
        else:
            # The same settings are applied by the generic staged developer.
            receipt['development']=develop(c,out,build,a.ranks,a.cpus,a.timeout,a.extensions,
                numerics,probe_builder=route_probes,stage_auditor=coverage,score=False,initial_unforced_flow_times=1,
                checkpoint_interval=a.checkpoint_interval)
            receipt['status']='DEVELOPED' if receipt['development']['complete'] else 'NOT_DEVELOPED'
    except Exception as error:
        receipt.update(status='FAIL',error=repr(error));out.mkdir(parents=True,exist_ok=True)
        state=out/'development.json'
        if state.exists():
            value=json.loads(state.read_text())
            value.update(complete=False,error=repr(error))
            if value['status']=='developing':value['status']='failed_input_or_audit'
            write_json(state,value)
        write_json(out/'operation.json',receipt);raise
    if receipt['binaries']!={n:sha(build/n) for n in receipt['binaries']}:raise ValueError('Binary changed during operation')
    write_json(out/'operation.json',receipt);print(a.surface,a.model,a.operation,receipt['status'],flush=True)


if __name__=='__main__':main()
