#!/usr/bin/env python3
"""Full-history/rank-change and map-mismatch checks for the R1 Newton map."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
from cao_tamura_experiment import (ROOT,THREADS,allocation,environment,execute,records,run,sha,write_json)
from cao_tamura_routes import coverage
from cao_tamura_cv3_budgets import flat_momentum_budget
from terrain_continuity import audit as strip_continuity
from validate_terrain import compare
from generate_m4 import info,parameter


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('build','prefix','recipe','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--cpus',required=True);a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    os.environ.update(THREADS);allocation(4,a.cpus);os.sched_setaffinity(0,set(map(int,a.cpus.split(','))))
    build=environment(a.build);out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(__file__,out/'producer.py')
    first=a.cpus.split(',')[0];helper=Path(__file__).parent/'cao_tamura_r1.py'
    recipe=json.loads(a.recipe.read_text());recipe['config'].update(steps=7,checkpoint=1)
    write_json(out/'whole-recipe.json',recipe)
    execute([sys.executable,helper,'pilot','--surface','R','--model','cls_fixed','--build',build,
        '--output',out/'whole','--cpus',first,'--timeout','900','--override',out/'whole-recipe.json'],out,'whole',1200)
    restart=ROOT/'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py'
    resumed=out/'resumed';resumed.mkdir()
    execute([sys.executable,restart,a.prefix/'case.xml',a.prefix/'case.fld',resumed/'case.xml','--steps','4'],out,'prepare')
    tree=ET.parse(resumed/'case.xml');parameter(tree.getroot(),'IO_CheckSteps',1)
    # The serial reference uses a direct solver. MPI continuation requires a
    # parallel-capable inverse; this changes the inverse, not the stage equations.
    info(tree.getroot(),'GlobalSysSoln','IterativeMultiLevelStaticCond')
    tree.write(resumed/'case.xml',encoding='utf-8',xml_declaration=True)
    manifest=json.loads((a.prefix/'manifest.json').read_text());manifest['case_sha256']=sha(resumed/'case.xml')
    manifest.setdefault('numerics',{})['GlobalSysSoln']='IterativeMultiLevelStaticCond'
    manifest['config'].update(steps=4,checkpoint=1);write_json(resumed/'manifest.json',manifest)
    receipt=run(resumed/'case.xml',build,4,a.cpus,900)
    s=receipt['stages'];free=max(float(r['free_coefficient_l2']) for r in records(resumed/'case.rans-residual.csv') if r['phase']=='1')
    if (free>1e-7 or min(s['minimum_k'],s['minimum_epsilon'])<=0 or s['maximum_guard_activations']
            or max(s[k] for k in ('maximum_weak_divergence','maximum_projected_normal',
                'maximum_normal_integral_defect','maximum_resolved_mass_defect'))>1e-7):
        raise ValueError('Continued native stage failed its physical gates')
    for name,folder in (('whole',out/'whole'),('resumed',resumed)):
        execute(allocation(1,first)+[build/'RANSWallFieldAudit',folder/'case.xml',folder/'case.fld',out/(name+'.csv')],out,'full-'+name)
    errors=compare(out/'whole.csv',out/'resumed.csv')
    if max(errors.values())>1e-6:raise ValueError('Full-history 1-to-4-rank fields differ: '+str(errors))
    mismatches={'map':('RANSNewtonMap','CurrentSolve')}
    if recipe.get('numerics',{}).get('RANSNewtonScalarScaling')=='Pointwise':
        mismatches['scalar-scaling']=('RANSNewtonScalarScaling','FieldMaximum')
    for name,(key,value) in mismatches.items():
        rejected=out/(name+'-mismatch');rejected.mkdir()
        tree=ET.parse(resumed/'case.xml');info(tree.getroot(),key,value)
        tree.write(rejected/'case.xml',encoding='utf-8',xml_declaration=True)
        try:execute(allocation(1,first)+[build/'IncNavierStokesSolver','case.xml'],rejected,'reject',120)
        except RuntimeError:
            if 'RANS restart model/configuration mismatch' not in (rejected/'reject.log').read_text():
                raise ValueError('Wrong '+name+'-mismatch rejection')
        else:raise ValueError('A mode-specific history silently changed '+name)
    result=dict(status='PASS',version='cao-tamura-r1-residual-map-history-v1',run=receipt,
        maximum_free_residual=free,scaled_field_errors=errors,map_mismatch_rejected=True,
        scalar_scaling_mismatch_rejected='scalar-scaling' in mismatches,
        coverage=coverage(resumed/'case.xml'),strip_continuity=strip_continuity(resumed,'case'),
        strip_momentum=flat_momentum_budget(resumed/'case.xml'),
        binaries={n:sha(build/n) for n in ('IncNavierStokesSolver','RANSWallFieldAudit')},
        cpus=a.cpus,threads=THREADS)
    result['producer_sha256']=sha(out/'producer.py')
    write_json(out/'evidence.json',result);print('residual-map history/MPI/mismatch PASS',errors,flush=True)


if __name__=='__main__':main()
