#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Short restart/cache/rollback/refinement controls, never production campaigns."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
from run_a93_a95_e94e95_matrix import case_xml,run_case
from run_cls_epawf_picard import THREADS,sha,state_hash,set_info,set_parameter,write
from continue_cls_epawf_pimple import header
from compare_cls_epawf_pimple import same_history


def load(path):return json.loads(Path(path).read_text())

def check(ok,message):
    if not ok:raise AssertionError(message)

def checkpoint(folder,steps):
    candidates=[p for p in Path(folder).glob('checkpoint*') if p.is_dir() and (p/'State.bin').exists() and header(p)['steps']==steps]
    check(bool(candidates),('checkpoint missing',steps,folder));return sorted(candidates)[0]

def physical_difference(a,b):
    result={}
    for k,v in a['physical_samples'].items():
        if not k.startswith('values_'):continue
        if k not in b['physical_samples']:continue
        w=a['physical_samples']['weights_'+k[7:]]
        diff=[x-y for x,y in zip(v,b['physical_samples'][k])]
        result[k]=math.sqrt(math.fsum(q*x*x for q,x in zip(w,diff))/math.fsum(w))
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--descriptor-probe',type=Path,required=True);p.add_argument('--matrix',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);os.environ.update(THREADS)
    results=dict(schema=1,checks=[],cases=[],status='running',threads=THREADS,ranks=1)
    def save():(a.output/'qualification.json').write_text(json.dumps(results,indent=2)+'\n')
    def control(name,model='E94E95',policy='Natural',wall='Accepted',**kwargs):
        row=run_case(a.solver,a.output/name,model,policy,wall,descriptor_probe=a.descriptor_probe,**kwargs);results['cases'].append(row);save()
        report=a.output/name/'case.evolution.json'
        return row,load(report) if report.exists() else None
    def record(name,fn):
        try:value=fn();results['checks'].append(dict(name=name,status='pass',result=value))
        except Exception as e:results['checks'].append(dict(name=name,status='fail',failure=repr(e)))
        save()
    matrix=load(a.matrix)
    passed=[r for r in matrix['cases'] if r['status']=='pass' and r['accepted_steps']>=16]
    for wall in ('Accepted','EP'):
        rows=[r for r in passed if r['model']=='E94E95' and r['kappa_policy']=='Natural' and r['wall']==wall]
        if not rows:
            results['checks'].append(dict(name='representative-'+wall,status='not_attempted',reason='16-step primary prerequisite failed'));continue
        base=rows[0];folder=Path(base['path']);whole=load(folder/'case.evolution.json')
        def restart(wall=wall,whole=whole,folder=folder):
            row,r=control('restart-'+wall,wall=wall,steps=8,restart=checkpoint(folder,8))
            check(row['status']=='pass','restart trajectory rejected');check(same_history(whole,r),'restart history differs')
            return dict(bitwise_equal=True,steps=r['steps'])
        record('compatible-restart-'+wall,restart)
        def rollback(wall=wall,folder=folder):
            row,r=control('rollback-'+wall,wall=wall,steps=3,reject=3)
            check(r and r['steps']==2 and 'deliberate PIMPLE rejection' in r['failure'],'rejection not exercised')
            receipt=header(checkpoint(folder,2),True)
            check(state_hash(receipt['coupled_state'])==state_hash(r['state']),'rollback coupled state differs')
            check([state_hash(v) for v in receipt['history']]==[state_hash(v) for v in r['native_history']],'rollback native history differs')
            check(receipt['history_times']==r['native_history_times'],'rollback history times differ')
            return dict(coupled_and_native_history_bitwise_equal=True)
        record('complete-rollback-'+wall,rollback)
        def cache(wall=wall,whole=whole):
            row,r=control('cache-'+wall,wall=wall,steps=16,cache='ExactStepV1')
            check(row['status']=='pass','cache trajectory rejected');check(same_history(whole,r),'cache changed numerical history')
            return dict(bitwise_equal=True)
        record('cache-on-off-'+wall,cache)
        def temporal(wall=wall,whole=whole):
            row,r=control('half-dt-'+wall,wall=wall,steps=32,dt=.00005)
            check(row['status']=='pass','dt-refined trajectory rejected')
            return dict(common_time=r['time'],physical_L2_difference=physical_difference(whole,r),drift=r['step_audits'][-1])
        record('time-refinement-'+wall,temporal)
        def spatial(wall=wall,whole=whole):
            row,r=control('spatial-'+wall,wall=wall,steps=16,refine=True)
            check(row['status']=='pass','spatial/quadrature-refined trajectory rejected')
            coarse=whole['raw_reference_norms'];fine=r['raw_reference_norms']
            return dict(coarse_initial_residual=coarse,fine_initial_residual=fine,
                residual_ratios=[y/x if x>1e-12 else None for x,y in zip(coarse,fine)],coarse_drift=whole['step_audits'][-1],fine_drift=r['step_audits'][-1])
        record('spatial-quadrature-refinement-'+wall,spatial)
        def perturb(wall=wall):
            row,r=control('perturbed-'+wall,wall=wall,steps=8,perturb=.005)
            check(row['status']=='pass','mild positive inner-k perturbation rejected')
            return dict(steps=r['steps'],amplitude_over_K0=.005,drift=r['step_audits'][-1])
        record('mild-admissible-perturbation-'+wall,perturb)
    # Signature and execution guards are checked using valid reference cases;
    # failures must identify the requested unsupported combination explicitly.
    source=next((r for r in passed if r['model']=='E94E95' and r['kappa_policy']=='Natural' and r['wall']=='Accepted'),None)
    if source:
        saved=checkpoint(source['path'],8)
        for label,model,policy,wall in [('changed-model','A95','Natural','Accepted'),('changed-calibration','E94E95','Prescribed','Accepted'),('changed-wall','E94E95','Natural','EP')]:
            def mismatch(label=label,model=model,policy=policy,wall=wall):
                row,r=control(label,model,policy,wall,steps=1,restart=saved)
                text=(a.output/label/'solver.log').read_text()
                check(row['exit_code']!=0 and 'checkpoint signature mismatch' in text,'incompatible restart not rejected by signature')
                return dict(rejected=True)
            record(label+'-restart-rejection',mismatch)
    for label,key,value,reason in [('nonlinear','RANSWallTemporalCoupling','ConvergedStage','one-sweep'),('inner-stress','RANSWallInnerStress','SaturatedR2','Linear'),('regularized','RANSWallVariant','Regularized','CLSHighRe')]:
        def guard(label=label,key=key,value=value,reason=reason):
            folder=a.output/('guard-'+label);xml=case_xml(folder,'E94E95','Natural','Accepted',steps=0,descriptor_probe=a.descriptor_probe)
            root=ET.parse(xml).getroot();set_info(root,key,value);write(root,xml)
            with (folder/'solver.log').open('w') as f:r=subprocess.run([str(a.solver),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,**THREADS),timeout=120)
            text=(folder/'solver.log').read_text();check(r.returncode!=0 and reason in text,('guard rejection reason',label,text[-400:]))
            return dict(rejected=True)
        record('unsupported-'+label,guard)
    results['status']='pass' if all(r['status']=='pass' for r in results['checks']) else 'fail';save()
    print(json.dumps(dict(status=results['status'],checks=[(r['name'],r['status']) for r in results['checks']])))
if __name__=='__main__':main()
