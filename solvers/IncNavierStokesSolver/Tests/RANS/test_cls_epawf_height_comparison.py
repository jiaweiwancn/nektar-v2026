#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused native tests of the declared height/sampler family, not old h0 waivers."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import gzip
import json
import os
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_height_comparison import case_xml,WORK
from run_cls_epawf_picard import load,require,run,write_json,sha,state_hash,THREADS,set_info,set_parameter,write
from run_cls_epawf_pimple import prepare,check_report,identity
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import freeze,header,atomic
from test_cls_epawf_k_bound import replay


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--screening',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--upstream-baseline',type=Path,required=True)
    a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    os.environ.update(THREADS)
    solver=a.solver.resolve();probe,_=freeze(a.probe.resolve(),a.output/'probe-runtime')
    os.environ['LD_LIBRARY_PATH']=str(a.output/'probe-runtime/lib')
    summary=dict(status='running',solver=identity(solver),probe=identity(probe),checks={})
    atomic(a.output/'progress.json',summary)
    def native(name,case,steps,restart=None,change=None,hard_reject=0):
        folder=a.output/name;xml=case_xml(case,folder,steps,restart=restart,hard_reject=hard_reject)
        if change:
            root=ET.parse(xml).getroot();change(root);write(root,xml)
        receipt=run([str(solver),str(xml)],folder,'solver.log',1800);write_json(folder/'receipt.json',receipt)
        path=folder/'case.evolution.json'
        return (load(path) if path.exists() else None),receipt,xml
    def probe_at(name,xml,state,operation='height_comparison',**extra):
        folder=a.output/name;folder.mkdir(parents=True)
        write_json(folder/'input.json',dict(operation=operation,state=state,history=state,time=0.,alpha=0.,require_history_h0=False,inner_stress_model='SaturatedR2',**extra))
        receipt=run([str(probe),str(xml), 'input.json','result.json'],folder,'probe.log',1800)
        write_json(folder/'receipt.json',receipt);require(receipt['exit_code']==0,'probe failed '+name)
        return load(folder/'result.json')
    def initial_checks():
        result={}
        for case in 'EFGH':
            initial,rc,xml=native('initial-'+case,case,0)
            require(rc['exit_code']==0 and not initial['failure'],'native initialization '+case)
            check_report(initial,True)
            require(initial['initial_k_bound_check']['activations']==0,'native h0 was repaired')
            check=probe_at('sampler-'+case,xml,initial['state'],'sampling_checks')
            result[case]=dict(initial_state_sha256=state_hash(initial['state']),sampling=check,
                              initialization=initial['height_comparison_final'])
            atomic(a.output/'initial-checks.json',result)
        # Canonical placeholder has no effect on the average functional, profiles,
        # boundary data or method identity. Use a valid but different location.
        base=load(a.output/'initial-G/case.evolution.json')
        changed,rc,_=native('unused-placeholder','G',0,restart=Path(base['checkpoint']),
                           change=lambda root:set_parameter(root,'RANSWallSampleHeight',.1))
        require(rc['exit_code']==0,'unused average placeholder changed restart signature')
        result['unused_placeholder_exact_history']=same_history(base,changed)
        # H, mesh, tail and selected functional are genuine incompatibilities.
        bad=[]
        for name,case,checkpoint in [('mesh','H',base['checkpoint']),('height','E',base['checkpoint']),
                 ('sampler','C',str(a.screening/'B/initial/checkpoint')),('tail','B',str(a.screening/'A/initial/checkpoint'))]:
            report,rc,xml=native('incompatible-'+name,case,0,restart=Path(checkpoint))
            text=(xml.parent/'solver.log').read_text(errors='replace')
            require(rc['exit_code']!=0 and ('signature' in text.lower() or 'checkpoint' in text.lower()),'incompatible '+name+' accepted')
            bad.append(dict(change=name,exit_code=rc['exit_code'],log_sha256=rc['log_sha256']))
        result['incompatible_cases']=bad
        # Between-node ground separation exercises native CG and essential data
        # with the selected mean functional present throughout the full audit.
        fixture=replay(probe,a.output/'initial-G/case.xml',base,a.output/'mean-between-nodes',fixture='between_nodes')
        require(fixture['k_bound']['success'] and not fixture['complete_audit_failure'],'mean ground fixture audit')
        result['mean_between_nodes']={k:fixture[k] for k in ('k_bound','non_k_change','essential_k_change','shared_cg_jump','quadrature_correction_mass_norm')}
        return result
    def restart_checks():
        result={}
        # These native four-step prefixes also seed the selected smaller-H pair.
        for case in 'GH':
            whole,rc,_=native(case+'/whole4',case,4)
            require(rc['exit_code']==0 and not whole['failure'],'short smaller-H trajectory failed '+case)
            check_report(whole,True)
            atomic(a.output/(case+'-prefix.json'),dict(checkpoint=whole['checkpoint'],report=str(a.output/case/'whole4/case.evolution.json'),steps=whole['steps'],time=whole['time']))
            if case=='H':
                half,rc,_=native('restart/half','H',2)
                require(rc['exit_code']==0 and not half['failure'],'restart half')
                final,rc,_=native('restart/final','H',2,restart=Path(half['checkpoint']))
                require(rc['exit_code']==0 and not final['failure'],'restart final')
                result['restart_exact_state_history_time']=same_history(whole,final)
                rejected,_,_=native('rollback','H',1,restart=Path(half['checkpoint']),hard_reject=3)
                require(rejected and 'deliberate' in rejected['failure'].lower(),'forced hard rejection missing')
                result['rejection_exact_state_history_time']=same_history(half,rejected)
            result[case]=dict(steps=whole['steps'],time=whole['time'],state_sha256=state_hash(whole['state']))
        return result
    def baseline_checks():
        folder=a.output/'upstream-active6'
        xml=prepare(folder,.000625,6,True,checkpoint=True,k_bound='MassProjection',floor_relative=1e-8,snapshots=True,
                    inventory='ConstitutiveV2',wall_policy='Upstream',bound_work=WORK)
        rc=run([str(solver),str(xml)],folder,'solver.log',1800);write_json(folder/'receipt.json',rc)
        current=load(folder/'case.evolution.json');baseline=load(a.upstream_baseline)
        require(rc['exit_code']==0 and any(v['k_bound']['activations'] for v in current['pimple_attempts']),'active baseline not exercised')
        return dict(point_upstream_active6_exact=same_history(current,baseline))
    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs={pool.submit(fn):name for name,fn in [('initial_and_functional',initial_checks),('restart',restart_checks),('baseline',baseline_checks)]}
        for future in as_completed(jobs):
            name=jobs[future]
            try:summary['checks'][name]=dict(passed=True,result=future.result())
            except BaseException as error:summary['checks'][name]=dict(passed=False,error=repr(error))
            atomic(a.output/'progress.json',summary)
    summary['status']='passed' if all(c['passed'] for c in summary['checks'].values()) else 'failed'
    atomic(a.output/'progress.json',summary);require(summary['status']=='passed','focused check failed; inspect retained artifacts')


if __name__=='__main__':main()
