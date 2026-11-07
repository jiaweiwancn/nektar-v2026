#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused native V3 tests, all fresh physical fixtures have z_m=H."""
import argparse
import os
import sys
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_height_comparison import case_xml
from run_cls_epawf_picard import load,write_json,require,run,THREADS,set_info,set_parameter,write,state_hash
from run_cls_epawf_pimple import check_report
from compare_cls_epawf_pimple import same_history
from continue_cls_epawf_pimple import freeze,header

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();require(sys.platform=='linux','Linux authoritative')
    a.output.mkdir(parents=True,exist_ok=False);os.environ.update(THREADS)
    solver,_=freeze(a.solver.resolve(),a.output/'runtime');probe,_=freeze(a.probe.resolve(),a.output/'probe-runtime')
    os.environ['LD_LIBRARY_PATH']=str(a.output/'runtime/lib');results={}
    def native(name,case,steps,inventory='ConstitutiveV3',restart=None,**kw):
        folder=a.output/name;xml=case_xml(case,folder,steps,inventory=inventory,restart=restart,**kw)
        rc=run([str(solver),str(xml)],folder,'solver.log',3600);write_json(folder/'receipt.json',rc)
        report=load(folder/'case.evolution.json');check_report(report,True)
        return report,xml,rc
    def query(name,xml,report,operation,**kw):
        folder=a.output/name;folder.mkdir();write_json(folder/'input.json',dict(operation=operation,state=report['state'],
            history=report['state'],time=report['time'],alpha=0.,require_history_h0=False,inner_stress_model='SaturatedR2',**kw))
        rc=run([str(probe),str(xml),'input.json','result.json'],folder,'probe.log',1800)
        write_json(folder/'receipt.json',rc);require(rc['exit_code']==0,'probe failed '+name);return load(folder/'result.json')
    # A real matched z_m=H control; V2's archived state is not relabelled V3.
    old,oldxml,rc=native('uniform-v2','G',4,inventory='ConstitutiveV2')
    require(not old['failure'] and rc['exit_code']==0,'matched V2 prefix failed')
    new,newxml,rc=native('uniform-v3','G',4)
    require(not new['failure'] and rc['exit_code']==0,'matched V3 prefix failed')
    results['uniform_exact_state_history_time']=same_history(old,new)
    results['uniform_inactive_bound']=all(v['bound_activations']==0 for v in new['pimple_attempts'])
    require(results['uniform_inactive_bound'],'matched uniform correction unexpectedly active')
    results['coverage_tripwire']=query('coverage',newxml,new,'coverage')
    results['mean_and_fixed_band_checks']=query('sampling',newxml,new,'sampling_checks')
    # Nonuniform synthetic polynomial is checked by sampling_checks, including
    # exact global mean rows, negative local k, and protected coefficients.
    whole,xml,rc=native('upstream-whole','U05-one',4)
    require(not whole['failure'] and rc['exit_code']==0,'Upstream prefix failed')
    results['upstream_mean_checks']=query('upstream-sampling',xml,whole,'sampling_checks')
    results['upstream_branches']=query('upstream-branches',xml,whole,'branch_coverage')
    half,_,rc=native('upstream-half','U05-one',2);require(not half['failure'] and rc['exit_code']==0,'half prefix failed')
    final,_,rc=native('upstream-restart','U05-one',2,restart=Path(half['checkpoint']))
    require(not final['failure'] and rc['exit_code']==0,'compatible restart failed')
    results['restart_exact_state_history_time']=same_history(whole,final)
    rejected,_,_=native('upstream-rollback','U05-one',1,restart=Path(half['checkpoint']),hard_reject=3)
    require('deliberate' in rejected['failure'].lower(),'forced hard rejection missing')
    results['rollback_exact_state_history_time']=same_history(half,rejected)
    # These are reader-negative tests: no incompatible physical step is run.
    incompatible={}
    for name in ('inventory','wall_policy','mesh','sampler'):
        folder=a.output/('incompatible-'+name)
        case='U05-graded' if name=='mesh' else 'U05-one'
        x=case_xml(case,folder,1,inventory='ConstitutiveV3',restart=Path(half['checkpoint']))
        tree=ET.parse(x).getroot()
        if name=='inventory':set_info(tree,'RANSWallPimpleKInventory','ConstitutiveV2')
        if name=='wall_policy':set_info(tree,'RANSWallTransition','UniformEPFlat')
        if name=='sampler':set_info(tree,'RANSWallKSampling','Point')
        write(tree,x);before=header(Path(half['checkpoint']))
        rc=run([str(solver),str(x)],folder,'solver.log',1800)
        require(rc['exit_code']!=0,'incompatible reader accepted '+name)
        require(header(Path(half['checkpoint']))==before,'reader changed input checkpoint')
        incompatible[name]=dict(exit_code=rc['exit_code'],input_checkpoint_unchanged=True,
            message=(folder/'solver.log').read_text()[-2000:])
    results['strict_reader_rejections']=incompatible
    results['one_sweep_attempts']=whole['pimple_attempts'];results['h0_sha256']=state_hash(whole['accepted_states'][0])
    write_json(a.output/'summary.json',results);print(a.output/'summary.json')

if __name__=='__main__':main()
