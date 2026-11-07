#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Active V3 rollback and Point/disabled/intrinsic branch controls; z_m=H."""
import argparse
import os
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_height_comparison import case_xml
from run_cls_epawf_picard import load,run,require,write_json,write,set_parameter,set_info,THREADS,state_hash
from run_cls_epawf_pimple import check_report
from continue_cls_epawf_pimple import header
from compare_cls_epawf_pimple import same_history

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--baseline-solver',type=Path)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    os.environ.update(THREADS);summary={}
    def native(name,case,steps,restart=None,parameters=None,**kw):
        folder=a.output/name;xml=case_xml(case,folder,steps,inventory='ConstitutiveV3',restart=restart,**kw)
        if parameters:
            tree=ET.parse(xml).getroot()
            for k,v in parameters.items():set_parameter(tree,k,v)
            write(tree,xml)
        rc=run([str(a.solver),str(xml)],folder,'solver.log',1800);write_json(folder/'receipt.json',rc)
        report=load(folder/'case.evolution.json');check_report(report,True);return report,xml,rc
    def query(name,xml,report,operation):
        folder=a.output/name;folder.mkdir();write_json(folder/'input.json',dict(operation=operation,state=report['state'],history=report['state'],
            time=report['time'],alpha=0.,require_history_h0=False,inner_stress_model='SaturatedR2'))
        rc=run([str(a.probe),str(xml),'input.json','result.json'],folder,'probe.log',1800);write_json(folder/'receipt.json',rc)
        require(rc['exit_code']==0,'diagnostic failed '+name);return load(folder/'result.json')
    cps=[p.parent for p in (a.root/'screening/U05-graded/trajectory').glob('checkpoint*/State.bin') if header(p.parent)['steps']==8]
    require(len(cps)==1,'need the immutable accepted graded step-8 checkpoint')
    cp=a.output/'checkpoint8';shutil.copytree(cps[0],cp)
    base,_,rc=native('restart-zero','U05-graded',0,restart=cp);require(rc['exit_code']==0 and not base['failure'],'restart startup rejected')
    require(base['initial_k_bound_check']['activations']==0,'restart startup changed retained k')
    rejected,_,_=native('active-reject','U05-graded',1,restart=cp,hard_reject=9)
    require('deliberate' in rejected['failure'].lower(),'post-bound forced rejection not reached')
    row=rejected['pimple_attempts'][-1];require(row['bound_activations']==1,'rollback did not exercise active bounding')
    summary['active_rollback_exact_state_history_time']=same_history(base,rejected)
    packet=load(next((a.output/'active-reject/pimple-states').glob('*.json')))
    lo,hi=map(int,rejected['offsets'][2:4]);raw=packet['raw'];corrected=packet['corrected']
    require(state_hash(raw[:lo]+raw[hi:])==state_hash(corrected[:lo]+corrected[hi:]),'bound modified non-k coefficients')
    free=set(map(int,packet['transport_rows']))
    require(state_hash([raw[i] for i in range(lo,hi) if i not in free])==
        state_hash([corrected[i] for i in range(lo,hi) if i not in free]),'bound modified essential k')
    summary['protected_coefficients_exact']=True;summary['active_attempt']=row;summary['checkpoint']=header(cp)
    # Same nonuniform evolved fields: inventory changes constraints, not the
    # selected native source/weak loads or dependent epsilon target.
    selected=base
    original_xml=a.output/'restart-zero/case.xml'
    v3=query('same-state-v3',original_xml,selected,'height_comparison')
    folder=a.output/'same-state-xml';folder.mkdir();tree=ET.parse(original_xml).getroot()
    set_info(tree,'RANSWallPimpleKInventory','ConstitutiveV2');write(tree,folder/'case.xml')
    v2=query('same-state-v2',folder/'case.xml',selected,'height_comparison')
    require(all(state_hash(v)==state_hash(v2['weak_k_component_loads'][k]) for k,v in v3['weak_k_component_loads'].items()),
        'inventory policy changed native source/weak loads')
    require(v3['matching_trace_defect']==v2['matching_trace_defect'],'inventory changed dependent matching')
    summary['same_state_all_native_weak_loads_exact']=True
    point,xml,rc=native('point-initial','U20-point',0);require(rc['exit_code']==0 and not point['failure'],'fresh Point bridge initialization failed')
    summary['point_sampling']=query('point-sampling',xml,point,'sampling_checks')
    branch=query('point-branches',xml,point,'branch_coverage');summary['point_branches']=branch
    reasons={v['anchor'][0]:v['prospective_reason'] for v in branch['state']['anchors']}
    require(reasons[-3.5]=='permitted_active' and reasons[-1.]=='spatially_excluded','Point prospective support is wrong')
    disabled,xml,rc=native('disabled-initial','U05-one',0,parameters={'RANSWallChi':0})
    require(rc['exit_code']==0 and not disabled['failure'],'disabled physical fixture initialization failed')
    d=query('disabled-branches',xml,disabled,'branch_coverage');summary['disabled']=d
    require(all(v['prospective_reason']=='disabled' for v in d['state']['anchors']),'disabled correction forced an active profile')
    gated,xml,rc=native('reynolds-gated-initial','U05-one',0,
        parameters={'RANSWallSensorReynoldsOff':20000,'RANSWallSensorReynoldsOn':23000})
    require(rc['exit_code']==0 and not gated['failure'],'gate test initialization failed')
    low=dict(gated);low['state']=list(gated['state']);lo,hi=map(int,gated['offsets'][2:4])
    scale=.25/gated['height_comparison_final']['reference']['K0']
    for i in range(lo,hi):low['state'][i]*=scale
    d=query('intrinsic-gate',xml,low,'branch_coverage');summary['intrinsic_gate']=d
    require(d['state']['anchors'][0]['prospective_reason']=='intrinsic_gate_inactive','inactive Reynolds gate constructed Q+')
    if a.baseline_solver:
        current,_,rc=native('diagnostic-map-new','U05-graded',2,restart=cp)
        require(rc['exit_code']==0 and not current['failure'],'new reporting-runtime map check failed')
        folder=a.output/'diagnostic-map-old';xml=case_xml('U05-graded',folder,2,inventory='ConstitutiveV3',restart=cp)
        rc=run([str(a.baseline_solver),str(xml)],folder,'solver.log',1800);write_json(folder/'receipt.json',rc)
        previous=load(folder/'case.evolution.json');require(rc['exit_code']==0 and not previous['failure'],'old frozen-runtime map check failed')
        summary['reporting_change_exact_active_state_history_time']=same_history(current,previous)
    write_json(a.output/'summary.json',summary);print(a.output/'summary.json')

if __name__=='__main__':main()
