#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact one-step revalidation of a diagnostic-only WP2 runtime successor."""
import argparse
from pathlib import Path
from wp2_native import Run
from wp2_qualification import verify_report
from run_cls_epawf_picard import require,sha,state_hash
from compare_cls_epawf_pimple import same_history
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser()
    for n in ['root','old','new','h1','preflight','parent_progress','output']:p.add_argument('--'+n.replace('_','-'),type=Path,required=True)
    p.add_argument('--model',choices=['CLS','standard'],required=True)
    a=p.parse_args();prefix=load(a.parent_progress);parent=Path(prefix['observations'][-1]['checkpoint']);before=sha(parent/'State.bin')
    a.output.mkdir(exist_ok=False,parents=True)
    reports=[];runs=[]
    K=.017598891430879275 if a.model=='CLS' else .01838575534504143
    for label,solver in [('old',a.old),('new',a.new)]:
        x=Run(a.root,solver,a.h1,a.preflight,a.output/label,a.model,4,16,'RadialEnergyV1')
        d=x.native('one',1,.02,parent);checks=verify_report(d,.07426793792419734,K)
        x.result.update(status='completed_bounded_trajectory',mode='exact-runtime-equivalence',dt=.02,accepted_steps=d['steps'],actual_time=d['time'],
                        observations=[dict(time=d['time'],accepted=d['steps'],checkpoint=d['checkpoint'],checkpoint_sha256=sha(Path(d['checkpoint'])/'State.bin'),
                                           report=x.result['results']['one']['report'],physical_budgets_passed=True,final_physical_checks=checks[-1])])
        x.save();reports.append(d);runs.append(x)
    require(same_history(*reports),'diagnostic runtime changed coupled fields or native history')
    for name in ['history_sha256','incoming_sha256','post_flow_sha256','pre_k_sha256','raw_k_sha256','corrected_k_sha256']:
        require(reports[0]['pimple_attempts'][0][name]==reports[1]['pimple_attempts'][0][name],'governing stage differs: '+name)
    for name in ['fixed_reference','fixed_defect']:
        require(state_hash(reports[0][name])==state_hash(reports[1][name]),name+' changed')
    for left,right in zip(reports[0]['state_audits'],reports[1]['state_audits']):
        require(left.keys()==right.keys() and all(state_hash(left[k])==state_hash(right[k]) for k in left),'native audit values changed')
    table=reports[1]['wp2_per_state_monitoring'][0]['physical_checks'];ids=[int(x) for x in table['physical_essential_rows']];values=table['physical_essential_values']
    prefix_errors=[]
    for dt,n in [('020',2),('010',4)]:
        d=load(a.root/('trial-dt'+dt+'-'+a.model)/('step'+str(n))/'case.evolution.json')
        error=max(abs(s[i]-v) for s in d['accepted_states'] for i,v in zip(ids,values));require(error==0,'earlier accepted prescribed coefficients changed');prefix_errors.append(error)
    require(sha(parent/'State.bin')==before,'runtime comparison rewrote parent checkpoint')
    atomic(a.output/'result.json',dict(passed=True,model=a.model,bitwise_governing_stages_state_history_audits=True,
        reference_defect_preserved=True,parent_unchanged=True,prefix_physical_essential_absolute_errors=prefix_errors,
        new_successor_progress=str(a.output/'new/progress.json'),new_progress_sha256=sha(a.output/'new/progress.json')))


if __name__=='__main__':main()
