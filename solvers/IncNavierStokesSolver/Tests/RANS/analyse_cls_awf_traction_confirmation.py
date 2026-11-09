#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resolve the specific early traction order with one shorter dt/4 control."""
import argparse
import math
from pathlib import Path
from analyse_cls_accepted_awf import report
from diagnose_cls_accepted_covariance import checkpoint
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import load,require,sha,state_hash


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--primary-analysis',type=Path,required=True)
    p.add_argument('--confirmation',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    primary=load(a.primary_analysis/'summary.json');policy=primary['policy'];progress=load(a.confirmation/'extra-fine/progress.json')
    require(progress['status']=='targets_completed','targeted confirmation is not complete')
    out=dict(status='recorded',purpose='traction-only early nonasymptotic discrimination; no dt/4 trajectory through .05 or T=5',
        primary_summary_sha256=sha(a.primary_analysis/'summary.json'),confirmation_progress_sha256=sha(a.confirmation/'extra-fine/progress.json'),
        desired_long_campaign_dt=.000625,triplet_dt=[.000625,.0003125,.00015625],rows={},accepted_steps=progress['accepted_steps'],accepted_time=progress['accepted_time'])
    for t in (.01,.02):
        key=str(t);entry=progress['observations'][key];cp=Path(entry['checkpoint']);c=checkpoint(cp);native=report(Path(entry['report']))
        require(sha(cp/'State.bin')==entry['checkpoint_sha256'] and state_hash(c['coupled'])==state_hash(native['state']) and c['time']==native['time'],'confirmation checkpoint/report identity differs')
        require([state_hash(v) for v in c['history']]==[state_hash(v) for v in native['native_history']] and state_hash(c['history_times'])==state_hash(native['native_history_times']),'confirmation GLM history differs')
        require(native['audit_revision']=='NativeElementSideV2' and not native['failure'],'confirmation hard audit/version failed')
        tau=[primary['observations'][label][key]['physical']['mean_wall_stress'] for label in ('middle','fine')]+[native['accepted_awf_final']['mean_wall_stress']]
        dc,df=abs(tau[0]-tau[1]),abs(tau[1]-tau[2]);order=math.log2(dc/df);coherent=policy['coherent_order_interval'][0]<=order<=policy['coherent_order_interval'][1]
        estimate=2**order*dc/(2**order-1)/.16 if coherent else None
        out['rows'][key]=dict(traction=tau,Dc=dc,Df=df,p=order,coherent=coherent,
            relative_error_estimate_for_dt_000625=estimate,within_target=estimate is not None and estimate<=policy['middle_step_relative_mean_traction_error_maximum'],
            checkpoint_sha256=sha(cp/'State.bin'),coupled_state_sha256=state_hash(c['coupled']),report_sha256=sha(Path(entry['report'])),
            reused_native_alpha_zero_diagnostic=True)
    fields=all(all(v['coherent_first_order'] and v['within_target'] for v in primary['triplets'][str(t)].values()) for t in policy['coherent_late_observations'])
    traction=all(v['coherent'] and v['within_target'] for v in out['rows'].values()) and all(primary['engineering_triplets'][str(t)]['traction_within_target'] for t in (.04,.05))
    out['combined_classification']='finite_horizon_fields_and_traction_consistent_first_order_within_targets' if fields and traction else 'traction_or_field_assessment_remains_inconclusive'
    atomic(a.output,out);print(a.output)

if __name__=='__main__':main()
