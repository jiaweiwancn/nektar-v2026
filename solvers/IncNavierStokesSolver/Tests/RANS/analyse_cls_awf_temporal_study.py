#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict immutable checkpoint triplet, common physical norms and one gauge."""
import argparse
from datetime import datetime,timezone
import json
import math
import os
from pathlib import Path
import shutil
import xml.etree.ElementTree as ET
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import load,require,sha,state_hash,set_info,write,THREADS
from analyse_cls_accepted_awf import differences,normalized,compact
from build_cls_accepted_awf_analysis_suite import compact_diagnostic


def middle(campaign,times):
    p=load(campaign/'progress.json');out={}
    ledger={}
    for line in (campaign/'steps.jsonl').read_text().splitlines():
        try:r=json.loads(line)
        except json.JSONDecodeError:break
        ledger[int(r['global_step'])]=r
    for t in times:
        n=round(t/.000625);segments=[s for s in p['segments'] if s['end_step']==n and not s['failure']]
        require(segments,'middle immutable native checkpoint unavailable: '+str(n))
        s=segments[0];cp=Path(s['checkpoint']);require(sha(cp/'State.bin')==ledger[n]['checkpoint_sha256'],'middle checkpoint ledger differs')
        out[str(t)]=dict(checkpoint=str(cp),xml=str(cp.parent/'case.xml'),attempt=ledger[n]['attempt'])
    return out


def observation(probe,row,folder):
    folder.mkdir(parents=True,exist_ok=False);cp=folder/'checkpoint';cp.mkdir();source=Path(row['checkpoint']);shutil.copy2(source/'State.bin',cp/'State.bin')
    q=checkpoint(cp);require(sha(cp/'State.bin')==sha(source/'State.bin'),'immutable copy differs')
    tree=ET.parse(row['xml']).getroot();set_info(tree,'RANSWallEvaluationReuse','ExactStepV1');write(tree,folder/'case.xml')
    packet=dict(common_packet(cp),operation='selected_original_operator',state=q['coupled'],history=q['coupled'],time=q['time'],alpha=0.,
        inner_stress_model='Linear',reporting_breakpoints=[.0025,.005,.02,.04],witnesses=[[-3,0,3.6880433085817326],[1.6,0,2.]])
    result,receipt=run_probe(probe,folder/'case.xml',packet,folder/'readback');d=result['diagnostics'];physical=d['physical']
    require(result['governing_solves']==result['physical_commits']==0 and not result['native_history_advanced'] and result['initialization_work']==[0,0],'readback advanced/lifted state')
    compacted=dict(checkpoint_sha256=sha(cp/'State.bin'),original_checkpoint=str(source),state_sha256=state_hash(q['coupled']),
        native_history_sha256=[state_hash(v) for v in q['history']],native_history_times_sha256=state_hash(q['history_times']),receipt=receipt,
        physical=compact_diagnostic(compact(physical)),complete_audit=result['complete_audits']['selected_state'],operators={k:v for k,v in d['operators'].items() if k!='blocks'},
        spatial_blocks=[{k:v for k,v in b.items() if k not in ('free_rows','selected_load','same_state_Original_load','difference_load')} for b in d['operators']['blocks']],
        covariance=d['selected_state']['minimum_outer'],common_covariance_witnesses=d['selected_state']['one_sided_witnesses'])
    return physical,compacted


def complete_engineering(out):
    policy=out['policy'];out['engineering_triplets']={}
    for t in policy['observation_times']:
        key=str(t)
        if not all(key in out['observations'].get(label,{}) for label in ('coarse','middle','fine')):continue
        rows=[out['observations'][label][key] for label in ('coarse','middle','fine')]
        traction=[r['physical']['mean_wall_stress'] for r in rows]
        dc,df=abs(traction[0]-traction[1]),abs(traction[1]-traction[2]);order=math.log2(dc/df) if min(dc,df)>.16e-11 else None
        coherent=order is not None and policy['coherent_order_interval'][0]<=order<=policy['coherent_order_interval'][1]
        relative=2**order*df/(2**order-1)/.16 if coherent else None
        out['engineering_triplets'][key]=dict(traction_values=traction,traction_Dc=dc,traction_Df=df,traction_p=order,
            relative_middle_traction_Richardson_error=relative,traction_target=policy['middle_step_relative_mean_traction_error_maximum'],
            traction_within_target=relative is not None and relative<=policy['middle_step_relative_mean_traction_error_maximum'],
            matching_defects=[r['physical']['matching_trace_defect'] for r in rows],
            strong_continuity=[r['complete_audit']['strong_continuity'] for r in rows],
            interface=[r['complete_audit']['interface'] for r in rows],
            covariance_minima=[r['complete_audit']['accepted_awf_covariance'] for r in rows],
            actual_momentum_splitting_load_defect=[r['frozen_step_diagnostics']['frozen_momentum_defect'] for r in rows],
            splitting_identity_error=[r['frozen_step_diagnostics']['momentum_splitting_identity_defect'] for r in rows])
    if out.get('classification')=='consistent_with_expected_order_and_within_declared_budget':
        if not all(out['engineering_triplets'].get(str(t),{}).get('traction_within_target',False) for t in policy['coherent_late_observations']):
            out['classification']='fields_within_budget_traction_assessment_inconclusive_or_outside_budget'
    return out


def main():
    os.environ.update(THREADS)
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',required=True,type=Path);p.add_argument('--middle-campaign',required=True,type=Path)
    p.add_argument('--probe',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    policy=load(a.study/'declared-policy.json');times=policy['observation_times'];sources={'middle':middle(a.middle_campaign,times)}
    for label in ('coarse','fine'):
        progress=load(a.study/label/'progress.json');sources[label]=progress['observations']
    out=dict(status='running',started_UTC=datetime.now(timezone.utc).isoformat(),policy=policy,policy_sha256=sha(a.study/'declared-policy.json'),
        common_pressure_gauge='one physical-volume mean shift across both pressure regions; no regional offset removal',observations={},triplets={})
    fields={}
    for label,rows in sources.items():
        fields[label]={};out['observations'][label]={}
        for t in times:
            if str(t) not in rows:continue
            physical,summary=observation(a.probe,rows[str(t)],a.output/label/str(t));fields[label][t]=physical
            attempt=rows[str(t)].get('attempt')
            if attempt is None:
                from analyse_cls_accepted_awf import report
                native=report(Path(rows[str(t)]['report']));attempt=native['pimple_attempts'][-1]
            summary['frozen_step_diagnostics']={k:attempt.get(k) for k in ('frozen_momentum_defect','momentum_splitting_identity_defect','refreshed_trace_defect',
                'corrected_continuity_gauge','bound_activations','mass_inverse_applications','k_added_frozen_residual','k_added_storage_residual')}
            summary['frozen_step_diagnostics']['momentum_defect_metric']='maximum absolute free-row frozen momentum load, distinct from the splitting identity error'
            out['observations'][label][str(t)]=summary;atomic(a.output/'summary.json',out)
    for t in times:
        if not all(t in fields[label] for label in ('coarse','middle','fine')):continue
        dc=differences(fields['coarse'][t],fields['middle'][t]);df=differences(fields['middle'][t],fields['fine'][t]);nc,nf=normalized(dc),normalized(df)
        values={}
        for field in ('u','v','k','epsilon','pressure'):
            tiny=min(nc[field],nf[field])<=policy['near_zero_normalized_difference']
            order=None if tiny else math.log2(dc[field]/df[field]);coherent=order is not None and policy['coherent_order_interval'][0]<=order<=policy['coherent_order_interval'][1]
            fine_error=df[field]/(2**order-1) if coherent else None;mid_error=2**order*fine_error if coherent else None
            normalized_middle=2**order*nf[field]/(2**order-1) if coherent else None
            target=policy['middle_step_Richardson_error_maximum'][field]
            values[field]=dict(Dc=dc[field],Df=df[field],normalized_Dc=nc[field],normalized_Df=nf[field],p=order,
                coherent_first_order=coherent,near_zero=tiny,fine_Richardson_error=fine_error,middle_Richardson_error=mid_error,
                normalized_middle_Richardson_error=normalized_middle,target=target,within_target=normalized_middle is not None and normalized_middle<=target)
        out['triplets'][str(t)]=values
    required=[out['triplets'].get(str(t),{}) for t in (.02,.04,.05)]
    if len(out['triplets'])<len(times):classification='rejected_or_incomplete_before_common_horizon'
    elif all(all(v.get('coherent_first_order') and v.get('within_target') for v in r.values()) and len(r)==5 for r in required):classification='consistent_with_expected_order_and_within_declared_budget'
    elif all(all(v.get('p') is not None and v['p']>0 for v in r.values()) and len(r)==5 for r in required):classification='ordered_but_accuracy_or_asymptotic_order_not_established'
    else:classification='inconclusive'
    out.update(status='recorded',classification=classification,completed_UTC=datetime.now(timezone.utc).isoformat());complete_engineering(out);atomic(a.output/'summary.json',out)
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        for field in ('u','v','k','epsilon','pressure'):
            ts=sorted(float(t) for t in out['triplets']);axes[0].plot(ts,[out['triplets'][str(t)][field]['p'] for t in ts],'o-',label=field)
            axes[1].semilogy(ts,[out['triplets'][str(t)][field]['normalized_Df'] for t in ts],'o-',label=field)
        axes[0].axhspan(.8,1.2,color='grey',alpha=.2);axes[0].set(xlabel='physical time',ylabel='observed temporal order')
        axes[1].set(xlabel='physical time',ylabel='normalized middle–fine difference');axes[0].legend();fig.savefig(a.output/'temporal-triplet.png',dpi=160);plt.close(fig)
    except ImportError:pass
    print(a.output/'summary.json')

if __name__=='__main__':main()
