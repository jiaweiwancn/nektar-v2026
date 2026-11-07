#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read completed height-case blocks and publish fixed-physical-window metrics."""
import argparse
import gzip
import json
import os
from pathlib import Path
import time
from analyse_cls_epawf_height_comparison import differences,compact_diagnostic,attempts_summary
from continue_cls_epawf_pimple import CampaignLock,atomic,process_identity,load_progress
from run_cls_epawf_picard import load,sha,require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--once',action='store_true')
    p.add_argument('--cases',nargs=2,default=['G','H'])
    p.add_argument('--initial-reports',nargs=2,type=Path)
    p.add_argument('--campaigns',nargs=2,type=Path)
    a=p.parse_args();require(bool(a.initial_reports)==bool(a.campaigns),'both case paths must be explicit')
    a.output.mkdir(parents=True,exist_ok=True);lock=CampaignLock(a.output);states={};initial={};digests={}
    reports=dict(zip(a.cases,a.initial_reports or [a.root/'focused'/c/'whole4/case.evolution.json' for c in a.cases]))
    campaigns=dict(zip(a.cases,a.campaigns or [a.root/('selected-'+c)/'trajectory' for c in a.cases]))
    for c in a.cases:
        initial[c]=load(reports[c])['height_comparison_initial'];states[c]={0:initial[c]}
        prefix=load(reports[c])
        if prefix.get('height_comparison_final'):states[c][prefix['steps']]=prefix['height_comparison_final']
    while True:
        result=dict(status='observing',observer=process_identity(os.getpid()),production_governing_solves=0,cases={},paired={})
        terminal=True
        for c in a.cases:
            root=campaigns[c];progress=load_progress(root/'progress.json')
            terminal=terminal and progress['status'] in ('numerically_rejected','targets_completed')
            for segment in progress['segments']:
                step=int(segment['end_step']);path=Path(segment['checkpoint']).parent/'case.evolution.json.gz'
                if step>160 and step%80 and step not in progress['requested_targets']:continue
                # receipt checkpoint's parent is the native segment directory.
                if not path.exists():continue
                if step in states[c]:continue
                try:
                    with gzip.open(path,'rt') as stream:d=json.load(stream,parse_int=float)
                except (EOFError,OSError,json.JSONDecodeError):
                    # The driver writes gzip before removing the complete
                    # original. Observe only completed immutable payloads.
                    continue
                require(d['steps']==step,'native block identity differs')
                q=d['height_comparison_final'];require(abs(q['time']-d['time'])<1e-13,'diagnostic time differs')
                states[c][step]=q;digests[str(path)]=sha(path)
            prefix=load(reports[c])['pimple_attempts']
            text=(root/'steps.jsonl').read_text() if (root/'steps.jsonl').exists() else ''
            if text and not text.endswith('\n'):text=text[:text.rfind('\n')+1]
            ledger=[json.loads(v) for v in text.splitlines()];attempts=prefix+[v['attempt'] for v in ledger if v['global_step']<=progress['accepted_steps']]
            q0=initial[c];s0=q0['common_samples'];bulk=sum(w*u for w,u in zip(s0['weights'],s0['u']))/sum(s0['weights']);flow=6/bulk
            case=dict(status=progress['status'],accepted_steps=progress['accepted_steps'],accepted_time=progress['accepted_time'],
                completed_targets=progress['completed_targets'],next_target=progress.get('next_target'),checkpoint=progress['checkpoint'],
                checkpoint_sha256=progress['checkpoint_sha256'],failure=progress.get('failure'),initial_bulk_u=bulk,flow_through_time=flow,H_over_ustar=q0['H']/.4,
                observations={},fixed_time_windows=[],work=attempts_summary(attempts,progress['accepted_time']))
            for n,q in sorted(states[c].items()):
                if not n:continue
                case['observations'][str(n)]=dict(time=q['time'],drift_from_actual_h0=differences(q,q0),
                    error_to_own_analytic=differences(q,q,True),diagnostic=compact_diagnostic(q),
                    elapsed_initial_flow_times=q['time']/flow,elapsed_wall_times=q['time']/(q0['H']/.4))
                if n-32 in states[c]:
                    old=states[c][n-32];require(abs(q['time']-old['time']-.02)<1e-12,'window physical duration differs')
                    case['fixed_time_windows'].append(dict(first_step=n-32,last_step=n,start=old['time'],end=q['time'],
                        physical_L2_change=differences(q,old),wall_stress_change=q['mean_wall_stress']-old['mean_wall_stress']))
                if n>160 and n-80 in states[c]:
                    old=states[c][n-80];require(abs(q['time']-old['time']-.05)<1e-11,'long window physical duration differs')
                    case['fixed_time_windows'].append(dict(first_step=n-80,last_step=n,start=old['time'],end=q['time'],
                        physical_L2_change=differences(q,old),wall_stress_change=q['mean_wall_stress']-old['mean_wall_stress']))
            result['cases'][c]=case
        left,right=a.cases
        for n in sorted(states[left].keys() & states[right].keys()):result['paired'][str(n)]=differences(states[left][n],states[right][n])
        result['immutable_report_sha256']=digests;result['status']='terminal_observed' if terminal else 'observing'
        atomic(a.output/'summary.json',result)
        if terminal or a.once:return
        time.sleep(15)


if __name__=='__main__':main()
