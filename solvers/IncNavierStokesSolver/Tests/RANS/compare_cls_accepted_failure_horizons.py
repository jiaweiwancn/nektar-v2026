#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable common-time checkpoint comparisons for the two faithful controls.

No checkpoint migration, field transfer, governing solve or history advance.
The strict native reader verifies each case with its own session and histories.
"""
import argparse
from datetime import datetime,timezone
from pathlib import Path
import shutil

from analyse_cls_accepted_awf import compact,differences,normalized,reference_errors
from build_cls_accepted_awf_analysis_suite import compact_diagnostic
from continue_cls_epawf_pimple import atomic,header
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from run_cls_epawf_picard import load,require,sha,state_hash


def native_checkpoints(directory):
    rows={}
    for file in directory.glob('checkpoint*/State.bin'):
        item=header(file.parent);n=item['steps']
        if n in rows:require(item['sha256']==rows[n][1]['sha256'],'conflicting native aliases')
        rows[n]=(file.parent,item)
    return rows


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--old-root',type=Path,required=True)
    p.add_argument('--controls-root',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--times',type=float,nargs='+',default=[0.,.02,.04,.06,.08,.098125,.1,.12]);a=p.parse_args()
    a.output.mkdir(exist_ok=False);old=a.old_root;root=a.controls_root
    progress=load(old/'selected-refined-2G/trajectory/progress.json')
    baseline=native_checkpoints(old/'refined-2G/Accepted-AWF-2G/initial')
    # The durable campaign inherited the completed 32-step short trajectory.
    # Its earlier checkpoints are genuine native data outside later segments.
    baseline.update(native_checkpoints(old/'refined-2G/Accepted-AWF-2G/trajectory'))
    for segment in progress['segments']:
        if not segment['failure']:
            baseline.update(native_checkpoints(Path(segment['checkpoint']).parent))
    baseline.update(native_checkpoints(Path(progress['checkpoint']).parent))
    folders={'baseline':None,'half-dt':root/'half-dt/native/Accepted-AWF-2G',
        'pressure-p2':root/'outer-pressure-p2-v2/native/Accepted-AWF-2G'}
    out=dict(schema=1,captured_UTC=datetime.now(timezone.utc).isoformat(),status='started',
        governing_solves=0,physical_commits=0,native_history_advanced=False,
        reporting_breakpoints=[.0025,.005,.02,.04],cases={},same_time_differences={},
        scope='strict alpha=0 native-checkpoint observations; common physical panels; one common pressure shift; own actual h0 distinct from analytical reference')
    physical={}
    for name,base in folders.items():
        available=baseline if base is None else {**native_checkpoints(base/'initial'),**native_checkpoints(base/'trajectory')}
        dt=.0003125 if name=='half-dt' else .000625
        case=out['cases'][name]=dict(dt=dt,observations={},unavailable_times=[]);physical[name]={}
        for t in a.times:
            n=round(t/dt);require(abs(n*dt-t)<1e-12,'comparison time is not an integer native step')
            if n not in available:case['unavailable_times'].append(t);continue
            source,h=available[n];folder=a.output/name/('step-%06d'%n);folder.mkdir(parents=True)
            cp=folder/'checkpoint';cp.mkdir();shutil.copy2(source/'State.bin',cp/'State.bin')
            require(sha(cp/'State.bin')==h['sha256'],'immutable comparison checkpoint differs')
            xml=source.parent/'case.xml';native=checkpoint(cp)
            packet=dict(common_packet(cp),operation='selected_original_operator',state=native['coupled'],
                history=native['coupled'],time=native['time'],alpha=0.,inner_stress_model='Linear',
                witnesses=[[-3.,0.,3.6880433085817326]],reporting_breakpoints=out['reporting_breakpoints'])
            result,receipt=run_probe(a.probe,xml,packet,folder/'readback');d=result['diagnostics'];q=d['physical'];op=d['operators']
            require(result['governing_solves']==result['physical_commits']==0 and result['initialization_work']==[0,0],
                'comparison performed a solve, physical commit or initialization Lift')
            require(op['selected_outer']=='Original' and all(max(map(abs,b['difference_load']),default=0)==0 for b in op['blocks']),
                'faithful same-state operator identity changed')
            physical[name][n]=q
            case['observations'][str(n)]=dict(time=native['time'],checkpoint=str(source),checkpoint_copy=str(cp),
                checkpoint_sha256=h['sha256'],state_sha256=state_hash(native['coupled']),
                native_history_sha256=[state_hash(v) for v in native['history']],
                native_history_times_sha256=state_hash(native['history_times']),receipt=receipt,
                physical=compact_diagnostic(compact(q)),error_to_own_reference=reference_errors(q),
                covariance_minimum_inner=d['selected_state']['minimum_inner'],
                covariance_minimum_outer=d['selected_state']['minimum_outer'],
                fixed_witnesses=d['selected_state']['one_sided_witnesses'])
            atomic(a.output/'summary.json',out)
        if 0 in physical[name]:
            for n,q in physical[name].items():
                change=differences(q,physical[name][0]);case['observations'][str(n)].update(
                    drift_from_actual_h0=change,common_normalized_drift=normalized(change))
        atomic(a.output/'summary.json',out)
    for name in ('half-dt','pressure-p2'):
        rows=out['same_time_differences'][name]={};dt=out['cases'][name]['dt']
        for n,q in physical[name].items():
            other=round(n*dt/.000625)
            if other not in physical['baseline']:continue
            change=differences(q,physical['baseline'][other]);rows[str(n)]=dict(time=q['time'],
                dimensional_field_difference=change,common_normalized_difference=normalized(change),
                mean_wall_stress_difference=q['mean_wall_stress']-physical['baseline'][other]['mean_wall_stress'])
    out.update(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat())
    atomic(a.output/'summary.json',out);print(a.output/'summary.json')


if __name__=='__main__':main()
