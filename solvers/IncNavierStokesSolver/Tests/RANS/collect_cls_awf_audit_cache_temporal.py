#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact receipts and checkpoint-verified live frontiers; never resume a job."""
import argparse
from datetime import datetime,timezone,timedelta
import hashlib
import json
from pathlib import Path
import subprocess
import time
import statistics
from continue_cls_epawf_pimple import atomic,header,process_identity,load_progress
from run_cls_epawf_picard import load,sha,state_hash,require


def witness(p):
    keys=('coordinate','geometry','element','side','inner','inventory','gradient_element_region','constitutive_region','gradient_region_matches_law',
          'k','epsilon','sampled_K','endpoint_K','velocity_gradient','divergence','epsilon_semantics','constitutive_branch')
    out={k:p[k] for k in keys if k in p}
    if 'selected' in p:out['eigenvalues']=p['selected']['eigenvalues']
    if p.get('same_state_Original') is not None:out['Original_eigenvalues']=p['same_state_Original']['eigenvalues']
    return out


def archived(root):
    d=load(root/'summary.json');out={k:v for k,v in d.items() if k!='cases'};out['cases']={}
    for name,row in d.get('cases',{}).items():
        full=load(root/name/'stages/result.json');audits=full['complete_audits']
        out['cases'][name]={k:v for k,v in row.items() if k not in ('stages','independent_tensor_checks')}
        out['cases'][name]['stages']={stage:dict(complete_audit=audits[stage],minimum_inner=witness(q['minimum_inner']),
            minimum_outer=witness(q['minimum_outer'])) for stage,q in full['diagnostics'].items()}
        checks=row.get('independent_tensor_checks',[])
        out['cases'][name]['independent_tensor_check_count']=len(checks)
        out['cases'][name]['maximum_independent_tensor_error']=max((r['tensor_error'] for r in checks),default=0.)
    return out


def campaign(path,snapshots):
    for _ in range(20):
        try:
            payload=(path/'progress.json').read_bytes();p=json.loads(payload);cp=Path(p['checkpoint']);h=header(cp,True)
            checkpoint_bytes=(cp/'State.bin').read_bytes()
            if not (h['sha256']==p['checkpoint_sha256']==hashlib.sha256(checkpoint_bytes).hexdigest() and h['steps']==p['accepted_steps'] and h['time']==p['accepted_time']):
                time.sleep(.05);continue
            break
        except (FileNotFoundError,json.JSONDecodeError):time.sleep(.05)
    else:raise RuntimeError('checkpoint rotation prevented a stable status capture')
    snapshot=snapshots/(str(h['steps'])+'-'+h['sha256'])
    snapshot.mkdir(parents=True,exist_ok=True);file=snapshot/'State.bin'
    if file.exists():require(sha(file)==h['sha256'],'immutable status snapshot differs')
    else:file.write_bytes(checkpoint_bytes)
    seal=path/'runtime-handoff.json';runtime=path/'runtime/identity.json';identity=load(runtime)['frozen']
    driver=process_identity(p['pid']);expected=p.get('native_identity');native=process_identity(expected['pid']) if expected else None
    status_file=path.parent/'status.json';status=load_progress(status_file) if status_file.exists() else {}
    expected_supervisor=status.get('supervisor_identity');supervisor=process_identity(expected_supervisor['pid']) if expected_supervisor else None
    return dict(path=str(path),status='handed_off' if seal.exists() else p['status'],accepted_steps=h['steps'],accepted_time=h['time'],
        checkpoint=str(cp),checkpoint_sha256=h['sha256'],coupled_state_sha256=state_hash(h['coupled_state']),
        immutable_checkpoint=str(snapshot),
        native_history_sha256=[state_hash(v) for v in h['history']],native_history_times_sha256=state_hash(h['history_times']),
        source_progress_sha256=hashlib.sha256(payload).hexdigest(),driver_identity=driver,native_identity=native,supervisor_identity=supervisor,
        supervisor_start_matches=bool(supervisor and supervisor==expected_supervisor),native_identity_matches=bool(native and native==expected),
        driver_start_matches=bool(driver and str(driver['start'])==str(p.get('process_start'))),
        executable=identity['path'],executable_sha256=identity['sha256'],library_manifest=str(runtime),library_manifest_sha256=sha(runtime),
        next_target_step=p.get('next_target'),next_target_time=p.get('next_target')*.000625 if p.get('next_target') else None,
        completed_targets=sorted(set(p.get('completed_targets',[]))|{n for n in p.get('requested_targets',[]) if n<=h['steps']}),
        recorded_successor_completed_targets=p.get('completed_targets'),failure=p.get('failure'),prefix_steps=p.get('prefix_steps',0),
        cumulative_bound_integral=p.get('cumulative_bound_integral'),audit_revision=p.get('audit_revision','LegacyInterfaceV1'),
        evaluation_reuse=p.get('evaluation_reuse','Off'),handoff=load(seal) if seal.exists() else None)


def timing_statistics(root):
    """Summarize existing alternating receipts; do not execute another replay."""
    out={}
    for case in ('BudgetH2','UpstreamOne'):
        modes={}
        for mode in ('Off','ExactStepV1'):
            folders=[root/case/f'pair-{i}-{mode}' for i in (1,2,3)]
            rows=[load(folder/'receipt.json') for folder in folders]
            median=lambda values:statistics.median(values)
            groups={name:{key:median(row['timing']['groups'][name][key] for row in rows)
                for key in ('cpu_per_step','exclusive_cpu_per_step','wall_per_step','exclusive_wall_per_step')}
                for name in rows[0]['timing']['groups']}
            first=[];continued=[]
            for folder in folders:
                steps=[json.loads(line) for line in (folder/'timing.jsonl').read_text().splitlines() if json.loads(line)['name']=='native.step']
                first.append(steps[0]['wall_seconds']);continued.append(statistics.mean(v['wall_seconds'] for v in steps[1:]))
            modes[mode]=dict(median_driver_seconds=median(v['driver_wall_seconds'] for v in rows),
                median_native_cpu_seconds=median(v['native_process_cpu_seconds'] for v in rows),
                median_native_wall_seconds=median(v['native_process_wall_seconds'] for v in rows),
                median_peak_rss_KiB=median(v['peak_rss_KiB'] for v in rows),
                median_first_physical_step_wall_seconds=median(first),median_continued_step_wall_seconds=median(continued),
                median_native_wall_outside_physical_steps=median(v['native_process_wall_seconds']-v['timing']['groups']['step/native.step']['wall'] for v in rows),
                groups=groups)
        out[case]=dict(steps_per_pair=8,alternating_pairs=3,modes=modes,
            interpretation='Per-group paired medians are descriptive, not additive. Each original timing tree separately passes exclusive partition checks. First-step and outside-step restart costs are distinct.')
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--campaign',action='append',default=[]);p.add_argument('--status-only',action='store_true',help='Refresh only dated live identities/checkpoint copies; do not repeat completed scientific analysis')
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    out=dict(schema=1,captured_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),linux_root=str(a.root),
        audit_revision='NativeElementSideV2',audit={},performance={},temporal={},long_windows={},handoffs={})
    if a.status_only:
        statuses={name:campaign(Path(path),a.output/'checkpoints'/name) for name,path in (item.split('=',1) for item in a.campaign)}
        atomic(a.output/'latest-status.json',dict(captured_Asia_Shanghai=out['captured_Asia_Shanghai'],campaigns=statuses));return
    for folder in ('archived-audit-v2','archived-half-dt-v2'):
        if (a.root/folder/'summary.json').exists():out['audit'][folder]=archived(a.root/folder)
    for folder in ('legacy-report-H2-Saturated-Original','legacy-report-H05-Linear-Regularized'):
        if (a.root/folder/'summary.json').exists():
            d=load(a.root/folder/'summary.json');d.pop('independent_tensor_checks',None);out['audit'][folder]=d
    for folder,file in [('interface-fixture','result.json'),('native-one-audit-repair','validation.json')]:
        if (a.root/folder/file).exists():out['audit'][folder]=load(a.root/folder/file)
    for file in ('rejection-rollback-v2.json','rejection-rollback-final.json','per-step-checkpoint-H2.json','per-step-checkpoint-final.json'):
        if (a.root/file).exists():out['audit'][file]=load(a.root/file)
    for folder in ('cache-contract-final','benchmarks-final','benchmarks-v2','legacy-graded-reconstruction-pair'):
        if (a.root/folder/'summary.json').exists():out['performance'][folder]=load(a.root/folder/'summary.json')
    if (a.root/'benchmarks-final/summary.json').exists():out['performance']['timing_statistics']=timing_statistics(a.root/'benchmarks-final')
    for root in a.root.glob('handoff-*'):
        if (root/'handoff.json').exists():out['handoffs'][root.name]=load(root/'handoff.json')
    operator_folder=next(name for name in ('operators-through-T2-verified','operators-through-T2','operators-long-final') if (a.root/name/'summary.json').exists())
    for folder,key in [('temporal-analysis','temporal'),(operator_folder,'long_windows')]:
        file=a.root/folder/'summary.json'
        if file.exists():
            d=load(file)
            if key=='temporal':
                for rows in d.get('observations',{}).values():
                    for q in rows.values():
                        q['common_covariance_witnesses']=[witness(x) for x in q.get('common_covariance_witnesses',[])]
                        q['covariance']=witness(q['covariance'])
            else:
                for q in d.get('cases',{}).values():
                    q['selected_covariance']=witness(q['selected_covariance']);q['Original_covariance']=witness(q['Original_covariance'])
            out[key]=d
    window_folder=next(name for name in ('long-window-through-T2-final','long-window-through-T2-verified','long-window-through-T2','long-window-analysis') if (a.root/name/'summary.json').exists())
    for file,key in [('temporal-validation-final.json','temporal_confirmation'),(window_folder+'/summary.json','window_classification')]:
        if (a.root/file).exists():out[key]=load(a.root/file)
    statuses={}
    for item in a.campaign:
        name,path=item.split('=',1);statuses[name]=campaign(Path(path),a.output/'checkpoints'/name)
    atomic(a.output/'latest-status.json',dict(captured_Asia_Shanghai=out['captured_Asia_Shanghai'],campaigns=statuses))
    indices=[]
    for pattern in ('**/summary.json','**/receipt.json','**/identity.json','**/handoff.json','*rollback*.json','temporal-policy.json',
                    'source-runtime.json','temporal-validation-final.json','per-step-checkpoint*.json','*.log'):
        for file in sorted(a.root.glob(pattern)):
            if file.is_relative_to(a.output):continue
            indices.append(dict(path=str(file),sha256=sha(file),bytes=file.stat().st_size))
    unique={r['path']:r for r in indices};atomic(a.output/'artifact-index.json',dict(root=str(a.root),
        indices=list(unique.values()),scope='Detailed receipts index every native packet/matrix, complete fields, checkpoint payload and library closure; these summaries do not replace native restart data.'))
    out['artifacts']=dict(index=str(a.output/'artifact-index.json'),sha256=sha(a.output/'artifact-index.json'))
    atomic(a.output/'native.json',out);print(a.output/'native.json')

if __name__=='__main__':main()
