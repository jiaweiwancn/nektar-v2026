#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H17 outcomes from receipts, never from an assumed successful matrix."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import shutil
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj3_measure import checkpoint_payload

HORIZON=505.02466521278404
WINDOW=[447.8483832973273,505.02466521278404]
STEPS=202010

def table(path,records):
    if not records:return
    with path.open('x',newline='')as f:
        w=csv.DictWriter(f,list(records[0]));w.writeheader();w.writerows(records)

def seal(directory):
    files=[identity(p)for p in sorted(directory.rglob('*'))if p.is_file()and p.name!='manifest.json']
    atomic(directory/'manifest.json',dict(schema='WJH17CompactManifestV1',files=files,scope='immutable publication snapshot; closure is separate'))

def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');entries=book['entries']
    def copy(src,name):
        src=Path(src)
        if not src.exists():return None
        dest=out/name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dest);assert identity(dest)['sha256']==identity(src)['sha256'];return identity(dest)
    def passed(label):return entries.get(label,{}).get('status')=='PASS'
    for name in ('accounting.json','source-intent.json','execution-config.json','audit-admission.json','target-admission.json','inputs/selected-bootstrap.json'):
        copy(root/name,'controls/'+name.replace('/','--'))
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','target-cache-compatibility.json','small-cache-compatibility.json','runtime-qualified-v1.json','audit-comparison.json'):
        p=root/'science'/name
        if p.exists():
            target=out/'controls'/(name+'.gz')
            target.parent.mkdir(parents=True,exist_ok=True)
            with p.open('rb')as f,gzip.open(target,'wb')as g:shutil.copyfileobj(f,g)
    operations=[]
    for label,v in entries.items():
        operations.append(dict(operation=label,phase=v['phase'],status=v['status'],reserved_seconds=v['reserved_seconds'],charged_seconds=v.get('seconds',v['reserved_seconds']),reserved_attempts=v.get('reserved_attempts',0),actual_or_unknown_attempts=v.get('charged_attempts',0),peak_tree_rss_bytes=v.get('peak_tree_rss_bytes',0)))
        if v['status']=='RESERVED':continue  # do not seal active logs/receipts
        for name in ('request.json','receipt.json','attempt-counts.json','native-exec-samples.json'):
            copy(root/'operations'/label/name,'operations/'+label+'/'+name)
        p=root/'operations'/label/'command.log'
        if p.exists():
            dest=out/'operations'/label/'command.log.gz';dest.parent.mkdir(parents=True,exist_ok=True)
            with p.open('rb')as f,gzip.open(dest,'wb')as g:shutil.copyfileobj(f,g)
    table(out/'operations.csv',operations)
    audit_records=[];audit_payloads=[];counter_summaries=[]
    for policy in ('Off','PatchExactV1'):
        folder=root/'science/audits'/policy
        copy(folder/'native.json.progress.json','audit-'+policy+'-progress.json')
        p=folder/'native.json'
        if not p.exists():continue
        v=load(p);audit_payloads.append(identity(p))
        counter_summaries.append(dict(policy=policy,immutable_catalog=v['immutable_rows'],scoped_dynamic_reuse=v['dynamic_reuse'],scope='whole batch, including strict imports and affected-output checks; audit-only visit commitments are separate',cold_role=v['results'][0]['operation_seconds'],return_to_A_role=v['results'][-1]['operation_seconds']))
        with p.open('rb')as f,gzip.open(out/('complete-'+policy+'.json.gz'),'wb')as g:shutil.copyfileobj(f,g)
        for n,r in enumerate(v['results']):
            c=r['complete'];s=r['rows_after']
            audit_records.append(dict(policy=policy,ordinal=n,absolute_step=r['step'],time=r['time'],alpha=r['alpha'],import_seconds=r['import_seconds'],audit_seconds=r['audit_seconds'],commitment_seconds=c['commitment_wall_seconds'],audit_less_commitment_seconds=r['audit_seconds']-c['commitment_wall_seconds'],operation_seconds=r['operation_seconds'],audit_row_requests=c['ordered_rows']['count'],observation_visits=c['ordered_observations']['count'],row_stream=c['ordered_rows']['root_sha256'],observation_stream=c['ordered_observations']['root_sha256'],audit_hash=c['audit_numerical_sha256'],catalog_constructs=s['constructs'],catalog_hits=s['hits'],catalog_misses=s['misses'],catalog_evictions=s['evictions'],catalog_bypasses=s['bypasses'],resident_entries=s['resident_entries'],catalog_live_bytes=s['live_bytes'],retained_bytes=s['retained_after_eviction_bytes']))
    table(out/'saved-audit-cost.csv',audit_records)
    atomic(out/'native-row-counters.json',dict(policies=counter_summaries,note='Off has no persistent catalog; its actual constructs/hits are scoped_dynamic_reuse counters, not the zero-valued catalog fields'))
    target=root/'science/target';stages=[];bindings=[]
    for n in range(1,5):
        directory=target/'monitor';p=directory/f'step-{n}.ready.json'
        if not p.exists():continue
        ready=load(p);r=load(verify(ready['report']));result_path=directory/f'step-{n}.ready.result.json';result=load(result_path)if result_path.exists()else None
        commit=directory/f'step-{n}.qualification-commit.json';s=r['immutable_row_reuse'];attempt=r['pimple_attempts'][0]
        stages.append(dict(absolute_step=n,time=r['time'],native_completed=True,qualified=commit.exists(),native_compute_seconds=ready['native_compute_seconds'],independent_seconds=result['measured_independent_seconds']if result else None,points=len(load(directory/f'step-{n}.ready.inventory.json')['points'])if result else None,owners=len(load(directory/f'step-{n}.ready.inventory.json')['owned'])if result else None,bound_calls=attempt['bound_calls'],bound_activations=attempt['bound_activations'],constructs=s['constructs'],hits=s['hits'],misses=s['misses'],evictions=s['evictions'],resident_entries=s['resident_entries'],live_bytes=s['live_bytes'],retained_bytes=s['retained_after_eviction_bytes']))
        for name in (f'step-{n}.ready.json',f'step-{n}.qualification-commit.json',f'step-{n}.native-replay-parity.json',f'step-{n}.independent-replay-parity.json',f'step-{n}.remaining.json'):
            copy(directory/name,'target/'+name)
        if result:
            atomic(out/f'target/step-{n}-physical-summary.json',dict(summary=result['summary'],all_regional_Q1_moments=result['all_regional_Q1_moments'],scope='same live H16 inventory',source=identity(result_path)))
        cp=checkpoint_payload(verify(ready['checkpoint']).parent,include_reference=True)
        assert vector_hash(cp['state'])==ready['binding']['state_sha256']
        bindings.append(dict(step=n,state=ready['binding']['state_sha256'],reference=ready['binding']['fixed_reference_sha256'],defect=ready['binding']['fixed_defect_sha256'],native_history=[vector_hash(v)for v in cp['native_history']],native_history_times=cp['native_history_times'],ready=identity(p),report=ready['report'],checkpoint=ready['checkpoint'],callback=ready['callback'],commit=identity(commit)if commit.exists()else None))
    table(out/'target-timeline.csv',stages)
    atomic(out/'heavy-input-bindings.json',dict(audits=audit_payloads,target=bindings,scope='files stay in original Linux roots; state-vector hashes differ from checkpoint file hashes'))
    if len(stages)==4 and all(v['qualified']for v in stages):
        # Reuse the already qualified final-union evaluation of the EXACT same
        # physical states under unchanged independent numerical functions.
        # This is copied H16 evidence, never a newly executed H17 scan.
        old=root.parent/'wj-h16-repaired-monitor-20261114'
        index=load(old/'authoritative-index-v1.json');manifest=load(verify(index['publication']))
        pins={v['path']:v for v in manifest['files']};reused=[]
        for name in ('summary.json','timeline.csv',*[f'step-{n}.json'for n in range(5)]):
            p=old/'publication-v1/analysis'/name;pin=pins[str(p)];verify(pin)
            copy(p,'reused-H16-final-union/'+name);reused.append(pin)
        atomic(out/'retrospective-reuse-binding.json',dict(status='REUSED_H16_PASS',new_H17_scans=0,new_scientific_steps=0,original_publication=index['publication'],copied=reused,per_stage_exact_native_and_independent_parity=[identity(target/f'monitor/step-{n}.independent-replay-parity.json')for n in range(1,5)],scope='H16 final union11456 points/17521 owners plus44800 physical Gauss rows; earlier online inventories were unchanged; original H16 progress counts in copied history remain historical'))
    for name in ('target-outcome.json','view/result.json','execution-diff.json','intent.json'):
        copy(target/name,'target/'+name.replace('/','--'))
    for case in ('A','B'):
        for name in ('qualification.json','execution-diff.json','view/result.json'):
            copy(root/f'science/small/{case}'/name,'small/'+case+'-'+name.replace('/','--'))
    copy(root/'science/row-fixture/native.json','controls/native-reuse-fixture.json')
    native_cases=sum(entries.get(label,{}).get('charged_attempts',0)for label in entries if entries[label].get('kind')in ('small','target'))
    audit_pass=(root/'science/audit-comparison.json').exists()and load(root/'science/audit-comparison.json')['status']=='PASS'
    small_pass=all((root/f'science/small/{c}/qualification.json').exists()and load(root/f'science/small/{c}/qualification.json')['status']=='PASS'for c in ('A','B'))
    outcome=load(target/'target-outcome.json')if(target/'target-outcome.json').exists()else dict(status='NOT_RUN',qualified=0,new_scientific_progress=0)
    strict=(target/'view/result.json').exists()and load(target/'view/result.json')['status']=='PASS'
    gates=dict(implementation_build='PASS'if(root/'science/runtime-build-v1.json').exists()else'UNRESOLVED',native_catalog_unit='PASS'if passed('native-row-catalog-unit-v1')else'NOT_RUN',native_reuse_fixture='PASS'if passed('native-reuse-readonly-fixture-v1')else'NOT_RUN',complete_saved_parity='PASS'if audit_pass else'UNRESOLVED'if book.get('full_audit_reserved')else'NOT_RUN',complete_full_audit_calls=book.get('full_audit_actual_calls',0),reserved_full_audit_calls=book.get('full_audit_reserved',0),small_governing='PASS'if small_pass else'NOT_RUN',governing_charged_attempts=native_cases,target=outcome,target_strict='PASS'if strict else'NOT_RUN',new_scientific_steps=0,performance='MEASURED'if len(stages)==4 else'NOT_RUN',affordability='NOT_QUALIFIED',H11_temporal='UNRESOLVED unchanged',physical_height='NOT_SELECTED',duration_preservation='NOT_RUN',experimental_validation='NOT_RUN')
    gates['scope_result']='PASS'if audit_pass and small_pass and strict and outcome['status']=='PASS'else'PARTIAL'
    atomic(out/'component-gates.json',gates)
    if len(stages)==4:
        native=sum(v['native_compute_seconds']for v in stages[1:])/3;assessment=sum(v['independent_seconds']for v in stages[1:])/3
        timings=[]
        for label,v in entries.items():
            if v.get('kind')=='target':
                p=root/'operations'/label/'timing.jsonl'
                if p.exists():timings.extend(json.loads(s)for s in p.read_text().splitlines())
        names=('native.step','startup.object','startup.initialize','native.acceptance_audit','audit.WJ.owned_rows_gradients','audit.WJ.velocity_accessor_validation','audit.WJ.scalar_row_accessor','audit.WJ.current_state_contractions','audit.WJ.principal_full_search','pimple.sweep','monitor.checkpoint_and_publish','monitor.independent_wait','monitor.qualification_commit')
        selected=[dict(step=v['step'],time=v['time'],name=v['name'],wall_seconds=v['wall_seconds'],cpu_seconds=v['process_cpu_seconds'],calls=v.get('calls',1),clock_scope=v.get('clock_scope','inclusive scope; do not add nested intervals'))for v in timings if v['name']in names]
        table(out/'target-native-timing.csv',selected)
        means={name:sum(v['wall_seconds']for v in selected if v['name']==name and v['step']in (2,3,4))/3 for name in names}
        inclusive=means['native.step']
        atomic(out/'conditional-cost.json',dict(horizon=HORIZON,terminal_window=WINDOW,steps=STEPS,warm_native_seconds=native,warm_assessment_seconds=assessment,warm_inclusive_loop_seconds=inclusive,warm_timing_means=means,native_only_days=native*STEPS/86400,observed_compute_plus_assessment_days=(native+assessment)*STEPS/86400,inclusive_loop_days=inclusive*STEPS/86400,H16_inclusive_days=211.896222547,scope='conditional arithmetic, not duration admission; inclusive loop includes publication/wait/commit; nested timers are not added; audit CPU0 is unmeasured; future witness growth/startup/terminal work unresolved'))
    figures(out,audit_records,stages)
    seal(out)
    atomic(root/'authoritative-index-v1.json',dict(publication=identity(out/'manifest.json'),component_gates=identity(out/'component-gates.json'),ledger_state='OPEN publication snapshot, separately closed later',new_scientific_steps=0))

def figures(out,audits,stages):
    if not audits and not stages:return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,3.4),layout='constrained')
    for policy in ('Off','PatchExactV1'):
        rows=[r for r in audits if r['policy']==policy]
        if rows:axes[0].plot([r['ordinal']for r in rows],[r['audit_less_commitment_seconds']for r in rows],'o-',label=policy)
    axes[0].set(xlabel='Saved audit ordinal (steps 0,1,2,3,4,0)',ylabel='Audit wall minus measured commitment time (s)',title='Identical saved-state operations');axes[0].legend()
    if stages:
        axes[1].plot([r['time']for r in stages],[r['native_compute_seconds']for r in stages],'o-',label='H17 native compute')
        axes[1].plot([r['time']for r in stages],[r['independent_seconds']for r in stages],'s-',label='Independent assessment')
        axes[1].axhline(83.0563815253,color='gray',linestyle='--',label='H16 warm native mean')
    axes[1].set(xlabel='T (Hh/U∞)',ylabel='Wall time (s)',title='H16 trajectory replay; fixed H_w=6 mm');axes[1].legend()
    fig.savefig(out/'row-reuse-cost.png',dpi=150);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
