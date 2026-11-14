#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H18 measured outcomes; no solver, projection or scientific reevaluation."""
import argparse,csv,gzip,hashlib,json,math,os,shutil,statistics,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj_h17_observer import native_load
from run_wj_h15 import size
from wj3_measure import checkpoint_payload

HORIZON=505.02466521278404
WINDOW=[447.8483832973273,HORIZON]
N=202010

def table(path,rows):
    if rows:
        with path.open('x',newline='')as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)

def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');entries=book['entries'];copied=[];heavy=[]
    def copy(p,name,compressed=False):
        p=Path(p)
        if not p.exists():return None
        q=out/name;q.parent.mkdir(parents=True,exist_ok=True)
        if compressed:
            with p.open('rb')as f,gzip.open(q,'wb',compresslevel=6)as g:shutil.copyfileobj(f,g)
            h=hashlib.sha256()
            with gzip.open(q,'rb')as g:
                for b in iter(lambda:g.read(1024**2),b''):h.update(b)
            assert h.hexdigest()==identity(p)['sha256']
        else:shutil.copyfile(p,q);assert identity(q)['sha256']==identity(p)['sha256']
        copied.append(dict(original=identity(p),stored=identity(q),compression='gzip lossless'if compressed else'none'));return q
    copy(root/'accounting.json','accounting-at-publication.json')
    for n in ('source-intent.json','execution-config.json','target-admission.json','inputs/selected-bootstrap.json'):
        copy(root/n,'controls/'+n.replace('/','--'))
    for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','runtime-qualified-v1.json','target-B-measured-admission-v1.json','pending-restored-native-map-parity-v1.json','component-checks.json','observer-small-preflight.json','target-cache-compatibility.json','small-cache-compatibility.json'):
        copy(root/'science'/n,'controls/'+n+'.gz',True)
    timeline=[];bindings=[];audits=[];moments=[];signed=[];timings=[];storage=[]
    for name in ('A','B'):
        folder=root/'science/target'/name;directory=folder/'monitor'
        if not folder.exists():continue
        for n in ('qualification.json','case.evolution.json','execution-diff.json','intent.json'):
            copy(folder/n,f'target-{name}/'+n)
        for n in ('startup.verification-commit.json','latest-native-completed.json','latest-durable-checkpoint.json','latest-fully-verified-checkpoint.json','qualified-cursor.json'):
            copy(directory/n,f'target-{name}/'+n)
        tr=root/f'operations/target-{name}-v1/timing.jsonl'
        scopes=[json.loads(s)for s in tr.read_text().splitlines()]if tr.exists()else[]
        for s in scopes:
            if s['name']in ('native.step','native.advance','native.acceptance_audit','native.physical_samples','monitor.checkpoint_and_publish','monitor.independent_wait','monitor.qualification_commit','startup.object','startup.initialize'):
                timings.append(dict(run=name,**{k:s.get(k)for k in ('step','time','name','wall_seconds','exclusive_wall_seconds','process_cpu_seconds')}))
        def duration(step,key):return sum(s['wall_seconds']for s in scopes if s['step']==step and s['name']==key)
        for n in range(1,5):
            p=directory/f'step-{n}.completion.json'
            if not p.exists():continue
            c=load(p);r=native_load(verify(c['report']));cp=checkpoint_payload(verify(c['checkpoint']).parent,include_reference=True);attempt=r['pimple_attempts'][0]
            full=r['native_full_audit_status']=='PASS';ready=directory/f'step-{n}.ready.json';result_path=directory/f'step-{n}.ready.result.json';commit=directory/f'step-{n}.qualification-commit.json'
            result=load(result_path)if result_path.exists()else None
            tele=next((json.loads(s)for s in (directory/'telemetry.jsonl').read_text().splitlines()if json.loads(s)['absolute_step']==n),{})if(directory/'telemetry.jsonl').exists()else{}
            row=dict(run=name,step=n,time=r['time'],native_completed=True,online=r['online_check_status'],native_full_audit=r['native_full_audit_status'],independent='PASS'if commit.exists()else'NOT_SCHEDULED'if not full else'UNRESOLVED',historical_H16_numerical_parity='PASS'if(folder/'qualification.json').exists()else'UNRESOLVED',
                bound_calls=attempt['bound_calls'],bound_activations=attempt['bound_activations'],material_corrections=attempt['additional_correction_or_restoration_solves'],native_compute_seconds=tele.get('native_compute_seconds'),
                inclusive_step_seconds=duration(n,'native.step'),native_full_audit_seconds=duration(n,'native.acceptance_audit'),physical_samples_seconds=duration(n,'native.physical_samples'),recovery_publication_seconds=duration(n,'monitor.checkpoint_and_publish'),
                observer_wait_seconds=duration(n,'monitor.independent_wait'),independent_seconds=result['measured_independent_seconds']if result else None,verification_commit_seconds=duration(n,'monitor.qualification_commit'),
                full_audit_timing_calls=sum(s['name']=='native.acceptance_audit'and s['step']==n for s in scopes),physical_samples_timing_calls=sum(s['name']=='native.physical_samples'and s['step']==n for s in scopes),observer_ready_files=int(ready.exists()),rich_k_snapshot=int('k_snapshot'in attempt),raw_corrected_k_equal=attempt['k_bound']['raw_sha256']==attempt['k_bound']['corrected_sha256'])
            timeline.append(row)
            bindings.append(dict(run=name,step=n,time=r['time'],contained_state_sha256=vector_hash(cp['state']),checkpoint_file=c['checkpoint'],callback=c['callback'],producer=c['report'],completion=identity(p),reference_sha256=vector_hash(cp['fixed_reference']),defect_sha256=vector_hash(cp['fixed_defect']),native_history_sha256=[vector_hash(h)for h in cp['native_history']],native_history_times=cp['native_history_times'],callback_history_sha256=c['binding']['callback_history_sha256'],original_qh0_sha256=c['binding']['original_qh0_sha256'],physical_signature_sha256=c['binding']['physical_signature_sha256'],alpha_bits=c['binding']['alpha_bits']))
            payload=[Path(c[k]['path'])for k in ('checkpoint','callback','report')]+[p]
            b=sum(q.stat().st_size for q in payload);storage.append(dict(run=name,step=n,role='scientific+recovery'if full else'recovery',native_recovery_bytes=b,native_recovery_files=len(payload),independent_result_bytes=result_path.stat().st_size if result_path.exists()else 0,inventory_bytes=(directory/f'step-{n}.ready.inventory.json').stat().st_size if result else 0))
            heavy.extend(identity(q)for q in payload)
            copy(p,f'target-{name}/step-{n}.completion.json');copy(commit,f'target-{name}/step-{n}.qualification-commit.json')
            signed.append(dict(run=name,step=n,stage_audits=r['stage_audits'],native_hard_audit=r['native_hard_audit'],pimple=attempt))
            if result:
                audits.append(dict(run=name,role='completed',step=n,time=r['time'],status=result['status'],summary=result['summary'],coverage=result['coverage'],cache=result['cache'],inventory_sha256=result['inventory_sha256']))
                for region,v in result['all_regional_Q1_moments'].items():
                    for i,value in enumerate(v['moments']):moments.append(dict(run=name,role='completed',step=n,region=region,test=i,signed_moment=value))
                heavy.extend((identity(result_path),identity(directory/f'step-{n}.ready.inventory.json')))
        p=directory/'startup.result.json'
        if p.exists():
            r=load(p);startup=native_load(directory/'startup.json');audits.append(dict(run=name,role='restored_input',step=startup['steps'],time=startup['time'],status=r['status'],summary=r['summary'],coverage=r['coverage'],cache=r['cache'],inventory_sha256=r['inventory_sha256']))
            for region,v in r['all_regional_Q1_moments'].items():
                for i,value in enumerate(v['moments']):moments.append(dict(run=name,role='restored_input',step=startup['steps'],region=region,test=i,signed_moment=value))
            heavy.extend((identity(p),identity(directory/'startup.inventory.json')))
        p=folder/'view/result.json';copy(p,'final-strict-result.json.gz',True)
    table(out/'target-timeline.csv',timeline);table(out/'disjoint-inclusive-timing.csv',timings);table(out/'all-signed-Q1-moments.csv',moments);table(out/'checkpoint-storage.csv',storage)
    atomic(out/'state-history-bindings.json',bindings)
    q=out/'audited-state-summaries.json';atomic(q,audits)
    # Lossless compact storage retains all signed numerical scientific maps.
    with q.open('rb')as f,gzip.open(str(q)+'.gz','wb')as g:shutil.copyfileobj(f,g)
    q.unlink()
    q=out/'signed-native-stage-data.json';atomic(q,signed)
    with q.open('rb')as f,gzip.open(str(q)+'.gz','wb')as g:shutil.copyfileobj(f,g)
    q.unlink()
    ops=[]
    for label,v in entries.items():
        ops.append(dict(operation=label,phase=v['phase'],status=v['status'],seconds=v.get('seconds',v['reserved_seconds']),reserved_seconds=v['reserved_seconds'],attempts=v.get('charged_attempts',0),peak_tree_rss_bytes=v.get('peak_tree_rss_bytes',0)))
        if v['status']!='RESERVED':
            for n in ('request.json','receipt.json','attempt-counts.json'):copy(root/'operations'/label/n,'operations/'+label+'/'+n)
    table(out/'operations.csv',ops)
    for name in ('legacy','A','B','fault'):
        for n in ('qualification.json','case.evolution.json','execution-diff.json'):copy(root/f'science/small/{name}'/n,'small/'+name+'-'+n)
    copy(root/'science/component-checks/native.json','controls/pending-online-native-checks.json')
    def passed(label):return entries.get(label,{}).get('status')=='PASS'
    full=len(timeline)==6 and all(r['online']=='PASS'and r['historical_H16_numerical_parity']=='PASS'for r in timeline)
    strict=any((root/f'science/target/{name}/view/result.json').exists()for name in ('B','A','final-original-CP0'))
    gates=dict(build='PASS',focused_python='PASS'if passed('focused-python-v1')else'NOT_RUN',small_guards_import_protocol='PASS'if passed('pending-online-components-v1')else'NOT_RUN',
        successful_small_attempts=sum(load(root/f'science/small/{n}/qualification.json')['actual_attempts']for n in ('legacy','A','B')if(root/f'science/small/{n}/qualification.json').exists()),
        deliberate_small_failed_attempts=load(root/'science/small/fault/qualification.json')['actual_attempts']if(root/'science/small/fault/qualification.json').exists()else 0,
        target_reserved_attempts=book.get('target_reserved_attempts',0),target_actual_attempts=book.get('target_charged_attempts',0),target_blocks=len(book.get('target_cases',[])),target_replay='PASS'if full else'PARTIAL',
        final_strict='PASS'if strict else'NOT_RUN',fresh_CP0=0,projections_initializations=0,new_scientific_steps=0,
        full_audited_states=[dict(run=a['run'],role=a['role'],step=a['step'],status=a['status'])for a in audits],
        intermediate_covariance_principal='NOT_SCHEDULED; exact H16 state reuse is separately recorded',quality_threshold=.002068763923454,Q1_absolute_gate=2e-15,
        temporal_adequacy='UNRESOLVED unchanged',duration_admission='CLOSED',experimental_validation='NOT_RUN')
    gates['scope_result']='PASS'if full and strict and all(a['status']=='PASS'for a in audits)and gates['small_guards_import_protocol']=='PASS'else'PARTIAL'
    atomic(out/'component-gates.json',gates)
    if full:
        warm=[r for r in timeline if r['run']=='A'and r['step']in (2,3)]
        core=statistics.mean(r['native_compute_seconds']for r in warm);recovery=statistics.mean(r['recovery_publication_seconds']for r in warm)
        terminal=next(r for r in timeline if r['run']=='A'and r['step']==4)
        scientific=terminal['native_full_audit_seconds']+terminal['physical_samples_seconds']+terminal['observer_wait_seconds']+terminal['verification_commit_seconds']
        first_terminal=math.ceil(WINDOW[0]/.0025);fulls=[s for s in range(1,N+1)if s==N or s%(50 if s>=first_terminal else 100)==0]
        forecast=(N*core+math.ceil(N/10)*recovery+len(fulls)*scientific)/86400
        rec=statistics.mean(v['native_recovery_bytes']for v in storage if v['role']=='recovery')
        sci=statistics.mean(v['native_recovery_bytes']+v['independent_result_bytes']+v['inventory_bytes']for v in storage if v['role']=='scientific+recovery')
        telemetry=sum((root/f'science/target/{n}/monitor/telemetry.jsonl').stat().st_size for n in ('A','B'))/6
        estimate=dict(horizon=HORIZON,window=WINDOW,steps=N,metadata_dt=.0025,ceiling_step_time=N*.0025,proposal_recovery_stride=10,proposal_full_stride=100,proposal_terminal_stride=50,
            proposed_scientific_checkpoints=len(fulls),proposed_recovery_checkpoints=math.ceil(N/10),measured_nondue_native_compute_seconds=core,measured_recovery_publication_seconds=recovery,
            measured_full_audit_plus_inclusive_wait_commit_seconds=scientific,independent_not_added_to_wait=True,conditional_duration_days=forecast,H16_every_step_inclusive_days=211.896222547,
            measured_recovery_bytes=rec,measured_scientific_package_bytes=sci,measured_telemetry_bytes_per_step=telemetry,
            archive_estimate_uncompressed_bytes=len(fulls)*sci+(3+10)*rec+N*telemetry,storage_scope='scientific archives plus newest3 recoveries and conservative10 pending-recovery segment; immutable configuration/source and startup/terminal extras separate; not an admitted archive',
            limitations='one short replay; CPU not hidden in observer; no established speedup statistics, scientifically accepted cadence, application tolerance or long-run affordability; checkpoint compression measured separately')
        atomic(out/'conditional-cost-storage.json',estimate)
        # Actual lossless compression benchmark on one recovery and one science
        # package; original checkpoint/report/callback files stay unchanged.
        compressed=[];dest=out/'compression-measurements';dest.mkdir()
        for name,n in (('recovery',2),('scientific',4)):
            c=load(root/f'science/target/A/monitor/step-{n}.completion.json')
            payload={role:c[role]for role in ('checkpoint','callback','report')}
            if name=='scientific':
                for role in ('inventory','result'):payload[role]=identity(root/f'science/target/A/monitor/step-{n}.ready.{role}.json')
            for role,pin in payload.items():
                p=verify(pin);q=dest/(name+'-'+role+'.gz');t=time.monotonic()
                with p.open('rb')as f,gzip.open(q,'wb',compresslevel=6)as g:shutil.copyfileobj(f,g)
                raw=gzip.open(q,'rb').read();assert hashlib.sha256(raw).hexdigest()==pin['sha256']
                compressed.append(dict(package=name,role=role,original=pin,compressed=identity(q),seconds=time.monotonic()-t,ratio=q.stat().st_size/p.stat().st_size))
        atomic(out/'compression-measurements.json',compressed)
        # Avoid copying optional large gzip checkpoint aliases into Git transport.
        heavy.extend(identity(p)for p in dest.iterdir());shutil.move(str(dest),str(root/'science/compression-measurements'))
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7,3.8),layout='constrained');x=list(range(len(timeline)));base=[r['native_compute_seconds']-r['native_full_audit_seconds']-r['physical_samples_seconds']for r in timeline]
        ax.bar(x,base,label='advance + online checks');b=base[:]
        for key,label in (('native_full_audit_seconds','full native audit'),('physical_samples_seconds','scientific physical samples'),('recovery_publication_seconds','recovery publication'),('observer_wait_seconds','observer barrier (includes assessment)'),('verification_commit_seconds','verification commit')):
            y=[r[key]for r in timeline];ax.bar(x,y,bottom=b,label=label);b=[u+v for u,v in zip(b,y)]
        ax.set(xticks=x,xticklabels=[f"{r['run']}:{r['step']}"for r in timeline],ylabel='Wall time (s)',xlabel='Replay run : absolute step',title='Fixed6 mm control, dt=.0025; all stages replay H16');ax.legend(fontsize=7,ncol=2)
        fig.savefig(out/'checkpoint-cost.png',dpi=160);plt.close(fig)
        selected=[v for v in storage if v['run']=='A'];fig,ax=plt.subplots(figsize=(6.2,3.4),layout='constrained')
        x=[v['step']for v in selected];native=[v['native_recovery_bytes']/1024**2 for v in selected];ind=[(v['independent_result_bytes']+v['inventory_bytes'])/1024**2 for v in selected]
        ax.bar(x,native,label='native state / callback / producer / manifest');ax.bar(x,ind,bottom=native,label='current scientific inventory and result')
        ax.set(xticks=x,xlabel='Absolute replay step',ylabel='Retained payload (MiB)',title='H18 pilot retains every recovery; full audit at step4');ax.legend(fontsize=7)
        fig.savefig(out/'checkpoint-storage.png',dpi=160);plt.close(fig)

    # Source/ABI provenance remains remote in full; transport precise compact pins.
    compiled=load(root/'science/compiled-source-objects.json')
    atomic(out/'compiled-source-summary.json',dict(native_source_commit=compiled['native_source_commit'],compiler=compiled['compiler'],objects=len(compiled['objects']),dependencies=compiled['local_dependencies'],full_record=identity(root/'science/compiled-source-objects.json'),executables=load(root/'science/runtime-build-v1.json')['executables']))
    # Reused H16 profile/fetch/error context is explicit, never another H18 scan.
    old=root.parent/'wj-h16-repaired-monitor-20261114';idx=load(old/'authoritative-index-v1.json');manifest=load(verify(idx['publication']));pins={v['path']:v for v in manifest['files']};reused=[]
    for n in ('summary.json','timeline.csv'):
        p=old/'publication-v1/analysis'/n
        if str(p)in pins:verify(pins[str(p)]);copy(p,'reused-H16/'+n);reused.append(pins[str(p)])
    atomic(out/'reused-context.json',dict(scope='sealed H16 final-union reference/drift/profile/fetch context; H18 exact state/history parity separate; no additional H18 target audit batch',files=reused,new_scientific_steps=0))
    atomic(out/'heavy-manifest.json',dict(files=heavy,scope='immutable selected inputs and outputs; files remain at original Linux paths'))
    atomic(out/'copied-evidence.json',copied)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='immutable publication; open accounting snapshot is not final closure'))
    print(json.dumps(gates),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
