#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Publish the bounded H19 study from retained packets; never launch a solver."""
import argparse,gzip,hashlib,json,math,os,shutil,statistics
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj_h17_observer import native_load
from wj_h3_analysis import algebra
from wj3_measure import checkpoint_payload
from publish_wj_h18 import table,HORIZON,WINDOW
from run_wj_h15 import size

OLD=Path('/home/jay/build/wj-h18-checkpoint-audits-20261114')

def witness_replay(report,arrays):
    output=[];physical=arrays['physical'];geometry=arrays['geometry']
    for name,w in report['wj_audit_packets'][0]['principal_witnesses'].items():
        point=np.asarray(w['position'],dtype='<f8');ids=np.flatnonzero(np.all(physical.view(np.uint64)==point.view(np.uint64),axis=1)&(geometry==w['geometry']))
        assert len(ids)==1,'native witness must select one exact physical owner'
        i=int(ids[0]);k=arrays['coarse_k'][i:i+1];e=arrays['coarse_epsilon'][i:i+1];G=arrays['coarse_G'][i:i+1];H=np.asarray(w['H_row_major']).reshape(3,3)
        R,nt,N=algebra(k,e,G);tangent=algebra(k,e,G.astype(complex)+1e-27j*H)[0].imag/1e-27
        dimensional=float(-np.sum(H*tangent[0]));normalized=dimensional/float(nt[0]);error=abs(normalized-w['replayed_over_nu_eff'])
        assert error<=2e-11 and abs(dimensional-w['minus_H_contract_dR'])<2e-13,'independent/native action replay differs'
        field_errors=dict(G=float(np.max(abs(G[0]-np.asarray(w['raw_G_row_major']).reshape(3,3)))),k=abs(float(k[0])-w['k']),epsilon=abs(float(e[0])-w['epsilon']),covariance=float(np.max(abs(np.linalg.eigvalsh(R[0])-w['covariance_eigenvalues']))))
        assert field_errors['G']<=2e-11 and field_errors['k']<=2e-11*.008809015145515192 and field_errors['epsilon']<=2e-11*.002315524054799254
        assert field_errors['covariance']<=2e-11*.008809015145515192
        output.append(dict(field_errors=field_errors,polarization=name,row=i,physical=w['position'],geometry=w['geometry'],normalized_native=w['replayed_over_nu_eff'],normalized_independent=normalized,
            normalized_error=error,dimensional_error=abs(dimensional-w['minus_H_contract_dR']),N_native=w['N'],N_independent=float(N[0]),nu_eff_native=w['nu_eff'],nu_eff_independent=float(nt[0])))
    return output

def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');heavy=[];copies=[]
    def copy(path,name,compress=False):
        p=Path(path)
        if not p.exists():return
        q=out/name;q.parent.mkdir(parents=True,exist_ok=True)
        if compress:
            with p.open('rb')as f,gzip.open(q,'wb',compresslevel=6)as g:shutil.copyfileobj(f,g)
            h=hashlib.sha256()
            with gzip.open(q,'rb')as g:
                for b in iter(lambda:g.read(1024**2),b''):h.update(b)
            assert h.hexdigest()==identity(p)['sha256']
        else:shutil.copyfile(p,q)
        copies.append(dict(original=identity(p),stored=identity(q),compression='gzip lossless'if compress else'none'))
    for p in ('source-intent.json','python-execution-source-v1.json','execution-config.json','execution-config-v1.json','target-admission.json','accounting.json','inputs/anchors-v2.json'):
        copy(root/p,'controls/'+p.replace('/','--'))
    for p in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json','runtime-qualified-v1.json','target-cache-compatibility.json','measured-admission-v1.json'):
        copy(root/'science'/p,'controls/'+p+'.gz',True)
    for p in ('lineage.json','qualification.json','execution-diff.json'):
        copy(root/'science/coarse-cp0'/p,'bootstrap/'+p+('.gz'if p=='qualification.json'else''),p=='qualification.json')
    copy(root/'science/component/qualification.json','component/qualification.json');copy(root/'science/component/result.json','component/native-result.json.gz',True)
    for p in ('qualification.json','outcome.json','execution-diff.json','intent.json'):
        copy(root/'science/coarse'/p,'coarse/'+p)
    moments=[];fields=[];errors=[];principals=[];pairs=[];summaries=[];witnesses=[];timeline=[];bindings=[];signed=[];storage=[]
    scopes_path=root/'operations/coarse-target-v1/timing.jsonl';scopes=[json.loads(s)for s in scopes_path.read_text().splitlines()]if scopes_path.exists()else[]
    def duration(step,name):return sum(s['wall_seconds']for s in scopes if s['step']==step and s['name']==name)
    for step in (0,1,2):
        folder=root/'science/coarse/monitor';stem='startup'if not step else f'step-{step}'
        p=folder/(stem+'.sensitivity.json')
        if not p.exists():
            pairs.append(dict(coarse_step=step,fine_step=2*step,time=step*.005,status='NOT_RUN',coarse_time_bits='',fine_time_bits=''));continue
        v=load(p);pair=v['pair'];pairs.append(dict(coarse_step=step,fine_step=2*step,time=pair['coarse_time'],status=v['status'],coarse_time_bits=pair['coarse_time_bits'],fine_time_bits=pair['fine_time_bits']))
        copy(p,stem+'-sensitivity.json.gz',True)
        inventory=load(verify(v['inventory']));arrays=np.load(verify(v['signed_arrays']))
        report=native_load(verify(pair['coarse_report']))
        witnesses.append(dict(step=step,replays=witness_replay(report,arrays)))
        heavy.extend((v['inventory'],v['signed_arrays'],identity(p)))
        for region,values in v['temporal'].items():
            for field,m in values.items():
                for kind in ('rms','maximum'):
                    f=m[kind];fields.append(dict(time=step*.005,region=region,field=field,metric=kind,scaled_difference=f['value'],limit=f['limit'],utilization=f['utilization'],floor=f['floor'],resolved=f['resolved'],status=f['status'],witness_row=m['row'],physical_measure=m['physical_measure']))
        for kind in ('representation','reference_error','drift'):
            branches={'initial':v[kind]}if kind=='representation'else v[kind]
            for branch,regions in branches.items():
                for region,values in regions.items():
                    for field,m in values.items():errors.append(dict(time=step*.005,kind=kind,branch=branch,region=region,field=field,RMS=m['rms']['value'],maximum=m['maximum']['value'],classification='DIAGNOSTIC_ONLY; provisional limits apply to temporal differences'))
        for branch,key in (('coarse','all_regional_Q1_moments'),('fine','fine_Q1')):
            for region,data in v[key].items():
                for i,value in enumerate(data['moments']):moments.append(dict(time=step*.005,branch=branch,region=region,test=i,signed_moment=value))
        for name,c in v['principal']['checks'].items():
            i=c['row'];ids=arrays['action_outer_ids'];j=int(np.flatnonzero(ids==i)[0]);rec=dict(time=step*.005,check=name,**c,physical=arrays['physical'][i].tolist(),geometry=int(arrays['geometry'][i]))
            if name!='minimized':
                rec['evaluations']={branch:{key:arrays[f'action_{name}_{branch}_{key}'][j].tolist()for key in ('dimensional','nu_eff','normalized','xi','a')}for branch in ('coarse','fine')}
            principals.append(rec)
        summaries.append(dict(step=step,pair=pair,status=v['status'],coarse=v['coarse_summary'],fine=v['fine_summary'],inventory=v['inventory'],signed_arrays=v['signed_arrays'],
            owned_rows=len(inventory['owned']),points=len(inventory['points']),quadrature_rows=inventory['quadrature_rows']))
        signed.append(dict(step=step,stage_audits=report['stage_audits'],native_hard_audit=report['native_hard_audit'],pimple=report['pimple_attempts'],
            instantaneous_spatial_residual='not retained by the ordinary per-step packet; final strict provides separately normalized alpha0 norms'))
        if step:
            completion=load(folder/f'step-{step}.completion.json');commit=folder/f'step-{step}.qualification-commit.json';cp=checkpoint_payload(verify(completion['checkpoint']).parent,include_reference=True)
            callback=native_load(verify(completion['callback']));b=completion['binding'];a=report['pimple_attempts'][0]
            timeline.append(dict(step=step,time=report['time'],dt=report['dt'],alpha=report['last_stage_alpha'],native='PASS',independent=v['physical_coarse_status'],provisional_screen=v['status'],qualification_commit=commit.exists(),
                U_solves=a['momentum_predictor_solves'],pressure_solves=a['pressure_system_solves'],k_solves=a['k_transport_solves'],epsilon_solves=a['epsilon_transport_solves'],bound_calls=a['bound_calls'],bound_activations=a['bound_activations'],material_corrections=a['additional_correction_or_restoration_solves'],
                native_compute_seconds=load(folder/f'step-{step}.ready.json')['native_compute_seconds'],inclusive_step_seconds=duration(step,'native.step'),full_audit_seconds=duration(step,'native.acceptance_audit'),physical_samples_seconds=duration(step,'native.physical_samples'),
                recovery_seconds=duration(step,'monitor.checkpoint_and_publish'),observer_wait_seconds=duration(step,'monitor.independent_wait'),verification_seconds=duration(step,'monitor.qualification_commit')))
            bindings.append(dict(step=step,time=cp['time'],dt=cp['dt'],checkpoint_file=completion['checkpoint'],contained_state_sha256=vector_hash(cp['state']),producer=completion['report'],callback=completion['callback'],callback_history_sha256=vector_hash(callback['callback_history']),
                reference_sha256=vector_hash(cp['fixed_reference']),defect_sha256=vector_hash(cp['fixed_defect']),native_history_sha256=[vector_hash(h)for h in cp['native_history']],native_history_times=cp['native_history_times'],physical_signature_sha256=hashlib.sha256(cp['signature'].encode()).hexdigest(),qualification_commit=identity(commit)if commit.exists()else None))
            storage.append(dict(step=step,recovery_bytes=sum(c['bytes']for c in (completion['checkpoint'],completion['callback'],completion['report']))+identity(folder/f'step-{step}.completion.json')['bytes'],
                comparison_bytes=p.stat().st_size+v['signed_arrays']['bytes'],inventory_bytes=v['inventory']['bytes']))
            heavy.extend(completion[k]for k in ('checkpoint','callback','report'))
            for q in (folder/f'step-{step}.completion.json',commit,folder/f'step-{step}.ack.consumed.json'):copy(q,'coarse/'+q.name)
    for name,values in (('time-pairs.csv',pairs),('field-sensitivity.csv',fields),('reference-drift.csv',errors),('all-signed-Q1-moments.csv',moments),('target-timeline.csv',timeline),('checkpoint-storage.csv',storage)):table(out/name,values)
    atomic(out/'principal-same-owner.json',principals);atomic(out/'native-independent-witness-replay.json',dict(status='PASS',records=witnesses,compiled_kernel_requests=0,independent_complex_direction_evaluations=sum(len(w['replays'])for w in witnesses)))
    atomic(out/'state-history-bindings.json',bindings);atomic(out/'audited-summaries.json',summaries)
    selected=load(root/'inputs/anchors-v2.json');coarse_cp=load(root/'science/coarse-cp0/lineage.json')
    history=[]
    for role,pin in (('fine_original',selected['original']['checkpoint']),('coarse_native',coarse_cp['coarse_checkpoint'])):
        cp=checkpoint_payload(verify(pin).parent,include_reference=True)
        history.append(dict(role=role,checkpoint_file=pin,state_sha256=vector_hash(cp['state']),fixed_reference_sha256=vector_hash(cp['fixed_reference']),fixed_defect_sha256=vector_hash(cp['fixed_defect']),eta=cp['eta'],
            signature_sha256=hashlib.sha256(cp['signature'].encode()).hexdigest(),dt=cp['dt'],time=cp['time'],step=cp['steps'],order=cp['order'],native_history_sha256=[vector_hash(h)for h in cp['native_history']],native_history_times=cp['native_history_times']))
    atomic(out/'bootstrap-native-history.json',dict(status='PASS',records=history,target_native_history_initializations=1,small_component_history_initializations=1,spatial_initializations=0,
        interpretation='native IMEX1 initialization, never copied/rescaled histories; zero explicit data may give equal derivative coefficients while dt labels/signatures legitimately differ'))
    q=out/'signed-native-stage.json';atomic(q,signed)
    with q.open('rb')as f,gzip.open(str(q)+'.gz','wb')as g:shutil.copyfileobj(f,g)
    q.unlink()
    strict=root/'science/final-strict/result.json';copy(strict,'final-strict-result.json.gz',True)
    # Reuse the existing fine terminal strict packet; no new fine native read.
    copy(OLD/'science/target/B/view/result.json','reused-fine-terminal-strict.json.gz',True)
    copy(root/'science/coarse/monitor/startup.verification-commit.json','coarse/startup.verification-commit.json')
    copy(root/'science/coarse/monitor/latest-fully-verified-checkpoint.json','coarse/latest-fully-verified-checkpoint.json')
    ops=[]
    for label,e in book['entries'].items():
        ops.append(dict(operation=label,phase=e['phase'],status=e['status'],charged_seconds=e.get('seconds',e['reserved_seconds']),reserved_seconds=e['reserved_seconds'],attempts=e.get('charged_attempts',0),peak_RSS=e.get('peak_tree_rss_bytes',0)))
        for n in ('receipt.json','request.json','attempt-counts.json'):copy(root/'operations'/label/n,'operations/'+label+'/'+n)
    table(out/'operations.csv',ops)
    for q in (root/'science/coarse-cp0/CP/State.bin',root/'science/final-strict/native.json'):
        if q.exists():heavy.append(identity(q))
    atomic(out/'heavy-manifest.json',dict(files=heavy))
    assert load(root/'science/coarse-cp0/lineage.json')['status']=='PASS'
    screen='PASS'if len(timeline)==2 and all(r['provisional_screen']=='PASS'for r in timeline)else next((r['provisional_screen']for r in timeline if r['provisional_screen']!='PASS'),'NOT_RUN')
    gates=dict(initial_spatial_equivalence='PASS',native_physical='PASS'if len(timeline)==2 else'PARTIAL',independent_physical='PASS'if len(timeline)==2 and all(v['independent']=='PASS'for v in timeline)else'PARTIAL_OR_OBSTRUCTED',provisional_sensitivity_screen=screen,
        coarse_CP0=book['bootstrap_reserved'],coarse_attempts=book.get('target_charged_attempts',0),coarse_completed=len(timeline),fine_attempts=0,
        final_strict='PASS'if strict.exists()else'NOT_RUN',focused_python_distinct_methods=18,native_component_history_initializations=1,target_history_initializations=1,
        projections_lifts_velocity_initializers=0,governing_builds='recorded separately from scientific attempts',quality_gate=.002068763923454,Q1_gate=2e-15,
        native_independent_action_replay='PASS',H11_temporal_adequacy='UNRESOLVED unchanged',production_dt='NOT_SELECTED',duration='NOT_ADMITTED')
    atomic(out/'component-gates.json',gates)
    if timeline:
        core=statistics.mean(r['native_compute_seconds']-r['full_audit_seconds']-r['physical_samples_seconds']for r in timeline)
        audit=statistics.mean(r['full_audit_seconds']+r['physical_samples_seconds']for r in timeline);wait=statistics.mean(r['observer_wait_seconds']+r['verification_seconds']for r in timeline);rec=statistics.mean(r['recovery_seconds']for r in timeline)
        n=math.ceil(HORIZON/.005);first=math.ceil(WINDOW[0]/.005);fulls=sum(s==n or s%(25 if s>=first else 50)==0 for s in range(1,n+1))
        old=load(OLD/'publication-v1/conditional-cost-storage.json');logs=(root/'science/coarse/monitor/telemetry.jsonl').stat().st_size/len(timeline)+(scopes_path.stat().st_size+(root/'science/coarse/native.log').stat().st_size)/len(timeline)
        forecast=dict(horizon=HORIZON,window=WINDOW,dt=.005,steps=n,ceiling_step_time=n*.005,recovery_stride=5,regular_stride=50,terminal_stride=25,
            recovery_delta_T=.025,scientific_delta_T=.25,terminal_delta_T=.125,scientific_checkpoints=fulls,recovery_checkpoints=math.ceil(n/5),
            measured_audited_step_online_compute_subscope_seconds=core,full_native_audit_seconds=audit,observed_comparison_wait_commit_seconds=wait,recovery_publication_seconds=rec,
            native_only_conditional_days=(n*core+fulls*audit)/86400,observed_H19_comparison_role_conditional_days=(n*core+fulls*(audit+wait)+math.ceil(n/5)*rec)/86400,
            measured_telemetry_timing_console_bytes_per_step=logs,H18_archive_bytes=old['archive_estimate_uncompressed_bytes'],H18_recurring_log_bytes=old['steps']*(old['measured_telemetry_bytes_per_step']+old['measured_full_block_timing_and_console_bytes_per_step']),
            hypothetical_archive_H18_science_packages_plus_H19_logs_bytes=fulls*old['measured_scientific_package_bytes']+13*old['measured_recovery_bytes']+n*logs,
            limitations='two fully audited coarse steps, no ordinary warm sample; conditional subscopes only. H19 wait includes fine-state sensitivity work, not hidden CPU. No future fine oracle, cadence/namespace/output/archive adequacy or production admission. Study NPZ comparisons are extra and not assumed routine duration payloads.')
        atomic(out/'conditional-cost-storage.json',forecast)
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(7,3.7),layout='constrained')
        names=list(dict.fromkeys(f['field']for f in fields));x=np.arange(len(names))
        for j,t in enumerate((.005,.01)):
            values=[max((f['utilization']for f in fields if f['time']==t and f['field']==name),default=np.nan)for name in names]
            ax.bar(x+(j-.5)*.35,values,.35,label=f'T={t:g}')
        ax.axhline(1,color='black',ls='--',label='provisional review limit');ax.set_yscale('log');ax.set(xticks=x,xticklabels=names,ylabel='Maximum screen utilization',title='Fixed 6 mm target: dt .005 versus .0025');ax.legend()
        fig.savefig(out/'sensitivity-screen.png',dpi=160);plt.close(fig)
        fig,ax=plt.subplots(figsize=(6,3.4),layout='constrained')
        for key in ('minimized','fine_direction','coarse_direction'):
            p=[r for r in principals if r['check']==key and r['time']>0];ax.plot([r['time']for r in p],[r['max_utilization']for r in p],'o-',label=key.replace('_',' '))
        ax.axhline(1,color='black',ls='--');ax.set(xlabel='Physical T',ylabel='Maximum same-owner budget utilization',title='Principal action and both minimizing directions');ax.legend(fontsize=8)
        fig.savefig(out/'principal-screen.png',dpi=160);plt.close(fig)
    atomic(out/'actual-storage.json',dict(science_bytes=size(root/'science'),source_build_runtime_bytes=sum(size(root/n)for n in ('source','build','runtime')),pilot_output_cap=2*1024**3,root=str(root),older_roots_unchanged=True))
    atomic(out/'copied-evidence.json',copies)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='publication snapshot; closed accounting sealed separately'))
    print(json.dumps(gates),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
