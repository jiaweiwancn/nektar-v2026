#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact three-level H20 publication from immutable saved packets only."""
import argparse,csv,gzip,hashlib,json,math,os,shutil,statistics
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash,bits
from wj_h17_observer import native_load
from wj_h20_analysis import DT,ADDITIONAL,START,END,STRIDE,TIMES,absolute_step,FLOOR,LIMITS
from wj_h20_observer import resolve,read_anchor
from publish_wj_h19 import witness_replay
from publish_wj_h18 import HORIZON,WINDOW
from wj3_measure import checkpoint_payload
from run_wj_h20 import usage
from run_wj_h15 import size

def table(path,records):
    if not records:return
    with Path(path).open('x',newline='')as f:
        w=csv.DictWriter(f,list(records[0]),lineterminator='\n');w.writeheader();w.writerows(records)

def publish(root):
    out=root/'publication-v1';out.mkdir();book=load(root/'accounting.json');inputs=load(root/'inputs/selected.json');copied=[];heavy=[]
    def copy(path,name,compressed=False):
        p=Path(path)
        if not p.exists():return
        q=out/name;q.parent.mkdir(parents=True,exist_ok=True)
        if compressed:
            with p.open('rb')as f,gzip.open(q,'wb',compresslevel=6)as g:shutil.copyfileobj(f,g)
        else:shutil.copyfile(p,q)
        copied.append(dict(original=identity(p),stored=identity(q),compression='gzip lossless'if compressed else'none'))
    for n in ('resumption.json','execution-config.json','execution-config-v2.json','active-execution-config.json','passive-resolution-v1.json','inputs/selected.json','accounting.json'):
        copy(root/n,'controls/'+n.replace('/','--'))
    for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','runtime-qualified-v1.json','target-cache-compatibility.json','initial-spatial-lineage.json'):
        copy(root/'science'/n,'controls/'+n+'.gz',True)
    copy(root/'science/target-cache-compatibility-v2-final.json','controls/target-cache-compatibility-v2.json.gz',True)
    copy(verify(inputs['compiled']),'controls/reused-compiled-source-summary.json.gz',True)
    timeline=[];bindings=[];moments=[];audits=[];timing=[];strict=[];storage=[]
    for z in ('M','C','F'):
        folder=root/'science'/z;directory=folder/'monitor';op=root/f'operations/target-{z}-v1/timing.jsonl'
        scopes=[json.loads(s)for s in op.read_text().splitlines()]if op.exists()else[]
        def duration(step,name):return sum(s['wall_seconds']for s in scopes if s['step']==step and s['name']==name)
        for s in scopes:
            if s['name']in ('native.step','native.advance','native.acceptance_audit','native.physical_samples','monitor.checkpoint_and_publish','monitor.independent_wait','monitor.qualification_commit','startup.object','startup.initialize'):
                timing.append(dict(level=z,**{k:s.get(k)for k in ('step','time','name','wall_seconds','exclusive_wall_seconds','process_cpu_seconds')}))
        for name in ('qualification.json','outcome.json','execution-diff.json','intent.json'):
            copy(folder/name,z+'/'+name)
        for name in ('startup.verification-commit.json','latest-native-completed.json','latest-durable-checkpoint.json','latest-fully-verified-checkpoint.json'):
            copy(directory/name,z+'/'+name)
        for n in range(START[z]+1,END[z]+1):
            p=directory/f'step-{n}.completion.json'
            if not p.exists():continue
            c=load(p);r=native_load(verify(c['report']));cp=checkpoint_payload(verify(c['checkpoint']).parent,include_reference=True);a=r['pimple_attempts'][0]
            rp=directory/f'step-{n}.result.json';v=load(rp)if rp.exists()else None;due=n%STRIDE[z]==0 or n==END[z];commit=directory/f'step-{n}.qualification-commit.json'
            tele=next((json.loads(s)for s in(directory/'telemetry.jsonl').read_text().splitlines()if json.loads(s)['absolute_step']==n),{})if(directory/'telemetry.jsonl').exists()else{}
            row=dict(level=z,step=n,time=r['time'],time_bits=bits(r['time']),dt=DT[z],alpha=r['last_stage_alpha'],native_completed=True,online=r['online_check_status'],native_full=r['native_full_audit_status'],independent=v['status']if v else'NOT_SCHEDULED',
                absolute_screen=v['comparison']['status']if v else'NOT_SCHEDULED',physical_commit=commit.exists(),contraction='DESCRIPTIVE_ONLY'if z=='F'and v else'NOT_YET_AVAILABLE',
                U_solves=a['momentum_predictor_solves'],pressure_solves=a['pressure_system_solves'],k_solves=a['k_transport_solves'],epsilon_solves=a['epsilon_transport_solves'],bound_calls=a['bound_calls'],bound_activations=a['bound_activations'],material_corrections=a['additional_correction_or_restoration_solves'],raw_corrected_k_equal=a['k_bound']['raw_sha256']==a['k_bound']['corrected_sha256'],
                native_compute_seconds=tele.get('native_compute_seconds'),inclusive_step_seconds=duration(n,'native.step'),native_audit_seconds=duration(n,'native.acceptance_audit'),physical_samples_seconds=duration(n,'native.physical_samples'),recovery_seconds=duration(n,'monitor.checkpoint_and_publish'),observer_wait_seconds=duration(n,'monitor.independent_wait'),verification_seconds=duration(n,'monitor.qualification_commit'),independent_seconds=v['measured_independent_seconds']if v else None)
            timeline.append(row)
            bindings.append(dict(level=z,step=n,time=cp['time'],dt=cp['dt'],checkpoint_file=c['checkpoint'],coupled_state_sha256=vector_hash(cp['state']),producer=c['report'],callback=c['callback'],
                callback_history_sha256=c['binding']['callback_history_sha256'],reference_sha256=vector_hash(cp['fixed_reference']),defect_sha256=vector_hash(cp['fixed_defect']),original_qh0_sha256=c['binding']['original_qh0_sha256'],
                native_history_sha256=[vector_hash(h)for h in cp['native_history']],native_history_times=cp['native_history_times'],signature_sha256=hashlib.sha256(cp['signature'].encode()).hexdigest(),completion=identity(p),qualified_commit=identity(commit)if commit.exists()else None))
            heavy.extend(c[k]for k in ('checkpoint','callback','report'));copy(p,z+'/'+p.name);copy(commit,z+'/'+commit.name)
            storage.append(dict(level=z,step=n,role='scientific'if due else'recovery',checkpoint_bytes=c['checkpoint']['bytes'],callback_bytes=c['callback']['bytes'],report_bytes=c['report']['bytes'],native_recovery_total=sum(c[k]['bytes']for k in ('checkpoint','callback','report'))+p.stat().st_size,independent_payload_bytes=v['signed_arrays']['bytes']if v else 0,inventory_bytes=v['inventory']['bytes']if v else 0,independent_result_bytes=rp.stat().st_size if v else 0))
        for p in sorted(directory.glob('*.result.json')):
            v=load(p);step=v['step'];copy(p,z+'/'+p.name+'.gz',True);heavy.extend((identity(p),v['inventory'],v['signed_arrays']))
            audits.append(dict(level=z,role='restored_startup'if p.name=='startup.result.json'else'completed',step=step,time=v['time'],status=v['status'],comparison=v['comparison']['status'],continuation=v['continuation_verdict'],inventory=v['inventory'],signed_arrays=v['signed_arrays']))
            for participant,physical in v['physical'].items():
                for region,data in physical['all_regional_Q1_moments'].items():
                    for i,value in enumerate(data['moments']):moments.append(dict(gate_level=z,gate_step=step,participant=participant,time=v['time'],region=region,test=i,signed_moment=value))
        sp=folder/'final-strict/result.json'
        if sp.exists():
            v=load(sp);interface=v['interface'];keys=('schema','state_sha256','gauge_unknowns','mass_flux','boundary_flux_correction','alpha_mass','time','history_sha256','scope','interface')
            atomic(out/(z+'-terminal-strict-compact.json'),dict(status=v['status'],source=identity(sp),selected_absolute_step=v['selected_absolute_step'],physical_mass_dual_norms=v['physical_mass_dual_norms'],weak=v['weak'],interface={k:interface[k]for k in keys if k in interface}))
            strict.append(dict(level=z,status=v['status'],selected_absolute_step=v['selected_absolute_step'],physical_time=DT[z]*v['selected_absolute_step'],physical_mass_dual_norms=v['physical_mass_dual_norms'],source=identity(sp)))
            heavy.extend((identity(sp),identity(folder/'final-strict/native.json')))
    table(out/'all-step-timeline.csv',timeline);table(out/'all-signed-Q1-moments.csv',moments);table(out/'disjoint-inclusive-timings.csv',timing);table(out/'checkpoint-storage.csv',storage)
    atomic(out/'state-history-bindings.json',bindings);atomic(out/'audit-role-status.json',audits);atomic(out/'terminal-strict-residuals.json',strict)
    pairs=[];fields=[];ratios=[];errors=[];principals=[];probes=[];replays=[];maxima=[];profiles=[];boundary=[];signed=[];initial=[];physical_extrema=[];native_packets=[]
    for t in TIMES:
        n=absolute_step('F',t);path=root/f'science/F/monitor/step-{n}.result.json'
        if not path.exists():pairs.append(dict(time=t,C_step=absolute_step('C',t),M_step=absolute_step('M',t),F_step=n,status='NOT_RUN',reason='required triplet not produced'));continue
        v=load(path);comparison=v['comparison'];diag=comparison['contraction'];inv=load(verify(v['inventory']));arrays=np.load(verify(v['signed_arrays']));owned=inv['owned'];count=len(owned)
        pairs.append(dict(time=t,C_step=absolute_step('C',t),M_step=absolute_step('M',t),F_step=n,status=comparison['status'],reason='one common support for all three actual states'))
        for pair,record in comparison['pairs'].items():
            for region,values in record['temporal'].items():
                for name,m in values.items():
                    for metric in ('rms','maximum'):
                        f=m[metric];fields.append(dict(time=t,pair=pair,region=region,field=name,metric=metric,scaled_difference=f['value'],limit=f['limit'],utilization=f['utilization'],floor=f['floor'],resolved=f['resolved'],status=f['status'],row=m['row'],physical_measure=m['physical_measure']))
            prefix=pair.replace('-','');ids=arrays['analysis_'+prefix+'_outer_ids']
            for name,check in record['principal']['checks'].items():
                i=check['row'];j=int(np.flatnonzero(ids==i)[0]);rec=dict(time=t,pair=pair,check=name,**check,physical=arrays['physical'][i].tolist(),geometry=int(arrays['geometry'][i]))
                if name!='minimized':rec['actions']={z:{k:arrays[f'analysis_{prefix}_{name}_{z}_{k}'][j].tolist()for k in ('dimensional','nu_eff','normalized','xi','a')}for z in ('coarse','fine')}
                principals.append(rec)
        for group,quantities in (('field',diag['fields']),('principal',diag['principal'])):
            for name,regions in quantities.items():
                for region,values in regions.items():
                    for metric in ('RMS','maximum'):
                        ratios.append(dict(time=t,group=group,quantity=name,region=region,metric=metric,**values[metric]))
        probes.append(dict(time=t,records=diag['fixed_probes'],inventory=v['inventory'],signed_arrays=v['signed_arrays']))
        for level,kinds in v['errors'].items():
            for kind,regions in kinds.items():
                for region,values in regions.items():
                    for name,m in values.items():errors.append(dict(time=t,level=level,kind=kind,region=region,field=name,RMS=m['rms']['value'],maximum=m['maximum']['value'],scope='diagnostic; temporal limits do not gate representation/drift'))
        for region,values in v['representation'].items():
            for name,m in values.items():initial.append(dict(time=t,region=region,field=name,RMS=m['rms']['value'],maximum=m['maximum']['value'],support=v['inventory']['sha256']))
        for level in ('C','M','F'):
            anchor=resolve(root,inputs,level,t);report,audit=read_anchor(anchor,level,t)
            pairs[-1].update({level+'_serialized_time':report['time'],level+'_time_bits':bits(report['time']),level+'_alpha':report['last_stage_alpha']})
            native_packets.append(dict(time=t,level=level,genuine_producer=anchor['report'],full_audit_source=anchor.get('audit_report',anchor['report']),
                native_hard_audit=audit['native_hard_audit'],genuine_stage_audits=report['stage_audits'],frozen_pimple_attempts=report['pimple_attempts'],
                scope='native physical/root/domain/linear facts; completed-stage under genuine callback; alpha0 spatial mass-dual norms only in final strict packets'))
            proxy={key:arrays[key]for key in ('physical','geometry')};proxy.update({'coarse_'+key:arrays[level+'_'+key]for key in ('k','epsilon','G')})
            replays.append(dict(time=t,level=level,producer=anchor['report'],native_audit_source=anchor.get('audit_report',anchor['report']),records=witness_replay(audit,proxy)))
            for component in ('U','G','k','R','pi','production','divergence','K','rc','destruction'):
                data=arrays[level+'_'+component]
                maxima.append(dict(time=t,level=level,quantity=component,minimum=np.min(data,axis=0).tolist(),maximum=np.max(data,axis=0).tolist()))
            outer_ids=arrays['analysis_outer_ids'];weights=arrays['weights'];outer=np.zeros(len(weights),bool);outer[outer_ids]=True
            def location(i):return dict(row=int(i),physical=arrays['physical'][i].tolist(),geometry=int(arrays['geometry'][i]),owner=owned[i]['owner']if i<count else dict(role='physical Gauss8 interior, unique element owner',geometry=int(arrays['geometry'][i])))
            extrema={}
            for name,values,ids in (
                ('in_plane',arrays[level+'_margin'],outer_ids),('transverse',arrays[level+'_transverse'],outer_ids),
                ('covariance',np.min(arrays[level+'_covariance'],axis=1),np.arange(len(weights))),('k',arrays[level+'_k'],np.arange(len(weights))),
                ('outer_epsilon',arrays[level+'_epsilon'],outer_ids),('root_N',arrays[level+'_root'],outer_ids),('outer_nu_eff',arrays[level+'_nu_eff'],outer_ids),('nu_d',arrays[level+'_nu_d'],np.arange(len(weights)))):
                i=int(ids[np.argmin(values[ids])]);extrema[name]=dict(value=float(values[i]),**location(i))
            div=arrays[level+'_divergence'];i=int(np.argmax(abs(div)));divergence=dict(absolute_maximum=float(abs(div[i])),signed_at_maximum=float(div[i]),**location(i))
            divergence['regions']={region:dict(RMS=math.sqrt(float(weights[mask]@div[mask]**2)/weights[mask].sum()),signed_mean=float(weights[mask]@div[mask]/weights[mask].sum()))for region,mask in (('full',np.ones(len(weights),bool)),('inner',~outer),('outer',outer))}
            physical_extrema.append(dict(time=t,level=level,extrema=extrema,strong_divergence=divergence,scope='finite common inventory plus physical Gauss8; distinct from signed Q1 moments'))
        # Signed physical means/components, fetch station ranges and one-sided
        # boundary traces are derived from saved arrays, never new PDE work.
        weights=arrays['weights'];outer_ids=arrays['analysis_outer_ids'];outer=np.zeros(len(weights),bool);outer[outer_ids]=True
        u0=arrays['physical'][:,2];x0=arrays['physical'][:,0]
        for pair in ('CM','MF','CF'):
            a,b=pair;delta={name:arrays[a+'_'+name]-arrays[b+'_'+name]for name in ('U','G','k','R','pi','production','K','rc','destruction')}
            for region,mask in (('full',np.ones(len(weights),bool)),('inner',~outer),('outer',outer)):
                w=weights[mask]
                signed.append(dict(time=t,pair=a+'-'+b,region=region,components={name:dict(mean=(np.tensordot(w,d[mask],axes=1)/w.sum()).tolist(),minimum=np.min(d[mask],axis=0).tolist(),maximum=np.max(d[mask],axis=0).tolist())for name,d in delta.items()}))
            for x in sorted({o['physical'][0]for o in owned}):
                ids=[i for i,o in enumerate(owned)if o['physical'][0]==x]
                profiles.append(dict(time=t,pair=a+'-'+b,x=x,owners=len(ids),max_U=float(np.max(np.linalg.norm(delta['U'][ids],axis=1))),max_G=float(np.max(np.linalg.norm(delta['G'][ids],axis=(1,2)))),max_k=float(np.max(abs(delta['k'][ids]))),max_R=float(np.max(np.linalg.norm(delta['R'][ids],axis=(1,2)))),scope='all signed profile components/heights/owners stay in NPZ'))
            for i,o in enumerate(owned):
                p=o['physical'];owner=o['owner'];roles=[]
                if p[0]==-10:roles.append(('inlet',0,-1))
                if p[0]==15:roles.append(('outlet',0,1))
                if owner['chart'][2]==0:roles.append(('ground',2,-1))
                if p[2]==20:roles.append(('top',2,1))
                for role,axis,sign in roles:
                    boundary.append(dict(time=t,pair=a+'-'+b,role=role,physical=p,owner=owner,scalar_convention='n dot D grad f; the outward Fick flux has opposite sign; represented trace, not integrated weak load',delta_mass_flux=float(sign*delta['U'][i,axis]),
                        delta_total_traction=(sign*(arrays[a+'_T'][i,:,axis]-arrays[b+'_T'][i,:,axis])).tolist(),
                        delta_k_diffusive_flux=float(sign*(arrays[a+'_nu_d'][i]*arrays[a+'_grad_k'][i,axis]-arrays[b+'_nu_d'][i]*arrays[b+'_grad_k'][i,axis])),
                        delta_outer_epsilon_flux=float(sign*(arrays[a+'_nu_d'][i]*arrays[a+'_grad_epsilon'][i,axis]-arrays[b+'_nu_d'][i]*arrays[b+'_grad_epsilon'][i,axis])/1.3)if outer[i]else None))
        arrays.close()
    for name,data in (('time-pairs.csv',pairs),('all-pair-regional-sensitivity.csv',fields),('contraction.csv',ratios),('reference-drift.csv',errors),('representation.csv',initial),('fetch-stations.csv',profiles)):table(out/name,data)
    atomic(out/'same-owner-principal-screens.json',principals);atomic(out/'fixed-probe-triplets.json',probes);atomic(out/'native-independent-replays.json',dict(status='PASS',records=replays,kernel_process_requests=0,method='independent implicit-root complex tangent against saved compiled native action'))
    atomic(out/'physical-extrema-and-divergence.json',physical_extrema)
    for name,data in (('signed-component-means.json',signed),('exterior-trace-differences.json',boundary),('field-component-extrema.json',maxima),('native-physical-frozen-stage-packets.json',native_packets)):
        p=out/name;atomic(p,data)
        with p.open('rb')as f,gzip.open(str(p)+'.gz','wb')as g:shutil.copyfileobj(f,g)
        p.unlink()
    for n in ('lineage.json','qualification.json','execution-diff.json'):
        copy(root/'science/F-cp0'/n,'F-bootstrap/'+n+'.gz',True)
    operations=[]
    for label,e in book['entries'].items():
        operations.append(dict(label=label,phase=e['phase'],status=e['status'],seconds=e.get('seconds',e['reserved_seconds']),reserved=e['reserved_seconds'],attempts=e.get('charged_attempts',0),peak_RSS=e.get('peak_tree_rss_bytes',0)))
        if e['status']!='RESERVED':
            for n in ('request.json','receipt.json','attempt-counts.json','native-exec-samples.json'):copy(root/'operations'/label/n,'operations/'+label+'/'+n)
    table(out/'operations.csv',operations)
    full=len(timeline)==22 and len(pairs)==4 and all(v['status']=='PASS'for v in pairs)and len(strict)==3
    gates=dict(result='PASS'if full else'PARTIAL_OR_STOPPED',initial_restart_integrity='PASS'if full else'RETAINED_PARTIAL',physical='PASS'if full and all(v['status']=='PASS'for v in audits)else'RETAINED_PARTIAL',absolute_provisional_screens='PASS'if full else'RETAINED_PARTIAL',
        contraction='DESCRIPTIVE_MIXED_RESULTS; see per-quantity table',new_attempts=sum(v.get('charged_attempts',0)for v in book['entries'].values()if v.get('kind')=='target'),new_completed=len(timeline),new_CP0=book['bootstrap_reserved'],new_strict_reads=len(strict),
        historical_configuration_launches=1,historical_unknown_charged_capacity=4,historical_strict_reads=1,spatial_initializations=0,native_builds=0,bound_activations=sum(v['bound_activations']for v in timeline),
        total_charged_so_far=usage(book),native_source='b2365e52b0d334d7b61c6be3a4a3bb9029491084',production_dt='NOT_SELECTED',duration='NOT_ADMITTED',H11_temporal_adequacy='UNRESOLVED historical scope')
    atomic(out/'component-gates.json',gates)
    # Conditional cost at unchanged PHYSICAL monitoring intervals; passive
    # triplet comparison is separated from a production monitor's cost.
    cost=[]
    monitor_rows=[r for r in timeline if r['level']=='M'and r['native_full']=='PASS']
    physical_wait=statistics.mean(r['observer_wait_seconds']+r['verification_seconds']for r in monitor_rows)if monitor_rows else None
    physical_independent=statistics.mean(r['independent_seconds']for r in monitor_rows)if monitor_rows else None
    physical_packages=[v for v in storage if v['level']=='M'and v['role']=='scientific']
    factor_path=root/'science/M/monitor/factors'
    factor_bytes=size(factor_path)if factor_path.exists()else 0
    package_bytes=statistics.mean(v['native_recovery_total']+v['independent_payload_bytes']+v['inventory_bytes']+v['independent_result_bytes']+factor_bytes for v in physical_packages)if physical_packages else None
    for z in ('C','M','F'):
        rows=[r for r in timeline if r['level']==z];ordinary=[r['native_compute_seconds']for r in rows if r['native_full']=='NOT_SCHEDULED'];audited=[r for r in rows if r['native_full']=='PASS']
        base=statistics.mean(ordinary)if ordinary else statistics.mean(r['native_compute_seconds']-r['native_audit_seconds']-r['physical_samples_seconds']for r in audited)if audited else None
        audit=statistics.mean(r['native_audit_seconds']+r['physical_samples_seconds']for r in audited)if audited else None
        rec=statistics.mean(r['recovery_seconds']for r in rows)if rows else None
        n=math.ceil(HORIZON/DT[z]);first=math.ceil(WINDOW[0]/DT[z]);rs,ss,ts=(round(v/DT[z])for v in (.025,.25,.125));ns=sum(s==n or s%(ts if s>=first else ss)==0 for s in range(1,n+1))
        log_paths=(root/f'science/{z}/monitor/telemetry.jsonl',root/f'science/{z}/native.log',root/f'operations/target-{z}-v1/timing.jsonl')
        log_bytes=sum(p.stat().st_size for p in log_paths if p.exists())/len(rows)if rows else None
        rec_bytes=statistics.mean(v['native_recovery_total']for v in storage if v['level']==z)if rows else None
        cost.append(dict(level=z,dt=DT[z],horizon=HORIZON,terminal_window=WINDOW,steps=n,ceiling_step_time=n*DT[z],recovery_stride=rs,regular_stride=ss,terminal_stride=ts,scientific_checkpoints=ns,
            ordinary_samples=len(ordinary),ordinary_or_subtracted_online_seconds=base,native_audit_seconds=audit,recovery_seconds=rec,native_only_conditional_days=(n*base+ns*audit+math.ceil(n/rs)*rec)/86400 if base is not None else None,
            physical_only_monitor_wait_commit_seconds_from_M=physical_wait,physical_only_independent_seconds_from_M=physical_independent,
            conditional_loop_days_without_future_temporal_oracle=(n*base+ns*(audit+physical_wait)+math.ceil(n/rs)*rec)/86400 if base is not None and physical_wait is not None else None,
            measured_log_bytes_per_step=log_bytes,measured_native_recovery_bytes=rec_bytes,measured_M_physical_science_package_bytes=package_bytes,
            hypothetical_archive_measured_encodings_bytes=ns*package_bytes+13*rec_bytes+n*log_bytes if package_bytes is not None and rec_bytes is not None else None,
            limitations='short initial-relaxation sample; no accepted production/error/cadence/archive budget; observed comparison work is not ordinary production monitoring; storage uses measured native JSON/lossless NPZ and full logs, not an uncompressed or optimized production format'))
    atomic(out/'conditional-cost.json',cost)
    atomic(out/'heavy-manifest.json',dict(files=heavy))
    atomic(out/'actual-storage.json',dict(science_bytes=size(root/'science'),all_new_output_before_archive=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if not p.name.startswith(('source','build','runtime'))),
        source_build_runtime_bytes=sum(size(p)if p.is_dir()else p.stat().st_size for p in root.iterdir()if p.name.startswith(('source','build','runtime'))),prior_output_bytes=25259673,pilot_combined_output_cap=2*1024**3))
    atomic(out/'copied-evidence.json',copied)
    if probes:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
        for quantity,ax in (('minimized',axes[0]),('theta_F',axes[1])):
            rows=[(v['time'],next(p for p in v['records']if p['name']=='outer_principal'and p['primary'])['principal'][quantity])for v in probes]
            ax.plot([t for t,v in rows],[v['signed_E1']for t,v in rows],'o-',label='C - M');ax.plot([t for t,v in rows],[v['signed_E2']for t,v in rows],'s-',label='M - F');ax.axhline(0,color='black',lw=.7)
            ax.set(xlabel='Physical T [H_h/U∞]',ylabel='Signed normalized action difference',title='Minimized'if quantity=='minimized'else'All states at F minimizing direction');ax.legend()
        fig.suptitle('Outlet x=15, geometry251 upper: physical height10.243007mm');fig.savefig(out/'signed-principal-differences.png',dpi=160);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
        for pair in ('C-M','M-F','C-F'):
            ts=sorted({r['time']for r in fields});axes[0].plot(ts,[max(r['utilization']for r in fields if r['time']==t and r['pair']==pair)for t in ts],'o-',label=pair)
        axes[0].axhline(1,color='black',ls='--');axes[0].set(xlabel='Physical T',ylabel='Largest field screen utilization',title='All regions, RMS and finite-inventory maximum');axes[0].legend()
        for field in ('U','G','R','pi'):
            rr=[v for v in ratios if v['group']=='field'and v['quantity']==field and v['region']=='full'and v['metric']=='RMS']
            axes[1].plot([r['time']for r in rr],[r['rho']if r['rho']is not None else np.nan for r in rr],'o-',label=field)
        axes[1].axhline(1,color='black',ls='--');axes[1].set(xlabel='Physical T',ylabel='D2 / D1, same-support physical RMS',title='Resolved contraction diagnostics');axes[1].legend()
        fig.savefig(out/'screens-and-field-contraction.png',dpi=160);plt.close(fig)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='publication snapshot; closed accounting sealed separately'))
    print(json.dumps(gates),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
