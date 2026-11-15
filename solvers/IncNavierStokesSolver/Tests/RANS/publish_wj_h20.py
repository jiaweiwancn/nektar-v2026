#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Seal H20's retained configuration obstruction; no scientific reevaluation."""
import argparse,gzip,hashlib,json,math,os,shutil
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits
from wj_h17_observer import native_load
from wj3_measure import checkpoint_payload
from publish_wj_h18 import table,HORIZON,WINDOW
from run_wj_h15 import size
from wj_h20_analysis import DT,ADDITIONAL,START,END,STRIDE,TIMES,absolute_step

def publish(root):
    out=root/'publication-v2';out.mkdir();book=load(root/'accounting.json');inputs=load(root/'inputs/selected.json');copied=[]
    assert book['target_levels']==['M']and book['bootstrap_reserved']==0 and book['target_reserved_attempts']==4
    assert book['target_halt']and not(root/'science/M/case.evolution.json').exists()
    assert not(root/'science/C').exists()and not(root/'science/F').exists()and not(root/'science/F-cp0').exists()
    def copy(path,name,compress=False):
        p=Path(path)
        if not p.exists():return
        q=out/name;q.parent.mkdir(parents=True,exist_ok=True)
        if compress:
            with p.open('rb')as f,gzip.open(q,'wb',compresslevel=6)as g:shutil.copyfileobj(f,g)
            with gzip.open(q,'rb')as g:assert hashlib.sha256(g.read()).hexdigest()==identity(p)['sha256']
        else:shutil.copyfile(p,q);assert identity(p)['sha256']==identity(q)['sha256']
        copied.append(dict(original=identity(p),stored=identity(q),compression='gzip lossless'if compress else'none'))
    for n in ('execution-config.json','admission-M.json','inputs/selected.json','accounting.json','source-correction-v1/manifest.json'):
        copy(root/n,'controls/'+n.replace('/','--'))
    for n in ('runtime-build-v1.json','runtime-dependencies-v1.json','runtime-qualified-v1.json','target-cache-compatibility.json','initial-spatial-lineage.json'):
        copy(root/'science'/n,'controls/'+n+'.gz',True)
    copy(verify(inputs['compiled']),'controls/reused-compiled-source-summary.json.gz',True)
    for n in ('outcome.json','intent.json','execution-diff.json','native-launch.json','observer.log','native.log','case.xml'):
        copy(root/'science/M'/n,'failed-M/'+n)
    intent=load(root/'science/M/intent.json');outcome=load(root/'science/M/outcome.json');log=(root/'science/M/native.log').read_text()
    assert intent['ack_timeout_seconds']==180 and outcome['native_exit']==1 and 'invalid bounded monitoring limits'in log
    source=root/'source/solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp';text=source.read_text()
    assert 'get<double>("ack_timeout_seconds")<=120'in text
    assert text.index('ConfigureMonitoring();')<text.index('void RANSEquilibriumWallEvolution::v_DoInitialise')
    reason=next(s.strip()for s in log.splitlines()if s.startswith('EP-AWF evolution initialisation:'))
    failure=dict(status='OBSTRUCTED_CONFIGURATION',first_native_reason=reason,requested_timeout_seconds=180,native_maximum_seconds=120,
        native_source=identity(source),native_function='ConfigureMonitoring in v_InitObject, before v_DoInitialise/restore and v_DoSolve/Advance',
        original_outcome=identity(root/'science/M/outcome.json'),original_log=identity(root/'science/M/native.log'),
        no_native_producer=True,no_new_checkpoint=True,no_startup_ready=True,
        actual_governing_attempts=0,actual_count_basis='explicit early configuration rejection; no integrator Advance reached',
        reservation=dict(maximum=4,charged=book['entries']['target-M-v1']['charged_attempts'],original_unknown_attempts=book['entries']['target-M-v1']['unknown_attempts'],
            original_receipt_unchanged=True,meaning='conservative reserved/unknown charge; never counted as accepted or observed governing attempts'),
        post_stop_correction=identity(root/'source-correction-v1/manifest.json'),governing_requalification='NOT_RUN; all target admission remains closed')
    assert not(root/'science/M/monitor/startup.ready.json').exists()
    atomic(out/'first-failure.json',failure)
    planned=[]
    for z in ('M','C','F'):
        planned.append(dict(level=z,dt=DT[z],restored_absolute_step=START[z],additional_attempt_ceiling=ADDITIONAL[z],final_absolute_step=END[z],final_time=.02,
            native_process_launches=int(z=='M'),observed_governing_attempts=0,completed=0,newly_qualified=0,reserved_attempts=4 if z=='M'else 0,
            conservatively_charged_attempt_capacity=4 if z=='M'else 0,status='OBSTRUCTED_CONFIGURATION'if z=='M'else'NOT_RUN'))
    table(out/'planned-versus-executed.csv',planned)
    pairs=[]
    for t in TIMES:
        pairs.append(dict(time=t,C_step=absolute_step('C',t),M_step=absolute_step('M',t),F_step=absolute_step('F',t),
            C_input_role='REUSED_H19'if t<=.01 else'NOT_RUN',M_input_role='REUSED_H18'if t<=.01 else'NOT_RUN',F_input_role='NOT_RUN',
            absolute_all_pair_screen='NOT_RUN',triplet_contraction='NOT_RUN',reason='first M configuration stop; no fabricated missing partners'))
    table(out/'physical-time-pairs.csv',pairs)
    table(out/'unavailable-comparisons.csv',[dict(time=t,pair=p,quantity=q,region=r,status='NOT_RUN',value='',reason='H20 required branch states absent')
        for t in TIMES for p in ('C-M','M-F','C-F')for q in ('U','G','k','R','pi','outer_epsilon','principal_minimized','principal_common_direction')for r in ('full','inner','outer')if not(q in ('outer_epsilon','principal_minimized','principal_common_direction')and r=='inner')])
    cp_records=[]
    for level,pin in (('common_original_M',inputs['original']['checkpoint']),('C_own_native_CP0',inputs['C_CP0'])):
        cp=checkpoint_payload(verify(pin).parent,include_reference=True)
        cp_records.append(dict(role=level,checkpoint=pin,dt=cp['dt'],step=cp['steps'],time=cp['time'],time_bits=bits(cp['time']),
            hashes={k:vector_hash(cp[k])for k in ('state','fixed_reference','fixed_defect')},eta=cp['eta'],signature_sha256=hashlib.sha256(cp['signature'].encode()).hexdigest(),
            native_history_sha256=[vector_hash(v)for v in cp['native_history']],native_history_times=cp['native_history_times']))
    atomic(out/'initial-and-anchor-lineage.json',dict(initial_spatial_equivalence='PASS for existing C/M; F NOT_RUN',records=cp_records,old_completed_anchors=inputs['states'],
        distinction='checkpoint file hashes are separate from contained coefficient-array hashes; M2 A producer and B restored audit remain distinct'))
    strict=root/'science/M/final-strict/result.json';s=load(strict);assert s['status']=='PASS'and s['selected_absolute_step']==4
    request=load(root/'science/M/final-strict/request.json');assert request['native_report']==inputs['states']['M']['4']['report']['path']
    copy(strict,'protected-original-M-strict.json.gz',True)
    keys=('schema','state_sha256','gauge_unknowns','mass_flux','boundary_flux_correction','alpha_mass','time','history_sha256','scope','interface')
    atomic(out/'protected-original-M-strict-compact.json',dict(status='PASS',role='new strict read of the OLD H18 M step4/T=.01 anchor; not an H20 terminal state',
        source=identity(strict),request=identity(root/'science/M/final-strict/request.json'),selected_absolute_step=4,physical_time=.01,physical_mass_dual_norms=s['physical_mass_dual_norms'],weak=s['weak'],interface={k:s['interface'][k]for k in keys if k in s['interface']}))
    copy(root/'science/M/final-strict/request.json','protected-original-M-strict-request.json')
    # Reuse only manifest-selected H19 context; never label it a new H20 result.
    oldpub=load(verify(inputs['H19_publication']));known={v['path']:v for v in oldpub['files']}
    for name in ('component-gates.json','conditional-cost-storage.json','all-signed-Q1-moments.csv'):
        pin=next(v for k,v in known.items()if k.endswith('/publication-v1/'+name));copy(verify(pin),'reused-H19/'+name)
    oldcost=load(out/'reused-H19/conditional-cost-storage.json')
    forecast=[]
    for z,dt in DT.items():
        n=math.ceil(HORIZON/dt);first=math.ceil(WINDOW[0]/dt);recovery=round(.025/dt);regular=round(.25/dt);terminal=round(.125/dt)
        science=sum(step==n or step%(terminal if step>=first else regular)==0 for step in range(1,n+1))
        forecast.append(dict(level=z,dt=dt,horizon=HORIZON,steps=n,ceiling_step_time=n*dt,recovery_stride=recovery,regular_stride=regular,terminal_stride=terminal,
            proposed_scientific_archives=science,proposed_recoveries=math.ceil(n/recovery),H20_timing_forecast='NOT_ESTABLISHED',prior_H19_comparison_role_days=oldcost['observed_H19_comparison_role_conditional_days']if z=='C'else None,
            scope='hypothetical unchanged physical cadences; no H20 governing timing or duration/archive admission'))
    atomic(out/'conditional-clock-storage-status.json',dict(records=forecast,H19_prior_only=identity(out/'reused-H19/conditional-cost-storage.json'),H20_measured_ordinary_compute='NOT_RUN',H20_storage_forecast='NOT_UPDATED; no complete trajectory payloads',
        horizon=HORIZON,terminal_window=WINDOW,output_issues='H19 approximate47.6-day/165.5GB remains conditional; production namespaces/log volume/retention and application budgets unvalidated'))
    ops=[]
    for label,e in book['entries'].items():
        ops.append(dict(operation=label,phase=e['phase'],status=e['status'],charged_seconds=e.get('seconds',e['reserved_seconds']),reserved_seconds=e['reserved_seconds'],reserved_attempts=e.get('reserved_attempts',0),charged_attempt_capacity=e.get('charged_attempts',0),unknown_attempts=e.get('unknown_attempts',False),peak_tree_RSS=e.get('peak_tree_rss_bytes',0)))
        if e['status']!='RESERVED':
            for n in ('request.json','receipt.json','attempt-counts.json','command.log','native-exec-samples.json'):copy(root/'operations'/label/n,'operations/'+label+'/'+n)
    table(out/'operations.csv',ops)
    atomic(out/'component-gates.json',dict(scope_result='OBSTRUCTED_CONFIGURATION',selected_provenance='PASS',existing_initial_C_M_spatial_identity='PASS',F_bootstrap='NOT_RUN',
        initial_focused_python='PASS',initial_methods=24,post_stop_passive_correction='PASS',corrected_methods=26,distinct_methods=26,total_method_executions=50,
        executed_native_source='b2365e52b0d334d7b61c6be3a4a3bb9029491084',executed_H20_Python='778b4b82',post_stop_corrected_Python='a9d9ab51; passive only',
        target_configuration='FAIL; ACK180 exceeds native120',new_native_physical='NOT_RUN',new_independent_physical='NOT_RUN',new_principal_quality='NOT_RUN',new_signed_Q1='NOT_RUN',
        absolute_provisional_sensitivity='NOT_RUN',triplet_contraction='NOT_RUN',protected_strict_original_M='PASS',maximum_historical_time=.01,new_scientific_attempts=0,new_CP0=0,native_builds=0,
        governor_process_launches=1,governing_attempts_observed=0,reserved_governing_capacity=4,charged_unknown_capacity=4,new_completed=0,new_qualified=0,new_final_strict_reads=1,
        H11_temporal_adequacy='UNRESOLVED unchanged',H19_provisional_screen='historical PASS unchanged',production_dt='NOT_SELECTED',duration='NOT_ADMITTED'))
    heavy=[identity(p)for p in (root/'science/M/final-strict/native.json',root/'science/M/final-strict/result.json',root/'operations/target-M-v1/command.log')]
    atomic(out/'heavy-manifest.json',dict(files=heavy,old_selected_heavy_inputs='in controls/inputs--selected.json; historical bytes unchanged'))
    atomic(out/'actual-storage.json',dict(science_bytes=size(root/'science'),source_build_runtime_bytes=sum(size(root/n)for n in ('source','source-correction-v1','build','runtime')if(root/n).exists()),pilot_output_cap=2*1024**3))
    atomic(out/'copied-evidence.json',copied)
    atomic(out/'manifest.json',dict(files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()],scope='publication snapshot; active publication operation log excluded; final accounting separately sealed'))
    print(json.dumps(failure),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();publish(a.root)
