#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Identity-bound WP2 gate ledger. Missing/failed evidence never closes WP2."""
import argparse,subprocess,time
from datetime import datetime
from zoneinfo import ZoneInfo
from pathlib import Path
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha,state_hash
from continue_cls_epawf_pimple import atomic


def binding(path):return dict(path=str(path),sha256=sha(path),bytes=path.stat().st_size)


def main():
    p=argparse.ArgumentParser()
    for n in ['root','checkpoints','output']:p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();r=a.root;read=lambda name:load(r/name)
    checkpoints=[load(a.checkpoints/f'checkpoint{i}.json') for i in [1,2,3]]
    static=dict(boundary=checkpoints[0]['boundary_consistency_both_models'],spatial=checkpoints[1]['passed_both_models'],
                short_timestep_cost=checkpoints[2]['short_timestep_and_cost_passed_both'])
    unperturbed=read('unperturbed-preservation-v2.json');dense=read('dense-quarter-budget-summary.json')
    boundary=read('native-boundary-comparison-bounds.json');windows=read('unperturbed-all-pairs-windows/summary.json')
    terminal=read('terminal-checks-after-transactions.json');models={}
    for model in ['CLS','standard']:
        cases=unperturbed['models'][model]['cases'];unpert_pass=True
        for phase,d in cases.items():
            require(sha(Path(d['evidence']['path']))==d['evidence']['sha256'],'unperturbed assessment changed')
            window=windows['results'][model][phase]
            require(sha(Path(window['path']))==window['sha256'],'all-pairs window receipt changed')
            for key in ['native_all_pass','full_tensor_all_pass','homogeneity_all_pass','independent_covariance_all_pass','dense_profile_all_pass']:
                unpert_pass=unpert_pass and d[key]
            unpert_pass=unpert_pass and d['actual_time']>=1.1942387693600194 and d['final_window']['transit_scaled_linear_trend']['passed']
            unpert_pass=unpert_pass and windows['results'][model][phase]['passed']
        time_mesh=all(v['all_passed'] for v in dense['models'][model]['comparisons'].values())
        time_mesh=time_mesh and all(v['passed'] for v in boundary['models'][model]['comparisons'].values())
        transaction=r/f'transit-transaction-{model}/progress.json';tr=load(transaction)
        require(sha(transaction)==terminal['results'][model]['transaction']['sha256'],'transit transaction receipt changed')
        transaction_pass=tr['status']=='passed_transit_transaction' and terminal['results'][model]['terminal_seen_grid_audit']['checks']['passed']
        recovery=r/f'recovery-{model}-followup-T4p2/assessment/result.json';recovery_pass=False;rec=dict(status='not_completed')
        if recovery.exists():
            d=load(recovery);report=Path(d['final_report']);require(sha(report)==d['final_report_sha256'],'recovery final report changed')
            trajectory=Path(d['native_progress']);require(sha(trajectory)==d['native_progress_sha256'],'recovery trajectory snapshot changed')
            endpoint=load(report);require(state_hash(endpoint['state'])==d['state_sha256'],'recovery endpoint state differs')
            rec=dict(receipt=binding(recovery),actual_time=d['actual_time'],transit_units=d['actual_time']/1.1942387693600194,
                     state_sha256=d['state_sha256'],checkpoint=endpoint['checkpoint'],checkpoint_sha256=sha(Path(endpoint['checkpoint'])/'State.bin'),
                     signature=endpoint['signature'],reference_sha256=state_hash(endpoint['fixed_reference']),defect_sha256=state_hash(endpoint['fixed_defect']),
                     native_histories_sha256=[state_hash(s) for s in endpoint['native_history']],native_history_times=endpoint['native_history_times'],
                     final_native_and_tensor_passed=d['final_physical_passed'],final_GL28_profiles_passed=d['final_GL28_profiles_passed'],
                     final_window_passed=d['GL28_window_passed'],complete_assessment_passed=d['passed'])
            envelope=r/f'recovery-{model}-final-all-pairs/{model}-recovery.json'
            if envelope.exists():
                e=load(envelope);require(e['state_sha256'][-1]==d['state_sha256'],'recovery window endpoint differs')
                rec['all_pair_window']=dict(receipt=binding(envelope),passed=e['checks']['passed'])
                recovery_pass=d['passed'] and e['checks']['passed'] and d['actual_time']>=1.1942387693600194
            rec['status']='passed' if recovery_pass else 'failed_or_incomplete_qualification'
        else:
            progress=read(f'preservation-recovery-{model}/progress.json');rec.update(native_status=progress['status'],
                verified_report_frontier=progress['observations'][-1]['time'],failure=progress.get('failure'))
            owner=read(f'recovery-{model}-followup-T4p2/progress.json')
            rec.update(qualification_owner_status=owner['status'],qualification_failure=owner.get('failure'))
        rejection=r/f'recovery-{model}-native-rejection.json'
        if rejection.exists():
            rejected=load(rejection)
            require(all(rejected['checks'].values()),'rejected-state replay or rollback check failed')
            for path,expected in rejected['bindings'].items():
                require(sha(Path(path))==expected,'rejection replay input changed')
            require(sha(Path(rejected['parent_checkpoint']['path']))==rejected['parent_checkpoint']['sha256'],
                    'accepted parent checkpoint changed after rejection')
            rec.update(status='failed_native_covariance_recovery',rejection_receipt=binding(rejection),
                actual_time=rejected['actual_accepted_time'],first_rejected_attempt_time=rejected['first_rejected_attempt_time'],
                accepted_steps=rejected['accepted_steps'],attempted_steps=rejected['attempted_steps'],
                independent_witness=rejected['witness'],failure_stage=rejected['failure_stage'])
            recovery_pass=False
        gates=dict(**static,unperturbed_physical_transit=bool(unpert_pass),time_mesh_sensitivity=bool(time_mesh),
                   restart_rollback_and_terminal_audit=bool(transaction_pass),perturbation_recovery=bool(recovery_pass))
        models[model]=dict(gates=gates,passed=all(gates.values()),recovery=rec,
                           unperturbed_cases={k:dict(actual_time=v['actual_time'],transit_units=v['transit_units'],unique_states=v['unique_states']) for k,v in cases.items()})
    executions={}
    for progress in r.glob('*/progress.json'):
        d=load(progress)
        for q in d.get('results',{}).values():
            if isinstance(q,dict) and 'observed_physical_attempts' in q:
                executions[q['report']]=q
    counts=dict(reported_native_attempts=sum(q['observed_physical_attempts'] for q in executions.values()),
                accepted_step_executions=sum(q['accepted_step_executions'] for q in executions.values()),
                deliberate_post_assignment_rejections=sum(q.get('failure')=='deliberate post-assignment explicit rejection' for q in executions.values()),
                native_invocations=len(executions),scope='unique native report paths; reused prefixes counted once, independent replay invocations counted separately')
    used=time.time()-read('origin.json')['started_epoch'];storage_kib=int(subprocess.check_output(['du','-sk',str(r)],text=True).split()[0])
    rss=max(int(q['usage'].split()[-1]) for q in executions.values())
    runtime=read('runtime-v12/identity.json')['frozen'];require(sha(Path(runtime['path']))==runtime['sha256'],'production executable changed')
    libraries=runtime['libraries']
    for path,expected in libraries.items():
        require(sha(Path(path))==expected,'frozen linked library changed: '+path)
    budget=dict(elapsed_wall_seconds=used,wall_limit_seconds=43200,storage_GiB=storage_kib/1024**2,storage_limit_GiB=96,
                maximum_native_RSS_GiB=rss/1024**2,process_limit_GiB=24,aggregate=read('aggregate.json'))
    within_budget=used<=43200 and storage_kib<=96*1024**2 and rss<=24*1024**2
    files=['unperturbed-preservation-v2.json','dense-quarter-budget-summary.json','native-boundary-comparison-bounds.json',
           'unperturbed-all-pairs-windows/summary.json','terminal-checks-after-transactions.json','recovery-initial-history-contract.json','runtime-v12/identity.json']
    result=dict(version='WP2ClosureLedgerV2',capture_Asia_Shanghai=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),
                closed=all(v['passed'] for v in models.values()) and within_budget,models=models,budget=budget,within_declared_budget=within_budget,
                physical_executions=counts,runtime=dict(executable=runtime['path'],sha256=runtime['sha256'],linked_closure_verified=True),
                evidence={n:binding(r/n) for n in files},checkpoint_evidence={str(i):binding(a.checkpoints/f'checkpoint{i}.json') for i in [1,2,3]},
                source_sha256=sha(Path(__file__)),hill_steps=0,
                interpretation='closure is only for the selected full-height neutral open flat problem; spatial/temporal/developed hill flow and experimental agreement remain separate')
    atomic(a.output,result)


if __name__=='__main__':main()
