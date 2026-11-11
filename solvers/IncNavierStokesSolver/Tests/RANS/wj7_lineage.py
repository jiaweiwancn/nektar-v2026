# SPDX-License-Identifier: MIT
"""Semantic WJ-6 correction and strict WJ-7 native endpoint import checks."""
import copy
from pathlib import Path
from wj_immutable import read, digest, artifacts, verify_artifacts, file_identity, write_once
from publish_wj6 import execution_stage_binding
from run_cls_epawf_picard import require, state_hash
from compare_cls_epawf_pimple import same_history
from wj3_measure import checkpoint_payload


def correct_stage_map(parent, output):
    old_path=parent/'compact-v3/accepted-stage-lineage.json'
    old=read(old_path);remaining=set(old);entries={};corrected=[]
    snapshots={}
    for policy in ('EP','Budget'):
        path=parent/'compact-v3'/(policy+'.json');snapshots[policy]=file_identity(path)
        for name,e in read(path)['execution'].items():
            report_identity=e.get('report')
            if not report_identity:continue
            verify_artifacts({e['receipt']['path']:e['receipt']})
            receipt=read(e['receipt']['path']);verify_artifacts(receipt['artifacts'])
            require(receipt['artifacts'].get(report_identity['path'])==report_identity,
                    'execution receipt does not bind producer report')
            d,stage,allocations=execution_stage_binding(report_identity)
            payload=lambda s:{k:v for k,v in s.items() if k!='report'}
            matches=[k for k in remaining if digest(payload(old[k]))==digest(payload(stage))]
            require(len(matches)==1,'stage packet not uniquely bound to actual report contents')
            key=matches[0];remaining.remove(key)
            fields={k:digest(v) for k,v in payload(stage).items()}
            callbacks=[dict(history_sha256=a['history_sha256'],incoming_sha256=a['incoming_sha256'],
                time=a['time'],alpha=a['alpha']) for a in d['pimple_attempts']]
            entry=dict(policy=policy,case=name,report=report_identity,original_key=key,
                checkpoint=artifacts(Path(d['checkpoint']).rglob('*')) if d.get('checkpoint') else {},
                stage_fields_sha256=fields,callbacks=callbacks,last_stage_time=d['last_stage_time'],
                last_stage_alpha=d['last_stage_alpha'],state_sha256=state_hash(d['state']),
                original_stage_packet_sha256=digest(old[key]),corrected=key!=report_identity['path'])
            entry['native_retained_packets_sha256']={k:digest(d[k]) for k in ('accepted_states','native_history','native_history_times',
                'stage_audits','step_audits','explicit_audits','pimple_attempts','audit_times','state_audits')}
            if entry['corrected']:
                require(allocations and old[key]['report']['path']==allocations[-1]['path'] and
                    old[key]['report']['sha256']==allocations[-1]['sha256'],
                    'incorrect producer is not the verified final arithmetic allocation')
                corrected.append(dict(original_key=key,correct_report=report_identity,policy=policy,case=name))
            else:require(old[key]['report']==report_identity,'unchanged report identity differs')
            entries[report_identity['path']]=entry
    require(not remaining and len(entries)==36 and len(corrected)==18,'unexpected WJ-6 correction coverage')
    result=dict(schema='WJ7StageLineageCorrectionV1',status='PASS',entries=entries,corrected_bindings=corrected,
        verified_entries=36,corrected_count=18,parent=artifacts([old_path,parent/'compact-v3/publication.json',
            parent/'export-v1/finalization.json',parent/'export-v1/transport.json']),sealed_snapshots=snapshots,
        code=artifacts([Path(__file__),Path(__file__).with_name('publish_wj6.py')]),
        semantic_scope='report role validated from native schema and sealed execution receipt; all original stage contents and callback contexts compared; no native rerun')
    write_once(output,result);return result


def verify_import(expected,actual,initial):
    """Native zero-step import; compare binary64 values, not printed paths."""
    require(actual['dt']==expected['dt'],'import_dt_mismatch')
    require(actual['signature']==expected['signature'],'import_signature_mismatch')
    require(actual['steps']==expected['steps'] and actual['time']==expected['time'],'import_time_step_mismatch')
    require(state_hash(actual['state'])==state_hash(expected['state']),'import_state_mismatch')
    same_history(expected,actual)
    require(state_hash(actual['fixed_reference'])==state_hash(expected['fixed_reference'])==state_hash(initial['fixed_reference']),
            'import_original_reference_mismatch')
    require(actual['initialization_lift_calls_iterations']==[0,0] and not actual['pimple_attempts'] and actual['implicit_callback_count']==0,
            'import_governing_solve_or_lift')
    require(actual['last_stage_alpha']==0,'readback must be startup, not fabricated preceding callback')
    return dict(status='PASS',state_sha256=state_hash(actual['state']),fixed_reference_sha256=state_hash(actual['fixed_reference']),
        original_primary_sha256=state_hash(initial['state']),native_time=actual['time'],absolute_steps=actual['steps'],
        preceding_completed_callback=dict(time=expected['last_stage_time'],alpha=expected['last_stage_alpha'],
            history_sha256=state_hash(expected['last_stage_history'])),zero_governing_solves=True,zero_lift=True)


def bind_producer(record,correction):
    identity=file_identity(record['report'])
    require(identity['sha256']==record['report_sha256'] and identity['path'] in correction['entries'],
            'import_producer_mismatch')
    require(correction['entries'][identity['path']]['report']==identity,'import_producer_role_mismatch')
    payload=checkpoint_payload(record['checkpoint']);report=read(identity['path'])
    require(state_hash(payload['state'])==state_hash(report['state'])==record['state_sha256'],'import_checkpoint_state_mismatch')
    require(payload['steps']==report['steps']==record['steps'] and payload['time']==report['time'],'import_checkpoint_time_mismatch')
    require(file_identity(Path(record['checkpoint'])/'State.bin')['sha256']==record['checkpoint_sha256'],'import_checkpoint_bytes_mismatch')
    return report,identity
