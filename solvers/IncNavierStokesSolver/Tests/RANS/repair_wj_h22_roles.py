#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only role correction and partial C/M attribution; no commit/ACK recovery."""
from pathlib import Path
import argparse,os
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from wj_h17_observer import native_load
from wj_h22_replay import read_completed,recovery
from wj_h22_audit_roles import native_full_map_parity
from wj_h22_observer import prefix
from wj_h22_phase import run

def repair(root):
    book=load(root/'accounting.json');assert book['target_levels']==['C','M']and book['target_reserved_attempts']==3
    assert book['target_halt']and book['integrity_halt']
    assert not(root/'science/F/native-launch.json').exists()
    inputs=load(root/'inputs/selected.json')
    index=load(verify(inputs['prior_indices']['20']['publication']))
    original_inputs=load(verify(next(v for v in index['files']if v['path'].endswith('/controls/inputs--selected.json'))))
    separate=original_inputs['states']['M']['2']['audit_report'];verify(separate)
    original_qualification=original_inputs['states']['M']['2']['qualified_commit'];verify(original_qualification)
    current_anchor=recovery(root/'science/M/monitor',2)
    current,cp=read_completed(current_anchor,'M',2)
    producer,_=read_completed(inputs['oracles']['M']['2'],'M',2)
    audited=native_load(verify(separate));proof=native_full_map_parity(current,producer,audited)
    cursor=root/'science/M/monitor/latest-fully-verified-checkpoint.json';before=identity(cursor)
    out=root/'science/passive-role-fix-v1';out.mkdir()
    corrected=dict(inputs,auxiliary_audits={'M':{'2':dict(audit_report=separate,original_verification=original_qualification,
        role='separate H18 B startup audit; A2 remains the completed producer')}})
    structure=load(verify(inputs['sparse_structure']));results={}
    for z,n in (('C',1),('M',2)):
        records=prefix(root,root/f'science/{z}/monitor',z,n,corrected,structure['offsets'])
        atomic(out/(z+'-retained-prefix.json'),dict(status='PASS_PASSIVE_ONLY',steps=records,native_qualification_written=False))
        results[z]=identity(out/(z+'-retained-prefix.json'))
    assert before==identity(cursor)and not(root/'science/M/monitor/step-2.qualification-commit.json').exists()
    atomic(out/'result.json',dict(status='PASS_PASSIVE_CORRECTION_ONLY',original_failed_gate=identity(root/'science/M/monitor/first-replay-capture-failure.json'),
        current=current_anchor,historical_producer=inputs['oracles']['M']['2'],separate_audit=separate,separate_verification=original_qualification,
        exact_full_map=proof,prefixes=results,qualified_cursor_unchanged=before,
        original_native_failure_unchanged=True,corrected_observer_native_qualification='NOT_RUN',F='NOT_RUN',native_work=0,
        source=[identity(Path(__file__).parent/n)for n in ('repair_wj_h22_roles.py','wj_h22_audit_roles.py','wj_h22_observer.py','wj_h22_phase.py')]))
    run(root,levels=('C','M'),out_name='partial-phase-analysis-v1')
    print('M2 separate full-map source is bitwise equal; completed callback remains A2; no ACK, commit or cursor promoted',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();repair(a.root)
