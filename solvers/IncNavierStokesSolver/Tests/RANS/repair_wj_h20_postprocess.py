#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resolve only M's saved string/int assertion; never invoke a native program."""
import argparse,os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import atomic,identity,verify,vector_hash
from wj_h17_observer import native_load,compatibility
from wj_h20_observer import read_anchor
from study_wj_h16 import same_cp
from wj3_measure import checkpoint_payload

def repair(root):
    book=load(root/'accounting.json');assert book['target_levels']==['M']and book['target_charged_attempts']==4
    assert not any(book.get(k)for k in ('coverage_halt','screen_halt','runtime_halt','qualification_halt','post_resolution_target_halt'))
    folder=root/'science/M';directory=folder/'monitor';outcome=load(folder/'outcome.json');summary=native_load(verify(outcome['producer']))
    assert outcome['native_exit']==0 and outcome['error']is None and summary['failure']==''
    assert summary['steps']==8 and summary['time']==.02 and summary['native_completed']==summary['implicit_callbacks']==4 and summary['independently_qualified']==2
    oldcfg=load(root/'execution-config.json')
    for pin in oldcfg['files']:verify(pin)
    evidence=[identity(folder/'outcome.json'),outcome['producer'],identity(root/'operations/target-M-v1/receipt.json'),identity(root/'operations/target-M-v1/command.log')]
    log=(root/'operations/target-M-v1/command.log').read_text()
    assert "['binding']['absolute_step']==END[level]"in log and 'AssertionError'in log
    for step in (5,6,7,8):
        completion=load(directory/f'step-{step}.completion.json');r=native_load(verify(completion['report']));cp=checkpoint_payload(verify(completion['checkpoint']).parent,include_reference=True);same_cp(r,cp)
        a=r['pimple_attempts'][0]
        assert all(a[k]==1 for k in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert a['bound_activations']==a['additional_correction_or_restoration_solves']==0 and a['k_bound']['raw_sha256']==a['k_bound']['corrected_sha256']
        evidence.extend(completion[k]for k in ('checkpoint','report','callback'))
        if step in (6,8):
            commit=directory/f'step-{step}.qualification-commit.json';q=load(commit);v=load(directory/f'step-{step}.result.json')
            assert int(q['binding']['absolute_step'])==step and v['status']=='PASS'and v['comparison']['status']=='NOT_YET_AVAILABLE'
            for p in v['physical'].values():
                assert p['summary']['status']=='PASS'
                assert max(abs(x)for data in p['all_regional_Q1_moments'].values()for x in data['moments'])<=2e-15
            evidence.extend((identity(commit),identity(directory/f'step-{step}.result.json')))
    q=load(directory/'latest-fully-verified-checkpoint.json');assert int(q['binding']['absolute_step'])==8
    assert isinstance(q['binding']['absolute_step'],str),'exact reproduced passive type mismatch'
    anchor={k:q[k]for k in ('checkpoint','callback','report')};anchor['qualified_commit']=identity(directory/'step-8.qualification-commit.json')
    read_anchor(anchor,'M',.02)
    strict=load(folder/'final-strict/result.json');assert strict['status']=='PASS'and strict['selected_absolute_step']==8
    evidence.append(identity(folder/'final-strict/result.json'))
    atomic(folder/'qualification.json',dict(status='PASS',native_exit=0,native_completed=4,independently_qualified=2,attempts=4,final_absolute_step=8,time=.02,
        full_audit_calls=summary['full_audit_calls'],physical_samples_calls=summary['physical_samples_calls'],producer=outcome['producer'],
        passive_correction='integer normalization of serialized absolute_step only; no native rerun or changed scientific decision'))
    evidence.append(identity(folder/'qualification.json'))
    atomic(root/'passive-resolution-v1.json',dict(status='VERIFIED_NATIVE_PASS',original_halt=book['target_halt'],evidence=evidence,
        scope='only the completed M string/int postcondition; first failed receipt and journal halt remain unchanged; later native/physical halts cannot be waived'))
    historical=load(root/'science/target-cache-compatibility.json')
    compatibility(root,verify(historical['historical_manifest']),root/'science/target-cache-compatibility-v2.json')
    c=load(root/'science/target-cache-compatibility-v2.json');c.update(schema=historical['schema'],dt_independence=historical['dt_independence'],permitted_difference='unchanged evaluator bytes imported from versioned source-v2; passive coordinator correction only')
    atomic(root/'science/target-cache-compatibility-v2-final.json',c)
    source=Path(__file__).resolve().parents[4];original_source=root/'source'
    files=[]
    for pin in oldcfg['files']:
        p=Path(pin['path']);replacement=source/p.relative_to(original_source)if p.is_relative_to(original_source)else p
        files.append(identity(replacement))
    files.extend((identity(root/'execution-config.json'),identity(root/'passive-resolution-v1.json'),identity(root/'science/target-cache-compatibility-v2-final.json'),identity(__file__)))
    cfg=dict(oldcfg,files=files,cache_compatibility=identity(root/'science/target-cache-compatibility-v2-final.json'),historical_configuration=identity(root/'execution-config.json'),
        corrected_scope='postcondition parse and explicit saved-evidence resolution only; no observer/equation change')
    atomic(root/'execution-config-v2.json',cfg)
    atomic(root/'active-execution-config.json',dict(configuration=identity(root/'execution-config-v2.json'),resolution=identity(root/'passive-resolution-v1.json')))
    print('M verified from saved bytes; no native rerun',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();repair(a.root)
