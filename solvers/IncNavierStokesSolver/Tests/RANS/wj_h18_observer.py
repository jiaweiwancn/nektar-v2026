#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Checkpoint barrier using the qualified evaluator and bounded current support."""
import argparse
from collections import OrderedDict
import os
from pathlib import Path
import time
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,pkey,RULE,add_witnesses,assess
from wj_h17_observer import load_compatible,native_load,equal
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp


def due(step,intent):
    stride=intent.get('terminal_stride',intent['scientific_stride'])if intent.get('terminal_start',0)and step>=intent['terminal_start']else intent['scientific_stride']
    return step==intent['maximum_absolute_step']or step%stride==0


def rotate_recoveries(records,protected,keep=3):
    """Selection only: callers may remove owned unprotected recovery payloads.

    Scientific/pending/first-failure anchors are protected by explicit identities.
    The H18 short qualification disables deletion and retains every checkpoint.
    """
    assert keep>=2
    protected=set(protected)
    old=sorted((v for v in records if v['id']not in protected and v['role']=='recovery'),key=lambda v:v['step'])
    return [v for v in old[:-keep]]


def rotate_owned(root,records,protected,keep=3):
    """Remove only explicitly eligible owned recovery payloads; keep manifests.

    Pending verification and first-failure/scientific/verified anchors must be
    protected. H18 qualification retains all short-run payloads instead.
    """
    root=Path(root).resolve();removed=[]
    for record in rotate_recoveries(records,protected,keep):
        assert record.get('verification')=='NOT_REQUIRED_SUPERSEDED','pending verification input cannot rotate'
        paths=[Path(pin['path']).resolve()for pin in record['payloads']]
        assert all(p.is_relative_to(root)and p!=root for p in paths),'rotation outside owned root'
        for pin,p in zip(record['payloads'],paths):
            verify(pin);p.unlink()
            fd=os.open(p.parent,os.O_RDONLY|os.O_DIRECTORY)
            try:os.fsync(fd)
            finally:os.close(fd)
        removed.append(record['id'])
    return removed


def parity(report,expected,full):
    old=native_load(verify(expected['report']))
    for k in ('state','fixed_reference','fixed_defect','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_alpha','last_stage_time'):
        equal(report[k],old[k],k)
    if full:equal(report['native_hard_audit'],old['native_hard_audit'],'complete native audit numerical map')


def check_current(report,cache,step):
    assert report['original_qh0_state_sha256']==cache.reference_identity['qh0']
    assert vector_hash(report['fixed_reference'])==cache.reference_identity['fixed_reference']
    assert vector_hash(report['fixed_defect'])==cache.reference_identity['fixed_defect']
    assert bits(report['time'])==bits(step*report['dt'])or abs(report['time']-step*report['dt'])<=128*2.220446049250313e-16
    if step:
        assert report['last_stage_alpha']==1/report['dt']and report['last_stage_time']==report['time']
        assert len(report['last_stage_history'])==len(report['state'])


def daemon(root,intent_path,compat):
    began=time.monotonic()
    while not intent_path.exists():
        assert time.monotonic()-began<20;time.sleep(.025)
    intent=load(intent_path);directory=Path(intent['directory'])
    assert intent['observer_pid']==os.getpid()
    # Validate immutable modules/runtime/configuration once per observer run.
    cache,manifest=load_compatible(compat)
    support=load(verify(intent['common_support']))
    fixed=support['points'];fixed_keys={pkey(p)for p in fixed}
    assert len(fixed)==len(fixed_keys)
    records=load(verify(intent['expected_manifest']))['states']
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=intent['nonce'],cache=manifest['cache'],setup_seconds=time.monotonic()-began,
        immutable_row_rebuilds=0,fixed_points=len(fixed),variable_limit=cache.variable_limit,history_union=False))

    def evaluate(report,step,stem,binding):
        start=time.monotonic();check_current(report,cache,step)
        required=OrderedDict((pkey(p),p)for p in fixed)
        if report.get('wj_audit_packets'):
            add_witnesses(required,report['wj_audit_packets'][0],cache.report,cache.config)
        rows,initial,stats=cache.view(list(required.values()))
        result,data=assess(rows,report['state'],initial,step)
        inventory=directory/(stem+'.inventory.json')
        atomic(inventory,dict(rule=RULE,points=list(required.values()),owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],
            physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count,fixed_support=intent['common_support'],current_additions=len(required)-len(fixed)))
        status=result['summary']['status']
        assert report.get('native_full_audit_status','PASS')=='PASS'
        result.update(status=status,binding=binding,inventory_sha256=identity(inventory)['sha256'],cache=stats,
            measured_independent_seconds=time.monotonic()-start,reference_identity=cache.reference_identity,
            original_qh0_sha256=cache.reference_identity['qh0'],coverage='fixed common support plus complete current native witness neighborhoods; no historical union')
        dest=directory/(stem+'.result.json');atomic(dest,result)
        # Same final sealed support/operations: numerical summaries and780 moments
        # must be exact. Additional current witnesses are reported separately.
        if str(step)in records and len(required)==len(fixed):
            pin=load(verify(records[str(step)]))
            expected=load(verify(pin['result']));oldinv=load(verify(pin['inventory']))
            if oldinv['points']==list(required.values()):
                equal(result['all_regional_Q1_moments'],expected['all_regional_Q1_moments'],'all780 signed Q1 moments')
                equal(result['summary'],expected['summary'],'same-support independent summary')
        return status,dest,inventory

    startup_path=directory/'startup.ready.json'
    checkpoint_mode=intent['policy']=='CheckpointAuditV1'
    if checkpoint_mode:
        while not startup_path.exists():
            if(directory/'observer-stop.json').exists():return
            time.sleep(.025)
        ready=load(startup_path);report=native_load(verify(ready['report']));structure=load(verify(ready['structure']))
        assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        step=report['steps'];status,dest,inventory=evaluate(report,step,'startup',ready['binding'])
        atomic(directory/'startup.verified.json',dict(status=status,binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
        if status!='PASS':return
    # Absolute monotone schedule; no archive glob or all-history seen set.
    selected=(n for n in range(intent['restored_absolute_step']+1,intent['maximum_absolute_step']+1)if not checkpoint_mode or due(n,intent))
    for step in selected:
        path=directory/f'step-{step}.ready.json'
        while not path.exists():
            if(directory/'observer-stop.json').exists():return
            time.sleep(.025)
        start=time.monotonic();ready=load(path);b=ready['binding']
        assert b['nonce']==intent['nonce']and b['source_sha256']==intent['source_sha256']and b['observer_sha256']==intent['observer_sha256']
        report=native_load(verify(ready['report']));cp=checkpoint_payload(verify(ready['checkpoint']).parent,include_reference=True)
        side=native_load(verify(ready['callback']));same_cp(report,cp)
        assert b['state_sha256']==vector_hash(report['state'])and b['callback_history_sha256']==vector_hash(side['callback_history'])
        equal(report['last_stage_history'],side['callback_history'],'genuine completed callback')
        parity(report,load(verify(records[str(step)])),True)
        status,dest,inventory=evaluate(report,step,f'step-{step}.ready',b)
        atomic(directory/f'step-{step}.parity.json',dict(status='PASS',baseline=records[str(step)],scope='bitwise state/history/callback/reference/defect/time/native audit; independent current support separately recorded'))
        cursor=directory/f'step-{step}.ready.qualified.json'
        if status=='PASS':atomic(cursor,dict(role='PREPARED_BEFORE_ACK; native commit required',binding=b,ready_sha256=identity(path)['sha256'],inventory_sha256=identity(inventory)['sha256'],checkpoint=ready['checkpoint'],callback=ready['callback'],report=ready['report']))
        continuation='TERMINAL'if step==intent['maximum_absolute_step']else'PROCEED'
        if continuation!='TERMINAL'and 'deadline_monotonic'in intent:
            forecast=1.25*((intent['maximum_absolute_step']-step)*100+intent.get('strict_reserve_seconds',120)+30)
            if intent['deadline_monotonic']-time.monotonic()<forecast:continuation='RESOURCE_LIMIT'
        atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status=status,ready_sha256=identity(path)['sha256'],binding=b,result=identity(dest),
            inventory=identity(inventory),qualified_cursor=identity(cursor)if status=='PASS'else identity(dest),continuation=continuation))
        atomic(directory/'observer-progress.json',dict(last_assessed_absolute_step=step,status=status,assessment_and_publication_seconds=time.monotonic()-start),False)
        if status!='PASS':return
    while not(directory/'observer-stop.json').exists():time.sleep(.025)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True);p.add_argument('--compat',type=Path,required=True);a=p.parse_args()
    daemon(a.root,a.intent,a.compat)
