#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""The existing H22 assessment with explicit retained-startup and missing-F roles."""
import argparse,os,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,assess,RULE
from wj_h17_observer import native_load,load_compatible
from wj_h18_observer import check_current
from wj_h20_observer import common_points
from wj_h22_observer import prefix,prepared_cursor
from wj_h22_completion import validate_role,restored_parity,capture_descriptor
from wj_h22_replay import recovery
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp
from run_wj_h22_completion import halt

def wait_stop(directory):
    # Zero NumSteps still performs final native budgets/residuals after commit.
    while not(directory/'observer-stop.json').exists():time.sleep(.025)

def daemon(root,intent_path):
    began=time.monotonic()
    while not intent_path.exists():assert time.monotonic()-began<20;time.sleep(.025)
    intent=load(intent_path);role=intent['completion_role'];is_m=role=='RetainedM2StartupV1'
    validate_role(intent,0 if is_m else 4)
    directory=Path(intent['directory']);assert intent['observer_pid']==os.getpid()
    inputs=load(verify(intent['selected_inputs']));level=intent['level']
    cache,m=load_compatible(root/'science/target-cache-compatibility.json')
    fixed=load(verify(inputs['common_support']))['points']
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=intent['nonce'],role=role,
        setup_seconds=time.monotonic()-began,immutable_row_rebuilds=0,cache=m['cache'],fixed_points=len(fixed)))

    def evaluate(ready,stem,step):
        start=time.monotonic();r=native_load(verify(ready['report']));structure=load(verify(ready['structure']))
        assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        check_current(r,cache,step);assert r['native_full_audit_status']=='PASS'
        points=common_points(cache,fixed,{level:r},{level:r});rows,initial,stats=cache.view(points)
        result,data=assess(rows,r['state'],initial,step)
        inventory=directory/(stem+'.inventory.json')
        atomic(inventory,dict(rule=RULE,points=points,owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],
            physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count,
            extension='unchanged H22 fixed support plus current native witnesses and cutoff neighborhoods'))
        result.update(status=result['summary']['status'],binding=ready['binding'],inventory_sha256=identity(inventory)['sha256'],
            measured_independent_seconds=time.monotonic()-start,cache=stats,actual_assessed_absolute_step=step)
        dest=directory/(stem+'.result.json');atomic(dest,result);return r,result,dest,inventory

    def failed(error,binding,stem):
        inv=directory/(stem+'.failure.inventory.json');dest=directory/(stem+'.failure.json')
        atomic(inv,dict(status='FAIL_NOT_A_SCIENTIFIC_INVENTORY',reason=repr(error)))
        atomic(dest,dict(status='FAIL',reason=repr(error),binding=binding,inventory_sha256=identity(inv)['sha256']))
        halt(root,'integrity_halt',dict(role=role,evidence=identity(dest),reason=repr(error)))
        return dest,inv

    while not(directory/'startup.ready.json').exists():
        if(directory/'observer-stop.json').exists():return
        time.sleep(.025)
    ready=load(directory/'startup.ready.json');step=2 if is_m else 0
    try:
        current=native_load(verify(ready['report']))
        if is_m:
            desc=inputs['retained_attempts']['M']['2']['anchor'];cp=checkpoint_payload(verify(desc['checkpoint']).parent,include_reference=True)
            proof=restored_parity(current,desc,inputs['oracles']['M']['2'],inputs['auxiliary_audits']['M']['2'],cp)
            atomic(directory/'startup-replay-parity.json',proof)
        else:same_cp(current,checkpoint_payload(verify(inputs['cp0']['F']['checkpoint']).parent,include_reference=True))
        _,result,dest,inventory=evaluate(ready,'startup',step)
        atomic(directory/'startup.verified.json',dict(status=result['status'],binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
        if result['status']!='PASS':halt(root,'coverage_halt',dict(role=role,step=step,evidence=identity(dest)))
    except Exception as error:
        dest,inventory=failed(error,ready['binding'],'startup')
        atomic(directory/'startup.verified.json',dict(status='FAIL',binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
        wait_stop(directory);return
    if is_m:
        while not(directory/'startup.verification-commit.json').exists():
            if(directory/'observer-stop.json').exists():return
            time.sleep(.025)
        q=load(directory/'startup.verification-commit.json')
        assert q['schema']=='WJStartupVerificationCommitV1'and q['completion']==intent['restart_completion']
        atomic(directory/'startup-observer-complete.json',dict(status='PASS',commit=identity(directory/'startup.verification-commit.json'),
            role='new verification only; original completion/callback remain producer',awaiting_native_exit=True))
        wait_stop(directory);return
    if result['status']!='PASS':wait_stop(directory);return
    path=directory/'step-4.ready.json'
    while not path.exists():
        if(directory/'observer-stop.json').exists():return
        time.sleep(.025)
    ready=load(path);b=ready['binding']
    try:
        assert int(b['absolute_step'])==4 and b['nonce']==intent['nonce']and b['observer_sha256']==intent['observer_sha256']and b['source_sha256']==intent['source_sha256']
        descriptors={str(n):capture_descriptor(recovery(directory,n),identity(root/'science/flow-vector-execution.json'))for n in range(1,5)}
        structure=load(verify(ready['structure']))
        records=prefix(root,directory,'F',4,inputs,structure['offsets'],descriptors)
        proof=directory/'prefix-replay-capture.json';atomic(proof,dict(status='PASS',steps=records,source_descriptors=descriptors,
            coverage='all four new F attempts, original H20 own producing completions'))
        _,result,dest,inventory=evaluate(ready,'step-4',4)
        status=result['status'];prepared=directory/'step-4.prepared.json'
        if status=='PASS':atomic(prepared,prepared_cursor(b,identity(path),identity(inventory),identity(proof)))
        atomic(directory/'step-4.ack.json',dict(schema='WJObserverAckV1',status=status,ready_sha256=identity(path)['sha256'],binding=b,
            result=identity(dest),inventory=identity(inventory),qualified_cursor=identity(prepared)if status=='PASS'else identity(dest),continuation='TERMINAL'))
        if status!='PASS':halt(root,'coverage_halt',dict(role=role,step=4,evidence=identity(dest)))
    except Exception as error:
        dest,inventory=failed(error,b,'step-4')
        atomic(directory/'step-4.ack.json',dict(schema='WJObserverAckV1',status='FAIL',ready_sha256=identity(path)['sha256'],binding=b,
            result=identity(dest),inventory=identity(inventory),qualified_cursor=identity(dest),continuation='TERMINAL'))
    wait_stop(directory)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True);a=p.parse_args();daemon(a.root,a.intent)
