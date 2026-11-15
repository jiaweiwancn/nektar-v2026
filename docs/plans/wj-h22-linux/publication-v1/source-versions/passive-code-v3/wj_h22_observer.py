#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H22 prefix replay/capture checks at existing scientific barriers only."""
import argparse,os,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,bits,vector_hash,assess,RULE
from wj_h17_observer import load_compatible,native_load,equal
from wj_h20_observer import common_points
from wj_h18_observer import check_current
from wj_h22_replay import DT,END,read_completed,recovery,parity,phases,flow_vectors
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp
from wj_h22_audit_roles import native_full_map_parity

def prepared_cursor(binding,ready_pin,inventory_pin,proof_pin):
    return dict(role='PREPARED_NOT_NATIVE_COMMIT',binding=binding,ready_sha256=ready_pin['sha256'],
                inventory_sha256=inventory_pin['sha256'],prefix_parity=proof_pin)

def prefix(root,directory,level,step,inputs,offsets):
    predecessor=checkpoint_payload(verify(inputs['cp0'][level]['checkpoint']).parent,include_reference=True)
    records=[]
    for n in range(1,step+1):
        anchor=recovery(directory,n);r,cp=read_completed(anchor,level,n)
        expected=inputs['oracles'][level][str(n)];old,oldcp=read_completed(expected,level,n)
        exact=parity(r,old,cp,oldcp)
        attempt=r['pimple_attempts'][0]
        assert attempt['scheduled_complete']and all(attempt[k]==1 for k in
            ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'))
        assert attempt['bound_activations']==attempt['additional_correction_or_restoration_solves']==0
        assert attempt['k_bound']['raw_sha256']==attempt['k_bound']['corrected_sha256']
        snapshot=dict(path=attempt['k_snapshot'],sha256=attempt['k_snapshot_sha256'])
        snapshot['bytes']=Path(snapshot['path']).stat().st_size
        packet=native_load(verify(snapshot));states,binding=phases(r,packet,predecessor,offsets)
        assert vector_hash(states['C'])==attempt['incoming_sha256']and vector_hash(packet['history'])==attempt['history_sha256']
        _,vectors=flow_vectors(attempt['flow_vectors'],states,r,offsets,identity(root/'science/flow-vector-execution.json'))
        matrix=identity(packet['matrix_file']);assert matrix['sha256']==packet['matrix_sha256']
        full_map=None
        if n==END[level]:
            auxiliary=inputs.get('auxiliary_audits',{}).get(level,{}).get(str(n))
            audit=native_load(verify(auxiliary['audit_report']))if auxiliary else None
            full_map=native_full_map_parity(r,old,audit)
        records.append(dict(step=n,exact=exact,current=anchor,historical=expected,phase_binding=binding,
             snapshot=snapshot,sparse_k_matrix=matrix,saved_k_matrix_operations=0,flow_vectors=vectors,
             historical_full_audit=old['native_full_audit_status'],new_full_audit=r['native_full_audit_status'],full_map_comparison=full_map))
        predecessor=cp
    return records

def daemon(root,intent_path,compat):
    start=time.monotonic()
    while not intent_path.exists():assert time.monotonic()-start<20;time.sleep(.025)
    intent=load(intent_path);directory=Path(intent['directory']);level=intent['level']
    assert intent['observer_pid']==os.getpid()
    cache,m=load_compatible(compat);inputs=load(verify(intent['selected_inputs']))
    fixed=load(verify(inputs['common_support']))['points']
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=intent['nonce'],setup_seconds=time.monotonic()-start,
        immutable_row_rebuilds=0,cache=m['cache'],fixed_points=len(fixed),role='exact H22 replay, not H20 sensitivity',compatibility=identity(compat)))
    def evaluate(ready,stem,step):
        began=time.monotonic();r=native_load(verify(ready['report']));structure=load(verify(ready['structure']))
        assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        check_current(r,cache,step)
        points=common_points(cache,fixed,{level:r},{level:r});rows,initial,stats=cache.view(points)
        result,data=assess(rows,r['state'],initial,step)
        inventory=directory/(stem+'.inventory.json')
        atomic(inventory,dict(rule=RULE,points=points,owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],
            physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count,
            extension='fixed common support, complete current native witnesses and current cutoff neighborhoods'))
        result.update(status=result['summary']['status'],binding=ready['binding'],inventory_sha256=identity(inventory)['sha256'],
            measured_independent_seconds=time.monotonic()-began,cache=stats)
        dest=directory/(stem+'.result.json');atomic(dest,result)
        return r,result,dest,inventory
    while not(directory/'startup.ready.json').exists():
        if(directory/'observer-stop.json').exists():return
        time.sleep(.025)
    ready=load(directory/'startup.ready.json');r,result,dest,inventory=evaluate(ready,'startup',0)
    same_cp(r,checkpoint_payload(verify(inputs['cp0'][level]['checkpoint']).parent,include_reference=True))
    atomic(directory/'startup.verified.json',dict(status=result['status'],binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
    if result['status']!='PASS':
        from run_wj_h22 import halt
        halt(root,'coverage_halt',dict(level=level,step=0,evidence=identity(dest)));return
    step=END[level];path=directory/f'step-{step}.ready.json'
    while not path.exists():
        if(directory/'observer-stop.json').exists():return
        time.sleep(.025)
    ready=load(path);b=ready['binding']
    assert int(b['absolute_step'])==step and b['nonce']==intent['nonce']and b['source_sha256']==intent['source_sha256']
    assert b['observer_sha256']==intent['observer_sha256']
    try:
        records=prefix(root,directory,level,step,inputs,cache.report['offsets'])
        proof=directory/'prefix-replay-capture.json';atomic(proof,dict(status='PASS',level=level,steps=records,
            scope='exact every-recovery parity and genuine selected phases/vectors before scientific ACK; nondue full audits remain NOT_SCHEDULED'))
        r,result,dest,inventory=evaluate(ready,f'step-{step}',step)
        status=result['status'];prepared=directory/f'step-{step}.prepared.json'
        if status=='PASS':atomic(prepared,prepared_cursor(b,identity(path),identity(inventory),identity(proof)))
        atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status=status,ready_sha256=identity(path)['sha256'],binding=b,
            result=identity(dest),inventory=identity(inventory),qualified_cursor=identity(prepared)if status=='PASS'else identity(dest),continuation='TERMINAL'))
        if status!='PASS':
            from run_wj_h22 import halt
            halt(root,'coverage_halt',dict(level=level,step=step,evidence=identity(dest)))
    except Exception as error:
        failed_inventory=directory/'first-replay-capture-failure.inventory.json'
        atomic(failed_inventory,dict(status='FAIL_NOT_A_SCIENTIFIC_INVENTORY',reason=repr(error)))
        failed=directory/'first-replay-capture-failure.json';atomic(failed,dict(status='FAIL',reason=repr(error),binding=b,
            inventory_sha256=identity(failed_inventory)['sha256']))
        from run_wj_h22 import halt
        halt(root,'integrity_halt',dict(level=level,step=step,evidence=identity(failed),reason=repr(error)))
        # A completed state remains unqualified; native consumes a genuine FAIL.
        atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status='FAIL',ready_sha256=identity(path)['sha256'],binding=b,
            result=identity(failed),inventory=identity(failed_inventory),qualified_cursor=identity(failed),continuation='TERMINAL'))
    while not(directory/'observer-stop.json').exists():time.sleep(.025)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True);p.add_argument('--compat',type=Path,required=True)
    a=p.parse_args();daemon(a.root,a.intent,a.compat)
