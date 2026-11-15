#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H23 synchronous scientific checkpoint barrier using qualified numerical loops."""
import argparse,os,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,assess,serial,RULE
from wj_h17_observer import native_load,load_compatible,equal
from wj_h18_observer import check_current
from wj_h20_observer import common_points
from wj_h20_analysis import PROBES,availability,continuation,all_pair_screens,triplet_diagnostics
from wj_h19_analysis import metrics,principal_screen
from wj_h23_continuation import DT,START,END,STRIDE,absolute_step,partners,resolve,qualified,read_completed,validate_intent,deadline_ok,progress_forecast
from wj_h22_replay import recovery,parity,integer
from wj_h22_audit_roles import native_full_map_parity
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp
from run_wj_h23 import halt

def evaluate(root,cache,fixed,inputs,report,phase,step,directory,stem,forced_points=None):
    t=report['time'];reports={phase:report};audits={phase:report};bindings={}
    for z in partners(phase):
        anchor=resolve(root,inputs,z,t);reports[z],_=qualified(anchor,z,absolute_step(z,t));audits[z]=reports[z];bindings[z]=anchor
    points=forced_points if forced_points is not None else common_points(cache,fixed,reports,audits);rows,initial,stats=cache.view(points)
    data={};physical={}
    physical_start=time.monotonic()
    for level,r in reports.items():
        n=round(r['time']/r['dt']);physical[level],data[level]=assess(rows,r['state'],initial,n)
    physical_seconds=time.monotonic()-physical_start;comparison_start=time.monotonic()
    status='PASS'if all(v['summary']['status']=='PASS'for v in physical.values())else'OBSTRUCTED'
    comparison=dict(status='NOT_YET_AVAILABLE'if phase=='M'else'NOT_APPLICABLE_INITIAL_CP0',availability=availability(phase),contraction=None)
    actions={};diagnostics=None
    if phase=='C':
        m=metrics(rows,data['C'],data['M']);p,actions=principal_screen(data['C'],data['M'],rows.outer)
        exceeded=any(v[k]['status']=='EXCEEDED'for region in m.values()for v in region.values()for k in ('rms','maximum'))
        comparison=dict(status='UNRESOLVED'if p['status']=='UNRESOLVED'else'EXCEEDED'if exceeded or p['status']=='EXCEEDED'else'PASS',availability=availability(phase),
            pairs={'C-M':dict(temporal=m,principal=p)},contraction='NOT_YET_AVAILABLE')
    if phase=='F':
        comparison,actions=all_pair_screens(rows,data);diagnostics,triplet_arrays=triplet_diagnostics(rows,data);actions.update(triplet_arrays)
        comparison.update(availability=availability(phase),contraction=diagnostics)
    comparison_seconds=time.monotonic()-comparison_start
    inventory=directory/(stem+'.inventory.json')
    atomic(inventory,dict(rule=RULE,extension='current participating witnesses, all current cutoffs +/-1e-6 at fixed H x; exact fixed probes/all owners',
        points=points,owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count))
    arrays={}
    for level,d in data.items():arrays.update({level+'_'+k:v for k,v in d.items()})
    arrays.update(initial_U=initial['U'],initial_G=initial['G'],initial_k=initial['k'],initial_epsilon=initial['epsilon'],initial_R=initial['R'],initial_pi=initial['pi'])
    arrays.update({'analysis_'+k:v for k,v in actions.items()});arrays.update(physical=np.asarray([p for p,o in rows.rows]),geometry=np.asarray([o['geometry']for p,o in rows.rows]),weights=rows.weights)
    payload=directory/(stem+'.fields.npz');np.savez_compressed(payload,**arrays)
    errors={level:dict(reference=metrics(rows,d,rows.exact),drift=metrics(rows,d,initial))for level,d in data.items()}
    for level,d in data.items():
        start_report,_=qualified(inputs['anchors'][level],level,START[level])
        start_data=d if vector_hash(reports[level]['state'])==vector_hash(start_report['state']) else rows.evaluate(start_report['state'],principal=False)
        errors[level]['change_since_02']=metrics(rows,d,start_data)
    scalar_summary={}
    for level,d in data.items():
        scalar_summary[level]={}
        for region,mask in (('full',np.ones(len(rows.rows),bool)),('inner',~rows.outer),('outer',rows.outer)):
            scalar_summary[level][region]={k:dict(minimum=float(np.min(d[k][mask])),maximum=float(np.max(d[k][mask])),integral=float(rows.weights[mask]@d[k][mask]))for k in ('production','turbulent_production','reference_production','destruction','K','rc','Lambda')}
    decision=comparison['status']
    if decision=='NOT_APPLICABLE_INITIAL_CP0':decision='NOT_YET_AVAILABLE'
    verdict=continuation(status,decision)
    result=dict(schema='WJH23ScientificCheckpointV1',phase=phase,step=step,time=t,time_bits=bits(t),status=status,physical=physical,comparison=comparison,
        continuation_verdict=verdict,partner_bindings=bindings,inventory=identity(inventory),signed_arrays=identity(payload),cache=stats,
        representation=metrics(rows,initial,rows.exact),errors=errors,interfaces={level:rows.interface_history(d)for level,d in data.items()},physical_scalar_summaries=scalar_summary,
        timings=dict(participating_physical_seconds=physical_seconds,comparison_seconds=comparison_seconds),participant_count=len(data),
        scope='finite current H23 support; .02 incoming differences retained; contraction descriptive only')
    # Fixed probes during physical-only/pair phases retain every adjacent owner.
    if diagnostics is None:
        probes=[]
        for probe in PROBES:
            ids=[i for i,(p,o)in enumerate(rows.rows[:rows.point_count])if np.asarray(p,dtype='<f8').tobytes()==np.asarray(probe['physical'],dtype='<f8').tobytes()]
            assert sum(all(rows.rows[i][1][k]==probe[k]for k in ('geometry','x_side','normal_side'))for i in ids)==1,'fixed owner missing'
            for i in ids:
                values={level:{k:np.asarray(d[k][i]).tolist()if np.all(np.isfinite(d[k][i]))and(k!='epsilon'or rows.outer[i])else None for k in ('U','G','k','epsilon','R','pi','production','margin','angle','nu_eff','K','rc','destruction')}for level,d in data.items()}
                probes.append(dict(name=probe['name'],row=i,physical=rows.rows[i][0],owner=rows.rows[i][1],states=values))
        result['fixed_probes']=probes
    return serial(result),inventory

def wait_stop(directory):
    while not(directory/'observer-stop.json').exists():time.sleep(.025)

def daemon(root,intent_path):
    began=time.monotonic()
    while not intent_path.exists():assert time.monotonic()-began<20;time.sleep(.025)
    v=load(intent_path);level=v['level'];validate_intent(v,END[level]-START[level])
    directory=Path(v['directory']);assert v['observer_pid']==os.getpid()
    inputs=load(verify(v['selected_inputs']));cache,manifest=load_compatible(root/'science/target-cache-compatibility.json')
    fixed=load(verify(inputs['common_support']))['points']
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=v['nonce'],setup_seconds=time.monotonic()-began,
        fixed_points=len(fixed),base_owned_rows=cache.base.point_count,variable_limit=cache.variable_limit,
        immutable_row_rebuilds=0,history_union=False,cache=manifest['cache']))
    def assess_ready(path,step,stem):
        start=time.monotonic();ready=load(path);r=native_load(verify(ready['report']));check_current(r,cache,step)
        structure=load(verify(ready['structure']));assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        assert r['dt']==DT[level]and integer(r['steps'])==step and r['native_full_audit_status']=='PASS'
        if stem=='startup':
            anchor=inputs['anchors'][level];old,cp=qualified(anchor,level,step)
            current_cp=checkpoint_payload(verify(v['restored_checkpoint']).parent,include_reference=True);same_cp(r,current_cp)
            proof=parity(r,old,current_cp,cp);proof['full_map']=native_full_map_parity(r,old)
            atomic(directory/'startup-restored-parity.json',dict(proof,genuine_producer=anchor,role='new audit, not a producing callback'))
        else:read_completed(recovery(directory,step),level,step)
        result,inventory=evaluate(root,cache,fixed,inputs,r,level,step,directory,stem)
        result.update(binding=ready['binding'],inventory_sha256=identity(inventory)['sha256'],measured_independent_seconds=time.monotonic()-start)
        dest=directory/(stem+'.result.json');atomic(dest,result)
        return ready,result,dest,inventory
    def failure(error,ready,stem):
        inv=directory/(stem+'.failure.inventory.json');dest=directory/(stem+'.failure.json')
        atomic(inv,dict(status='FAIL_NOT_SCIENTIFIC_INVENTORY',reason=repr(error)))
        atomic(dest,dict(status='FAIL',reason=repr(error),binding=ready['binding'],inventory_sha256=identity(inv)['sha256']))
        halt(root,'integrity_halt',dict(level=level,step=integer(ready['binding']['absolute_step']),reason=repr(error),evidence=identity(dest)))
        return dest,inv
    def await_path(path):
        while not path.exists():
            if(directory/'observer-stop.json').exists():return False
            time.sleep(.025)
        return True
    path=directory/'startup.ready.json'
    if not await_path(path):return
    ready=load(path)
    try:
        ready,result,dest,inventory=assess_ready(path,START[level],'startup')
        verdict=result['continuation_verdict']
        atomic(directory/'startup.verified.json',dict(status='PASS'if verdict=='PROCEED'else'FAIL',binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
        if verdict!='PROCEED':
            halt(root,'screen_halt'if verdict=='REVIEW_STOP'else'coverage_halt',dict(level=level,step=START[level],reason=verdict,evidence=identity(dest)))
            wait_stop(directory);return
    except Exception as error:
        dest,inv=failure(error,ready,'startup')
        atomic(directory/'startup.verified.json',dict(status='FAIL',binding=ready['binding'],result=identity(dest),inventory=identity(inv)))
        wait_stop(directory);return
    for step in range(START[level]+1,END[level]+1):
        if step%STRIDE[level]and step!=END[level]:continue
        path=directory/f'step-{step}.ready.json'
        if not await_path(path):return
        ready=load(path);binding=ready['binding']
        try:
            assert integer(binding['absolute_step'])==step and binding['nonce']==v['nonce']
            assert binding['observer_sha256']==v['observer_sha256']and binding['source_sha256']==v['source_sha256']
            ready,result,dest,inv=assess_ready(path,step,f'step-{step}')
            status=result['status'];verdict=result['continuation_verdict'];cursor=directory/f'step-{step}.prepared.json'
            if status=='PASS':atomic(cursor,dict(role='PREPARED_PHYSICAL_NOT_NATIVE_COMMIT',binding=binding,
                ready_sha256=identity(path)['sha256'],inventory_sha256=identity(inv)['sha256']))
            next_step='REVIEW_STOP'if verdict=='REVIEW_STOP'else'TERMINAL'if step==END[level]else'PROCEED'
            if next_step=='PROCEED':
                costs=dict(ordinary=v['ordinary_step_forecast_seconds'],audit=v['audit_forecast_seconds'],observer=v['observer_forecast_seconds'])
                now=time.monotonic();forecast=progress_forecast(load(root/'accounting.json'),level,step,now-v['operation_start_monotonic'],**costs)
                forecast['governor_deadline_ok']=deadline_ok(v['deadline_monotonic'],now,level,step,**costs)
                atomic(directory/f'step-{step}.resource-admission.json',forecast)
                if forecast['status']!='PASS'or not forecast['governor_deadline_ok']:next_step='RESOURCE_LIMIT'
            atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status=status,ready_sha256=identity(path)['sha256'],binding=binding,
                result=identity(dest),inventory=identity(inv),qualified_cursor=identity(cursor)if status=='PASS'else identity(dest),continuation=next_step))
            if verdict!='PROCEED'or next_step=='RESOURCE_LIMIT':
                halt(root,'resource_halt'if next_step=='RESOURCE_LIMIT'else'screen_halt'if verdict=='REVIEW_STOP'else'coverage_halt',
                    dict(level=level,step=step,reason=next_step if next_step=='RESOURCE_LIMIT'else verdict,evidence=identity(dest)))
                break
        except Exception as error:
            dest,inv=failure(error,ready,f'step-{step}')
            atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status='FAIL',ready_sha256=identity(path)['sha256'],binding=binding,
                result=identity(dest),inventory=identity(inv),qualified_cursor=identity(dest),continuation='TERMINAL'))
            break
    wait_stop(directory)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True)
    a=p.parse_args();daemon(a.root,a.intent)
