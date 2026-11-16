#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H25 physical-only checkpoint barrier; qualified numerical evaluator unchanged."""
import argparse,os,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,RULE,assess
from wj_h17_observer import native_load,load_compatible,equal
from wj_h18_observer import check_current,due
from wj_h20_observer import common_points
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp

def daemon(root,ip):
    began=time.monotonic()
    while not ip.exists():assert time.monotonic()-began<20;time.sleep(.025)
    v=load(ip);d=Path(v['directory']);assert v['observer_pid']==os.getpid()
    cache,manifest=load_compatible(verify(v['cache_compatibility']))
    fixed=load(verify(v['common_support']))['points']
    atomic(d/'observer-started.json',dict(pid=os.getpid(),nonce=v['nonce'],cache=manifest['cache'],setup_seconds=time.monotonic()-began,
        fixed_points=len(fixed),variable_limit=cache.variable_limit,history_union=False))
    def await_path(p):
        while not p.exists():
            if(d/'observer-stop.json').exists():return False
            time.sleep(.025)
        return True
    def evaluate(path,step,stem):
        start=time.monotonic();ready=load(path);r=native_load(verify(ready['report']));check_current(r,cache,step)
        structure=load(verify(ready['structure']));assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        assert r['native_full_audit_status']=='PASS'and int(r['steps'])==step
        if stem=='startup':
            cp=checkpoint_payload(verify(v['restored_checkpoint']).parent,include_reference=True);same_cp(r,cp)
            old=native_load(verify(v['completed_anchor']['report']))
            for k in ('state','fixed_reference','fixed_defect','native_history','native_history_times','time','steps','dt','signature','last_stage_history','last_stage_time','last_stage_alpha'):
                equal(r[k],old[k],'restored '+k)
            atomic(d/'startup-restored-parity.json',dict(status='PASS',anchor=v['completed_anchor'],history_role='genuine completed, not alpha0'))
        else:
            cp=checkpoint_payload(verify(ready['checkpoint']).parent,include_reference=True);same_cp(r,cp)
            side=native_load(verify(ready['callback']));equal(r['last_stage_history'],side['callback_history'],'genuine producing callback')
            assert ready['binding']['callback_history_sha256']==vector_hash(side['callback_history'])
        points=common_points(cache,fixed,{'current':r},{'current':r});rows,initial,stats=cache.view(points)
        result,data=assess(rows,r['state'],initial,step)
        inv=d/(stem+'.inventory.json')
        atomic(inv,dict(rule=RULE,points=points,owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],
            physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count,current_only=True))
        result.update(status=result['summary']['status'],binding=ready['binding'],inventory_sha256=identity(inv)['sha256'],
            cache=stats,measured_independent_seconds=time.monotonic()-start,scope='full current physical-only assessment; no temporal oracle')
        dest=d/(stem+'.result.json');atomic(dest,result)
        return ready,result,dest,inv
    def evaluate_safe(path,step,stem):
        ready=load(path)
        try:return evaluate(path,step,stem)
        except Exception as e:
            inv=d/(stem+'.failure.inventory.json');dest=d/(stem+'.failure.json')
            atomic(inv,dict(status='FAIL_NOT_PHYSICAL_INVENTORY',reason=repr(e)))
            result=dict(status='FAIL',reason=repr(e),binding=ready['binding'],inventory_sha256=identity(inv)['sha256'])
            atomic(dest,result);return ready,result,dest,inv
    path=d/'startup.ready.json'
    if not await_path(path):return
    ready,result,dest,inv=evaluate_safe(path,v['restored_absolute_step'],'startup')
    atomic(d/'startup.verified.json',dict(status=result['status'],binding=ready['binding'],result=identity(dest),inventory=identity(inv)))
    if result['status']=='PASS':
        for step in range(v['restored_absolute_step']+1,v['maximum_absolute_step']+1):
            if not due(step,v):continue
            path=d/f'step-{step}.ready.json'
            if not await_path(path):return
            ready,result,dest,inv=evaluate_safe(path,step,f'step-{step}');b=ready['binding']
            assert b['nonce']==v['nonce']and b['observer_sha256']==v['observer_sha256']and b['source_sha256']==v['source_sha256']
            cursor=d/f'step-{step}.prepared.json'
            if result['status']=='PASS':atomic(cursor,dict(role='PREPARED_PHYSICAL_NOT_NATIVE_COMMIT',binding=b,
                ready_sha256=identity(path)['sha256'],inventory_sha256=identity(inv)['sha256']))
            atomic(d/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status=result['status'],binding=b,ready_sha256=identity(path)['sha256'],
                result=identity(dest),inventory=identity(inv),qualified_cursor=identity(cursor)if result['status']=='PASS'else identity(dest),
                continuation='TERMINAL'if step==v['maximum_absolute_step']else'PROCEED'))
            if result['status']!='PASS':break
    # Native NumSteps=0/final diagnostics still run after a startup commit.
    while not(d/'observer-stop.json').exists():time.sleep(.025)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True)
    a=p.parse_args();daemon(a.root,a.intent)
