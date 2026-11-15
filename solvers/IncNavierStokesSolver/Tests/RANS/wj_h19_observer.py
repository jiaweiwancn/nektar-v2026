#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H19 physical-time observer; H18 replay/import behavior is unchanged."""
import argparse
from collections import OrderedDict
import os,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,pkey,RULE,add_witnesses,assess,serial
from wj_h17_observer import load_compatible,native_load,equal
from wj_h18_observer import check_current
from wj_h19_analysis import physical_pair,compare
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp

def common_points(cache,fixed,coarse,fine,fine_audit):
    required=OrderedDict((pkey(p),p)for p in fixed)
    for report in (coarse,fine_audit):
        for packet in report.get('wj_audit_packets',[]):
            if packet.get('principal_witnesses'):add_witnesses(required,packet,cache.report,cache.config)
    # Every fixed matching-outer x, each current state's actual cutoff and its
    # two physical inward neighborhoods. Exact rows/owners are never merged.
    ids=[i for i,(p,o)in enumerate(cache.base.rows[:cache.base.point_count])if o['interface_role']=='matching_outer']
    xs=sorted({float(cache.base.rows[i][0][0])for i in ids})
    for r in (coarse,fine):
        d=cache.base.evaluate(r['state'],principal=False)
        for x in xs:
            i=next(i for i in ids if cache.base.rows[i][0][0]==x)
            for offset in (-1e-6,0.,1e-6):
                p=[x,0.,cache.config.d0+float(d['rc'][i])+offset]
                assert cache.config.d0<p[2]<cache.config.d0+cache.config.H
                required.setdefault(pkey(p),p)
    return list(required.values())

def assess_pair(cache,fixed,coarse,fine,fine_audit,step,directory,stem):
    if step:physical_pair(coarse,fine,step)
    points=common_points(cache,fixed,coarse,fine,fine_audit)
    rows,initial,stats=cache.view(points)
    physical,a=assess(rows,coarse['state'],initial,step)
    fine_physical,b=assess(rows,fine['state'],initial,2*step)
    comparison,actions=compare(rows,a,b,initial)
    inventory=directory/(stem+'.inventory.json')
    atomic(inventory,dict(rule=RULE,extension_rule='both current native witness neighborhoods and each fixed matching-outer x/current-state cutoff +/-1e-6',
        points=points,owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count))
    # Large signed arrays remain outside Git; actual row order is hash bound.
    arrays={}
    for name,d in (('coarse',a),('fine',b),('initial',initial),('exact',rows.exact)):
        for field,value in d.items():arrays[name+'_'+field]=value
    arrays.update({'action_'+k:v for k,v in actions.items()})
    arrays['physical']=np.asarray([p for p,o in rows.rows]);arrays['geometry']=np.asarray([o['geometry']for p,o in rows.rows]);arrays['weights']=rows.weights
    payload=directory/(stem+'.signed-fields.npz');np.savez_compressed(payload,**arrays)
    # Useful signed station/profile and fixed-owner records remain compact.
    witnesses={v['row']for v in comparison['principal']['checks'].values()}
    for r in comparison['temporal'].values():witnesses.update(v['row']for v in r.values())
    detail=[]
    for i in sorted(witnesses):
        detail.append(dict(row=i,physical=rows.rows[i][0],owner=rows.rows[i][1],
            coarse={k:np.asarray(a[k][i]).tolist()if np.all(np.isfinite(a[k][i]))and(k!='epsilon'or rows.outer[i])else None for k in ('U','G','k','epsilon','R','pi','production','margin','angle','nu_eff','K','rc')},
            fine={k:np.asarray(b[k][i]).tolist()if np.all(np.isfinite(b[k][i]))and(k!='epsilon'or rows.outer[i])else None for k in ('U','G','k','epsilon','R','pi','production','margin','angle','nu_eff','K','rc')}))
    comparison.update(physical_coarse_status=physical['summary']['status'],physical_fine_status=fine_physical['summary']['status'],
        coarse_summary=physical['summary'],fine_summary=fine_physical['summary'],all_regional_Q1_moments=physical['all_regional_Q1_moments'],fine_Q1=fine_physical['all_regional_Q1_moments'],
        inventory=identity(inventory),signed_arrays=identity(payload),extreme_details=detail,cache=stats)
    return serial(comparison),physical,inventory

def daemon(root,intent_path,compat):
    began=time.monotonic()
    while not intent_path.exists():assert time.monotonic()-began<20;time.sleep(.025)
    intent=load(intent_path);directory=Path(intent['directory']);assert intent['observer_pid']==os.getpid()
    c=load(compat);assert c['dt_independence']['status']=='PASS'
    cache,manifest=load_compatible(compat);fixed=load(verify(intent['common_support']))['points']
    selected=load(verify(intent['fine_anchors']))
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=intent['nonce'],setup_seconds=time.monotonic()-began,cache=manifest['cache'],
        compatibility=identity(compat),immutable_row_rebuilds=0,fixed_points=len(fixed)))
    def evaluate(ready,step,stem):
        start=time.monotonic();report=native_load(verify(ready['report']));check_current(report,cache,step)
        structure=load(verify(ready['structure']));assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        if not step:
            original=selected['original'];fine=native_load(verify(original['producer']));fine_audit=fine
            equal(report['state'],fine['state'],'spatial qh0 identity, not full cross-dt artifact')
        else:
            anchor=selected['fine'][str(step*2)];fine=native_load(verify(anchor['report']));fine_audit=native_load(verify(anchor.get('audit_report',anchor['report'])))
            cp=checkpoint_payload(verify(ready['checkpoint']).parent,include_reference=True);same_cp(report,cp)
            side=native_load(verify(ready['callback']));equal(report['last_stage_history'],side['callback_history'],'coarse producing callback')
            fcp=checkpoint_payload(verify(anchor['checkpoint']).parent,include_reference=True);same_cp(fine,fcp)
            fside=native_load(verify(anchor['callback']));equal(fine['last_stage_history'],fside['callback_history'],'fine producing callback')
            for r,s in ((report,side),(fine,fside)):
                assert s['state_sha256']==vector_hash(r['state'])and s['callback_history_sha256']==vector_hash(r['last_stage_history'])
                assert s['alpha_mass']==1/r['dt']and bits(s['time'])==bits(r['time'])
        comparison,physical,inventory=assess_pair(cache,fixed,report,fine,fine_audit,step,directory,stem)
        comparison['pair']=dict(coarse_step=step,fine_step=2*step,coarse_time=report['time'],fine_time=fine['time'],coarse_time_bits=bits(report['time']),fine_time_bits=bits(fine['time']),coarse_dt=report['dt'],fine_dt=fine['dt'],coarse_alpha=200 if step else 0,fine_alpha=400 if step else 0,coarse_report=ready['report'],fine_report=selected['fine'][str(2*step)]['report']if step else selected['original']['producer'])
        temporal=directory/(stem+'.sensitivity.json');atomic(temporal,comparison)
        status=physical['summary']['status']
        if comparison['physical_fine_status']!='PASS':status='OBSTRUCTED'
        physical.update(status=status,binding=ready['binding'],inventory_sha256=identity(inventory)['sha256'],
            provisional_screen=dict(status=comparison['status'],record=identity(temporal)),measured_independent_seconds=time.monotonic()-start)
        dest=directory/(stem+'.result.json');atomic(dest,physical)
        return status,comparison['status'],dest,inventory
    while not(directory/'startup.ready.json').exists():
        if(directory/'observer-stop.json').exists():return
        time.sleep(.025)
    ready=load(directory/'startup.ready.json');status,screen,dest,inventory=evaluate(ready,0,'startup')
    atomic(directory/'startup.verified.json',dict(status=status,binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
    if status!='PASS' or screen!='PASS':return
    for step in (1,2):
        path=directory/f'step-{step}.ready.json'
        while not path.exists():
            if(directory/'observer-stop.json').exists():return
            time.sleep(.025)
        ready=load(path);b=ready['binding']
        assert b['nonce']==intent['nonce']and b['source_sha256']==intent['source_sha256']and b['observer_sha256']==intent['observer_sha256']
        status,screen,dest,inventory=evaluate(ready,step,f'step-{step}')
        cursor=directory/f'step-{step}.prepared.json'
        if status=='PASS':atomic(cursor,dict(role='PREPARED_BEFORE_ACK; physical verification only',binding=b,ready_sha256=identity(path)['sha256'],inventory_sha256=identity(inventory)['sha256']))
        continuation='REVIEW_STOP'if screen!='PASS'else'TERMINAL'if step==2 else'PROCEED'
        atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status=status,ready_sha256=identity(path)['sha256'],binding=b,result=identity(dest),inventory=identity(inventory),
            qualified_cursor=identity(cursor)if status=='PASS'else identity(dest),continuation=continuation))
        if status!='PASS' or screen!='PASS':
            from run_wj_h19 import halt
            halt(root,'screen_halt'if screen!='PASS'else'coverage_halt',dict(step=step,status=screen if screen!='PASS'else status,first_evidence=identity(dest)))
            # Stay alive until native consumes the stop/records its commit.
            break
    while not(directory/'observer-stop.json').exists():time.sleep(.025)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True);p.add_argument('--compat',type=Path,required=True);a=p.parse_args();daemon(a.root,a.intent,a.compat)
