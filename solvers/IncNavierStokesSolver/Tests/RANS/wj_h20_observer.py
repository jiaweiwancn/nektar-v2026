#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Phase-aware H20 checkpoint barrier; no future-oracle waits."""
import argparse,os,time
from pathlib import Path
from collections import OrderedDict
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,pkey,add_witnesses,assess,serial,RULE
from wj_h17_observer import load_compatible,native_load,equal
from wj_h18_observer import check_current
from wj_h19_analysis import metrics,principal_screen
from wj_h20_analysis import DT,STRIDE,END,PROBES,absolute_step,selected,availability,continuation,all_pair_screens,triplet_diagnostics
from wj3_measure import checkpoint_payload
from study_wj_h16 import same_cp

def read_anchor(anchor,level,t,qualified=True):
    r=native_load(verify(anchor['report']));cp=checkpoint_payload(verify(anchor['checkpoint']).parent,include_reference=True);same_cp(r,cp)
    c=selected(r,level,t);side=native_load(verify(anchor['callback']))
    assert side['state_sha256']==vector_hash(r['state'])and side['callback_history_sha256']==vector_hash(c['history'])
    assert bits(side['time'])==bits(c['time'])and side['alpha_mass']==c['alpha']
    equal(side['callback_history'],c['history'],'genuine own callback, not cross-dt history')
    if qualified:
        q=load(verify(anchor['qualified_commit']));assert q['schema']in ('WJQualificationCommitV1','WJStartupVerificationCommitV1')
        assert int(q['binding']['absolute_step'])==c['step']and q['binding']['state_sha256']==vector_hash(r['state'])
        assert q['binding']['callback_history_sha256']==vector_hash(c['history'])
    audit=native_load(verify(anchor.get('audit_report',anchor['report'])))
    assert vector_hash(audit['state'])==vector_hash(r['state'])and audit['time']==r['time']
    return r,audit

def resolve(root,inputs,level,t):
    step=absolute_step(level,t)
    if str(step)in inputs['states'].get(level,{}):return inputs['states'][level][str(step)]
    directory=Path(root)/'science'/level/'monitor'
    commit=directory/f'step-{step}.qualification-commit.json'
    assert commit.exists(),'required computed partner absent: '+level+str(step)
    q=load(commit);assert q['schema']=='WJQualificationCommitV1'
    return dict(**{k:q[k]for k in ('checkpoint','callback','report')},qualified_commit=identity(commit),audit_report=q['report'],completion=identity(directory/f'step-{step}.completion.json'))

def common_points(cache,fixed,reports,audits):
    required=OrderedDict((pkey(p),p)for p in fixed)
    for p in PROBES:required.setdefault(pkey(p['physical']),p['physical'])
    for r in audits.values():
        for packet in r.get('wj_audit_packets',[]):
            if packet.get('principal_witnesses'):add_witnesses(required,packet,cache.report,cache.config)
    indices={}
    for i,(p,o)in enumerate(cache.base.rows[:cache.base.point_count]):
        if o['interface_role']=='matching_outer':indices.setdefault(float(p[0]),i)
    for r in reports.values():
        d=cache.base.evaluate(r['state'],principal=False)
        for x,i in sorted(indices.items()):
            for offset in (-1e-6,0.,1e-6):
                p=[x,0.,cache.config.d0+float(d['rc'][i])+offset]
                assert cache.config.d0<p[2]<cache.config.d0+cache.config.H
                required.setdefault(pkey(p),p)
    return list(required.values())

def evaluate(root,cache,fixed,inputs,report,phase,step,directory,stem):
    t=report['time'];reports={phase:report};audits={phase:report};bindings={}
    if phase=='C':
        anchor=resolve(root,inputs,'M',t);reports['M'],audits['M']=read_anchor(anchor,'M',t);bindings['M']=anchor
    elif phase=='F'and step:
        for level in ('C','M'):
            anchor=resolve(root,inputs,level,t);reports[level],audits[level]=read_anchor(anchor,level,t);bindings[level]=anchor
    points=common_points(cache,fixed,reports,audits);rows,initial,stats=cache.view(points)
    data={};physical={}
    for level,r in reports.items():
        n=round(r['time']/r['dt']);physical[level],data[level]=assess(rows,r['state'],initial,n)
    status='PASS'if all(v['summary']['status']=='PASS'for v in physical.values())else'OBSTRUCTED'
    comparison=dict(status='NOT_YET_AVAILABLE'if phase=='M'else'NOT_APPLICABLE_INITIAL_CP0',availability=availability(phase),contraction=None)
    actions={};diagnostics=None
    if phase=='C':
        m=metrics(rows,data['C'],data['M']);p,actions=principal_screen(data['C'],data['M'],rows.outer)
        exceeded=any(v[k]['status']=='EXCEEDED'for region in m.values()for v in region.values()for k in ('rms','maximum'))
        comparison=dict(status='UNRESOLVED'if p['status']=='UNRESOLVED'else'EXCEEDED'if exceeded or p['status']=='EXCEEDED'else'PASS',availability=availability(phase),
            pairs={'C-M':dict(temporal=m,principal=p)},contraction='NOT_YET_AVAILABLE')
    if phase=='F'and step:
        comparison,actions=all_pair_screens(rows,data);diagnostics,triplet_arrays=triplet_diagnostics(rows,data);actions.update(triplet_arrays)
        comparison.update(availability=availability(phase),contraction=diagnostics)
    inventory=directory/(stem+'.inventory.json')
    atomic(inventory,dict(rule=RULE,extension='current participating witnesses, all current cutoffs +/-1e-6 at fixed H x; exact fixed probes/all owners',
        points=points,owned=[dict(physical=p,owner=o)for p,o in rows.rows[:rows.point_count]],physical_rows=len(rows.rows),quadrature_rows=len(rows.rows)-rows.point_count))
    arrays={}
    for level,d in data.items():arrays.update({level+'_'+k:v for k,v in d.items()})
    arrays.update(initial_U=initial['U'],initial_G=initial['G'],initial_k=initial['k'],initial_epsilon=initial['epsilon'],initial_R=initial['R'],initial_pi=initial['pi'])
    arrays.update({'analysis_'+k:v for k,v in actions.items()});arrays.update(physical=np.asarray([p for p,o in rows.rows]),geometry=np.asarray([o['geometry']for p,o in rows.rows]),weights=rows.weights)
    payload=directory/(stem+'.fields.npz');np.savez_compressed(payload,**arrays)
    errors={level:dict(reference=metrics(rows,d,rows.exact),drift=metrics(rows,d,initial))for level,d in data.items()}
    scalar_summary={}
    for level,d in data.items():
        scalar_summary[level]={}
        for region,mask in (('full',np.ones(len(rows.rows),bool)),('inner',~rows.outer),('outer',rows.outer)):
            scalar_summary[level][region]={k:dict(minimum=float(np.min(d[k][mask])),maximum=float(np.max(d[k][mask])),integral=float(rows.weights[mask]@d[k][mask]))for k in ('production','turbulent_production','reference_production','destruction','K','rc','Lambda')}
    decision=comparison['status']
    if decision=='NOT_APPLICABLE_INITIAL_CP0':decision='NOT_YET_AVAILABLE'
    verdict=continuation(status,decision)
    result=dict(schema='WJH20ScientificCheckpointV1',phase=phase,step=step,time=t,time_bits=bits(t),status=status,physical=physical,comparison=comparison,
        continuation_verdict=verdict,partner_bindings=bindings,inventory=identity(inventory),signed_arrays=identity(payload),cache=stats,
        representation=metrics(rows,initial,rows.exact),errors=errors,interfaces={level:rows.interface_history(d)for level,d in data.items()},physical_scalar_summaries=scalar_summary,
        scope='finite current support; planned comparison absence explicit; contraction alone never an admission gate')
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

def daemon(root,intent_path,compat):
    start=time.monotonic()
    while not intent_path.exists():assert time.monotonic()-start<20;time.sleep(.025)
    intent=load(intent_path);directory=Path(intent['directory']);phase=intent['level'];assert intent['observer_pid']==os.getpid()
    cache,manifest=load_compatible(compat);inputs=load(verify(intent['selected_inputs']));fixed=load(verify(inputs['common_support']))['points']
    atomic(directory/'observer-started.json',dict(pid=os.getpid(),nonce=intent['nonce'],setup_seconds=time.monotonic()-start,immutable_row_rebuilds=0,cache=manifest['cache'],compatibility=identity(compat)))
    def assess_ready(path,step,stem):
        began=time.monotonic();ready=load(path);r=native_load(verify(ready['report']));check_current(r,cache,step)
        structure=load(verify(ready['structure']));assert structure['fields']==cache.report['fields']and structure['offsets']==cache.report['offsets']
        assert r['dt']==DT[phase]and r['steps']==step
        if step and 'callback'in ready:
            read_anchor(ready,phase,r['time'],False)
        result,inventory=evaluate(root,cache,fixed,inputs,r,phase,step,directory,stem)
        result.update(binding=ready['binding'],inventory_sha256=identity(inventory)['sha256'],measured_independent_seconds=time.monotonic()-began)
        dest=directory/(stem+'.result.json');atomic(dest,result)
        return ready,result,dest,inventory
    while not(directory/'startup.ready.json').exists():
        if(directory/'observer-stop.json').exists():return
        time.sleep(.025)
    ready=load(directory/'startup.ready.json');r=native_load(verify(ready['report']));step=r['steps']
    # Restore proof preserves original accepted producer/callback roles separately.
    cp=checkpoint_payload(verify(intent['restored_checkpoint']).parent,include_reference=True);same_cp(r,cp)
    ready,result,dest,inventory=assess_ready(directory/'startup.ready.json',step,'startup')
    verdict=result['continuation_verdict']
    atomic(directory/'startup.verified.json',dict(status=result['status']if verdict=='PROCEED'else'REVIEW_STOP'if verdict=='REVIEW_STOP'else'FAIL',binding=ready['binding'],result=identity(dest),inventory=identity(inventory)))
    if verdict!='PROCEED':
        from run_wj_h20 import halt
        halt(root,'screen_halt'if verdict=='REVIEW_STOP'else'coverage_halt',dict(phase=phase,step=step,reason=verdict,evidence=identity(dest)));return
    for step in range(intent['restored_absolute_step']+1,END[phase]+1):
        if step%STRIDE[phase]and step!=END[phase]:continue
        path=directory/f'step-{step}.ready.json'
        while not path.exists():
            if(directory/'observer-stop.json').exists():return
            time.sleep(.025)
        ready,result,dest,inventory=assess_ready(path,step,f'step-{step}')
        b=ready['binding'];assert b['nonce']==intent['nonce']and b['observer_sha256']==intent['observer_sha256']and b['source_sha256']==intent['source_sha256']
        cursor=directory/f'step-{step}.prepared.json'
        if result['status']=='PASS':atomic(cursor,dict(role='PREPARED_PHYSICAL_VERIFICATION',binding=b,ready_sha256=identity(path)['sha256'],inventory_sha256=identity(inventory)['sha256']))
        verdict=result['continuation_verdict'];next_step='REVIEW_STOP'if verdict=='REVIEW_STOP'else'TERMINAL'if step==END[phase]else'PROCEED'
        if next_step=='PROCEED':
            remaining=END[phase]-step;sciences=sum(n%STRIDE[phase]==0 for n in range(step+1,END[phase]+1))
            forecast=1.25*(remaining*intent['ordinary_step_forecast_seconds']+sciences*(44+intent['observer_forecast_seconds'])+120+20)
            if intent['deadline_monotonic']-time.monotonic()<forecast:next_step='RESOURCE_LIMIT'
        atomic(directory/f'step-{step}.ack.json',dict(schema='WJObserverAckV1',status=result['status'],ready_sha256=identity(path)['sha256'],binding=b,result=identity(dest),inventory=identity(inventory),
            qualified_cursor=identity(cursor)if result['status']=='PASS'else identity(dest),continuation=next_step))
        if verdict!='PROCEED':
            from run_wj_h20 import halt
            halt(root,'screen_halt'if verdict=='REVIEW_STOP'else'coverage_halt',dict(phase=phase,step=step,reason=verdict,evidence=identity(dest)));break
    while not(directory/'observer-stop.json').exists():time.sleep(.025)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--intent',type=Path,required=True);p.add_argument('--compat',type=Path,required=True);a=p.parse_args();daemon(a.root,a.intent,a.compat)
