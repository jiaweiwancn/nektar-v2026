#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Versioned passive H4 native readbacks and owned same-point comparisons."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h4 import load,save,sha,THREADS
from wj_h3_analysis import actual,identity,state_context,pressure_complement
from wj_h4_points import pair_rows,aligned,vector_hash
from wj_h4_study import CASES,OBSERVER,active_inventory
from measure_wj_h2 import checkpoint_for
from wj_h2_weak import independent_weak


def native_readback(root,name,rp,folder,role):
    folder.mkdir(parents=True,exist_ok=False);r=load(rp)
    try:cp,payload=checkpoint_for(dict(xml=str(rp.with_name('case.xml')),report=str(rp)),r['steps'])
    except AssertionError:
        if r['steps']!=0:raise
        # No accepted evolving checkpoint was produced. Bind the real startup
        # producer, never the failed candidate or a fabricated checkpoint.
        rp=root/'cases'/name/'cp0/case.evolution.json';r=load(rp)
        cp,payload=checkpoint_for(dict(xml=str(rp.with_name('case.xml')),report=str(rp)),0);role='startup_only:no_completed_stage'
    request=dict(operation='open_checks',inner_stress_model='Linear',native_report=str(rp),native_report_sha256=sha(rp),
        native_checkpoint=cp,checkpoint_sha256=sha(Path(cp)/'State.bin'),allow_failed_producer=bool(r['failure']),
        physical_observations=[[x,0,z]for x in (0.,3.,4.,8.)for z in (.5875,.6720111978440669,.8818157002137355,.9080394908593612,1.1947077409162705,2.,12.5)],
        weak=True,weak_outer_order=32,assembly=False,reuse=False)
    save(folder/'request.json',request);argv=[str(OBSERVER),str(rp.with_name('case.xml')),str(folder/'request.json'),str(folder/'native.json')]
    start=time.monotonic()
    with (folder/'command.log').open('x')as f:child=subprocess.run(argv,cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=180)
    if child.returncode:raise RuntimeError('native readback failed; retain original log')
    save(folder/'native-completion.json',dict(exit=child.returncode,seconds=time.monotonic()-start,argv=argv))
    return finish_readback(folder,role)


def finish_readback(folder,role,receipt_path=None):
    """Seal a complete passive packet; a missing process exit stays unknown."""
    request=load(folder/'request.json');rp=Path(request['native_report']);r=load(rp);cp=Path(request['native_checkpoint'])
    assert sha(rp)==request['native_report_sha256'] and sha(cp/'State.bin')==request['checkpoint_sha256']
    p=load(folder/'native.json');assert p['initialization_work']==[0,0] and p['governing_solves']==0
    assert p['report_sha256']==sha(rp) and p['checkpoint_sha256']==sha(cp/'State.bin')
    step=p['selected_absolute_step'];state,h,alpha,t,packet=state_context(r,step)
    assert p['state_sha256']==vector_hash(state)==packet['state_sha256']
    assert p['callback_history_sha256']==vector_hash(h)
    if step:assert p['callback_history_sha256']==packet['stage_history_sha256']
    completed=load(folder/'native-completion.json')if (folder/'native-completion.json').exists()else None
    result=dict(role=role,argv=[str(OBSERVER),str(rp.with_name('case.xml')),str(folder/'request.json'),str(folder/'native.json')],
        seconds=None if completed is None else completed['seconds'],process_exit='NOT_RECORDED'if completed is None else completed['exit'],
        recovered_complete_packet=completed is None,producer=identity(rp),checkpoint=identity(cp/'State.bin'),
        request=identity(folder/'request.json'),native=identity(folder/'native.json'),log=identity(folder/'command.log'),
        callback_history_sha256=p['callback_history_sha256'],stage_norms=p['physical_mass_dual_norms'],
        offline_pressure_projection=pressure_complement(r,r['state']),
        independent_weak=independent_weak(p['weak_checks'],dict(r,open_boundary=p['open_boundary']),enforce=False))
    save(receipt_path or folder/'receipt.json',result);return result


def observe(root):
    out=root/'native-observations-v1';out.mkdir();results={}
    for name in CASES:
        ready=root/'cases'/name/'prepared.json'
        if not ready.exists():continue
        if name in ('base-x1-half','local-x3'):
            rp=root/'cases'/name/'cp0/case.evolution.json'
            results[name+'-cp0']=native_readback(root,name,rp,out/(name+'-cp0'),'represented startup')
        path=root/'cases'/name/'trajectory.json'
        if path.exists():
            t=load(path);rp=Path(t['calls'][-1]['producer']['path'])
            results[name]=native_readback(root,name,rp,out/name,'terminal accepted state or accepted predecessor')
        save(out/'index.json',results)
    print('Native readbacks',[(k,v['independent_weak']['status'])for k,v in results.items()],flush=True)


def resume_observe(root):
    """New immutable view of completed v1 packets plus only missing readbacks."""
    old=root/'native-observations-v1';out=root/'native-observations-v2';out.mkdir();results={}
    tasks=[]
    for name in CASES:
        if not (root/'cases'/name/'prepared.json').exists():continue
        if name in ('base-x1-half','local-x3'):tasks.append((name+'-cp0',name,root/'cases'/name/'cp0/case.evolution.json','represented startup'))
        f=root/'cases'/name/'trajectory.json'
        if f.exists():tasks.append((name,name,Path(load(f)['calls'][-1]['producer']['path']),'terminal accepted state or predecessor'))
    for key,name,rp,role in tasks:
        previous=old/key
        if (previous/'receipt.json').exists():
            result=load(previous/'receipt.json')
            for field in ('producer','checkpoint','request','native','log'):
                assert sha(result[field]['path'])==result[field]['sha256'],'changed completed passive input/output'
        elif (previous/'native.json').exists():
            # Full native packet exists; the parent timed out during independent
            # algebra. Verify it against original inputs; do not invent its exit.
            # Keep the original folder untouched: the recovered receipt is v2.
            result=finish_readback(previous,role,out/(key+'-recovery.json'))
        else:result=native_readback(root,name,rp,out/key,role)
        results[key]=result;save(out/'index.json',results)
    save(root/'native-observation-view.json',identity(out/'index.json'))
    print('Versioned readbacks complete',[(k,v['independent_weak']['status'])for k,v in results.items()],flush=True)


def group_name(case):
    aliases={'retained-B-zero-x1':'base-x1','retained-B-zero-x3':'base-x3','retained-A-dt02':'A02','retained-A-dt01':'A01'}
    if case in aliases:return aliases[case]
    return next((n for n in sorted(CASES,key=len,reverse=True)if case.startswith(n+'-')),case)


def samples(root):
    inv=active_inventory(root);catalog=load(root/'catalog.json');groups={}
    for entry in catalog:
        assert sha(entry['path'])==entry['sha256'],'scan bytes differ from sealed catalog'
        p=load(entry['path']);assert sha(p['producer']['path'])==p['producer']['sha256'],'wrong producer bytes'
        r=load(p['producer']['path']);name=group_name(entry['case'])
        step=entry['step'] if entry['step'] is not None else r['steps']+1
        intended=aligned(step,r['dt'],entry['time']);role='rejected' if entry['native_status']=='rejected candidate' else 'accepted'
        state,history,alpha,t,packet=state_context(r,step,role=='rejected')
        assert vector_hash(state)==p['state_sha256']==packet['state_sha256']
        assert p['callback_history_sha256']==packet['stage_history_sha256']
        if step:assert vector_hash(history)==p['callback_history_sha256']
        assert p['inventory_sha256']==inv['sha256'],'stale scan in final catalog'
        key=(int(step),role)
        if key in groups.setdefault(name,{}):raise ValueError('duplicate native stage binding')
        groups[name][key]=dict(entry,dt=r['dt'],step=step,intended=intended,role=role,producer=p['producer'])
    return inv,groups


def flatten(pair,a,b):
    v=dict(point_id=pair['point_id'],x=pair['physical'][0],z=pair['physical'][2],region=pair['region'],
        intended_time=pair['intended_time'],recorded_time_a=pair['recorded_time_a'],recorded_time_b=pair['recorded_time_b'])
    for side,row in (('a',a),('b',b)):
        for key in ('geometry','element','x_side','normal_side','interface_role'):v[key+'_'+side]=row['owner'][key]
        for key in ('k','epsilon','epsilon_scope','tau','pi','sampled_K','N','F_N','root_residual','nu_eff','nu_d','Praw','Pdev','raw_divergence'):v[key+'_'+side]=row.get(key)
        for i,axis in enumerate(('x','y','z')):
            v['U'+axis+'_'+side]=row['U'][i]
            v['covariance_'+str(i)+'_'+side]=row['covariance'][i]
            for j,axis2 in enumerate(('x','y','z')):v['G'+axis+axis2+'_'+side]=row['G'][i][j]
        v['principal_'+side]=None if row['principal']is None else row['principal']['in_plane']
        v['transverse_'+side]=None if row['principal']is None else row['principal']['transverse']
        for vector in ('xi','polarization'):
            for j in range(3):v[vector+str(j)+'_'+side]=None if row['principal']is None else row['principal'][vector][j]
    v['delta_Gxz']=pair['delta_G'][0][2];v['delta_tau']=pair['delta_tau'];v['delta_pi']=pair['delta_pi']
    v['delta_principal']=None if pair['margin_a']is None else pair['margin_b']-pair['margin_a']
    return v


def compare(root):
    out=root/'comparisons-v1';out.mkdir();inv,groups=samples(root);save(out/'inventory.json',inv)
    pairs=[('time-base-x1','base-x1','base-x1-half'),('time-refined-x1','refined-x1','refined-x1-half'),
        ('time-base-x3','base-x3','base-x3-half'),('normal-dt01','base-x1','refined-x1'),
        ('normal-dt005','base-x1-half','refined-x1-half'),('horizontal-dt01','base-x1','base-x3'),
        ('horizontal-dt005','base-x1-half','base-x3-half'),('local-dt01','base-x3','local-x3'),
        ('local-dt005','base-x3-half','local-x3-half'),('retained-A-time','A02','A01')]
    summary={};minima=[];initial={}
    for name,g in groups.items():
        for entry in g.values():
            p=load(entry['path']);m=p['minimum'];w=m['row']
            minima.append(dict(case=name,time=entry['time'],step=entry['step'],role=entry['role'],
                native_status=entry['native_status'],expanded_status=entry['expanded_status'],prefix_qualified=entry['prefix_qualified'],
                margin=m['margin'],x=w['physical'][0],z=w['physical'][2],geometry=w['owner']['geometry'],
                x_side=w['owner']['x_side'],normal_side=w['owner']['normal_side'],Gxz=w['G'][0][2],tau=w['tau'],
                k=w['k'],epsilon=w['epsilon'],nu_eff=w['nu_eff'],nu_d=w['nu_d'],
                state_sha256=p['state_sha256'],callback_history_sha256=p['callback_history_sha256'],scan=entry['path']))
    def matching_times(g):return {e['intended']:e for e in g.values()}
    for label,na,nb in pairs:
        if na not in groups or nb not in groups:summary[label]=dict(status='NOT_RUN',reason='missing trajectory');continue
        ga,gb=matching_times(groups[na]),matching_times(groups[nb]);common=sorted(ga.keys()&gb.keys());records=[]
        full=out/(label+'.csv.gz');critical=out/(label+'-limiting.csv.gz')
        with gzip.open(full,'wt',newline='')as ff,gzip.open(critical,'wt',newline='')as cf:
            fw=cw=None
            for target in common:
                # Only the declared scientific cadence, retaining rejected
                # candidates separately where an actual common target exists.
                if (target/ __import__('fractions').Fraction('.01')).denominator!=1:continue
                ea,eb=ga[target],gb[target];a,b=load(ea['path']),load(eb['path'])
                rows=pair_rows(a,b,ea['step'],ea['dt'],ea['time'],eb['step'],eb['dt'],eb['time'])
                lookup=lambda p:{(v['point_id'],v['owner']['region'],v['owner']['geometry']):v for v in p['rows']}
                aa,bb=lookup(a),lookup(b);limiting={a['minimum']['row']['point_id'],b['minimum']['row']['point_id']}
                extrema={};locations={};limiting_rows=[]
                for pair in rows:
                    x=aa[(pair['point_id'],pair['region'],pair['owner_a']['geometry'])];y=bb[(pair['point_id'],pair['region'],pair['owner_b']['geometry'])]
                    flat=flatten(pair,x,y);flat.update(case_a=na,case_b=nb,native_role_a=ea['role'],native_role_b=eb['role'],
                        expanded_a=ea['expanded_status'],expanded_b=eb['expanded_status'])
                    if fw is None:
                        fw=csv.DictWriter(ff,fieldnames=list(flat));cw=csv.DictWriter(cf,fieldnames=list(flat));fw.writeheader();cw.writeheader()
                    fw.writerow(flat)
                    tags=inv['points'][pair['point_id']]['tags']
                    if pair['point_id']in limiting or any(t.startswith(('native_witness','declared_'))for t in tags):cw.writerow(flat)
                    if pair['point_id']in limiting:limiting_rows.append(flat)
                    for key in ('delta_Gxz','delta_tau','delta_pi','delta_principal'):
                        if flat[key]is None:continue
                        if key not in extrema or abs(flat[key])>extrema[key]:
                            extrema[key]=abs(flat[key]);locations[key]=dict(physical=pair['physical'],owner_a=pair['owner_a'],owner_b=pair['owner_b'],signed=flat[key])
                records.append(dict(intended_time=float(target),recorded_times=[ea['time'],eb['time']],
                    native_roles=[ea['role'],eb['role']],expanded_statuses=[ea['expanded_status'],eb['expanded_status']],
                    row_count=len(rows),minima=[a['minimum'],b['minimum']],same_point_at_limiting_locations=limiting_rows,
                    extrema=extrema,extrema_locations=locations))
                if target==0 and (label.startswith('time-') or label=='retained-A-time'):
                    initial[label]=dict(equal_primary_coefficients=a['state_sha256']==b['state_sha256'],state_sha256_a=a['state_sha256'],state_sha256_b=b['state_sha256'])
        expected=[i*.01 for i in range(9)]
        observed=[v['intended_time']for v in records]
        missing=[t for t in expected if not any(abs(t-u)<1e-14 for u in observed)]
        summary[label]=dict(status='COMPLETE_DATA' if not missing else 'PARTIAL_DATA',records=records,
            missing_actual_common_times=missing,
            rejected_or_unqualified_times=[v['intended_time']for v in records if v['native_roles']!=['accepted','accepted'] or v['expanded_statuses']!=['PASS','PASS']],
            qualification_scope='data coverage is separate from native and expanded physical admissibility',full=identity(full),limiting=identity(critical))
    save(out/'summary.json',summary);save(out/'minima.json',minima);save(out/'initial-physical-equality.json',initial)
    print('Common-point comparisons complete',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('observe','resume_observe','compare'),required=True)
    a=p.parse_args();globals()[a.mode](a.root)
