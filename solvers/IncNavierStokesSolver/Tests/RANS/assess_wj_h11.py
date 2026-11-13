#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H11 final physical-union comparisons and only the two selected saved systems."""
import argparse
from pathlib import Path
import math
import numpy as np
from run_wj_h11 import H10,load,save,sha
from study_wj_h11 import parent,audit,halt
from study_wj_h9 import identity,check,serial,norm_metrics
from wj_h11_analysis import MATRIX,FACE18,ratio,selected,common_step
from wj_h10_analysis import arithmetic,matrix,gates
from wj_h7_analysis import Rows,context,classify
from wj_h8_analysis import extend_witnesses,FACE,local_rows
from wj_h9_local import MomentAudit,phase_bindings,actions
from wj_h4_points import vector_hash
from study_wj_h6 import configuration
from validate_wj1 import Probe


def arithmetic_views(root):
    out=root/'science/arithmetic-v1';out.mkdir();results={};count=0
    for branch,dt,n,target in MATRIX:
        case=parent(root,branch)
        if branch=='coarse':
            r=load(case['parent']['producer']['path']);step=1
        else:
            if not Path(case['report']).is_file():results[branch]=dict(status='NOT_RUN');continue
            r=load(case['report']);step=2 if branch=='fine'else 4
            if r['steps']<step:results[branch]=dict(status='NOT_RUN');continue
        cx,a,index=selected(r,step);rp=Path(case['parent']['producer']['path']if branch=='coarse'else case['report'])
        check(dict(path=a['k_snapshot'],bytes=Path(a['k_snapshot']).stat().st_size,sha256=a['k_snapshot_sha256']));snap=load(a['k_snapshot'])
        if branch=='coarse':
            old=case['parent']['arithmetic'];check(old);ar=load(old['path']);ar['reuse_identity']=old
        else:
            m=load(a['momentum_arithmetic']['path']);p=load(a['pressure_arithmetic']['path'])
            for d in (a['momentum_arithmetic'],a['pressure_arithmetic']):assert sha(d['path'])==d['sha256']
            for d in (m,p):assert d['time']==cx['time']and d['alpha_mass']==cx['alpha']and d['history_sha256']==vector_hash(cx['history'])
            inv=load(root/'science/inventory-parent.json');rows=Rows(r,[v['physical']for v in inv['points']],configuration(6.));mom=MomentAudit(rows)
            ar=arithmetic(r,snap,m,p,matrix(m['mass_or_correction_map']['path'],m['mass_or_correction_map']),matrix(p['mass_or_correction_map']['path'],p['mass_or_correction_map']),matrix(p['matrix']['path'],p['matrix']),matrix(m['matrix']['path'],m['matrix']),mom,step=step)
            delta=np.asarray(snap['post_flow'])-snap['committed'];pi=rows.field(delta,4)[0]+rows.field(delta,5)[0]
            ar['pressure_physical']=norm_metrics(rows,pi);ar['impulse_physical']=norm_metrics(rows,dt*pi);profiles=[]
            for label,values in [('predictor',ar['predictor']),('pressure_correction',ar['pressure_correction_velocity']),('total',ar['total_velocity_increment']),('predictor_per_dt',ar['predictor_per_dt'])]:
                vector=np.zeros(len(delta));vector[np.array(ar['velocity_rows'],int)]=values;u=rows.field(vector,0)[0];w=rows.field(vector,1)[0];ar[label+'_physical_L2']=float(np.sqrt(rows.weights@(u*u+w*w)))
            for j,(pt,owner)in enumerate(rows.rows[:rows.point_count]):
                if pt[0]in(0.,4.,8.)and min(abs(pt[2]-z)for z in(.2375,FACE,FACE18))<=1.01e-6:profiles.append(dict(physical=pt,owner=owner,delta_pi=float(pi[j]),j_dt=float(dt*pi[j])))
            ar['pressure_profiles']=profiles;ar['diagnostic_mass_factorizations']=mom.projection.factorizations;ar['diagnostic_mass_solves']=mom.solves;count+=1
        # Full H10 profile rows remain at their original sealed identity.
        if branch=='coarse':ar['pressure_profiles']=[v for v in ar['pressure_profiles']if min(abs(v['physical'][2]-z)for z in(.2375,FACE,FACE18))<=1.01e-6]
        ar.update(producer=identity(rp),absolute_step=step,time=cx['time'],alpha=cx['alpha'],different_increment_durations=True)
        save(out/(branch+'.json'),serial(ar));results[branch]=dict(status=ar['status'],packet=identity(out/(branch+'.json')))
    save(out/'result.json',dict(status='PASS'if len(results)==3 and all(v['status']=='PASS'for v in results.values())else 'NOT_RUN',branches=results,new_attempt_level_arithmetic_sets=count,
        scope='coarse H10 step1 reused, fine H11 step2 and quarter H11 step4; all end at.01 with different incoming times/durations; matvecs only'))


def compare_fields(rows,data,time,c):
    result=[]
    scales=dict(U=1.,G=1.,k=c.K,epsilon=c.A/c.H,R=c.K,pi=c.u_tau**2,mechanical_pressure=c.u_tau**2,margin=1.)
    for scope,mask in [('full',np.ones(len(rows.rows),bool)),('inner',~rows.outer),('outer',rows.outer)]:
        for field in ('U','G','k','epsilon','R','pi','mechanical_pressure','margin'):
            if field in ('epsilon','margin')and scope=='inner':continue
            active=mask&(rows.outer if field in('epsilon','margin')else True);measure=float(np.sum(rows.weights[active]));a,b,d=[v[field][active]for v in data]
            for norm in ('L2','maximum'):
                def value(v):
                    v=v.reshape(len(v),-1);square=np.sum(v*v,axis=1)
                    return float(np.sqrt(rows.weights[active]@square))if norm=='L2'else float(np.sqrt(max(square)))
                floor=4e-11*scales[field]*(math.sqrt(measure)if norm=='L2'else 1.)
                result.append(dict(time=time,scope=scope,field=field,norm=norm,scale=scales[field],measure=measure,**ratio(value(a-b),value(b-d),floor)))
    return result


def common(root):
    out=root/'science/common-v1';out.mkdir();c=configuration(6.);inv=load(root/'science/inventory-parent.json');reports={};available={};base=None
    for branch,dt,n,target in MATRIX:
        case=parent(root,branch);old=load(case['parent']['producer']['path']);base=base or old
        available[branch]={s:(old,context(old,s),case['parent']['producer'])for s in(0,1)}
        if Path(case['report']).is_file():
            r=load(case['report']);reports[branch]=r;inv=extend_witnesses(r,inv,c)
            for s in range(2,r['steps']+1):available[branch][s]=(r,context(r,s),identity(Path(case['report'])))
        for r0 in [old,*([reports[branch]]if branch in reports else [])]:assert r0['fields']==base['fields']and r0['parameters']==base['parameters']
    save(out/'inventory.json',inv);rows=Rows(base,[v['physical']for v in inv['points']],c);mom=MomentAudit(rows);initial=rows.evaluate(base['accepted_states'][0]);records=[];projections=[];pairs=[];primary=[];fixed=[];seen={};kernel=Probe(str(H10/'runtime/WJKernelProbe'));queries=0;pressure_solves=0
    primary_ids=[j for j,(p,o)in enumerate(rows.rows[:rows.point_count])if p[0]==8. and p[2]in(FACE,FACE18)and o['region']=='outer']
    assert len(primary_ids)==4,'both owners of both exact physical faces required'
    local_keys={tuple(v)for v in inv['local_keys']};local_ids=[j for j,(p,o)in enumerate(rows.rows[:rows.point_count])if tuple(float(x).hex()for x in p)in local_keys or(p[0]==8. and min(abs(p[2]-z)for z in(FACE,FACE18))<=1.01e-6)]
    try:
        with(out/'local-all-states.jsonl').open('x')as local_file:
            for quarter_step in range(17):
                time=quarter_step*.0025;group={};group_cache={}
                for branch,dt,n,target in MATRIX:
                    if quarter_step%round(dt/.0025):continue
                    step=quarter_step//round(dt/.0025)
                    if step not in available[branch]:continue
                    r,cx,producer=available[branch][step];h=vector_hash(cx['state'])
                    if h in group_cache:data=group_cache[h]
                    else:data=rows.evaluate(cx['state']);group_cache[h]=data
                    group[branch]=data
                    if h not in seen:
                        summary=rows.summary(data,initial);summary['divergence_scaled_max']=float(np.max(abs(data['divergence'])/np.maximum(1.,np.linalg.norm(data['G'],axis=(1,2)))))
                        pp=mom.projection.measure(data);pressure_solves+=pp['mass_solves'];summary['moments']={k:dict(moments=v['pressure_test_moments'],projection_L2=v['continuous_Q1_projection_L2'])for k,v in pp['regions'].items()}
                        ii=np.flatnonzero(rows.outer);j=ii[np.argmin(data['margin'][ii])];summary['minimum_volume_inclusive']=float(data['margin'][j]);summary['volume_inclusive_witness']=rows.witness(data,j)
                        status=gates(summary,summary['moments'],startup=False);all_status=classify(data['margin'][j],float(np.min(data['transverse'][ii])),float(np.min(data['covariance'])),float(min(np.min(data['k']),np.min(data['epsilon'][ii]))))
                        if all_status!='PASS':status=all_status
                        H=np.outer([-np.sin(data['angle'][j]),0,np.cos(data['angle'][j])],[np.cos(data['angle'][j]),0,np.sin(data['angle'][j])])
                        v=kernel.closure(data['k'][j],data['epsilon'][j],data['G'][j]);dv=kernel.closure(data['k'][j],data['epsilon'][j],data['G'][j],(0,0,H));queries+=2
                        compiled=-np.sum(H*dv['R'])/v['scalar'][1];assert abs(compiled-data['margin'][j])<2e-11
                        seen[h]=dict(status=status,summary=serial(summary),projection=serial(pp),compiled_replay=float(compiled))
                    record=dict(branch=branch,step=step,time=cx['time'],comparison_time=time,alpha=cx['alpha'],producer=producer,state_sha256=h,history_sha256=vector_hash(cx['history']),**seen[h]);records.append(record)
                    projections.append(dict(branch=branch,step=step,time=cx['time'],**seen[h]['projection']))
                    loc=local_rows(rows,data,local_ids);local_file.write(__import__('json').dumps(serial(dict(branch=branch,step=step,time=cx['time'],records=loc)),allow_nan=False)+'\n')
                    for j in primary_ids:
                        fixed.append(dict(branch=branch,step=step,time=cx['time'],record=serial(rows.witness(data,j))))
                    if record['status']!='PASS':
                        if record['status']in('OBSTRUCTED','UNRESOLVED'):
                            from run_wj_h7 import record_coverage_halt
                            record_coverage_halt(root/'accounting.json',dict(status=record['status'],role='accepted',trajectory=branch,branch=branch,step=step,time=cx['time'],witness=seen[h]['summary']['witness']))
                        halt(root,record['status'],'final common-union retrospective qualification failed',branch=branch,step=step,time=cx['time'])
                if quarter_step%4==0:
                    if set(group)!=set(v[0]for v in MATRIX):
                        primary.append(dict(time=time,status='NOT_RUN',reason='required same-time branch state missing'));continue
                    data=[group[b]for b,dt,n,target in MATRIX];pairs.extend(compare_fields(rows,data,time,c))
                    for j in primary_ids:
                        a,b,d=[v['margin'][j]for v in data];primary.append(dict(time=time,physical=rows.rows[j][0],owner=rows.rows[j][1],values=[float(a),float(b),float(d)],signed_D1=float(a-b),signed_D2=float(b-d),**ratio(abs(a-b),abs(b-d),4e-11)))
                    active=np.flatnonzero(rows.outer);d1=data[0]['margin'][active]-data[1]['margin'][active];d2=data[1]['margin'][active]-data[2]['margin'][active];j1=active[np.argmax(abs(d1))];j2=active[np.argmax(abs(d2))]
                    primary.append(dict(time=time,scope='maximum same-owner outer support including volume nodes',D1_witness=dict(physical=rows.rows[j1][0],owner=rows.rows[j1][1],signed=float(d1[np.argmax(abs(d1))])),D2_witness=dict(physical=rows.rows[j2][0],owner=rows.rows[j2][1],signed=float(d2[np.argmax(abs(d2))])),**ratio(float(max(abs(d1))),float(max(abs(d2))),4e-11)))
    finally:kernel.close()
    save(out/'states.json',serial(records));save(out/'pressure-projections.json',serial(projections));save(out/'field-pairs.json',serial(pairs));save(out/'principal-pairs.json',serial(primary));save(out/'fixed-owner-history.json',serial(fixed))
    required={(b,s)for b,dt,n,t in MATRIX for s in range(t+1)};present={(v['branch'],v['step'])for v in records};complete=present==required
    positive=[v for v in primary if v['time']>0];temporal='NOT_RUN'if not complete else 'UNRESOLVED'if any(v['status']in('UNRESOLVED','DOES_NOT_CONTRACT','NOT_RUN')for v in positive)else 'CONTRACTS'
    status='PASS'if complete and all(v['status']=='PASS'for v in records)else 'OBSTRUCTED'if any(v['status']!='PASS'for v in records)else 'NOT_RUN'
    save(out/'result.json',dict(status=status,temporal=temporal,complete=complete,roles=len(records),unique_physical_states=len(seen),missing=[list(v)for v in sorted(required-present)],
        inventory=identity(out/'inventory.json'),owner_rows=rows.point_count,physical_volume_rows=len(rows.rows)-rows.point_count,diagnostic_mass_factorizations=mom.projection.factorizations,diagnostic_mass_solves=pressure_solves,kernel_queries=queries,
        scope='all roles use original H10 step0/1 producers or genuine new H11 callbacks; comparisons use final identical physical rows/owners; finite inventory only'))
    print('H11 final common',status,'temporal',temporal,'roles',len(records),'owners',rows.point_count,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('arithmetic','common'),required=True);a=p.parse_args();(arithmetic_views if a.mode=='arithmetic'else common)(a.root)
