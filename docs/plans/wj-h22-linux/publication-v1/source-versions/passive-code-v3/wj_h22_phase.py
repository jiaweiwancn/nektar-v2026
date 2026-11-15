#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""All-seven actual-stage telescopes on exactly the two H21 owners."""
from pathlib import Path
import numpy as np
from run_wj_h4 import load
from wj_h15_observer import identity,verify,atomic,vector_hash
from wj_h17_observer import native_load
from wj_h21_outlet import linear,local,plain,scale,comparison,OWNERS,PAIRS
from wj_h12_phase import decompositions,identity as algebra_identity,cancellation
from wj_h22_replay import DT,END,PHASES,read_completed,recovery,phases,flow_vectors
from wj3_measure import checkpoint_payload
from wj_h9_local import same

NAMES=('U','G','devG','trace_G','grad_k','grad_epsilon','grad_pi','k','epsilon','pi','mechanical_pressure',
       'R','covariance','tau','B','d','nu_eff','nu_d','a','g','mhat','Praw','N','F_N','transverse')
TRANSITIONS=(('lift','P','C'),('flow','C','F'),('k','F','K'),('epsilon','K','A'))

def write(path,v):atomic(path,plain(v))

def branch_window(rows,name):
    if any(name not in row[p]for row in rows for p in PHASES):
        return dict(status='NOT_EVALUABLE',reason='an actual provisional constitutive observation is undefined')
    pieces=[{tag:np.asarray(row[b][name])-row[a][name]for tag,a,b in TRANSITIONS}for row in rows]
    terms={tag:sum(v[tag]for v in pieces)for tag,_,_ in TRANSITIONS}
    start=np.asarray(rows[0]['P'][name]);end=np.asarray(rows[-1]['A'][name])
    return dict(status='PASS',start=start,end=end,pieces=pieces,cumulative=terms,
        identity=algebra_identity(end-start-sum(terms.values()),[start,end,*terms.values()],scale(name)))

def cross_window(coarse,finer,name):
    if coarse['status']!='PASS'or finer['status']!='PASS':return dict(status='NOT_EVALUABLE')
    incoming=coarse['start']-finer['start'];end=coarse['end']-finer['end']
    assert np.all(incoming==0),'early window must not inherit H21 late-window differences'
    terms={tag:coarse['cumulative'][tag]-finer['cumulative'][tag]for tag,_,_ in TRANSITIONS}
    return dict(status='PASS',incoming=incoming,endpoint=end,contributions=terms,
        identity=algebra_identity(end-incoming-sum(terms.values()),[end,incoming,*terms.values()],scale(name)),
        cancellation=cancellation(incoming,list(terms.values()),end,scale(name)))

def flow_local(structure,states,values,increments,owners):
    records=[];dev=lambda G:G-np.trace(G)*np.eye(3)/3
    for j,owner in enumerate(owners):
        p=linear(structure,increments['predictor'],owner,j)
        c=linear(structure,increments['pressure'],owner,j)
        total=linear(structure,increments['total'],owner,j)
        U=np.asarray(values['F'][j]['U'])-values['C'][j]['U'];G=np.asarray(values['F'][j]['G'])-values['C'][j]['G']
        base=values['C'][j];tau=base['k']/base['epsilon']
        assert values['F'][j]['k']==base['k']and values['F'][j]['epsilon']==base['epsilon']
        Bp=tau*dev(p['G']);Bc=tau*dev(c['G']);Bg=tau*dev(G)
        r=dict(owner=owner,predictor_U=p['U'],pressure_U=c['U'],total_U=U,
             predictor_G=p['G'],pressure_G=c['G'],total_G=G,
             predictor_B=Bp,pressure_B=Bc,total_B=Bg,tau_C=tau,
             U_identity=algebra_identity(U-p['U']-c['U'],[values['F'][j]['U'],base['U'],p['U'],c['U']],1.),
             G_identity=algebra_identity(G-p['G']-c['G'],[values['F'][j]['G'],base['G'],p['G'],c['G']],1.),
             B_linear_identity=algebra_identity(Bg-Bp-Bc,[Bg,Bp,Bc],1.),
             coefficient_to_gradient_identity=algebra_identity(G-total['G'],[values['F'][j]['G'],base['G'],total['G']],1.),
             scope='genuine linear increments, not predictor-only or pressure-only constitutive states')
        if 'B'in values['F'][j]and 'B'in base:
            r['actual_flow_B_identity']=algebra_identity(values['F'][j]['B']-base['B']-Bg,[values['F'][j]['B'],base['B'],Bg],1.)
        records.append(r)
    return records

def run(root,levels=('C','M','F'),out_name='phase-analysis-v1'):
    assert tuple(levels)in (('C','M','F'),('C','M'))
    inputs=load(root/'inputs/selected.json');structure=load(verify(inputs['sparse_structure']))
    manifest=load(verify(inputs['H21_endpoints']));files={Path(v['path']).name:v for v in manifest['files']}
    owners=load(verify(files['control--selected.json']))['owners']
    expected=load(verify(files['endpoint-values.json']))
    out=root/'science'/out_name;out.mkdir()
    cache={};buffers={};role_records=[];attempts=[];flows=[];transitions=[];windows={};endpoint_checks=[]
    for z in levels:
        predecessor=checkpoint_payload(verify(inputs['cp0'][z]['checkpoint']).parent,include_reference=True)
        predecessor_binding=dict(role='original own-dt CP0',anchor=inputs['cp0'][z],step=0,time=0.,alpha=0.,
            callback_history_sha256=None,completed_callback='NONE_INITIAL_ROLE',native_history_sha256=[vector_hash(h)for h in predecessor['native_history']])
        branch_rows=[[],[]]
        for n in range(1,END[z]+1):
            anchor=recovery(root/f'science/{z}/monitor',n);report,cp=read_completed(anchor,z,n)
            a=report['pimple_attempts'][0];sp=dict(path=a['k_snapshot'],bytes=Path(a['k_snapshot']).stat().st_size,sha256=a['k_snapshot_sha256'])
            packet=native_load(verify(sp));states,binding=phases(report,packet,predecessor,structure['offsets'])
            inc,vectors=flow_vectors(a['flow_vectors'],states,report,structure['offsets'],identity(root/'science/flow-vector-execution.json'))
            values={}
            for phase in PHASES:
                state=states[phase];digest=vector_hash(state)
                if digest in cache:assert same(buffers[digest],state),'alias hash collision or source corruption'
                else:
                    buffers[digest]=state.copy();cache[digest]=[local(linear(structure,state,o,j))for j,o in enumerate(owners)]
                values[phase]=cache[digest]
                for j,owner in enumerate(owners):
                    role_records.append(dict(level=z,step=n,phase=phase,owner=owner,state_sha256=digest,values=values[phase][j],
                        role='native accepted predecessor'if phase=='P'else'native completed endpoint'if phase=='A'else'actual provisional phase',
                        native_full_audit_status=report['native_full_audit_status']if phase=='A'else'NOT_AN_ACCEPTED_STATE'if phase!='P'else'predecessor own role',
                        predecessor_binding=predecessor_binding if phase=='P'else None,
                        independent_qualification=('COMMITTED'if (root/f'science/{z}/monitor/step-{n}.qualification-commit.json').exists()else'NOT_SCHEDULED'if report['native_full_audit_status']=='NOT_SCHEDULED'else'NOT_COMMITTED')if phase=='A'else'NOT_A_QUALIFICATION_ROLE',
                        attempt_callback_sha256=vector_hash(report['last_stage_history']),attempt_time=report['time'],alpha=report['last_stage_alpha']))
            for j,(geo,side)in enumerate(OWNERS):
                row={p:values[p][j]for p in PHASES};branch_rows[j].append(row)
                for tag,before,after in TRANSITIONS:
                    transitions.append(dict(level=z,step=n,geometry=geo,phase=tag,
                        values=decompositions(row[before],row[after]),
                        signed={name:(np.asarray(row[after][name])-row[before][name])for name in NAMES if name in row[before]and name in row[after]}))
                if n==END[z]:
                    old=next(v['values']for v in expected if v['level']==z and v['time']==.005 and v['geometry']==geo)
                    checks={name:comparison(row['A'][name],old[name],scale(name))for name in NAMES if name in row['A']and name in old}
                    endpoint_checks.append(dict(level=z,geometry=geo,status='PASS',checks=checks))
            flows.append(dict(level=z,step=n,time=report['time'],owners=flow_local(structure,states,values,inc,owners)))
            attempts.append(dict(level=z,step=n,current=anchor,phase_snapshot=sp,binding=binding,vectors=vectors,
                 matrix=dict(path=packet['matrix_file'],sha256=packet['matrix_sha256']),matrix_operations=0))
            predecessor=cp
            predecessor_binding=dict(role='genuine native completed predecessor',anchor=anchor,step=n,time=report['time'],alpha=report['last_stage_alpha'],
                callback_history_sha256=vector_hash(report['last_stage_history']),completed_callback=anchor['callback'],
                native_history_sha256=[vector_hash(h)for h in cp['native_history']])
        for j,(geo,side)in enumerate(OWNERS):
            windows[z,geo]={name:branch_window(branch_rows[j],name)for name in NAMES}
    paired=[];branch=[];flow_windows=[]
    for (z,geo),quantities in windows.items():branch.append(dict(level=z,geometry=geo,quantities=quantities))
    for geo,_ in OWNERS:
        for a,b in PAIRS:
            if a not in levels or b not in levels:continue
            paired.append(dict(pair=a+'-'+b,geometry=geo,quantities={name:cross_window(windows[a,geo][name],windows[b,geo][name],name)for name in NAMES}))
        for z in levels:
            selected=[next(v for v in f['owners']if v['owner']['owner']['geometry']==geo)for f in flows if f['level']==z]
            parts={name:sum(np.asarray(v[name])for v in selected)for name in ('predictor_U','pressure_U','total_U','predictor_G','pressure_G','total_G','predictor_B','pressure_B','total_B')}
            checks={name:algebra_identity(parts['total_'+name]-parts['predictor_'+name]-parts['pressure_'+name],
                    [parts['total_'+name],parts['predictor_'+name],parts['pressure_'+name]],1.)for name in ('U','G','B')}
            if windows[z,geo]['B']['status']=='PASS':
                checks['actual_cumulative_flow_B']=algebra_identity(windows[z,geo]['B']['cumulative']['flow']-parts['total_B'],
                    [windows[z,geo]['B']['cumulative']['flow'],parts['total_B']],1.)
            flow_windows.append(dict(level=z,geometry=geo,contributions=parts,checks=checks))
    number=sum(END[z]for z in levels)
    assert len(role_records)==10*number and len(attempts)==len(flows)==number
    write(out/'roles.json',role_records);write(out/'attempt-bindings.json',attempts);write(out/'transition-decompositions.json',transitions)
    write(out/'branch-telescopes.json',branch);write(out/'cross-level-telescopes.json',paired)
    write(out/'genuine-flow-increments.json',flows);write(out/'cumulative-flow-split.json',flow_windows);write(out/'endpoint-H21-replay.json',endpoint_checks)
    write(out/'counts.json',dict(requested_attempts=7,attempts=number,phase_owner_roles=8*number,predecessor_owner_roles=2*number,total_roles=10*number,
          distinct_state_vectors=len(cache),distinct_constitutive_owner_evaluations=2*len(cache),predictor_pressure_packets=number,
          endpoint_replay_owners=2*len(levels),branch_windows=2*len(levels),cross_level_windows=len(paired),
          saved_matrix_operations=0,new_governing_solves=0,new_strict_reads=0,new_kernel_processes=0,
          undefined_phase_owner_roles=sum(r['values']['constitutive_status']!='DEFINED'for r in role_records),
          coverage='COMPLETE'if number==7 else'PARTIAL_C_M_ONLY; F_NOT_RUN_AFTER_OBSERVER_HALT'))
    print('H22 retained actual-phase/flow-vector/window identities PASS; coverage',number,'of seven attempts',flush=True)
