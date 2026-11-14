#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Passive common-support views of H15's immutable completed packets.

No launch, projection, lift, constitutive cache, quadrature change or solve.
The live observer's original decisions remain distinct from this final union.
"""
import argparse
from collections import OrderedDict
import csv
import json
import math
import os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import load_cache,verify,identity,atomic,pkey,okey,assess
from study_wj_h15 import qualified_step
from study_wj_h9 import serial
from wj_h5_neutral import KAPPA

FIELDS=('U','G','k','epsilon','R','pi','mechanical_pressure')

def scaled_norms(norms,c):
    scales=dict(U=1.,G=1.,k=c.K,epsilon=c.A/c.H,R=c.u_tau**2,pi=c.u_tau**2,mechanical_pressure=c.u_tau**2)
    return {region:{field:{key:value/scales[field]for key,value in values.items()}for field,values in fields.items()}for region,fields in norms.items()}

def fetch_rows(rows,data):
    """Signed sampled ranges at common physical heights/normal owners."""
    groups={}
    for i,(p,o)in enumerate(rows.rows[:rows.point_count]):
        key=(float(p[2]).hex(),float(o['lower'][2]).hex(),float(o['upper'][2]).hex(),o['normal_side'],o['region'])
        groups.setdefault(key,[]).append(i)
    result=[]
    for key,ii in groups.items():
        if len({float(rows.rows[i][0][0]).hex()for i in ii})<2:continue
        record=dict(z_code=rows.rows[ii[0]][0][2],z_mm=40*rows.rows[ii[0]][0][2],region=key[-1],normal_side=key[-2],normal_chart_interval=[float.fromhex(key[1]),float.fromhex(key[2])],sampled_x_range=[min(rows.rows[i][0][0]for i in ii),max(rows.rows[i][0][0]for i in ii)])
        for name in ('U','G','k','pi','K','rc'):
            values=data[name][ii];record[name]=dict(minimum=np.min(values,axis=0),maximum=np.max(values,axis=0),range=np.ptp(values,axis=0))
        if key[-1]=='outer':
            values=data['epsilon'][ii];record['epsilon']=dict(minimum=float(np.min(values)),maximum=float(np.max(values)),range=float(np.ptp(values)))
        result.append(record)
    return serial(result)

def exterior_rows(rows,data):
    """Actual stress/flux with fixed prescribed data separately labelled."""
    c=rows.config;out=[];xmin=min(e['lower'][0]for e in rows.report['fields'][0]['mesh']);xmax=max(e['upper'][0]for e in rows.report['fields'][0]['mesh'])
    for i,(p,o)in enumerate(rows.rows[:rows.point_count]):
        faces=[]
        if p[0]==xmin:faces.append(('inlet',[-1.,0,0]))
        if p[0]==xmax:faces.append(('outlet',[1.,0,0]))
        if p[2]==c.d0:faces.append(('ground',[0.,0,-1.]))
        if p[2]==c.top:faces.append(('top',[0.,0,1.]))
        for role,normal in faces:
            n=np.array(normal);un=float(data['U'][i]@n);traction=data['T'][i]@n
            value=dict(role=role,physical=p,owner=o,normal=normal,normal_velocity=un,actual_traction=traction,k_advective=float(data['k'][i]*un),k_diffusive_outward=float(-data['nu_d'][i]*data['grad_k'][i]@n))
            if o['region']=='outer':
                value['epsilon_advective']=float(data['epsilon'][i]*un);value['epsilon_diffusive_outward']=float(-data['nu_d'][i]/1.3*data['grad_epsilon'][i]@n)
            if role=='outlet':
                k=rows.exact['k'][i];pi=rows.exact['pi'][i]
                if o['region']=='outer':T=-rows.exact['R'][i]+(2*k/3-pi)*np.eye(3)
                else:
                    G=rows.exact['G'][i];T=c.D(p[2]-c.d0)*(G+G.T)-pi*np.eye(3)
                value['prescribed_complete_traction']=T@n
                value['prescribed_scalar_diffusive_outward']=0.
            elif role=='top':
                value['prescribed_tangential_traction']=[c.u_tau**2,0.,0.]
                value['normal_traction_role']='reaction to prescribed impermeable normal velocity'
                value['prescribed_k_diffusive_outward']=0.
                value['prescribed_epsilon_diffusive_outward']=c.u_tau*c.A*KAPPA/(1.3*(c.top-c.d0))
            else:value['traction_role']='state-dependent Dirichlet reaction; no copied final-state load'
            out.append(value)
    return serial(out)

def csv_rows(path,records):
    with Path(path).open('x',newline='')as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)

def analysis(root):
    folder=root/'science/target';end=qualified_step(folder);assert end>0
    out=root/'science/analysis-v1';out.mkdir()
    cache,manifest=load_cache(root/'science/observer-cache/manifest.json')
    inventory=load(folder/f'monitor/step-{end}.ready.inventory.json')
    rows,initial,stats=cache.view(inventory['points']);witness_points=set();reports={}
    for step in range(1,end+1):
        ready=load(folder/f'monitor/step-{step}.ready.json');r=load(verify(ready['report']));reports[step]=r
        assert r['original_qh0_state_sha256']==manifest['reference_identity']['qh0']
        for w in r['wj_audit_packets'][0]['principal_witnesses'].values():witness_points.add(pkey(w['position']))
    # Retained CP0/CP1 witnesses are already in the pinned base inventory.
    for name in ('cp0','pilot1'):
        p=load(root/f'science/observer-cache/{name}-parity.json')['original']['summary']['witness']['physical'];witness_points.add(pkey(p))
    selected=[i for i,(p,o)in enumerate(rows.rows[:rows.point_count])if pkey(p)in witness_points]
    initial_errors=rows.norms(initial,rows.exact);timeline=[];histories=[];step_records=[]
    for step in range(end+1):
        state=cache.qh0 if step==0 else reports[step]['state'];result,data=assess(rows,state,initial,step);s=result['summary']
        assert s['status']=='PASS','final-union retrospective qualification failed'
        current=dict(step=step,time=step*.0025,reused_initial_state=step==0,inventory=identity(folder/f'monitor/step-{end}.ready.inventory.json'),quality=s['status'],
            E_initial=initial_errors,E_eq=s['reference_error'],E_drift=s['initial_drift'],scaled_E_initial=scaled_norms(initial_errors,rows.config),scaled_E_eq=scaled_norms(s['reference_error'],rows.config),scaled_E_drift=scaled_norms(s['initial_drift'],rows.config),summary=s,all_regional_Q1_moments=result['all_regional_Q1_moments'])
        atomic(out/f'step-{step}.json',current)
        atomic(out/f'step-{step}-fetch.json',fetch_rows(rows,data))
        atomic(out/f'step-{step}-exterior.json',exterior_rows(rows,data))
        for i in selected:
            point=rows.witness(data,i)
            point.update(step=step,time=step*.0025,role='previous/current native or independent witness; all legitimate owners')
            if not rows.outer[i]:
                for k in ('margin','transverse','angle'):point[k]=None
            else:
                G=data['G'][i];B=data['k'][i]/data['epsilon'][i]*(G-np.trace(G)*np.eye(3)/3);S=(B+B.T)/2;O=(B-B.T)/2
                ss=float(np.trace(S@S));ww=float(np.trace(O@O));N=float(data['root'][i])
                point['root_derivative']=3*N*N-3.6*N-2*ww-2.7*ss
                point['root_residual']=(N-1.8)*(N*N-2*ww)-2.7*ss*N
            histories.append(point)
        row=dict(step=step,time=step*.0025,new_scientific_progress=step>=2,min_in_plane=s['min_in_plane'],min_transverse=s['min_transverse'],min_covariance=s['min_covariance'],max_strong_divergence=s['raw_divergence_max'],max_inner_Q1=result['all_regional_Q1_moments']['4']['maximum'],max_outer_Q1=result['all_regional_Q1_moments']['5']['maximum'],U_eq_L2=s['reference_error']['full']['U']['L2'],U_drift_L2=s['initial_drift']['full']['U']['L2'],E_G=s['errors']['E_G'],E_epsilon=s['errors']['E_epsilon'],E_U=s['errors']['E_U'],E_k=s['errors']['E_k'])
        timeline.append(row)
        if step:
            live=load(folder/f'monitor/step-{step}.ready.result.json');r=reports[step]
            step_records.append(dict(step=step,time=r['time'],native_hard_status='PASS',independent_live_status=live['status'],retrospective_final_union_status=s['status'],source_report=identity(folder/f'monitor/step-{step}.report.json'),native_compute_seconds=live['native_compute_seconds'],independent_seconds=live['measured_independent_seconds'],cache=live['cache'],replay=live['replay'],frozen_subsolve=r['pimple_attempts'][0],completed_stage=r['stage_audits']))
    atomic(out/'witness-history.json',histories);atomic(out/'native-and-independent.json',step_records);csv_rows(out/'timeline.csv',timeline)
    atomic(out/'summary.json',dict(status='PASS',qualified=end,new_scientific_progress=max(0,end-1),final_union_owner_rows=rows.point_count,physical_quadrature_rows=len(rows.rows)-rows.point_count,cache=stats,
        exact_reference_role='analytical continuum control; never a state reset',qh0_role='original H14 actual initialized state',internal_reference_role='original native fixed_reference/defect, separately hash-bound',comparison_scales=dict(U=1,G=1,k=rows.config.K,epsilon=rows.config.A/rows.config.H,R=rows.config.u_tau**2,pi=rows.config.u_tau**2),
        clocks=dict(horizon=505.02466521278404,terminal_window=[447.8483832973273,505.02466521278404],dt=.0025,required_steps=202010),references=manifest['reference_identity']))
    figures(out,timeline,step_records)

def figures(out,timeline,steps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=np.array([v['time']for v in timeline]);fig,axes=plt.subplots(1,2,figsize=(9,3.4),constrained_layout=True)
    axes[0].plot(t,[v['min_in_plane']for v in timeline],'o-',label='finite-inventory minimum')
    axes[0].axhline(.002068763923454,color='tab:red',ls='--',label='half-reference quality gate');axes[0].set(xlabel='T [Hh/U∞]',ylabel='m / νeff');axes[0].legend(fontsize=8)
    axes[1].plot(t,[v['U_eq_L2']for v in timeline],'o-',label='current − analytical reference')
    axes[1].plot(t,[v['U_drift_L2']for v in timeline],'s-',label='current − original initialized qh0');axes[1].set(xlabel='T [Hh/U∞]',ylabel='Physical full-domain velocity L2');axes[1].legend(fontsize=8)
    fig.suptitle('Fixed 6 mm full-domain neutral control; stage 1 is H14 replay');fig.savefig(out/'margins-and-drift.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.5,3.3),constrained_layout=True)
    ax.plot([v['step']for v in steps],[v['native_compute_seconds']for v in steps],'o-',label='native compute before publication/wait')
    ax.plot([v['step']for v in steps],[v['independent_seconds']for v in steps],'s-',label='independent assessment/publication')
    ax.set(xlabel='Absolute completed step (1 = replay)',ylabel='Wall seconds',xticks=[v['step']for v in steps]);ax.legend(fontsize=8);fig.savefig(out/'per-stage-cost.png',dpi=160);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();analysis(a.root)
