#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Final common-union replay of already completed H16 packets; no native work."""
import argparse
import os
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,verify,atomic,pkey,assess,H14
from wj_h16_observer import load_compatible
from study_wj_h15 import qualified_step
from finish_wj_h15 import scaled_norms,fetch_rows,exterior_rows,csv_rows

def analysis(root):
    folder=root/'science/target';summary=load(folder/'case.evolution.json');end=summary['native_completed'];qualified=qualified_step(folder)
    out=root/'science/analysis-v1';out.mkdir();cache,manifest=load_compatible(root/'science/target-cache-compatibility.json')
    assert end>0,'no completed target packet to analyze'
    inventory=load(folder/f'monitor/step-{end}.ready.inventory.json');rows,initial,stats=cache.view(inventory['points']);reports={};witness_points=set()
    for step in range(1,end+1):
        ready=load(folder/f'monitor/step-{step}.ready.json');r=load(verify(ready['report']));reports[step]=r
        assert r['original_qh0_state_sha256']==manifest['reference_identity']['qh0']
        for w in r['wj_audit_packets'][0]['principal_witnesses'].values():witness_points.add(pkey(w['position']))
    for name in ('cp0','pilot1'):witness_points.add(pkey(load(root/f'science/retained-parity/{name}.json')['summary']['witness']['physical']))
    selected=[i for i,(p,o)in enumerate(rows.rows[:rows.point_count])if pkey(p)in witness_points]
    initial_errors=rows.norms(initial,rows.exact);timeline=[];history=[];steps=[]
    for step in range(end+1):
        state=cache.qh0 if step==0 else reports[step]['state'];result,data=assess(rows,state,initial,step);s=result['summary']
        record=dict(step=step,time=step*.0025,qualified_restart_eligible=step<=qualified,final_union_status=s['status'],E_initial=initial_errors,E_eq=s['reference_error'],E_drift=s['initial_drift'],scaled_E_initial=scaled_norms(initial_errors,rows.config),scaled_E_eq=scaled_norms(s['reference_error'],rows.config),scaled_E_drift=scaled_norms(s['initial_drift'],rows.config),summary=s,all_regional_Q1_moments=result['all_regional_Q1_moments'])
        atomic(out/f'step-{step}.json',record);atomic(out/f'step-{step}-fetch.json',fetch_rows(rows,data));atomic(out/f'step-{step}-exterior.json',exterior_rows(rows,data))
        for i in selected:
            v=rows.witness(data,i);v.update(step=step,time=step*.0025,role='previous/current native or independent witness; full legitimate owners')
            if rows.outer[i]:
                G=data['G'][i];B=data['k'][i]/data['epsilon'][i]*(G-np.trace(G)*np.eye(3)/3);S=(B+B.T)/2;O=(B-B.T)/2;ss=float(np.trace(S@S));ww=float(np.trace(O@O));N=float(data['root'][i]);v['root_derivative']=3*N*N-3.6*N-2*ww-2.7*ss;v['root_residual']=(N-1.8)*(N*N-2*ww)-2.7*ss*N
            else:
                for k in ('margin','transverse','angle'):v[k]=None
            history.append(v)
        timeline.append(dict(step=step,time=step*.0025,new_scientific_progress=step>=2 and step<=qualified,qualified_restart_eligible=step<=qualified,quality=s['status'],min_in_plane=s['min_in_plane'],min_transverse=s['min_transverse'],min_covariance=s['min_covariance'],max_strong_divergence=s['raw_divergence_max'],max_inner_Q1=result['all_regional_Q1_moments']['4']['maximum'],max_outer_Q1=result['all_regional_Q1_moments']['5']['maximum'],U_eq_L2=s['reference_error']['full']['U']['L2'],U_drift_L2=s['initial_drift']['full']['U']['L2'],E_G=s['errors']['E_G'],E_epsilon=s['errors']['E_epsilon'],E_U=s['errors']['E_U'],E_k=s['errors']['E_k']))
        if step:
            live=load(folder/f'monitor/step-{step}.ready.result.json');ready=load(folder/f'monitor/step-{step}.ready.json');r=reports[step]
            steps.append(dict(step=step,time=r['time'],native_hard='PASS',independent_live=live['status'],final_union=s['status'],report=ready['report'],checkpoint=ready['checkpoint'],callback=ready['callback'],commit=identity(folder/f'monitor/step-{step}.qualification-commit.json')if step<=qualified else None,native_compute_seconds=ready['native_compute_seconds'],independent_seconds=live['measured_independent_seconds'],cache=live['cache'],inventory=identity(folder/f'monitor/step-{step}.ready.inventory.json'),replay=live['replay'],frozen_subsolve=r['pimple_attempts'][0],completed_stage=r['stage_audits'],step_audits=r['step_audits'],explicit_calls=r['actual_explicit_calls']))
    atomic(out/'witness-history.json',history);atomic(out/'native-and-independent.json',steps);csv_rows(out/'timeline.csv',timeline)
    warm=[v for v in steps if v['step']>=2 and v['step']<=qualified];native=float(np.mean([v['native_compute_seconds']for v in warm]))if warm else None
    observed=float(np.mean([v['native_compute_seconds']+v['independent_seconds']for v in warm]))if warm else None
    atomic(out/'summary.json',dict(status='PASS'if qualified==end==4 and all(v['quality']=='PASS'for v in timeline)else'OBSTRUCTED_OR_INCOMPLETE',native_completed=end,qualified=qualified,new_scientific_progress=max(0,qualified-1),final_union_owner_rows=rows.point_count,quadrature_rows=len(rows.rows)-rows.point_count,final_required_points=len(inventory['points']),cache=stats,reference_identity=manifest['reference_identity'],clocks=dict(horizon=505.02466521278404,terminal_window=[447.8483832973273,505.02466521278404],dt=.0025,required_steps=202010),costs=dict(warm_samples=len(warm),native_mean_seconds=native,observed_native_plus_assessment_mean_seconds=observed,native_only_conditional_days=None if native is None else native*202010/86400,observed_loop_lower_bound_days=None if observed is None else observed*202010/86400,additional_costs='checkpoint/ready/wait/native commit, startup, strict read, future witness growth not added to inclusive nested timers; no duration admission'),scope='finite four-stage support, not unrestricted duration-bounded monitoring; application/experimental drift budgets unresolved'))
    figures(out,timeline,steps)
    failed=next((v for v in timeline if v['quality']!='PASS'),None)
    if failed:
        from run_wj_h16 import halt
        halt(root,'coverage_halt',dict(reason='retrospective final-union qualification failed',earliest_affected_step=failed['step'],earliest_affected_time=failed['time'],classification=failed['quality'],native_and_live_decisions_preserved=True))
    assert all(v['quality']=='PASS'for v in timeline),'final union retrospectively failed; preserve native/live decisions'

def figures(out,timeline,steps):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t=[v['time']for v in timeline];fig,axes=plt.subplots(1,2,figsize=(9,3.3),constrained_layout=True)
    axes[0].plot(t,[v['min_in_plane']for v in timeline],'o-',label='all-owner / Gauss8 minimum');axes[0].axhline(.002068763923454,color='tab:red',ls='--',label='half-reference quality gate');axes[0].set(xlabel='T [Hh/U∞]',ylabel='m / νeff');axes[0].legend(fontsize=8)
    axes[1].plot(t,[v['U_eq_L2']for v in timeline],'o-',label='current − analytical reference');axes[1].plot(t,[v['U_drift_L2']for v in timeline],'s-',label='current − original qh0');axes[1].set(xlabel='T [Hh/U∞]',ylabel='Physical full-domain U L2');axes[1].legend(fontsize=8)
    fig.suptitle('Fixed 6 mm neutral control: first step is H14 replay');fig.savefig(out/'margins-and-drift.png',dpi=160);plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.5,3.2),constrained_layout=True)
    ax.plot([v['step']for v in steps],[v['native_compute_seconds']for v in steps],'o-',label='native compute before publication / wait');ax.plot([v['step']for v in steps],[v['independent_seconds']for v in steps],'s-',label='independent assessment / publication');ax.set(xlabel='Absolute step (1 = replay)',ylabel='Wall seconds',xticks=[v['step']for v in steps]);ax.legend(fontsize=8);fig.savefig(out/'per-stage-cost.png',dpi=160);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();analysis(a.root)
