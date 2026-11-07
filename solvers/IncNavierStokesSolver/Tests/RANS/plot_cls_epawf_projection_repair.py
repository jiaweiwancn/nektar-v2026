#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Static figures from immutable projection and continuation readbacks."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_cls_epawf_picard import load


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--uniform',type=Path,required=True)
    p.add_argument('--readback',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=True)
    r=load(a.readback/'summary.json');verification=load(a.root/'repair-verification.json')
    trace=[json.loads(x) for x in (a.root/'replay-total20000/trace.jsonl').read_text().splitlines()]
    separated=[v for v in trace if v['kind']=='separation'];solved=[v for v in trace if v['kind']=='solved']
    floor=verification['larger_total']['bound']['floor']
    fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    # The plotting threshold is only a log-axis convention, never a solver floor.
    ax[0,0].semilogy([v['pass'] for v in separated],[max(1e-16,floor-v['ground_minimum']) for v in separated],lw=.8)
    ax[0,0].axvline(1025,color='r',ls='--',label='old budget stop');ax[0,0].legend(fontsize=8)
    ax[0,0].set(xlabel='Physical separation pass',ylabel='Ground floor violation (log display >= 1e-16)')
    ax[0,1].plot([v['total_work'] for v in solved],[1e6*v['last_completed_primal_objective'] for v in solved],lw=.8)
    ax[0,1].axvline(2000,color='r',ls='--');ax[0,1].set(xlabel='Completed KKT solves',ylabel='Mass objective / 1e-6')
    budgets=r['budgets'];t=[v['time'] for v in budgets]
    ax[1,0].plot(t,[v['correction_raw_ratio'] for v in budgets]);ax[1,0].set(xlabel='Physical time',ylabel='Correction / raw step mass norm')
    for name,label in (('signed_integral','Net'),('positive_integral','Positive'),('negative_integral','Negative')):
        ax[1,1].plot(t,np.cumsum([v[name] for v in budgets]),label=label)
    ax[1,1].set(xlabel='Physical time',ylabel='Accumulated correction integral');ax[1,1].legend(fontsize=8)
    for axis in ax.flat:axis.grid(alpha=.25)
    fig.suptitle('Attempt-73 budget repair and retained Upstream evolution')
    fig.savefig(a.output/'projection-and-budget.png',dpi=180);plt.close(fig)

    milestones=[]
    for step in (160,800,1024):
        f=a.uniform/f'uniform-monitor/step-{step:06d}/summary.json'
        if f.exists():
            d=load(f);milestones.append((d['states']['step'+str(step)],d))
    fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    states=[('h0',milestones[0][1]['states']['initial'])]+[(f"T={s['time']:.3g}",s) for s,_ in milestones]
    for label,s in states:
        profile=[v for v in s['profiles'] if abs(v['coordinate'][0]-1.6)<1e-10]
        z=[v['coordinate'][2] for v in profile]
        ax[0,0].plot([v['u_v_k_outer_epsilon'][2] for v in profile],z,'o-',ms=3,label=label)
        ax[0,1].plot([v['u_v_k_outer_epsilon'][0] for v in profile],z,'o-',ms=3,label=label)
    for axis in ax[0]:axis.set_ylim(0,.05);axis.legend(fontsize=8);axis.set_ylabel('Height z');axis.grid(alpha=.25)
    ax[0,0].set_xlabel('k at x=1.60');ax[0,1].set_xlabel('Velocity u at x=1.60')
    t=[s['time'] for _,s in states];near0=states[0][1]['near_k_integral'];stress0=states[0][1]['mean_wall_stress']
    ax[1,0].plot(t,[100*(s['near_k_integral']/near0-1) for _,s in states],'o-',label='First-cell integrated k')
    ax[1,0].plot(t,[100*(s['mean_wall_stress']/stress0-1) for _,s in states],'o-',label='Mean wall stress')
    ax[1,0].set(xlabel='Physical time',ylabel='Change from own initialized h0 (%)');ax[1,0].legend(fontsize=8)
    for i,label in ((0,'u'),(2,'k'),(3,'epsilon')):
        ax[1,1].semilogy(t,[s['spatial_block_norms'][i] for _,s in states],'o-',label=label)
    ax[1,1].set(xlabel='Physical time',ylabel='Original spatial residual block (native units)');ax[1,1].legend(fontsize=8)
    for axis in ax[1]:axis.grid(alpha=.25)
    fig.suptitle('Independent UniformEPFlat control: unchanged frozen runtime, no bounding')
    fig.savefig(a.output/'uniform-continuation.png',dpi=180);plt.close(fig)
    print(a.output)


if __name__=='__main__':main()
