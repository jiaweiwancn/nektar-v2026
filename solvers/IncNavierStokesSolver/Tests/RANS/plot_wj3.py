#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Standalone scientific figures from compact WJ-3 native receipts."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('summaries',nargs='+',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=True,parents=True)
    data=[json.loads(s.read_text()) for s in a.summaries]
    fig,axes=plt.subplots(1,3,figsize=(10,3.2),layout='constrained')
    colors=['#1769aa','#c45b16','#2a8040']
    for s in data:
        style='-' if s['policy']=='EP' else '--';rows=s.get('spatial',[])
        for j,(name,index) in enumerate([('u',0),('k',2),('epsilon',3)]):
            axes[0].semilogy([r['level'] for r in rows],[r['raw_reference_norms'][index] for r in rows],style+'o',color=colors[j],label=s['policy']+' '+name)
        for f,c in [('u',colors[0]),('k',colors[1]),('outer_epsilon',colors[2])]:
            main=s.get('main',[])
            if main:
                axes[1].loglog([r['dt'] for r in main],[r['comparison']['field_errors']['initial_drift']['fields'][f]['L2'] for r in main],style+'o',color=c,label=s['policy']+' '+f)
            rate=s.get('temporal',{}).get('rates',{}).get(f)
            if rate:
                axes[2].semilogy([0,1],[rate['A'],rate['B']],style+'o',color=c,label=f"{s['policy']} {f} p={rate['p']:.2f}" if rate['p'] is not None else f"{s['policy']} {f} unresolved")
    axes[0].set(xlabel='Normal bisection level (p=3)',ylabel='Spatial residual, scaled mass dual',title='Startup: alpha_mass=0, eta=0',xticks=[0,1,2])
    axes[1].set(xlabel='dt (s)',ylabel='Drift from represented CP0, physical L2',title='Common-time drift; not stage defect')
    durations='; '.join(f"{s['policy']}: T={s.get('frozen_design',{}).get('T',0):.6g}s, T u*/H={s.get('frozen_design',{}).get('T_ustar_over_H',0):.6g}" for s in data)
    axes[2].set(xlabel='Difference pair',ylabel='Physical L2 timestep difference',xticks=[0,1],xticklabels=['D vs D/2','D/2 vs D/4'],title='Resolved short-time comparison')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=6)
    fig.suptitle(durations,fontsize=9)
    for ext in ('svg','pdf','png'):fig.savefig(a.output/('wj3-spatial-temporal.'+ext),dpi=170)
    fig,axes=plt.subplots(1,2,figsize=(8,3.1),layout='constrained')
    for s in data:
        main=s.get('main',[]);style='-o' if s['policy']=='EP' else '--s'
        for field in ('u','k','outer_epsilon'):
            axes[0].loglog([r['dt'] for r in main],[r['comparison']['field_errors']['continuum_error']['fields'][field]['L2'] for r in main],style,label=s['policy']+' '+field)
        for name in ('traction_max','k_flux_max','pressure_jump_error'):
            axes[1].loglog([r['dt'] for r in main],[max(r['comparison']['field_errors']['interface'][name],1e-18) for r in main],style,label=s['policy']+' '+name)
    axes[0].set(xlabel='dt (s)',ylabel='Physical L2 total continuum error',title='Spatial representation error retained')
    axes[1].set(xlabel='dt (s)',ylabel='Interface sampled maximum',title='Separate physical interface diagnostics')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=6)
    fig.suptitle(durations,fontsize=8)
    for ext in ('svg','pdf','png'):fig.savefig(a.output/('wj3-errors-interface.'+ext),dpi=170)

    fig,axes=plt.subplots(2,2,figsize=(8.4,5.5),layout='constrained')
    for axis,index,label in zip(axes.flat[:3],(0,2,3),('u momentum','k','Entire outer epsilon')):
        for s in data:
            color='#1769aa' if s['policy']=='EP' else '#c45b16'
            rows=s.get('main',[])
            axis.loglog([r['dt'] for r in rows],[r['comparison']['spatial_mass_dual_norms'][index] for r in rows],'-o',color=color,label=s['policy']+' spatial alpha=0')
            axis.loglog([r['dt'] for r in rows],[r['comparison']['completed_stage'][index] for r in rows],'--s',color=color,label=s['policy']+' actual completed stage')
        axis.set(xlabel='dt (s)',ylabel='Scaled consistent-mass dual norm',title=label)
        axis.grid(alpha=.2);axis.legend(fontsize=6)
    for s in data:
        rows=s.get('spatial',[])
        axes[1,1].loglog([r['unknowns'] for r in rows],[r['first_step_pressure_adjustment']['fields']['reduced_pressure']['L2'] for r in rows],'-o',label=s['policy'])
    axes[1,1].set(xlabel='Native unknowns',ylabel='First-step reduced-pressure drift, L2',title='Algebraic pressure adjustment; dt=1e-6 s')
    axes[1,1].grid(alpha=.2);axes[1,1].legend(fontsize=6)
    fig.suptitle(durations+'; eta=0',fontsize=8)
    for ext in ('svg','pdf','png'):fig.savefig(a.output/('wj3-residual-stage-pressure.'+ext),dpi=170)

if __name__=='__main__':main()
