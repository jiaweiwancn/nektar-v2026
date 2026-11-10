#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Standalone scientific WJ-4 drift, sensitivity and separate residual figures."""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wj_immutable import read


def main():
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    summaries={policy:read(a.study/policy/'summary.json') for policy in ('EP','Budget')}
    fig,axes=plt.subplots(2,3,figsize=(10,6),layout='constrained')
    for i,policy in enumerate(('EP','Budget')):
        for tag,t in summaries[policy].get('longer',{}).items():
            rows=t['trajectory'];style='-o' if tag=='coarse' else '--s'
            for f in ('velocity','k','outer_epsilon'):
                axes[i,0].semilogy([r['time'] for r in rows],[r['errors']['initial_drift']['full']['fields'][f]['normalized_L2'] for r in rows],style,label=tag+' '+f)
            for name in ('traction','k_flux','pressure_offset'):
                axes[i,1].semilogy([r['time'] for r in rows],[max(r['errors']['interface'][name],1e-18) for r in rows],style,label=tag+' '+name)
            axes[i,2].semilogy([r['time'] for r in rows],[r['spatial_residual'][0] for r in rows],style,label=tag+' spatial u alpha=0')
            axes[i,2].semilogy([r['time'] for r in rows],[r['stage_defect'][0] for r in rows],style,alpha=.5,label=tag+' actual stage u')
        axes[i,0].axhline(1e-4,color='k',ls=':',lw=.8)
        for j,label in enumerate(('Normalized initial-field drift','Normalized interface jumps','Consistent-mass dual residual')):
            axes[i,j].set(xlabel='Physical time (s)',ylabel=label,title=policy);axes[i,j].grid(alpha=.2);axes[i,j].legend(fontsize=6)
    fig.suptitle('Finite reference-Dirichlet strip; frozen T4=0.008192 s, T4 u*/H=0.004096; spatial residual and stage defect are distinct',fontsize=9)
    for ext in ('png','pdf','svg'):fig.savefig(a.output/('wj4-drift-interface-residual.'+ext),dpi=170)
    plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(10,3.5),layout='constrained')
    for policy,s in summaries.items():
        for nx,t in s.get('horizontal',{}).items():
            rows=t['trajectory']
            for scope in ('full','central'):
                axes[0].semilogy([r['time'] for r in rows],[r['errors']['initial_drift'][scope]['fields']['velocity']['normalized_L2'] for r in rows],'-o',label=policy+' Nx'+nx+' '+scope)
            axes[1].plot([r['time'] for r in rows],[r['errors']['strong_divergence'] for r in rows],'-o',label=policy+' Nx'+nx)
        for scope in ('full','central'):
            rows=list(s.get('timestep_comparisons',{}).values())
            axes[2].semilogy([r['time'] for r in rows],[max(r['comparison'][scope]['fields']['velocity']['normalized_L2'],1e-18) for r in rows],'-o',label=policy+' '+scope)
    axes[0].set(ylabel='Normalized velocity drift',title='Horizontal sensitivity, L1')
    axes[1].set(ylabel='max |div U| H/u*',title='Finite observation inventory');axes[1].axhline(1e-3,color='k',ls=':')
    axes[2].set(ylabel='Normalized velocity timestep difference',title='Same-T pair, L2; no rate claim');axes[2].axhline(1e-5,color='k',ls=':')
    for ax in axes:ax.set(xlabel='Physical time (s)');ax.grid(alpha=.2);ax.legend(fontsize=6)
    fig.suptitle('T3=0.004096 s (T3 u*/H=0.002048); T4=0.008192 s (T4 u*/H=0.004096)',fontsize=9)
    for ext in ('png','pdf','svg'):fig.savefig(a.output/('wj4-horizontal-timestep.'+ext),dpi=170)
    plt.close(fig)

if __name__=='__main__':main()
