#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Scientific figures for an immutable V2 campaign-prefix analysis."""
import argparse
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_cls_epawf_picard import load


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--summary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    s=load(a.summary);states=sorted(s['states'].values(),key=lambda x:x['time']);b=s['budgets']
    plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25})
    def save(fig,name):
        fig.tight_layout();fig.savefig(a.output/(name+'.png'),dpi=170);fig.savefig(a.output/(name+'.pdf'));plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(10,4))
    for v in states:
        rows=[x for x in v['profiles'] if abs(x['coordinate'][0]-1.6)<1e-12]
        label='h0 / neutral' if v['time']==0 else 'T='+format(v['time'],'.5g')
        axes[0].plot([x['u_v_k_outer_epsilon'][0] for x in rows],[x['coordinate'][2] for x in rows],'.-',label=label)
        near=[x for x in rows if x['coordinate'][2]<=.05]
        axes[1].plot([x['u_v_k_outer_epsilon'][2] for x in near],[x['coordinate'][2] for x in near],'.-',label=label)
        outer=[x for x in rows if len(x['u_v_k_outer_epsilon'])==4]
        axes[2].plot([x['u_v_k_outer_epsilon'][3] for x in outer],[x['coordinate'][2] for x in outer],'.-',label=label)
    for ax,title,x in zip(axes,('Streamwise velocity','Near-ground k','Outer dissipation'),('u','k','epsilon')):
        ax.set(title=title,xlabel=x,ylabel='z')
    axes[0].legend(fontsize=7);fig.suptitle('Retained fields at x=1.60; V2 bounded one-sweep evolution',y=1.03)
    save(fig,'profiles')
    t=[x['time'] for x in states];bt=[x['time'] for x in b]
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    for field,label in (('global_k_integral','global'),('near_k_integral','z <= .01')):
        axes[0,0].plot(t,[x[field]/states[0][field] for x in states],'.-',label=label)
    axes[0,0].set(title='Integrated k',ylabel='ratio to initial integral');axes[0,0].legend()
    axes[0,1].plot(t,[x['mean_wall_stress']/.16 for x in states],'.-');axes[0,1].set(title='Mean ground stress',ylabel='tau / u_star^2')
    for idx,label in ((0,'u'),(2,'k'),(3,'epsilon')):
        initial=states[0]['spatial_block_norms'][idx]
        axes[0,2].semilogy(t,[x['spatial_block_norms'][idx]/initial for x in states],'.-',label=label)
    axes[0,2].set(title='Original spatial residual, alpha=0',ylabel='native block norm / initial norm');axes[0,2].legend()
    for key,label in (('raw_increment_mass_norm','raw'),('correction_mass_norm','bound'),('corrected_increment_mass_norm','accepted')):
        axes[1,0].semilogy(bt,[x.get(key,math.nan) or math.nan for x in b],'.-',label=label)
    axes[1,0].set(title='k increments (inactive bounds omitted)',ylabel='physical mass norm');axes[1,0].legend()
    for key,label in (('signed_integral','signed'),('positive_integral','positive'),('negative_integral','negative'),('absolute_integral','absolute')):
        axes[1,1].plot(bt,[x[key]/.000625 for x in b],'.-',label=label)
    axes[1,1].set(title='Correction integrals per unit time',ylabel='integral(k) / time, case units');axes[1,1].legend(fontsize=7)
    axes[1,2].plot(bt,[x['minimum']['value']/x['floor'] for x in b],'.-',label='inventory')
    axes[1,2].plot(bt,[x['ground']['minimum']/x['floor'] for x in b],'--',label='ground')
    axes[1,2].set(title='Protected minima',ylabel='k / imposed floor',yscale='log');axes[1,2].legend()
    for ax in axes.flat:ax.set_xlabel('physical time')
    save(fig,'evolution-budget')


if __name__=='__main__':main()
