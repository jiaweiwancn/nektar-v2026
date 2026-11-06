#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Scientific source/coupled-field plots; uses immutable readbacks only."""
import argparse
import math
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_cls_epawf_picard import load


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sources',type=Path,required=True);p.add_argument('--uniform-analysis',type=Path,required=True)
    p.add_argument('--upstream',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False);plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25})
    def save(fig,name):
        fig.tight_layout();fig.savefig(a.output/(name+'.png'),dpi=160,bbox_inches='tight')
        fig.savefig(a.output/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
    def points(name,variant=-1,x=1.6,top=.01):
        s=load(a.sources/name/'result.json')
        return sorted([v for v in s['points'] if v.get('valid') and v['variant']==variant and abs(v['x']-x)<1e-12 and v['z']<=top],key=lambda v:v['z'])
    fig,axes=plt.subplots(1,3,figsize=(12,4.4))
    q=points('upstream-h0');ep=points('upstream-h0',1);later=points('upstream-T.02')
    for key,label in [('production','represented P'),('reference_production','analytical P_A0'),('epsilon_A','original epsilon_A')]:
        axes[0].plot([v[key] for v in q],[v['z'] for v in q],'.-',label=label)
    axes[0].set(xscale='log',xlabel='production / destruction [k / time]',title='Initial source balance, x=1.60')
    for rows,label in [(q,'h0, selected chi=0'),(ep,'same h0, chi=1 source only'),(later,'T=.02, selected chi=0')]:
        axes[1].plot([v['source'] for v in rows],[v['z'] for v in rows],'.-',label=label)
    axes[1].set_xscale('symlog',linthresh=1.);axes[1].set(xlabel='P - D [k / time]',title='Source-only counterfactual')
    for rows,label,key in [(q,'h0: represented / analytic','production'),(later,'T=.02: linear / analytic','linear_production_same_gradient'),(later,'T=.02: saturated / analytic','production')]:
        axes[2].plot([v[key]/v['reference_production'] for v in rows],[v['z'] for v in rows],'.-',label=label)
    axes[2].set(xscale='log',xlabel='production / analytical P_A0',title='Gradient error and saturation')
    for ax in axes:
        ax.set_yscale('symlog',linthresh=1e-6);ax.set_ylim(0,.01);ax.set_ylabel('z [length]')
        ax.axhline(q[0]['z_epsilon'],color='black',ls=':',lw=.8);ax.legend(fontsize=7)
    save(fig,'source-balance')
    u=load(a.uniform_analysis/'summary.json');names=[n for n in ('initial','step8','step16','step32','step64','step160') if n in u['states']]
    old_paths={'initial':a.upstream/'analysis16/initial','step8':a.upstream/'analysis16/step8','step16':a.upstream/'analysis16/step16','step32':a.upstream/'monitor/step-000032/step32'}
    fig,axes=plt.subplots(1,3,figsize=(12,4.5));colors=plt.get_cmap('viridis')
    for i,name in enumerate(names):
        time=u['states'][name]['time'];color=colors(i/max(1,len(names)-1));uniform_source='uniform-h0' if name=='initial' else 'uniform-'+name
        if (a.sources/uniform_source/'result.json').exists():
            rows=points(uniform_source)
            axes[0].plot([v['local_k'] for v in rows],[v['z'] for v in rows],'.-',color=color,label=f'uniform T={time:g}')
        if name in old_paths:
            src='upstream-h0' if name=='initial' else ('upstream-T.02' if name=='step32' else None)
            if src and (a.sources/src/'result.json').exists():
                rows=points(src);axes[0].plot([v['local_k'] for v in rows],[v['z'] for v in rows],'--',color=color,label=f'upstream T={time:g}')
        rows=[v for v in u['states'][name]['profiles'] if abs(v['coordinate'][0]-1.6)<1e-12]
        axes[1].plot([v['u_v_k_outer_epsilon'][0] for v in rows],[v['coordinate'][2] for v in rows],'.-',color=color,label=f'T={time:g}')
        ep=[v for v in rows if len(v['u_v_k_outer_epsilon'])==4]
        axes[2].plot([v['u_v_k_outer_epsilon'][3] for v in ep],[v['coordinate'][2] for v in ep],'.-',color=color)
    axes[0].set(xlabel='transported k [velocity squared]',ylabel='z',title='Near-ground k at x=1.60')
    axes[0].set_yscale('symlog',linthresh=1e-6);axes[0].set_ylim(0,.01);axes[0].legend(fontsize=6)
    axes[1].set(xlabel='u [velocity]',ylabel='z',title='Uniform EP retained velocity');axes[1].legend(fontsize=7)
    axes[2].set(xlabel='epsilon [k / time]',ylabel='z',title='Uniform EP outer dissipation')
    save(fig,'control-profiles')
    original=[load(a.upstream/'analysis16/summary.json'),load(a.upstream/'monitor/step-000032/summary.json')]
    old_states={};old_budget={}
    for s in original:
        old_states.update(s['states']);old_budget.update({v['step']:v for v in s['budgets']})
    horizon=max(v['time'] for v in u['states'].values());old_states={k:v for k,v in old_states.items() if v['time']<=horizon+1e-13}
    old_budget=[v for k,v in sorted(old_budget.items()) if v['time']<=horizon+1e-13]
    fig,axes=plt.subplots(2,3,figsize=(12,7))
    for states,budget,label,style in [(old_states,old_budget,'Upstream','--'),(u['states'],u['budgets'],'Uniform EP','-')]:
        ss=sorted(states.values(),key=lambda v:v['time']);t=[v['time'] for v in ss];bt=[v['time'] for v in budget]
        axes[0,0].plot([0.]+bt,[ss[0]['scales_U_K_E_H_L_volume'][1]]+[v['ground']['minimum'] for v in budget],style,label=label+' ground')
        if label=='Uniform EP':axes[0,0].plot(bt,[v['minimum']['value'] for v in budget],':',label=label+' retained bound rows')
        for key,tag in [('global_k_integral','global'),('near_k_integral','z<=.01')]:
            axes[0,1].plot(t,[v[key]/ss[0][key] for v in ss],style,marker='.',label=label+' '+tag)
        axes[0,2].plot(t,[v['mean_wall_stress']/.16 for v in ss],style,marker='.',label=label)
        for index,tag in [(0,'U'),(2,'k')]:
            axes[1,0].semilogy(t,[v['spatial_block_norms'][index] for v in ss],style,marker='.',label=label+' '+tag)
        axes[1,1].plot(bt,[v['signed_integral_rate'] for v in budget],style,label=label+' net')
        axes[1,1].plot(bt,[v['absolute_integral']/.000625 for v in budget],style,label=label+' absolute',alpha=.6)
        axes[1,2].plot(bt,[v['raw_increment_mass_norm'] for v in budget],style,label=label+' raw')
        axes[1,2].plot(bt,[v['correction_mass_norm'] for v in budget],style,label=label+' repair',alpha=.6)
    for ax,title,y in zip(axes.flat,['Ground and retained-row minima','Integrated k / own h0','Mean ground stress','Original spatial residual, alpha=0','Correction supply rate','Raw k update and correction'],
                          ['k','ratio','tau / u_star squared','native scaled block norm','integral(k) / time','physical mass norm']):
        ax.set(xlabel='physical time',ylabel=y,title=title);ax.legend(fontsize=6)
    save(fig,'comparison-budgets')


if __name__=='__main__':main()
