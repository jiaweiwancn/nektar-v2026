#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete physical balance readback and export compact scientific figures."""
import argparse
import json
import math
from pathlib import Path
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from run_cls_epawf_picard import load,run,write_json,require,state_hash
from compare_cls_epawf_pimple import differences


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.campaign=a.campaign.resolve();a.prior=a.prior.resolve();a.probe=a.probe.resolve()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    summary=load(a.campaign/'analysis-final/summary.json');budget=[v for v in summary['budgets'] if v['accepted']]
    state_metrics={0:summary['states']['h0'],4:summary['states']['contact0025'],8:summary['states']['contact005'],12:summary['states']['later']}
    prior=load(a.prior/'fixed-horizon/bound-0/case.evolution.json');progress=summary['progress']
    primary=load(Path(progress['segment'])/'case.evolution.json');finer=load(a.campaign/'finer-from-h0/case.evolution.json')
    rows={round(v['time']/.000625):v for v in prior['pimple_attempts'] if v['time']<=.0025+1e-13}
    for line in (a.campaign/'main/steps.jsonl').read_text().splitlines():
        v=json.loads(line);rows[v['global_step']]=v['attempt']
    def probe(name,state,time,operation='contact',xml=None):
        folder=a.output/name;folder.mkdir();shutil.copy2(xml or a.prior/'fixed-horizon/bound-0/case.xml',folder/'case.xml')
        write_json(folder/'input.json',dict(operation=operation,state=state,history=state,time=time,alpha=0.,
            inner_stress_model='SaturatedR2',require_history_h0=False))
        receipt=run([str(a.probe),'case.xml','input.json','result.json'],folder,'probe.log',1800)
        require(receipt['exit_code']==0,'balance readback failed');write_json(folder/'receipt.json',receipt)
        return load(folder/'result.json')
    for step in range(1,13):
        if step not in state_metrics:
            packet=load(Path(rows[step]['k_snapshot']))
            state_metrics[step]=probe('step-'+str(step),packet['corrected'],packet['time'])
    scal=summary['states']['h0']['scales_U_K_E_H_L_volume'];K0=scal[1];volume=scal[-1];reference_norm=K0*math.sqrt(volume)
    for v in budget:
        v['correction_over_reference_field_norm']=v['bound_norm']/reference_norm
        v['correction_over_current_integrated_k']=v['signed_integral']/v['global_k']
        v['near_correction_over_current_near_k']=v['near_integral']/v['near_k']
    windows=[]
    for lo,hi in ((0,4),(4,8),(8,12)):
        b=[v for v in budget if lo<v['step']<=hi];duration=(hi-lo)*.000625
        physical=[0.]*4;flux=[0.,0.];rate=[0.,0.]
        for step in range(lo+1,hi+1):
            for source,weight in ((state_metrics[step-1],.000625/2),(state_metrics[step],.000625/2)):
                for i,x in enumerate(source['physical_production_destruction_global_near']):physical[i]+=weight*x
                flow=source['physical_outward_advection_diffusion_near']
                require(flow['success'],'near control flux diagnostic outside domain')
                for i,x in enumerate(flow['outward_advection_diffusion']):flux[i]+=weight*x
                for term in source['component_integral_rates']:
                    rate[0]+=weight*term['global'];rate[1]+=weight*term['near']
        net=sum(v['signed_integral'] for v in b);near=sum(v['near_integral'] for v in b)
        windows.append(dict(start=lo*.000625,end=hi*.000625,correction=net,near_correction=near,
            positive=sum(v['positive_negative'][0] for v in b),negative=sum(v['positive_negative'][1] for v in b),
            correction_rate=net/duration,near_correction_rate=near/duration,
            global_delta_k=state_metrics[hi]['global_k_integral']-state_metrics[lo]['global_k_integral'],
            near_delta_k=state_metrics[hi]['near_k_integral']-state_metrics[lo]['near_k_integral'],
            integrated_physical_P_D_global_near=physical,integrated_gradient_flux_advection_diffusion_near=flux,
            integrated_mass_inverted_spatial_rate_global_near=rate,
            time_quadrature='trapezoid at every accepted endpoint; flux is represented-gradient flux, not an exact locally conservative CG numerical flux'))
    # The two full trajectories are compared on matching accepted times.
    comparisons=[]
    for coarse_step in (8,11):
        coarse_packet=load(Path(rows[coarse_step]['k_snapshot']));time=coarse_step*.000625
        fine_state=finer['accepted_states'][2*coarse_step]
        x=probe('coarse-samples-'+str(coarse_step),coarse_packet['corrected'],time,'samples')
        y=probe('fine-samples-'+str(coarse_step),fine_state,time,'samples',a.campaign/'finer-from-h0/case.xml')
        x['initial']=y['initial']=prior['accepted_states'][0]
        comparisons.append(dict(time=time,physical_L2=differences(x,y)))
    fine_end=probe('fine-end',finer['state'],finer['time'],xml=a.campaign/'finer-from-h0/case.xml')
    physical_initial=probe('initial-samples',prior['accepted_states'][0],0.,'samples')
    physical_ref=probe('reference-samples',prior['fixed_reference'],0.,'samples')
    physical_initial['initial']=physical_ref['initial']=prior['accepted_states'][0]
    references=dict(h0_sha256=state_hash(prior['accepted_states'][0]),reference_sha256=state_hash(prior['fixed_reference']),
        maximum_coefficient_difference=max(abs(x-y) for x,y in zip(prior['accepted_states'][0],prior['fixed_reference'])),
        physical_L2=differences(physical_initial,physical_ref))
    result=dict(status='completed',windows=windows,matched_finer_comparisons=comparisons,reference_comparison=references,
        finer_final=fine_end,budgets=budget,state_metrics=state_metrics,
        finer_qualification_events=[{k:v for k,v in e.items() if k not in ('state','history','audit')} for e in finer['physical_range_events']],
        finer_final_audit=finer['state_audits'][-1])
    write_json(a.output/'summary.json',result)
    plt.rcParams.update({'font.size':9,'axes.grid':True,'grid.alpha':.25})
    def save(fig,name):
        fig.tight_layout();fig.savefig(a.output/(name+'.png'),dpi=170);fig.savefig(a.output/(name+'.pdf'));plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(10,4))
    for name,color in (('reference','0.65'),('h0','black'),('contact0025','tab:blue'),('contact005','tab:orange'),('later','tab:red')):
        rows_profile=[v for v in summary['states'][name]['profiles'] if v['coordinate'][0]==1.6]
        z=[v['coordinate'][2] for v in rows_profile];u=[v['u_v_k_outer_epsilon'][0] for v in rows_profile]
        k=[v['u_v_k_outer_epsilon'][2] for v in rows_profile]
        label=name if name in ('h0','reference') else 't='+format(summary['states'][name]['time'],'.4g')
        style='--' if name=='h0' else '-'
        axes[0].plot(u,z,style,color=color,label=label)
        axes[1].plot([x/K0 for x,r in zip(k,z) if r<=.01],[r for r in z if r<=.01],style,color=color)
        ep=[(v['u_v_k_outer_epsilon'][3],v['coordinate'][2]) for v in rows_profile if len(v['u_v_k_outer_epsilon'])==4]
        axes[2].plot([x for x,z in ep],[z for x,z in ep],style,color=color)
    for ax,label in zip(axes,('u','k / K0','outer epsilon')):ax.set_xlabel(label);ax.set_ylabel('z')
    axes[0].legend(fontsize=8);axes[1].set_title('First wall-normal cell');fig.suptitle('Retained profiles at x=1.60; primary run rejected at t=.008125',fontsize=10)
    save(fig,'profiles')
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    steps=sorted(state_metrics);tt=[n*.000625 for n in steps];g0=state_metrics[0]['global_k_integral'];n0=state_metrics[0]['near_k_integral']
    axes[0,0].plot(tt,[state_metrics[n]['global_k_integral']/g0 for n in steps],label='global k / initial')
    axes[0,0].plot(tt,[state_metrics[n]['near_k_integral']/n0 for n in steps],label='near k / initial');axes[0,0].legend()
    t=[v['time'] for v in budget]
    for key,label in (('raw_norm','raw step'),('corrected_norm','corrected step'),('bound_norm','bound')):
        axes[0,1].plot(t,[v[key]/reference_norm for v in budget],label=label)
    axes[0,1].set_ylabel('physical mass norm / (K0 sqrt(volume))');axes[0,1].legend()
    positive=negative=0.;pp=[];nn=[]
    for v in budget:positive+=v['positive_negative'][0];negative+=v['positive_negative'][1];pp.append(positive/n0);nn.append(negative/n0)
    axes[1,0].plot(t,[v['cumulative_accepted_integral']/n0 for v in budget],label='net')
    axes[1,0].plot(t,pp,label='positive');axes[1,0].plot(t,nn,label='negative');axes[1,0].set_ylabel('cumulative correction / initial near k');axes[1,0].legend()
    for i,label in ((0,'u'),(2,'k'),(3,'epsilon')):
        axes[1,1].semilogy([summary['states'][n]['time'] for n in ('h0','contact0025','contact005','later')],
            [summary['states'][n]['spatial_block_norms'][i] for n in ('h0','contact0025','contact005','later')],marker='o',label=label)
    axes[1,1].set_ylabel('original alpha=0 scaled residual');axes[1,1].legend()
    for ax in axes.flat:ax.set_xlabel('physical time')
    save(fig,'evolution-budget')
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    for ax,name in zip(list(axes.flat)[:3],('pre_k_first','contact0025','contact005')):
        rows_scan=[v for v in summary['states'][name]['scan'] if v['variant']=='complete'];hs=[v['h'] for v in rows_scan]
        for key,label in (('raw_rate','raw / h'),('bound_rate','bound / h'),('net_rate','net / h')):
            ax.semilogx(hs,[v['points'][0][key] for v in rows_scan],marker='.',label=label)
        ax.axhline(summary['states'][name]['contacts'][0]['mass_inverted_rate'],ls=':',color='black',label='independent f_k')
        ax.set_title(name);ax.set_xlabel('diagnostic h');ax.set_ylabel('rate at selected ground minimum');ax.legend(fontsize=7)
        axes[1,1].loglog(hs,[v['raw_minus_full_spatial_rate_mass_norm'] for v in rows_scan],marker='.',label=name)
    axes[1,1].set_xlabel('diagnostic h');axes[1,1].set_ylabel('mass norm: raw rate minus f_k');axes[1,1].legend(fontsize=8)
    save(fig,'contact-rates')
    print(a.output/'summary.json')


if __name__=='__main__':main()
