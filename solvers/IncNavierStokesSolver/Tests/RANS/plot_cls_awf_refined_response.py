#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Scientific figures from verified immutable refined-response readbacks."""
import argparse,os
from pathlib import Path
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import THREADS,load,require,sha


def physical(row):
    f=Path(row['receipt']['output']);require(sha(f)==row['receipt']['output_sha256'],'plot input changed')
    return load(f)['diagnostics']['physical']


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--retained',type=Path,required=True);p.add_argument('--response',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False)
    import numpy as np
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    d=load(a.retained/'summary.json');fig,axes=plt.subplots(1,3,figsize=(12,4.2),layout='constrained');inputs=[]
    for label in ('P1-h0','half-first-h0','P1-late'):
        q=physical(d['cases'][label]);s=q['common_samples'];x=np.array(s['x']);z=np.array(s['z']);station=min(x,key=lambda v:abs(v+3.5));mask=(x==station)&(z<=.5)
        axes[0].plot(100*(np.array(s['du_dz'])[mask]/np.array(s['reference_du_dz'])[mask]-1),z[mask],label=label)
        r=d['cases'][label]['ground_reaction']['retained_state'];axes[1].plot([label,label],[100*(r['sampled_weighted_mean_shear']/.16-1),100*(r['inferred_weighted_mean_shear']/.16-1)],'-',color='lightgray')
        axes[1].plot(label,100*(r['sampled_weighted_mean_shear']/.16-1),'o',color='tab:blue',label='weighted sampled' if label=='P1-h0' else None)
        axes[1].plot(label,100*(r['inferred_weighted_mean_shear']/.16-1),'s',color='tab:orange',label='weighted spatial reaction' if label=='P1-h0' else None)
        axes[2].semilogy(range(4),q['common_normalized_spatial_first4'],'o-',label=label)
        inputs.append(d['cases'][label]['receipt'])
    axes[0].set(xlabel='dU/dz reference error (%)',ylabel='physical z',ylim=(0,.5),title=f'Native free column x={station:.5g}')
    axes[0].set_yscale('symlog',linthresh=.00125);axes[0].axvline(0,color='k',lw=.5)
    axes[1].set(ylabel='Weighted shear reference error (%)');axes[1].axhline(0,color='k',lw=.5);axes[1].tick_params(axis='x',labelrotation=18)
    axes[2].set_xticks(range(4),['Ux','W','k','epsilon']);axes[2].set(ylabel='Common-normalized selected spatial residual')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(a.output/'representation-reaction.png',dpi=170);plt.close(fig)
    if a.response:
        d=load(a.response/'summary.json');times=sorted(map(float,d['response']));rows=[d['response'][str(t)] for t in times]
        fig,axes=plt.subplots(2,3,figsize=(13,8),layout='constrained');trace=d.get('stepwise_k_response',{}).get('rows',[])
        if trace:axes[0,0].plot([v['elapsed_time'] for v in trace],[v['k_amplification'] for v in trace],lw=1)
        axes[0,0].plot(times,[v['k_amplification'] for v in rows],'o');axes[0,0].axhline(1,color='k',ls=':');axes[0,0].set(ylabel='k amplification relative to actual seed')
        for k in ('u','v','epsilon','pressure'):axes[0,1].semilogy(times[1:],[v['normalized_difference'][k] for v in rows[1:]],'o-',label=k)
        axes[0,1].set(ylabel='Induced field L2 / fixed common scale')
        axes[0,2].semilogy(times[1:],[v['normalized_difference']['k'] for v in rows[1:]],'o-',label='perturbed minus control')
        axes[0,2].semilogy(times[1:],[v['control_drift_from_parent']['k'] for v in rows[1:]],'o-',label='control drift from parent');axes[0,2].set(ylabel='Normalized k L2')
        for t in times:
            c=physical(d['observations']['control'][str(t)]);q=physical(d['observations']['perturbed'][str(t)])
            s=c['common_samples'];x=np.array(s['x']);z=np.array(s['z']);station=min(x,key=lambda v:abs(v+1));mask=(x==station)&(z<=.5)
            delta=(np.array(q['common_samples']['k'])-np.array(s['k']))/q['reference']['K0']
            axes[1,0].plot(delta[mask],z[mask],label=f'elapsed {t:g}')
        axes[1,0].set(xlabel='perturbed-control k / K0',ylabel='physical z',ylim=(0,.5),title=f'Native column x={station:.5g}');axes[1,0].set_yscale('symlog',linthresh=.0025)
        xx=sorted(set(x));zz=sorted(set(z[z<.5]));field=delta.reshape(len(xx),-1)[:,:len(zz)]
        im=axes[1,1].contourf(xx,zz,field.T,levels=20,cmap='RdBu_r');fig.colorbar(im,ax=axes[1,1],label='k difference / K0')
        axes[1,1].set(xlabel='physical x',ylabel='physical z',title=f'Elapsed {times[-1]:g}');axes[1,1].set_yscale('symlog',linthresh=.0025)
        for r,name in [(0,'inner'),(1,'outer')]:axes[1,2].semilogy(times[1:],[v['perturbation_localization']['strong_divergence_difference_L2_inner_outer'][r] for v in rows[1:]],'o-',label=name)
        axes[1,2].set(ylabel='Strong-divergence difference physical L2')
        for ax in (axes[0,0],axes[0,1],axes[0,2],axes[1,2]):ax.set_xlabel('Elapsed physical time')
        for ax in axes.flat:ax.grid(alpha=.2)
        for ax in (axes[0,1],axes[0,2],axes[1,0],axes[1,2]):ax.legend(fontsize=7)
        fig.suptitle('Refined P1 finite-horizon response: fixed reference, retained epsilon, zero bound supply')
        fig.savefig(a.output/'late-response.png',dpi=170);plt.close(fig)
    atomic(a.output/'manifest.json',dict(retained_summary=str(a.retained/'summary.json'),retained_summary_sha256=sha(a.retained/'summary.json'),
        response_summary=str(a.response/'summary.json') if a.response else None,script_sha256=sha(Path(__file__)),
        images={f.name:sha(f) for f in a.output.glob('*.png')},governing_solves=0,physical_commits=0))


if __name__=='__main__':main()
