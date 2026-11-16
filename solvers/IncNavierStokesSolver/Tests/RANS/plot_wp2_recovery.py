#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Plot actual native recovery observations; no interpolated numerical states."""
import argparse,time,resource
from pathlib import Path
from cao_tamura_h2d_state import load
from run_cls_epawf_picard import require,sha
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(exist_ok=False);started=time.monotonic();data={}
    for model in ['CLS','standard']:
        progress_path=a.root/f'preservation-recovery-{model}/progress.json';progress=load(progress_path);rows={};sources=[]
        for entry in progress['results'].values():
            path=Path(entry['report'])
            if not path.exists():continue
            d=load(path);sources.append(dict(path=str(path),sha256=sha(path)))
            for row in d.get('wp2_per_state_monitoring',[]):
                physical=row['physical_checks'];errors=physical['regional_errors_inner_outer'];t=row['time']
                value=dict(actual_time=t,state_sha256=physical['state_sha256'],
                    maximum_pointwise_k_error=max(e['pointwise_normalized_maxima'][2] for e in errors),
                    maximum_unforced_block_norm=max(physical['unanchored_free_block_norms'][:4]),
                    strong_divergence_RMS_over_ustar_H=physical['strong_divergence_RMS_over_ustar_H'])
                key=round(t,12)
                if key in rows:require(rows[key]['state_sha256']==value['state_sha256'],'plot restart boundary differs')
                else:rows[key]=value
        data[model]=dict(trajectory=str(progress_path),trajectory_sha256=sha(progress_path),sources=sources,history=[rows[k] for k in sorted(rows)])
    atomic(a.output/'data.json',data)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(11,3.5),layout='constrained')
    fields=[('maximum_pointwise_k_error','Maximum k reference error',.01),
            ('maximum_unforced_block_norm','Unforced spatial block norm',.001),
            ('strong_divergence_RMS_over_ustar_H','Strong divergence RMS / (u_tau/H_A)',.001)]
    for ax,(field,label,limit) in zip(axes,fields):
        for model,color in [('CLS','#2563a8'),('standard','#c45b12')]:
            values=data[model]['history'];values=[v for v in values if v[field]>0]
            ax.semilogy([v['actual_time'] for v in values],[v[field] for v in values],label='Regularized CLS' if model=='CLS' else 'Calibrated standard',color=color,lw=1.5)
        ax.axhline(limit,color='#333333',ls='--',lw=1,label='Final-state accuracy limit')
        ax.axvline(1.1942387693600194,color='#888888',ls=':',lw=1)
        ax.set(xlabel='Actual physical time (s)',ylabel=label);ax.grid(alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('Flat perturbation recovery, dt = 0.01 s; no bound correction',fontsize=11)
    fig.savefig(a.output/'recovery.png',dpi=180);plt.close(fig)
    atomic(a.output/'provenance.json',dict(source_sha256=sha(Path(__file__)),data_sha256=sha(a.output/'data.json'),image_sha256=sha(a.output/'recovery.png'),
        seconds=time.monotonic()-started,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='native per-step metrics, both regions; connecting lines visualize executed observations and do not define interpolated physical states'))


if __name__=='__main__':main()
