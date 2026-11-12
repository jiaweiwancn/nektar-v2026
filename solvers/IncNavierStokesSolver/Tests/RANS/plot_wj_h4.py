#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H4 finite-inventory minima and both physical normal traces; no solves."""
import argparse
import csv
import gzip
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_wj_h4 import save
from wj_h3_analysis import identity


def main(root):
    source=root/'publication-v1';out=root/'figures-v2';out.mkdir()
    with (source/'principal-history.csv').open()as f:rows=list(csv.DictReader(f))
    names=('base-x1','base-x1-half','refined-x1','refined-x1-half','base-x3','base-x3-half','local-x3','local-x3-half')
    assert all(np.isfinite(float(v['margin']))for v in rows)
    plt.rcParams.update({'font.size':8,'figure.dpi':160})
    fig,axes=plt.subplots(1,2,figsize=(9,3.7))
    for name in names:
        rr=sorted((v for v in rows if v['case']==name),key=lambda v:float(v['time']))
        accepted=[v for v in rr if v['role']=='accepted'];rejected=[v for v in rr if v['role']=='rejected']
        if not accepted:continue
        for ax in axes:
            line=ax.plot([float(v['time'])for v in accepted],[float(v['margin'])for v in accepted],'-o',ms=2,label=name)[0]
            ax.plot([float(v['time'])for v in rejected],[float(v['margin'])for v in rejected],'x',color=line.get_color(),ms=8,mew=1.5)
    for ax in axes:ax.axhline(0,color='k',lw=.7);ax.set(xlabel='T=t U∞/Hh',ylabel='principal / νeff at each state’s own minimum')
    axes[0].legend(ncol=2,fontsize=6);axes[1].set_ylim(-.0005,.002)
    fig.suptitle('Expanded finite inventory; crosses are rejected candidates; locations can differ')
    fig.tight_layout();fig.savefig(out/'principal-history.png');plt.close(fig)
    with gzip.open(source/'paired-declared.csv.gz','rt')as f:paired=list(csv.DictReader(f))
    traces={};wanted={'refined-x1-half','base-x3-half','local-x3-half'}
    for v in paired:
        if float(v['x'])!=3. or float(v['z'])!=1.1947077409162705:continue
        for side in ('a','b'):
            case=v['case_'+side]
            if case not in wanted:continue
            key=(case,float(v['intended_time']),v['normal_side_'+side])
            value=(float(v['principal_'+side]),float(v['Gxz_'+side]),float(v['tau_'+side]),v['native_role_'+side])
            if key in traces:assert traces[key]==value,'inconsistent duplicate physical observation'
            traces[key]=value
    fig,axes=plt.subplots(1,3,figsize=(10,3.8));colors=dict(zip(sorted(wanted),('C0','C1','C2')))
    for case in sorted(wanted):
        for side in ('interior','upper','lower'):
            rr=sorted((t,v)for(c,t,s),v in traces.items()if c==case and s==side)
            if not rr:continue
            label=case+': '+{'interior':'interior','upper':'below face','lower':'above face'}[side]
            style='--'if side=='lower'else '-'
            for j,ax in enumerate(axes):
                ax.plot([t for t,v in rr],[v[j]for t,v in rr],style,color=colors[case],lw=1,label=label)
                ax.plot([t for t,v in rr if v[3]=='rejected'],[v[j]for t,v in rr if v[3]=='rejected'],'x',color=colors[case],ms=7)
    for ax in axes:ax.set_xlabel('T=t U∞/Hh')
    axes[0].axhline(0,color='k',lw=.6);axes[0].set_ylabel('principal / νeff')
    axes[1].set_ylabel('Gxz Hh/U∞');axes[2].set_ylabel('(k/epsilon) U∞/Hh')
    axes[0].legend(fontsize=5.5);fig.suptitle('x=3, physical z=1.1947077409162705; dt=.005; all normal owners')
    fig.tight_layout();fig.savefig(out/'face-traces.png');plt.close(fig)
    save(out/'manifest.json',dict(source=identity(__file__),inputs=[identity(source/'principal-history.csv'),identity(source/'paired-declared.csv.gz')],
        outputs=[identity(p)for p in out.glob('*.png')],role='versioned figures explicitly showing rejected candidates; no governing or diagnostic solve'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);main(p.parse_args().root)
