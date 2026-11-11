#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Scientific WJ-7 finite-horizon figures from compact physical tables."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wj_immutable import write_once,artifacts


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--tables',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    def rows(name):return list(csv.DictReader((a.tables/(name+'.csv')).open()))
    def save(fig,name):
        fig.tight_layout()
        for ext in ('png','pdf','svg'):fig.savefig(a.output/(name+'.'+ext),dpi=180)
        plt.close(fig)
    data=rows('interfaces');end=max(float(v['time']) for v in data)
    title='Retained Tx requalification; T3 NOT_RUN' if end<.004095 else 'Finite T3 strip qualification'
    fig,ax=plt.subplots(figsize=(7,4));fig.suptitle(title)
    for policy in ('EP','Budget'):
        for tag in ('coarse','fine','baseline'):
            r=sorted([v for v in data if v['policy']==policy and v['tag']==tag and float(v['order'])==8],key=lambda v:float(v['time']))
            ax.plot([float(v['time']) for v in r],[float(v['normalized_strong_divergence']) for v in r],marker='.',ls='-' if policy=='EP' else '--',label=policy+' '+tag)
    ax.axhline(.001,color='black',ls=':',label='fixed budget');ax.set(xlabel='Physical time [s] (T3=0.004096; T3 u*/H=0.002048)',ylabel='Sampled max |div U| H/u*');ax.legend(fontsize=7,ncol=2);ax.grid(alpha=.2);save(fig,'divergence')
    data=rows('sensitivity');fig,ax=plt.subplots(figsize=(7,4));fig.suptitle(title)
    for policy in ('EP','Budget'):
        for group in ('horizontal','timestep'):
            endpoints=sorted({int(float(v['qualification_end'])) for v in data if v['policy']==policy})
            values=[max((float(v['value'])+float(v['uncertainty']))/float(v['limit']) for v in data if v['policy']==policy and v['group']==group and int(float(v['qualification_end']))==n) for n in endpoints]
            ax.plot([n*3.2e-5 for n in endpoints],values,marker='o',label=policy+' '+group)
    ax.axhline(1,color='black',ls=':');ax.set(xlabel='Qualified prefix endpoint [s]',ylabel='Worst (difference + both uncertainties) / limit');ax.grid(alpha=.2);ax.legend(fontsize=8);save(fig,'matched-sensitivity')
    data=rows('fields');fig,axes=plt.subplots(1,2,figsize=(10,4));fig.suptitle(title)
    for policy in ('EP','Budget'):
        for kind,ax in zip(('continuum_error','initial_drift'),axes):
            r=sorted([v for v in data if v['policy']==policy and v['tag']=='fine' and float(v['order'])==8 and v['kind']==kind and v['scope']=='full' and v['field']=='velocity'],key=lambda v:float(v['time']))
            ax.plot([float(v['time']) for v in r],[float(v['normalized_L2']) for v in r],marker='.',label=policy)
            ax.axhline(1e-4,color='black',ls=':');ax.set(xlabel='Physical time [s]',ylabel='Normalized physical L2',title=kind.replace('_',' '));ax.grid(alpha=.2);ax.legend()
    save(fig,'velocity-error-drift');write_once(a.output/'figures.json',dict(files=artifacts(a.output.rglob('*')),scope='finite Dirichlet strip; no stationarity claim'))


if __name__=='__main__':main()
