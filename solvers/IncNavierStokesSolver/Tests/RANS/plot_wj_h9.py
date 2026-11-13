#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Final figure QA from sealed scalar records, charged to H9's closure reserve."""
import argparse
import csv
import json
from pathlib import Path
import time


def plot(root,out):
    start=time.monotonic()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from study_wj_h9 import identity
    from run_wj_h9 import load,save
    source=root/'publication-v3';book=load(root/'closure/accounting-final.json')
    assert book['closed']and book['closure_charge_seconds']==90
    out.mkdir();impulses=load(source/'signed-impulses.json');moments=load(source/'all-Q1-moments.json')
    with(source/'flow-pressure-increments.csv').open()as f:flow=list(csv.DictReader(f))
    fig,axs=plt.subplots(1,3,figsize=(11.5,3.6),constrained_layout=True)
    for b in ('coarse','fine','quarter'):
        dt=impulses[b]['dt'];label=f'dt={dt:g}'
        # Preserve owner samples, but connect only in physical-height order.
        samples=sorted((v['physical'][2]*40,v['value'])for v in impulses[b]['profiles']if v['physical'][0]==8 and v['owner']['region']=='outer'and v['physical'][2]<.7)
        axs[0].plot([v for z,v in samples],[z for z,v in samples],'.-',label=label)
        r=[v for v in flow if v['branch']==b];axs[1].semilogy([float(v['time'])for v in r],[float(v['pressure_increment_L2'])for v in r],'.-',label=label)
        m=[v for v in moments if v['branch']==b]
        color=axs[2].semilogy([v['time']for v in m],[v['regions']['5']['C_projection_L2']for v in m],'.-',label=label+' C')[0].get_color()
        axs[2].semilogy([v['time']for v in m],[v['regions']['5']['F_projection_L2']for v in m],'x--',color=color,label=label+' F')
    axs[0].set(xlabel='first-step j = dt*(piF-piC), at x=8',ylabel='physical z (mm)')
    axs[1].set(xlabel='attempt time (Hh/Uinf)',ylabel='pressure-increment physical L2')
    axs[2].set(xlabel='attempt time (Hh/Uinf)',ylabel='outer continuous-Q1 projection L2')
    for ax in axs:ax.grid(alpha=.2);ax.legend(fontsize=7)
    fig.savefig(out/'startup-impulse.png',dpi=160);plt.close(fig)
    save(out/'receipt.json',dict(role='scalar-only final figure QA; sorted physical coordinates, logarithmic pressure-increment axis',
        original_figure=identity(source/'startup-impulse.png'),source_records=[identity(source/n)for n in ('signed-impulses.json','all-Q1-moments.json','flow-pressure-increments.csv')],
        closed_accounting=identity(root/'closure/accounting-final.json'),charge='included in fully charged90s closure reserve; no new science reservation or ledger reset',
        seconds=time.monotonic()-start,new_governing=0,new_CP0=0,new_strict=0,new_builds=0,source=identity(__file__)))
    save(out/'seal.json',dict(files=[identity(p)for p in sorted(out.iterdir())if p.is_file()]))
    print('PASS final scalar-only figure QA',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();plot(a.root,a.output)
