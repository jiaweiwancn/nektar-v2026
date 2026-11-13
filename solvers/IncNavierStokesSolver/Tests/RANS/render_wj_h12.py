#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Presentation-only layout repair from sealed CSV; no constitutive evaluation."""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def render(publication,out):
    with(publication/'signed-window-contributions.csv').open()as f:
        rows=[v for v in csv.DictReader(f)if v['observable']=='mhat'and v['owner_index']=='0']
    assert len(rows)==4
    out.mkdir(exist_ok=False)
    fig,axes=plt.subplots(2,2,figsize=(9,6),sharey=True,layout='constrained')
    for w in range(2):
        for j,pair in enumerate(('coarse-fine','fine-quarter')):
            row=next(v for v in rows if int(v['window'])==w and v['pair']==pair)
            values=[float(row[k])*1e6 for k in ('incoming','lift','flow','k','epsilon','endpoint')]
            ax=axes[w,j];ax.bar(range(6),values,color=['.5','.8','#4477aa','#228833','#ee7733','black'])
            ax.axhline(0,color='.4',lw=.6);ax.set_xticks(range(6),['incoming','lift','flow','k','ε','end'])
            ax.tick_params(labelsize=9);ax.set_title(pair+'; '+('[0,.01]'if w==0 else '[.03,.04]'),fontsize=11)
    fig.supylabel('signed normalized principal difference × 10⁶',fontsize=11)
    fig.suptitle('Geometry 17 upper, x=8, physical z=10.243007 mm\ncoarser minus finer; inherited history remains separate',fontsize=12)
    fig.savefig(out/'signed-window-contributions.png',dpi=170);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--publication',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();render(a.publication,a.out)
