#!/usr/bin/env python3
"""Render qualified M4 data; this script does not run numerical solvers."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('verification',type=Path)
p.add_argument('output',type=Path)
a=p.parse_args()
r=json.loads(a.verification.read_text())
if r['verification_status']!='passed' or not r['milestone_closure']:
    raise SystemExit('A complete passed M4 qualification is required')
a.output.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
labels=('U','W','k','epsilon')
colors=('#1565a8','#c46a15','#2d8054','#924e9d')
markers=('o','s','^','D')
fig,axes=plt.subplots(1,2,figsize=(10.4,4.0),constrained_layout=True)
for ax,key,xkey,title,xlabel in ((axes[0],'unanchored_p_refinement','degree','All-region p refinement','Polynomial degree p'),
                                (axes[1],'unanchored_h_refinement','unknowns','Inner h refinement; outer mesh fixed','Total coupled coefficients')):
    rows=r[key]
    for b,(label,color,marker) in enumerate(zip(labels,colors,markers)):
        ax.semilogy([v[xkey] for v in rows],[max(v['raw_norms'][b],1e-18) for v in rows],marker=marker,color=color,label=label)
    ax.set(title=title,xlabel=xlabel,ylabel='Scaled physical mass-dual defect')
    ax.grid(alpha=.23,which='both');ax.legend(frameon=False)
fig.suptitle('Unanchored EP-AWF reference: chi=1, eta=0')
for extension in ('png','svg'):fig.savefig(a.output/f'reference-refinement.{extension}',dpi=170)
plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(10.4,4.0),constrained_layout=True)
for ax,row in zip(axes,r['manufactured_time_refinement']):
    dt=np.asarray(row['dt']);err=np.asarray(row['scaled_physical_field_errors'])
    for b,(label,color,marker) in enumerate(zip(labels,colors,markers)):
        ax.loglog(dt,err[:,b],marker=marker,color=color,label=label)
    guide=max(err[0,:4])*(dt/dt[0])**row['order']
    ax.loglog(dt,guide,'--',color='0.3',alpha=.7,label=f"dt^{row['order']} guide")
    ax.set(title=f"Native IMEX{row['order']}",xlabel='Timestep (s)',ylabel='Scaled physical L2 state error')
    ax.grid(alpha=.23,which='both');ax.legend(frameon=False)
fig.suptitle('Unanchored manufactured evolution: t=0.3 to 0.7 s')
for extension in ('png','svg'):fig.savefig(a.output/f'temporal-order.{extension}',dpi=170)
plt.close(fig)
print('Rendered M4 plots from source',r['source_commit'])
