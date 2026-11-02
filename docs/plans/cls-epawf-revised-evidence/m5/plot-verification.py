#!/usr/bin/env python3
"""Plot small scientific M5 evidence from a completed Linux qualification."""
import argparse,json,os
os.environ.update({name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser();p.add_argument('verification',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
data=json.loads(a.verification.read_text());assert data['verification_status']=='passed' and data['milestone_closure']
a.output.mkdir(parents=True,exist_ok=True);sections=data['sections']
fig,axs=plt.subplots(2,2,figsize=(10.5,7.6),constrained_layout=True)
for case in sections['manufactured']:
 x=[r['degree'] for r in case['sequence']]
 axs[0,0].semilogy(x,[max(r['error']['epsilon_rms'],1e-18) for r in case['sequence']],'o-',label=case['family'])
 if case['family']!='feet':axs[0,1].semilogy(x,[max(r['error']['pressure_rms'],1e-18) for r in case['sequence']],'o-',label=case['family'])
axs[0,0].set(title='Outer dissipation: fixed physical geometry',xlabel='Native velocity/scalar degree p',ylabel='Physical RMS epsilon error (m²/s³)')
axs[0,1].set(title='Normal-momentum pressure solution',xlabel='Native velocity/scalar degree p',ylabel='Regional pressure RMS norm (m²/s²)')
levels=sections['patch_partition']['quadrature_refinement']
axs[1,0].semilogy([r['outer_order'] for r in levels],[max(r['maximum_weak_moment_difference'],1e-18) for r in levels],'o-')
axs[1,0].axhline(1e-10,color='k',linestyle='--',linewidth=.8,label='Acceptance 1e-10')
axs[1,0].set(title='Same field on one and two patch columns',xlabel='Outer Gauss order per axis',ylabel='Max. scaled common weak-load error')
for case in sections['manufactured']:
 x=[r['degree'] for r in case['sequence']]
 axs[1,1].semilogy(x,[max(r['error']['strong_divergence_rms'],1e-18) for r in case['sequence']],'o-',label=case['family'])
axs[1,1].set(title='Strong divergence diagnostic',xlabel='Native velocity/scalar degree p',ylabel='Physical RMS div(U) (1/s)')
for ax in (axs[0,0],axs[0,1],axs[1,1]):ax.set_xticks([3,4,5,6])
axs[1,0].set_xticks([12,24,48]);axs[1,0].set_xticklabels(['12','24','48'])
for ax in axs.flat:ax.grid(True,which='both',alpha=.2);ax.legend(fontsize=8)
fig.suptitle('M5 conservative mapped CG verification · '+data['source_commit'][:12],fontsize=13)
for ext in ('png','pdf'):fig.savefig(a.output/('curved-verification.'+ext),dpi=180)
print(a.output/'curved-verification.png')
