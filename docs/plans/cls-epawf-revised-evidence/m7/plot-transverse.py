import hashlib,json,platform
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
source=root/'m7-strict-transverse-002/validation.json'
v=json.loads(source.read_text());assert v['success']
out=root/'m7-figures-003';out.mkdir(exist_ok=False)
s=v['sections'];fig,axes=plt.subplots(1,3,figsize=(12.6,4.1),layout='constrained')
plt.rcParams.update({'pdf.fonttype':42,'font.size':10})
trials=[r for r in s['trials'] if r['operator'].startswith('p7-n8-') and r['order']==2 and r['dt']==.00015625 and r['amplitude']==.0001]
trials=sorted(trials,key=lambda r:float(r['operator'].split('xi')[1]))
x=[float(r['operator'].split('xi')[1]) for r in trials];g=[r['gain'] for r in trials]
axes[0].plot(x,g,'o-',color='#a83131',lw=2);axes[0].axhline(1,color='.6',ls='--',lw=.9)
axes[0].set(xlabel=r'Wavenumber $\xi_x$ [m$^{-1}$]',ylabel=r'Energy ratio $E(T)/E(0)$',title='a  Adverse growth at positive covariance',xticks=x)
axes[0].set_yscale('log');axes[0].set_yticks([1,2,4,8],labels=['1','2','4','8']);axes[0].text(.04,.89,'T = 0.05 s\nOriginal CLS coefficients',transform=axes[0].transAxes,va='top',fontsize=9)
for rec,color in zip(s['temporal'],['#31768d','#94612c']):
 rows=rec['records'];d=np.array([r['dt'] for r in rows]);err=np.array([r['relative_error_to_semidiscrete_exponential'] for r in rows]);order=rec['order']
 axes[1].loglog(d,err,'o-',color=color,label=f'Native IMEX{order}',lw=1.8)
 axes[1].loglog(d,err[-1]*(d/d[-1])**order,':',color=color,lw=1)
axes[1].set(xlabel=r'Time step $\Delta t$ [s]',ylabel='Relative field error against exp(tL)',title='b  Time refinement resolves the growth')
axes[1].legend(frameon=False,fontsize=9)
sp=next(a for a in s['spatial'] if a['xi']==4)['records']
for rows,mark,color,label in [([r for r in sp if r['cells']==2],'o','#6a4c93','p = 3, 5, 7; 2 cells'),([r for r in sp if r['degree']==3],'s','#2a8064','2, 4, 8 cells; p = 3')]:
 rows=sorted(rows,key=lambda r:r['degree']*r['cells']);dof=[r['degree']*r['cells']-1 for r in rows]
 axes[2].semilogy(dof,[r['actual_native_relative_field_error'] for r in rows],mark+'-',color=color,label=label,lw=1.8)
axes[2].set(xlabel='Free wall-normal coefficients',ylabel='Relative native field error',title=r'c  Spatial refinement at $\xi_x=4$')
axes[2].legend(frameon=False,fontsize=8.5)
for ax in axes:
 ax.grid(alpha=.18,which='both');ax.spines[['top','right']].set_visible(False)
fig.suptitle('Unchanged CLS: reduced linearized transverse experiment above the wall layer',fontsize=13)
fig.savefig(out/'m7-transverse.png',dpi=190)
fig.savefig(out/'m7-transverse.pdf')
manifest=dict(schema=1,data_source=str(source),data_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),source_commit=v['source_commit'],compiled_physics_commit=json.loads((root/'m7-attempt-004/runtime-fingerprint.json').read_text())['source_commit'],script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),python=platform.python_version(),matplotlib=matplotlib.__version__,scope='reduced linearized three-component transverse Fourier/native-CG system, not nonlinear terrain flow',figures={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('m7-transverse.*')})
(out/'receipt.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(manifest,indent=2))
