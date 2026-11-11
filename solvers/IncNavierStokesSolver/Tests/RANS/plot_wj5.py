#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Export scientific WJ-5 finite-strip diagnosis/sensitivity figures."""
import argparse
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def save(fig,directory,name):
    fig.tight_layout()
    for suffix in ('png','pdf','svg'):fig.savefig(directory/(name+'.'+suffix),dpi=180,bbox_inches='tight')
    plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--compact',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    load=lambda n:list(csv.DictReader((a.compact/n).open()))
    div=load('divergence.csv');projection=load('pressure-projection.csv');loads=load('near-h-loads.csv');fields=load('fields.csv');diagnoses=json.loads((a.compact/'diagnoses.json').read_text())
    fig,axes=plt.subplots(1,2,figsize=(10,3.7),sharey=True)
    for ax,policy in zip(axes,('EP','Budget')):
        values=[r for r in div if r['policy']==policy and r['field']=='divergence' and r['scope']=='full' and float(r['steps'])>0 and (r['case'].startswith('retained') or r['case'].startswith('control/') or r['case'].startswith('screen/'))]
        labels=[r['case'].replace('retained-','').replace('control/terminal','L1 Nx2 half-dt').replace('screen/','').replace('/terminal','') for r in values]
        ax.bar(range(len(values)),[2*float(r['sampled_Linf']) for r in values],color='#315c8c')
        ax.axhline(1e-3,color='#a33232',ls='--',label='fixed divergence budget')
        ax.set_xticks(range(len(values)),labels,rotation=25,ha='right');ax.set_title(policy+' at Tx = .001024 s');ax.set_ylabel(r'sampled $|\nabla\cdot U|H/u_*$');ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Unchanged reference-Dirichlet strip; Tx u*/H = .000512',fontsize=11)
    save(fig,out,'divergence-screen')
    faces=load('face-traces.csv');fig,axes=plt.subplots(1,2,figsize=(10,3.8),sharey=True)
    for ax,policy in zip(axes,('EP','Budget')):
        names=list(dict.fromkeys(r['case'] for r in faces if r['policy']==policy and float(r['steps'])>0 and (r['case']=='retained-Nx2-32' or r['case']=='control/terminal' or r['case'].startswith('screen/') and r['case'].endswith('/terminal'))))
        for name in names:
            q=[r for r in faces if r['policy']==policy and r['case']==name and float(r['x'])==1. and .075<=float(r['z'])<=.205]
            q.sort(key=lambda r:(float(r['z']),float(r['normal_side'])))
            ax.plot([float(r['z']) for r in q],[2*float(r['divergence']) for r in q],'.-',ms=3,label=name.replace('/terminal',''))
        ax.axvline(.1,color='0.6',ls=':');ax.axhline(1e-3,color='#a33232',ls='--');ax.axhline(-1e-3,color='#a33232',ls='--');ax.set_xlabel('z [m], both actual normal traces');ax.set_ylabel(r'$(G_{xx}+G_{zz})H/u_*$ at x=1');ax.set_title(policy);ax.grid(alpha=.2);ax.legend(fontsize=7)
    save(fig,out,'normal-face-divergence')
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for ax,policy in zip(axes,('EP','Budget')):
        q=[r for r in projection if r['policy']==policy and float(r['steps'])==32 and float(r['region'])==0 and (r['case'].startswith('retained') or r['case'].startswith('screen/'))]
        x=np.arange(len(q));ax.bar(x-.17,[float(r['Q1_projection_L2']) for r in q],.34,label='continuous Q1 projection');ax.bar(x+.17,[float(r['complement_L2']) for r in q],.34,label='orthogonal complement')
        ax.set_yscale('log');ax.set_xticks(x,[r['case'].replace('/terminal','') for r in q],rotation=22,ha='right');ax.set_ylabel('physical divergence L2 [m/s in 2D]');ax.set_title(policy+' inner region at Tx');ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
    save(fig,out,'pressure-test-complement')
    fig,axes=plt.subplots(2,2,figsize=(8.5,6.4));image=None
    for row,policy in enumerate(('EP','Budget')):
        candidates=[key for key in diagnoses if key.startswith(policy+'/screen/') and key.endswith('/terminal')]
        chosen=[policy+'/retained-Nx2-32',candidates[-1] if candidates else policy+'/retained-Nx2-32']
        for col,key in enumerate(chosen):
            energy=np.zeros((4,4))
            for e in diagnoses[key]['element_modes']:
                if e['lower'][2]<.08-1e-13 or e['upper'][2]>.2+1e-13:continue
                vol=(e['upper'][0]-e['lower'][0])*(e['upper'][2]-e['lower'][2])
                energy+=vol*np.array(e['Legendre_x_z'])[2]**2/np.outer([1,3,5,7],[1,3,5,7])
            ax=axes[row,col];image=ax.imshow(np.log10(np.maximum(energy,1e-30)),origin='lower',vmin=-28,vmax=-7,cmap='viridis')
            ax.axhline(1.5,color='white',ls=':');ax.axvline(1.5,color='white',ls=':');ax.set_xticks(range(4));ax.set_yticks(range(4));ax.set_xlabel('normal Legendre degree');ax.set_ylabel('horizontal Legendre degree');ax.set_title(key.replace('/terminal',''),fontsize=9)
    fig.colorbar(image,ax=axes.ravel().tolist(),label='log10 elementwise divergence L2 energy',shrink=.8)
    # A global continuous Q1 projection differs from independent element Q1
    # projections. Low element degrees may also cancel against global tests.
    fig.subplots_adjust(right=.83,wspace=.45,hspace=.5)
    for suffix in ('png','pdf','svg'):fig.savefig(out/('element-divergence-modes.'+suffix),dpi=180,bbox_inches='tight')
    plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for ax,policy in zip(axes,('EP','Budget')):
        keys=[policy+'/retained-Nx2-32']+[key for key in diagnoses if key.startswith(policy+'/screen/') and key.endswith('/terminal')]
        for j,key in enumerate(keys):
            q=sorted(diagnoses[key]['forcing']['samples'],key=lambda r:r['position'][0]);x=[r['position'][0] for r in q]
            for term,style in (('actual_F_x','-'),('mechanical_pressure_gradient_x','--'),('conservative_tangential_convection_x',':')):
                ax.plot(x,[r[term] for r in q],style,color=plt.cm.tab10(j),label=key.split('/',1)[1].replace('/terminal','')+' '+term)
        ax.set_xlabel('x [m], outer one-sided H sampler');ax.set_ylabel('acceleration [m/s²]');ax.set_title(policy+' evolved forcing at Tx');ax.grid(alpha=.2);ax.legend(fontsize=6)
    save(fig,out,'evolved-roughness-forcing')
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    terms=['storage','conservative_advection','viscous_stress','reduced_pressure','roughness_force']
    for ax,policy in zip(axes,('EP','Budget')):
        names=list(dict.fromkeys(r['case'] for r in loads if r['policy']==policy and r['field']=='u' and float(r['steps'])==32 and (r['case'].startswith('retained') or r['case'].startswith('screen/'))))
        x=np.arange(len(terms));width=.8/max(1,len(names))
        for j,name in enumerate(names):
            vals=[next(float(r['value']) for r in loads if r['policy']==policy and r['case']==name and r['field']=='u' and r['term']==term) for term in terms]
            ax.bar(x+(j-(len(names)-1)/2)*width,np.array(vals)*1e6,width,label=name.replace('/terminal',''))
        ax.axhline(0,color='0.4',lw=.7);ax.set_xticks(x,['storage','advection','viscous stress','reduced pressure','roughness force'],rotation=25,ha='right');ax.set_ylabel(r'signed weak load [$10^{-6}$ m$^3$/s$^2$]');ax.set_title(policy+'; continuous test on [.08,.2]');ax.legend(fontsize=7);ax.grid(axis='y',alpha=.2)
    save(fig,out,'near-h-momentum-loads')
    fig,axes=plt.subplots(2,2,figsize=(10,6.5),sharex=True)
    for ax,field in zip(axes.flat,('velocity','k','outer_epsilon','reduced_pressure')):
        for policy,color in (('EP','#315c8c'),('Budget','#b16a28')):
            for tag,style in (('Nx2','-'),('Nx2half','--')):
                for kind,alpha in (('continuum_error',1.),('initial_drift',.45)):
                    q=[r for r in fields if r['policy']==policy and r['scope']=='full' and r['field']==field and r['kind']==kind and r['case'].startswith('selected/'+tag+'/')]
                    q.sort(key=lambda r:float(r['time']))
                    if q:ax.semilogy([float(r['time']) for r in q],[max(float(r['normalized_L2']),1e-18) for r in q],style,color=color,alpha=alpha,label=policy+' '+tag+' '+('error' if kind=='continuum_error' else 'drift'))
        ax.axhline(1e-4,color='#a33232',ls=':');ax.set_title(field);ax.set_ylabel('normalized physical L2');ax.grid(alpha=.2);ax.set_xlabel('t [s]');ax.legend(fontsize=6,ncol=2)
    fig.suptitle('Selected local mesh; frozen T3=.004096 s, T3 u*/H=.002048',fontsize=11)
    save(fig,out,'selected-error-drift')
    if (a.compact/'sensitivity.csv').exists():
        sensitivity=load('sensitivity.csv');fig,axes=plt.subplots(1,2,figsize=(10,3.7),sharey=True)
        for ax,group in zip(axes,('selected_horizontal','selected_timestep')):
            for policy,color in (('EP','#315c8c'),('Budget','#b16a28')):
                q=[r for r in sensitivity if r['policy']==policy and r['group']==group];times=sorted(set(float(r['cadence_index']) for r in q))
                ax.plot(np.array(times)*.001024,[max(float(r['fraction']) for r in q if float(r['cadence_index'])==t) for t in times],'o-',color=color,label=policy)
            ax.axhline(1,color='#a33232',ls='--');ax.set_yscale('log');ax.set_xlabel('common physical time [s]');ax.set_ylabel('(difference + uncertainty) / 10% budget');ax.set_title(group.replace('selected_','')+' sensitivity');ax.grid(alpha=.2);ax.legend()
        save(fig,out,'complete-sensitivity')
    manifest={f.name:__import__('hashlib').sha256(f.read_bytes()).hexdigest() for f in out.iterdir() if f.is_file()}
    (out/'figures.json').write_text(json.dumps(dict(schema='WJ5ScientificFiguresV1',files=manifest,scope='serial flat2D reference-Dirichlet finite strip; spatial residual and actual-stage defect remain separate tables'),indent=2)+'\n')

if __name__=='__main__':main()
