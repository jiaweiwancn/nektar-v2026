#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Draw retained native Q3 polynomials, not lines between sparse profile points."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run_cls_epawf_picard import load,require,write_json
from run_cls_epawf_height_comparison import UPSTREAM_CASES

def basis(x,degree=3,derivative=False):
    if derivative:return np.array([-.5,.5,-x/2,.5-1.5*x*x][:degree+1])
    return np.array([(1-x)/2,(1+x)/2,(1-x*x)/4,x*(1-x*x)/2][:degree+1])

def value(report,state,block,x,z,derivative=False):
    mesh=report['fields'][block]['mesh']
    e=next(e for e in mesh if e['lower'][0]-1e-13<=x<=e['upper'][0]+1e-13 and e['lower'][2]-1e-13<=z<=e['upper'][2]+1e-13)
    n=int(e['degree'])+1;lx=e['upper'][0]-e['lower'][0];lz=e['upper'][2]-e['lower'][2]
    a=2*(x-e['lower'][0])/lx-1;b=2*(z-e['lower'][2])/lz-1
    coefficients=np.array([state[int(report['offsets'][block])+int(i)]*s for i,s in zip(e['ids'],e['sign'])]).reshape((n,n))
    return float(basis(b,n-1,derivative)@coefficients@basis(a,n-1))*(2/lz if derivative else 1)

def primitive(x):return np.array([x/2-x*x/4,x/2+x*x/4,x/4-x**3/12,x*x/4-x**4/8])

def integral(report,state,top,xbounds=(-4.,2.)):
    total=0.;b=2
    for e in report['fields'][b]['mesh']:
        x0,x1=max(e['lower'][0],xbounds[0]),min(e['upper'][0],xbounds[1]);z0,z1=e['lower'][2],min(e['upper'][2],top)
        if x1<=x0 or z1<=z0:continue
        lx=e['upper'][0]-e['lower'][0];lz=e['upper'][2]-e['lower'][2]
        px=primitive(2*(x1-e['lower'][0])/lx-1)-primitive(2*(x0-e['lower'][0])/lx-1)
        pz=primitive(2*(z1-e['lower'][2])/lz-1)-primitive(2*(z0-e['lower'][2])/lz-1)
        c=np.array([state[int(report['offsets'][b])+int(i)]*s for i,s in zip(e['ids'],e['sign'])]).reshape((4,4))
        total+=lx*lz/4*float(pz@c@px)
    return total

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--screening',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    data={c:load(a.screening/c/'trajectory/case.evolution.json') for c in UPSTREAM_CASES}
    initial={c:load(a.screening/c/'initial/case.evolution.json') for c in UPSTREAM_CASES}
    n=min(32,min(int(d['steps']) for d in data.values()));t=n*.000625;validation={}
    for c,d in data.items():
        q=d['height_comparison_final'];error=0.
        for profile in q['profiles']:
            for row in profile['profile']:
                error=max(error,abs(value(d,d['state'],2,profile['x'],row[0])-row[6]),
                    abs(value(d,d['state'],0,profile['x'],row[0])-row[2]))
        integral_error=0.
        for band in q['fixed_bands']['regions']:
            if band.get('x_region','all')!='all':continue
            integral_error=max(integral_error,abs(integral(d,d['state'],band['upper'])-band['k_integral']))
        require(error<1e-10 and integral_error<1e-10,'native coefficient rendering/integration disagrees with probe')
        validation[c]=dict(native_point_error=error,physical_integral_error=integral_error)
    xs=[-3.5,-2.625,-1.,1.6];pairs=[('U05-graded','U05-one'),('U20-graded','U20-one')]
    for normalized in [False,True]:
        fig,axes=plt.subplots(2,4,figsize=(14,7),sharex=True)
        for row,(g,o) in enumerate(pairs):
            H=initial[g]['height_comparison_final']['H'];K0=initial[g]['height_comparison_final']['reference']['K0']
            z=np.r_[0.,np.geomspace(1e-7,H,400)]
            for col,x in enumerate(xs):
                ax=axes[row,col]
                for c,color in [(g,'#1764ab'),(o,'#ce4b24')]:
                    d=data[c];s=d['accepted_states'][n]
                    ax.plot([value(d,s,2,x,r)/K0 for r in z],z/H if normalized else z,color=color,label=c)
                    ax.plot([value(d,initial[c]['state'],2,x,r)/K0 for r in z],z/H if normalized else z,color=color,ls=':',alpha=.65)
                ax.axvline(1,color='gray',lw=.7);ax.set_yscale('symlog',linthresh=1e-4 if normalized else 1e-5)
                ax.set_ylim(0,1 if normalized else H)
                ax.grid(alpha=.2);ax.set_title(f'H={H:g}, x={x:g}');ax.set_xlabel('k / K0')
                if col==0:ax.set_ylabel('z / H' if normalized else 'z')
                if col==0:ax.legend(fontsize=8)
        fig.suptitle(f'Actual native k at T={t:g}; dotted: each case native h0')
        fig.tight_layout();fig.savefig(a.output/('k-normalized.png' if normalized else 'k-physical.png'),dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,4,figsize=(14,7))
    for row,(g,o) in enumerate(pairs):
        q=initial[g]['height_comparison_final'];H=q['H'];ref=q['reference'];ell=ref['velocity_scale'];z=np.r_[0.,np.geomspace(1e-7,H,400)]
        gradient=ref['U_H']/np.log1p(H/ell)/(z+ell)
        for col,x in enumerate(xs):
            ax=axes[row,col];ax.plot(z,gradient,color='black',lw=1,label='analytical AWF')
            for c,color in [(g,'#1764ab'),(o,'#ce4b24')]:
                d=data[c];ax.plot(z,[value(d,d['accepted_states'][n],0,x,r,True) for r in z],color=color,label=c)
                ax.plot(z,[value(d,initial[c]['state'],0,x,r,True) for r in z],color=color,ls=':',alpha=.65)
            ax.set_xscale('symlog',linthresh=1e-5);ax.set_yscale('symlog',linthresh=1)
            ax.set_xlim(0,H)
            ax.set_title(f'H={H:g}, x={x:g}');ax.set_xlabel('z');ax.set_ylabel('du / dz');ax.grid(alpha=.2)
            if col==0:ax.legend(fontsize=8)
    fig.suptitle(f'Native gradients at T={t:g}; dotted: native h0');fig.tight_layout();fig.savefig(a.output/'gradients.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(11,7))
    for c,d in data.items():
        times=[j*d['dt'] for j in range(n+1)];states=d['accepted_states'][:n+1]
        for ax,top,title in [(axes[0,0],.01,'Integrated k, 0<z<=.01'),(axes[0,1],.1,'Integrated k, 0<z<=.1')]:
            ax.plot(times,[integral(d,s,top) for s in states],label=c);ax.set_title(title)
        attempts=[v for v in d['pimple_attempts'] if v['time']<=t+1e-13 and v['scheduled_complete'] and not v['failure']]
        ts=[v['time'] for v in attempts];budgets=[v['k_bound'] for v in attempts]
        axes[1,0].plot(ts,np.cumsum([b['integral_total'] for b in budgets]),label=c)
        axes[1,1].plot(ts,[b['correction_to_step_ratio'] for b in budgets],label=c)
    axes[1,0].set_title('Cumulative signed bound supply (integrated k)');axes[1,1].set_title('Correction / raw-step physical mass norm')
    for ax in axes.flat:ax.set_xlabel('Physical time');ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(a.output/'bands-and-correction.png',dpi=150);plt.close(fig)
    write_json(a.output/'plot-validation.json',dict(common_step=n,common_time=t,native_rendering=validation,
        note='exact tensor Modified_A Q3 coefficient evaluation and polynomial band integration; no clipping or state repair'))

if __name__=='__main__':main()
