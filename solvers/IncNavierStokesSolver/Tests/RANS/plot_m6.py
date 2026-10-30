#!/usr/bin/env python3
"""Reduce native M6 CSVs and generate standalone scientific study figures."""
import argparse
import csv
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def rows(path):
    return list(csv.DictReader(path.open())) if path.exists() else []


def parameters(path):
    result={}
    for node in ET.parse(path).getroot().findall('CONDITIONS/PARAMETERS/P'):
        key,value=node.text.split('=',1)
        try:result[key.strip()]=float(value)
        except ValueError:pass
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--study',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True);study=args.study;reports={}
    plt.rcParams.update({'font.size':10,'axes.grid':True,'grid.alpha':.25,'savefig.dpi':160})
    def save(fig,name):
        fig.tight_layout();fig.savefig(args.output/(name+'.png'));fig.savefig(args.output/(name+'.pdf'));plt.close(fig)
    for directory in sorted(study.iterdir()):
        if not directory.is_dir():continue
        name=directory.name;session=directory/(name+'.xml')
        if not session.exists():continue
        p=parameters(session);errors=rows(directory/(name+'.rans-errors.csv'))
        if not errors:continue
        last=max(float(r['time']) for r in errors);fields={}
        for field in sorted({r['field'] for r in errors}):
            group=[r for r in errors if r['field']==field];end=next(r for r in reversed(group) if float(r['time'])==last)
            fields[field]={key:float(end[key]) for key in end if key not in ('time','field')}
            fields[field]['peak_l2']=max(float(r['normalized_l2']) for r in group)
            fields[field]['peak_linf']=max(float(r['pointwise_relative_linf']) for r in group)
        stages=rows(directory/(name+'.rans.csv'));walls=rows(directory/(name+'.awf.csv'))
        required=['u','k','eps' if 'eps' in fields else 'eps_tilde','rans_nut','rans_R_xx','rans_R_yy','rans_R_zz','rans_R_xz']
        record=dict(time=last,flow_times=last*p['Uref']/p['L'],parameters=p,fields=fields,
            target_met=all(fields[f]['peak_l2']<.001 and fields[f]['peak_linf']<.001 for f in required))
        if stages:
            record['stages']=dict(count=len(stages),last={k:float(v) for k,v in stages[-1].items()},
                max_weak_divergence=max(float(r['weak_divergence_l2']) for r in stages),
                max_live_systems=max(float(r['live_global_systems']) for r in stages),
                guards=sum(float(r['guard_activations']) for r in stages))
        if walls:
            record['wall']=dict(last={k:float(v) for k,v in walls[-1].items()},
                chi_min=min(float(r['chi_min']) for r in walls),chi_max=max(float(r['chi_max']) for r in walls),
                traction_relative_error=max(abs(float(r['traction_x'])/float(r['area'])/p['Ustar']**2-1) for r in walls),
                max_k_budget=max(abs(float(r['k_budget'])) for r in walls),
                max_momentum_budget=max(abs(float(r['momentum_budget_x'])) for r in walls),
                guards=sum(float(r['guard_activations']) for r in walls),
                fallbacks=sum(float(r[k]) for r in walls for k in ('cutoff_fallback_points','regime_fallback_points','range_fallback_points')))
        record['completed_five_flows']=record['flow_times']>=5
        reports[name]=record
    fig,axes=plt.subplots(2,4,figsize=(14,6),sharex=True)
    styles=[('standard','Standard','#7653a3','-'),('cls_original','CLS original','#c44238','-'),
            ('cls_fixed','CLS fixed','#e1912a','--'),('cls_switched','CLS switched','#15826d',':')]
    for dim in (2,3):
        for model,label,color,line in styles:
            name=f'five_d{dim}_{model}';directory=study/name
            error=rows(directory/(name+'.rans-errors.csv'));p=parameters(directory/(name+'.xml'))
            for j,field in enumerate(('u','k','eps' if model=='standard' else 'eps_tilde','rans_nut')):
                group=[r for r in error if r['field']==field]
                axes[dim-2,j].semilogy([float(r['time'])*p['Uref']/p['L'] for r in group],
                    [max(1e-12,100*float(r['normalized_l2'])) for r in group],label=label,color=color,ls=line,lw=1.5)
                if group and float(group[-1]['time'])*p['Uref']/p['L']<5:
                    axes[dim-2,j].scatter([float(group[-1]['time'])*p['Uref']/p['L']],
                        [max(1e-12,100*float(group[-1]['normalized_l2']))],marker='x',color=color,s=30)
        for j,title in enumerate(('Velocity','Turbulent energy','Dissipation','Eddy viscosity')):
            ax=axes[dim-2,j];ax.axhline(.1,color='black',lw=.8,ls='--');ax.set_title(f'{dim}D {title}');ax.set_ylim(1e-7,100)
            if j==0:ax.set_ylabel('Relative L2 error (%)')
            if dim==3:ax.set_xlabel('Flow-through times')
    axes[0,0].legend(fontsize=8);save(fig,'preservation')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    hnames=['study_h6','five_d2_cls_switched','study_h24'];pnames=['study_p3','five_d2_cls_switched','study_p5']
    rates={}
    for field,label in (('eps_tilde','Dissipation'),('rans_nut','Eddy viscosity'),('u','Velocity')):
        eh=[reports[n]['fields'][field]['normalized_l2'] for n in hnames]
        ep=[reports[n]['fields'][field]['normalized_l2'] for n in pnames]
        axes[0].loglog([6,12,24],eh,'o-',label=label);axes[1].semilogy([3,4,5],ep,'o-',label=label)
        rates[field]=[math.log(a/b,2) for a,b in zip(eh,eh[1:])]
    axes[0].set_xlabel('Vertical geometric layers (degree 4)');axes[1].set_xlabel('Polynomial degree (12 layers)')
    axes[0].set_ylabel('Final normalized L2 error');axes[1].set_ylabel('Final normalized L2 error')
    axes[0].set_title('Fixed H=2 m, P=1 m, N=4 m');axes[1].set_title('Fixed physical matching locations');axes[1].legend()
    save(fig,'refinement')
    fig,axes=plt.subplots(1,3,figsize=(11,4.6),sharey=True)
    for name,label in (('five_d2_cls_switched','2D'),('five_d3_cls_switched','3D')):
        values=rows(study/name/'profiles.csv')
        for station,color in zip((75,150,225,285),plt.cm.viridis(np.linspace(.1,.9,4))):
            for ax,field in zip(axes,('u','eps_tilde','rans_nut')):
                group=[r for r in values if r['field']==field and float(r['x'])==station]
                z=[float(r['height']) for r in group]
                error=[100*(float(r['value'])-float(r['reference']))/float(r['reference']) for r in group]
                ax.plot(error,z,ls='-' if label=='2D' else '--',color=color,label=f'{label} x={station:g} m')
    for ax,label in zip(axes,('Velocity','Dissipation','Eddy viscosity')):
        ax.set_xlabel('Profile error (%)');ax.set_title(label);ax.axvline(0,color='black',lw=.7)
    axes[0].set_ylabel('Physical height (m)');axes[2].legend(fontsize=7)
    save(fig,'profiles')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for name,label in (('five_d2_cls_original','Original'),('five_d2_cls_switched','Switched'),('study_perturbed','Perturbed')):
        data=rows(study/name/(name+'.awf.csv'));p=parameters(study/name/(name+'.xml'))
        t=[float(r['time'])*p['Uref']/p['L'] for r in data]
        axes[0].plot(t,[float(r['mean_production'])-float(r['mean_dissipation']) for r in data],label=label)
        axes[1].plot(t,[float(r['chi_min']) for r in data],label=label)
    axes[0].set_ylabel('Mean strip P - epsilon (m²/s³)');axes[1].set_ylabel('Minimum effective correction weight')
    for ax in axes:ax.set_xlabel('Flow-through times');ax.legend()
    save(fig,'wall-response')
    result=dict(cases=reports,h_rates=rates)
    (args.output/'metrics.json').write_text(json.dumps(result,indent=2)+'\n')
    with (args.output/'metrics.csv').open('w') as file:
        writer=csv.writer(file,lineterminator='\n');writer.writerow(['case','flows','field','final_L2','peak_L2','peak_Linf','first_layer_L2','upper_layer_L2'])
        for name,record in reports.items():
            for field,value in record['fields'].items():writer.writerow([name,record['flow_times'],field,value['normalized_l2'],value['peak_l2'],value['peak_linf'],value['first_layer_relative_l2'],value['upper_layer_relative_l2']])
    print(json.dumps(dict(cases=len(reports),h_rates=rates),indent=2))

if __name__=='__main__':main()
