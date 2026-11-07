#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Common physical norms/bands for the declared height experiment; no coefficient comparisons."""
import argparse
import datetime
import math
from pathlib import Path
import json
from run_cls_epawf_picard import load,write_json,require,sha
from run_cls_epawf_height_comparison import CASES


def norm(weights,values):return math.sqrt(max(0.,math.fsum(w*x*x for w,x in zip(weights,values))))
def differences(a,b,reference=False):
    a,b=a['common_samples'],b['common_samples'];w=a['weights']
    require(len(w)==len(b['weights']) and max(abs(x-y) for x,y in zip(w,b['weights']))<1e-14,'physical comparison rules differ')
    fields=('u','k','epsilon') if reference else ('u','v','k','epsilon','pressure')
    out={}
    for key in fields:
        weights=a['epsilon_weights'] if key=='epsilon' else w
        other=b[('reference_'+key) if reference else key]
        delta=[x-y for x,y in zip(a[key],other)]
        if key=='pressure':
            shift=math.fsum(v*x for v,x in zip(weights,delta))/sum(weights);delta=[x-shift for x in delta];out['common_gauge_shift']=shift
        out[key]=norm(weights,delta)
    return out


def independent_reference(native):
    out=[];us=.4;nu=1.5e-5;h=.1
    for row in native['rows']:
        K=row['K0'];root=math.sqrt(K);x=h*root/(70*nu)
        n=max(.5-.4*x**.7,1-.79*x**(-.28));Zv=nu*10.7*(1-x**n)/root
        slope=.09*2.55*root;ell=nu/slope-Zv;H=row['H']
        j0=math.log1p(H/ell)/slope;UH=us*us*j0;z0=H*math.exp(-row['kappa']*UH/us)
        require(abs(j0-row['J0'])<1e-11 and abs(z0-row['z0_c'])<1e-14,'independent aerodynamic reference differs')
        out.append(dict(H=H,z_m=row['z_m'],J0=j0,U_H=UH,z0_c=z0,velocity_scale=ell,
                        fixed_z0_traction=(us/row['kappa']*math.log(H/native['z0_baseline']))/j0))
    return out


def compact_diagnostic(q):
    return {k:v for k,v in q.items() if k not in ('common_samples','weak_k_component_loads')}


def attempts_summary(attempts,through):
    attempts=[a for a in attempts if a.get('numerically_accepted',False) and a['time']<=through+1e-13]
    bands={};increments=[];solves={}
    for a in attempts:
        for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves'):
            solves[k]=solves.get(k,0)+a.get(k,0)
        require(all(a[k]==1 for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls')) and a['additional_governing_equation_solves']==0,'production schedule changed')
        b=a.get('k_bound',{})
        for r in (b.get('physical_bands') or {}).get('regions',[]):
            acc=bands.setdefault(r['name'],dict(lower=r['lower'],upper=r['upper'],signed=0.,positive=0.,negative=0.,absolute=0.))
            for k in ('signed','positive','negative','absolute'):acc[k]+=r['correction_'+k]
        increments.append({k:a[k] for k in a if k in ('time','raw_state_sha256','corrected_state_sha256','k_raw_step_mass_norm','k_corrected_step_mass_norm','k_raw_linear_residual','k_added_frozen_residual','k_bound_storage_residual','refreshed_trace_defect','refreshed_original_block_norms','seconds')})
    return dict(accepted_steps=len(attempts),solves=solves,activations=sum(a.get('bound_activations',0) for a in attempts),
        seconds=sum(a['seconds'] for a in attempts),correction_integral=sum(a.get('k_bound',{}).get('integral_total',0.) for a in attempts),
        correction_mass_norm_sum=sum(a.get('k_bound',{}).get('mass_norm_correction',0.) for a in attempts),
        maximum_correction_raw_increment_ratio=max((a.get('k_bound',{}).get('correction_over_raw_step',0.) for a in attempts),default=0.),
        fixed_band_correction_totals=bands,increments=increments)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--screening',type=Path);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();r=a.root.resolve();s=(a.screening or r/'screening').resolve();a.output.mkdir(parents=True,exist_ok=True)
    out=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),root=str(r),
        common_norms='dimensional physical L2; pressure difference has one common weighted-mean gauge removed; epsilon only z>=2',
        common_scales=dict(U=.4,K=.51050801305988613,H=2.,volume=24.,outer_volume=12.),
        reference=load(r/'reference-heights.json'),cases={},paired_differences={},cross_height_differences={})
    out['independent_reference']=independent_reference(out['reference'])
    states={}
    for case in CASES:
        path=s/case/'result.json'
        if not path.exists():continue
        result=load(path);initial=result.get('initialization')
        if not initial:
            out['cases'][case]=result;continue
        states[case]={0:initial};q=out['cases'][case]=dict(definition=result['definition'],status=result['status'],steps=result['steps'],time=result['time'],failure=result['failure'],
            actual_h0_sha256=result['actual_h0_sha256'],initial_checkpoint=result['native_h0_checkpoint'],accepted_checkpoint=result['accepted_checkpoint'],
            initial=compact_diagnostic(initial),initial_error_to_own_analytic= differences(initial,initial,True),observations={},qualification_events=result['qualification_events'],
            total_execution_seconds=result['trajectory_receipt']['seconds'])
        for step,packet in result['observations'].items():states[case][int(step)]=load(Path(packet['path']))
        for step,d in states[case].items():
            if step==0:continue
            row=dict(diagnostic=compact_diagnostic(d),drift_from_actual_h0=differences(d,initial),error_to_own_analytic=differences(d,d,True),
                     work=attempts_summary(result['attempts'],d['time']))
            # Source estimates are time integrals of sampled physical rates.
            # They are explicitly separate from exact discrete correction sums.
            estimates={};old=initial
            for n,current in sorted(states[case].items()):
                if n==0 or n>step:continue
                dt=current['time']-old['time']
                for x,y in zip(old['fixed_bands']['regions'],current['fixed_bands']['regions']):
                    budget=estimates.setdefault(x['name'],dict(production=0.,destruction=0.,valid=True))
                    for k in ('production','destruction'):
                        if x[k] is None or y[k] is None:budget['valid']=False
                        else:budget[k]+=.5*dt*(x[k]+y[k])
                old=current
            row['sampled_continuum_time_integrals']=estimates;q['observations'][str(step)]=row
        out['cases'][case]['sampling_checks']={k:v for k,v in result['sampling_checks'].items() if k!='samples'}
    for first,second in [('A','B'),('B','C'),('C','D'),('E','F'),('G','H')]:
        if first in states and second in states:
            out['paired_differences'][first+'-'+second]={str(step):differences(states[first][step],states[second][step]) for step in states[first].keys() & states[second].keys()}
    for first,second in [('C','E'),('E','G'),('D','F'),('F','H')]:
        if first in states and second in states:
            out['cross_height_differences'][first+'-'+second]=dict(
                reference_velocity_offset=norm(states[first][0]['common_samples']['weights'],[x-y for x,y in zip(states[first][0]['common_samples']['reference_u'],states[second][0]['common_samples']['reference_u'])]),
                represented={str(step):differences(states[first][step],states[second][step]) for step in states[first].keys() & states[second].keys()})
    write_json(a.output/'native.json',out)
    try:
        import matplotlib;matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        for pair in [('C','D'),('E','F'),('G','H')]:
            if not all(c in states for c in pair):continue
            fig,axes=plt.subplots(1,3,figsize=(12,4.5),layout='constrained')
            for c in pair:
                for step in (0,32):
                    if step not in states[c]:continue
                    q=states[c][step];p=min(q['profiles'],key=lambda p:abs(p['x']-1.6));rows=p['profile'];H=q['H']
                    for ax,index,scale in zip(axes,[6,2,4],[1,1,H]):
                        ax.plot([v[index] for v in rows if v[0]<=H],[v[0]/scale for v in rows if v[0]<=H],label=f'{c} t={q["time"]:.3g}')
                p=states[c][0]['profiles'][-1];rows=p['profile'];H=states[c][0]['H']
                axes[1].plot([v[3] for v in rows if v[0]<=H],[v[0] for v in rows if v[0]<=H],':',color='black')
                axes[2].plot([v[5] for v in rows if v[0]<=H],[v[0]/H for v in rows if v[0]<=H],':',color='black')
            axes[0].set(xlabel='transported k',ylabel='physical z',ylim=(0,.1));axes[1].set(xlabel='velocity u',ylabel='physical z')
            axes[2].set(xlabel='represented du/dz',ylabel='z/H',xscale='symlog');axes[0].legend(fontsize=8)
            fig.savefig(a.output/('profiles-'+''.join(pair)+'.png'),dpi=160);plt.close(fig)
        fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
        for case,by_step in states.items():
            times=[];integ=[];tau=[]
            for step,q in sorted(by_step.items()):
                times.append(q['time']);integ.append(next(v['k_integral'] for v in q['fixed_bands']['regions'] if v['name']=='z_le_0.01'));tau.append(q['mean_wall_stress'])
            axes[0].plot(times,integ,'o-',label=case);axes[1].plot(times,tau,'o-',label=case)
        axes[0].set(xlabel='physical time',ylabel='integral k, 0<z<=.01');axes[1].set(xlabel='physical time',ylabel='mean wall stress');axes[0].legend(ncol=2)
        fig.savefig(a.output/'fixed-band-wall-stress.png',dpi=160);plt.close(fig)
    except ImportError:pass
    print(a.output/'native.json')


if __name__=='__main__':main()
