#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Apply the existing fixed-window policy to immutable selected-operator readbacks."""
import argparse
import json
import math
import os
from pathlib import Path
from run_cls_epawf_picard import load,sha,THREADS,require
from continue_cls_epawf_pimple import atomic


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--operators',type=Path,required=True)
    p.add_argument('--policy',type=Path,required=True);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False)
    source=load(a.operators/'summary.json');policy=load(a.policy);cases={int(n.split('-')[-1]):c for n,c in source['cases'].items() if n.startswith('Regularized-')}
    full={n:load(Path(c['full_free_loads_retained']))['diagnostics']['physical'] for n,c in cases.items()}
    h0=full[0];samples=h0['common_samples'];bulk=sum(w*u for w,u in zip(samples['weights'],samples['u']))/sum(samples['weights'])
    flow=6/bulk;wall=2/.4;out=dict(policy=policy,policy_sha256=sha(a.policy),source_summary_sha256=sha(a.operators/'summary.json'),
        wall_time_scale=wall,flow_through_time=flow,bulk_initial_velocity=bulk,checkpoints={},windows={},classification='continuing_transient')
    records={}
    for line in (a.campaign/'steps.jsonl').read_text().splitlines():
        try:r=json.loads(line)
        except json.JSONDecodeError:break
        records[int(r['global_step'])]=r
    for n,c in sorted(cases.items()):
        q=full[n];audit=c['complete_audit'];face=audit['interface'];area=face[0]
        out['checkpoints'][str(n)]=dict(time=c['time'],elapsed_wall_scales=c['time']/wall,elapsed_flow_throughs=c['time']/flow,
            common_normalized_spatial_first4=q['common_normalized_spatial_first4'],native_scaled_spatial_blocks=q['scaled_native_spatial_blocks'],
            covariance=audit['accepted_awf_covariance'],strong_continuity=audit['strong_continuity'],interface=face,
            mean_pressure_jump=face[1]/area,mean_normal_traction_jump=face[7]/area,
            recovered_mean_nonpressure_normal_stress_jump=(face[1]+face[7])/area,
            mean_traction=q['mean_wall_stress'],traction_relative_reference_error=q['mean_wall_stress']/.16-1,
            matching_defect=q['matching_trace_defect'],constitutive_modification=audit.get('constitutive_modification'),
            scalar_domain=audit['accepted_awf_scalar_domain'],
            continuous_ground=records[n]['attempt']['k_bound']['corrected_ground'] if n in records else None,
            profile_summaries=c['physical']['profile_summaries'],physical_integrated_source_rates=audit['accepted_awf_sources'],
            inner_constitutive_modification=audit.get('inner_constitutive_modification'),inner_stress=audit.get('inner_stress'),
            fixed_bands=q['fixed_bands'],drift_from_h0=c.get('drift_from_actual_native_h0'),
            counterfactual_Original_minimum=c['Original_covariance']['same_state_Original']['eigenvalues'][0])
    for name,window in source['windows'].items():
        first,last=map(int,name.split('-'));duration=window['end']-window['start'];end=cases[last];field=window['common_normalized_field_change']
        q=full[last];rows=[records[i] for i in range(first+1,last+1) if i in records]
        require(len(rows)==last-first,'native correction window is incomplete')
        b=[r['attempt']['k_bound'] for r in rows];signed=sum(x['integral_total'] for x in b)
        integrals=[x['integrals_inner_outer_near_positive_negative'] for x in b]
        positive=sum(x[3] for x in integrals);negative=sum(x[4] for x in integrals)
        fields_ok=all(field[k]<=limit for k,limit in policy['normalized_window_L2_change_maximum'].items())
        residual_ok=max(q['scaled_native_spatial_blocks'])<=policy['normalized_selected_spatial_residual_maximum_each_block']
        supply_rate=abs(signed)/duration/(policy['common_normalization']['K']*24/wall)
        out['windows'][name]=dict(window,primary_window=abs(duration-policy['primary_physical_window'])<1e-12,
            fields_pass=fields_ok,residuals_pass=residual_ok,traction_pass=abs(window['traction_change'])/.16<=policy['relative_wall_stress_window_change_maximum'],
            matching_pass=q['matching_trace_defect']/policy['common_normalization']['epsilon']<=policy['relative_matching_defect_maximum'],
            correction_signed=signed,correction_positive=positive,correction_negative=negative,correction_absolute=positive-negative,
            correction_rate_signed=signed/duration,normalized_bound_supply_rate=supply_rate,bound_activations=sum(x['activations'] for x in b))
    primary=[w for w in out['windows'].values() if w['primary_window']]
    passed=lambda w:all(w[k] for k in ('fields_pass','residuals_pass','traction_pass','matching_pass')) and w['normalized_bound_supply_rate']<=policy['bound_supply_per_time_relative_to_K0_volume_per_wall_time_maximum']
    if len(primary)>=3 and all(passed(w) for w in primary[-3:]):out['classification']='candidate_steady_selected_operator_at_declared_screen'
    out['caveat']='Not a temporal-error or spatial-accuracy certificate. Instantaneous source rates and time-integrated discrete correction amounts have distinct units. No statistically stationary claim follows from endpoint samples.'
    atomic(a.output/'summary.json',out)
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    names=('u','v','k','epsilon','pressure')
    for key in names:
        ws=[w for w in out['windows'].values() if w['primary_window']]
        axes[0,0].semilogy([w['end'] for w in ws],[w['common_normalized_field_change'][key] for w in ws],'o-',label=key)
    times=[v['time'] for v in out['checkpoints'].values()];rows=list(out['checkpoints'].values())
    for i,key in enumerate(('U','W','k','epsilon')):axes[0,1].semilogy(times,[r['common_normalized_spatial_first4'][i] for r in rows],'o-',label=key)
    axes[0,0].axhline(1e-4,color='grey',ls=':');axes[0,0].axhline(1e-3,color='grey',ls='--');axes[0,1].axhline(1e-3,color='grey',ls=':')
    axes[0,0].set(xlabel='window endpoint T (window .5)',ylabel='normalized physical L2 change');axes[0,1].set(xlabel='T',ylabel='normalized alpha=0 selected residual')
    axes[0,0].legend();axes[0,1].legend()
    axes[1,0].plot(times,[r['mean_traction'] for r in rows],'o-');axes[1,0].axhline(.16,color='grey',ls=':');axes[1,0].set(xlabel='T',ylabel='mean wall traction')
    axes[1,1].plot(times,[r['strong_continuity'][0] for r in rows],'o-',label='inner divergence L2')
    axes[1,1].plot(times,[r['strong_continuity'][1] for r in rows],'o-',label='outer divergence L2');axes[1,1].set(xlabel='T',ylabel='dimensional divergence norm');axes[1,1].legend()
    fig.savefig(a.output/'long-window-trends.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(3,3,figsize=(12,10),layout='constrained')
    for column,initial in enumerate(h0['profiles']):
        for n,q in sorted(full.items()):
            profile=next(v for v in q['profiles'] if v['x']==initial['x']);rows=profile['profile'];label=f"T={cases[n]['time']:.2f}"
            axes[0,column].plot([v[2] for v in rows],[v[0] for v in rows],'.-',label=label)
            axes[1,column].plot([v[6] for v in rows],[v[0] for v in rows],'.-',label=label)
            band=[v for v in rows if v[0]<=.1]
            axes[2,column].plot([v[6] for v in band],[v[0] for v in band],'.-',label=label)
        ref=initial['profile'];axes[0,column].plot([v[3] for v in ref],[v[0] for v in ref],'k--',label='analytical reference')
        axes[1,column].plot([v[7] for v in ref],[v[0] for v in ref],'k--')
        band=[v for v in ref if v[0]<=.1];axes[2,column].plot([v[7] for v in band],[v[0] for v in band],'k--')
        for row,field in enumerate(('U','k','k (z<=.1)')):axes[row,column].set(xlabel=field,ylabel='physical z')
        axes[0,column].set_title(f"x={initial['x']}")
    axes[0,0].legend(fontsize=8);fig.savefig(a.output/'long-run-profiles.png',dpi=160);plt.close(fig)
    print(a.output/'summary.json')

if __name__=='__main__':main()
