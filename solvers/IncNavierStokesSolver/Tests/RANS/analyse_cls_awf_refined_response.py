#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only finite-horizon response about the actual refined P1 candidate."""
import argparse
from datetime import datetime,timezone
import math,os
from pathlib import Path
from analyse_cls_accepted_awf import normalized
from analyse_cls_epawf_height_comparison import differences
from cls_awf_native_metrics import flat_partition,gl4_differences,difference_localization
from continue_cls_epawf_pimple import atomic
from diagnose_cls_accepted_covariance import checkpoint,read_report
from qualify_cls_awf_neutral import observe
from run_cls_epawf_picard import THREADS,load,require,sha


def stepwise_k(study,mass_readback,initial_norm,through):
    """All accepted states, using M_ff (never the frozen transport A_k)."""
    import numpy as np
    from scipy.sparse import csr_matrix
    raw=load(mass_readback);require(raw['governing_solves']==raw['physical_commits']==0,'invalid native mass readback')
    b=raw['diagnostics']['neutral_reference_balance']['mass_metrics'][2]
    free=np.array(b['free_local_rows'],dtype=int);M=csr_matrix(np.array(b['consistent_free_mass']));offset=int(b['offset'])
    size=len(checkpoint(study/'parent')['coupled']);n=raw['diagnostics']['physical']['native_degrees_of_freedom'][2]
    essential=np.array(sorted(set(range(int(n)))-set(map(int,free))),dtype=int)
    progress=[load(study/name/'progress.json') for name in ('control','perturbed')]
    segments=[{int(v['end_step']):v for v in q['segments'] if v['end_step']<=through} for q in progress];rows=[]
    require(set(segments[0])==set(segments[1]),'accepted response blocks differ')
    for end in sorted(segments[0]):
        ss=[segments[j][end] for j in (0,1)];reports=[read_report(Path(v['receipt']['report']+'.gz')) for v in ss]
        first=int(ss[0]['start_step']);require(first==int(ss[1]['start_step']),'response block starts differ')
        states=[np.array(q['accepted_states']) for q in reports]
        require(states[0].shape==states[1].shape==(end-first+1,size),'accepted state inventory differs')
        delta=states[1][:,offset:offset+int(n)]-states[0][:,offset:offset+int(n)]
        require(np.max(abs(delta[:,essential]),initial=0)==0,'perturbed/control essential k data differ')
        norms=np.sqrt(np.maximum(0,np.sum(delta[:,free]*(M@delta[:,free].T).T,axis=1)))
        for j in range(1,end-first+1):rows.append(dict(step=first+j,elapsed_time=(first+j)*progress[0]['dt'],physical_k_L2=float(norms[j]),k_amplification=float(norms[j]/initial_norm)))
    require([v['step'] for v in rows]==list(range(1,through+1)),'stepwise response has a gap')
    return dict(mass_readback=str(mass_readback),mass_readback_sha256=sha(mass_readback),metric='physical consistent free k mass; exact zero essential differences verified',
        maximum_accepted_step_amplification=max(v['k_amplification'] for v in rows),rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('study','probe','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--mass-readback',type=Path,help='strict native balanced readback containing the consistent k mass')
    p.add_argument('--times',nargs='+',type=float,default=[.05,.1,.25,.5]);a=p.parse_args();os.environ.update(THREADS)
    a.output.mkdir(parents=True,exist_ok=False);study=load(a.study/'summary.json');require(study['validation']['status']=='passed','unvalidated seed contract')
    policy=load(a.study/'declared-policy.json');out=dict(status='running',study_summary_sha256=sha(a.study/'summary.json'),policy=policy,
        probe_sha256=sha(a.probe),started_UTC=datetime.now(timezone.utc).isoformat(),observations={},response={},governing_solves=0,physical_commits=0)
    data={};progress={name:load(a.study/name/'progress.json') for name in ('control','perturbed')}
    for name in ('control','perturbed'):
        out['observations'][name]={};data[name]={};q=progress[name]
        selections={'0.0':dict(checkpoint=q['seed_checkpoint'],xml=str(a.study/name/'initial/case.xml'),report=str(a.study/name/'initial/case.evolution.json'))}
        for t in a.times:
            require(str(t) in q['observations'],'requested response endpoint is not yet accepted: '+str(t))
            selections[str(t)]=q['observations'][str(t)]
        for t,row in selections.items():
            cp=Path(row['checkpoint']);saved=checkpoint(cp);selection=dict(checkpoint=str(cp),checkpoint_sha256=sha(cp/'State.bin'),xml=row['xml'],
                native_step=saved['steps'],native_time=saved['time'],source_report=row['report'],H=.5,balance=False)
            summary,physical=observe(a.probe,selection,a.output/name/t);out['observations'][name][t]=summary;data[name][t]=physical
            atomic(a.output/'summary.json',out)
    seed=differences(data['perturbed']['0.0'],data['control']['0.0'])
    expected=progress['perturbed']['seed_identity']['inner_k_seed']['physical_L2'];require(abs(seed['k']-expected)<2e-14,'actual seed mass norm disagrees')
    require(all(seed[k]==0 for k in ('u','v','epsilon','pressure')),'unintended initial induced field')
    out['initial_difference']=seed;out['initial_scaled_difference']=normalized(seed)
    for key in ['0.0']+list(map(str,a.times)):
        c,q=data['control'][key],data['perturbed'][key];delta=differences(q,c);norm=normalized(delta);x,_=flat_partition(Path(out['observations']['control'][key]['selection']['xml']))
        check=gl4_differences(q,c,x);require(max(abs(delta[k]-check['dimensional'][k]) for k in norm)<2e-12,'response GL4/16 norms disagree')
        out['response'][key]=dict(physical_time=c['time'],dimensional_difference=delta,normalized_difference=norm,
            k_amplification=delta['k']/seed['k'],combined_scaled_field_norm=math.sqrt(sum(v*v for v in norm.values())),
            perturbation_localization=difference_localization(q,c,x),independent_gl4=check,
            mean_traction_control=c['mean_wall_stress'],mean_traction_perturbed=q['mean_wall_stress'],
            traction_difference=q['mean_wall_stress']-c['mean_wall_stress'],
            control_drift_from_parent=normalized(differences(c,data['control']['0.0'])))
    values=[out['response'][str(t)]['k_amplification'] for t in a.times];last=values[-1]
    if a.mass_readback:
        out['stepwise_k_response']=stepwise_k(a.study,a.mass_readback,seed['k'],round(max(a.times)/progress['control']['dt']))
        for t in a.times:
            row=out['stepwise_k_response']['rows'][round(t/progress['control']['dt'])-1]
            require(abs(row['physical_k_L2']-out['response'][str(t)]['dimensional_difference']['k'])<2e-14,'native mass/common-panel response norm differs')
    if max(a.times)<.5:classification='interim_response_bounded_horizon_not_completed'
    elif len(values)>=2 and last<.9 and last<values[-2]:classification='decay_of_tested_k_disturbance_over_observed_horizon'
    elif len(values)>=2 and last>1.1 and last>values[-2]:classification='sustained_growth_over_observed_horizon'
    else:classification='bounded_or_unresolved_at_observed_horizon'
    peak=max(values+[out.get('stepwise_k_response',{}).get('maximum_accepted_step_amplification',0.)])
    out.update(status='passed',classification=classification,transient_amplification=peak>1.05,
        elapsed_horizon=max(a.times),conditional_extension_needed=None if max(a.times)<.5 else classification=='bounded_or_unresolved_at_observed_horizon',
        scope='one finite-amplitude shape about a selected discrete candidate; no general stability, faithful-Original evolution or timestep-independent-root claim',
        completed_UTC=datetime.now(timezone.utc).isoformat());atomic(a.output/'summary.json',out)


if __name__=='__main__':main()
