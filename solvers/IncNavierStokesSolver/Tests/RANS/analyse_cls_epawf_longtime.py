#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Actual-state contact/h scans and physical-window readback of the campaign.

All k-only solves below are diagnostic replays, never production callbacks.
Exploratory settling criteria are recorded before computing statistics.
"""
import argparse
import json
import math
from pathlib import Path
import shutil
import struct

from run_cls_epawf_picard import load,run,write_json,require,state_hash,sha
from run_cls_epawf_pimple import identity
from compare_cls_epawf_pimple import differences,same_history
from continue_cls_epawf_pimple import header


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--prior',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.campaign=a.campaign.resolve();a.prior=a.prior.resolve();a.probe=a.probe.resolve()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    criteria=dict(primary_window=.5,sensitivity_window=1.,minimum_adjacent_windows=3,
        minimum_elapsed_flow_through_times=1.,normalized_field_change=1e-3,
        normalized_spatial_rate=1e-3,normalized_correction_rate=1e-3,
        normalized_matching_defect=1e-3,relative_mean_and_RMS_change=.01,
        normalization='Ustar, K0, epsilon(H), Ustar^2; physical volume/area weighted norms; H/Ustar time')
    write_json(a.output/'criteria.json',criteria)
    prior=load(a.prior/'fixed-horizon/bound-0/case.evolution.json')
    progress=load(a.campaign/'main/progress.json');last=load(Path(progress['segment'])/'case.evolution.json')
    def packet(report,step):
        row=next(v for v in report['pimple_attempts'] if round(v['time']/report['dt'])==step)
        return load(Path(row['k_snapshot']))
    p4,p8=packet(prior,4),packet(prior,8)
    states=dict(h0=(prior['accepted_states'][0],0.,None),reference=(prior['fixed_reference'],0.,None),
        pre_k_first=(p4['pre_k'],p4['time'],p4),contact0025=(p4['corrected'],p4['time'],p4),
        contact005=(p8['corrected'],p8['time'],p8),later=(last['state'],last['time'],packet(last,int(last['steps']))))
    summary=dict(schema=1,criteria=criteria,progress=progress,probe=identity(a.probe),states={},checks={})
    def native(name,state,time,source=None,operation='contact',h=(),isolate=False,export_mass=False,xml=None,alpha=0.):
        folder=a.output/name;folder.mkdir();shutil.copy2(xml or a.prior/'fixed-horizon/bound-0/case.xml',folder/'case.xml')
        witnesses=[v['coordinate'] for v in (source or {}).get('k_bound',{}).get('active_constraints',[])]
        request=dict(operation=operation,state=state,history=(source or {}).get('history',state),time=time,alpha=alpha,
            h=list(h),witnesses=witnesses,inner_stress_model='SaturatedR2',require_history_h0=False,
            isolate=isolate,export_mass=export_mass)
        write_json(folder/'input.json',request)
        receipt=run([str(a.probe),'case.xml','input.json','result.json'],folder,'probe.log',1800)
        write_json(folder/'receipt.json',receipt);require(receipt['exit_code']==0,'diagnostic replay failed: '+name)
        return load(folder/'result.json')
    hs=(.000625,.00015625,1e-5,1e-6,1e-7,1e-8,1e-9)
    for name,(state,t,source) in states.items():
        result=native(name,state,t,source,h=hs if source else (),isolate=name=='contact005',export_mass=name=='contact0025')
        require(result['retained_state_unchanged'] and result['state_sha256']==state_hash(state),'selected field changed')
        compact={k:v for k,v in result.items() if k not in ('free_k_mass','k_integral_load','near_k_integral_load','free_k_rows','offsets','scan')}
        compact['scan']=[]
        for v in result['scan']:
            require(v['raw_linear_residual']<1e-9,'diagnostic linear solve failed')
            compact['scan'].append({k:x for k,x in v.items() if k not in ('storage_reaction','added_frozen_residual','bound')})
            compact['scan'][-1].update(storage_reaction_max=max(map(abs,v['storage_reaction'])),
                added_frozen_residual_max=max(map(abs,v['added_frozen_residual'])),
                bound={k:x for k,x in v['bound'].items() if k not in ('active_constraints',)})
        summary['states'][name]=compact;write_json(a.output/'summary.partial.json',summary)
    # This audit replay identifies the actual first uncaught constitutive k
    # witness, rather than inferring it from the already positive ground.
    failure=packet(last,int(last['steps'])+1)
    summary['rejection_witness']=native('rejection-witness',last['rejected_state'],failure['time'],failure,
        operation='audit_witness',alpha=failure['alpha'])
    cp=header(Path(progress['checkpoint']),True)
    require(state_hash(cp['coupled_state'])==state_hash(last['state']),'last accepted checkpoint state differs')
    require([state_hash(v) for v in cp['history']]==[state_hash(v) for v in last['native_history']],'last accepted checkpoint history differs')
    require(cp['history_times']==last['native_history_times'],'checkpoint history times differ')
    summary['checks']['last_accepted_checkpoint_matches_all_native_history']=True
    # Physical fixed windows from actual retained coefficients, no interpolation.
    sampled={}
    for name in ('h0','reference','contact0025','contact005','later'):
        state,t,source=states[name];r=native('samples-'+name,state,t,source,operation='samples')
        r['initial']=prior['accepted_states'][0];sampled[name]=r
    summary['fixed_window_differences']=[]
    for x,y in (('h0','contact0025'),('contact0025','contact005'),('contact005','later')):
        first=dict(sampled[x],time=sampled[y]['time'])
        summary['fixed_window_differences'].append(dict(start=sampled[x]['time'],end=sampled[y]['time'],
            duration=sampled[y]['time']-sampled[x]['time'],physical_L2=differences(first,sampled[y])))
    reference=dict(sampled['reference'],time=sampled['later']['time'])
    initial=dict(sampled['h0'],time=sampled['later']['time'])
    summary['later_minus_neutral_reference']=differences(reference,sampled['later'])
    summary['later_minus_actual_h0']=differences(initial,sampled['later'])
    mass_packet=load(a.output/'contact0025/result.json');mass=mass_packet['free_k_mass'];rows=list(map(int,mass_packet['free_k_rows']))
    offset=int(mass_packet['offsets'][2]);integral=mass_packet['k_integral_load'];near=mass_packet['near_k_integral_load']
    dot=lambda x,y:math.fsum(a*b for a,b in zip(x,y))
    def mv(m,x):return [dot(r,x) for r in m]
    def norm(x):return math.sqrt(max(0.,dot(x,mv(mass,x))))
    attempts={round(v['time']/prior['dt']):v for v in prior['pimple_attempts'] if v['time']<=.0025+1e-13}
    for line in (a.campaign/'main/steps.jsonl').read_text().splitlines():
        r=json.loads(line);attempts[r['global_step']]=r['attempt']
    attempts[round(failure['time']/prior['dt'])]=last['pimple_attempts'][-1]
    budgets=[];cumulative=0.
    for step,row in sorted(attempts.items()):
        q=load(Path(row['k_snapshot']));raw=[q['raw'][i]-q['pre_k'][i] for i in rows]
        correction=[q['corrected'][i]-q['raw'][i] for i in rows];net=[q['corrected'][i]-q['pre_k'][i] for i in rows]
        defect=max(abs(x+y-z) for x,y,z in zip(raw,correction,net));require(defect<1e-15,'increment cancellation identity differs')
        accepted=step<=last['steps'];b=row['k_bound'];cumulative+=b['integral_total'] if accepted else 0.
        payload=Path(q['matrix_file']).read_bytes();n=len(rows);matrix=struct.unpack('<'+str(n*n)+'d',payload)
        added=[dot(matrix[i*n:(i+1)*n],correction) for i in range(n)];storage=[v/prior['dt'] for v in mv(mass,correction)]
        with (a.output/f'budget-vectors-{step:05d}.json').open('w') as stream:
            json.dump(dict(step=step,A_delta=added,M_delta_per_dt=storage),stream)
        full_k=q['corrected'][offset:offset+len(integral)];delta_full=[q['corrected'][offset+i]-q['raw'][offset+i] for i in range(len(integral))]
        budgets.append(dict(step=step,time=q['time'],accepted=accepted,raw_norm=norm(raw),corrected_norm=norm(net),bound_norm=norm(correction),
            correction_raw_ratio=norm(correction)/norm(raw),bound_rate_norm=norm(correction)/prior['dt'],
            signed_integral=b['integral_total'],absolute_integral=b['integrals_inner_outer_near_positive_negative'][3]-b['integrals_inner_outer_near_positive_negative'][4],
            positive_negative=b['integrals_inner_outer_near_positive_negative'][3:],near_integral=dot(delta_full,near),
            global_k=dot(full_k,integral),near_k=dot(full_k,near),cumulative_accepted_integral=cumulative,
            ground=b['corrected_ground'],active_locations=b['active_constraints'],
            A_delta_max=max(map(abs,added)),M_delta_per_dt_max=max(map(abs,storage)),cancellation_identity=defect,
            original_stage_block_norms=row['refreshed_original_block_norms'],matching=row['refreshed_trace_defect']))
    require(abs(cumulative-progress['cumulative_bound_integral'])<1e-17,'joined prefix/segment budget differs')
    summary['budgets']=budgets
    summary['checks']['joined_correction_budget_no_double_count']=True
    summary['elapsed_H_over_ustar']=last['time']/summary['states']['h0']['H_over_ustar']
    summary['elapsed_initial_flow_through']=last['time']/summary['states']['h0']['flow_through_time']
    summary['classification']='unresolved evolution / hard constitutive-evaluation-domain rejection before any assessment window'
    summary['status']='diagnostics_completed';write_json(a.output/'summary.json',summary)
    print(a.output/'summary.json')


if __name__=='__main__':main()
