#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read an immutable prefix of the live V2 campaign without advancing it."""
import argparse
import gzip
import json
import math
from pathlib import Path
import shutil
import time
from run_cls_epawf_picard import load,run,require,write_json,state_hash,H0,sha
from run_cls_epawf_pimple import identity
from continue_cls_epawf_pimple import header
from compare_cls_epawf_pimple import differences


def read_report(folder):
    p=folder/'case.evolution.json'
    for attempt in range(5):
        try:
            if p.exists():return load(p)
            with gzip.open(str(p)+'.gz','rt') as f:return json.load(f,parse_int=float)
        except (FileNotFoundError,EOFError,json.JSONDecodeError):
            if attempt==4:raise
            time.sleep(.1)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign',type=Path,required=True);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--through-step',type=int)
    p.add_argument('--steps',type=int,nargs='+',default=[4,8,12,16,32])
    a=p.parse_args();a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    root=a.campaign.resolve()/'trajectory';progress=load(root/'progress.json')
    text=(root/'steps.jsonl').read_text();lines=text.splitlines();records={}
    if text and not text.endswith('\n'):lines=lines[:-1]
    for line in lines:
        r=json.loads(line);require(r['global_step'] not in records,'duplicate accepted ledger entry');records[r['global_step']]=r
    cutoff=a.through_step or int(progress['accepted_steps']);require(cutoff<=progress['accepted_steps'],'requested prefix is not accepted')
    records={k:v for k,v in records.items() if k<=cutoff}
    checkpoint=Path(records[cutoff]['checkpoint'])
    shutil.copytree(checkpoint,a.output/'checkpoint-at-readback')
    saved=header(a.output/'checkpoint-at-readback',True)
    require(saved['steps']==cutoff and saved['sha256']==records[cutoff]['checkpoint_sha256'],'recorded checkpoint no longer retained; choose a fixed target')
    require(len(records)==cutoff and set(records)==set(range(1,cutoff+1)),'noncontiguous fresh-h0 prefix')
    rows=[r['attempt'] for r in records.values()]
    require(all(v['scheduled_complete'] and v['governing_solves_complete'] for v in rows),'accepted incomplete schedule')
    for v in rows:
        require(all(v[k]==1 for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','bound_calls')),'changed schedule')
        require(v['additional_governing_equation_solves']==0,'additional governing solve')
    initial=read_report(Path(records[1]['segment']))['accepted_states'][0];require(state_hash(initial)==H0,'initial h0 differs')
    summary=dict(schema=1,status_at_readback=progress['status'],last_verified_step=cutoff,last_verified_time=records[cutoff]['time'],
        progress=progress,probe=identity(a.probe.resolve()),states={},windows=[],budgets=[],
        classification='unresolved evolution; T=.02 is only a first milestone, not a long-time assessment',
        ranks=1,threads=1,governing_solves={k:sum(v[k] for v in rows) for k in ('momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves','additional_governing_equation_solves')})
    summary['completed_native_segments']=[]
    for segment in sorted({v['segment'] for v in records.values()}):
        folder=Path(segment)
        if not (folder/'case.evolution.json').exists() and not (folder/'case.evolution.json.gz').exists():continue
        report=read_report(folder)
        attempts=report['pimple_attempts'];accepted_here=len(report['accepted_states'])-1
        require(report['implicit_callback_count']==len(attempts),'native callback counter mismatch')
        if not report['failure']:require(len(attempts)==accepted_here,'hidden implicit callback in completed segment')
        summary['completed_native_segments'].append(dict(path=segment,callbacks=report['implicit_callback_count'],
            accepted=accepted_here,global_end_step=report['steps'],failure=report['failure'],
            initialization_k_bound=report.get('initial_k_bound_check'),initialization_lift_calls_iterations=report.get('initialization_lift_calls_iterations'),
            qualification_events=[{k:v for k,v in e.items() if k not in ('state','history','audit')} for e in report.get('physical_range_events',[])]))
    sampled={}
    def replay(name,state,time,xml,witnesses=()):
        folder=a.output/name;folder.mkdir();shutil.copy2(xml,folder/'case.xml')
        request=dict(operation='contact',state=state,history=state,time=time,alpha=0.,h=[],
            inner_stress_model='SaturatedR2',require_history_h0=False,export_mass=name=='initial',isolate=False,witnesses=list(witnesses))
        write_json(folder/'input.json',request)
        receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],folder,'probe.log',3600)
        write_json(folder/'receipt.json',receipt);require(receipt['exit_code']==0,'readback failed '+name)
        contact=load(folder/'result.json');require(contact.get('retained_state_unchanged'),'readback changed state')
        request['operation']='samples';write_json(folder/'samples-input.json',request)
        receipt=run([str(a.probe.resolve()),'case.xml','samples-input.json','samples.json'],folder,'samples.log',3600)
        require(receipt['exit_code']==0,'sample readback failed '+name)
        samples=load(folder/'samples.json');samples['initial']=initial;sampled[name]=samples
        summary['states'][name]={k:v for k,v in contact.items() if k not in ('free_k_mass','free_k_rows','offsets','k_integral_load','near_k_integral_load','scan')}
        weights=samples['physical_samples']['wall_weights'];stress=samples['physical_samples']['wall_stress']
        summary['states'][name]['mean_wall_stress']=math.fsum(w*t for w,t in zip(weights,stress))/math.fsum(weights)
        summary['states'][name]['rms_wall_stress']=math.sqrt(math.fsum(w*t*t for w,t in zip(weights,stress))/math.fsum(weights))
        write_json(a.output/'summary.partial.json',summary)
        return contact
    xml=Path(records[1]['segment'])/'case.xml';mass_spec=replay('initial',initial,0.,xml)
    selected=[n for n in sorted(set(a.steps)) if n<=cutoff]
    for n in selected:
        row=records[n];packet=Path(row['attempt']['k_snapshot'])
        if packet.exists():q=load(packet)
        else:
            report=read_report(Path(row['segment']));start=int(report['steps'])-len(report['accepted_states'])+1
            q=dict(corrected=report['accepted_states'][n-start],time=row['time'])
        witnesses=sorted({tuple(v['coordinate']) for v in row['attempt']['k_bound']['active_constraints'] if v['coordinate'][2]<=mass_spec['scales_U_K_E_H_L_volume'][3]})
        replay('step'+str(n),q['corrected'],q['time'],Path(row['segment'])/'case.xml',witnesses)
    names=['initial']+['step'+str(n) for n in selected if 'step'+str(n) in sampled]
    for x,y in zip(names,names[1:]):
        first=dict(sampled[x],time=sampled[y]['time'])
        summary['windows'].append(dict(start=sampled[x]['time'],end=sampled[y]['time'],
            physical_L2=differences(first,sampled[y])))
    # Actual physical consistent-mass functionals exported by the same native spaces.
    import numpy as np
    mass=np.array(mass_spec['free_k_mass']);free=np.array(mass_spec['free_k_rows'],dtype=int);off=int(mass_spec['offsets'][2])
    integral=np.array(mass_spec['k_integral_load']);near=np.array(mass_spec['near_k_integral_load'])
    dual=np.linalg.solve(mass.T,np.array([integral[free-off],near[free-off]]).T).T
    native_net_norms={}
    for entry in summary['completed_native_segments']:
        report=read_report(Path(entry['path']));states=report['accepted_states'];start=int(report['steps'])-len(states)+1
        for i,(before,after) in enumerate(zip(states,states[1:]),1):
            delta=np.array(after)-np.array(before);step=start+i
            if step>cutoff:continue
            essential=set(range(off,off+len(integral)))-set(free)
            require(all(delta[j]==0 for j in essential),'prescribed k changed; free-mass increment norm requires full lifting accounting')
            d=delta[free];native_net_norms[step]=math.sqrt(d@(mass@d))
    cumulative=0.;positive=0.;negative=0.;near_total=0.
    for n,row in sorted(records.items()):
        v=row['attempt'];b=v['k_bound'];cumulative+=b['integral_total']
        ints=b['integrals_inner_outer_near_positive_negative'];positive+=ints[3];negative+=ints[4];near_total+=ints[2]
        out=dict(step=n,time=row['time'],correction_mass_norm=b['mass_norm_correction'],raw_increment_mass_norm=b['mass_norm_raw_step'],
            correction_raw_ratio=b['correction_to_step_ratio'],correction_rate_mass_norm=b['mass_norm_correction']/.000625,
            signed_integral=b['integral_total'],signed_integral_rate=b['integral_total']/.000625,
            near_integral=ints[2],positive_integral=ints[3],negative_integral=ints[4],absolute_integral=ints[3]-ints[4],
            minimum=b['corrected_minimum'],raw_minimum=b['raw_minimum'],floor=b['floor'],ground=b['corrected_ground'],separations=b['separation_passes'],
            bound_calls=b['calls'],activations=b['activations'],mass_applications=b['mass_factor_applications'],kkt_solves=b['kkt_solves'],
            active_set_iterations=b['active_set_iterations'],active_constraints=b['active_constraints'],
            inventory_counts=b['inventory_counts'],constraint_families=b['discovered_constraint_families'],
            raw_frozen_residual=v['k_frozen_linear_residual'],corrected_frozen_residual=v['k_corrected_frozen_linear_residual'],
            added_frozen_residual=v['k_added_frozen_residual'],storage_reaction=v['k_added_storage_residual'],
            cumulative_integral=cumulative,matching=v['refreshed_trace_defect'],continuity=v['corrected_continuity_gauge'],
            B_min=v['k_domain']['minimum_b'],lambda_min=v['k_domain']['minimum_lambda'])
        if n in native_net_norms:out['corrected_increment_mass_norm']=native_net_norms[n]
        packet=Path(v['k_snapshot'])
        if packet.exists():
            q=load(packet);raw=np.array(q['raw']);corrected=np.array(q['corrected']);before=np.array(q['pre_k']);dk=(corrected-raw)[free];net=(corrected-before)[free]
            value=math.sqrt(net@(mass@net))
            if n in native_net_norms:require(abs(value-native_net_norms[n])<1e-12,'native accepted k increment differs from saved correction')
            out['corrected_increment_mass_norm']=value
            out['global_k']=float(integral@corrected[off:off+len(integral)]);out['near_k']=float(near@corrected[off:off+len(near)])
            out['cancellation_identity']=float(np.max(np.abs(net-(raw-before)[free]-dk)))
            out['pre_k_minus_native_history_max']=float(np.max(np.abs((before-np.array(q['history']))[off:off+len(integral)])))
            matrix=np.frombuffer(Path(q['matrix_file']).read_bytes(),dtype='<f8').reshape(len(free),len(free));h=1/q['alpha']
            rd=np.array(q['raw_increment']);input_load=-h*(dual@np.array(q['transport_residual']))
            implicit=-h*(dual@(matrix@rd-mass@rd/h));repair=np.array([integral[free-off]@dk,near[free-off]@dk]);actual=np.array([integral[free-off]@net,near[free-off]@net])
            out['exact_discrete_global_near_budget']=dict(input_load=input_load.tolist(),implicit_increment=implicit.tolist(),bound=repair.tolist(),actual=actual.tolist(),closure=(actual-input_load-implicit-repair).tolist())
            require(out['cancellation_identity']<1e-15 and max(abs(x) for x in out['exact_discrete_global_near_budget']['closure'])<1e-12,'increment budget identity failed')
        summary['budgets'].append(out)
    if cutoff==progress['accepted_steps']:require(abs(cumulative-progress['cumulative_bound_integral'])<1e-15,'live prefix cumulative correction differs')
    summary['cumulative_correction']=dict(net=cumulative,positive=positive,negative=negative,absolute=positive-negative,near=near_total)
    summary['window_integrated_budgets']=[]
    for x,y in zip(names,names[1:]):
        start=sampled[x]['time'];end=sampled[y]['time'];selected_budget=[v for v in summary['budgets'] if start+1e-13<v['time']<=end+1e-13]
        physical=[(end-start)*(a+b)/2 for a,b in zip(summary['states'][x]['physical_production_destruction_global_near'],summary['states'][y]['physical_production_destruction_global_near'])]
        row=dict(start=start,end=end,correction_net=sum(v['signed_integral'] for v in selected_budget),
            correction_positive=sum(v['positive_integral'] for v in selected_budget),correction_negative=sum(v['negative_integral'] for v in selected_budget),
            correction_absolute=sum(v['absolute_integral'] for v in selected_budget),correction_near=sum(v['near_integral'] for v in selected_budget),
            physical_production_destruction_global_near=physical,
            physical_time_quadrature='trapezoid at retained endpoints; not an exactly closed CG control-volume balance')
        if all('exact_discrete_global_near_budget' in v for v in selected_budget):
            row['exact_frozen_discrete']={k:np.sum([v['exact_discrete_global_near_budget'][k] for v in selected_budget],axis=0).tolist() for k in ('input_load','implicit_increment','bound','actual','closure')}
        summary['window_integrated_budgets'].append(row)

    summary['checkpoint_identity_at_readback']=dict(path=str(checkpoint),snapshot=str(a.output/'checkpoint-at-readback'),sha256=saved['sha256'],steps=saved['steps'],time=saved['time'])
    summary['status']='analysis_completed';write_json(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
