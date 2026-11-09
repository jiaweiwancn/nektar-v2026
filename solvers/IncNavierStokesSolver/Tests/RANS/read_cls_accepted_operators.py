#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Strict accepted-checkpoint spatial operators on retained fields, no evolution."""
import argparse
from datetime import datetime,timezone
from pathlib import Path
import shutil
import os

from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from continue_cls_epawf_pimple import header,atomic
from run_cls_epawf_picard import load,require,sha,state_hash,THREADS
from analyse_cls_accepted_awf import compact,differences,normalized
from build_cls_accepted_awf_analysis_suite import compact_diagnostic


def component_summary(folder):
    """Use the exported physical M_ff, never the frozen transport matrix A_k."""
    import numpy as np
    from scipy.linalg import cho_factor,cho_solve
    summary=load(folder/'summary.json');mass_file=Path(summary['cases']['Regularized-0']['full_free_loads_retained']);first=load(mass_file)['diagnostics']
    mass=first['native_free_k_mass'];rows=np.array(mass['rows'],dtype=int)
    matrix=np.array(mass['matrix']);factor=cho_factor(matrix,lower=True)
    terms=('conservative_advection','k_diffusion','inner_k_production','inner_k_destruction',
        'outer_k_production','outer_k_destruction','prescribed_scalar_source','prescribed_body',
        'roughness_force','inner_momentum_stress','outer_momentum_stress','epsilon_diffusion',
        'epsilon_production','epsilon_destruction','epsilon_yap','continuity','gauge')
    out=dict(metric='actual physical consistent free-k M_ff from PicardKMass; essential rows eliminated',
        mass_source_sha256=sha(mass_file),
        free_rows=len(rows),governing_solves=0,cases={},
        signs='residual convention: production loads negative, destruction positive; norms do not retain sign',
        scope='full free-mode component loads and consistent-mass dual norms, not a conservative physical-band flux balance')
    for label,c in summary['cases'].items():
        if not label.startswith('Regularized-'):continue
        file=Path(c['full_free_loads_retained']);d=load(file)['diagnostics'];koffset=int(d['selected_state']['offsets'][2])
        require(list(rows)==d['operators']['blocks'][2]['free_rows'],'free scalar rows changed between immutable snapshots')
        values=d['physical']['weak_k_component_loads'];indices=rows-koffset
        selected={name:np.array(values[name])[indices] for name in terms}
        selected['spatial_sum']=np.array(values['spatial_sum'])[indices]
        selected['total_constrained']=np.array(values['total_constrained'])[indices]
        combined=sum(selected[name] for name in terms)
        require(np.max(abs(combined-selected['spatial_sum']))<1e-12,'physical k decomposition counted/omitted a term')
        names=list(selected);loads=np.column_stack([selected[n] for n in names]);rates=cho_solve(factor,loads)
        require(np.max(abs(matrix@rates-loads))<1e-11,'physical mass solve failed')
        norms={name:float(max(0.,loads[:,i]@rates[:,i])**.5) for i,name in enumerate(names)}
        native=d['operators']['blocks'][2]['mass_dual_dimensional'][0]
        require(abs(norms['total_constrained']-native)<1e-10*max(native,1e-6),'component metric differs from native physical residual norm')
        out['cases'][label]=dict(time=c['time'],result_sha256=sha(file),mass_dual_dimensional=norms,
            component_sum_maximum_error=float(np.max(abs(combined-selected['spatial_sum']))),
            native_total_norm_error=abs(norms['total_constrained']-native))
    atomic(folder/'component-mass-dual.json',out)
    return out


def main():
    os.environ.update(THREADS)
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',required=True,type=Path)
    p.add_argument('--old-root',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--reuse',type=Path,help='reuse immutable readbacks with identical probe and checkpoint hashes')
    p.add_argument('--handoff',type=Path,help='validated completed-boundary receipt when the old driver retired before appending its segment summary')
    p.add_argument('--steps',type=int,nargs='+',default=[0,320,640,800]);a=p.parse_args();a.output.mkdir(exist_ok=False)
    campaign=a.old_root/'selected-refined-regularized-2G/trajectory';progress=load(campaign/'progress.json')
    # The committed append frontier may include an in-flight later line. Only
    # chosen completed native identities are consumed, never repaired here.
    ledger={}
    with (campaign/'steps.jsonl').open() as stream:
        import json
        for line in stream:
            try:r=json.loads(line)
            except json.JSONDecodeError:break
            if r['global_step'] in a.steps:ledger[int(r['global_step'])]=r
    previous=load(a.reuse/'summary.json') if a.reuse else {}
    rows={};data={};out=dict(schema=1,started_UTC=datetime.now(timezone.utc).isoformat(),status='started',
        probe_sha256=sha(a.probe),governing_solves=0,physical_commits=0,native_history_advanced=False,
        reporting_breakpoints=[.0025,.005,.02,.04],cases={},windows={},unavailable_steps=[])
    originals=a.old_root/'refined-2G/Accepted-AWF-2G/initial'
    selections=[('Original-h0',originals/'checkpoint',originals/'case.xml')]
    handoff=load(a.handoff/'handoff.json') if a.handoff else None
    if handoff:
        require(handoff['status']=='handoff_completed_with_observed_step' and Path(handoff['parent_campaign'])==campaign and
            handoff['complete_corrected_audit']['accepted_awf_hard_admissible'][0]==1 and
            handoff['isolated_restart_equivalence']['bitwise_checkpoint_payload'],'incomplete or unrelated handoff readback receipt')
        out['handoff_receipt']=dict(path=str(a.handoff/'handoff.json'),sha256=sha(a.handoff/'handoff.json'))
    for n in a.steps:
        if n==0:
            base=a.old_root/'initial-history-replication/refined-regularized-2G'
            selections.append(('Regularized-0',base/'checkpoint',base/'case.xml'));continue
        segments=[v for v in progress['segments'] if int(v['end_step'])==n and not v['failure']]
        if not segments and handoff and n==handoff['parent_step']:
            chosen=a.handoff/'parent/checkpoint';h=header(chosen,True);prefix=load(a.handoff/'parent/prefix.evolution.json')
            require(h['steps']==n==progress['accepted_steps']==prefix['steps'] and h['time']==prefix['time'] and
                h['sha256']==handoff['parent_checkpoint_sha256']==progress['checkpoint_sha256']==ledger[n]['checkpoint_sha256'],
                'handoff checkpoint/ledger identity differs')
            require(state_hash(h['coupled_state'])==state_hash(prefix['state']) and
                [state_hash(v) for v in h['history']]==[state_hash(v) for v in prefix['native_history']] and
                state_hash(h['history_times'])==state_hash(prefix['native_history_times']),
                'handoff native reference/history lineage differs')
            selections.append(('Regularized-'+str(n),chosen,Path(progress['segment'])/'case.xml'));continue
        if not segments:out['unavailable_steps'].append(n);continue
        base=Path(segments[0]['checkpoint']).parent;candidates=[]
        for file in base.glob('checkpoint*/State.bin'):
            if header(file.parent)['steps']==n:candidates.append(file.parent)
        # Native finalization can publish two byte-identical aliases for the
        # block endpoint.  Identity is the accepted checkpoint hash, not the
        # number of directory names that refer to those same bytes.
        require(candidates and len({sha(c/'State.bin') for c in candidates})==1,
            'missing or conflicting retained native checkpoints for step '+str(n))
        chosen=Path(segments[0]['checkpoint'])
        require(chosen in candidates and n in ledger and sha(chosen/'State.bin')==ledger[n]['checkpoint_sha256'],
            'native checkpoint is not bound to accepted segment and ledger')
        selections.append(('Regularized-'+str(n),chosen,base/'case.xml'))
    for label,source,xml in selections:
        folder=a.output/label;folder.mkdir();cp=folder/'checkpoint';cp.mkdir();shutil.copy2(source/'State.bin',cp/'State.bin')
        require(sha(cp/'State.bin')==sha(source/'State.bin'),'checkpoint copy differs');native=checkpoint(cp)
        packet=dict(common_packet(cp),operation='selected_original_operator',state=native['coupled'],history=native['coupled'],
            time=native['time'],alpha=0.,inner_stress_model='Linear',witnesses=[[-3.,0.,3.6880433085817326]],
            reporting_breakpoints=out['reporting_breakpoints'],export_k_mass=label=='Regularized-0')
        reused=previous.get('cases',{}).get(label)
        if reused and previous.get('probe_sha256')==sha(a.probe) and reused['checkpoint_sha256']==sha(cp/'State.bin'):
            file=Path(reused['full_free_loads_retained']);result=load(file)
            receipt=dict(output=str(file),output_sha256=sha(file),reuse='immutable identical-probe/checkpoint readback',
                source_summary=str(a.reuse/'summary.json'),source_summary_sha256=sha(a.reuse/'summary.json'))
        else:result,receipt=run_probe(a.probe,xml,packet,folder/'readback')
        d=result['diagnostics'];physical=d['physical'];op=d['operators']
        require(result['governing_solves']==0 and result['initialization_work']==[0,0] and not result['native_history_advanced'],'operator observation advanced/restored fields')
        require(op['selected_native_residual_maximum_difference']<1e-9 and op['unaffected_physical_term_maximum_difference']==0,'inconsistent operator substitution')
        require(op['selected_continuity_gauge_loads']==op['counterfactual_continuity_gauge_loads'],'counterfactual changed continuity/gauge')
        if label=='Original-h0':require(all(max(map(abs,b['difference_load']),default=0)==0 for b in op['blocks']),'Original counterfactual not exact no-op')
        cov=d['selected_state'];c=dict(checkpoint=str(source),checkpoint_copy=str(cp),checkpoint_sha256=sha(cp/'State.bin'),
            time=native['time'],step=native['steps'],state_sha256=state_hash(native['coupled']),
            native_history_sha256=[state_hash(v) for v in native['history']],native_history_times_sha256=state_hash(native['history_times']),
            receipt=receipt,physical=compact_diagnostic(compact(physical)),complete_audit=result.get('complete_audits',{}).get('selected_state'),
            selected_covariance=cov['minimum_outer'],Original_covariance=cov['minimum_same_state_Original_outer'],
            operators={k:v for k,v in op.items() if k!='blocks'},
            blocks=[{k:v for k,v in b.items() if k not in ('free_rows','selected_load','same_state_Original_load','difference_load')} for b in op['blocks']],
            full_free_loads_retained=receipt['output'])
        out['cases'][label]=c;data[label]=physical;atomic(a.output/'summary.json',out)
    if 'Regularized-0' in data:
        for label,d in data.items():
            if label.startswith('Regularized-'):
                out['cases'][label]['drift_from_actual_native_h0']=differences(d,data['Regularized-0'])
        ordered=sorted(int(label.split('-')[-1]) for label in data if label.startswith('Regularized-'))
        windows=set(zip(ordered,ordered[1:])) | {(0,320),(320,640),(0,800)}
        # Fixed requested endpoints yield adjacent physical windows; include
        # every second endpoint to expose window-length sensitivity as well.
        windows.update(zip(ordered,ordered[2:]))
        for first,last in sorted(windows):
            left,right='Regularized-'+str(first),'Regularized-'+str(last)
            if left in data and right in data:
                change=differences(data[right],data[left]);out['windows'][str(first)+'-'+str(last)]=dict(start=data[left]['time'],end=data[right]['time'],
                    dimensional_field_change=change,common_normalized_field_change=normalized(change),
                    traction_change=data[right]['mean_wall_stress']-data[left]['mean_wall_stress'],
                    interpretation='fixed physical common-panel endpoints; instantaneous source rates are not time-integrated budgets')
    out.update(status='passed',completed_UTC=datetime.now(timezone.utc).isoformat());atomic(a.output/'summary.json',out)
    component_summary(a.output)
    print(a.output/'summary.json')


if __name__=='__main__':main()
