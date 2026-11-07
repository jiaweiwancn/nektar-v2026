#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable status/checkpoint copies and compact projection-repair evidence."""
import argparse
import datetime
import hashlib
import gzip
import json
from pathlib import Path
from capture_cls_epawf_upstream_status import process_command,lock_record
from continue_cls_epawf_pimple import load_progress,process_identity,header
from run_cls_epawf_picard import load,write_json,require,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--old-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);r=a.root.resolve();index=[]
    def pin(path,role):
        q=dict(path=str(path.resolve()),sha256=sha(path),bytes=path.stat().st_size,role=role);index.append(q);return q
    out=dict(captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),root=str(r),
        base='a0ccc1c7febc3ed18efb5f4e5fc6d1116778f22f',source=load(r/'source-validation.json'),
        checkpoint_compatibility='unchanged physical mass metric, floor/reserve, V3 constraint policy, equations, dt and history; versioned finite-QP numerical backend separately recorded',
        scope='declared V3 physical inventory plus continuous ground; not global volume positivity',replays={},tests={},campaigns={})
    analysis=load(r/'replay-analysis.json')
    for name,d in analysis['replays'].items():
        b=d.pop('bound');d['work_and_feasibility']={k:b[k] for k in
            ('success','failure','work_controls','separation_passes','inventory_rebuilds','active_set_iterations','kkt_solves',
             'maximum_work_one_add','finite_pool_scans','finite_pool_reentries','exact_cut_reuses','finite_pool_complete',
             'mass_factorizations','mass_factor_applications','seconds_total','seconds_algebra','seconds_inventory',
             'floor','arithmetic_reserve','raw_minimum','corrected_minimum','raw_ground','corrected_ground',
             'primal_floor_defect','ground_floor_defect','stationarity_mass_dual','complementarity','dual_feasibility_defect',
             'primal_dual_gap_working_set','finite_pool_normalized_violation','raw_sha256','corrected_sha256') if k in b}
        out['replays'][name]=d
        if (r/name/'manifest.json').exists():pin(r/name/'manifest.json','strict saved-candidate replay manifest')
        for f in ('result.json','pool.json','pool.json.Mff.bin','pool.json.Cwhite.bin','trace.jsonl','input.json','case.xml'):
            pin(r/name/f,'immutable exact-candidate diagnostic '+name)
    out['cross_pools']=analysis['cross_pools'];out['between_replay_mass_difference']=analysis['between_replay_mass_difference']
    out['cut_pool_interpretation']=analysis['interpretation']
    for name in ('native-first','native-earlier','native-transactions'):
        f=r/name/'summary.json'
        if not f.exists():out['tests'][name]=dict(status='running_or_not_published');continue
        q=load(f);out['tests'][name]={k:v for k,v in q.items() if k not in ('solver','native_attempt')}
        if 'native_attempt' in q:
            t=q['native_attempt'];out['tests'][name]['step65_counts']={k:t[k] for k in
                ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves',
                 'bound_calls','additional_governing_equation_solves','scheduled_complete')}
        pin(f,'focused native verification receipt')
        for f in (r/name).glob('*/case.evolution.json'):pin(f,'native verification state and all histories')
    f=r/'replay-budget-stop/coherence.json'
    if f.exists():out['tests']['interrupted_finite_pool']=load(f);pin(f,'budget-stop coherent failed-primal/dual verification')
    out['tests']['ctest']=(r/'ctest-pool-v1.log').read_text();pin(r/'ctest-pool-v1.log','focused Linux CTest results')
    paths={'repaired_graded':r/'repaired-campaign/trajectory','original_graded_terminal':a.old_root/'selected-U05-graded/trajectory',
        'upstream_one':a.old_root/'selected-U05-one/trajectory',
        'uniform_G':r.parent/'cls-epawf-height-sampler-20261107/selected-G/trajectory',
        'uniform_H':r.parent/'cls-epawf-height-sampler-20261107/selected-H/trajectory'}
    for name,path in paths.items():
        if not (path/'progress.json').exists():continue
        for retry in range(20):
            q=load_progress(path/'progress.json');cp=Path(q['checkpoint'])
            try:payload=(cp/'State.bin').read_bytes()
            except FileNotFoundError:continue
            if hashlib.sha256(payload).hexdigest()==q['checkpoint_sha256']:break
        else:raise RuntimeError('no coherent published checkpoint for '+name)
        folder=a.output/name;folder.mkdir();copy=folder/'checkpoint';copy.mkdir();(copy/'State.bin').write_bytes(payload)
        write_json(folder/'progress.json',q);pin(copy/'State.bin','immutable dated accepted checkpoint: '+name)
        driver=process_identity(q['pid']) if q.get('pid') else None
        expected=q.get('native_identity');native=process_identity(expected['pid']) if expected else None
        identity=load(path/'runtime/identity.json');pin(path/'runtime/identity.json','frozen runtime closure: '+name)
        require(sha(Path(q['executable']))==q['executable_sha256'],'campaign executable changed')
        for library,digest in identity['frozen']['libraries'].items():require(sha(Path(library))==digest,'frozen linked library changed')
        c=dict(status=q['status'],accepted_steps=q['accepted_steps'],accepted_time=q['accepted_time'],checkpoint_copy=header(copy),
            original_checkpoint=str(cp),checkpoint_sha256=q['checkpoint_sha256'],completed_targets=q['completed_targets'],
            next_target=q.get('next_target'),failure=q.get('failure'),cumulative_bound_integral=q['cumulative_bound_integral'],
            bound_solver=q.get('bound_solver','WorkingSet'),backend_recorded_in_progress='bound_solver' in q,
            driver=driver,native=native,native_matches=bool(native and native==expected),driver_command=process_command(driver),
            native_command=process_command(native),campaign_lock=lock_record(path/'.campaign.lock',[driver,native]),
            executable_sha256=q['executable_sha256'],source_revision=q.get('source_revision'))
        status=path.parent/'status.json'
        if status.exists():
            supervisor=load(status);expected=supervisor.get('supervisor_identity');actual=process_identity(expected['pid']) if expected else None
            c['supervisor_matches']=bool(actual and actual==expected);c['supervisor_command']=process_command(actual)
        if name=='repaired_graded':
            text=(path/'steps.jsonl').read_text()
            if text and not text.endswith('\n'):text=text[:text.rfind('\n')+1]
            records=[json.loads(v) for v in text.splitlines()]
            records=[v for v in records if v['global_step']<=q['accepted_steps']]
            ids=[int(v['global_step']) for v in records]
            require(ids==list(range(int(q['prefix_steps'])+1,int(q['accepted_steps'])+1)),'accepted ledger duplicated or has a gap')
            c['accepted_step_records']=[]
            for record in records:
                t=record['attempt'];b=t['k_bound'];require(t['bound_calls']==1 and t['additional_governing_equation_solves']==0,'hidden equation solve')
                require(all(t[k]==1 for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves','epsilon_transport_solves','k_transport_solves')),'governing solve schedule changed')
                c['accepted_step_records'].append(dict(global_step=record['global_step'],time=record['time'],
                    solve_counts={k:t[k] for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves',
                        'epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')},
                    k_linear_residual=t['k_frozen_linear_residual'],k_corrected_residual=t['k_corrected_frozen_linear_residual'],
                    added_frozen_residual=t['k_added_frozen_residual'],storage_reaction=t['k_added_storage_residual'],
                    raw_mass_norm=b['mass_norm_raw_step'],bound_mass_norm=b['mass_norm_correction'],
                    correction_raw_ratio=b['correction_to_step_ratio'],integral=b['integral_total'],
                    corrected_minimum=b['corrected_minimum'],ground=b['corrected_ground'],
                    kkt_solves=b['kkt_solves'],pool_reentries=b['finite_pool_reentries'],inventory_rebuilds=b['inventory_rebuilds'],
                    mass_factor_applications=b['mass_factor_applications'],finite_pool_complete=b['finite_pool_complete'],
                    physical_floor_defect=b['primal_floor_defect'],ground_floor_defect=b['ground_floor_defect'],
                    physical_correction_bands=b['physical_bands'],
                    correction_seconds=b['seconds_total'],refreshed_trace_defect=t['refreshed_trace_defect']))
                c['accepted_step_records'][-1]['refreshed_stage_residual_block_norms']=t['refreshed_original_block_norms']
            require(abs(q['prefix_integral']+sum(v['attempt']['k_bound']['integral_total'] for v in records)-q['cumulative_bound_integral'])<1e-14,'prefix budget counted twice or lost')
            c['prefix_integral']=q['prefix_integral'];c['prefix_counted_once']=True
            if q['segments']:
                segment=Path(q['segments'][-1]['checkpoint']).parent;report=segment/'case.evolution.json.gz'
                if report.exists():
                    with gzip.open(report,'rt') as stream:complete=json.load(stream)
                    d=complete['height_comparison_final'];pin(report,'completed native block with retained epsilon and alpha=0 readback')
                    c['completed_block_alpha_zero']=dict(time=d['time'],mean_wall_stress=d['mean_wall_stress'],
                        dimensional_spatial_blocks=d['dimensional_mass_dual_spatial_first4'],
                        common_normalized_spatial_blocks=d['common_normalized_spatial_first4'],
                        refreshed_matching_defect=d['matching_trace_defect'],qualification=complete['state_audits'][-1]['r2_qualified'][0]==1)
        out['campaigns'][name]=c
    for n in ('pool-v1/probe/identity.json','pool-v1/solver/identity.json','pool-v1/source.json',
              'source-validation.json','initial-provenance.json','build-pool-v1.log','replay-analysis.json'):
        pin(r/n,'source, build or immutable runtime identity')
    pin(a.old_root/'terminal-U05-graded/attempt.json','original failed attempt remains unchanged')
    require(sha(a.old_root/'terminal-U05-graded/checkpoint/State.bin')=='4a953092ce4e85c6a6fec9c172ce5eee83bf820c8b0600762637afd106a9abe0','original accepted checkpoint changed')
    write_json(a.output/'native.json',out);write_json(a.output/'artifact-index.json',dict(root=str(r),artifacts=index))
    print(a.output/'native.json')


if __name__=='__main__':main()
