#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Index completed accepted-AWF artifacts without advancing or repairing fields.

The durable adapter exports real immutable diagnostic objects and a deliberately
incomplete, nonrestartable analysis-only attempt view. It is not a native solver
report, checkpoint, restart migration, or accepted physical time step.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import os

from analyse_cls_accepted_awf import report
from continue_cls_epawf_pimple import atomic
from run_cls_accepted_awf import CASES, FAMILY
from run_cls_epawf_picard import load, require, sha, state_hash


def receipt(path):
    return dict(path=str(path.resolve()),sha256=sha(path),bytes=path.stat().st_size)


def accepted_records(payload, through):
    records = {}
    lines=payload.decode().splitlines()
    if payload and not payload.endswith(b'\n'):
        lines=lines[:-1]  # An in-flight later append is not an accepted record.
    for line in lines:
        row = json.loads(line)
        n = int(row['global_step'])
        if n > through:
            continue
        require(n not in records,'duplicate source accepted-step identity')
        require(row['attempt']['scheduled_complete'] and not row['attempt']['failure'],
            'unaccepted attempt in source native checkpoint ledger')
        records[n] = row
    require(sorted(records)==list(range(1,through+1)),
        'durable adapter requires a contiguous native trajectory from actual h0')
    return [records[n] for n in sorted(records)]


def durable_case(case, campaign, target):
    progress_file=campaign/'progress.json';progress_bytes=progress_file.read_bytes()
    progress=json.loads(progress_bytes,parse_int=float)
    require(progress['case_family']==FAMILY and progress['accepted_awf_case']==case and
        progress.get('prefix_steps',0)==0,'adapter requires the declared fresh-h0 native campaign')
    completed=[]
    for segment in progress['segments']:
        path=Path(segment['checkpoint']).parent/'case.evolution.json.gz'
        if path.exists() and not segment['failure']:
            completed.append((int(segment['end_step']),path))
    require(completed,'no immutable completed native block to observe')
    completed.sort();last=completed[-1][0]
    chosen={8,16,32,80,160,last}
    source_initial=completed[0][1];native_initial=report(source_initial)
    require(native_initial['wall_policy']=='CLSCalibratedBudgetAWF' and
        native_initial['accepted_states'],'wrong durable native formulation/h0')
    initial=native_initial['accepted_awf_initial']
    h0=state_hash(native_initial['accepted_states'][0])
    definition=load(source_initial.parent/'definition.json')
    require(definition['case']==case and definition['fresh_native_initialization'],
        'first durable block is not actual fresh native initialization')
    target.mkdir(exist_ok=False);(target/'snapshots').mkdir();(target/'trajectory').mkdir()
    (target/'source-progress.json').write_bytes(progress_bytes)
    snapshots={};sources=[receipt(source_initial)]
    final=None
    for n,path in completed:
        if n not in chosen:
            continue
        native=report(path)
        require(native['steps']==n and not native['failure'],'immutable native block identity differs')
        diagnostic=native.get('accepted_awf_final')
        extracted=target/'snapshots'/('native-step-%06d.json'%n)
        atomic(extracted,diagnostic or dict(valid=False,failure='native optional diagnostic unavailable'))
        snapshots[str(n)]=dict(time=native['time'],exit_code=0,path=str(extracted),sha256=sha(extracted),
            interpretation='exact extraction of a completed immutable native report; zero PDE solves',
            original_native_report=receipt(path))
        sources.append(receipt(path))
        if n==last:
            final=native
    require(final is not None,'latest native diagnostic missing')
    ledger_file=campaign/'steps.jsonl';ledger_bytes=ledger_file.read_bytes()
    (target/'trajectory/source-native-checkpoint-ledger.jsonl').write_bytes(ledger_bytes)
    records=accepted_records(ledger_bytes,last)
    # Only analyse_cls_accepted_awf consumes pimple_attempts from this view.
    # State, histories and checkpoint fields are intentionally absent.
    view=dict(schema='accepted-awf-analysis-attempt-view-v1',is_native_solver_report=False,
        may_restart=False,production_governing_solves_added=0,
        interpretation='analysis-only view of exact accepted native checkpoint-ledger records',
        source_ledger=receipt(target/'trajectory/source-native-checkpoint-ledger.jsonl'),
        original_live_ledger=str(ledger_file),through_native_step=last,
        native_record_identities=[dict(global_step=r['global_step'],checkpoint_sha256=r['checkpoint_sha256']) for r in records],
        pimple_attempts=[r['attempt'] for r in records])
    atomic(target/'trajectory/case.evolution.json',view)
    result=dict(case=case,status='durable_native_campaign_observation',
        analysis_adapter=True,may_restart=False,definition=definition,actual_h0_sha256=h0,
        initialization=initial,initialization_work=native_initial['initialization_lift_calls_iterations'],
        final=final.get('accepted_awf_final'),steps=last,time=final['time'],failure=None,snapshots=snapshots,
        campaign_status=progress['status'],campaign_failure=progress.get('failure'),
        campaign_frontier=dict(steps=progress['accepted_steps'],
            time=progress['accepted_time'],checkpoint=progress['checkpoint'],checkpoint_sha256=progress['checkpoint_sha256']),
        source_revision=progress['source_revision'],frozen_executable_sha256=progress['executable_sha256'],
        native_sources=sources,source_progress=receipt(target/'source-progress.json'),
        original_live_progress=str(progress_file),
        cumulative_correction_integral=sum(r['attempt'].get('k_bound',{}).get('integral_total',0.) for r in records))
    atomic(target/'result.json',result)
    return dict(kind='durable-native-analysis-adapter',case=case,path=str(target),
        result=receipt(target/'result.json'),attempt_view=receipt(target/'trajectory/case.evolution.json'),
        native_sources=sources,through_step=last,through_time=final['time'],may_restart=False)


def mapping(text):
    name,path=text.split('=',1)
    require(name in CASES,'unknown declared accepted-AWF case')
    return name,Path(path).resolve()


def compact_diagnostic(diagnostic):
    out={k:v for k,v in diagnostic.items() if k not in ('profiles','independent_vertical_polynomial_minima')}
    if 'profiles' in diagnostic:
        out['profile_summaries']=[]
        for profile in diagnostic['profiles']:
            row={k:v for k,v in profile.items() if k not in ('profile','functional_ids','functional_phi')}
            points=profile.get('profile',[])
            row['physical_profile_rows_retained_in_full_analysis']=len(points)
            if points:row.update(ground_row=points[0],highest_row=points[-1])
            out['profile_summaries'].append(row)
    inventory=diagnostic.get('independent_vertical_polynomial_minima')
    if inventory:
        selected=dict(inventory['selected_inventory']);selected.pop('profile_branch_pattern',None)
        panels=[];ray_rows=[]
        for anchor in inventory['anchors']:
            for panel in anchor['vertical_cubic_panels']:
                panels.append(dict(panel,x=anchor['anchor'][0]))
            for profile in anchor['profiles']:
                ray_rows.extend(profile['rows'])
        out['independent_vertical_polynomial_minima_summary']=dict(
            scope=inventory['scope'],state_sha256=inventory['state_sha256'],selected_inventory=selected,
            anchors=len(inventory['anchors']),vertical_panels=len(panels),
            minimum_vertical_panel=min(panels,key=lambda x:x['minimum']) if panels else None,
            profile_ray_rows=len(ray_rows),minimum_profile_ray=min(ray_rows,key=lambda x:x['k']) if ray_rows else None,
            unregistered_profile_ray_rows=sum(not x['visited_by_selected_inventory'] for x in ray_rows))
    return out


def compact_analysis(path,output):
    full=load(path);out=dict(full);out['full_analysis_artifact']=receipt(path)
    out['compaction']='full fields/profiles, per-step hashes and inventories remain indexed on Linux; this file retains numerical summaries'
    out['cases']={}
    for case,row in full['cases'].items():
        compact=dict(row)
        if 'initial' in row:compact['initial']=compact_diagnostic(row['initial'])
        compact['observations']={}
        for n,observation in row.get('observations',{}).items():
            item=dict(observation)
            if 'diagnostic' in item:item['diagnostic']=compact_diagnostic(item['diagnostic'])
            if 'exact_discrete_work_and_correction' in item:
                work=dict(item['exact_discrete_work_and_correction']);increments=work.pop('increments',[])
                work['individual_increments_retained_in_full_analysis']=len(increments)
                if increments:
                    work['last_increment']=increments[-1]
                    work['maximum_raw_step_mass_norm']=max(x['raw_step_mass_norm'] for x in increments)
                    work['maximum_bound_mass_norm']=max(x['correction_mass_norm'] for x in increments)
                item['exact_discrete_work_and_correction']=work
            compact['observations'][n]=item
        out['cases'][case]=compact
    out['plots']=[receipt(p) for p in sorted(path.parent.glob('*.png'))]
    output.parent.mkdir(parents=True,exist_ok=True);atomic(output,out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--case-source',type=mapping,action='append',default=[],metavar='CASE=COMPLETED_CASE_DIRECTORY')
    p.add_argument('--durable-source',type=mapping,action='append',default=[],metavar='CASE=CAMPAIGN_TRAJECTORY')
    p.add_argument('--compact-analysis',type=Path,
        help='Summarize a completed analysis native.json instead; --output names the compact JSON file')
    a=p.parse_args();a.output=a.output.resolve()
    require(not a.output.exists(),'analysis suite already exists; use a new dated snapshot, never overwrite evidence')
    if a.compact_analysis:
        require(not a.case_source and not a.durable_source,'analysis compaction is separate from suite construction')
        compact_analysis(a.compact_analysis.resolve(),a.output);print(a.output);return
    names=[c for c,_ in a.case_source+a.durable_source]
    require(names and len(set(names))==len(names),'duplicate or empty analysis case inventory')
    a.output.mkdir(parents=True)
    index=dict(schema=1,family=FAMILY,captured_UTC=datetime.now(timezone.utc).isoformat(),
        adapter='read-only accepted-AWF analysis suite',production_governing_solves_added=0,may_restart=False,
        script_sha256=sha(Path(__file__)),cases={},
        legacy_companion_metadata_caveat='Some pre-fix definitions label outer_equations Original for every variant; actual model_choices and native model_variant remain authoritative')
    for case,path in a.case_source:
        result_file=path/'result.json';result=load(result_file)
        require(result['case']==case and result['definition']['family']==FAMILY,
            'completed result family/case identity differs')
        require(result['status'] not in ('running','prepared'),'unfinished wrapper result cannot be indexed as completed')
        os.symlink(path,a.output/case,target_is_directory=True)
        index['cases'][case]=dict(kind='completed-native-wrapper-symlink',source=str(path),result=receipt(result_file),
            status=result['status'],steps=result.get('steps'),time=result.get('time'))
    for case,path in a.durable_source:
        index['cases'][case]=durable_case(case,path,a.output/case)
    atomic(a.output/'adapter-index.json',index)
    print(a.output/'adapter-index.json')


if __name__=='__main__':
    main()
