#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read current campaign pointers and retain verified immutable checkpoint copies."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

from continue_cls_epawf_pimple import atomic, header, process_identity
from run_cls_epawf_picard import load, require, sha, state_hash


def latest_record(path, step):
    if not path.exists():return None
    with path.open('rb') as stream:
        stream.seek(max(0,path.stat().st_size-2*1024*1024))
        lines=stream.read().decode().splitlines()
    for line in reversed(lines):
        try:record=json.loads(line)
        except json.JSONDecodeError:continue
        if record.get('global_step')==step:return record
    return None


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root=a.root.resolve();base=root.parent
    require(platform.system()=='Linux' and base==Path('/media/jay/新加卷/nektar-v2026-runs'),
        'authoritative Linux run volume required')
    roots={
        'RepairedUpstreamGraded':base/'cls-epawf-u05-graded-projection-repair-20261108/repaired-campaign/trajectory',
        'UpstreamOne':base/'cls-epawf-upstream-branches-20261107/selected-U05-one/trajectory',
        'UniformG':base/'cls-epawf-height-sampler-20261107/selected-G/trajectory',
        'UniformH':base/'cls-epawf-height-sampler-20261107/selected-H/trajectory',
        'BudgetFaithful2G':root/'faithful-2G-campaign/trajectory',
        'BudgetFaithfulRefined2G':root/'selected-refined-2G/trajectory',
        'BudgetRegularized05G':root/'selected-regularized-05G/trajectory',
        'BudgetRegularizedSaturated05G':root/'selected-regularized-saturated-05G/trajectory',
        'BudgetRegularizedSaturated05O':root/'selected-regularized-saturated-05O/trajectory',
        'BudgetRefinedRegularized2G':root/'selected-refined-regularized-2G/trajectory'}
    stamp=datetime.now(timezone.utc);directory=root/('status-capture-'+stamp.strftime('%Y%m%dT%H%M%S%fZ'))
    directory.mkdir();out=dict(schema=1,captured_UTC=stamp.isoformat(),root=str(root),
        production_governing_solves_added=0,scope='read-only current pointers; immutable verified checkpoint bytes; no retry or migration',
        library_verification='Complete closures verified in the separately dated preexisting/new campaign snapshots; this capture rechecks executable bytes.',cases={})
    for name,campaign in roots.items():
        captured=(campaign/'progress.json').read_bytes();progress=json.loads(captured,parse_int=float)
        saved=directory/name;saved.mkdir();(saved/'progress.json').write_bytes(captured)
        checkpoint=Path(progress['checkpoint']);payload=(checkpoint/'State.bin').read_bytes()
        require(hashlib.sha256(payload).hexdigest()==progress['checkpoint_sha256'],'checkpoint publication/hash differs: '+name)
        (saved/'checkpoint').mkdir();(saved/'checkpoint/State.bin').write_bytes(payload)
        receipt=header(saved/'checkpoint',True)
        require(receipt['steps']==progress['accepted_steps'] and receipt['time']==progress['accepted_time'],
            'checkpoint/frontier differs: '+name)
        driver=process_identity(progress['pid']);expected=progress.get('native_identity')
        native=process_identity(expected['pid']) if expected else None
        supervisor_status=load(campaign.parent/'status.json');supervisor_expected=supervisor_status.get('supervisor_identity')
        supervisor=process_identity(supervisor_expected['pid']) if supervisor_expected else None
        image=Path(progress['executable']);require(sha(image)==progress['executable_sha256'],'runtime image changed: '+name)
        current=dict(status=progress['status'],accepted_steps=receipt['steps'],accepted_time=receipt['time'],
            next_target=progress.get('next_target'),completed_targets=progress['completed_targets'],failure=progress.get('failure'),
            source_progress=str(campaign/'progress.json'),immutable_progress=str(saved/'progress.json'),
            progress_sha256=sha(saved/'progress.json'),source_checkpoint=str(checkpoint),checkpoint=str(saved/'checkpoint'),
            checkpoint_sha256=receipt['sha256'],coupled_state_sha256=state_hash(receipt['coupled_state']),
            native_history_sha256=[state_hash(v) for v in receipt['history']],native_history_times_sha256=state_hash(receipt['history_times']),
            executable=str(image),executable_sha256=progress['executable_sha256'],
            driver_identity=driver,driver_matches=bool(driver and str(driver['start'])==str(progress.get('process_start'))),
            native_identity=native,native_matches=bool(native and native==expected),
            supervisor_identity=supervisor,supervisor_matches=bool(supervisor and supervisor==supervisor_expected),
            cumulative_bound_integral=progress.get('cumulative_bound_integral'),
            case_family=progress.get('case_family'),case_definition=progress.get('case_definition'),
            dt=receipt['dt'],physical_time_in_H_over_ustar=receipt['time']/(progress.get('case_definition',{}).get('H',.5)/.4))
        record=latest_record(campaign/'steps.jsonl',receipt['steps'])
        if record:
            require(record['checkpoint_sha256']==receipt['sha256'],'last accepted record/checkpoint differs: '+name)
            attempt=record['attempt'];bound=attempt.get('k_bound',{})
            current['last_accepted_step']=dict(
                solves={k:attempt.get(k) for k in ('outer_sweeps','momentum_predictor_solves','pressure_system_solves',
                    'epsilon_transport_solves','k_transport_solves','bound_calls','additional_governing_equation_solves')},
                bound_activations=bound.get('activations'),minimum_k=bound.get('corrected_minimum'),
                continuous_ground=bound.get('corrected_ground'),raw_mass_norm=bound.get('mass_norm_raw_step'),
                bound_mass_norm=bound.get('mass_norm_correction'),accepted_mass_norm=bound.get('mass_norm_accepted_step'),
                seconds=attempt.get('seconds'),refreshed_matching_defect=attempt.get('refreshed_trace_defect'))
        if progress['status']=='running':require(current['driver_matches'],'live driver identity differs: '+name)
        out['cases'][name]=current
    a.output.parent.mkdir(parents=True,exist_ok=True);atomic(a.output,out)
    atomic(directory/'status.json',out);print(a.output)


if __name__=='__main__':main()
