#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Controlled accepted-boundary handoff; preserve old runtime/checkpoint lineage.

Pause the identified supervisor, then stop its driver only after the native
child has ended at its requested integer boundary. Revalidate and independently
replay the exact checkpoint before launching a separately locked successor.
"""
import argparse
from datetime import datetime,timezone
import gzip
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from benchmark_cls_awf_evaluation_reuse import replay,compare
from continue_cls_epawf_pimple import atomic,header,process_identity,CampaignLock,prefix_budget,load_progress
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from follow_cls_accepted_awf_pair import supervisor_command
from run_cls_epawf_picard import THREADS,load,require,sha,state_hash


def matches(actual,expected):
    return bool(actual and expected and int(actual['pid'])==int(expected['pid']) and str(actual['start'])==str(expected['start']))


def command_line(expected):
    return Path(f"/proc/{int(expected['pid'])}/cmdline").read_bytes().decode().replace('\0',' ')


def checked_signal(expected,kind):
    require(matches(process_identity(expected['pid']),expected),'process identity changed before signal')
    os.kill(int(expected['pid']),kind)


def live_children(pid):
    path=Path(f'/proc/{pid}/task/{pid}/children')
    return [v for p in path.read_text().split() if (v:=process_identity(int(p))) is not None] if path.exists() else []


def boundary(progress):
    return progress.get('status') not in ('numerically_rejected','targets_completed') and \
        progress.get('accepted_steps')==progress.get('segment_target') and \
        not matches(process_identity(progress.get('native_pid') or -1),progress.get('native_identity'))


def wait_parent_lock(path,timeout=30.):
    deadline=time.monotonic()+timeout
    while True:
        try:return CampaignLock(path)
        except RuntimeError as error:
            if 'campaign lock held' not in str(error) or time.monotonic()>=deadline:raise
            time.sleep(.05)


def observe_successor(a,h,first_sha,save):
    """Observe an already launched successor; a read gap never launches work."""
    deadline=time.monotonic()+a.timeout
    pointer=a.output/'trajectory/progress.json'
    while time.monotonic()<deadline:
        try:nextp=load_progress(pointer)
        except FileNotFoundError:
            time.sleep(.2);continue
        require(nextp['status']!='numerically_rejected','successor numerically rejected; preserved terminal campaign')
        if nextp['accepted_steps']>h['steps']:
            records=(a.output/'trajectory/steps.jsonl').read_text().splitlines()
            row=next(json.loads(line) for line in records if json.loads(line)['global_step']==h['steps']+1)
            require(row['checkpoint_sha256']==first_sha,'first production successor checkpoint differs from isolated restart')
            save(status='handoff_completed_with_observed_step',observed_successor_step=nextp['accepted_steps'],observed_successor_time=nextp['accepted_time'],
                observed_first_checkpoint_sha256=row['checkpoint_sha256'],completed_UTC=datetime.now(timezone.utc).isoformat());return
        time.sleep(.2)
    raise RuntimeError('successor launched but no accepted step observed before explicit timeout')


def launch_validated_successor(a,old,h,inputs,first_sha,save):
    parent_lock=wait_parent_lock(a.campaign)
    current=load(a.campaign/'progress.json')
    require(current['accepted_steps']==h['steps'] and current['checkpoint_sha256']==h['sha256'],'parent advanced after validation')
    command=supervisor_command(a.solver,a.output,old['accepted_awf_case'],inputs,a.targets,a.timeout,
        old['accepted_awf_companion'],old['accepted_awf_companion_reason'])
    command[command.index('--block-steps')+1]=str(a.block_steps)
    command+=['--prefix-campaign',str(a.campaign),'--inner-refinement',old.get('accepted_awf_inner_refinement','none'),
        '--evaluation-reuse','ExactStepV1','--audit-revision','NativeElementSideV2']
    (a.output/'tmp').mkdir(exist_ok=True)
    env=dict(os.environ,**THREADS,LD_LIBRARY_PATH=str(a.solver.parent/'lib'),TMPDIR=str(a.output/'tmp'))
    with (a.output/'supervisor-launch.log').open('w') as log:
        child=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,pass_fds=(parent_lock.stream.fileno(),))
    new_identity=process_identity(child.pid);atomic(a.campaign/'runtime-handoff.json',dict(successor=str(a.output),
        successor_supervisor_identity=new_identity,parent_checkpoint_sha256=h['sha256'],receipt=str(a.output/'handoff.json')))
    parent_lock.close();save(status='successor_launched',successor_supervisor_identity=new_identity,command=command)
    observe_successor(a,h,first_sha,save)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--solver',type=Path,required=True);p.add_argument('--reference-solver',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--block-steps',type=int,choices=(8,32),default=8);p.add_argument('--timeout',type=float,default=7200)
    p.add_argument('--targets',type=int,nargs='+',default=[32,160,800,1600,3200,8000]);a=p.parse_args()
    require(sys.platform=='linux','Linux process and native checkpoint environment required')
    a.output.mkdir(parents=True,exist_ok=False);out=dict(status='preflight',started_UTC=datetime.now(timezone.utc).isoformat(),
        parent_campaign=str(a.campaign),old_audit_revision='LegacyInterfaceV1',new_audit_revision='NativeElementSideV2',
        evaluation_reuse='ExactStepV1',solver_sha256=sha(a.solver),probe_sha256=sha(a.probe),verification_governing_sweeps=0)
    def save(**changes):out.update(changes);atomic(a.output/'handoff.json',out)
    save()
    old=load(a.campaign/'progress.json');status=load(a.campaign.parent/'status.json')
    require(old['wall_policy']=='CLSCalibratedBudgetAWF' and old['case_family']=='AcceptedAWFV1','handoff is restricted to existing accepted-AWF campaigns')
    require(old['status'] not in ('numerically_rejected','targets_completed'),'terminal campaign cannot be handed off')
    require(not (a.campaign/'runtime-handoff.json').exists(),'parent already handed off')
    driver=process_identity(old['pid']);supervisor=status['supervisor_identity']
    require(matches(driver,dict(pid=old['pid'],start=old['process_start'])),'old driver is not live')
    require(matches(process_identity(supervisor['pid']),supervisor),'old supervisor identity differs')
    for identity,script in ((driver,'continue_cls_epawf_pimple.py'),(supervisor,'supervise_cls_epawf_pimple.py')):
        command=command_line(identity)
        require(script in command and str(a.campaign.parent) in command,'process command is not the selected campaign')
    save(old_driver_identity=driver,old_supervisor_identity=supervisor,old_executable_sha256=old['executable_sha256'],
        old_runtime_identity=load(a.campaign/'runtime/identity.json'),status='waiting_for_native_boundary')
    checked_signal(supervisor,signal.SIGSTOP);paused=False;terminated=False
    try:
        deadline=time.monotonic()+a.timeout
        while time.monotonic()<deadline:
            current=load(a.campaign/'progress.json')
            require(current['status']!='numerically_rejected','parent numerically rejected before handoff')
            if boundary(current):
                checked_signal(driver,signal.SIGSTOP)
                if not live_children(driver['pid']):paused=True;break
                checked_signal(driver,signal.SIGCONT)
            time.sleep(.02)
        require(paused,'no completed native boundary before finite handoff timeout')
        old=load(a.campaign/'progress.json');require(boundary(old) and not live_children(driver['pid']),'boundary changed')
        source=Path(old['segment'])/'case.evolution.json'
        if source.exists():payload=source.read_bytes()
        else:
            source=source.with_suffix('.json.gz')
            with gzip.open(source,'rb') as f:payload=f.read()
        native=json.loads(payload,parse_int=float);cp=Path(old['checkpoint']);h=header(cp,True)
        require(not native['failure'] and native['steps']==old['segment_target']==h['steps'] and native['time']==h['time'],'native boundary/report/checkpoint differs')
        require(state_hash(native['state'])==state_hash(h['coupled_state']) and
            [state_hash(v) for v in native['native_history']]==[state_hash(v) for v in h['history']] and
            state_hash(native['native_history_times'])==state_hash(h['history_times']),'native boundary state/history identity differs')
        budget=prefix_budget(a.campaign,h);inputs=a.output/'parent';inputs.mkdir();(inputs/'checkpoint').mkdir()
        shutil.copy2(cp/'State.bin',inputs/'checkpoint/State.bin');(inputs/'prefix.evolution.json').write_bytes(payload)
        shutil.copy2(a.campaign/'progress.json',inputs/'progress.json');shutil.copy2(a.campaign/'steps.jsonl',inputs/'steps.jsonl')
        # The immutable prefix copy identifies the actual parent CP path/hash.
        # Its progress pointer keeps original paths; prefix_budget uses hashes.
        save(status='revalidating_checkpoint',parent_checkpoint=str(cp),parent_checkpoint_sha256=h['sha256'],
            parent_step=h['steps'],parent_time=h['time'],parent_history_sha256=[state_hash(v) for v in h['history']],
            parent_history_times_sha256=state_hash(h['history_times']),prefix_correction_integral=budget,
            completed_targets_before_handoff=[n for n in a.targets if n<=h['steps']])
        q=checkpoint(inputs/'checkpoint');packet=dict(common_packet(inputs/'checkpoint'),operation='selected_original_operator',
            state=q['coupled'],history=q['coupled'],alpha=0.,time=q['time'],inner_stress_model=native['inner_stress_model'],
            witnesses=[[-3.,0.,3.6880433085817326]],reporting_breakpoints=[.0025,.005,.02,.04])
        diagnostic,receipt=run_probe(a.probe,Path(old['segment'])/'case.xml',packet,a.output/'checkpoint-reaudit')
        audit=diagnostic['complete_audits']['selected_state']
        require(audit['admissible'][0]==audit['accepted_awf_hard_admissible'][0]==1,'corrected audit rejected the parent checkpoint')
        fixture=a.output/'verification-input';fixture.mkdir();(fixture/'start').mkdir();shutil.copy2(cp/'State.bin',fixture/'start/State.bin')
        shutil.copy2(Path(old['segment'])/'case.xml',fixture/'production.xml')
        atomic(fixture/'capture.json',dict(segment=dict(start_step=h['steps'],end_step=h['steps']+1)))
        ref=replay(fixture,a.output/'verification-uncached',a.reference_solver,'Off',a.timeout)
        fast=replay(fixture,a.output/'verification-cached',a.solver,'ExactStepV1',a.timeout)
        proof=compare(ref,fast);save(status='validated_boundary',checkpoint_reaudit=receipt,complete_corrected_audit=audit,
            isolated_restart_equivalence=proof,verification_governing_sweeps=2,
            expected_first_successor_checkpoint_sha256=fast['checkpoint_sha256'])
        # Retire exactly the identified paused owners. No native child exists.
        for identity in (supervisor,driver):
            checked_signal(identity,signal.SIGTERM)
            try:os.kill(int(identity['pid']),signal.SIGCONT)
            except ProcessLookupError:pass
        terminated=True
        for _ in range(300):
            if not matches(process_identity(driver['pid']),driver) and not matches(process_identity(supervisor['pid']),supervisor):break
            time.sleep(.05)
        require(not matches(process_identity(driver['pid']),driver) and not matches(process_identity(supervisor['pid']),supervisor),'old owners did not stop')
        # A successor supervisor retains the old lock too, preventing a stale
        # old driver from resuming this lineage concurrently.
        launch_validated_successor(a,old,h,inputs,fast['checkpoint_sha256'],save);return

    except BaseException as error:
        save(status='handoff_failed',failure=repr(error),old_owners_terminated=terminated)
        if not terminated:
            # Infrastructure/setup failures leave the old physical run usable;
            # never resume after a demonstrated corrected hard rejection.
            hard='rejected' in str(error).lower()
            if not hard:
                if paused and matches(process_identity(driver['pid']),driver):checked_signal(driver,signal.SIGCONT)
                if matches(process_identity(supervisor['pid']),supervisor):checked_signal(supervisor,signal.SIGCONT)
            else:save(old_owners_paused_for_hard_failure=True)
        raise

if __name__=='__main__':main()
