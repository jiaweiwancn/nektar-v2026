#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable, fixed-dt continuation in checkpointed blocks; no numerical retries.

Run under nohup on Linux. The progress pointer is replaced atomically only
after a published native checkpoint exists. Every restart is validated by
the native reader. Detailed step packets are retained at selected contacts
and failures; other accepted packets/matrices are hashed then discarded.
"""
import argparse
import fcntl
import gzip
import json
import os
from pathlib import Path
import queue
import shutil
import socket
import struct
import subprocess
import threading
import time
import xml.etree.ElementTree as ET

from run_cls_epawf_picard import THREADS, load, require, sha, state_hash, write_json, set_parameter, set_info, write
from run_cls_epawf_pimple import prepare, identity, check_report, add_bound_work_arguments, bound_work_arguments


def atomic(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, path)


def load_progress(path,attempts=20):
    """Bounded reader retry for observed progress-replacement visibility gaps.

    A permanent absence/corruption remains an infrastructure error. Never
    interpret a failed read as permission to launch a replacement process.
    """
    for attempt in range(attempts):
        try:return load(path)
        except (FileNotFoundError,json.JSONDecodeError):
            if attempt==attempts-1:raise
            time.sleep(.05)


class CampaignLock:
    """The native child inherits the same flock: an orphan keeps ownership."""
    def __init__(self, output):
        self.stream=(output/'.campaign.lock').open('a+')
        try:fcntl.flock(self.stream,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            self.stream.close()
            raise RuntimeError('campaign lock held by a live driver or native child')
    def close(self):
        # Never LOCK_UN: a surviving child shares this open file description.
        self.stream.close()


def process_identity(pid):
    root=Path('/proc')/str(int(pid))
    try:
        fields=(root/'stat').read_text().rsplit(')',1)[1].split()
        if fields[0]=='Z':return None
        return dict(pid=int(pid),start=fields[19],executable=str((root/'exe').resolve(strict=True)),
                    cwd=str((root/'cwd').resolve(strict=True)))
    except (FileNotFoundError,ProcessLookupError):return None


def assert_native_stopped(progress):
    expected=progress.get('native_identity')
    if expected:
        actual=process_identity(expected['pid'])
        require(actual is None or actual!=expected,
                'recorded native solver is still alive; replacement/reconciliation refused')
    elif progress.get('native_pid'):
        # Old progress lacks a start token. Do not guess that a live PID is
        # safe to replace; retain it for explicit process/checkpoint inspection.
        require(process_identity(progress['native_pid']) is None,
                'legacy recorded native PID is alive; inspect identity before recovery')


class StepLedger:
    """One durable record per (global step, checkpoint SHA); idempotent replay."""
    def __init__(self,path):
        self.path=path;self.records={};lines=path.read_text().splitlines() if path.exists() else []
        repair=False
        for i,line in enumerate(lines):
            try:record=json.loads(line)
            except json.JSONDecodeError:
                require(i==len(lines)-1,'corrupt interior step record');repair=True;break
            step=record['global_step']
            if step in self.records:
                require(self.records[step]['checkpoint_sha256']==record['checkpoint_sha256'],
                        'conflicting checkpoint hashes for one accepted step')
                repair=True
            else:self.records[step]=record
        if repair:
            temporary=path.with_suffix('.reconcile.tmp')
            with temporary.open('w') as stream:
                for r in self.records.values():stream.write(json.dumps(r,separators=(',',':'))+'\n')
                stream.flush();os.fsync(stream.fileno())
            os.replace(temporary,path)
    def append(self,record):
        step=record['global_step']
        if step in self.records:
            require(self.records[step]['checkpoint_sha256']==record['checkpoint_sha256'],
                    'accepted-step identity conflict')
            return False
        with self.path.open('a') as stream:
            stream.write(json.dumps(record,separators=(',',':'))+'\n');stream.flush();os.fsync(stream.fileno())
        self.records[step]=record;return True
    def budget(self,progress,through):
        return progress.get('prefix_integral',0.)+sum(
            r['attempt'].get('k_bound',{}).get('integral_total',0.)
            for n,r in self.records.items() if progress.get('prefix_steps',0)<n<=through)


def header(checkpoint, full=False):
    b = (checkpoint/'State.bin').read_bytes(); n = struct.unpack_from('<Q', b, 24)[0]
    order, steps, t, dt = struct.unpack_from('<QQdd', b, 32+n)
    result=dict(order=order, steps=steps, time=t, dt=dt, sha256=sha(checkpoint/'State.bin'))
    if full:
        cursor=32+n+32+8  # signature, order/steps/time/dt, eta
        def vector():
            nonlocal cursor
            count=struct.unpack_from('<Q',b,cursor)[0];cursor+=8
            values=list(struct.unpack_from('<'+str(count)+'d',b,cursor));cursor+=8*count
            return values
        vector();slots=struct.unpack_from('<Q',b,cursor)[0];cursor+=8
        result['history']=[vector() for _ in range(slots)];result['history_times']=vector()
        result['coupled_state']=vector()
    return result


def prefix_budget(campaign, receipt):
    """Read an immutable prefix ledger; never repair or rewrite its source."""
    progress=load(campaign/'progress.json')
    require(progress['accepted_steps']==receipt['steps'] and progress['checkpoint_sha256']==receipt['sha256'],
            'prefix ledger/checkpoint identity differs')
    records={}
    for line in (campaign/'steps.jsonl').read_text().splitlines():
        row=json.loads(line);step=row['global_step']
        require(step not in records,'duplicate prefix ledger record')
        records[step]=row
    start=int(progress.get('prefix_steps',0));end=int(receipt['steps'])
    require(sorted(records)==list(range(start+1,end+1)),'prefix ledger is not contiguous')
    require(records[end]['checkpoint_sha256']==receipt['sha256'],'prefix last record differs')
    amount=progress.get('prefix_integral',0.)+sum(records[n]['attempt'].get('k_bound',{}).get('integral_total',0.) for n in sorted(records))
    require(abs(amount-progress['cumulative_bound_integral'])<1e-15,'prefix correction budget differs')
    return amount


def freeze(executable, output):
    output.mkdir(parents=True, exist_ok=False); lib = output/'lib'; lib.mkdir()
    before = identity(executable); shutil.copy2(executable, output/executable.name)
    # Preserve both the loader/SONAME path and its resolved basename.
    linked = subprocess.check_output(['ldd', str(executable)], text=True)
    for word in linked.split():
        if not word.startswith('/'): continue
        source = Path(word); names = {source.name, source.resolve().name}
        for name in names:
            target = lib/name
            if target.exists(): require(sha(target) == sha(source), 'library basename collision')
            else: shutil.copy2(source, target)
    os.environ['LD_LIBRARY_PATH'] = str(lib)
    after = identity(output/executable.name)
    require(after['sha256'] == before['sha256'], 'executable copy differs')
    require(set(after['libraries'].values()) == set(before['libraries'].values()), 'frozen library closure differs')
    require(all(Path(x).parent == lib for x in after['libraries'] if 'ld-linux' not in x), 'unfrozen linked library')
    atomic(output/'identity.json', dict(original=before, frozen=after))
    return output/executable.name, after


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--solver', type=Path, required=True)
    p.add_argument('--evaluation-reuse',choices=('Off','ExactStepV1'),default='Off')
    p.add_argument('--audit-revision',choices=('NativeElementSideV2',))
    p.add_argument('--checkpoint', type=Path)
    p.add_argument('--prefix-report', type=Path)
    p.add_argument('--prefix-campaign',type=Path)
    p.add_argument('--from-h0',action='store_true')
    p.add_argument('--wall-policy',choices=('Upstream','UniformEPFlat','CLSCalibratedBudgetAWF'),default='Upstream')
    p.add_argument('--height-case',choices=list('ABCDEFGH'))
    p.add_argument('--upstream-case',choices=('U05-graded','U05-one','U20-graded','U20-one'))
    p.add_argument('--accepted-awf-case',choices=('Accepted-AWF-2G','Accepted-AWF-05G','Accepted-AWF-05O','Accepted-AWF-2O'))
    p.add_argument('--accepted-awf-companion',choices=('faithful','regularized','saturated','regularized-saturated'),default='faithful')
    p.add_argument('--accepted-awf-companion-reason',default='')
    p.add_argument('--accepted-awf-perturbation',type=float,default=0.)
    p.add_argument('--inner-refinement','--accepted-awf-inner-refinement',dest='accepted_awf_inner_refinement',choices=('none','cutoff-v1'),default='none')
    p.add_argument('--pressure-space',choices=('baseline-v1','inner-p2-v1'),default='baseline-v1')
    p.add_argument('--inventory',choices=('LegacyV1','ConstitutiveV2','ConstitutiveV3','AcceptedBudgetV1'),default='LegacyV1')
    p.add_argument('--bound-solver',choices=('WorkingSet','RetainedPool'),default='WorkingSet')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--targets', type=int, nargs='+', default=[32,160,800,1600,3200,8000])
    p.add_argument('--block-steps', type=int, default=16)
    p.add_argument('--timeout', type=float, default=7200)
    p.add_argument('--resume', action='store_true')
    add_bound_work_arguments(p)
    a = p.parse_args(); a.output = a.output.resolve()
    a.checkpoint = a.checkpoint.resolve() if a.checkpoint else None
    require(1 <= a.block_steps <= 64 and a.timeout > 0, 'invalid block controls')
    require(sum(bool(x) for x in (a.height_case,a.upstream_case,a.accepted_awf_case))<=1,'choose one declared comparison family')
    if a.upstream_case:require(a.wall_policy=='Upstream' and a.inventory=='ConstitutiveV3','Upstream family policy mismatch')
    if a.accepted_awf_case:
        from run_cls_accepted_awf import WORK,check_report as accepted_check_report
        require(a.wall_policy=='CLSCalibratedBudgetAWF' and a.inventory=='AcceptedBudgetV1' and
            a.bound_solver=='RetainedPool' and bound_work_arguments(a)==WORK,'accepted-AWF family/method controls mismatch')
        require(a.accepted_awf_companion=='faithful' or a.accepted_awf_companion_reason.strip(),'companion needs its preserved faithful blocker')
        require(a.accepted_awf_inner_refinement=='none' or a.accepted_awf_case.endswith('G'),'focused cutoff refinement requires a graded accepted-AWF case')
        require(a.pressure_space=='baseline-v1' or (a.accepted_awf_case=='Accepted-AWF-05G' and
            a.accepted_awf_inner_refinement=='cutoff-v1'),'inner-p2-v1 requires the matched cutoff-v1 H=.5 graded case')
    else:
        require(a.wall_policy!='CLSCalibratedBudgetAWF' and a.inventory!='AcceptedBudgetV1','accepted model requires its declared case family')
        require(a.accepted_awf_companion=='faithful' and a.accepted_awf_perturbation==0. and a.accepted_awf_inner_refinement=='none','accepted model controls on another case family')
        require(a.pressure_space=='baseline-v1','pressure-space control requires its declared accepted-AWF case')
    os.environ.update(THREADS); os.environ['TMPDIR'] = str(a.output/'tmp')
    if not a.resume:a.output.mkdir(parents=True,exist_ok=False)
    require(not (a.output/'runtime-handoff.json').exists(),'campaign sealed by a validated runtime handoff; use its recorded successor')
    lock=CampaignLock(a.output)
    ledger=StepLedger(a.output/'steps.jsonl')
    if a.resume:
        progress = load(a.output/'progress.json'); a.checkpoint = Path(progress['checkpoint']) if progress['checkpoint'] else None
        require(progress['status']!='numerically_rejected','numerical rejection is terminal; use a separately identified follow-up')
        require(progress.get('evaluation_reuse','Off')==a.evaluation_reuse,'evaluation reuse differs from campaign')
        require(progress.get('audit_revision')==a.audit_revision,'audit revision differs from campaign')
        require(progress.get('inventory','LegacyV1')==a.inventory,'inventory differs from campaign')
        require(progress.get('bound_solver','WorkingSet')==a.bound_solver,'finite-QP backend differs from campaign')
        require(progress.get('wall_policy','Upstream')==a.wall_policy,'wall policy differs from campaign')
        require(progress.get('height_case')==a.height_case,'height case differs from campaign')
        require(progress.get('upstream_case')==a.upstream_case,'Upstream case differs from campaign')
        require(progress.get('accepted_awf_case')==a.accepted_awf_case,'accepted-AWF case differs from campaign')
        require(progress.get('accepted_awf_companion','faithful')==a.accepted_awf_companion,'accepted-AWF companion differs from campaign')
        require(progress.get('accepted_awf_perturbation',0.)==a.accepted_awf_perturbation,'accepted-AWF initialization differs from campaign')
        require(progress.get('accepted_awf_inner_refinement','none')==a.accepted_awf_inner_refinement,'accepted-AWF inner mesh policy differs from campaign')
        require(progress.get('accepted_awf_pressure_space','baseline-v1')==a.pressure_space,'accepted-AWF pressure space differs from campaign')
        if a.accepted_awf_case:require(progress.get('accepted_awf_companion_reason','')==a.accepted_awf_companion_reason,'companion provenance differs from campaign')
        require(progress.get('projection_work_controls',{})==bound_work_arguments(a),'projection work controls differ from recorded campaign')
        assert_native_stopped(progress)
        old_stat=Path('/proc')/str(int(progress['pid']))/'stat'
        if old_stat.exists():
            stat=old_stat.read_text().split()
            require(stat[21]!=progress.get('process_start') or stat[2]=='Z','campaign is already running')
        if progress.get('failure'):
            progress.setdefault('interruption_history',[]).append(dict(reason=progress.pop('failure'),time=progress.get('updated_unix')))
        executable = a.output/'runtime'/a.solver.name
        os.environ['LD_LIBRARY_PATH'] = str(a.output/'runtime/lib')
        pinned = load(a.output/'runtime/identity.json')['frozen']
        # Reconcile only genuine native publications after an orchestration
        # interruption. Never manufacture or edit a checkpoint payload.
        segment=Path(progress.get('segment','/nonexistent'))
        candidates=sorted((x.parent for x in segment.glob('checkpoint*/State.bin')),key=lambda x:header(x)['steps'])
        for checkpoint in candidates:
            receipt=header(checkpoint,True);step=receipt['steps']
            if step<=progress['accepted_steps']:continue
            require(step==progress['accepted_steps']+1,'noncontiguous checkpoint recovery')
            packets=[x for x in segment.glob('pimple-states/*.json') if abs(load(x)['time']-receipt['time'])<1e-13]
            require(len(packets)==1,'missing recovery packet');packet=load(packets[0])
            require(state_hash(receipt['coupled_state'])==state_hash(packet['corrected']),'publication/candidate mismatch')
            records=[json.loads(x[len('EP-AWF PIMPLE '):]) for x in (segment/'solver.log').read_text().splitlines() if x.startswith('EP-AWF PIMPLE ')]
            attempt=next(x for x in records if abs(x['time']-receipt['time'])<1e-13)
            record=dict(global_step=step,time=receipt['time'],checkpoint=str(checkpoint),
                checkpoint_sha256=receipt['sha256'],attempt=attempt,recovered_publication=True)
            ledger.append(record)
            progress.update(accepted_steps=step,accepted_time=receipt['time'],checkpoint=str(checkpoint),
                checkpoint_sha256=receipt['sha256'],cumulative_bound_integral=ledger.budget(progress,step))
        progress['cumulative_bound_integral']=ledger.budget(progress,progress['accepted_steps'])
        # A driver can die after its child records a numerical failure but
        # before the terminal progress replace. That is still not a retry.
        orphan_report=segment/'case.evolution.json'
        if orphan_report.exists():
            report=load(orphan_report)
            if report.get('failure'):
                progress.update(status='numerically_rejected',failure=report['failure'],updated_unix=time.time())
                atomic(a.output/'progress.json',progress)
                raise ValueError('native report contains terminal numerical rejection; no replacement launched')

    else:
        (a.output/'tmp').mkdir()
        executable, pinned = freeze(a.solver.resolve(), a.output/'runtime')
        if a.from_h0:
            require(a.checkpoint is None and a.prefix_report is None,'h0 is not checkpoint migration')
            progress=dict(schema=2,status='starting',checkpoint=None,accepted_steps=0,accepted_time=0.,dt=.000625,
                cumulative_bound_integral=0.,prefix_steps=0,prefix_integral=0.,inventory=a.inventory,wall_policy=a.wall_policy,
                requested_targets=a.targets,completed_targets=[],segments=[])
        else:
            require(a.checkpoint and a.prefix_report,'restart needs real checkpoint and prefix report')
            prefix = load(a.prefix_report); info = header(a.checkpoint)
            require(info['order'] == 1 and info['dt'] == .000625 and prefix['steps'] == info['steps'], 'incompatible primary restart')
            expected_policy='CLSCalibratedBudgetAWF' if a.accepted_awf_case else ('UniformEPFlat-v1' if a.wall_policy=='UniformEPFlat' else 'Upstream')
            require(prefix.get('wall_policy','Upstream')==expected_policy,'wrong prefix wall policy')
            if a.accepted_awf_case:
                from run_cls_accepted_awf import check_report as accepted_check_report,model_choices
                accepted_check_report(prefix,True)
                require(prefix['model_variant']==model_choices(a.accepted_awf_companion)['outer_variant'] and
                    prefix['inner_stress_model']==model_choices(a.accepted_awf_companion)['inner_stress'],'wrong prefix companion')
            require(prefix['k_bound_mode'] == 'MassProjection', 'wrong prefix mode')
            integral = sum(v.get('k_bound',{}).get('integral_total',0.) for v in prefix['pimple_attempts'] if v['time']<=info['time']+1e-13)
            if a.prefix_campaign:integral=prefix_budget(a.prefix_campaign.resolve(),info)
            progress = dict(schema=2, status='starting', checkpoint=str(a.checkpoint), accepted_steps=info['steps'],
                accepted_time=info['time'], dt=info['dt'], cumulative_bound_integral=integral,inventory=a.inventory,wall_policy=a.wall_policy,
                prefix_report=str(a.prefix_report.resolve()), prefix_report_sha256=sha(a.prefix_report),
                prefix_campaign=str(a.prefix_campaign.resolve()) if a.prefix_campaign else None,
                prefix_steps=info['steps'], prefix_integral=integral, requested_targets=a.targets, completed_targets=[], segments=[])
    progress.update(pid=os.getpid(), host=socket.gethostname(), started_unix=time.time(),
                    evaluation_reuse=a.evaluation_reuse,audit_revision=a.audit_revision,
                    projection_work_controls=bound_work_arguments(a),bound_solver=a.bound_solver,height_case=a.height_case,upstream_case=a.upstream_case,
                    accepted_awf_case=a.accepted_awf_case,accepted_awf_companion=a.accepted_awf_companion,
                    accepted_awf_companion_reason=a.accepted_awf_companion_reason,accepted_awf_perturbation=a.accepted_awf_perturbation,
                    accepted_awf_inner_refinement=a.accepted_awf_inner_refinement,
                    accepted_awf_pressure_space=a.pressure_space,
                    case_family='AcceptedAWFV1' if a.accepted_awf_case else ('UpstreamHeightV1' if a.upstream_case else ('HeightSamplerV1' if a.height_case else 'LegacyOuter6')),
                    process_start=Path('/proc/self/stat').read_text().split()[21], executable=str(executable),
                    executable_sha256=pinned['sha256'], source_revision=subprocess.check_output(['git','-C',str(Path(__file__).resolve().parents[4]),'rev-parse','HEAD'],text=True).strip())
    def save(**changes):
        progress.update(changes); progress['updated_unix'] = time.time(); atomic(a.output/'progress.json', progress)
    def verify_identity():
        current = identity(executable)
        require(current['sha256'] == pinned['sha256'] and current['libraries'] == pinned['libraries'], 'runtime identity changed')
    save()
    try:
        while progress['accepted_steps'] < max(a.targets):
            if (a.output/'stop-after-block.json').exists():
                save(status='stopped_at_checkpoint',stop_request_sha256=sha(a.output/'stop-after-block.json'));return
            verify_identity()
            workers = subprocess.check_output(['ps','-eo','comm=,nlwp='],text=True)
            count = sum(int(x.split()[-1]) for x in workers.splitlines() if x.split()[0] in ('IncNavierStokes','cc1plus','RANSWallResolut'))
            require(count+3 <= 72, 'worker budget unavailable; no native process launched')
            start = int(progress['accepted_steps']); stage = next(v for v in a.targets if v > start)
            end = min(start+a.block_steps, stage)
            directory = a.output/f'steps-{start+1:06d}-{end:06d}'
            if directory.exists():
                # An interrupted attempt is immutable; retry from the committed
                # native checkpoint in a distinct attempt directory.
                suffix=1
                while Path(str(directory)+f'-attempt{suffix}').exists(): suffix+=1
                directory=Path(str(directory)+f'-attempt{suffix}')
            if a.accepted_awf_case:
                from run_cls_accepted_awf import case_xml
                xml=case_xml(a.accepted_awf_case,directory,end-start,.000625,
                    restart=Path(progress['checkpoint']) if progress['checkpoint'] else None,
                    companion=a.accepted_awf_companion,companion_reason=a.accepted_awf_companion_reason,
                    perturbation=a.accepted_awf_perturbation,inner_refinement=a.accepted_awf_inner_refinement,
                    pressure_space=a.pressure_space)
            elif a.height_case or a.upstream_case:
                from run_cls_epawf_height_comparison import case_xml,WORK
                require((a.upstream_case and a.wall_policy=='Upstream' and a.inventory=='ConstitutiveV3') or
                    (a.height_case and a.wall_policy=='UniformEPFlat' and a.inventory in ('ConstitutiveV2','ConstitutiveV3')),
                    'height continuation requires its declared policy')
                require(not bound_work_arguments(a) or bound_work_arguments(a)==WORK,'height continuation uses its validated explicit controls')
                if a.upstream_case:require(bound_work_arguments(a)==WORK,'new Upstream continuation requires explicit validated projection controls')
                xml=case_xml(a.upstream_case or a.height_case,directory,end-start,.000625,
                    restart=Path(progress['checkpoint']) if progress['checkpoint'] else None,inventory=a.inventory,wall_policy=a.wall_policy,bound_solver=a.bound_solver)
            else:
                xml = prepare(directory, .000625, end-start, True, restart=Path(progress['checkpoint']) if progress['checkpoint'] else None,
                              checkpoint=True, k_bound='MassProjection', floor_relative=1e-8, snapshots=True,inventory=a.inventory,wall_policy=a.wall_policy,
                              bound_work=bound_work_arguments(a),bound_solver=a.bound_solver)
            tree=ET.parse(xml).getroot();set_info(tree,'RANSWallEvaluationReuse',a.evaluation_reuse)
            set_parameter(tree,'IO_CheckSteps',1);write(tree,xml)
            if a.height_case or a.upstream_case or a.accepted_awf_case:
                progress['case_definition']=load(directory/'definition.json')
            command=[str(executable),'case.xml']; started=time.monotonic()
            child=subprocess.Popen(command,cwd=directory,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                   text=True,bufsize=1,env=os.environ.copy(),pass_fds=(lock.stream.fileno(),))
            save(status='running',native_pid=child.pid,native_identity=process_identity(child.pid),segment=str(directory),segment_target=end,next_target=stage,
                 segment_start_steps=start,command=command,xml_sha256=sha(xml),existing_worker_threads=count)
            lines=queue.Queue()
            def reader():
                for line in child.stdout: lines.put(line)
                lines.put(None)
            threading.Thread(target=reader,daemon=True).start()
            pending=None; local=[]; infra=None
            with (directory/'solver.log').open('w') as log:
                while True:
                    if time.monotonic()-started>a.timeout:
                        infra='block_timeout';child.terminate();child.wait(timeout=30);break
                    try: line=lines.get(timeout=5)
                    except queue.Empty:
                        save();continue
                    if line is None:
                        child.wait();break
                    for line in (line,):
                        log.write(line);log.flush()
                        if line.startswith('EP-AWF PIMPLE '): pending=json.loads(line[len('EP-AWF PIMPLE '):])
                        if line.startswith('EP-AWF IMEX1 step '):
                            fields=line.split(); step=int(fields[3])
                            require(pending is not None and pending['scheduled_complete'], 'missing completed schedule')
                            deadline=time.monotonic()+30
                            candidates=[]
                            while not candidates and time.monotonic()<deadline:
                                candidates=[x.parent for x in directory.glob('checkpoint*/State.bin') if header(x.parent)['steps']==step]
                                if not candidates:time.sleep(.05)
                            require(candidates, 'accepted step has no published checkpoint')
                            checkpoint=max(candidates,key=lambda x:(x/'State.bin').stat().st_mtime_ns)
                            receipt=header(checkpoint); require(receipt['dt']==.000625, 'checkpoint dt changed')
                            bound=pending.get('k_bound',{}); amount=bound.get('integral_total',0.)
                            record=dict(global_step=step,time=receipt['time'],segment=str(directory),checkpoint=str(checkpoint),
                                checkpoint_sha256=receipt['sha256'],attempt=pending)
                            ledger.append(record)
                            local.append(record)
                            save(accepted_steps=step,accepted_time=receipt['time'],checkpoint=str(checkpoint),
                                 checkpoint_sha256=receipt['sha256'],cumulative_bound_integral=ledger.budget(progress,step))
                            # Keep early contact/failure evidence and every
                            # 32nd step. Delete only *accepted* unselected packets.
                            packet=Path(pending['k_snapshot']) if pending.get('k_snapshot') else None
                            if packet and step>32 and step%32 and step!=start+1 and step not in a.targets:
                                saved=load(packet)
                                with (a.output/'discarded-packets.jsonl').open('a') as tomb:
                                    tomb.write(json.dumps(dict(step=step,path=str(packet),sha256=sha(packet),matrix_sha256=saved['matrix_sha256']))+'\n')
                                Path(saved['matrix_file']).unlink();packet.unlink()
                            # Keep two most recent checkpoints and fixed-time
                            # samples. The pointer is already atomically safe.
                            for old in directory.glob('checkpoint*/State.bin'):
                                n=header(old.parent)['steps']
                                if n<step-1 and n%32 and n not in a.targets: shutil.rmtree(old.parent)
                            pending=None
                    save()
            verify_identity()
            report_file=directory/'case.evolution.json'
            if report_file.exists():
                report=load(report_file)
                if a.audit_revision:require(report.get('audit_revision')==a.audit_revision,'native audit revision differs')
                require(report.get('evaluation_reuse',{}).get('mode','Off')==a.evaluation_reuse,'native evaluation reuse differs')
                if report.get('initialization_complete') is False:
                    save(status='numerically_rejected',failure=report['failure'],exit_code=child.returncode);return
                if a.accepted_awf_case:accepted_check_report(report,start>0,progress['case_definition'])
                else:check_report(report,start>0)
                require(report['initial_k_bound_check'].get('work_controls',{}).get('finite_qp_solver','WorkingSet')==a.bound_solver,
                    'native runtime did not use the recorded finite-QP backend')
                require(report['steps']==progress['accepted_steps'], 'report/checkpoint count differs')
                for record,audit in zip(local,report['state_audits'][1:]):record['audit']=audit
                with gzip.open(directory/'accepted-metrics.json.gz','wt') as f: json.dump(local,f,separators=(',',':'))
                summary=dict(start_step=start,end_step=int(report['steps']),time=report['time'],requested_end_step=end,
                    seconds=time.monotonic()-started,exit_code=child.returncode,failure=report['failure'],
                    report_sha256=sha(report_file),checkpoint=progress['checkpoint'])
                atomic(directory/'receipt.json',summary); progress['segments'].append(summary)
                if not infra and child.returncode==0:
                    require(report['steps']==end,'successful block missed integer target')
                    if end in a.targets: progress['completed_targets'].append(end)
                    # Large accepted-state/audit arrays are compressed per
                    # bounded-size segment; all hard checks ran before commit.
                    with report_file.open('rb') as src,gzip.open(str(report_file)+'.gz','wb') as dst: shutil.copyfileobj(src,dst)
                    report_file.unlink()
                    save(status='between_blocks');continue
                save(status='infrastructure_interrupted' if infra else 'numerically_rejected',failure=infra or report['failure']);return
            save(status='infrastructure_interrupted',failure=infra or 'process_exit_without_native_report',exit_code=child.returncode);return
        save(status='targets_completed',next_target=None)
    except BaseException as error:
        if 'child' in locals() and child.poll() is None:
            child.terminate();child.wait(timeout=30)
        save(status='orchestration_interrupted',failure=repr(error));raise


if __name__=='__main__': main()
