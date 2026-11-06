#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Durable, fixed-dt continuation in checkpointed blocks; no numerical retries.

Run under nohup on Linux. The progress pointer is replaced atomically only
after a published native checkpoint exists. Every restart is validated by
the native reader. Detailed step packets are retained at selected contacts
and failures; other accepted packets/matrices are hashed then discarded.
"""
import argparse
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

from run_cls_epawf_picard import THREADS, load, require, sha, state_hash, write_json, set_parameter, write
from run_cls_epawf_pimple import prepare, identity, check_report


def atomic(path, value):
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as f:
        json.dump(value, f, indent=2, allow_nan=False); f.write('\n'); f.flush(); os.fsync(f.fileno())
    os.replace(temporary, path)


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
    p.add_argument('--checkpoint', type=Path, required=True)
    p.add_argument('--prefix-report', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--targets', type=int, nargs='+', default=[32,160,800,1600,3200,8000])
    p.add_argument('--block-steps', type=int, default=16)
    p.add_argument('--timeout', type=float, default=7200)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args(); a.output = a.output.resolve(); a.checkpoint = a.checkpoint.resolve()
    require(1 <= a.block_steps <= 64 and a.timeout > 0, 'invalid block controls')
    os.environ.update(THREADS); os.environ['TMPDIR'] = str(a.output/'tmp')
    if a.resume:
        progress = load(a.output/'progress.json'); a.checkpoint = Path(progress['checkpoint'])
        require(progress['status']!='numerically_rejected','numerical rejection is terminal; use a separately identified follow-up')
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
            with (a.output/'steps.jsonl').open('a') as f:
                f.write(json.dumps(dict(global_step=step,time=receipt['time'],checkpoint=str(checkpoint),
                    checkpoint_sha256=receipt['sha256'],attempt=attempt,recovered_publication=True))+'\n')
            progress.update(accepted_steps=step,accepted_time=receipt['time'],checkpoint=str(checkpoint),
                checkpoint_sha256=receipt['sha256'],cumulative_bound_integral=progress['cumulative_bound_integral']+packet['k_bound']['integral_total'])
    else:
        a.output.mkdir(parents=True, exist_ok=False); (a.output/'tmp').mkdir()
        executable, pinned = freeze(a.solver.resolve(), a.output/'runtime')
        prefix = load(a.prefix_report); info = header(a.checkpoint)
        require(info['order'] == 1 and info['dt'] == .000625 and prefix['steps'] == info['steps'], 'incompatible primary restart')
        require(prefix['k_bound_mode'] == 'MassProjection', 'wrong prefix mode')
        integral = sum(v.get('k_bound',{}).get('integral_total',0.) for v in prefix['pimple_attempts'])
        progress = dict(schema=1, status='starting', checkpoint=str(a.checkpoint), accepted_steps=info['steps'],
            accepted_time=info['time'], dt=info['dt'], cumulative_bound_integral=integral,
            prefix_report=str(a.prefix_report.resolve()), prefix_report_sha256=sha(a.prefix_report),
            prefix_steps=info['steps'], prefix_integral=integral, requested_targets=a.targets, completed_targets=[], segments=[])
    progress.update(pid=os.getpid(), host=socket.gethostname(), started_unix=time.time(),
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
            xml = prepare(directory, .000625, end-start, True, restart=Path(progress['checkpoint']),
                          checkpoint=True, k_bound='MassProjection', floor_relative=1e-8, snapshots=True)
            tree=ET.parse(xml).getroot(); set_parameter(tree,'IO_CheckSteps',1); write(tree,xml)
            command=[str(executable),'case.xml']; started=time.monotonic()
            child=subprocess.Popen(command,cwd=directory,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                                   text=True,bufsize=1,env=os.environ.copy())
            save(status='running',native_pid=child.pid,segment=str(directory),segment_target=end,next_target=stage,
                 segment_start_steps=start,command=command,xml_sha256=sha(xml),existing_worker_threads=count)
            lines=queue.Queue()
            def reader():
                for line in child.stdout: lines.put(line)
                lines.put(None)
            threading.Thread(target=reader,daemon=True).start()
            pending=None; local=[]; infra=None
            with (directory/'solver.log').open('w') as log, (a.output/'steps.jsonl').open('a') as metrics:
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
                            metrics.write(json.dumps(record,separators=(',',':'))+'\n');metrics.flush();os.fsync(metrics.fileno())
                            local.append(record)
                            save(accepted_steps=step,accepted_time=receipt['time'],checkpoint=str(checkpoint),
                                 checkpoint_sha256=receipt['sha256'],cumulative_bound_integral=progress['cumulative_bound_integral']+amount)
                            # Keep early contact/failure evidence and every
                            # 32nd step. Delete only *accepted* unselected packets.
                            packet=Path(pending['k_snapshot']) if pending.get('k_snapshot') else None
                            if packet and step>32 and step%32:
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
                report=load(report_file);check_report(report,True)
                require(report['steps']==progress['accepted_steps'], 'report/checkpoint count differs')
                for i,audit in enumerate(report['state_audits'][1:]):
                    local[i]['audit']=audit
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
