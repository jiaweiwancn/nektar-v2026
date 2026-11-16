#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Separate H25R1 finite account. Historical H25 counters are never enlarged."""
import argparse
import datetime
import json
import os
from pathlib import Path
import signal
import subprocess
import time
from run_wj_h4 import THREADS, load, save, sha, input_files
from run_wj_h7 import update_journal, complete_operation
from run_wj_h25 import usage
from wj6_resources import scientific_workers, tree_usage
from wj7_command import owned_groups
from run_wj_h15 import size

os.environ.update(THREADS)


def initialize(parent, revision):
    root = parent/'resumption-r1'
    root.mkdir(exist_ok=False)
    for n in ('science', 'operations', 'tmp', 'launches'): (root/n).mkdir()
    save(root/'accounting.json', dict(schema='WJH25R1FiniteContinuationV1', closed=False,
        ceilings=dict(total=21600, native=14400, build=2400, passive=3600),
        setup_seconds=300, closure_seconds=900, entries={}, reserved_attempts=0,
        maximum_attempts=384, governor_blocks=0, maximum_governor_blocks=24,
        source_revision=revision, original_account=dict(path=str(parent/'accounting.json'),
            sha256=sha(parent/'accounting.json'), usage=usage(load(parent/'accounting.json'))),
        scope='prospective WJ-H25R1 allocation; original charges/ceilings preserved'))


def admission(parent, phase):
    allocation = scientific_workers()
    extra = 16 if phase == 'build' else 10 if phase == 'native' else 2
    free = os.statvfs(parent).f_bavail * os.statvfs(parent).f_frsize
    return allocation, extra, free, allocation['workers'] + extra <= 72 and free >= 15*1024**3


def execute(parent, label, argv, seconds, phase, attempts=0, final=False):
    root = parent/'resumption-r1'
    assert parent.name == 'wj-h25-production-20261116' and root.is_dir()
    assert label.replace('-','').isalnum() and phase in ('native','build','passive')
    assert 0 < seconds <= dict(native=2400,build=900,passive=300)[phase]
    assert 0 <= attempts <= 128 and (not attempts or phase == 'native')
    allocation, extra, free, okay = admission(parent, phase)
    folder = root/'operations'/label
    assert not folder.exists(), 'immutable operation label already used'
    folder.mkdir()
    request = dict(argv=list(map(str,argv)), phase=phase, reserved_seconds=seconds,
        attempts=attempts, threads=THREADS, allocation=allocation,
        additional_workers_reserved=extra, free_bytes=free, final_operation=final,
        inputs=input_files(argv), runner_sha256=sha(__file__),
        utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    save(folder/'request.json', request)
    if not okay:
        event = dict(status='RESOURCE_NOT_ADMITTED',phase='passive',seconds=1.,
            reserved_seconds=1.,native_launches=0, reason='host workers or free-space reserve',
            reserved_attempts=0, actual_attempts=0, known=True)
        def refuse(b):
            assert not b['closed'] and label not in b['entries']
            assert usage(b)['passive']+1 <= b['ceilings']['passive']
            b['entries'][label] = event
        update_journal(root/'accounting.json',refuse)
        save(folder/'receipt.json',dict(request,**event))
        print(label,event['status'],allocation['workers'],'workers',flush=True)
        return False

    def reserve(b):
        assert not b['closed'] and label not in b['entries']
        u = usage(b)
        assert u['total']+seconds <= b['ceilings']['total']
        assert u[phase]+seconds+(900 if phase=='native' and not final else 0) <= b['ceilings'][phase]
        if attempts:
            assert b['reserved_attempts']+attempts <= b['maximum_attempts']
            assert b['governor_blocks'] < b['maximum_governor_blocks']
            b['reserved_attempts'] += attempts
            b['governor_blocks'] += 1
        b['entries'][label] = dict(status='RESERVED',phase=phase,reserved_seconds=seconds,
            reserved_attempts=attempts)
    update_journal(root/'accounting.json',reserve)
    env = dict(os.environ,**THREADS,PYTHONPATH=os.pathsep.join((str(Path(__file__).resolve().parent),
        str(parent/'source-v8/solvers/IncNavierStokesSolver/Tests/RANS'))),
        TMPDIR=str(root/'tmp'),MPLCONFIGDIR=str(root/'tmp/matplotlib'),
        NEKTAR_RANS_TIMING_FILE=str(folder/'timing.jsonl'),WJ_H25_OPERATION=str(folder),
        WJ_H25_DEADLINE=str(time.monotonic()+seconds),
        WJ_H14_PARENT='/home/jay/build/wj-h13-fullheight-paired-domain-20261113/publication-v1/paired-case')
    start=time.monotonic(); peak=0; groups=set(); reason=None; maps=[]; seen=set(); last_check=-100.

    def native_maps(pid):
        nonlocal reason
        try:
            paths=sorted({str(Path(s.split(maxsplit=5)[5]).resolve())
                for s in Path(f'/proc/{pid}/maps').read_text().splitlines()
                if len(s.split(maxsplit=5))==6 and s.split(maxsplit=5)[5].startswith('/') and not s.endswith(' (deleted)')})
            if any(p.endswith('.elf') for p in paths) and any(p.endswith('/libc.so.6') for p in paths):
                key=(pid,tuple(paths))
                if key not in seen:
                    assert len(maps)<32, 'bounded native map catalog'
                    seen.add(key)
                    okay=[p for p in paths if p.endswith('/libc.so.6')]==['/home/jay/build/wj-h22-early-phase-replay-20261115/runtime/glibc/libc.so.6']
                    maps.append(dict(pid=pid,status='PASS' if okay else 'FAIL',files=paths,elapsed=time.monotonic()-start))
                    if not okay: reason='MIXED_NATIVE_RUNTIME'
            for p in Path(f'/proc/{pid}/task/{pid}/children').read_text().split(): native_maps(int(p))
        except (FileNotFoundError,ProcessLookupError): pass

    with (folder/'command.log').open('x') as log:
        child=subprocess.Popen(['/usr/bin/time','-v',*map(str,argv)],cwd=parent,env=env,
            stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while child.poll() is None:
            try:
                rss,_=tree_usage(os.getpid());peak=max(peak,rss);groups |= owned_groups(child.pid);native_maps(child.pid)
                if rss>(16 if phase=='build' else 8)*1024**3: reason='RSS_LIMIT'
                if time.monotonic()-start>seconds: reason='TIME_LIMIT'
                if time.monotonic()-start-last_check>3:
                    last_check=time.monotonic()-start
                    if scientific_workers()['workers']>72: reason='HOST_WORKER_LIMIT'
                    if os.statvfs(root).f_bavail*os.statvfs(root).f_frsize<15*1024**3: reason='FREE_SPACE_LIMIT'
                    if size(root/'science')>8*1024**3: reason='OUTPUT_LIMIT'
                    code_size=sum(size(p) for p in root.iterdir() if p.is_dir() and p.name.startswith(('source','build','runtime')))
                    if code_size>8*1024**3: reason='CODE_STORAGE_LIMIT'
            except Exception as e: reason='MONITOR_ERROR '+repr(e)
            if reason:
                for g in groups|{child.pid}:
                    try: os.killpg(g,signal.SIGTERM)
                    except ProcessLookupError: pass
                break
            time.sleep(.25)
        try: code=child.wait(timeout=5)
        except subprocess.TimeoutExpired:
            for g in groups|{child.pid}:
                try: os.killpg(g,signal.SIGKILL)
                except ProcessLookupError: pass
            code=child.wait()
    event=dict(status='PASS' if code==0 and reason is None else 'FAIL',phase=phase,
        reserved_seconds=seconds,seconds=time.monotonic()-start,exit=code,reason=reason,
        peak_tree_rss_bytes=peak,reserved_attempts=attempts,unknown_attempts=bool(attempts),
        log_sha256=sha(folder/'command.log'))
    if (folder/'attempt-counts.json').exists():
        actual=load(folder/'attempt-counts.json')
        assert type(actual['actual_attempts']) is int and 0<=actual['actual_attempts']<=attempts
        event.update(actual);event['unknown_attempts']=not actual.get('known',True)
    save(folder/'native-exec-maps.json',maps)
    complete_operation(root/'accounting.json',label,event)
    save(folder/'receipt.json',dict(request,**event))
    print(label,event['status'],event['seconds'],flush=True)
    return code==0 and reason is None


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    p.add_argument('--initialize',action='store_true');p.add_argument('--revision')
    p.add_argument('--label');p.add_argument('--seconds',type=float);p.add_argument('--phase',choices=('build','native','passive'))
    p.add_argument('--attempts',type=int,default=0);p.add_argument('--final-operation',action='store_true')
    p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    if a.initialize: initialize(a.root,a.revision)
    else:
        cmd=a.command[1:] if a.command[:1]==['--'] else a.command
        raise SystemExit(0 if execute(a.root,a.label,cmd,a.seconds,a.phase,a.attempts,a.final_operation) else 1)
