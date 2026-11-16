#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Synthetic paused-family fixture; never a scientific validation result."""
import json,os,signal,sys
from pathlib import Path


def parent(pid):return int((Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()[1])


control=Path(sys.argv[1]).resolve();owned=parent(os.getpid());outer=parent(owned)
commands={}
for pid in (Path('/proc')/str(outer)/'task'/str(outer)/'children').read_text().split():
    try:commands[int(pid)]=(Path('/proc')/pid/'cmdline').read_bytes().decode().split('\0')[:-1]
    except (FileNotFoundError,ProcessLookupError):pass
# Match an argv token, never a substring in another process's Python -c body.
guards=[pid for pid,cmd in commands.items() if any(Path(word).name=='guard_rans_worker_budget.py' for word in cmd)]
assert len(guards)==1 and str(control) in commands[guards[0]],'foreign or ambiguous worker guard'
assert str(control) in commands[owned] and '--owned-child' in commands[owned],'foreign parent'
(control/'paused-test-marker.json').write_text(json.dumps(dict(leaf=os.getpid(),parent=owned,guard=guards[0],outer=outer)))
for pid in (guards[0],owned):
    fd=os.pidfd_open(pid);signal.pidfd_send_signal(fd,signal.SIGSTOP);os.close(fd)
os.kill(os.getpid(),signal.SIGSTOP)
raise RuntimeError('paused fixture unexpectedly resumed instead of package termination')
