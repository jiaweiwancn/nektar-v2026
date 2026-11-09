#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Real Linux process identities and boundary guards for runtime adoption."""
from pathlib import Path
import os
import signal
import subprocess
import sys
import time
import tempfile
from unittest.mock import patch
from handoff_cls_awf_runtime import matches,checked_signal,boundary,live_children,command_line,wait_parent_lock,observe_successor
from continue_cls_epawf_pimple import process_identity
from run_cls_epawf_picard import require


def main():
    require(sys.platform=='linux','Linux authoritative process test required')
    child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'])
    try:
        expected=process_identity(child.pid);require(matches(expected,expected),'live identity missing')
        require(command_line(dict(expected,pid=float(expected['pid'])))==command_line(expected),'JSON float PID changed /proc command lookup')
        wrong=dict(expected,start='wrong-start')
        try:checked_signal(wrong,signal.SIGSTOP)
        except ValueError as e:require('identity changed' in str(e),'wrong guard reason')
        else:raise AssertionError('PID start mismatch signalled a process')
        checked_signal(expected,signal.SIGSTOP)
        for _ in range(100):
            status=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')',1)[1].split()[0]
            if status=='T':break
            time.sleep(.01)
        require(status=='T','identified child was not paused')
        checked_signal(expected,signal.SIGCONT)
        require(any(v['pid']==child.pid for v in live_children(os.getpid())),'live child omitted')
        p=dict(status='running',accepted_steps=8,segment_target=8,native_pid=child.pid,native_identity=expected)
        require(not boundary(p),'live native child was treated as a handoff boundary')
        child.terminate();child.wait(timeout=10)
        require(boundary(p),'completed boundary with ended child missing')
        require(not boundary(dict(p,accepted_steps=7)),'partial block was treated as a completed boundary')
        require(not boundary(dict(p,status='numerically_rejected')),'numerical rejection was eligible for handoff')
        require(not matches(process_identity(child.pid),expected),'dead PID retained ownership')
    finally:
        if child.poll() is None:os.kill(child.pid,signal.SIGCONT);child.terminate();child.wait(timeout=10)
    marker=object()
    with patch('handoff_cls_awf_runtime.CampaignLock',side_effect=[RuntimeError('campaign lock held'),RuntimeError('campaign lock held'),marker]) as factory:
        require(wait_parent_lock(Path('/unused'),1) is marker and factory.call_count==3,'positive lock reacquisition race not handled')
    with patch('handoff_cls_awf_runtime.CampaignLock',side_effect=RuntimeError('campaign lock held')):
        try:wait_parent_lock(Path('/unused'),0)
        except RuntimeError:pass
        else:raise AssertionError('lock timeout bypassed ownership')
    from types import SimpleNamespace
    with tempfile.TemporaryDirectory() as folder:
        output=Path(folder);(output/'trajectory').mkdir();(output/'trajectory/steps.jsonl').write_text('{"global_step": 9, "checkpoint_sha256": "expected"}\n')
        completed={}
        with patch('handoff_cls_awf_runtime.load_progress',side_effect=[FileNotFoundError('startup replacement'),dict(status='running',accepted_steps=9,accepted_time=.005625)]):
            observe_successor(SimpleNamespace(output=output,timeout=5),dict(steps=8),'expected',lambda **kw:completed.update(kw))
        require(completed['status']=='handoff_completed_with_observed_step','startup visibility gap lost a healthy successor receipt')
    print('runtime handoff: start-token guard, real STOP/CONT, live-child exclusion, complete-boundary guard, terminal rejection passed')

if __name__=='__main__':main()
