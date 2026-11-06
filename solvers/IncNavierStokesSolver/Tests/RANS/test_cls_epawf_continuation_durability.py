#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused Linux campaign ownership and crash-window recovery checks."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from continue_cls_epawf_pimple import CampaignLock,StepLedger,process_identity,assert_native_stopped
from run_cls_epawf_picard import require


def fails(fn,word):
    try:fn()
    except (RuntimeError,ValueError) as error:
        require(word in str(error),str(error));return
    raise AssertionError('expected rejection: '+word)


def main():
    require(sys.platform=='linux','authoritative durability checks run on Linux')
    with tempfile.TemporaryDirectory(prefix='pimple-lock-') as folder:
        root=Path(folder);lock=CampaignLock(root)
        fails(lambda:CampaignLock(root),'campaign lock held')
        child=subprocess.Popen([sys.executable,'-c','import time;time.sleep(60)'],pass_fds=(lock.stream.fileno(),))
        try:
            expected=process_identity(child.pid);require(expected is not None,'native identity missing')
            lock.close() # simulate the driver's death, preserving its child
            fails(lambda:CampaignLock(root),'campaign lock held')
            fails(lambda:assert_native_stopped(dict(native_identity=expected)),'still alive')
            # PID reuse with a different start token must not target a process.
            assert_native_stopped(dict(native_identity=dict(expected,start='different-token')))
        finally:child.terminate();child.wait(timeout=10)
        replacement=CampaignLock(root);replacement.close()
        assert_native_stopped(dict(native_identity=expected))
        path=root/'steps.jsonl';ledger=StepLedger(path)
        prefix=dict(prefix_steps=4,prefix_integral=.25)
        record=dict(global_step=5,checkpoint_sha256='real-checkpoint-sha',attempt=dict(k_bound=dict(integral_total=.125)))
        require(ledger.append(record),'first publication missing')
        # Crash after durable append, before progress replace. Recover twice.
        for _ in range(2):
            ledger=StepLedger(path);require(not ledger.append(record),'recovery duplicated record')
            require(ledger.budget(prefix,5)==.375,'recovery duplicated correction budget')
        require(len(path.read_text().splitlines())==1,'duplicate accepted-step record')
        fails(lambda:ledger.append(dict(record,checkpoint_sha256='conflict')),'identity conflict')
        # A duplicate old journal line and incomplete final write are repaired
        # deterministically, without admitting an incomplete accepted step.
        with path.open('a') as stream:stream.write(json.dumps(record)+'\n'+ '{"global_step":6')
        ledger=StepLedger(path)
        require(len(ledger.records)==1 and ledger.budget(prefix,5)==.375,'journal recovery differs')
        require(len(path.read_text().splitlines())==1,'journal not canonical after reconciliation')
        terminal=root/'terminal';terminal.mkdir()
        progress=terminal/'progress.json';progress.write_text(json.dumps(dict(status='numerically_rejected',checkpoint=None)))
        original=progress.read_bytes()
        command=[sys.executable,'-B',str(Path(__file__).with_name('continue_cls_epawf_pimple.py')),
            '--solver',str(root/'never-execute'),'--output',str(terminal),'--resume']
        rejected=subprocess.run(command,capture_output=True,text=True,timeout=10)
        require(rejected.returncode and 'numerical rejection is terminal' in rejected.stderr,'ordinary resume retried numerical rejection')
        require(progress.read_bytes()==original,'terminal progress changed')
    print(json.dumps(dict(simultaneous_resume_excluded=True,orphan_native_excluded=True,
        native_identity_checked=True,idempotent_step_and_budget_recovery=True,partial_final_record_repaired=True,numerical_rejection_terminal=True)))


if __name__=='__main__':main()
