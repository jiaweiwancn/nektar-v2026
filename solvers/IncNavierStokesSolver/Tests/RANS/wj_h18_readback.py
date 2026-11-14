# SPDX-License-Identifier: MIT
"""Select a committed verification, including a verified pending-restart input."""
from pathlib import Path
from run_wj_h4 import load
from wj_h15_observer import identity,verify


def selection(root):
    for name in ('B','A'):
        folder=Path(root)/'science/target'/name;directory=folder/'monitor';cursor=directory/'latest-fully-verified-checkpoint.json'
        if not cursor.exists():continue
        q=load(cursor);kind=q['schema'];assert kind in ('WJQualificationCommitV1','WJStartupVerificationCommitV1'),'prepared/uncommitted pointer is not verified'
        if kind=='WJStartupVerificationCommitV1':
            if int(q['binding']['absolute_step'])==0:return None  # genuine original CP0 producer, alpha0
            commit=directory/'startup.verification-commit.json';c=load(verify(q['completion']))
        else:
            commit=directory/f"step-{int(q['binding']['absolute_step'])}.qualification-commit.json";c=q
        assert identity(commit)['sha256']==identity(cursor)['sha256'],'published verified cursor differs from native commit'
        for k in ('checkpoint','callback','report'):verify(c[k])
        return folder,dict(operation='open_checks',inner_stress_model='Linear',native_report=c['report']['path'],native_report_sha256=c['report']['sha256'],
            native_checkpoint=str(Path(c['checkpoint']['path']).parent),checkpoint_sha256=c['checkpoint']['sha256'],completed_sidecar=c['callback'],qualified_cursor=identity(commit),
            weak=True,weak_compact=True,assembly=False,reuse=False,retain_scalar_observations=False)
    return None
