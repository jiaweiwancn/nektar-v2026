#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Create one isolated H19 allocation and select immutable prior artifacts."""
import os,shutil,subprocess,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,H14

ROOT=Path('/home/jay/build/wj-h19-target-timestep-20261115')
OLD=Path('/home/jay/build/wj-h18-checkpoint-audits-20261114')

def main():
    ROOT.mkdir()
    for n in ('science','operations','inputs','tmp'): (ROOT/n).mkdir()
    book=dict(schema='WJH19LedgerV1',closed=False,ceiling_seconds=5400,setup_charged_seconds=60,closure_reserve_seconds=300,
        starting_revision='a8fc3f50eab7b50fcdec235a22f82a140ea7ee06',target_reserved_attempts=0,bootstrap_reserved=0,
        entries={'selected-bootstrap-v1':dict(status='RESERVED',phase='operation',reserved_seconds=120)})
    atomic(ROOT/'accounting.json',book);start=time.monotonic()
    try:
        fs=subprocess.check_output(['findmnt','-n','-o','FSTYPE,SOURCE,TARGET','-T',str(ROOT)],text=True).strip()
        assert fs.startswith('ext4 ')and shutil.disk_usage(ROOT).free>12*1024**3
        data=bytes(range(256))*4096;p=ROOT/'tmp/fsync-readback'
        with p.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
        assert p.read_bytes()==data
        index=load(OLD/'authoritative-index-v1.json')
        for v in index.values():verify(v)
        accounting=load(verify(index['accounting']));assert accounting['closed']and accounting['total_charged_seconds']==2305.3583236689738
        bindings=load(H14/'publication-v2/state-bindings.json')
        original=next(v for v in bindings['states']if v['case']=='cp0')
        for k in ('xml','checkpoint','producer'):verify(original[k])
        fine={}
        for step,folder,commit in ((2,OLD/'science/target/B/monitor','startup.verification-commit.json'),(4,OLD/'science/target/B/monitor','step-4.qualification-commit.json')):
            q=load(folder/commit)
            c=load(verify(q['completion']))if step==2 else q
            fine[str(step)]={k:c[k]for k in ('checkpoint','report','callback')}
            fine[str(step)]['qualified_commit']=identity(folder/commit)
            for k in ('checkpoint','report','callback'):verify(fine[str(step)][k])
        support=identity(OLD/'inputs/target/common-support.json');verify(support)
        atomic(ROOT/'inputs/selected.json',dict(original=original,fine=fine,common_support=support,H18_closed_accounting=index['accounting'],
            cache_manifest=identity(OLD/'science/target-cache-compatibility.json'),filesystem=fs,free_bytes=shutil.disk_usage(ROOT).free))
        shutil.copytree(OLD/'source',ROOT/'source')
        shutil.copytree(OLD/'build',ROOT/'build')
        e=dict(status='PASS',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,native_work=0)
    except Exception as error:
        e=dict(status='OBSTRUCTED',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,reason=repr(error));raise
    finally:
        book['entries']['selected-bootstrap-v1']=e
        # Initial bootstrap is single-writer, before any process admission.
        atomic(ROOT/'accounting.json',book,False);atomic(ROOT/'operations/bootstrap-receipt.json',e)
    print(e,flush=True)
if __name__=='__main__':main()
