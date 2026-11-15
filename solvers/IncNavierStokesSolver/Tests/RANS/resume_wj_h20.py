#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Explicitly authorized continuation, preserving the sealed H20 prefix."""
import os,shutil,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic

OLD=Path('/home/jay/build/wj-h20-three-level-20261115')
ROOT=Path('/home/jay/build/wj-h20-three-level-resume-20261115')
AUTHORIZATION='Resume and close **obstructed H20 by doing whatever necessary.**'

def initialize():
    index=load(OLD/'authoritative-index-v1.json')
    for pin in index.values():verify(pin)
    prior=load(verify(index['accounting']));assert prior['closed']
    assert prior['total_charged_seconds']==714.0562545959547
    assert prior['target_levels']==['M']and prior['bootstrap_reserved']==0
    outcome=load(OLD/'science/M/outcome.json');intent=load(OLD/'science/M/intent.json')
    assert outcome['native_exit']==1 and outcome['producer']is None and intent['ack_timeout_seconds']==180
    assert not(OLD/'science/M/case.evolution.json').exists()and not(OLD/'science/M/monitor/startup.ready.json').exists()
    assert not any(prior.get(k)for k in ('coverage_halt','screen_halt','runtime_halt'))
    ROOT.mkdir()
    for name in ('inputs','science','operations','tmp'):(ROOT/name).mkdir()
    native=sum(e.get('seconds',e['reserved_seconds'])for e in prior['entries'].values()if e['phase']in ('qualification','target','target_readback'))
    book=dict(schema='WJH20ResumptionLedgerV1',closed=False,ceiling_seconds=7200,setup_charged_seconds=60,closure_reserve_seconds=300,
        carryforward_seconds=prior['total_charged_seconds'],carryforward_native_seconds=native,carryforward_accounting=index['accounting'],historical_charged_attempt_capacity=4,historical_strict_reads=1,
        target_reserved_attempts=0,bootstrap_reserved=0,strict_reservations=[],entries={
            'resumption-bootstrap-v1':dict(status='RESERVED',phase='operation',reserved_seconds=120),
            'coordination-transport-v1':dict(status='FULLY_CHARGED',phase='operation',reserved_seconds=120,seconds=120)})
    atomic(ROOT/'accounting.json',book);start=time.monotonic()
    try:
        assert shutil.disk_usage(ROOT).free>10*1024**3
        p=ROOT/'tmp/fsync-readback';data=bytes(range(256))*4096
        with p.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
        assert p.read_bytes()==data
        selected=load(OLD/'inputs/selected.json')
        for v in selected['states'].values():
            for state in v.values():
                for key in ('checkpoint','report','callback','qualified_commit','audit_report'):verify(state[key])
        atomic(ROOT/'inputs/selected.json',selected)
        atomic(ROOT/'resumption.json',dict(schema='WJH20ExplicitResumptionV1',authorization=AUTHORIZATION,starting_revision='e9fc43a5cb35debc0d512138d5230ee5c368b140',
            sealed_prefix=index,carryforward_seconds=prior['total_charged_seconds'],historical_unknown_capacity=4,
            authorization_scope='one corrected M launch followed by C and F; maximum22 actual governing attempts, one CP0 and three terminal strict reads; no scientific input changes',
            counters='old launch/unknown-capacity/read retained separately; no accepted governing state is rerun',
            historical_halt='retained immutable configuration rejection; new authorization does not override a physical/coverage halt'))
        shutil.copytree(OLD/'source',ROOT/'source')
        entry=dict(status='PASS',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start)
    except Exception as error:
        entry=dict(status='OBSTRUCTED',phase='operation',reserved_seconds=120,seconds=time.monotonic()-start,reason=repr(error));raise
    finally:
        book['entries']['resumption-bootstrap-v1']=entry;atomic(ROOT/'accounting.json',book,False);atomic(ROOT/'operations/bootstrap-receipt.json',entry)
    print(entry,flush=True)

if __name__=='__main__':initialize()
