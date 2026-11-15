#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""New finite H22 ledger and selected immutable inputs; never a CP0 producer."""
import base64,io,json,os,shutil,subprocess,tarfile,time
from pathlib import Path
from run_wj_h4 import THREADS,load,save
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from wj6_resources import scientific_workers

ROOT=Path('/home/jay/build/wj-h22-early-phase-replay-20261115')
OLD=Path('/home/jay/build/wj-h19-target-timestep-20261115')

def bootstrap(bundle):
    ROOT.mkdir()
    for name in ('science','operations','inputs','tmp'): (ROOT/name).mkdir()
    book=dict(schema='WJH22LedgerV1',closed=False,ceiling_seconds=4500,setup_charged_seconds=60.,closure_reserve_seconds=300.,
              target_reserved_attempts=0,target_levels=[],strict_reservations=[],entries={
              'bootstrap-selected':dict(status='RESERVED',phase='operation',reserved_seconds=120.)})
    save(ROOT/'accounting.json',book);start=time.monotonic()
    try:
        allocation=scientific_workers();assert allocation['workers']+8<=72
        fs=subprocess.check_output(['findmnt','-n','-o','FSTYPE,SOURCE,TARGET','-T',str(ROOT)],text=True).strip()
        assert fs.startswith('ext4 ')and shutil.disk_usage(ROOT).free>12*1024**3
        probe=ROOT/'tmp/fsync-readback';data=bytes(range(256))*4096
        with probe.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
        assert probe.read_bytes()==data
        shutil.copytree(OLD/'source',ROOT/'source',symlinks=True)
        with tarfile.open(fileobj=io.BytesIO(base64.b64decode(bundle)),mode='r:gz')as tar:
            for member in tar.getmembers():
                assert member.isfile()and not Path(member.name).is_absolute()and '..'not in Path(member.name).parts
                path=ROOT/member.name;path.parent.mkdir(parents=True,exist_ok=True)
                # Source is a new owned clone, never the historical capsule.
                path.write_bytes(tar.extractfile(member).read())
        selected=load(ROOT/'inputs/selected.json')
        for v in selected['prior_indices'].values():
            for key in ('publication','accounting'):verify(v[key])
            assert load(v['accounting']['path'])['closed']
        for anchor in selected['cp0'].values():
            for value in anchor.values():
                if isinstance(value,dict)and 'path'in value:verify(value)
        for level,states in selected['oracles'].items():
            for anchor in states.values():
                for value in anchor.values():
                    if isinstance(value,dict)and 'path'in value:verify(value)
        atomic(ROOT/'inputs/setup.json',dict(status='PASS',filesystem=fs,free_bytes=shutil.disk_usage(ROOT).free,
               allocation=allocation,threads=THREADS,spatial_or_history_initializations=0,bootstrap_driver=identity(__file__)))
        entry=dict(status='PASS',phase='operation',reserved_seconds=120.,seconds=time.monotonic()-start)
    except Exception as error:
        entry=dict(status='OBSTRUCTED',phase='operation',reserved_seconds=120.,seconds=time.monotonic()-start,reason=repr(error));raise
    finally:
        book['entries']['bootstrap-selected']=entry;save(ROOT/'accounting.json',book);save(ROOT/'operations/bootstrap-receipt.json',entry)
    print(entry,flush=True)

if __name__=='__main__':
    import sys
    bootstrap(sys.stdin.read())
