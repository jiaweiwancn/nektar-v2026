#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One isolated H20 ledger; selected sealed states and in-place H19 runtime."""
import os,shutil,subprocess,time
from pathlib import Path
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import atomic,identity,verify
ROOT=Path('/home/jay/build/wj-h20-three-level-20261115')
OLD=Path('/home/jay/build/wj-h19-target-timestep-20261115')

def main():
    ROOT.mkdir()
    for n in ('inputs','science','operations','tmp'): (ROOT/n).mkdir()
    book=dict(schema='WJH20LedgerV1',closed=False,ceiling_seconds=7200,setup_charged_seconds=60,closure_reserve_seconds=300,
        starting_revision='fd009626de46e164c0f645b9879a68d2102ae11f',target_reserved_attempts=0,bootstrap_reserved=0,strict_reservations=[],
        entries={'selected-bootstrap-v1':dict(status='RESERVED',phase='operation',reserved_seconds=120),
                 'network-coordination-unknown-v1':dict(status='UNKNOWN_FULL_CHARGE',phase='operation',reserved_seconds=180,seconds=180,reason='conservative unmetered source/network/receipt coordination; never reclaimed')})
    atomic(ROOT/'accounting.json',book);began=time.monotonic()
    try:
        fs=subprocess.check_output(['findmnt','-n','-o','FSTYPE,SOURCE,TARGET','-T',str(ROOT)],text=True).strip()
        assert fs.startswith('ext4 ')and shutil.disk_usage(ROOT).free>10*1024**3
        p=ROOT/'tmp/fsync-readback';data=bytes(range(256))*4096
        with p.open('xb')as f:f.write(data);f.flush();os.fsync(f.fileno())
        assert p.read_bytes()==data
        index=load(OLD/'authoritative-index-v1.json')
        for v in index.values():verify(v)
        oldbook=load(verify(index['accounting']));assert oldbook['closed']and oldbook['total_charged_seconds']==1259.9079674030363
        manifest=load(verify(index['publication']));published={v['path']:v for v in manifest['files']}
        def sealed(relative):
            p=OLD/relative;verify(published[str(p)]);return load(p)
        records=sealed('publication-v1/state-history-bindings.json');fine=sealed('publication-v1/controls/inputs--anchors-v2.json')
        states={'C':{},'M':{}}
        for v in records:
            n=v['step'];folder=OLD/'science/coarse/monitor';completion=folder/f'step-{n}.completion.json'
            c=load(verify(published[str(OLD/f'publication-v1/coarse/step-{n}.completion.json')]))
            # The copied publication pins the original completion's payloads.
            for key in ('checkpoint','callback','report'):verify(c[key])
            states['C'][str(n)]={k:c[k]for k in ('checkpoint','callback','report')}
            states['C'][str(n)].update(completion=identity(completion),qualified_commit=v['qualification_commit'],audit_report=c['report'])
            assert c['checkpoint']==v['checkpoint_file']and c['report']==v['producer']
            verify(v['qualification_commit'])
        states['M']={n:dict(v)for n,v in fine['fine'].items()}
        for n,v in states['M'].items():
            for key in ('checkpoint','callback','report','audit_report','qualified_commit'):verify(v[key])
            if n=='4':
                directory=Path(v['report']['path']).parent;v['completion']=identity(directory/'step-4.completion.json')
        original=fine['original']
        for k in ('xml','producer','checkpoint'):verify(original[k])
        assert original['checkpoint']['sha256']=='d27ee0e2f9b9bf8639d47449c5e8e338212363417faca7bf014d1008f44f69a3'
        lineage=sealed('publication-v1/bootstrap/lineage.json')
        assert lineage['coarse_checkpoint']['sha256']=='9d51909af6a22d28a5b1ade258b982f19e491cbc8057fd93304c11e681934889';verify(lineage['coarse_checkpoint'])
        runtime=load(OLD/'science/runtime-build-v1.json')
        assert [v['ELF']['sha256']for v in runtime['executables']]==['067ccdf0954efb6ea9921887edd911c57da910be1c82ffddaadba4f7f544d3e8','637656771bdeb4b94d1487748489357ab58e0030a80cd52d5f55289018ed375c']
        for v in runtime['executables']:verify(v['ELF']);verify(v['wrapper'])
        for v in load(OLD/'science/runtime-dependencies-v1.json')['libraries']:verify(v['actual'])
        compiled=sealed('publication-v1/compiled-source-summary.json')
        for v in compiled['dependencies']:verify(v)
        atomic(ROOT/'inputs/selected.json',dict(status='PASS',filesystem=fs,free_bytes=shutil.disk_usage(ROOT).free,states=states,original=original,
            C_CP0=lineage['coarse_checkpoint'],C_CP0_report=lineage['producer'],C_CP0_xml=lineage['coarse_xml'],common_support=fine['common_support'],
            H19_closed_accounting=index['accounting'],H19_publication=index['publication'],runtime_build=identity(OLD/'science/runtime-build-v1.json'),runtime_dependencies=identity(OLD/'science/runtime-dependencies-v1.json'),
            compiled=published[str(OLD/'publication-v1/compiled-source-summary.json')],runtime_root=str(OLD),cache_manifest=identity(OLD/'science/target-cache-compatibility.json')))
        shutil.copytree(OLD/'source',ROOT/'source')
        e=dict(status='PASS',phase='operation',reserved_seconds=120,seconds=time.monotonic()-began,native_work=0)
    except Exception as error:
        e=dict(status='OBSTRUCTED',phase='operation',reserved_seconds=120,seconds=time.monotonic()-began,reason=repr(error));raise
    finally:
        book['entries']['selected-bootstrap-v1']=e;atomic(ROOT/'accounting.json',book,False);atomic(ROOT/'operations/bootstrap-receipt.json',e)
    print(e,flush=True)
if __name__=='__main__':main()
