#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Export read-only Linux metadata for CT15 publication on the source Mac.

No solver, probe, FE calculation or scientific reclassification is performed.
Closed Linux receipts remain authoritative. This administrative helper belongs
to the fixed transport/publication reserve; it does not open a numerical family
or supply worker admission. In particular an unclosed queue stays unclosed.
"""
import argparse, json, os, resource, subprocess, time
from pathlib import Path
from ct_accepted15_collect import START,digest,historical_ledger,read,write


def main():
    start=time.monotonic();own=time.process_time();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    resource.setrlimit(resource.RLIMIT_CPU,(60,65))
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','ct14','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();a.root=a.root.resolve();a.output.mkdir(exist_ok=False)
    names={'request-manifest.json','completed-artifacts.json','results.json','progress.json',
        'execution.json','package-budget-watch.json','live-budget.json','resource-final.json',
        'worker-guard-status.json','queue-deadline-watch.json','queue-cancellation.json',
        'retention-controls.json','retention-controls-seal.json','heldout-admission.json',
        'frozen-model-before-heldout.json','frozen-complete-model-before-heldout.json',
        'audit-replay-2.json','audit-replay-4.json','replay-proof-2.json','replay-proof-4.json'}
    artifacts={};resources={};unclosed=[];watches={};directory_bytes={'.':0}
    for folder in sorted(a.root.iterdir()):
        if not folder.is_dir() or folder.name in ('source','build'):continue
        for path in sorted(folder.iterdir()):
            if path.is_file() and path.name in names:
                artifacts[str(path.relative_to(a.root))]=dict(bytes=path.stat().st_size,sha256=digest(path))
        if folder.name.startswith('resources-'):
            value=read(folder/'execution.json')
            if value is None:unclosed.append(folder.name)
            else:resources[folder.name]=value
            for name in ('package-budget-watch.json','live-budget.json','worker-guard-status.json','resource-final.json','queue-deadline-watch.json','queue-cancellation.json'):
                q=read(folder/name)
                if q is not None:watches[str((folder/name).relative_to(a.root))]=q
    # One walk measures actual logical live occupancy, including all attempts,
    # source/build/runtime and staging. Symlink targets are not counted twice.
    for folder,dirs,files in os.walk(a.root):
        for name in files:
            path=Path(folder)/name
            try:n=path.lstat().st_size
            except FileNotFoundError:continue
            ancestor=path.parent
            while True:
                key=str(ancestor.relative_to(a.root));directory_bytes[key]=directory_bytes.get(key,0)+n
                if ancestor==a.root:break
                ancestor=ancestor.parent
    source=a.root/'source'
    revision=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
    tree=subprocess.check_output(['git','-C',str(source),'status','--short'],text=True)
    historical=historical_ledger(a.root,a.ct14)
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result=dict(schema='CT15ReadOnlyLinuxPublicationSnapshotV1',external_root=str(a.root),observed_unix=time.time(),
        original_allocation_start_unix=START,whole_elapsed_wall=time.time()-START,
        source_revision=revision,source_working_tree=tree,
        artifacts=artifacts,closed_resources=resources,unclosed_resources=unclosed,watches=watches,
        directory_logical_bytes=directory_bytes,exception_reserves=read(a.root/'exception-reserves.json'),
        historical_ledger=historical,
        helper_usage=dict(self_CPU=time.process_time()-own,reaped_children_CPU=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,wall=time.monotonic()-start,
            accounting='inside the existing 240-CPU-second source/transport/publication reserve; do not add a second reserve'),
        scope='read-only metadata export only; no numerical execution or scientific/worker-admission PASS; missing/failed receipts are preserved explicitly; logical bytes differ from physical allocated disk and cumulative writes')
    write(a.output/'snapshot.json',result)
    write(a.output/'seal.json',dict(files={'snapshot.json':digest(a.output/'snapshot.json')}))
    print(json.dumps(dict(status='metadata_exported',closed_count=len(resources),unclosed=unclosed,CPU=sum(r['process_cpu_seconds'] for r in resources.values()),logical_live_bytes=directory_bytes['.'])))


if __name__=='__main__':raise SystemExit(main())
