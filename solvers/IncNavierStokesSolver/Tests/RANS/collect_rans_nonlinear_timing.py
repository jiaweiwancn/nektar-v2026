#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Completed-stage CPU costs; suspended wall time is reported separately."""
import argparse
import collections
import json
from pathlib import Path
from rans_request_manifest import RequestManifest,sha
from rans_qualification import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();rows=[]
    for folder in sorted((a.root/'campaigns').iterdir()):
        if not (folder/'request-manifest.json').exists():continue
        q=json.loads((folder/'request-manifest.json').read_text())['request'];seal=RequestManifest(folder,q,True)
        for endpoint in sorted(folder.glob('T-*')):
            if not seal.completed(endpoint) or not (endpoint/'timing.jsonl').exists():continue
            metrics=json.loads((endpoint/'metrics.json').read_text());native=json.loads((endpoint/'case.evolution.json').read_text())
            n=len(native.get('pimple_attempts',[]));cpu=collections.Counter();wall=collections.Counter();counts=collections.Counter()
            for line in (endpoint/'timing.jsonl').read_text().splitlines():
                item=json.loads(line);name=item['name'];cpu[name]+=item['exclusive_process_cpu_seconds'];wall[name]+=item['exclusive_wall_seconds'];counts[name]+=1
            rows.append(dict(label=folder.name,configuration={k:q[k] for k in ('model','policy','wall')},time=metrics['time'],dt=q['dt'],accepted_steps_in_segment=n,wall_seconds=metrics['wall_seconds'],instrumented_exclusive_process_cpu_seconds=sum(cpu.values()),instrumented_CPU_per_step=sum(cpu.values())/n if n else None,exclusive_CPU_by_scope=dict(cpu),exclusive_wall_by_scope=dict(wall),calls_by_scope=dict(counts),timing_sha256=sha(endpoint/'timing.jsonl'),source_native_sha256=sha(endpoint/'case.evolution.json')))
    result=dict(status='pass' if rows else 'not_completed',cases=rows,definition='sum exclusive process-CPU intervals only; inclusive intervals are not added. Instrumented CPU covers initialization, assembly, diagnostics, solves and output where instrumented. Wall time includes resource waits/suspension.',scope='cost evidence only, not a numerical accuracy or conservation gate')
    if a.output.exists():raise FileExistsError(a.output)
    write_json(a.output,result);return 0 if rows else 1


if __name__=='__main__':raise SystemExit(main())
