#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Immutable progress/checkpoint evidence; read-only toward all live campaigns."""
import argparse
import datetime
import gzip
import json
from pathlib import Path
import shutil
import subprocess
from run_cls_epawf_picard import load,sha,write_json,THREADS
from capture_cls_epawf_projection_repair import snapshot


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    r=a.root.resolve();runs=r.parent
    out=dict(captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),threads=THREADS,ranks_per_native=1,
        root=str(r),campaigns={},artifacts=[],source_revision=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip())
    for name,campaign in [('upstream',runs/'cls-epawf-projection-repair-20261107/repaired-campaign'),
                          ('uniform',runs/'cls-epawf-uniform-ep-flat-20261106/uniform-campaign'),
                          ('selected-G',r/'selected-G'),('selected-H',r/'selected-H')]:
        if not (campaign/'trajectory/progress.json').exists():continue
        record=snapshot(campaign,a.output/(name+'-checkpoint'));out['campaigns'][name]=record
        if record['progress']['status']=='numerically_rejected':
            segment=Path(record['progress']['segment']);path=segment/'case.evolution.json'
            if not path.exists():path=path.with_suffix(path.suffix+'.gz')
            target=a.output/(name+'-terminal'+''.join(path.suffixes));shutil.copy2(path,target)
            d=json.load(gzip.open(path,'rt') if path.suffix=='.gz' else path.open());attempt=d['pimple_attempts'][-1]
            record['first_failure']={k:attempt[k] for k in ('time','failure','scheduled_complete','flow_domain','flow_ground') if k in attempt}
            record['terminal_report']=dict(path=str(target),sha256=sha(target),contains='actual rejected_state and last_stage_history; no retry')
    # Index native checkpoints/packets/reports and logs, leaving dense matrices
    # and fields on Linux. Do not hash mutable active files as immutable evidence.
    for family in ('screening','focused','dt-half'):
        root=r/family
        if not root.exists():continue
        for path in sorted(root.rglob('*')):
            if not path.is_file() or any(x in path.parts for x in ('runtime','probe-runtime','lib','tmp')):continue
            if path.name in ('progress.json','solver.log') or path.suffix not in ('.json','.bin','.gz','.xml','.log'):continue
            out['artifacts'].append(dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path)))
    for path in sorted(a.output.rglob('State.bin')):
        out['artifacts'].append(dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path)))
    write_json(a.output/'status-and-manifest.json',out);print(a.output/'status-and-manifest.json')


if __name__=='__main__':main()
