#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read campaign ledgers/processes and archive verified immutable checkpoint copies."""
import argparse
import datetime
import json
import os
import shutil
import subprocess
from pathlib import Path
from continue_cls_epawf_pimple import header,load_progress,process_identity,atomic
from run_cls_epawf_picard import load,sha

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    root=a.root.resolve();base=root.parent;cases={}
    paths={'Upstream05graded':root/'selected-U05-graded/trajectory', 'Upstream05one':root/'selected-U05-one/trajectory',
        'UniformG':base/'cls-epawf-height-sampler-20261107/selected-G/trajectory',
        'UniformH':base/'cls-epawf-height-sampler-20261107/selected-H/trajectory',
        'historicalUniform':base/'cls-epawf-uniform-ep-flat-20261106/uniform-campaign/trajectory',
        'terminal103':base/'cls-epawf-projection-repair-20261107/repaired-campaign/trajectory'}
    for name,path in paths.items():
        if not (path/'progress.json').exists():continue
        progress=load_progress(path/'progress.json');copy=a.output/name;copy.mkdir()
        atomic(copy/'progress.json',progress)
        # A live checkpoint can rotate. Open the payload once, verify the
        # recorded hash, then retain the exact opened bytes atomically.
        checkpoint=Path(progress['checkpoint']);payload=None
        try:payload=(checkpoint/'State.bin').read_bytes()
        except FileNotFoundError:pass
        if payload is not None:
            import hashlib
            digest=hashlib.sha256(payload).hexdigest()
            if digest==progress['checkpoint_sha256']:
                cp=copy/'checkpoint';cp.mkdir();(cp/'State.bin').write_bytes(payload)
                receipt=header(cp)
            else:receipt=dict(status='checkpoint_rotated_before_copy',observed_sha256=digest)
        else:receipt=dict(status='checkpoint_rotated_before_copy')
        driver=process_identity(progress['pid']) if progress.get('pid') else None
        native=process_identity(progress['native_pid']) if progress.get('native_pid') else None
        identity=progress.get('native_identity')
        record=dict(progress=progress,checkpoint_copy=receipt,driver_identity=driver,native_identity=native,
            driver_matches=bool(driver and str(driver['start'])==str(progress.get('process_start'))),
            native_matches=(native==identity if identity else bool(native and native['executable']==progress.get('executable'))))
        if native:
            try:
                env=(Path('/proc')/str(native['pid'])/'environ').read_bytes().split(b'\0')
                record['computational_threads']={v.split(b'=',1)[0].decode():v.split(b'=',1)[1].decode() for v in env
                    if v.split(b'=',1)[0].decode() in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
            except FileNotFoundError:record['native_finished_during_capture']=True
        identity_file=path/'runtime/identity.json'
        if identity_file.exists():
            record['frozen_runtime_manifest_sha256']=sha(identity_file)
            shutil.copy2(identity_file,copy/'runtime-identity.json')
        cases[name]=record
    screening={}
    for c in ('U05-graded','U05-one','U20-graded','U20-one'):
        path=root/'screening'/c/'result.json'
        if path.exists():
            d=load(path);screening[c]={k:d.get(k) for k in ('status','steps','time','failure','actual_h0_sha256','accepted_checkpoint')}
    workers=[]
    for line in subprocess.check_output(['ps','-eo','pid,nlwp,pcpu,args'],text=True).splitlines()[1:]:
        if any(t in line for t in ('/IncNavierStokesSolver ','/RANSWallResolutionProbe ','cc1plus','cmake --build')) and 'capture_cls' not in line:
            workers.append(line.strip())
    record=dict(schema=1,captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),cases=cases,screening=screening,
        worker_processes=workers,production_rank_policy=1,computational_threads_per_runtime=1,
        historical_experiments_modified=False)
    atomic(a.output/'status.json',record);print(a.output/'status.json')

if __name__=='__main__':main()
