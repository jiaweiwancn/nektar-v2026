#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh independent implicit-root replay of saved completed native actions."""
import argparse,os
from pathlib import Path
from run_wj_h4 import THREADS
os.environ.update(THREADS)
import numpy as np
from wj_h15_observer import identity,atomic
from wj_h17_observer import native_load
from publish_wj_h19 import witness_replay


def run(root):
    base=root/'resumption-r1';d=base/'science/C-wide-v1';records=[]
    for nc,nm,folder in ((10,20,'W1-lookup-v1'),(12,24,'W0')):
        p=d/f'temporal-{nc}-signed-fields.npz'
        with np.load(p,allow_pickle=False)as z:arrays={k:z[k]for k in z.files}
        for prefix,report in (('C',d/f'monitor/step-{nc}.report.json'),('M',root/f'science/{folder}/monitor/step-{nm}.report.json')):
            proxy={k:arrays[k]for k in ('physical','geometry')}
            proxy.update({'coarse_'+k:arrays[prefix+'_'+k]for k in ('k','epsilon','G')})
            r=native_load(report);assert r['native_full_audit_status']=='PASS'
            records.append(dict(level=prefix,step=r['steps'],time=r['time'],native_report=identity(report),
                independent_fields=identity(p),records=witness_replay(r,proxy)))
    out=base/'science/native-independent-replays-v1.json';assert not out.exists()
    atomic(out,dict(status='PASS',records=records,normalized_action_tolerance=2e-11,
        dimensional_action_tolerance=2e-13,field_reconstruction='original H19 scaled 2e-11 native witness rules',
        kernel_process_requests=0,method='independent implicit-root complex tangent versus saved compiled native action; exact physical point and geometry owner',
        scope='new completed unperturbed C/M .05/.06 states only; no finite seed or production-duration admission'))
    print('native-independent saved-action replays PASS',max(v['normalized_error']for r in records for v in r['records']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
