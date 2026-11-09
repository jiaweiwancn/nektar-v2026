#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from profile_cls_awf_campaign import timing_summary

with tempfile.TemporaryDirectory() as temporary:
    trace=Path(temporary)/'trace.jsonl'
    environment=dict(os.environ);environment.pop('NEKTAR_RANS_TIMING_FILE',None)
    subprocess.run([sys.argv[1]],env=environment,check=True)
    assert not trace.exists()
    subprocess.run([sys.argv[1]],env=dict(environment,NEKTAR_RANS_TIMING_FILE=str(trace)),check=True)
    result=timing_summary(trace)
    assert result['physical_steps']==1
    assert set(result['groups'])=={'step/native.step','step/normal','step/manual_stop','step/exception','block_overhead/process'}
    assert all(g['calls']==1 for g in result['groups'].values())
print('opt-in, nesting, context, manual stop, exception unwind and time partition passed')
