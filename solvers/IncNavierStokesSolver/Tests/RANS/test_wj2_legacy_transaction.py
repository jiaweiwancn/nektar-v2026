#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded existing CLS budget control for the shared WJ-2 native paths."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import uuid
from run_cls_epawf_picard import THREADS,require

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--solver',type=Path,required=True)
p.add_argument('--probe',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();os.environ.update(THREADS)
a.output.mkdir(parents=True,exist_ok=True)
output=a.output/('run-'+uuid.uuid4().hex[:12])
subprocess.run([sys.executable,'-B',str(Path(__file__).with_name('test_cls_accepted_awf_native.py')),
    '--solver',str(a.solver),'--probe',str(a.probe),'--output',str(output),
    '--case','Accepted-AWF-2G','--companion','regularized',
    '--companion-reason','Explicit legacy Regularized budget control for WJ-2 shared source, PIMPLE and checkpoint paths'],check=True)
result=json.loads((output/'summary.json').read_text())
require(result['status']=='passed' and result['restart_exact'] and result['rollback_exact'],
        'legacy control must execute and pass its complete native transactions')
print(output/'summary.json')
