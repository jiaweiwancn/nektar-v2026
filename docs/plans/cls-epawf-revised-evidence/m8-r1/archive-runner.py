#!/usr/bin/env python3
"""Exact one-worker archive reconstruction wrapper; no solver invocation."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True, type=Path)
parser.add_argument('--attempt', required=True)
args = parser.parse_args()
root, attempt = args.root, args.attempt
script, source, out = root/'cls_epawf_m8_r1_archive.py', root/'archive-source', root/attempt
threads = {key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}
os.environ.update(threads)
os.environ['TMPDIR'] = str(root/'tmp')
command = ['timeout','1800','taskset','-c','48','python3','-B',str(script),'--source',str(source),
    '--archive','/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8/flow-001/qualification/control-coarse',
    '--output',str(out)]
receipt = dict(schema=1, started_utc=datetime.now(timezone.utc).isoformat(), command=command,
    driver_sha256=hashlib.sha256(script.read_bytes()).hexdigest(),
    wrapper_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    threads=threads, worker_limit=1, CPU='48', solver_or_trajectory_launched=False)
start = time.monotonic()
with (root/(attempt+'.log')).open('w') as log:
    result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=os.environ.copy())
if out.exists():
    shutil.copy2(script,out/'archive-driver.py')
receipt.update(exit_code=result.returncode, seconds=time.monotonic()-start,
    finished_utc=datetime.now(timezone.utc).isoformat(),
    log_sha256=hashlib.sha256((root/(attempt+'.log')).read_bytes()).hexdigest())
(root/(attempt+'-runner.json')).write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
print((root/(attempt+'.log')).read_text()[-12000:])
raise SystemExit(result.returncode)
