"""Execute the one frozen R1 kernel proof/witness suite on Linux CPU50."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

base=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised');root=base/'m8-r1'
tools=root/'kernel-tools-001';out=root/'kernel-runner-001';out.mkdir(exist_ok=False)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
expected={'cls_epawf_m8_r1_kernel.py':'43214f63aa5a96400a85dd8ded0416055ed87d668a051baf5cd125595e9fe293',
    'kernel-plan.json':'fdb75bce68563ebac786d4e4fc8109fc635183469d360ee0edfc94fec17b2da0',
    'candidate-contract.md':'082a791fb086a08bd83285e68ae5e05923fdb688c7919410263bb894a6f22313'}
for name,value in expected.items():assert sha(tools/name)==value,name
env=dict(os.environ);env.update({key:'1' for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
    'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')});env['PYTHONDONTWRITEBYTECODE']='1';env['TMPDIR']=str(root/'tmp')
argv=['taskset','-c','50','python3','-B',str(tools/'cls_epawf_m8_r1_kernel.py'),
    '--archive-source',str(root/'archive-source'),'--runs',str(base),
    '--ledger',str(root/'archive-attempt-003/archive-ledger.json'),
    '--plan',str(tools/'kernel-plan.json'),'--contract',str(tools/'candidate-contract.md'),
    '--output',str(root/'kernel-attempt-001')]
receipt=dict(schema=1,status='running',command=argv,input_hashes=expected,runner_sha256=sha(__file__),
    numerical_workers=1,numerical_threads=1,cpu=50,production_adopted=False,trajectory_launched=False,
    connection='tailscale:100.105.247.110',started_utc=datetime.now(timezone.utc).isoformat())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');start=time.monotonic()
try:
    with (out/'run.log').open('w') as log:code=subprocess.run(argv,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=1200).returncode
except subprocess.TimeoutExpired:code=124
receipt.update(status='passed' if code==0 else 'failed',exit_code=code,seconds=time.monotonic()-start,
    log_sha256=sha(out/'run.log'),finished_utc=datetime.now(timezone.utc).isoformat())
(out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print((out/'run.log').read_text()[-3000:]);print(json.dumps(receipt));raise SystemExit(code)
