#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Finish the recorded prototype run without repeating its completed tests.

A fresh reproduction uses qualify-linux.py instead. This helper reuses that
committed script's legacy replay checks and repairs only the prototype's final
identity-receipt lookup. It never edits a native state or checkpoint.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
sys.dont_write_bytecode=True
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
r=a.root.resolve(strict=True);s=r/'source';t=s/'solvers/IncNavierStokesSolver/Tests/RANS'
sys.path.insert(0,str(t))
from run_cls_epawf_picard import readback,sha
from run_cls_epawf_pimple import identity
os.environ.update(LD_LIBRARY_PATH=str(r/'runtime-qualified/lib'),TMPDIR=str(r/'tmp'),PYTHONPATH=str(t),PYTHONDONTWRITEBYTECODE='1')
qpath=r/'qualification-qualified/qualification.json';q=json.loads(qpath.read_text())
assert q['status']=='pass' and len(q['checks'])==24
def extra(name,fn):
    try:value=fn();q['checks'].append(dict(name=name,status='pass',result=value))
    except Exception as error:q['checks'].append(dict(name=name,status='fail',failure=repr(error)))
    q['status']='pass' if all(x['status']=='pass' for x in q['checks']) else 'fail'
    qpath.write_text(json.dumps(q,indent=2)+'\n')

driver=s/'docs/plans/a93-a95-e94e95-awf-evidence/qualify-linux.py'
text=driver.read_text()
begin=text.index("resolution=r/'runtime-qualified/RANSWallResolutionProbe'")
end=text.index("execute(['python3','-B',t/'collect_a93_a95_e94e95.py'",begin)
exec(compile(text[begin:end],str(driver)+':legacy-replay-block','exec'))

receipt=r/'runtime-qualified/source-build.json';record=json.loads(receipt.read_text())
current=identity(r/'runtime-qualified/IncNavierStokesSolver')
record['runtime_unchanged']=all(current[key]==record['solver'][key] for key in ('sha256','libraries'))
record['final_content_status']=subprocess.check_output(['git','-c','core.filemode=false','-C',str(s),'status','--porcelain'],text=True)
record['completion_helper_sha256']=sha(Path(__file__))
record['replay_driver_sha256']=sha(driver)
record['replay_build_provenance_sha256']=sha(r/'replay-build-provenance.json')
record['prototype_receipt_error']='KeyError frozen after all scheduled tests; corrected identity comparison without rerunning or altering trajectories'
record['prototype_log_sha256']=sha(r/'qualified-driver.log')
assert record['runtime_unchanged'] and not record['final_content_status']
receipt.write_text(json.dumps(record,indent=2)+'\n')
command=['python3','-B',str(t/'collect_a93_a95_e94e95.py'),'--matrix',str(r/'matrix-qualified/matrix.json'),
         '--qualification',str(qpath),'--output',str(r/'matrix-qualified/final-results.json')]
subprocess.run(command,check=True)
assert (r/'runtime-qualified/actual-loaded-libraries.json').is_file()
print(json.dumps(dict(status=q['status'],checks=len(q['checks']),runtime_unchanged=record['runtime_unchanged'])),flush=True)
