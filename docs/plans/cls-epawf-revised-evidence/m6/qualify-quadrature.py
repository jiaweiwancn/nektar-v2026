#!/usr/bin/env python3
"""Supplement M6 with fixed-state-family tangential quadrature refinement.

Run concurrently only with an explicit worker allocation. No compilation or
source checkout is performed; source and executable hashes must stay fixed.
"""
import argparse
from datetime import datetime,timezone
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
os.environ.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
os.environ['PYTHONDONTWRITEBYTECODE']='1';sys.dont_write_bytecode=True
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
p.add_argument('--qualification',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
source=a.root/'source';build=a.root/'build';sys.path.insert(0,str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
import numpy as np
from validate_cls_epawf_m6 import Runner,manufactured_error,stage_audit
from cls_epawf_m6_cases import transition_session
from validate_cls_epawf_m5 import spaces
def git(*args):return subprocess.check_output(['git','-C',str(source),*args],text=True).strip()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
assert not git('status','--porcelain')
revision=git('rev-parse','HEAD');binary=digest(build/'IncNavierStokesSolver')
run=Runner(build,a.output);started=time.monotonic()
receipt={'schema':1,'started_utc':datetime.now(timezone.utc).isoformat(),'source_commit':revision,
         'solver_sha256':binary,'script_sha256':digest(__file__),'status':'running',
         'ranks':1,'threads':1,'allocation':'one numerical worker on CPU 4, alongside at most four main qualification ranks; no compiler job'}
try:
    refined=run.evolution('q24',generator=transition_session,family='flat',degree=3,steps=2,dt=.0025,quadrature=24)
    baseline=json.loads((a.qualification/'mms-flat-p3/case.evolution.json').read_text())
    assert baseline['time']==refined['time'] and baseline['offsets']==refined['offsets']
    stage=stage_audit(refined);error=manufactured_error(refined);base_error=manufactured_error(baseline)
    offsets=baseline['offsets'];differences=[]
    for block,field in enumerate(spaces(baseline)):
        delta=np.asarray(refined['state'][offsets[block]:offsets[block+1]])-baseline['state'][offsets[block]:offsets[block+1]]
        differences.append(float(np.sqrt(max(0.,delta@field.mass()@delta)/field.volume)))
    assert differences[baseline['dimension']+1]<.1*base_error['epsilon']
    assert error['epsilon']<1.01*base_error['epsilon']
    receipt.update(status='passed',quadrature_orders=[12,24],fixed_degree=3,fixed_mesh=True,
                   baseline_error=base_error,refined_error=error,physical_field_rms_differences=differences,stage=stage,
                   epsilon_quadrature_change_over_spatial_error=differences[baseline['dimension']+1]/base_error['epsilon'])
except Exception as error:
    receipt.update(status='failed',error=str(error));raise
finally:
    receipt.update(seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat(),commands=run.commands,
                   source_end_commit=git('rev-parse','HEAD'),source_end_status=git('status','--porcelain'),solver_end_sha256=digest(build/'IncNavierStokesSolver'))
    assert receipt['source_end_commit']==revision and not receipt['source_end_status'] and receipt['solver_end_sha256']==binary
    (a.output/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'seconds':receipt['seconds'],'epsilon_quadrature_change_over_spatial_error':receipt.get('epsilon_quadrature_change_over_spatial_error')}),flush=True)
