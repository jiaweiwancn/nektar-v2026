#!/usr/bin/env python3
"""Recheck archived neutral states at the frozen isolated-identity tolerance."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('source',type=Path);p.add_argument('qualification',type=Path);p.add_argument('output',type=Path)
a=p.parse_args()
sys.path.insert(0,str(a.source/'solvers/IncNavierStokesSolver/Tests/RANS'))
import validate_cls_epawf_m4 as m
labels={'constrained_physical_L2_projection','physical_global_pressure_gauge','native_exact_horizontal_invariance',
        'shared_matching_velocity_trace','full_matching_traction','required_relative_pressure_offset','CLS_reference_pressure_jump',
        'matching_k_and_dependent_epsilon','zero_inner_k_flux','zero_outer_k_flux','strong_ground_no_slip','compatible_top_epsilon_trace_and_flux'}
original=m.near

def isolated(actual,expected,name,scale=1.,tolerance=1e-10):
    return original(actual,expected,name,scale,min(tolerance,1e-10) if name in labels else tolerance)
m.near=isolated
revision=subprocess.check_output(['git','-C',str(a.source),'rev-parse','HEAD'],text=True).strip()
records=[]
for path in sorted(a.qualification.glob('*/case.evolution.json')):
    raw=json.loads(path.read_text())
    if raw['neutral_reference']:
        result=m.check_reference(raw,projection=path.parent.name=='unanchored-coarse')
        records.append(dict(case=path.parent.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),maximum_divergence=result['maximum_divergence']))
if not records:raise RuntimeError('No completed neutral cases')
result=dict(schema=1,verification_status='passed',source_commit=revision,isolated_identity_tolerance=1e-10,
            checks=m.checks,neutral_case_count=len(records),cases=records,
            maximum_isolated_scaled_errors={k:v for k,v in m.errors.items() if k in labels},
            runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='cases'},indent=2),flush=True)
