#!/usr/bin/env python3
"""Audit the completed M4 run and remaining registered Linux regressions."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import socket
import subprocess

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--root',type=Path,default=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised'))
p.add_argument('--regression',default='m4-attempt-013')
p.add_argument('--output',type=Path)
a=p.parse_args();root=a.root;source=root/'source';build=root/'build'
main=Path('/media/jay/新加卷/nektar-v2026')
def command(args,env=None):return subprocess.check_output(list(map(str,args)),text=True,env=env).strip()
def git(path,*args):return command(['git','-C',path,*args])
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
qualification=root/'m4-attempt-012/oracle/validation.json'
q=json.loads(qualification.read_text());revision=git(source,'rev-parse','HEAD')
assert q['verification_status']=='passed' and q['milestone_closure'] and q['source_commit']==revision
assert not git(source,'status','--porcelain') and not git(main,'status','--porcelain')
assert q['source_end_commit']==revision and not q['source_end_status']
for name,digest in q['source_sha256'].items():assert sha(source/name)==digest,name
for name,digest in q['binary_sha256'].items():assert sha(build/name)==digest,name
receipt=json.loads((root/a.regression/'receipt.json').read_text())
assert receipt['status']=='passed' and receipt['source_commit']==revision
assert receipt['source_end_commit']==revision and not receipt['source_end_status']
assert receipt['main_end_commit']==receipt['main_commit'] and not receipt['main_end_status']
log=root/a.regression/'suite.log';logtext=log.read_text()
assert '100% tests passed, 0 tests failed out of 25' in logtext
isolated=json.loads((root/'m4-attempt-012/isolated-interfaces.json').read_text())
assert isolated['verification_status']=='passed' and isolated['source_commit']==revision
assert isolated['isolated_identity_tolerance']==1e-10 and max(isolated['maximum_isolated_scaled_errors'].values())<=1e-10
history=json.loads((root/'m4-attempt-012/history-norms.json').read_text())
assert history['status']=='passed' and history['source_commit']==revision and history['maximum_scaled_history_error']<=1e-8
out=a.output or root/'m4-final-audit';out.mkdir(exist_ok=False)
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
ldd=command(['ldd',build/'IncNavierStokesSolver'],env);(out/'ldd.txt').write_text(ldd+'\n')
assert 'not found' not in ldd
paths={str(Path(name).resolve()) for name in re.findall(r'(?:=>\s+|^\s*)(/[^\s]+)',ldd,re.M)}
libraries={name:sha(name) for name in sorted(paths)}
assert any('geometry-lib/libSpatialDomains' in name for name in paths)
prior=json.loads((source/'docs/plans/cls-epawf-revised-evidence/m3/integration.json').read_text())
assert libraries==prior['dynamic_libraries_sha256'],'Shared library set or content changed from M3'
expected={'docs/plans/cls-epawf-revised-contract-v1.json':'4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a',
          'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex':'d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f',
          'library/SolverUtils/Turbulence/RANSModels.cpp':'86a7154ba3cb27e4d3319aa42d156a228fd1b180f35f9b217159822ec9d66b05',
          'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp':'5baf29819fed6b81ad1b042bc09258a3e0779d15de13076b688bf32f80ddcf7e'}
for name,digest in expected.items():assert sha(source/name)==digest,name
inventory=command(['ctest','--test-dir',build,'-N'])
assert 'Total Tests: 26' in inventory
(out/'test-inventory.txt').write_text(inventory+'\n')
flags={}
for target in ('RANSEquilibriumTimeProbe','RANSEquilibriumEvolutionProbe','IncNavierStokesSolver'):
 path=build/'CMakeFiles'/(target+'.dir')/'flags.make'
 flags[target]=path.read_text();assert '-Werror' in flags[target],target
m3=[]
for path in build.glob('cls-epawf-m3-*/validation.json'):
 data=json.loads(path.read_text())
 if data['source_commit']==revision and data['verification_status']=='passed' and data['milestone_closure']:m3.append((path,data))
assert len(m3)==1
attempts=[]
for directory in sorted(root.glob('m4-attempt-*')):
 path=directory/'receipt.json'
 if path.exists():
  data=json.loads(path.read_text())
  attempts.append(dict(path=str(path),sha256=sha(path),source_commit=data['source_commit'],status=data['status'],
                       commands=[dict(command=c['command'],exit_code=c['exit_code'],seconds=c['seconds']) for c in data['commands']]))
raw_index=[];maximum_stage=0.;maximum_stage_constraint=0.
for path in sorted(qualification.parent.rglob('*')):
 if path.is_file() and (path.suffix in ('.xml','.json','.log','.bin') or path.name=='State.bin.tmp'):
  raw_index.append(dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path)))
  if path.name=='case.evolution.json':
   data=json.loads(path.read_text())
   for stage in data['stage_audits']:
    maximum_stage=max(maximum_stage,stage[4])
    maximum_stage_constraint=max(maximum_stage_constraint,max(stage[6+data['dimension']+2:]))
(out/'raw-artifact-index.json').write_text(json.dumps(raw_index,indent=2)+'\n')
report=dict(schema=1,status='passed',milestone='M4',source_commit=revision,source_status=git(source,'status','--porcelain'),
            completed_utc=datetime.now(timezone.utc).isoformat(),machine=socket.gethostname(),platform=platform.platform(),
            source_version=command([build/'IncNavierStokesSolver','--version'],env),
            connection='Tailscale:100.105.247.110',compiler_workers=16,maximum_test_MPI_ranks=4,
            thread_environment={name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')},
            registered_test_inventory_count=26,CTest_regression_passed=25,full_registered_M4_oracle_passed=True,
            combined_registered_coverage=26,qualification_execution='Full M4 oracle once; remaining 25 registered CTests at the same source and binary, without duplicate M4 solves',
            M4_recorded_checks=q['checks'],M4_checks=q['checks']+1,final_source_guard_passed=True,native_time_checks=q['native_time_probe']['checks'],
            qualification_path=str(qualification),qualification_sha256=sha(qualification),
            regression_receipt_path=str(root/a.regression/'receipt.json'),regression_log_sha256=sha(log),
            isolated_interface_checks=isolated['checks'],isolated_interface_tolerance=1e-10,physical_history_checks=history['checks'],maximum_physical_history_error=history['maximum_scaled_history_error'],
            m3_checks=m3[0][1]['checks'],m3_qualification_path=str(m3[0][0]),
            frozen_source_sha256=expected,binary_sha256=q['binary_sha256'],static_kernel_archive_sha256=sha(build/'libRANSKernels.a'),
            dynamic_libraries_sha256=libraries,dynamic_libraries_unchanged_from_M3=True,runtime_LD_LIBRARY_PATH=env['LD_LIBRARY_PATH'],
            ldd_sha256=sha(out/'ldd.txt'),compile_flags=flags,compiler=command(['c++','--version']).splitlines()[0],
            cmake=command(['cmake','--version']).splitlines()[0],mpi=command(['mpirun','--version']).splitlines()[0],
            unchanged_main_checkout=dict(commit=git(main,'rev-parse','HEAD'),status=git(main,'status','--porcelain'),expected=receipt['main_commit']),
            source_manifest_verified=True,maximum_accepted_stage_residual=maximum_stage,maximum_accepted_stage_constraint=maximum_stage_constraint,raw_artifact_index=str(out/'raw-artifact-index.json'),raw_artifact_index_sha256=sha(out/'raw-artifact-index.json'),raw_artifact_count=len(raw_index),attempts=attempts,collector_sha256=sha(Path(__file__)))
assert report['unchanged_main_checkout']['commit']==report['unchanged_main_checkout']['expected']
(out/'integration.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:report[k] for k in ('status','source_commit','CTest_regression_passed','M4_checks','native_time_checks','m3_checks')},indent=2))
