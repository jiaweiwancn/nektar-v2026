#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze and qualify an already-built isolated Linux experiment; no retries.

Requires ROOT/source and ROOT/build. Run after the documented eight-worker
build. Every destination is new; existing evidence is never overwritten.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

sys.dont_write_bytecode = True
# Set numerical-library limits before any imported test module loads NumPy.
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
                                'BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
a = p.parse_args()
r = a.root.resolve(strict=True)
assert platform.system() == 'Linux'
s = r/'source'
t = s/'solvers/IncNavierStokesSolver/Tests/RANS'
sys.path.insert(0, str(t))
from run_cls_epawf_picard import THREADS, sha, state_hash, readback
from run_cls_epawf_pimple import identity
from continue_cls_epawf_pimple import freeze, header
from run_a93_a95_e94e95_matrix import run_case
from compare_cls_epawf_pimple import same_history
from rans_qualification import coverage,dependencies
dependencies("numpy","scipy")

os.environ.update(THREADS)
os.environ.update(TMPDIR=str(r/'tmp'), PYTHONDONTWRITEBYTECODE='1', PYTHONPATH=str(t))
(r/'tmp').mkdir(exist_ok=True)
git = lambda *args: subprocess.check_output(['git', '-c', 'core.filemode=false', '-C', str(s), *args], text=True).strip()
assert not git('status', '--porcelain'), 'source content must be clean'
solver, solver_identity = freeze(r/'build/IncNavierStokesSolver', r/'runtime-qualified')
probe, probe_identity = freeze(r/'build/RANSNewModelsProbe', r/'runtime-probe-qualified')
os.environ['LD_LIBRARY_PATH'] = str(r/'runtime-qualified/lib')
sources = {}
for directory in ('library/SolverUtils/Turbulence', 'library/UnitTests/SolverUtils', 'solvers/IncNavierStokesSolver'):
    for f in (s/directory).rglob('*'):
        if f.is_file() and (f.suffix in ('.h', '.cpp') or f.name == 'CMakeLists.txt'):
            sources[str(f.relative_to(s))] = sha(f)
record = dict(revision=git('rev-parse', 'HEAD'), content_status=git('status', '--porcelain'),
              compiler=subprocess.check_output(['c++', '--version'], text=True), source_files=sources,
              source_manifest_sha256=hashlib.sha256(json.dumps(sources, sort_keys=True).encode()).hexdigest(),
              cmake_cache_sha256=sha(r/'build/CMakeCache.txt'), solver=solver_identity, probe=probe_identity,
              ranks=1, threads=THREADS, tmpdir=os.environ['TMPDIR'], commands=[], required_failures=[])
receipt = r/'runtime-qualified/source-build.json'
def save(): receipt.write_text(json.dumps(record, indent=2)+'\n')
save()

def execute(command, log, env=None, witness=False):
    command = list(map(str, command))
    with (r/log).open('w') as stream:
        child = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT, env=env)
        if witness:
            # Read only this private native executable. Capture dynamically
            # loaded MPI plugins as well as libraries reported by ldd.
            target = str(solver.resolve())
            found = False
            while child.poll() is None and not found:
                for proc in Path('/proc').iterdir():
                    if not proc.name.isdigit(): continue
                    try:
                        if os.readlink(proc/'exe') != target: continue
                        time.sleep(.5)
                        mapped = sorted({line.split()[-1] for line in (proc/'maps').read_text().splitlines()
                                         if '/' in line and '.so' in line.split()[-1]})
                        data = dict(pid=int(proc.name), executable=target, cwd=os.readlink(proc/'cwd'),
                                    libraries={f:sha(Path(f)) for f in mapped},
                                    environment={k:v for item in (proc/'environ').read_bytes().split(b'\0') if b'=' in item
                                                 for k,v in [item.decode().split('=',1)] if k in THREADS or k in ('LD_LIBRARY_PATH','TMPDIR')})
                        (r/'runtime-qualified/actual-loaded-libraries.json').write_text(json.dumps(data, indent=2)+'\n')
                        found = True
                        break
                    except (FileNotFoundError, ProcessLookupError, PermissionError): pass
                if not found: time.sleep(.25)
        code = child.wait()
    record['commands'].append(dict(command=command, exit_code=code, log=str(r/log)))
    if code:record['required_failures'].append(dict(stage=log,exit_code=code))
    save()
    print(log, code, flush=True)
    return code

execute(['python3','-B',t/'validate_a93_a95_e94e95.py','--probe',probe,'--output',r/'kernel-results-qualified.json'], 'kernel-qualified.log')
execute(['python3','-B',t/'run_a93_a95_e94e95_matrix.py','--solver',solver,'--descriptor-probe',probe,
         '--output',r/'matrix-qualified','--steps','16'], 'matrix-qualified-driver.log', witness=True)
execute(['python3','-B',t/'qualify_a93_a95_e94e95_native.py','--solver',solver,'--descriptor-probe',probe,
         '--matrix',r/'matrix-qualified/matrix.json','--output',r/'qualification-qualified'], 'qualification-qualified-driver.log')

# The gradient-dependent A95 admissibility cache deserves its own comparison.
qpath = r/'qualification-qualified/qualification.json'
q = json.loads(qpath.read_text())
def extra(name, fn):
    if name not in q.setdefault('expected_check_names',[]):q['expected_check_names'].append(name)
    try: value=fn(); q['checks'].append(dict(name=name,status='pass',result=value))
    except Exception as error: q['checks'].append(dict(name=name,status='fail',failure=repr(error)))
    q['coverage']=coverage(q['checks'],set(q['expected_check_names']),key=lambda x:x['name'])
    q['status']='pass' if q['coverage']['status']=='pass' and all(x['status']=='pass' for x in q['checks']) else 'fail'
    qpath.write_text(json.dumps(q,indent=2)+'\n')
def a95_cache():
    whole=json.loads((r/'matrix-qualified/A95-Natural-Accepted/case.evolution.json').read_text())
    folder=r/'qualification-qualified/cache-A95-Accepted'
    row=run_case(solver,folder,'A95','Natural','Accepted',steps=16,cache='ExactStepV1',descriptor_probe=probe)
    q['cases'].append(row)
    assert row['status']=='pass'
    assert same_history(whole,json.loads((folder/'case.evolution.json').read_text()))
    return dict(bitwise_equal=True, model='A95', wall='Accepted', steps=16)
extra('A95-current-gradient-cache-on-off',a95_cache)

execute(['python3','-B',t/'capture_a93_baseline.py','--source',s,'--build',r/'build','--output',r/'regression-qualified'], 'regression-qualified-driver.log')
execute(['ctest','--test-dir',r/'build','--parallel','1','--output-on-failure','-R',
         '^(RANSStressTangentOracle|RANSHistoryRejectionChecks|RANSStandardTransportSmoke|RANSCLSTransportSmoke|RANSAWFOriginalSmoke|RANSAWFEquilibriumSmoke|RANSStandardAWFEquilibriumSmoke)$'], 'legacy-regression-ctest.log')

def legacy_restart():
    base_solver=r/'baseline/runtime/IncNavierStokesSolver'
    base_env=dict(os.environ,LD_LIBRARY_PATH=str(base_solver.parent/'lib'))
    old=r/'cls-reviewed-whole4';new=r/'cls-reviewed-restart2'
    flags=['--dt','.000625','--diagnostic','--wall-policy','UniformEPFlat','--checkpoint','--checkpoint-interval','1']
    assert execute(['python3','-B',t/'run_cls_epawf_pimple.py','--solver',base_solver,'--output',old,'--steps','4',*flags], 'cls-reviewed-whole4-driver.log',env=base_env)==0
    cp=next(x for x in old.glob('checkpoint*') if x.is_dir() and (x/'State.bin').exists() and header(x)['steps']==2)
    assert execute(['python3','-B',t/'run_cls_epawf_pimple.py','--solver',solver,'--output',new,'--steps','2','--restart',cp,*flags], 'cls-reviewed-restart2-driver.log')==0
    payload=[]
    for folder in (old,new):
        cp=next(x for x in folder.glob('checkpoint*') if x.is_dir() and (x/'State.bin').exists() and header(x)['steps']==4)
        h=header(cp,True)
        payload.append(dict(state=state_hash(h['coupled_state']),history=[state_hash(v) for v in h['history']],times=h['history_times'],time=h['time'],steps=h['steps']))
    assert payload[0]==payload[1], 'reviewed CLS continuation changed native binary64 history'
    return dict(bitwise_equal=True, reviewed_binary_source='99112a1cbb9c5448dd3de69f50600bf91465062b',
                reviewed_runtime=identity(base_solver), payload=payload[0])
extra('reviewed-CLS-checkpoint-continuation',legacy_restart)

# This separate replay tool retains its legacy CLS scope. Native new-model
# residuals above use the actual model-aware forward audit implementation.
resolution=r/'runtime-qualified/RANSWallResolutionProbe'
shutil.copy2(r/'build/RANSWallResolutionProbe',resolution)
(r/'resolution-runtime-identity.json').write_text(json.dumps(identity(resolution),indent=2)+'\n')
for label in ('cls-reviewed-whole4','cls-reviewed-restart2'):
    def resolution_check(label=label):
        folder=r/label;report=json.loads((folder/'case.evolution.json').read_text());manifest={}
        from types import SimpleNamespace
        norms=readback(SimpleNamespace(resolution_probe=resolution,timeout=120),folder,report['state'],report['last_stage_history'],report['last_stage_time'],report['last_stage_alpha'],manifest)
        error=max(abs(x-y) for x,y in zip(norms,report['stage_audits'][-1][6:]))
        assert error<1e-8, 'independent native CLS residual replay differs'
        return dict(maximum_absolute_norm_difference=error,receipt_sha256=sha(folder/'readback.json'))
    extra('legacy-resolution-'+label,resolution_check)
for wall in ('EP','Accepted'):
    def resolution_guard(wall=wall):
        folder=r/'qualification-qualified'/('guard-resolution-'+wall);folder.mkdir(exist_ok=False)
        shutil.copy2(r/'matrix-qualified'/('E94E95-Natural-'+wall)/'case.xml',folder/'case.xml')
        (folder/'input.json').write_text('{}\n')
        with (folder/'solver.log').open('w') as stream:
            status=subprocess.run([str(resolution),'case.xml','input.json','output.json'],cwd=folder,stdout=stream,stderr=subprocess.STDOUT,timeout=120)
        assert status.returncode!=0 and 'standalone CLS resolution replay does not support' in (folder/'solver.log').read_text()
        return dict(rejected=True)
    extra('unsupported-standalone-resolution-'+wall,resolution_guard)

execute(['python3','-B',t/'collect_a93_a95_e94e95.py','--matrix',r/'matrix-qualified/matrix.json','--qualification',qpath,
         '--output',r/'matrix-qualified/compact-results.json'], 'collect-qualified.log')
current_identity=identity(solver)
record['runtime_unchanged'] = all(current_identity[key]==solver_identity[key] for key in ('sha256','libraries'))
record['final_content_status']=git('status','--porcelain')
save()
for label,path,key,wanted in (
    ('kernel',r/'kernel-results-qualified.json','status','pass'),
    ('matrix',r/'matrix-qualified/matrix.json','status','pass'),
    ('continuation',qpath,'status','pass'),
    ('regression',r/'regression-qualified/baseline.json','status','pass'),
    ('collection',r/'matrix-qualified/compact-results.json','bounded_matrix_pass',True)):
    try:
        packet=json.loads(path.read_text())
        if packet.get(key)!=wanted:record['required_failures'].append(dict(stage=label,reason='structured status failed'))
    except Exception as error:record['required_failures'].append(dict(stage=label,reason=repr(error)))
if len(q['checks'])!=28 or q.get('coverage',{}).get('status')!='pass':record['required_failures'].append(dict(stage='coverage',reason='28 required continuation checks missing/failed'))
if not record['runtime_unchanged'] or record['final_content_status']:record['required_failures'].append(dict(stage='provenance',reason='source/runtime changed'))
record['status']='pass' if not record['required_failures'] else 'fail';save()
print('QUALIFICATION_DRIVER_'+record['status'].upper(),flush=True)
sys.exit(0 if record['status']=='pass' else 1)
