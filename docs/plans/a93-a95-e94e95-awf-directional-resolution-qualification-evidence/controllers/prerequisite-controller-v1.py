from pathlib import Path
import os,sys,subprocess,json,time,resource,hashlib,shutil
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS'
while not (r/'build-v3-execution.json').exists():time.sleep(5)
assert json.loads((r/'build-v3-execution.json').read_text())['status']=='pass'
binaries=['IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks','RANSWallResolutionProbe']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for n in binaries:
 dest=r/'runtime'/n
 assert not dest.exists()
 shutil.copy2(r/'build'/n,dest)
files={str(p.relative_to(r/'source')):sha(p) for p in (r/'source').rglob('*') if p.is_file() and p.suffix in ('.cpp','.h','.hpp') and '.git' not in p.parts}
identity=dict(source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r/'source',text=True).strip(),source_dirty=subprocess.check_output(['git','status','--porcelain'],cwd=r/'source',text=True),compiled_source_sha256=files,binaries={n:sha(r/'runtime'/n) for n in binaries},libraries={str(p.name):sha(p) for p in (r/'runtime/lib').glob('*') if p.is_file()},ldd={n:subprocess.check_output(['ldd',str(r/'runtime'/n)],text=True) for n in binaries},thread_environment={k:v for k,v in os.environ.items() if k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
(r/'runtime-identity.json').write_text(json.dumps(identity,indent=2)+'\n')
commands=[
 ['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R','^(RANSNativeDiagnosticChecks|RANSOwnershipRequiredScopeContract|RANSBudgetAccuracyContract|RANSNonlinearEvidenceContract)$'],
 [sys.executable,'-B',str(s/'check_rans_directional_layout.py'),'--runtime',str(r/'runtime'),'--output',str(r/'directional-layout-v1')],
 [sys.executable,'-B',str(s/'reconcile_rans_frozen_scalar_budget.py'),'--root',str(r),'--ownership-root','/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3','--runtime',str(r/'runtime')]]
events=[]
for i,cmd in enumerate(commands):
 before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.monotonic()
 with (r/'logs'/('prerequisite-v1-'+str(i)+'.log')).open('x') as log:code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
 after=resource.getrusage(resource.RUSAGE_CHILDREN)
 events.append(dict(command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,user_CPU_seconds=after.ru_utime-before.ru_utime,system_CPU_seconds=after.ru_stime-before.ru_stime))
 (r/'prerequisite-v1-execution.json').write_text(json.dumps(dict(status='running',events=events),indent=2)+'\n')
(r/'prerequisite-v1-execution.json').write_text(json.dumps(dict(status='pass' if all(e['exit_code']==0 for e in events) else 'failed_or_unresolved',events=events),indent=2)+'\n')
