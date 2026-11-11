from pathlib import Path
import os,sys,json,subprocess,time,resource,hashlib,shutil
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';events=[]
def execute(name,cmd):
 t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 with (r/'logs'/(name+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
 after=resource.getrusage(resource.RUSAGE_CHILDREN)
 events.append(dict(name=name,command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
 (r/'directional-science-execution.json').write_text(json.dumps(dict(status='running',events=events),indent=2)+'\n')
 return code
code=execute('native-export-aux-build',['cmake','--build',str(r/'build'),'--parallel','2','--target','RANSPreservationProbe'])
if code==0:
 aux=r/'runtime-export-check-v2';aux.mkdir();(aux/'lib').symlink_to(r/'runtime/lib',target_is_directory=True)
 shutil.copy2(r/'build/RANSPreservationProbe',aux/'RANSPreservationProbe');(aux/'RANSNewModelsProbe').symlink_to(r/'runtime/RANSNewModelsProbe')
 sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
 pinned=json.loads((r/'runtime-identity.json').read_text());changes=[]
 for name,digest in pinned['compiled_source_sha256'].items():
  path=r/'source'/name
  if sha(path)!=digest:changes.append(name)
 assert changes==['library/UnitTests/SolverUtils/RANSTensorValidation.h'],changes
 (r/'export-check-runtime-identity.json').write_text(json.dumps(dict(source_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=r/'source',text=True).strip(),read_only_test_changes=changes,solver_and_governing_objects='unchanged pinned 11ace15 snapshot; scientific runtime unmodified',probe_sha256=sha(aux/'RANSPreservationProbe'),header_sha256=sha(r/'source/library/UnitTests/SolverUtils/RANSTensorValidation.h')),indent=2)+'\n')
 code=execute('native-export-layout-v1',[sys.executable,'-B',str(s/'check_rans_directional_layout.py'),'--runtime',str(aux),'--output',str(r/'directional-export-layout-v1')])
if code==0:
 while not (r/'directional-history-v3/results.json').exists():time.sleep(5)
 if json.loads((r/'directional-history-v3/results.json').read_text())['status']=='pass':
  code=execute('A95-directional-control-v1',[sys.executable,'-B',str(s/'run_rans_directional_resolution.py'),'--root',str(r),'--old-root','/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3','--contract',str(r/'source/docs/plans/a93-a95-e94e95-awf-directional-resolution-qualification-evidence/experiment-contract.json')])
 else:code=1
(r/'directional-science-execution.json').write_text(json.dumps(dict(status='pass' if code==0 else 'failed_or_unresolved',events=events),indent=2)+'\n')
