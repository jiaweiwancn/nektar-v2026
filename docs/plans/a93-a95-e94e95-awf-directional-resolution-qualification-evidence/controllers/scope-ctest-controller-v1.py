from pathlib import Path
import subprocess,json,time,resource
r=Path(__file__).parent
commands=[['cmake','-S',str(r/'source/solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(r/'build'),'-DRANS_BUILD_SOLVER=ON','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(r/'dependency/libSpatialDomains.so.5.9.0'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'],['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R','^(RANSDirectionalRequiredScopeContract|RANSOwnershipRequiredScopeContract|RANSBudgetAccuracyContract)$']]
events=[]
for i,cmd in enumerate(commands):
 t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 with (r/'logs'/('scope-ctest-v1-'+str(i)+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
 after=resource.getrusage(resource.RUSAGE_CHILDREN);events.append(dict(command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
 if code:break
(r/'scope-ctest-execution.json').write_text(json.dumps(dict(status='pass' if len(events)==2 and all(e['exit_code']==0 for e in events) else 'fail',events=events),indent=2)+'\n')
