from pathlib import Path
import subprocess,os,json,time,resource
r=Path(__file__).parent
commands=[['cmake','-S',str(r/'source/solvers/IncNavierStokesSolver/Tests/RANS'),'-B',str(r/'build'),'-DRANS_BUILD_SOLVER=ON','-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake','-DRANS_PREBUILT_SPATIAL_LIBRARY='+str(r/'dependency/libSpatialDomains.so.5.9.0'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'],['cmake','--build',str(r/'build'),'--parallel','4','--target','IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks','RANSWallResolutionProbe']]
start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN);events=[]
for cmd in commands:
 with (r/'logs/build-v1.log').open('a') as log:a=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 events.append(dict(command=cmd,exit_code=a.returncode))
 if a.returncode:break
after=resource.getrusage(resource.RUSAGE_CHILDREN)
(r/'build-v1-execution.json').write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in events) else 'fail',events=events,wall_seconds=time.monotonic()-start,user_CPU_seconds=after.ru_utime-before.ru_utime,system_CPU_seconds=after.ru_stime-before.ru_stime),indent=2)+'\n')
