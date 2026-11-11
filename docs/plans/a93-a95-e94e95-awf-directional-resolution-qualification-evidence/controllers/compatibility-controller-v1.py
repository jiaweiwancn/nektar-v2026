from pathlib import Path
import subprocess,sys,json,time,resource,os
r=Path(__file__).parent;o=Path('/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3');s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS'
commands=[
 [sys.executable,'-B',str(s/'check_rans_directional_compatibility.py'),'--old-runtime',str(o/'runtime'),'--runtime',str(r/'runtime'),'--output',str(r/'isotropic-compatibility-v1')],
 [sys.executable,'-B',str(s/'check_rans_directional_history.py'),'--root',str(r),'--mesh',str(r/'meshes/A95-N-Accepted/Balanced5V1.json'),'--old-root',str(o)],
 [sys.executable,'-B',str(s/'check_rans_native_ownership.py'),'--probe',str(r/'runtime/RANSPreservationProbe'),'--descriptor-probe',str(r/'runtime/RANSNewModelsProbe'),'--directional','--output',str(r/'directional-geometry-v1')],
 ['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R','^(RANSNativeOwnershipGeometry|RANSNativeAdmissibilityRegression|RANSFrozenBoundaryInterfaceRegression|RANSForwardInterfaceTolerance)$']]
events=[]
for i,cmd in enumerate(commands):
 before=resource.getrusage(resource.RUSAGE_CHILDREN);t=time.monotonic()
 with (r/'logs'/('compatibility-v1-'+str(i)+'.log')).open('x') as log:code=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
 after=resource.getrusage(resource.RUSAGE_CHILDREN);events.append(dict(command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
 (r/'compatibility-v1-execution.json').write_text(json.dumps(dict(status='running',events=events),indent=2)+'\n')
(r/'compatibility-v1-execution.json').write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in events) else 'failed_or_unresolved',events=events),indent=2)+'\n')
