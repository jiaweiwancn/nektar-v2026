from pathlib import Path
import os,sys,subprocess,json,time,resource
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS'
cmd=[sys.executable,'-B',str(s/'supplement_rans_directional_resolution.py'),'--root',str(r),'--old-root','/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3']
t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with (r/'logs/A95-directional-supplement-v1.log').open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
after=resource.getrusage(resource.RUSAGE_CHILDREN)
(r/'directional-supplement-execution.json').write_text(json.dumps(dict(status='pass' if code==0 else 'failed_or_unresolved',command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime),indent=2)+'\n')
