import pathlib,subprocess,sys,time,resource,os,json
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_qualification import write_json,THREADS
os.environ.update(THREADS,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'))
command=[sys.executable,'-B',str(s/'collect_rans_directional_resolution.py'),'--root',str(r),'--old-root','/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-ownership-20261111-corrected-v3','--output',str(r/'compact-export-v1')]
start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
with (r/'logs/directional-collection-v1.log').open('x') as f:code=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT).returncode
after=resource.getrusage(resource.RUSAGE_CHILDREN)
write_json(r/'collection-process-v1.json',dict(command=command,exit_code=code,wall_seconds=time.monotonic()-start,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,expected_nonzero_for_failed_scientific_gates=True))
raise SystemExit(code)
