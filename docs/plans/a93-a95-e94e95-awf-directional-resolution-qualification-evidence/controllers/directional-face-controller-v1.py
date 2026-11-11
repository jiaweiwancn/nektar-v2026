import pathlib,sys,time,resource,subprocess,os,json
r=pathlib.Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_qualification import THREADS,write_json
os.environ.update(THREADS,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'))
wait=time.monotonic()
while not (r/'directional-resolution-supplement-v1/completed-artifacts.json').exists():
 if time.monotonic()-wait>1200:raise RuntimeError('sealed primary budget evidence not completed within diagnostic wait')
 time.sleep(1)
waiting=time.monotonic()-wait;start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
command=[sys.executable,'-B',str(s/'diagnose_rans_directional_budget_faces.py'),'--root',str(r),'--cohort','A95-N-Accepted-directional-dt-0.002','--runtime-name','runtime-face-decomposition-v2','--output',str(r/'directional-physical-face-decomposition-v1')]
with (r/'logs/directional-physical-face-decomposition-v1.log').open('x') as f:code=subprocess.run(command,stdout=f,stderr=subprocess.STDOUT).returncode
after=resource.getrusage(resource.RUSAGE_CHILDREN)
write_json(r/'directional-face-execution.json',dict(command=command,exit_code=code,waiting_for_sealed_inputs_seconds=waiting,wall_seconds=time.monotonic()-start,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
raise SystemExit(code)
