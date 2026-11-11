from pathlib import Path
import subprocess,sys,json,time,resource
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';events=[]
for name in ('test_rans_directional_scope.py','test_rans_ownership_contract.py','test_rans_budget_accuracy_contract.py'):
 cmd=[sys.executable,'-B',str(s/name)];t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 with (r/'logs'/(name+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
 after=resource.getrusage(resource.RUSAGE_CHILDREN);events.append(dict(command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
(r/'scope-negative-execution.json').write_text(json.dumps(dict(status='pass' if all(e['exit_code']==0 for e in events) else 'fail',events=events),indent=2)+'\n')
