from pathlib import Path
import subprocess,sys,json,time,resource,os
r=Path(__file__).parent;c=r/'a93-redistribution-v2';s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS'
ad=json.loads((c/'admissions/A93-N-Accepted.json').read_text());common=ad['common_directory'];tag='A93-N-Accepted-duration-0.002'
commands=[]
for m in ('B5-dt2','B11-dt2'):
 commands.append([sys.executable,'-B',str(s/'audit_rans_balanced_evolved.py'),'--root',str(c),'--analysis-runtime',str(c/'runtime-resolution-v3'),'--campaign-root',str(c),'--label',tag+'-'+m,'--output-label',tag+'-0.016-'+m,'--times','.016','--common-directory',common])
commands.append([sys.executable,'-B',str(s/'diagnose_rans_duration_residual.py'),'--root',str(c),'--admission',str(c/'admissions/A93-N-Accepted.json'),'--cohort',tag,'--time','.016','--output',str(c/'signed-load-supplement')])
events=[]
for i,cmd in enumerate(commands):
 t=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
 with (r/'logs'/('a93-decomposition-v1-'+str(i)+'.log')).open('x') as f:code=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
 after=resource.getrusage(resource.RUSAGE_CHILDREN);events.append(dict(command=cmd,exit_code=code,wall_seconds=time.monotonic()-t,process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime))
 (r/'a93-decomposition-execution.json').write_text(json.dumps(dict(events=events,status='running'),indent=2)+'\n')
 if code:break
(r/'a93-decomposition-execution.json').write_text(json.dumps(dict(events=events,status='pass' if len(events)==len(commands) and all(e['exit_code']==0 for e in events) else 'failed_or_unresolved'),indent=2)+'\n')
