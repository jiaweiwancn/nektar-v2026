from pathlib import Path
import sys,os,time,json,subprocess
r=Path(__file__).parent;s=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(s))
from rans_worker_budget import spawn
from rans_qualification import THREADS,write_json
from qualify_rans_preservation import load
os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime-resolution-v3/lib'),TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
ad=load(r/'admissions/A93-N-Accepted.json');events=[]
for m in ('B5-dt2','B11-dt2'):
 label='A93-N-Accepted-duration-0.002-'+m;outlabel='A93-N-Accepted-duration-0.002-0.5-'+m
 cmd=['python3','-B',str(s/'audit_rans_balanced_evolved.py'),'--root',str(r),'--analysis-runtime',str(r/'runtime-resolution-v3'),'--campaign-root',str(r),'--label',label,'--output-label',outlabel,'--times','.5','--common-directory',ad['common_directory']]
 start=time.monotonic();p=None
 with (r/('A93-half-spatial-'+m+'.log')).open('x') as f:
  while p is None:
   p=spawn(r,cmd,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
   if p is None:time.sleep(3)
  pid,status,usage=os.wait4(p.pid,0);code=os.waitstatus_to_exitcode(status);p.returncode=code
 events.append(dict(member=m,command=cmd,exit_code=code,wall_seconds=time.monotonic()-start,process_cpu_seconds=usage.ru_utime+usage.ru_stime))
 write_json(r/'A93-half-spatial-execution.json',dict(status='pass' if all(x['exit_code']==0 for x in events) else 'fail',checks=events))
