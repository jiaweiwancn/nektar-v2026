from pathlib import Path
import os,sys,subprocess,time,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
model='A93';results=[];env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',NEKTAR_RANS_PROBE_TIMEOUT='14400',LD_LIBRARY_PATH=str(r/'runtime-resolution-v3/lib'))
for member in ('B5-dt2','B11-dt2'):
 label=model+'-N-Accepted-duration-0.002-'+member;output=model+'-N-Accepted-duration-0.002-0.5-'+member
 cmd=['python3','-B',str(t/'audit_rans_balanced_evolved.py'),'--root',str(r),'--campaign-root',str(r),'--analysis-runtime',str(r/'runtime-resolution-v3'),'--label',label,'--output-label',output,'--times','0.5','--common-directory','/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-balanced-20261110/diagnosis/'+model+'-Natural-Accepted/common','--constitutive-samples']
 with (r/(output+'-diagnostic.log')).open('x') as f:
  p=None
  while p is None:
   p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
   if p is None:time.sleep(5)
  code=p.wait()
 results.append(dict(member=member,exit_code=code,command=cmd))
 (r/(model+'-half-spatial-diagnostic.json')).write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in results) else 'fail',results=results),indent=2)+'\n')
