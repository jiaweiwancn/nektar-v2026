from pathlib import Path
import os,subprocess,json,sys,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';old=r.parent/'a93-a95-e94e95-balanced-20261110';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),LD_LIBRARY_PATH=str(r/'runtime/lib'));checks=[]
def run(name,args):
 with (r/(name+'.log')).open('x') as f:
  p=None
  while p is None:
   p=spawn(r,['python3','-B',*map(str,args)],env=env,stdout=f,stderr=subprocess.STDOUT)
   if p is None:time.sleep(5)
  code=p.wait()
 checks.append(dict(name=name,exit_code=code,command=list(map(str,args))));print(name,code,flush=True)
 assert code==0,(name,code)
for member in ('B5-dt2','B11-dt2'):
 label='A93-N-EP-'+member;run('ep-short-observer-'+member,[t/'audit_rans_pilot_samples.py','--root',old,'--label',label,'--through-time','.016','--output-directory',r/'reused-short-observers'/label,'--tmp-directory',r/'tmp'])
run('ep-frozen-quadrature',[t/'check_rans_balanced_state_quadrature.py','--root',r,'--campaign-root',old,'--model','A93','--walls','EP','--output-label','frozen-quadrature-A93-N-EP','--time','.016'])
run('admit-A93-N-EP-loader-v2',[t/'admit_rans_balanced_duration.py','--root',r,'--previous-root',old,'--model','A93','--wall','EP','--output',r/'admissions/A93-N-EP-loader-v2.json'])
(r/'ep-control-preparation.json').write_text(json.dumps(dict(status='pass',checks=checks,scope='read-only admission preparation, no perturbation or duration advancement'),indent=2)+'\n')
