from pathlib import Path
import os,sys,json,subprocess,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400');checks=[];bodies={}
for name,runtime,face,resume,reason in [('v4-repeat','runtime-budget-v4',16,False,None),('v5-default','runtime-budget-v5',16,False,None),('v5-default','runtime-budget-v5',16,True,None),('v5-face32','runtime-budget-v5',32,False,None),('v4-unsupported-face32','runtime-budget-v4',32,False,'physical_face_order_not_consumed'),('v5-default','runtime-budget-v5',32,True,'immutable request mismatch')]:
 output='test-neutral-A93-B5-0.016-'+name;destination=r/'physical-budgets'/output;before=sha(destination/'packet.json.gz') if (destination/'packet.json.gz').exists() else None
 cmd=['python3','-B',str(t/'audit_rans_physical_scalar_budget.py'),'--root',str(r),'--label','A93-N-Accepted-duration-0.002-B5','--time','.016','--output-label',output,'--runtime-name',runtime,'--face-order',str(face)]+(['--resume'] if resume else [])
 log=r/('face-control-'+name+'-'+str(resume)+'-'+str(face)+'.log')
 with log.open('x') as f:
  p=None
  while p is None:
   p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
   if p is None:time.sleep(5)
  code=p.wait()
 ok=code==0 if reason is None else code!=0 and reason in log.read_text()
 if before:ok=ok and before==sha(destination/'packet.json.gz')
 checks.append(dict(name=name,resume=resume,face=face,exit_code=code,expected_reason=reason,status='pass' if ok else 'fail',command=cmd,log_sha256=sha(log)));assert ok,checks[-1]
 if reason is None:bodies[name]=json.loads((destination/'results.json').read_text())
exact=bodies['v4-repeat']['physical_budget']==bodies['v5-default']['physical_budget'];checks.append(dict(name='default-face16-bitwise-physical-budget',status='pass' if exact else 'fail'))
change=max(abs(v-bodies['v5-face32']['scaled'][eq][key]) for eq,row in bodies['v5-default']['scaled'].items() for key,v in row.items());checks.append(dict(name='independent-face-order-control',maximum_scaled_change=change,limit=1e-9,status='pass' if change<=1e-9 else 'unresolved'))
(r/'physical-budget-face-regression.json').write_text(json.dumps(dict(status='pass' if all(x['status']=='pass' for x in checks) else 'fail',checks=checks,scope='physical budget probe only; exact default runtime bridge and frozen-coefficient face quadrature'),indent=2)+'\n')
