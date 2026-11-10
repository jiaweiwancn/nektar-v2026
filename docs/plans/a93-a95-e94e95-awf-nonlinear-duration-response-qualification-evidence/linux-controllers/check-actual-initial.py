from pathlib import Path
import os,subprocess,json,sys,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import RequestManifest,sha
from continue_cls_epawf_pimple import header
from rans_qualification import exact_vector_hash
label='A95-actual-initial-binding';mesh=r/'meshes/A95-N-Accepted/Balanced5V1.json';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));checks=[]
for resume in (False,True):
 cmd=['python3','-B',str(t/'run_rans_neutral_pilot.py'),'--root',str(r),'--label',label,'--model','A95','--wall','Accepted','--mesh-file',str(mesh),'--dt','.002','--velocity-amplitude','.00002','--milestones','.016','.1','--initialize-only']+(['--resume'] if resume else [])
 with (r/('actual-initial-'+str(resume)+'.log')).open('x') as f:
  p=None
  while p is None:
   p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
   if p is None:time.sleep(5)
  code=p.wait()
 assert code==0,code
 folder=r/'campaigns'/label;req=json.loads((folder/'request-manifest.json').read_text())['request'];actual=req['actual_initial'];h=header(folder/'actual-initial',True);assert h['steps']==0 and h['time']==0 and h['dt']==.002;assert exact_vector_hash(h['coupled_state'])==actual['state_sha256']!=req['raw_neutral_sha256'];assert actual['receipt']['raw_effective_max_difference']==0;assert RequestManifest(folder,req,True).completed(folder/'actual-initial')
 checks.append(dict(resume=resume,exit_code=code,actual_state_sha256=actual['state_sha256'],neutral_reference_sha256=req['reference_sha256'],native_history_sha256=actual['history_sha256'],performed_steps=0))
(r/'actual-initial-regression.json').write_text(json.dumps(dict(status='pass',checks=checks,scope='actual perturbed coefficients and native zero-step history; identical exact-request resume; no time advancement'),indent=2)+'\n')
