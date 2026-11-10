from pathlib import Path
import argparse,subprocess,os,json,time,sys
p=argparse.ArgumentParser();p.add_argument('--model',required=True);p.add_argument('--controller-pid',type=int,required=True);a=p.parse_args();r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import sha
model=a.model;name=model+'-N-Accepted';ad=r/'admissions'/(name+'-loader-v2.json');proof=r/'duration-audits'/(name+'-duration-0.002-0.5')/'results.json';dest=r/('conditional-response-'+model+'.json');events=[]
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
def save(status,**kw):dest.write_text(json.dumps(dict(status=status,model=model,events=events,**kw),indent=2)+'\n')
while not proof.exists():
 try:cmd=Path('/proc/'+str(a.controller_pid)+'/cmdline').read_bytes()
 except FileNotFoundError:cmd=b''
 if b'run_rans_balanced_duration.py' not in cmd or model.encode() not in cmd:
  save('not_admitted',reason='neutral controller stopped before a passing final audit');raise SystemExit(1)
 save('waiting_for_neutral_duration');time.sleep(10)
proof_data=json.loads(proof.read_text())
if proof_data['status']!='pass':save('not_admitted',reason='neutral final audit failed');raise SystemExit(1)
assert proof_data['configuration']['model']==model and proof_data['configuration']['wall']=='Accepted'
cmd=['python3','-B',str(t/'run_rans_nonlinear_response.py'),'--root',str(r),'--admission',str(ad),'--accepted-neutral-proof',str(proof),'--dt','.002','--max-time','.1','--amplitude','2e-5','--concurrency','2']
with (r/(name+'-response-full.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
events.append(dict(name='matched_response',exit_code=code,command=cmd));save('response_finished',neutral_proof_sha256=sha(proof))
cohort=name+'-response-2e-05-dt-0.002';receipt=r/'response-cohorts'/cohort/'results.json';data=json.loads(receipt.read_text()) if receipt.exists() else dict(milestones=[])
for event in data['milestones']:
 T=event['time'];comparison=Path(event['comparison'])
 if not comparison.exists() or json.loads(comparison.read_text()).get('status')!='pass':continue
 for kind,head in [('perturbed',cohort),('baseline',name+'-duration-0.002')]:
  for member in ('B5','B11'):
   for q in (1,2):
    output=cohort+'-'+kind+'-'+member+'-'+str(T)+'-Q'+str(q)
    command=['python3','-B',str(t/'audit_rans_physical_scalar_budget.py'),'--root',str(r),'--label',head+'-'+member,'--time',str(T),'--output-label',output,'--quadrature-index',str(q)]
    with (r/(output+'.log')).open('x') as f:
     process=None
     while process is None:
      process=spawn(r,command,env=env,stdout=f,stderr=subprocess.STDOUT)
      if process is None:time.sleep(5)
     status=process.wait()
    events.append(dict(name=output,exit_code=status,command=command));save('physical_budget_reconstruction')
save('pass' if code==0 and all(x['exit_code']==0 for x in events) else 'failed_or_unresolved',neutral_proof_sha256=sha(proof),scope='neutral admission required before any perturbation advancement; independent physical flux defects reported without a forced conservation identity')
