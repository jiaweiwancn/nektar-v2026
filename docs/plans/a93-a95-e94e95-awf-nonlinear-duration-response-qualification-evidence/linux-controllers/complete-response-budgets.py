from pathlib import Path
import argparse,subprocess,os,sys,json,time
p=argparse.ArgumentParser();p.add_argument('--model',required=True);a=p.parse_args();r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
from rans_request_manifest import RequestManifest,sha
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400');name=a.model+'-N-Accepted';cohort=name+'-response-2e-05-dt-0.002';events=[];dest=r/('complete-response-budgets-'+a.model+'.json')
def save(status,**kw):dest.write_text(json.dumps(dict(status=status,model=a.model,events=events,**kw),indent=2)+'\n')
while True:
 q=r/('conditional-response-'+a.model+'.json');ready=json.loads(q.read_text()) if q.exists() else {}
 if ready.get('status') in ('pass','failed_or_unresolved','not_admitted'):break
 save('waiting_for_response_pipeline');time.sleep(10)
if ready.get('status')=='not_admitted':save('not_admitted',reason='neutral duration gate did not pass');raise SystemExit(1)
receipt=r/'response-cohorts'/cohort/'results.json';body=json.loads(receipt.read_text()) if receipt.exists() else dict(milestones=[])
for milestone in body['milestones']:
 T=milestone['time'];path=Path(milestone['comparison']);c=json.loads(path.read_text()) if path.exists() else {}
 if c.get('status')!='pass':continue
 for kind,head in [('perturbed',cohort),('baseline',name+'-duration-0.002')]:
  for member in ('B5','B11','B5-dt2','B11-dt2'):
   output=cohort+'-'+kind+'-'+member+'-'+str(T)+'-Q2-F32';command=['python3','-B',str(t/'audit_rans_physical_scalar_budget.py'),'--root',str(r),'--label',head+'-'+member,'--time',str(T),'--output-label',output,'--quadrature-index','2','--face-order','32']
   with (r/(output+'.log')).open('x') as f:
    proc=None
    while proc is None:
     proc=spawn(r,command,env=env,stdout=f,stderr=subprocess.STDOUT)
     if proc is None:time.sleep(5)
    code=proc.wait()
   events.append(dict(name=output,exit_code=code,command=command));save('physical_face_and_timestep_controls')
 output=r/'physical-budget-analyses'/(cohort+'-'+str(T));command=['python3','-B',str(t/'analyze_rans_physical_budgets.py'),'--root',str(r),'--cohort',cohort,'--time',str(T),'--output',str(output)]
 with (r/(cohort+'-physical-analysis-'+str(T)+'.log')).open('x') as f:code=subprocess.run(command,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 events.append(dict(name='analysis-'+str(T),exit_code=code,command=command));save('analyzed')
save('pass' if events and all(x['exit_code']==0 for x in events) else 'failed_or_unresolved',scope='complete frozen volume/face/space/timestep coverage; physical budget accuracy classifications remain distinct')
