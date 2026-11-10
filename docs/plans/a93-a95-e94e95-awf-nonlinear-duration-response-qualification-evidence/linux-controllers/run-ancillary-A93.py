from pathlib import Path
import subprocess,os,json,time,sys
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';proof=r/'duration-audits/A93-N-Accepted-duration-0.002-0.5/results.json';assert json.loads(proof.read_text())['status']=='pass';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400');events=[]
def run(name,args):
 cmd=['python3','-B',*map(str,args)]
 with (r/(name+'.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 events.append(dict(name=name,exit_code=code,command=cmd));(r/'ancillary-A93-progress.json').write_text(json.dumps(dict(status='running',events=events),indent=2)+'\n');return code
accepted=r/'admissions/A93-N-Accepted-loader-v2.json';ep=r/'admissions/A93-N-EP-loader-v2.json'
neutral=run('A93-EP-short-neutral',[t/'run_rans_balanced_duration.py','--root',r,'--admission',ep,'--dt','.002','--max-time','.016','--concurrency','2'])
if neutral==0:run('A93-EP-short-response',[t/'run_rans_nonlinear_response.py','--root',r,'--admission',ep,'--accepted-neutral-proof',proof,'--dt','.002','--max-time','.016','--amplitude','2e-5','--concurrency','2'])
run('A93-Accepted-half-response',[t/'run_rans_nonlinear_response.py','--root',r,'--admission',accepted,'--accepted-neutral-proof',proof,'--dt','.002','--max-time','.016','--amplitude','1e-5','--concurrency','2'])
full=r/'response-cohorts/A93-N-Accepted-response-2e-05-dt-0.002/T-0.016';half=r/'response-cohorts/A93-N-Accepted-response-1e-05-dt-0.002/T-0.016';epbody=r/'response-cohorts/A93-N-EP-response-2e-05-dt-0.002/T-0.016/results.json'
while not (full/'results.json').exists():
 status=json.loads((r/'conditional-response-A93.json').read_text())['status']
 if status in ('not_admitted','failed_or_unresolved'):break
 time.sleep(5)
if (full/'results.json').exists() and (half/'results.json').exists() and json.loads((full/'results.json').read_text())['status']=='pass' and json.loads((half/'results.json').read_text())['status']=='pass':
 output=r/'half-amplitude-comparison';output.mkdir();run('A93-half-amplitude-comparison',[t/'compare_rans_half_response.py','--full',full,'--half',half,'--output',output/'results.json'])
if epbody.exists() and json.loads(epbody.read_text())['status']=='pass' and (full/'results.json').exists() and json.loads((full/'results.json').read_text())['status']=='pass':
 # This read-only multi-replay comparison reserves one complete native family.
 sys.path.insert(0,str(t));from rans_worker_budget import spawn
 cmd=['python3','-B',str(t/'compare_rans_wall_responses.py'),'--root',str(r),'--model','A93','--accepted-prefix','A93-N-Accepted-response-2e-05-dt-0.002','--ep-prefix','A93-N-EP-response-2e-05-dt-0.002','--dt','.002','--time','.016','--output',str(r/'wall-response-comparison')]
 with (r/'A93-wall-response-comparison.log').open('x') as f:
  process=None
  while process is None:
   process=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
   if process is None:time.sleep(5)
  code=process.wait()
 events.append(dict(name='A93-wall-response-comparison',exit_code=code,command=cmd))
(r/'ancillary-A93-results.json').write_text(json.dumps(dict(status='pass' if all(x['exit_code']==0 for x in events) else 'failed_or_unresolved',events=events,scope='selected EP only through .016; half amplitude at .016; accepted model neutral .5 admission required'),indent=2)+'\n')
