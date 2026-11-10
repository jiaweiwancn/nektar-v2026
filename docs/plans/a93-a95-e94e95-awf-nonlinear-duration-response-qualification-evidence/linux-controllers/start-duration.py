from pathlib import Path
import os,subprocess,json,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';old=r.parent/'a93-a95-e94e95-balanced-20261110';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));checks=[]
def run(name,cmd):
 with (r/(name+'.log')).open('x') as f:code=subprocess.run(list(map(str,cmd)),env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 checks.append(dict(name=name,exit_code=code,command=list(map(str,cmd))));(r/'launch-results.json').write_text(json.dumps(dict(status='running',checks=checks),indent=2)+'\n');print(name,code,flush=True);return code
assert run('build-probe-v2',['cmake','--build',r/'build','--parallel','2','--target','RANSPreservationProbe'])==0
assert run('freeze-probe',['python3','-B',t/'freeze_rans_diagnostic_runtime.py','--root',r])==0
assert run('targeted-ctest',['ctest','--test-dir',r/'build','--parallel','1','--output-on-failure','-R','^(RANSDurationCoverageContract|RANSPilotEvidenceContract|RANSQualificationFailureContract)$'])==0
for model in ('A93','A95'):
 tag=model+'-N-Accepted'
 if run('admit-v3-'+tag,['python3','-B',t/'admit_rans_balanced_duration.py','--root',r,'--previous-root',old,'--model',model,'--wall','Accepted']):continue
 cmd=['python3','-B',str(t/'run_rans_balanced_duration.py'),'--root',str(r),'--admission',str(r/'admissions'/(tag+'.json')),'--dt','.002','--max-time','.5','--concurrency','2']
 log=(r/(tag+'-duration.log')).open('x');p=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);print('DURATION',model,p.pid,flush=True)
 checks.append(dict(name=model+'-duration-launched',pid=p.pid,command=cmd))
(r/'launch-results.json').write_text(json.dumps(dict(status='pass' if all(x.get('exit_code',0)==0 for x in checks) else 'failed_preparation',checks=checks),indent=2)+'\n')
