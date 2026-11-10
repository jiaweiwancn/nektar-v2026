from pathlib import Path
import subprocess,json,os,sys
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));checks=[]
for name,cmd in [('configure-ctest-current',['cmake','-S',str(t),'-B',str(r/'build')]),('ctest-current',['ctest','--test-dir',str(r/'build'),'-j','1','--output-on-failure','-R','^(RANSDurationCoverageContract|RANSNonlinearEvidenceContract|RANSPilotEvidenceContract|RANSQualificationFailureContract)$'])]:
 with (r/(name+'.log')).open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 checks.append(dict(name=name,command=cmd,exit_code=code));
 if code!=0:break
(r/'ctest-current.json').write_text(json.dumps(dict(status='pass' if len(checks)==2 and all(x['exit_code']==0 for x in checks) else 'fail',checks=checks),indent=2)+'\n')
