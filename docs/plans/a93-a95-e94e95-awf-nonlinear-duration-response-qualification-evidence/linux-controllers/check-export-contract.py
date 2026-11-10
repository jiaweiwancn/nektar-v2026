from pathlib import Path
import subprocess,os,json,time
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));checks=[]
for name in ('test_rans_duration_contract.py','test_rans_pilot_evidence.py','test_rans_qualification.py','test_rans_nonlinear_evidence.py'):
 with (r/(name+'.log')).open('x') as f:code=subprocess.run(['python3','-B',str(t/name)],env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 checks.append(dict(name=name,exit_code=code,status='pass' if code==0 else 'fail'))
code=subprocess.run(['python3','-B',str(t/'check_rans_nonlinear_evidence_failures.py'),'--output',str(r/'complete-export-negative-controls')],env=env).returncode
checks.append(dict(name='complete-export-cli-controls',exit_code=code,status='pass' if code==0 else 'fail'))
(r/'targeted-regressions.json').write_text(json.dumps(dict(status='pass' if all(x['status']=='pass' for x in checks) else 'fail',checks=checks,source_revision=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip()),indent=2)+'\n')
