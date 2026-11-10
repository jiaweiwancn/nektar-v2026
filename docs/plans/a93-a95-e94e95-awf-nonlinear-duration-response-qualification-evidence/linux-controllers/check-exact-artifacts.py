from pathlib import Path
import subprocess,os,json,sys
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_request_manifest import sha
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'));checks=[]
for name in ('test_rans_duration_contract.py','test_rans_pilot_evidence.py','test_rans_qualification.py','test_rans_nonlinear_evidence.py'):
 with (r/(name+'-exact-v2.log')).open('x') as f:code=subprocess.run(['python3','-B',str(t/name)],env=env,stdout=f,stderr=subprocess.STDOUT).returncode
 checks.append(dict(name=name,exit_code=code,status='pass' if code==0 else 'fail'))
files=['test_rans_duration_contract.py','test_rans_pilot_evidence.py','test_rans_qualification.py','test_rans_nonlinear_evidence.py','rans_duration_contract.py','rans_nonlinear_evidence_contract.py','check_rans_duration_failures.py','check_rans_nonlinear_evidence_failures.py']
result=dict(status='pass' if all(x['status']=='pass' for x in checks) else 'fail',checks=checks,python=sys.version,source_file_sha256={n:sha(t/n) for n in files},source_revision=subprocess.check_output(['git','-C',str(r/'source'),'rev-parse','HEAD'],text=True).strip())
with (r/'targeted-regressions-current.json').open('x') as f:f.write(json.dumps(result,indent=2)+'\n')
