from pathlib import Path
import os,sys,subprocess,time,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
cmd=['python3','-B',str(t/'audit_rans_response.py'),'--root',str(r),'--admission',str(r/'admissions/A95-N-Accepted-loader-v2.json'),'--cohort-prefix','A95-N-Accepted-response-2e-05-dt-0.002','--time','.016','--receipt-tag','complete-supplement-v2','--supplement-from',str(r/'response-audits/A95-N-Accepted-response-2e-05-dt-0.002-0.016/results.json')]
with (r/'A95-Accepted-complete-response-audit.log').open('x') as f:code=subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
(r/'A95-Accepted-complete-response-audit-process.json').write_text(json.dumps(dict(exit_code=code,command=cmd,scope='new independent-check supplement preserves the original sealed failure and original trajectories'),indent=2)+'\n')
