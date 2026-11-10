from pathlib import Path
import os,sys,subprocess,time,json
r=Path(__file__).resolve().parent;t=r/'source/solvers/IncNavierStokesSolver/Tests/RANS';sys.path.insert(0,str(t));from rans_worker_budget import spawn
cmd=['python3','-B',str(t/'diagnose_rans_volume_ownership.py'),'--packet',str(r/'evolved-audits/response-A95-N-Accepted-response-2e-05-dt-0.002-0.016-B5/T-0.0160/packet.json.gz'),'--runtime',str(r/'runtime-resolution-v3'),'--model','A95','--wall','Accepted','--output',str(r/'volume-ownership/A95-Accepted-full-B5-0.016')];env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(r/'tmp'),NEKTAR_RANS_PROBE_TIMEOUT='14400')
with (r/'A95-response-volume-ownership.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,cmd,env=env,stdout=f,stderr=subprocess.STDOUT)
  if p is None:time.sleep(5)
 code=p.wait()
(r/'A95-response-volume-ownership-process.json').write_text(json.dumps(dict(exit_code=code,command=cmd,scope='diagnostic counterfactual only; original qualification remains failed'),indent=2)+'\n')
