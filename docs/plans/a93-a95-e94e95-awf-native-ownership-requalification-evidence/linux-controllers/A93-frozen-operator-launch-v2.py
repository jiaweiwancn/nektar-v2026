import pathlib,sys,subprocess,os,time
r=pathlib.Path(__file__).parent;sys.path.insert(0,str(r/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
from rans_worker_budget import spawn
with (r/'A93-frozen-operator-readonly-v2.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,['python3','-B',str(r/'A93-frozen-operator-readonly-v2.py')],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
  if p is None:time.sleep(3)
 code=p.wait()
(r/'A93-frozen-operator-readonly-v2-exit.txt').write_text(str(code))
