import pathlib,os,sys,time,subprocess
r=pathlib.Path(__file__).parent;sys.path.insert(0,str(r/'source/solvers/IncNavierStokesSolver/Tests/RANS'))
from rans_worker_budget import spawn
while not (r/'diagnosis/A93-Natural-Accepted/common/packet.json.gz').exists():time.sleep(5)
with (r/'A93-frozen-operator-readonly.log').open('x') as f:
 p=None
 while p is None:
  p=spawn(r,['python3','-B',str(r/'A93-frozen-operator-readonly.py')],stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ))
  if p is None:time.sleep(3)
 code=p.wait()
(r/'A93-frozen-operator-readonly-exit.txt').write_text(str(code))
