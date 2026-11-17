from pathlib import Path
import subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
run('materialize-v19','analysis',30,['python3','-B',r/'inputs/materialize-v19.py'])
t=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27'
run('cutoff-panels-build-v19','build',360,['python3','-B',t/'build.py','--root',r,'--version','19','--targets','IncNavierStokesSolver,WJH27CutoffPanelProbe'])
print('SOURCE19_BUILD_FINISHED',flush=True)
