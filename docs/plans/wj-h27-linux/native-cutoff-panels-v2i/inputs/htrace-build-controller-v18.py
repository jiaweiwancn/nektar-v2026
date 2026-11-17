from pathlib import Path
import subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
run('htrace-materialize-v18','analysis',30,['python3','-B',r/'inputs/htrace-materialize-v18.py'])
t=r/'source-v18/solvers/IncNavierStokesSolver/Tests/RANS/h27'
run('htrace-build-v18','build',360,['python3','-B',t/'build.py','--root',r,'--version','18','--targets','IncNavierStokesSolver,WJH27TraceGeometryProbe'])
print('SOURCE18_HTRACE_BUILD_FINISHED',flush=True)
