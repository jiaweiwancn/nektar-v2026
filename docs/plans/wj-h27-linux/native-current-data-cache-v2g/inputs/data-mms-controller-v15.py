from pathlib import Path
import subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
run('coupled-fixtures-v15','analysis',30,['python3','-B',r/'inputs/data-fixtures-v15.py'])
o=r/'science/coupled-fixtures-v15';t=r/'tools-data-validation-v15/h27';table=r/'science/reference-tables-v1.json'
def native(name,cap):
    case=o/(name+'.xml');run('coupled-'+name+'-native-v15','science',cap,[r/'runtime-v15/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
for family in ['poly','hist']:
    for policy in ['fresh','cache']:native(family+'-'+policy,20)
    run(family+'-cache-comparison-v15','analysis',30,['python3','-B',t/'implicit_algorithm_check.py','--comparison','NewtonCacheV1',
        '--new-preconditioner','CurrentDataDriftV1','--old',o/(family+'-fresh.result.json'),'--new',o/(family+'-cache.result.json'),
        '--table',table,'--out',o/(family+'-cache-comparison.json'),'--seconds','25'])
    run(family+'-cache-current-v15','analysis',30,['python3','-B',t/'coupled_oracle.py','--native',o/(family+'-cache.result.json'),
        '--table',table,'--out',o/(family+'-cache-current.json'),'--seconds','25'])
native('poly-driftzero',10)
run('data-driftzero-event-v15','analysis',20,['python3','-B',r/'inputs/data-zero-events-v15.py'])
run('poly-driftzero-current-v15','analysis',30,['python3','-B',t/'coupled_oracle.py','--native',o/'poly-driftzero.result.json',
    '--table',table,'--out',o/'poly-driftzero-current.json','--seconds','25'])
print('SOURCE15_CURRENT_DATA_MMS_PHASE_FINISHED',flush=True)
