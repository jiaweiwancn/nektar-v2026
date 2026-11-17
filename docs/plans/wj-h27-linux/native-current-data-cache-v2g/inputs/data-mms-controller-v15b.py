from pathlib import Path
import subprocess,json
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
o=r/'science/coupled-fixtures-v15';t=r/'tools-data-validation-v15/h27';table=r/'science/reference-tables-v1.json'
assert json.loads((o/'poly-cache-comparison.json').read_text())['status']=='PASS_IMPLICIT_ALGORITHM_COMPARISON_ONLY'
def oracle(name,label,out):
    run(label,'analysis',30,['python3','-B',t/'coupled_oracle.py','--native',o/(name+'.result.json'),'--table',table,
        '--out',o/out,'--seconds','25','--direct-seconds','20'])
oracle('poly-cache','poly-cache-current-v15b','poly-cache-current-v2.json')
def native(name,cap):
    case=o/(name+'.xml');run('coupled-'+name+'-native-v15','science',cap,[r/'runtime-v15/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
for policy in ['fresh','cache']:native('hist-'+policy,20)
run('hist-cache-comparison-v15','analysis',30,['python3','-B',t/'implicit_algorithm_check.py','--comparison','NewtonCacheV1',
    '--new-preconditioner','CurrentDataDriftV1','--old',o/'hist-fresh.result.json','--new',o/'hist-cache.result.json',
    '--table',table,'--out',o/'hist-cache-comparison.json','--seconds','25'])
oracle('hist-cache','hist-cache-current-v15','hist-cache-current.json')
native('poly-driftzero',10)
run('data-driftzero-event-v15','analysis',20,['python3','-B',r/'inputs/data-zero-events-v15.py'])
oracle('poly-driftzero','poly-driftzero-current-v15','poly-driftzero-current.json')
print('SOURCE15_CURRENT_DATA_MMS_PHASE_FINISHED',flush=True)
