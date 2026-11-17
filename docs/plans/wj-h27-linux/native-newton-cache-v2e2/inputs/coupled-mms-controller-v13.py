from pathlib import Path
import subprocess,json
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args,check=True):return subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=check).returncode
assert json.loads((r/'operations/coupled-build-v13/receipt.json').read_text())['status']=='PASS'
run('coupled-fixtures-v13','analysis',30,['python3','-B',r/'inputs/coupled-fixtures-v13.py'])
o=r/'science/coupled-fixtures-v13';t=r/'tools-cache-validation-v13/h27';table=r/'science/reference-tables-v1.json'
def native(name,cap):
 case=o/(name+'.xml');run('coupled-'+name+'-native-v13','science',cap,[r/'runtime-v13/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
for family,cap in [('poly',30),('hist',60)]:
 native(family+'-fresh',cap);native(family+'-cache',cap)
 run('coupled-'+family+'-cache-comparison-v13','analysis',30,['python3','-B',t/'implicit_algorithm_check.py','--comparison','NewtonCacheV1','--old',o/(family+'-fresh.result.json'),'--new',o/(family+'-cache.result.json'),'--table',table,'--out',o/(family+'-cache-comparison.json'),'--seconds','25'])
 run('coupled-'+family+'-cache-current-v13','analysis',90 if family=='poly'else 120,['python3','-B',t/'coupled_oracle.py','--native',o/(family+'-cache.result.json'),'--table',table,'--out',o/(family+'-cache-current.json'),'--seconds','80'if family=='poly'else'110','--direct-seconds','60'])
native('poly-driftzero',60)
run('coupled-zero-drift-events-v13','analysis',30,['python3','-B',r/'inputs/coupled-zero-drift-events-v13.py'])
run('coupled-zero-drift-current-v13','analysis',60,['python3','-B',t/'coupled_oracle.py','--native',o/'poly-driftzero.result.json','--table',table,'--out',o/'poly-driftzero-current.json','--seconds','55','--direct-seconds','45'])
print('SOURCE13_MMS_CACHE_PHASE_FINISHED',flush=True)
