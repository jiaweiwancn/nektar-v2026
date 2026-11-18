from pathlib import Path
import json,hashlib,subprocess
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def write(p,j):
 with p.open('x')as f:json.dump(j,f,indent=2);f.write('\n')
assert not any(v['status']=='RESERVED'for v in json.loads((r/'accounting.json').read_text())['entries'].values())
run('owner-transition-fixtures-v19','analysis',30,['python3','-B',r/'inputs/owner-transition-fixtures-v19.py'])
o=r/'science/cutoff-owner-transition-v19';t=r/'tools-cutoff-owner-transition-v19/h27';table=r/'science/reference-tables-v1.json';cases=json.loads((o/'cases.json').read_text());prelaunch={}
for name,cap in [('fresh',80),('cached',80),('prefix',50),('cold',60),('reload15',20),('reload18',20)]:
 case=o/(name+'.xml');incoming=cases['cases'][name]['restart']
 if incoming:
  prelaunch[name]=bound(Path(incoming));write(o/(name+'-prelaunch.json'),dict(input=prelaunch[name],session=cases['cases'][name]['session']))
 run('owner-transition-'+name+'-native-v19','science',cap,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
for name in('fresh','cached','cold'):
 run('owner-transition-'+name+'-current-v19','analysis',80,['python3','-B',t/'coupled_oracle.py','--native',o/(name+'.result.json'),'--table',table,'--out',o/(name+'-current.json'),'--seconds','75','--direct-seconds','40'])
 assert json.loads((o/(name+'-current.json')).read_text())['current_coupled_step_gate_pass']is True
run('owner-transition-cache-comparison-v19','analysis',90,['python3','-B',t/'implicit_algorithm_check.py','--old',o/'fresh.result.json','--new',o/'cached.result.json','--table',table,'--comparison','NewtonCacheV1','--new-preconditioner','CurrentDataDriftV1','--out',o/'cache-comparison.json','--seconds','85'])
assert json.loads((o/'cache-comparison.json').read_text())['comparison_gate_pass']is True
p=o/'manifest.json';write(p,dict(schema='WJH27CutoffOwnerTransitionV1',root=str(o),source_identity=cases['source_identity'],geometry=cases['geometry'],
 runs={name:bound(o/(name+'.result.json'))for name in cases['cases']},cold_input=prelaunch['cold'],reload_inputs={str(step):prelaunch['reload'+str(step)]for step in(15,18)},
 current_equations={name:bound(o/(name+'-current.json'))for name in('fresh','cached','cold')},cache_comparison=bound(o/'cache-comparison.json')))
run('owner-transition-independent-v19','analysis',120,['python3','-B',t/'cutoff_owner_transition_check.py','--manifest',p,'--table',table,'--out',o/'independent.json','--seconds','115'])
assert json.loads((o/'independent.json').read_text())['foundation_gate_pass']is True
print('SOURCE19_OWNER_TRANSITION_QUALIFICATION_FINISHED',flush=True)
