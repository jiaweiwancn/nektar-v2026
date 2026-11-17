from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args,check=True):return subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=check)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
assert json.loads((r/'operations/cutoff-panels-publication-v19/receipt.json').read_text())['status']=='PASS'
assert not any(v['status']=='RESERVED'for v in json.loads((r/'accounting.json').read_text())['entries'].values())
run('cutoff-recovery-fixtures-v19','analysis',30,['python3','-B',r/'inputs/cutoff-recovery-fixtures-v19.py'])
o=r/'science/cutoff-recovery-v19';t=r/'tools-cutoff-recovery-v19/h27';table=r/'science/reference-tables-v1.json';pre=json.loads((o/'prelaunch.json').read_text())
for name,cap in [('cold9',60),('zero12',20),('prior-source',15)]:
 case=o/(name+'.xml');run('cutoff-recovery-'+name+'-native-v19','science',cap,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case],check=name!='prior-source')
 if name=='prior-source':assert json.loads((r/'operations/cutoff-recovery-prior-source-native-v19/receipt.json').read_text())['status']=='FAIL'
run('cutoff-recovery-cold-current-v19','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',o/'cold9.result.json','--table',table,'--out',o/'cold-current.json','--seconds','140'])
manifest=dict(schema='WJH27RetainedRecoveryContinuationV1',root=str(o),source_identity=pre['source_identity'],control=pre['control'],inputs=pre['inputs'],
 cold=dict(result=bound(o/'cold9.result.json'),session=pre['cases']['cold9']['session']),zero=dict(result=bound(o/'zero12.result.json'),session=pre['cases']['zero12']['session']),
 current_equations=bound(o/'cold-current.json'),negative=dict(input=pre['cases']['prior-source']['input'],output_prefix=pre['cases']['prior-source']['output_prefix'],
  receipt=bound(r/'operations/cutoff-recovery-prior-source-native-v19/receipt.json'),log=bound(r/'operations/cutoff-recovery-prior-source-native-v19/command.log')))
p=o/'manifest.json'
with p.open('x')as f:json.dump(manifest,f,indent=2);f.write('\n')
run('cutoff-recovery-history-independent-v19','analysis',65,['python3','-B',t/'retained_recovery_check.py','--manifest',p,'--table',table,'--out',o/'independent.json','--seconds','60'])
assert json.loads((o/'independent.json').read_text())['foundation_gate_pass']is True
run('cutoff-recovery-publication-v19','evidence',60,['python3','-B',r/'inputs/cutoff-recovery-publication-v19.py'])
print('SOURCE19_RETAINED_RECOVERY_FINISHED',flush=True)
