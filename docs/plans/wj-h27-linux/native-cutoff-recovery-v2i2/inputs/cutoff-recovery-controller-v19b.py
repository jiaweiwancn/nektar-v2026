from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
for label,status in [('cutoff-recovery-cold9-native-v19','PASS'),('cutoff-recovery-zero12-native-v19','PASS'),('cutoff-recovery-prior-source-native-v19','FAIL'),('cutoff-recovery-cold-current-v19','FAIL')]:assert json.loads((r/'operations'/label/'receipt.json').read_text())['status']==status
assert not any(v['status']=='RESERVED'for v in json.loads((r/'accounting.json').read_text())['entries'].values())
run('cutoff-recovery-fixture-repair-v19b','analysis',15,['python3','-B',r/'inputs/cutoff-recovery-fixture-repair-v19b.py'])
o=r/'science/cutoff-recovery-v19';t=r/'tools-cutoff-recovery-v19b/h27';table=r/'science/reference-tables-v1.json';pre=json.loads((o/'prelaunch.json').read_text())
run('cutoff-recovery-cold-current-v19b','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',o/'cold9.result.json','--table',table,'--out',o/'cold-current-v19b.json','--seconds','140'])
manifest=dict(schema='WJH27RetainedRecoveryContinuationV1',root=str(o),source_identity=pre['source_identity'],control=pre['control'],inputs=pre['inputs'],
 cold=dict(result=bound(o/'cold9.result.json'),session=pre['cases']['cold9']['session']),zero=dict(result=bound(o/'zero12.result.json'),session=pre['cases']['zero12']['session']),
 current_equations=bound(o/'cold-current-v19b.json'),negative=dict(input=pre['cases']['prior-source']['input'],output_prefix=pre['cases']['prior-source']['output_prefix'],
  receipt=bound(r/'operations/cutoff-recovery-prior-source-native-v19/receipt.json'),log=bound(r/'operations/cutoff-recovery-prior-source-native-v19/command.log')))
p=o/'manifest.json'
with p.open('x')as f:json.dump(manifest,f,indent=2);f.write('\n')
run('cutoff-recovery-history-independent-v19','analysis',65,['python3','-B',t/'retained_recovery_check.py','--manifest',p,'--table',table,'--out',o/'independent.json','--seconds','60'])
assert json.loads((o/'independent.json').read_text())['foundation_gate_pass']is True
run('cutoff-recovery-publication-v19','evidence',60,['python3','-B',r/'inputs/cutoff-recovery-publication-v19b.py'])
print('SOURCE19_RETAINED_RECOVERY_FINISHED',flush=True)
