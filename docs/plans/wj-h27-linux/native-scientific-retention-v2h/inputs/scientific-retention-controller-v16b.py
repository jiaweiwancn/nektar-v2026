from pathlib import Path
import json,subprocess,hashlib,shutil
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/scientific-retention-v16';t=r/'tools-scientific-retention-v16b'
assert not any(v.get('status')=='RESERVED' for v in json.loads((r/'accounting.json').read_text())['entries'].values())
shutil.copytree(r/'tools-scientific-retention-v16',t,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['scientific_accumulator.py','scientific_retention_check.py','scientific_recovery_check.py']:
 shutil.copy2(r/'inputs'/('v16b-'+name),t/'h27'/name)
shutil.copy2(r/'inputs/scientific-retention-admission-v16b.md',t/'scientific-retention-admission-v16b.md')
base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
(t/'scientific-overlay-v2.json').write_text(json.dumps(dict(scope='Source15 unchanged, immutable V2 per-batch metadata and owned scientific retention',files=[bound(p) for p in sorted(t.rglob('*')) if p.is_file() and p.suffix in ['.py','.csv','.md']]),indent=2)+'\n')
a=o/'accumulated-v2';a.mkdir();previous=None
m=json.loads((o/'retention-manifest.json').read_text())
for i in [1,2,3]:
 packet=json.loads((o/('epoch'+str(i)+'-accumulator-manifest.json')).read_text())
 packet['resume']=None if previous is None else bound(previous)
 p=a/('epoch'+str(i)+'-manifest.json')
 with p.open('x') as f:json.dump(packet,f,indent=2);f.write('\n')
 prefix=a/('generation'+str(i))
 run('scientific-accumulator-'+str(i)+'-v16b','analysis',30,['python3','-B',t/'h27/scientific_accumulator.py','--manifest',p,'--table',r/'science/reference-tables-v1.json','--output',prefix,'--seconds','25'])
 previous=Path(str(prefix)+'.json');m['epochs'][i-1]['generation']=bound(previous)
p=o/'retention-manifest-v2.json'
with p.open('x') as f:json.dump(m,f,indent=2);f.write('\n')
assessment=o/'retention-assessment-v2.json'
run('scientific-retention-assessment-v16b','analysis',120,['python3','-B',t/'h27/scientific_retention_check.py','--manifest',p,'--table',r/'science/reference-tables-v1.json','--out',assessment,'--seconds','110','--commit-retention'])
reloads={}
for step in [9,12]:
 name='reload'+str(step);case=o/(name+'.xml');incoming=bound(o/('epoch3.recovery-'+str(step)+'.bin'))
 run('scientific-retained-'+name+'-native-v16b','science',30,[r/'runtime-v15/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
 reloads[str(step)]=dict(result=bound(o/(name+'.result.json')),session=bound(case),input=incoming)
proof=dict(schema='WJH27PostRetentionColdImportsV1',root=str(o),assessment=bound(assessment),control=bound(o/'control.result.json'),final_generation=bound(previous),reloads=reloads)
p=o/'post-retention-manifest-v2.json'
with p.open('x') as f:json.dump(proof,f,indent=2);f.write('\n')
run('scientific-post-retention-check-v16b','analysis',30,['python3','-B',t/'h27/scientific_recovery_check.py','--manifest',p,'--out',o/'post-retention-check-v2.json','--seconds','25'])
print('SCIENTIFIC_RETENTION_COMPLETE',flush=True)
