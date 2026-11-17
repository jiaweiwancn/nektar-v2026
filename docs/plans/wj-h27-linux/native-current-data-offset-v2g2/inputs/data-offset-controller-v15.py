from pathlib import Path
import subprocess,json,hashlib,shutil,xml.etree.ElementTree as ET
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
run('data-offset-fixtures-v15','analysis',20,['python3','-B',r/'inputs/data-offset-fixtures-v15.py'])
o=r/'science/physical-data-cold-offset-v15';old=r/'science/coupled-fixtures-v15';t=r/'tools-data-offset-v15/h27';table=r/'science/reference-tables-v1.json'
for name,cap in [('prefix',80),('cold',150)]:
    if name=='cold':checkpoint=bound(o/'prefix.final.bin')
    case=o/(name+'.xml');run('data-offset-'+name+'-native-v15','science',cap,[r/'runtime-v15/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
report=o/'cold-current-independent.json'
run('data-offset-cold-current-v15','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',o/'cold.result.json','--table',table,'--out',report,'--seconds','140'])
manifest=dict(schema='WJH27ColdCacheOffsetV1',root=str(o),source_identity=json.loads((old/'physical-cache.result.json').read_text())['state']['source_identity'],
    control=dict(result=bound(old/'physical-cache.result.json'),session=bound(old/'physical-cache.xml')),
    prefix=dict(result=bound(o/'prefix.result.json'),session=bound(o/'prefix.xml')),
    cold=dict(result=bound(o/'cold.result.json'),session=bound(o/'cold.xml')),cold_input=checkpoint,current_equation_report=bound(report))
path=o/'manifest.json'
with path.open('x')as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
run('data-offset-history-independent-v15','analysis',120,['python3','-B',t/'cold_cache_offset_check.py','--manifest',path,'--table',table,'--out',o/'independent.json','--seconds','110'])
print('SOURCE15_OFFSET_COLD_PHASE_FINISHED',flush=True)
