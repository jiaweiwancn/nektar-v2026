from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
run('scientific-retention-fixtures-v16','analysis',30,['python3','-B',r/'inputs/scientific-fixtures-v16.py'])
o=r/'science/scientific-retention-v16';t=r/'tools-scientific-retention-v16/h27';table=r/'science/reference-tables-v1.json'
def bound(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
def native(name,cap):
    case=o/(name+'.xml');run('scientific-'+name+'-native-v16','science',cap,[r/'runtime-v15/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
native('control',180);epochs=[];previous=None
for i in [1,2,3]:
    name='epoch'+str(i)
    restart=None if i==1 else bound(o/('epoch'+str(i-1)+'.final.bin'))
    resume=None if previous is None else bound(previous)
    native(name,80)
    manifest=dict(schema='WJH27ScientificBatchV1',root=str(o),window=[.0173,.0537],
        run=dict(result=bound(o/(name+'.result.json')),index=bound(o/(name+'.science-index.json')),session=bound(o/(name+'.xml'))),
        resume=resume,restart=restart)
    path=o/(name+'-accumulator-manifest.json')
    with path.open('x')as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
    output=o/('generation'+str(i))
    run('scientific-accumulator-'+str(i)+'-v16','analysis',60,['python3','-B',t/'scientific_accumulator.py','--manifest',path,
        '--table',table,'--output',output,'--seconds','55'])
    previous=Path(str(output)+'.json');epochs.append(dict(name=name,native=manifest['run'],generation=bound(previous),restart=restart))
report=o/'epoch3-current-independent.json'
run('scientific-final-current-v16','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',o/'epoch3.result.json',
    '--table',table,'--out',report,'--seconds','140'])
manifest=dict(schema='WJH27ScientificRetentionV1',root=str(o),window=[.0173,.0537],
    control=dict(result=bound(o/'control.result.json'),index=bound(o/'control.science-index.json'),session=bound(o/'control.xml')),
    epochs=epochs,current_equation_report=bound(report),source_identity=json.loads((o/'control.result.json').read_text())['state']['source_identity'])
with(o/'retention-manifest.json').open('x')as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
print('SCIENTIFIC_RETENTION_DATA_PHASE_FINISHED',flush=True)
