from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args,check=True):return subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=check).returncode
run('terminal-cadence2-fixtures-v14','analysis',20,['python3','-B',r/'inputs/terminal-cadence2-fixtures-v14.py'])
o=r/'science/coupled-terminal-cadence2-v14';old=r/'science/coupled-terminal-v14';t=r/'tools-terminal-validation-v14b/h27';table=r/'science/reference-tables-v1.json'
def native(name,cap):
    case=o/(name+'.xml');run('terminal-cadence2-'+name+'-native-v14','science',cap,[r/'runtime-v14/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
def bound(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
native('scheduled',240);native('prefix',140);cold_input=bound(o/'prefix.final.bin');native('cold',140)
report=o/'cold-current-independent.json'
run('terminal-cadence2-cold-current-v14','analysis',180,['python3','-B',t/'physical_transient_oracle.py','--native',o/'cold.result.json','--table',table,'--out',report,'--seconds','170'])
runs={}
for name in ['scheduled','dense','prefix','cold','aligned','initial']:
    folder=o if name in ['scheduled','prefix','cold']else old
    runs[name]=dict(result=bound(folder/(name+'.result.json')),index=bound(folder/(name+'.science-index.json')),session=bound(folder/(name+'.xml')))
manifest=dict(schema='WJH27TerminalFoundationV1',root=str(o),window=[.0173,.0537],scientific_stride=2,
    source_identity=json.loads((o/'scheduled.result.json').read_text())['state']['source_identity'],runs=runs,cold_input=cold_input,
    segment_ownership=[dict(run='prefix',first_step=0,last_step=6,interval='[0,6]'),dict(run='cold',first_step=6,last_step=12,interval='(6,12]')],
    failed_index=bound(old/'failed.science-index.json'),failed_run=bound(old/'failed.failure.json'),
    current_equation_report=bound(report),early_MMS_report=bound(old/'mms-scheduling-independent.json'),
    prior_failed_report=bound(old/'terminal-independent.json'),prior_failure_detail=bound(old/'terminal-failure-detail.json'))
path=o/'terminal-manifest.json'
with path.open('x')as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
run('terminal-cadence2-independent-v14','analysis',120,['python3','-B',t/'terminal_output_check.py','--manifest',path,'--table',table,'--out',o/'terminal-independent.json','--seconds','110'],False)
print('SOURCE14_CADENCE2_PHASE_FINISHED',flush=True)
