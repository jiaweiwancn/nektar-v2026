from pathlib import Path
import json,subprocess,hashlib

r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args,check=True):
    return subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=check).returncode
run('terminal-fixtures-v14','analysis',30,['python3','-B',r/'inputs/terminal-fixtures-v14.py'])
o=r/'science/coupled-terminal-v14';t=r/'tools-terminal-validation-v14/h27';table=r/'science/reference-tables-v1.json'
def native(name,cap,success=True):
    case=o/(name+'.xml');code=run('terminal-'+name+'-native-v14','science',cap,[r/'runtime-v14/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case],success)
    if not success:assert code!=0 and (o/(name+'.failure.json')).exists()and not(o/(name+'.result.json')).exists()
def bound(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
native('aligned',10);native('initial',10);native('failed',10,False)
run('terminal-mms-scheduling-independent-v14','analysis',20,['python3','-B',r/'inputs/terminal-mms-check-v14.py'])
native('scheduled',240);native('dense',240);native('prefix',140)
cold_input=bound(o/'prefix.final.bin');native('cold',140)
report=o/'cold-current-independent.json'
run('terminal-cold-current-v14','analysis',180,['python3','-B',t/'physical_transient_oracle.py','--native',o/'cold.result.json','--table',table,'--out',report,'--seconds','170'])
manifest=dict(schema='WJH27TerminalFoundationV1',root=str(o),window=[.0173,.0537],
    source_identity=json.loads((o/'scheduled.result.json').read_text())['state']['source_identity'],
    runs={name:dict(result=bound(o/(name+'.result.json')),index=bound(o/(name+'.science-index.json')),session=bound(o/(name+'.xml')))
        for name in ['scheduled','dense','prefix','cold','aligned','initial']},cold_input=cold_input,
    segment_ownership=[dict(run='prefix',first_step=0,last_step=6,interval='[0,6]'),dict(run='cold',first_step=6,last_step=12,interval='(6,12]')],
    failed_index=bound(o/'failed.science-index.json'),failed_run=bound(o/'failed.failure.json'),current_equation_report=bound(report),
    early_MMS_report=bound(o/'mms-scheduling-independent.json'))
path=o/'terminal-manifest.json'
with path.open('x')as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
run('terminal-output-independent-v14','analysis',600,['python3','-B',t/'terminal_output_check.py','--manifest',path,'--table',table,'--out',o/'terminal-independent.json','--seconds','580'],False)
print('SOURCE14_TERMINAL_PHASE_FINISHED',flush=True)
