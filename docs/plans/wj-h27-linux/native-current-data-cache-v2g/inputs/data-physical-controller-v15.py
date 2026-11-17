from pathlib import Path
import subprocess,json,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args,check=True):return subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=check).returncode
o=r/'science/coupled-fixtures-v15';t=r/'tools-data-validation-v15/h27';table=r/'science/reference-tables-v1.json'
for name,status in {'poly-cache-comparison.json':'PASS_IMPLICIT_ALGORITHM_COMPARISON_ONLY','poly-cache-current-v2.json':'PASS_CURRENT_COUPLED_STEP_ONLY',
    'hist-cache-comparison.json':'PASS_IMPLICIT_ALGORITHM_COMPARISON_ONLY','hist-cache-current.json':'PASS_CURRENT_COUPLED_STEP_ONLY',
    'driftzero-event-independent.json':'PASS_NUMERICAL_REFRESH_EVENT_CONTROL_ONLY','poly-driftzero-current.json':'PASS_CURRENT_COUPLED_STEP_ONLY'}.items():
    assert json.loads((o/name).read_text())['status']==status
def native(name,cap):
    case=o/(name+'.xml');run('coupled-'+name+'-native-v15','science',cap,[r/'runtime-v15/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
def bound(path):return dict(path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest())
native('physical-fresh',240);native('physical-operator',180);native('physical-cache',180)
for name,policy in [('physical-operator','OperatorDriftV1'),('physical-cache','CurrentDataDriftV1')]:
    run(name+'-comparison-v15','analysis',60,['python3','-B',t/'implicit_algorithm_check.py','--comparison','NewtonCacheV1','--new-preconditioner',policy,
        '--old',o/'physical-fresh.result.json','--new',o/(name+'.result.json'),'--table',table,'--out',o/(name+'-comparison.json'),'--seconds','55'])
native('physical-prefix',80);inputs=dict(cold=bound(o/'physical-prefix.final.bin'));native('physical-cold',150)
inputs['recovery9']=bound(o/'physical-cache.recovery-9.bin');native('physical-recovery9',80)
inputs['reload12']=bound(o/'physical-cache.recovery-12.bin');native('physical-reload12',30);native('physical-retry',240)
runs={key:str(o/(name+'.result.json'))for key,name in dict(control='physical-cache',prefix='physical-prefix',cold='physical-cold',
    recovery9='physical-recovery9',reload12='physical-reload12',retry='physical-retry').items()};reports={}
for key in ['control','cold','recovery9','retry']:
    target=o/('physical-'+key+'-current.json')
    run('physical-'+key+'-current-v15','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',runs[key],
        '--table',table,'--out',target,'--seconds','140'])
    reports[key]=bound(target)
manifest=dict(schema='WJH27PhysicalCacheHistoryV1',root=str(o),source_identity=json.loads((o/'physical-cache.result.json').read_text())['state']['source_identity'],
    preconditioner_policy='CurrentDataDriftV1',terminal_zero_export=True,runs=runs,
    sessions={key:str(Path(path).with_name(Path(path).name.replace('.result.json','.xml')))for key,path in runs.items()},
    input_checkpoints=inputs,current_equation_reports=reports)
path=o/'physical-history-manifest.json'
with path.open('x')as stream:json.dump(manifest,stream,indent=2);stream.write('\n')
run('physical-data-cache-history-independent-v15','analysis',120,['python3','-B',t/'physical_cache_history_check.py','--manifest',path,
    '--table',table,'--out',o/'physical-history-independent.json','--seconds','110'],False)
print('SOURCE15_CURRENT_DATA_PHYSICAL_PHASE_FINISHED',flush=True)
