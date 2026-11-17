from pathlib import Path
import subprocess,json,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args,check=True):return subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=check).returncode
o=r/'science/coupled-fixtures-v13';t=r/'tools-cache-validation-v13/h27';table=r/'science/reference-tables-v1.json'
for family in ['poly','hist']:
 assert json.loads((o/(family+'-cache-comparison.json')).read_text())['status']=='PASS_IMPLICIT_ALGORITHM_COMPARISON_ONLY'
 assert json.loads((o/(family+'-cache-current.json')).read_text())['status']=='PASS_CURRENT_COUPLED_STEP_ONLY'
assert json.loads((o/'driftzero-event-independent.json').read_text())['status']=='PASS_NUMERICAL_REFRESH_EVENT_CONTROL_ONLY'
def native(name,cap):
 case=o/(name+'.xml');run('coupled-'+name+'-native-v13','science',cap,[r/'runtime-v13/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
def file(p):return dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
native('physical-fresh',300);native('physical-cache',240)
run('coupled-physical-cache-comparison-v13','analysis',180,['python3','-B',t/'implicit_algorithm_check.py','--comparison','NewtonCacheV1','--old',o/'physical-fresh.result.json','--new',o/'physical-cache.result.json','--table',table,'--out',o/'physical-cache-comparison.json','--seconds','170'])
native('physical-prefix',80);inputs=dict(cold=file(o/'physical-prefix.final.bin'))
native('physical-cold',180);inputs['recovery9']=file(o/'physical-cache.recovery-9.bin');native('physical-recovery9',100)
inputs['reload12']=file(o/'physical-cache.recovery-12.bin');native('physical-reload12',30);native('physical-retry',300)
runs={key:str(o/(name+'.result.json'))for key,name in dict(control='physical-cache',prefix='physical-prefix',cold='physical-cold',recovery9='physical-recovery9',reload12='physical-reload12',retry='physical-retry').items()}
reports={}
for key in ['control','cold','recovery9','retry']:
 target=o/('physical-'+key+'-current.json')
 run('coupled-physical-'+key+'-current-v13','analysis',180,['python3','-B',t/'physical_transient_oracle.py','--native',runs[key],'--table',table,'--out',target,'--seconds','170'])
 reports[key]=file(target)
manifest=o/'physical-history-manifest.json';manifest.write_text(json.dumps(dict(schema='WJH27PhysicalCacheHistoryV1',root=str(o),source_identity=json.loads((o/'physical-cache.result.json').read_text())['state']['source_identity'],runs=runs,sessions={key:str(Path(path).with_name(Path(path).name.replace('.result.json','.xml')))for key,path in runs.items()},input_checkpoints=inputs,current_equation_reports=reports),indent=2)+'\n')
run('physical-cache-history-independent-v13','analysis',300,['python3','-B',t/'physical_cache_history_check.py','--manifest',manifest,'--table',table,'--out',o/'physical-history-independent.json','--seconds','290'],False)
print('SOURCE13_PHYSICAL_CACHE_PHASE_FINISHED',flush=True)
