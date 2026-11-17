from pathlib import Path
import subprocess,json,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
assert not any(v['status']=='RESERVED'for v in json.loads((r/'accounting.json').read_text())['entries'].values())
run('mms-requal-fixtures-v19','analysis',30,['python3','-B',r/'inputs/mms-requalification-fixtures-v19.py'])
o=r/'science/mms-requalification-v19';t=r/'tools-mms-requalification-v19/h27';table=r/'science/reference-tables-v1.json';source=json.loads((o/'cases.json').read_text())['source_identity']
for name,cap in [('dt08',30),('dt04',45),('dt02',90),('hfine',180),('qfine',150),('tighter',100)]:
 case=o/'temporal'/(name+'.xml');run('mms-requal-'+name+'-native-v19','science',cap,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
run('mms-requal-temporal-independent-v19','analysis',450,['python3','-B',t/'mms_study.py','--manifest',o/'temporal/manifest.json','--table',table,'--out',o/'temporal/independent.json','--seconds','430'])
assert json.loads((o/'temporal/independent.json').read_text())['study_gate_pass']is True
print('SOURCE19_JFNK_TEMPORAL_STUDY_PASSED',flush=True)
for name,cap in [('h2',30),('h4',40),('h8',100),('h8q',150),('h8tol',120)]:
 case=o/'spatial'/(name+'.xml');run('mms-requal-'+name+'-native-v19','science',cap,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
run('mms-requal-spatial-independent-v19','analysis',450,['python3','-B',t/'steady_mms_spatial.py','--manifest',o/'spatial/manifest.json','--table',table,'--out',o/'spatial/independent.json','--seconds','430'])
assert json.loads((o/'spatial/independent.json').read_text())['spatial_MMS_gate_pass']is True
print('SOURCE19_STEADY_SPATIAL_STUDY_PASSED',flush=True)
cases={}
for name in ('poly8','poly12','hist8','hist12'):
 case=o/'crossing'/(name+'.xml');native=o/'crossing'/(name+'.result.json');report=o/'crossing'/(name+'-independent.json')
 run('mms-cross-'+name+'-native-v19','science',30,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
 run('mms-cross-'+name+'-current-v19','analysis',60,['python3','-B',t/'coupled_oracle.py','--native',native,'--table',table,'--out',report,'--seconds','55','--direct','--direct-seconds','40'])
 j=json.loads(report.read_text());assert j['current_coupled_step_gate_pass']is True
 cases[name]=dict(native=bound(native),independent=bound(report))
p=o/'crossing/manifest.json'
with p.open('x')as f:json.dump(dict(schema='WJH27MovingCutoffMMSV1',root=str(o/'crossing'),source_identity=source,geometry=bound(o/'crossing/geometry.json'),cases=cases),f,indent=2);f.write('\n')
run('mms-cross-coverage-independent-v19','analysis',60,['python3','-B',t/'mms_cutoff_crossing_check.py','--manifest',p,'--table',table,'--out',o/'crossing/independent.json','--seconds','55'])
assert json.loads((o/'crossing/independent.json').read_text())['foundation_gate_pass']is True
print('SOURCE19_MMS_REQUALIFICATION_FINISHED',flush=True)
