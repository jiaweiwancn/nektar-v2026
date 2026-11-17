from pathlib import Path
import subprocess,json,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
for label in ('mms-requal-temporal-independent-v19','mms-requal-spatial-independent-v19','mms-cross-poly8-native-v19'):assert json.loads((r/'operations'/label/'receipt.json').read_text())['status']=='PASS'
assert json.loads((r/'operations/mms-cross-poly8-current-v19/receipt.json').read_text())['status']=='FAIL'
assert not any(v['status']=='RESERVED'for v in json.loads((r/'accounting.json').read_text())['entries'].values())
o=r/'science/mms-requalification-v19';t=r/'tools-mms-requalification-v19/h27';table=r/'science/reference-tables-v1.json';source=json.loads((o/'cases.json').read_text())['source_identity'];cases={}
for name in ('poly8','poly12','hist8','hist12'):
 case=o/'crossing'/(name+'.xml');native=o/'crossing'/(name+'.result.json');direct=name.endswith('8')
 if name!='poly8':run('mms-cross-'+name+'-native-v19','science',30,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
 report=o/'crossing'/(name+'-independent-v19b.json');args=['python3','-B',t/'coupled_oracle.py','--native',native,'--table',table,'--out',report,'--seconds','190'if direct else'55','--direct-seconds','160'if direct else'40']
 if direct:args.append('--direct')
 run('mms-cross-'+name+'-current-v19b','analysis',200 if direct else 60,args)
 j=json.loads(report.read_text());assert j['current_coupled_step_gate_pass']is True
 if direct:assert j['direct_reference_gate_pass']is True
 cases[name]=dict(native=bound(native),independent=bound(report))
p=o/'crossing/manifest.json'
with p.open('x')as f:json.dump(dict(schema='WJH27MovingCutoffMMSV1',root=str(o/'crossing'),source_identity=source,geometry=bound(o/'crossing/geometry.json'),cases=cases),f,indent=2);f.write('\n')
run('mms-cross-coverage-independent-v19','analysis',60,['python3','-B',t/'mms_cutoff_crossing_check.py','--manifest',p,'--table',table,'--out',o/'crossing/independent.json','--seconds','55'])
assert json.loads((o/'crossing/independent.json').read_text())['foundation_gate_pass']is True
print('SOURCE19_MMS_REQUALIFICATION_FINISHED',flush=True)
