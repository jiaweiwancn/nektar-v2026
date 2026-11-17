from pathlib import Path
import json,subprocess,hashlib,shutil,xml.etree.ElementTree as ET
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
def write(p,j):
 with p.open('x')as f:json.dump(j,f,indent=2);f.write('\n')
assert json.loads((r/'operations/cutoff-panels-build-v19/receipt.json').read_text())['status']=='PASS'
assert not any(v['status']=='RESERVED'for v in json.loads((r/'accounting.json').read_text())['entries'].values())
run('cutoff-panels-fixtures-v19','analysis',30,['python3','-B',r/'inputs/cutoff-panels-fixtures-v19.py'])
o=r/'science/cutoff-panels-v19';t=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27';table=r/'science/reference-tables-v1.json';source=json.loads((r/'source-v19.json').read_text())['expected_source_identity'];oldsource='1c9b1fe7027814da2b575babb055cf988833e69cf11d63cf6eda34e71a2000f4'
run('cutoff-panels-operator-native-v19','science',40,[r/'runtime-v19/WJH27CutoffPanelProbe',r/'science/htrace-regression-v18/reference-r2.xml',table,o/'operator-request.json',o/'operator-native.json'])
p=o/'operator-manifest.json';write(p,dict(schema='WJH27CutoffPanelOperatorV1',root=str(o),source_identity=source,native=bound(o/'operator-native.json'),
 witness=bound(r/'science/radial-quadrature-v18/cutoff-crossing-witness-v1.json'),request=bound(o/'operator-request.json')))
run('cutoff-panels-operator-independent-v19','analysis',30,['python3','-B',t/'cutoff_panel_check.py','--manifest',p,'--table',table,'--out',o/'operator-independent.json'])
assert json.loads((o/'operator-independent.json').read_text())['operator_gate_pass'] is True
cases={}
for name in ('poly','hist','physical'):
 case=o/(name+'.xml');run('cutoff-panels-'+name+'-native-v19','science',120 if name=='physical'else 30,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
 report=o/(name+'-current.json')
 if name=='physical':run('cutoff-panels-physical-current-v19','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',o/(name+'.result.json'),'--table',table,'--out',report,'--seconds','140'])
 else:run('cutoff-panels-'+name+'-current-v19','analysis',35,['python3','-B',t/'coupled_oracle.py','--native',o/(name+'.result.json'),'--table',table,'--out',report,'--seconds','30','--direct-seconds','25'])
 cases[name]=dict(old=bound(r/'science/htrace-regression-v18'/(name+'.result.json')),new=bound(o/(name+'.result.json')),current_equations=bound(report))
p=o/'regression-manifest.json';write(p,dict(schema='WJH27HTraceRegressionV1',repair_scope='MovingCutoffPanelsV1',root=str(o),source_identities=dict(old=oldsource,new=source),cases=cases))
run('cutoff-panels-field-regression-v19','analysis',90,['python3','-B',t/'htrace_regression_check.py','--manifest',p,'--table',table,'--out',o/'regression-independent.json','--seconds','85'])
assert json.loads((o/'regression-independent.json').read_text())['foundation_gate_pass'] is True
print('CUTOFF_PANEL_REGRESSION_PREREQUISITES_PASSED',flush=True)
case=o/'reference-r2.xml';run('cutoff-panels-radial-native-v19','science',900,[r/'runtime-v19/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
for order,cap,seconds in [(32,400,380),(48,900,850)]:
 run('cutoff-panels-radial-order'+str(order)+'-v19','analysis',cap,['python3','-B',t/'physical_oracle.py','--native',o/'reference-r2.result.json','--table',table,'--out',o/('independent'+str(order)+'.json'),'--order',str(order),'--seconds',str(seconds)])
p=o/'q16-panels-comparison-manifest.json';write(p,dict(schema='WJH27SteadyQuadratureComparisonV1',comparison='MovingCutoffPanelsV1',root=str(o),old_source_identity=oldsource,source_identity=source,
 old=bound(r/'science/radial-quadrature-v18/reference-q16.result.json'),new=bound(o/'reference-r2.result.json'),old_quadratureX=16,new_quadratureX=8,
 independent={str(order):bound(o/('independent'+str(order)+'.json'))for order in (32,48)}))
run('cutoff-panels-radial-comparison-v19','analysis',120,['python3','-B',t/'steady_quadrature_compare.py','--manifest',p,'--table',table,'--out',o/'q16-panels-independent.json','--seconds','110'])
assert json.loads((o/'q16-panels-independent.json').read_text())['quadrature_foundation_gate_pass'] is True
print('CUTOFF_PANEL_RADIAL_EQUATIONS_PASSED',flush=True)
run('cutoff-panels-radial-zero-v19','science',60,['python3','-B',r/'inputs/cutoff-panels-zero-v19.py'])
print('SOURCE19_CUTOFF_PANEL_QUALIFICATION_FINISHED',flush=True)
