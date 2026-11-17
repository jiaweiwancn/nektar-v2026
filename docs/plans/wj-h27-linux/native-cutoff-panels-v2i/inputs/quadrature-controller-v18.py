from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
assert json.loads((r/'operations/radial-reference-comparison-x2-v18/receipt.json').read_text())['status']=='PASS'
assert not any(value['status']=='RESERVED'for value in json.loads((r/'accounting.json').read_text())['entries'].values())
run('radial-quadrature-fixtures-v18','analysis',30,['python3','-B',r/'inputs/quadrature-fixtures-v18.py'])
o=r/'science/radial-quadrature-v18';t=r/'tools-radial-quadrature-v18/h27';table=r/'science/reference-tables-v1.json';old=r/'science/htrace-regression-v18/reference-r2.result.json'
run('radial-q8-order32-v18','analysis',400,['python3','-B',t/'physical_oracle.py','--native',old,'--table',table,'--out',o/'q8-independent32.json','--order','32','--seconds','380'])
case=o/'reference-q16.xml'
run('radial-q16-native-v18','science',1400,[r/'runtime-v18/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
for order,cap,seconds in [(32,400,380),(48,900,850)]:
    run('radial-q16-order'+str(order)+'-v18','analysis',cap,['python3','-B',t/'physical_oracle.py','--native',o/'reference-q16.result.json','--table',table,'--out',o/('q16-independent'+str(order)+'.json'),'--order',str(order),'--seconds',str(seconds)])
new=o/'reference-q16.result.json';source=json.loads((r/'source-v18.json').read_text())['expected_source_identity']
manifest=dict(schema='WJH27SteadyQuadratureComparisonV1',root=str(o),source_identity=source,old=bound(old),new=bound(new),old_quadratureX=8,new_quadratureX=16,
    independent={str(order):bound(o/('q16-independent'+str(order)+'.json'))for order in [32,48]})
p=o/'q16-comparison-manifest.json'
with p.open('x')as f:json.dump(manifest,f,indent=2);f.write('\n')
run('radial-q16-comparison-v18','analysis',120,['python3','-B',t/'steady_quadrature_compare.py','--manifest',p,'--table',table,'--out',o/'q8-q16-independent.json','--seconds','110'])
print('RADIAL_QUADRATURE16_FINISHED',flush=True)
