from pathlib import Path
import json,subprocess,hashlib
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=['python3','-B',str(r/'tools/run.py'),'--root',str(r)]
def run(label,phase,cap,args):subprocess.run(base+['--label',label,'--phase',phase,'--seconds',str(cap),'--']+list(map(str,args)),check=True)
def bound(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
assert json.loads((r/'operations/htrace-build-v18/receipt.json').read_text())['status']=='PASS'
t=r/'source-v18/solvers/IncNavierStokesSolver/Tests/RANS/h27';table=r/'science/reference-tables-v1.json'
for name,case in [('r2',r/'science/physical-reference-r2-v17/reference-r2.xml'),('h2',r/'science/physical-reference-h2-v1/reference-h2.xml'),('coarse',r/'science/coupled-fixtures-v15/physical-cache.xml')]:
    out=r/'science'/('trace-geometry-'+name+'-v18.json')
    run('htrace-'+name+'-probe-v18','science',40,[r/'runtime-v18/WJH27TraceGeometryProbe',case,out])
    j=json.loads(out.read_text());assert j['horizontal_trace_gate_pass'] is True and j['native_face_gate_pass'] is True
run('htrace-fixtures-v18','analysis',30,['python3','-B',r/'inputs/htrace-fixtures-v18.py'])
o=r/'science/htrace-regression-v18';cases={}
for name in ['poly','hist','physical']:
    case=o/(name+'.xml')
    run('htrace-'+name+'-native-v18','science',120 if name=='physical' else 30,[r/'runtime-v18/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
    report=o/(name+'-current.json')
    if name=='physical':
        run('htrace-physical-current-v18','analysis',150,['python3','-B',t/'physical_transient_oracle.py','--native',o/(name+'.result.json'),'--table',table,'--out',report,'--seconds','140'])
    else:
        run('htrace-'+name+'-current-v18','analysis',35,['python3','-B',t/'coupled_oracle.py','--native',o/(name+'.result.json'),'--table',table,'--out',report,'--seconds','30','--direct-seconds','25'])
    oldname={'poly':'poly-cache','hist':'hist-cache','physical':'physical-cache'}[name]
    cases[name]=dict(old=bound(r/'science/coupled-fixtures-v15'/(oldname+'.result.json')),new=bound(o/(name+'.result.json')),current_equations=bound(report))
source=json.loads((r/'source-v18.json').read_text())['expected_source_identity']
manifest=dict(schema='WJH27HTraceRegressionV1',root=str(o),source_identities=dict(old='47fb3597c4e56348361cfab6eadf9a34ba8d89027e57ed01153ca4d9c5ab10d1',new=source),cases=cases)
p=o/'regression-manifest.json'
with p.open('x')as f:json.dump(manifest,f,indent=2);f.write('\n')
run('htrace-field-regression-v18','analysis',90,['python3','-B',t/'htrace_regression_check.py','--manifest',p,'--table',table,'--out',o/'regression-independent.json','--seconds','85'])
print('HTRACE_REGRESSION_PREREQUISITES_PASSED',flush=True)
case=o/'reference-r2.xml'
run('radial-reference-native-v18','science',2000,[r/'runtime-v18/IncNavierStokesSolver','--use-opt-file',str(case)+'.opt',case])
run('radial-reference-independent-v18','analysis',400,['python3','-B',t/'physical_oracle.py','--native',o/'reference-r2.result.json','--table',table,'--out',o/'independent-r2.json','--seconds','380'])
run('radial-reference-comparison-h2-v18','analysis',450,['python3','-B',t/'physical_oracle.py','--native',r/'science/physical-reference-h2-v1/reference-h2.result.json','--compare-finer',o/'reference-r2.result.json','--table',table,'--out',o/'comparison-h2-r2.json','--seconds','425'])
run('radial-reference-comparison-x2-v18','analysis',600,['python3','-B',t/'physical_oracle.py','--native',r/'science/physical-reference-x2-v1/reference-x2.result.json','--compare-finer',o/'reference-r2.result.json','--table',table,'--out',o/'comparison-x2-r2.json','--seconds','570'])
print('SOURCE18_RADIAL_REFINEMENT_FINISHED',flush=True)
