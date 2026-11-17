from pathlib import Path
import json,hashlib,shutil,xml.etree.ElementTree as ET
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/radial-quadrature-v18';o.mkdir()
t=r/'tools-radial-quadrature-v18';shutil.copytree(r/'tools-scientific-retention-v16b',t,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/steady_quadrature_compare-v18.py',t/'h27/steady_quadrature_compare.py')
def identity(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
for name in ['physical_oracle.py','implicit_algorithm_check.py','physical_cache_history_check.py','physical_transient_oracle.py','history_recovery_check.py']:
    assert(t/'h27'/name).read_bytes()==(r/'source-v18/solvers/IncNavierStokesSolver/Tests/RANS/h27'/name).read_bytes(),name
old=r/'science/htrace-regression-v18/reference-r2.xml';seed=r/'science/htrace-regression-v18/reference-r2.result.json';packet=json.loads(seed.read_text())
assert packet['steady_solver']['converged'] and packet['state']['quadratureX']==8
xml=ET.parse(old);x=xml.getroot();opts={i.get('PROPERTY'):i for i in x.findall('.//SOLVERINFO/I')};parameters={p.text.split('=',1)[0].strip():p for p in x.findall('.//PARAMETERS/P')}
opts['WJOutputPrefix'].set('VALUE',str(o/'reference-q16'));opts['WJInitialFields'].set('VALUE',str(seed));opts['WJInitialSourceIdentity'].set('VALUE',packet['state']['source_identity'])
parameters['WJQuadratureX'].text='WJQuadratureX = 16';parameters['WJSteadySeconds'].text='WJSteadySeconds = 1200'
p=o/'reference-q16.xml';xml.write(p,encoding='unicode');shutil.copy2(str(old)+'.opt',str(p)+'.opt')
with(o/'cases-manifest.json').open('x')as f:
    json.dump(dict(schema='WJH27RadialQuadratureCasesV1',session=identity(p),options=identity(Path(str(p)+'.opt')),seed=identity(seed),template=identity(old),source=identity(r/'source-v18.json'),runtime=identity(r/'science/runtime-v18.json'),
        helper_sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv']],initial_fields_are_only_fresh_steady_guess=True,TASK3_READY=False),f,indent=2);f.write('\n')
print('RADIAL_QUADRATURE_FIXTURES_READY')
