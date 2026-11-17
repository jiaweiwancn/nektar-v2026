from pathlib import Path
import json,hashlib,shutil,xml.etree.ElementTree as ET
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/htrace-regression-v18';o.mkdir();cases={}
def identity(p):return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
for name,template in [('poly','poly-cache'),('hist','hist-cache'),('physical','physical-cache')]:
    old=r/'science/coupled-fixtures-v15'/(template+'.xml');xml=ET.parse(old);x=xml.getroot();si=x.find('.//SOLVERINFO');opts={i.get('PROPERTY'):i for i in si}
    assert 'WJRestart' not in opts and 'WJInitialFields' not in opts and opts['WJStepPreconditioner'].get('VALUE')=='CurrentDataDriftV1'
    opts['WJOutputPrefix'].set('VALUE',str(o/name));path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(str(old)+'.opt',str(path)+'.opt')
    cases[name]=dict(template=identity(old),session=identity(path),options=identity(Path(str(path)+'.opt')),old=identity(old.with_suffix('.result.json')))
old=r/'science/physical-reference-r2-v17/reference-r2.xml';xml=ET.parse(old)
xml.getroot().find(".//SOLVERINFO/I[@PROPERTY='WJOutputPrefix']").set('VALUE',str(o/'reference-r2'))
path=o/'reference-r2.xml';xml.write(path,encoding='unicode');shutil.copy2(str(old)+'.opt',str(path)+'.opt')
with(o/'cases-manifest.json').open('x')as f:
    json.dump(dict(schema='WJH27HTraceCasesV1',cases=cases,radial_template=identity(old),radial_case=identity(path),source=identity(r/'source-v18.json'),runtime=identity(r/'science/runtime-v18.json'),TASK3_READY=False),f,indent=2);f.write('\n')
print('HTRACE_FIXTURES_READY')
