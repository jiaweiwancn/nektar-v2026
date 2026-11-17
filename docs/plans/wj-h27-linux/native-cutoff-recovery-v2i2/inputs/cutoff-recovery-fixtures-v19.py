from pathlib import Path
import json,shutil,sys,xml.etree.ElementTree as ET
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/cutoff-recovery-v19';o.mkdir();tools=r/'tools-cutoff-recovery-v19';tools.mkdir()
base=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS'
for name in ('h26','h27'):shutil.copytree(base/name,tools/name,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/retained_recovery_check-v19.py',tools/'h27/retained_recovery_check.py')
shutil.copy2(r/'inputs/cutoff-recovery-admission-v19.md',tools/'admission.md')
atomic(tools/'overlay.json',dict(base=identity(r/'source-v19.json'),checker=identity(tools/'h27/retained_recovery_check.py'),files=[identity(p)for p in sorted(tools.rglob('*'))if p.is_file()]))
control=r/'science/cutoff-panels-v19/physical.result.json';template=control.with_name('physical.xml');packet=json.loads(control.read_text());assert packet['rotated_recovery_count']==2 and len(packet['recovery_seconds'])==4
inputs={str(step):identity(r/'science/cutoff-panels-v19'/('physical.recovery-'+str(step)+'.bin'))for step in (9,12)}
old=r/'science/htrace-regression-v18/physical.final.bin';cases={}
for name,steps,path in [('cold9',3,Path(inputs['9']['path'])),('zero12',0,Path(inputs['12']['path'])),('prior-source',0,old)]:
 xml=ET.parse(template);x=xml.getroot();si=x.find('.//SOLVERINFO');opts={i.get('PROPERTY'):i for i in si};assert 'WJRestart'not in opts and 'WJInitialFields'not in opts
 opts['WJOutputPrefix'].set('VALUE',str(o/name));ET.SubElement(si,'I',PROPERTY='WJRestart',VALUE=str(path))
 for item in x.findall('.//PARAMETERS/P'):
  key=item.text.split('=')[0].strip()
  if key in ('NumSteps','IO_CheckSteps','WJScientificSteps'):item.text=key+' = '+str(dict(NumSteps=steps,IO_CheckSteps=0,WJScientificSteps=3 if name=='cold9'else 0)[key])
 case=o/(name+'.xml');xml.write(case,encoding='unicode');shutil.copy2(str(template)+'.opt',str(case)+'.opt')
 cases[name]=dict(session=identity(case),options=identity(Path(str(case)+'.opt')),input=identity(path),output_prefix=str(o/name))
atomic(o/'prelaunch.json',dict(schema='WJH27RetainedRecoveryInputsV1',source_identity=packet['state']['source_identity'],control=dict(result=identity(control),session=identity(template)),
 inputs=inputs,cases=cases,source=identity(r/'source-v19.json'),runtime=identity(r/'science/runtime-v19.json'),helpers=identity(tools/'overlay.json'),TASK3_READY=False))
print('CUTOFF_RECOVERY_FIXTURES_READY')
