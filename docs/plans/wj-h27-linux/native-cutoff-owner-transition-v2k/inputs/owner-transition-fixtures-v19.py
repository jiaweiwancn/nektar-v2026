from pathlib import Path
import json,shutil,xml.etree.ElementTree as ET
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/cutoff-owner-transition-v19';o.mkdir();t=r/'tools-cutoff-owner-transition-v19';shutil.copytree(r/'tools-mms-requalification-v19',t)
shutil.copy2(r/'inputs/cutoff_owner_transition_check-v19.py',t/'h27/cutoff_owner_transition_check.py');shutil.copy2(r/'inputs/owner-transition-admission-v19.md',t/'owner-transition-admission.md')
atomic(t/'transition-overlay.json',dict(parent=identity(r/'tools-mms-requalification-v19/overlay.json'),checker=identity(t/'h27/cutoff_owner_transition_check.py'),source=identity(r/'source-v19.json'),files=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()]))
template=r/'science/mms-requalification-v19/crossing/hist8.xml';geometry=r/'science/mms-requalification-v19/crossing/geometry.json';cases={}
for name,count,policy,recovery,science,restart in [('fresh',20,'RebuildEveryNewtonV1',0,1,''),('cached',20,'CurrentDataDriftV1',3,1,''),('prefix',10,'CurrentDataDriftV1',0,1,''),('cold',10,'CurrentDataDriftV1',0,1,str(o/'prefix.final.bin')),('reload15',0,'CurrentDataDriftV1',0,0,str(o/'cached.recovery-15.bin')),('reload18',0,'CurrentDataDriftV1',0,0,str(o/'cached.recovery-18.bin'))]:
 xml=ET.parse(template);x=xml.getroot();si=x.find('.//SOLVERINFO');opts={item.get('PROPERTY'):item for item in si};assert 'WJInitialFields'not in opts and 'WJRestart'not in opts
 opts['WJOutputPrefix'].set('VALUE',str(o/name));opts['WJStepPreconditioner'].set('VALUE',policy)
 if restart:ET.SubElement(si,'I',PROPERTY='WJRestart',VALUE=restart)
 for item in x.findall('.//PARAMETERS/P'):
  key=item.text.split('=')[0].strip()
  if key in('NumSteps','IO_CheckSteps','WJScientificSteps'):item.text=key+' = '+str(dict(NumSteps=count,IO_CheckSteps=recovery,WJScientificSteps=science)[key])
 case=o/(name+'.xml');xml.write(case,encoding='unicode');shutil.copy2(str(template)+'.opt',str(case)+'.opt')
 cases[name]=dict(session=identity(case),options=identity(Path(str(case)+'.opt')),steps=count,policy=policy,recovery_stride=recovery,scientific_stride=science,restart=restart)
atomic(o/'cases.json',dict(schema='WJH27CutoffOwnerTransitionCasesV1',source_identity=json.loads((r/'source-v19.json').read_text())['expected_source_identity'],template=identity(template),geometry=identity(geometry),
 source=identity(r/'source-v19.json'),runtime=identity(r/'science/runtime-v19.json'),helpers=identity(t/'transition-overlay.json'),cases=cases,TASK3_READY=False))
print('OWNER_TRANSITION_FIXTURES_READY')
