from pathlib import Path
import shutil,xml.etree.ElementTree as ET
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'tools-data-offset-v15'
shutil.copytree(r/'tools-data-validation-v15b',t,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/cold_cache_offset_check-v15.py',t/'h27/cold_cache_offset_check.py')
shutil.copy2(r/'inputs/data-offset-admission-v15.md',t/'offset-admission-v15.md')
atomic(t/'offset-overlay.json',dict(scope='Additional unscheduled-cache-refresh cold import qualification; no executable/physical changes',
    parent=identity(r/'tools-data-validation-v15b/history-overlay.json'),sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv','.md']]))
o=r/'science/physical-data-cold-offset-v15';o.mkdir();old=r/'science/coupled-fixtures-v15';cases={}
for name,count,template in [('prefix',4,'physical-prefix'),('cold',8,'physical-cold')]:
    xml=ET.parse(old/(template+'.xml'));x=xml.getroot()
    for p in x.findall('.//PARAMETERS/P'):
        if p.text.strip().split('=',1)[0].strip()=='NumSteps':p.text='NumSteps = '+str(count)
    for item in x.findall('.//SOLVERINFO/I'):
        if item.get('PROPERTY')=='WJOutputPrefix':item.set('VALUE',str(o/name))
        if item.get('PROPERTY')=='WJRestart':item.set('VALUE',str(o/'prefix.final.bin'))
    path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(old/(template+'.xml.opt'),str(path)+'.opt')
    cases[name]=dict(session=identity(path),options=identity(str(path)+'.opt'),template=identity(old/(template+'.xml')))
atomic(o/'cases-manifest.json',dict(schema='WJH27ColdOffsetCasesV1',cases=cases,source=identity(r/'source-v15.json'),runtime=identity(r/'science/runtime-v15.json'),overlay=identity(t/'offset-overlay.json')))
print('SOURCE15_OFFSET_FIXTURES_READY')
