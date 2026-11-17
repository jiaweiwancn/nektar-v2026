from pathlib import Path
import json,shutil,xml.etree.ElementTree as ET
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');old=r/'science/coupled-terminal-v14'
assert json.loads((old/'terminal-independent.json').read_text())['status']=='FAIL'
assert (old/'terminal-failure-detail.json').is_file()
t=r/'tools-terminal-validation-v14b';shutil.copytree(r/'tools-terminal-validation-v14',t,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/terminal_output_check-v14b.py',t/'h27/terminal_output_check.py')
shutil.copy2(r/'inputs/terminal-admission-v14b.md',t/'admission-v14b.md')
atomic(t/'cadence2-overlay.json',dict(scope='Output-cadence2 numerical refinement; unchanged source14 executable and fixed gates',
    parent=identity(r/'tools-terminal-validation-v14/terminal-overlay.json'),sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv','.md']]))
o=r/'science/coupled-terminal-cadence2-v14';o.mkdir();cases={}
source=old/'scheduled.xml'
for name in ['scheduled','prefix','cold']:
    xml=ET.parse(source);x=xml.getroot();p=x.find('.//PARAMETERS');s=x.find('.//SOLVERINFO')
    for item in p:
        key=item.text.strip().split('=',1)[0].strip()
        if key=='WJScientificSteps':item.text=key+' = 2'
        if key=='NumSteps':item.text=key+' = '+str(12 if name=='scheduled'else 6)
    for item in s:
        if item.get('PROPERTY')=='WJOutputPrefix':item.set('VALUE',str(o/name))
    if name=='cold':ET.SubElement(s,'I',PROPERTY='WJRestart',VALUE=str(o/'prefix.final.bin'))
    path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(str(source)+'.opt',str(path)+'.opt')
    cases[name]=dict(session=identity(path),options=identity(str(path)+'.opt'),template=identity(source))
atomic(o/'cases-manifest.json',dict(schema='WJH27TerminalCadenceRefinementV1',scientific_stride=2,source=identity(r/'source-v14.json'),
    runtime=identity(r/'science/runtime-v14.json'),overlay=identity(t/'cadence2-overlay.json'),cases=cases,
    reused_cases={name:identity(old/(name+'.result.json'))for name in ['dense','aligned','initial']},
    unchanged_limits=True,unchanged_physics=True,prior_failed_report=identity(old/'terminal-independent.json')))
print('SOURCE14_CADENCE2_FIXTURES_READY')
