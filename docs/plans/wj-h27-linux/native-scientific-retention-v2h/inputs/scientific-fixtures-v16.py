from pathlib import Path
import shutil,xml.etree.ElementTree as ET
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'tools-scientific-retention-v16'
shutil.copytree(r/'tools-data-offset-v15',t,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/scientific_accumulator-v16.py',t/'h27/scientific_accumulator.py')
shutil.copy2(r/'inputs/scientific-retention-admission-v16.md',t/'scientific-retention-admission-v16.md')
atomic(t/'scientific-overlay.json',dict(scope='Source15 executable unchanged; bounded persisted scientific accumulation and owned retention',
    sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv','.md']]))
o=r/'science/scientific-retention-v16';o.mkdir();old=r/'science/coupled-fixtures-v15';cases={}
for name in ['control','epoch1','epoch2','epoch3','reload9','reload12']:
    template='physical-cache'if name=='control'else'physical-prefix'if name=='epoch1'else'physical-reload12'if name.startswith('reload')else'physical-cold'
    xml=ET.parse(old/(template+'.xml'));x=xml.getroot();p=x.find('.//PARAMETERS');s=x.find('.//SOLVERINFO')
    params={item.text.strip().split('=',1)[0].strip():item for item in p};infos={item.get('PROPERTY'):item for item in s}
    def parameter(key,value):
        if key not in params:params[key]=ET.SubElement(p,'P')
        params[key].text=key+' = '+str(value)
    def info(key,value):
        if key not in infos:infos[key]=ET.SubElement(s,'I',PROPERTY=key)
        infos[key].set('VALUE',str(value))
    reload=name.startswith('reload');count=0 if reload else 12 if name=='control'else 4
    for key,value in dict(NumSteps=count,WJScientificSteps=0 if reload else 1,IO_InfoSteps=1,
            IO_CheckSteps=0 if reload else 3,WJTerminalBegin=.0173,WJTerminalEnd=.0537,WJInjectRetryStep=-1).items():parameter(key,value)
    info('WJOutputPrefix',o/name);info('WJTerminalOutput','False'if reload else'True')
    if name in ['epoch2','epoch3']:info('WJRestart',o/('epoch'+str(int(name[-1])-1)+'.final.bin'))
    if reload:info('WJRestart',o/('epoch3.recovery-'+name[len('reload'):]+'.bin'))
    path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(old/(template+'.xml.opt'),str(path)+'.opt')
    cases[name]=dict(session=identity(path),options=identity(str(path)+'.opt'),template=identity(old/(template+'.xml')))
atomic(o/'cases-manifest.json',dict(schema='WJH27ScientificRetentionCasesV1',cases=cases,runtime=identity(r/'science/runtime-v15.json'),
    executable_source=identity(r/'source-v15.json'),overlay=identity(t/'scientific-overlay.json'),window=[.0173,.0537],no_new_build=True))
print('SCIENTIFIC_RETENTION_FIXTURES_READY')
