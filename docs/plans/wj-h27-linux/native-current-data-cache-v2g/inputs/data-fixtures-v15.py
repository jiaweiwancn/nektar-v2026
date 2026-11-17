from pathlib import Path
import json,shutil,xml.etree.ElementTree as ET
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117')
assert json.loads((r/'operations/coupled-build-v15/receipt.json').read_text())['status']=='PASS'
t=r/'tools-data-validation-v15';shutil.copytree(r/'tools-terminal-validation-v14b',t,ignore=shutil.ignore_patterns('__pycache__'))
for name in ['implicit_algorithm_check','physical_cache_history_check']:
    shutil.copy2(r/'inputs'/(name+'-v15.py'),t/'h27'/(name+'.py'))
for name in ['admission','implicit-comparison','physical-history']:
    shutil.copy2(r/'inputs'/('data-'+name+'-v15.md'),t/(name+'-v15.md'))
atomic(t/'data-policy-overlay.json',dict(scope='Prospective source15 independent validation; compiled capsule remains unchanged',
    sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv','.md']]))
o=r/'science/coupled-fixtures-v15';o.mkdir();old=r/'science/coupled-fixtures-v13'
names=['poly-fresh','poly-cache','hist-fresh','hist-cache','poly-driftzero','physical-fresh','physical-operator',
       'physical-cache','physical-prefix','physical-cold','physical-recovery9','physical-reload12','physical-retry']
cases={}
for name in names:
    template=old/(('physical-cache'if name=='physical-operator'else name)+'.xml')
    xml=ET.parse(template);x=xml.getroot();p=x.find('.//PARAMETERS');s=x.find('.//SOLVERINFO')
    parameters={item.text.strip().split('=',1)[0].strip():item for item in p};solver={item.get('PROPERTY'):item for item in s}
    def parameter(key,value):
        if key not in parameters:parameters[key]=ET.SubElement(p,'P')
        parameters[key].text=key+' = '+str(value)
    def info(key,value):
        if key not in solver:solver[key]=ET.SubElement(s,'I',PROPERTY=key)
        solver[key].set('VALUE',str(value))
    policy='RebuildEveryNewtonV1'if name.endswith('-fresh')else'OperatorDriftV1'if name=='physical-operator'else'CurrentDataDriftV1'
    info('WJOutputPrefix',o/name);info('WJStepPreconditioner',policy);info('WJTerminalOutput','False')
    restart={'physical-cold':'physical-prefix.final.bin','physical-recovery9':'physical-cache.recovery-9.bin',
             'physical-reload12':'physical-cache.recovery-12.bin'}.get(name)
    if restart:info('WJRestart',o/restart)
    if name=='physical-reload12':
        info('WJTerminalOutput','True');parameter('WJTerminalBegin',.0173);parameter('WJTerminalEnd',.0537)
    path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(str(template)+'.opt',str(path)+'.opt')
    cases[name]=dict(template=identity(template),session=identity(path),options=identity(str(path)+'.opt'),preconditioner=policy)
atomic(o/'cases-manifest.json',dict(schema='WJH27CurrentDataCacheCasesV1',source=identity(r/'source-v15.json'),
    runtime=identity(r/'science/runtime-v15.json'),overlay=identity(t/'data-policy-overlay.json'),cases=cases,
    active_model='FiniteNuIsotropicWallFullWJCollarV2',unapproved_v4_implemented=False))
print('SOURCE15_FIXTURES_READY',len(cases))
