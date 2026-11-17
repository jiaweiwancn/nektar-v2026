from pathlib import Path
import json, shutil, xml.etree.ElementTree as ET
from resources import atomic, identity

r=Path('/home/jay/build/wj-h27-coupling-r1-20261117')
assert json.loads((r/'operations/coupled-build-v14/receipt.json').read_text())['status']=='PASS'
t=r/'tools-terminal-validation-v14'
shutil.copytree(r/'tools-cache-validation-v13',t,ignore=shutil.ignore_patterns('__pycache__'))
shutil.copy2(r/'inputs/terminal_output_check-v14.py',t/'h27/terminal_output_check.py')
shutil.copy2(r/'inputs/terminal-admission-v14.md',t/'admission-v14.md')
atomic(t/'terminal-overlay.json',dict(scope='Prospective test-only source14 terminal validation overlay; compiled source14 unchanged',
    sources=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()and p.suffix in ['.py','.csv','.md']]))
o=r/'science/coupled-terminal-v14';o.mkdir()
templates=r/'science/coupled-fixtures-v13';cases={}
for name in ['aligned','initial','failed','scheduled','dense','prefix','cold']:
    physical=name in ['scheduled','dense','prefix','cold'];source=templates/('physical-cache.xml'if physical else'poly-cache.xml')
    xml=ET.parse(source);root=xml.getroot();p=root.find('.//PARAMETERS');s=root.find('.//SOLVERINFO')
    parameters={item.text.strip().split('=',1)[0].strip():item for item in p}
    solver={item.get('PROPERTY'):item for item in s}
    def parameter(key,value):
        if key not in parameters:parameters[key]=ET.SubElement(p,'P')
        parameters[key].text=key+' = '+str(value)
    def info(key,value):
        if key not in solver:solver[key]=ET.SubElement(s,'I',PROPERTY=key)
        solver[key].set('VALUE',str(value))
    count=6 if name in ['prefix','cold'] else 12 if physical else 1 if name=='failed' else 4
    interval=(.0173,.0537)if physical else(.02,.06)if name=='aligned'else(0.,.02)if name=='failed'else(0.,.06)
    for key,value in dict(NumSteps=count,WJScientificSteps=1 if name=='dense' else 3 if physical else 0,
            WJTerminalBegin=interval[0],WJTerminalEnd=interval[1],IO_InfoSteps=1,WJInjectRetryStep=-1).items():parameter(key,value)
    for key,value in dict(WJOutputPrefix=o/name,WJTerminalOutput='True',WJOutputCoefficients='True').items():info(key,value)
    if name=='cold':info('WJRestart',o/'prefix.final.bin')
    if name=='failed':
        parameter('WJStepNewtonIterations',1);parameter('WJCoupledTolerance',1e-14)
    path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(str(source)+'.opt',str(path)+'.opt')
    cases[name]=dict(template=identity(source),session=identity(path),options=identity(str(path)+'.opt'),
        steps=count,dt=.005 if physical else.02,window=list(interval))
atomic(o/'cases-manifest.json',dict(schema='WJH27TerminalCasesV1',source=identity(r/'source-v14.json'),
    runtime=identity(r/'science/runtime-v14.json'),overlay=identity(t/'terminal-overlay.json'),cases=cases,
    admission=identity(t/'admission-v14.md'),scientific_model='FiniteNuIsotropicWallFullWJCollarV2',v4_implemented=False))
print(json.dumps(dict(cases=list(cases),source14_built=True,overlay=identity(t/'terminal-overlay.json'))))
