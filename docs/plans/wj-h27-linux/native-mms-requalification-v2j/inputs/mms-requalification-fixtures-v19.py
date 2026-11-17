from pathlib import Path
import os,sys,json,shutil,xml.etree.ElementTree as ET
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');base=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS';t=r/'tools-mms-requalification-v19';t.mkdir();o=r/'science/mms-requalification-v19';o.mkdir()
for name in ('h26','h27'):shutil.copytree(base/name,t/name,ignore=shutil.ignore_patterns('__pycache__'))
for name in ('coupled_oracle.py','mms_study.py','mms_cutoff_crossing_check.py'):shutil.copy2(r/'inputs'/('mms-v19-'+name),t/'h27'/name)
rough=t/'cao-tamura-fig3';rough.mkdir();shutil.copy2(base/'cao-tamura-fig3/rough-fig3.csv',rough/'rough-fig3.csv')
shutil.copy2(r/'inputs/mms-requalification-admission-v19.md',t/'admission.md')
atomic(t/'overlay.json',dict(source=identity(r/'source-v19.json'),model='approved v2 unchanged; v4 inactive',files=[identity(p)for p in sorted(t.rglob('*'))if p.is_file()]))
source=json.loads((r/'source-v19.json').read_text())['expected_source_identity'];assert source=='56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf';cases={}
def coupled(xml):
 x=xml.getroot();si=x.find('.//SOLVERINFO');params=x.find('.//PARAMETERS')
 for key,value in dict(WJImplicitSolver='CoupledJFNKV1',WJStepPreconditioner='CurrentDataDriftV1',WJNonlinearTrace='True').items():
  item=si.find("I[@PROPERTY='"+key+"']")
  if item is None:item=ET.SubElement(si,'I',PROPERTY=key)
  item.set('VALUE',value)
 for key,value in dict(WJStepPCRelativeDrift=.02,WJStepPCMaximumUses=8,WJStepNewtonSeconds=30,WJStepNewtonIterations=20,WJStepKrylovDimension=40,WJStepLinearCycles=8,WJStepFDStep=1e-6).items():
  old=next((p for p in params if p.text.split('=')[0].strip()==key),None)
  if old is None:old=ET.SubElement(params,'P')
  old.text=key+' = '+str(value)
 return x
for group,folder,names in [('temporal',r/'science/mms-study-v11-v3',('dt08','dt04','dt02','hfine','qfine','tighter')),('spatial',r/'science/steady-mms-spatial-v1',('h2','h4','h8','h8q','h8tol'))]:
 d=o/group;d.mkdir()
 for name in names:
  old=folder/(name+'.xml');xml=ET.parse(old);x=coupled(xml);si=x.find('.//SOLVERINFO');assert si.find("I[@PROPERTY='WJRestart']")is None and si.find("I[@PROPERTY='WJInitialFields']")is None
  si.find("I[@PROPERTY='WJOutputPrefix']").set('VALUE',str(d/name));path=d/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(str(old)+'.opt',str(path)+'.opt')
  cases[group+'/'+name]=dict(template=identity(old),session=identity(path),options=identity(Path(str(path)+'.opt')))
 if group=='temporal':
  manifest=json.loads((folder/'temporal-manifest.json').read_text())
 else:manifest=json.loads((folder/'manifest.json').read_text());manifest['source_identity']=source
 for row in manifest['runs']:row['native']=str(d/(row['label']+'.result.json'))
 atomic(d/'manifest.json',manifest)
sys.path.insert(0,str(t/'h27'));import numpy as np
import coupled_mms_check as mms
import physical_oracle as po
from coupled_fixtures import generate
table=r/'science/reference-tables-v1.json';target=mms.IndependentMMS(json.loads(table.read_text()),1,'PolynomialScalarsV1');ref=po.Reference(json.loads(table.read_text()))
value=float(target.functions['K'](.5,target.H,.2)[0]);face=float(ref.cutoff(np.asarray([value]))[0]);cross=o/'crossing';cross.mkdir()
atomic(cross/'geometry.json',dict(schema='WJH27FrozenCrossingGeometryV1',cutoff_face=face,choice=dict(x=.5,time=.2,K=value),native_fields_inspected=False,table=identity(table),source=identity(Path(__file__)),axes=dict(x=[0,.5,1],r=[0,face,.15,.325,.5])))
for name,variant,qx,qr in [('poly8','PolynomialScalarsV1',8,12),('poly12','PolynomialScalarsV1',12,18),('hist8','Historical',8,12),('hist12','Historical',12,18)]:
 path=cross/(name+'.xml');generate(path,table,cross/name,nx=2,nr=2,dt=.02,steps=2,quadrature_x=qx,quadrature_r=qr,linear=3e-13,nonlinear=3e-12,backend='FastFull3DV1',mms_body_amplitude=1,mms_variant=variant,nonlinear_trace=True)
 xml=ET.parse(path);x=coupled(xml);changed=0
 for vertex in x.find('GEOMETRY/VERTEX'):
  xx,zz,dd=map(float,vertex.text.split())
  if abs(zz-(.0875+.075))<1e-13:vertex.text=f'{xx:.17g} {(.0875+face):.17g} {dd:.17g}';changed+=1
 assert changed==3;xml.write(path,encoding='unicode')
 cases['crossing/'+name]=dict(session=identity(path),options=identity(Path(str(path)+'.opt')),variant=variant,quadrature=[qx,qr])
atomic(o/'cases.json',dict(schema='WJH27Source19MMSCasesV1',source_identity=source,source=identity(r/'source-v19.json'),runtime=identity(r/'science/runtime-v19.json'),helpers=identity(t/'overlay.json'),cases=cases,TASK3_READY=False))
print('SOURCE19_MMS_FIXTURES_READY')
