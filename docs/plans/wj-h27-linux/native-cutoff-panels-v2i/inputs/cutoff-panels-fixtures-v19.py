from pathlib import Path
import os,sys,json,time,shutil,xml.etree.ElementTree as ET
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):os.environ[key]='1'
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'source-v19/solvers/IncNavierStokesSolver/Tests/RANS/h27';sys.path.insert(0,str(t));sys.path.insert(0,str(r/'tools'))
from resources import atomic,identity
import numpy as np
import physical_oracle as po
o=r/'science/cutoff-panels-v19';o.mkdir();cases={}
oldfolder=r/'science/htrace-regression-v18';source=json.loads((r/'source-v19.json').read_text())['expected_source_identity']
for name in ('poly','hist','physical'):
 old=oldfolder/(name+'.xml');xml=ET.parse(old);si=xml.getroot().find('.//SOLVERINFO');opts={i.get('PROPERTY'):i for i in si}
 assert 'WJRestart'not in opts and 'WJInitialFields'not in opts and opts['WJStepPreconditioner'].get('VALUE')=='CurrentDataDriftV1'
 opts['WJOutputPrefix'].set('VALUE',str(o/name));path=o/(name+'.xml');xml.write(path,encoding='unicode');shutil.copy2(str(old)+'.opt',str(path)+'.opt')
 cases[name]=dict(template=identity(old),session=identity(path),options=identity(Path(str(path)+'.opt')),old=identity(old.with_suffix('.result.json')))
old=r/'science/radial-quadrature-v18/reference-q16.xml';xml=ET.parse(old);x=xml.getroot();opts={i.get('PROPERTY'):i for i in x.findall('.//SOLVERINFO/I')}
opts['WJOutputPrefix'].set('VALUE',str(o/'reference-r2'));opts['WJInitialFields'].set('VALUE',str(r/'science/radial-quadrature-v18/reference-q16.result.json'))
assert opts['WJInitialSourceIdentity'].get('VALUE')=='1c9b1fe7027814da2b575babb055cf988833e69cf11d63cf6eda34e71a2000f4'
for item in x.findall('.//PARAMETERS/P'):
 key=item.text.split('=')[0].strip()
 if key in ('WJQuadratureX','WJSteadySeconds'):item.text=key+' = '+dict(WJQuadratureX='8',WJSteadySeconds='700')[key]
path=o/'reference-r2.xml';xml.write(path,encoding='unicode');shutil.copy2(str(old)+'.opt',str(path)+'.opt')
raw=r/'science/htrace-regression-v18/reference-r2.result.json';witness=r/'science/radial-quadrature-v18/cutoff-crossing-witness-v1.json';data,_=po.load_native(raw)
ref=po.Reference(json.loads((r/'science/reference-tables-v1.json').read_text()));oracle=po.PhysicalOracle(data,ref,20,time.monotonic()+20)
w=json.loads(witness.read_text());assert identity(raw)['sha256']==w['native']['sha256'];traces=[];cells=[]
for e in oracle.elements:
 if abs(e['bounds'][2]-ref.H)<1e-12:
  h=oracle.owner[e['geometry']];assert h['geometry']==e['geometry'];xl,xr,_,_=e['bounds']
  traces.append(dict(owner=e['geometry'],left=xl,right=xr,coefficients=(oracle.s['k'][np.asarray(e['k'][:4],int)]*np.asarray(e['ks'][:4])).tolist()))
 if e['geometry']in {row['owner']for row in w['owners']}:
  ids=np.asarray(e['u'],int);signs=np.asarray(e['us']);kid=np.asarray(e['k'],int);ks=np.asarray(e['ks']);n=data['sizes'][0]
  cells.append(dict(owner=e['geometry'],bounds=e['bounds'],signed_u=(oracle.s['flow'][ids]*signs).tolist(),
    signed_v=(oracle.s['flow'][n+ids]*signs).tolist(),signed_k=(oracle.s['k'][kid]*ks).tolist(),k_signs=ks.tolist()))
assert len(traces)==28 and len(cells)==2
atomic(o/'operator-request.json',dict(schema='WJH27CutoffPanelRequestV1',source=identity(raw),traces=traces,cells=cells))
atomic(o/'cases-manifest.json',dict(schema='WJH27CutoffPanelCasesV1',cases=cases,radial_template=identity(old),radial_case=identity(path),
 source=identity(r/'source-v19.json'),runtime=identity(r/'science/runtime-v19.json'),source_identity=source,witness=identity(witness),request=identity(o/'operator-request.json'),TASK3_READY=False))
print('CUTOFF_PANEL_FIXTURES_READY')
