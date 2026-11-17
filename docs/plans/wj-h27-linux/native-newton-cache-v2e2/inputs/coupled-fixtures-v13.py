from pathlib import Path
import sys,json,hashlib,xml.etree.ElementTree as E
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');t=r/'source-v13/solvers/IncNavierStokesSolver/Tests/RANS/h27';sys.path.insert(0,str(t));from coupled_fixtures import generate
out=r/'science/coupled-fixtures-v13';out.mkdir();table=r/'science/reference-tables-v1.json';rows={}
template='/home/jay/build/wj-h26-native-flow-20261117/science/profile-v3/case.xml'
def case(name,cache=True,uses=8,drift=.02,step_seconds=5,**kw):
 p=out/(name+'.xml');row=generate(p,table,out/name,backend='FastFull3DV1',anderson_depth=0,nonlinear_trace=True,**kw)
 xml=E.parse(p).getroot();si=xml.find('CONDITIONS/SOLVERINFO');pars=xml.find('CONDITIONS/PARAMETERS')
 for key,value in dict(WJImplicitSolver='CoupledJFNKV1',WJStepPreconditioner='OperatorDriftV1'if cache else'RebuildEveryNewtonV1').items():E.SubElement(si,'I',PROPERTY=key,VALUE=value)
 for key,value in dict(WJStepNewtonSeconds=step_seconds,WJStepNewtonIterations=20,WJStepKrylovDimension=40,WJStepLinearCycles=8,WJStepFDStep=1e-6,WJStepPCRelativeDrift=drift,WJStepPCMaximumUses=uses).items():E.SubElement(pars,'P').text=f'{key} = {value}'
 E.indent(xml);p.write_text(E.tostring(xml,encoding='unicode')+'\n');row.update(sha256=hashlib.sha256(p.read_bytes()).hexdigest(),preconditioner='OperatorDriftV1'if cache else'RebuildEveryNewtonV1');rows[name]=row
poly=dict(mms_variant='PolynomialScalarsV1',mms_body_amplitude=1,dt=.02,steps=4,linear=3e-13,nonlinear=3e-12,check_current_actions=True)
for name,cache in [('poly-fresh',False),('poly-cache',True)]:case(name,cache,uses=3,**poly)
for name,cache in [('hist-fresh',False),('hist-cache',True)]:case(name,cache,uses=3,steps=8,mms_variant='Historical',dt=.0025,linear=1e-11,nonlinear=3e-11,quadrature_x=12,quadrature_r=18,check_current_actions=True)
case('poly-driftzero',drift=0,**poly)
physical=dict(mode='Physical',dt=.005,linear=1e-12,nonlinear=3e-12,maximum_iterations=180,template=template,recovery=3,science=3)
for name,cache in [('physical-fresh',False),('physical-cache',True)]:case(name,cache,step_seconds=50,steps=12,**physical)
case('physical-prefix',step_seconds=50,steps=3,**physical)
case('physical-cold',step_seconds=50,steps=9,restart=out/'physical-prefix.final.bin',**physical)
case('physical-recovery9',step_seconds=50,steps=3,restart=out/'physical-cache.recovery-9.bin',**physical)
case('physical-reload12',step_seconds=50,steps=0,restart=out/'physical-cache.recovery-12.bin',**physical)
case('physical-retry',step_seconds=50,steps=12,retry=4,**physical)
(out/'cases-manifest.json').write_text(json.dumps(dict(schema='WJH27Source13CacheMatrixV1',scope='Equivalent inverse reuse on approved v2; physical reference/full duration/dt/terminal/production gates remain open',cases=rows),indent=2)+'\n')
