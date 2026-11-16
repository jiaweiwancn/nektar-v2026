#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H25 fair physical-mesh comparison and narrow production candidates."""
import argparse,hashlib,json,math,os,re,subprocess,sys,time,uuid
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic,vector_hash,bits,RULE
from wj_h17_observer import native_load
from cls_epawf_m3_cases import set_info,set_parameter,write

HERE=Path(__file__).resolve().parent
OLD=Path('/home/jay/build/wj-h24-flat-continuation-20261115')

def measure_matches(actual,expected):
    # Derived coordinate products are not serialized payload comparisons.
    # Same 128-epsilon arithmetic bound used for reconstructed identities.
    return math.isfinite(actual)and math.isfinite(expected)and abs(actual-expected)<=128*sys.float_info.epsilon*(1+abs(actual)+abs(expected))

def put_same_or_new(path,value):
    if path.exists():assert load(path)==value,'immutable prior setup differs: '+str(path)
    else:atomic(path,value)

def stock_case(original,folder,steps=8,dt=.0025):
    """Materialize physical coordinates, never feed composite chart to VCS."""
    folder.mkdir(parents=True,exist_ok=False);x=ET.parse(original).getroot();co=x.find('CONDITIONS')
    parameters={p.text.split('=')[0].strip():float(p.text.split('=')[1])for p in co.findall('./PARAMETERS/P')}
    H=parameters['RANSWallHeight'];d0=parameters['RANSWallMapDisplacement'];Z=parameters['RANSWallMapTunnelTop']
    points={};coordinates=[]
    for v in x.findall('./GEOMETRY/VERTEX/V'):
        a,b,c=map(float,v.text.split());z=d0+(b if b<=H else H+(b-H)*(Z-d0-H)/(Z-H))
        v.text=f'{a:.17g} {z:.17g} {c:.17g}';points[int(v.get('ID'))]=(a,z);coordinates.append((a,z))
    xs=sorted({p[0]for p in coordinates});zs=sorted({p[1]for p in coordinates})
    assert len(xs)==15 and len(zs)==51 and len(x.findall('./GEOMETRY/ELEMENT/Q'))==700
    assert min(zs)==d0 and max(zs)==Z and any(z==d0+H for z in zs)
    measures=dict(full=(xs[-1]-xs[0])*(Z-d0),inner=(xs[-1]-xs[0])*H,outer=(xs[-1]-xs[0])*(Z-d0-H))
    assert all(measure_matches(measures[k],v)for k,v in dict(full=497.8125,inner=3.75,outer=494.0625).items())
    geometry=x.find('GEOMETRY');composite=geometry.find('COMPOSITE')
    for e in list(composite):
        if e.get('ID')!='0':composite.remove(e)
    edges={n:[]for n in ('ground','inlet','outlet','top')}
    for e in geometry.findall('./EDGE/E'):
        a,b=(points[int(i)]for i in e.text.split());eid=int(e.get('ID'))
        for name,axis,value in [('ground',1,d0),('inlet',0,xs[0]),('outlet',0,xs[-1]),('top',1,Z)]:
            if a[axis]==b[axis]==value:edges[name].append(eid)
    for n,(name,ids)in enumerate(edges.items(),1):ET.SubElement(composite,'C',ID=str(n)).text='E['+','.join(map(str,ids))+']'
    inner=co.find("FUNCTION[@NAME='InitialConditionsInner']/E[@VAR='u']").get('VALUE')
    outer=co.find("FUNCTION[@NAME='InitialConditionsOuter']/E[@VAR='u']").get('VALUE')
    # Inactive outer logarithm is finite even at the no-slip ground.
    r=f'(y-({d0:.17g}))';outer=outer.replace(r,f'max({r},{H:.17g})')
    U=f'(y<{d0+H:.17g})*({inner})+(y>={d0+H:.17g})*({outer})'
    topU=outer.replace('y',f'({Z:.17g})')
    for child in list(co):co.remove(child)
    si=ET.SubElement(co,'SOLVERINFO')
    for k,v in dict(SolverType='VelocityCorrectionScheme',EQTYPE='UnsteadyNavierStokes',AdvectionForm='Convective',Projection='Galerkin',GlobalSysSoln='DirectMultiLevelStaticCond').items():ET.SubElement(si,'I',PROPERTY=k,VALUE=v)
    scheme=ET.SubElement(co,'TIMEINTEGRATIONSCHEME');ET.SubElement(scheme,'METHOD').text='IMEX';ET.SubElement(scheme,'ORDER').text='1'
    ps=ET.SubElement(co,'PARAMETERS')
    for k,v in dict(Time=.04,TimeStep=dt,NumSteps=steps,Kinvis=parameters['Kinvis'],IO_InfoSteps=1,IO_CheckSteps=0,IO_Timer_Level=4).items():ET.SubElement(ps,'P').text=f'{k} = {v}'
    vs=ET.SubElement(co,'VARIABLES')
    for n,v in enumerate(('u','v','p')):ET.SubElement(vs,'V',ID=str(n)).text=v
    es=x.find('EXPANSIONS');es.clear()
    for fields,nm in [('u,v',4),('p',3)]:ET.SubElement(es,'E',COMPOSITE='C[0]',FIELDS=fields,NUMMODES=f'{nm},{nm}',BASISTYPE='Modified_A,Modified_A',POINTSTYPE='GaussLobattoLegendre,GaussLobattoLegendre',NUMPOINTS='12,16')
    br=ET.SubElement(co,'BOUNDARYREGIONS');bc=ET.SubElement(co,'BOUNDARYCONDITIONS')
    for i,name in enumerate(edges):
        ET.SubElement(br,'B',ID=str(i)).text=f'C[{i+1}]';reg=ET.SubElement(bc,'REGION',REF=str(i))
        if name=='outlet':
            ET.SubElement(reg,'N',VAR='u',VALUE='0');ET.SubElement(reg,'N',VAR='v',VALUE='0');ET.SubElement(reg,'D',VAR='p',VALUE='0')
        else:
            ET.SubElement(reg,'D',VAR='u',VALUE='0'if name=='ground'else U if name=='inlet'else topU)
            ET.SubElement(reg,'D',VAR='v',VALUE='0');ET.SubElement(reg,'N',VAR='p',USERDEFINEDTYPE='H',VALUE='0')
    f=ET.SubElement(co,'FUNCTION',NAME='InitialConditions')
    for name,value in [('u',U),('v','0'),('p','0')]:ET.SubElement(f,'E',VAR=name,VALUE=value)
    dest=folder/'case.xml';write(x,dest)
    atomic(folder/'physical-mesh.json',dict(original=identity(original),stock=identity(dest),vertices=points,cell_count=700,physical_measures=measures,
        physical_interface=d0+H,chart_input=False,velocity_global_dofs=2*(3*14+1)*(3*50+1),pressure_global_dofs=(2*14+1)*(2*50+1),
        pressure='global continuous Q2, not WJ regional duplicated Q1',governing_points=[12,16],native_WJ_representation_points=[8,8],
        differences=['positive molecular viscosity throughout stock, high-Re WJ outer omits molecular viscosity','stock fixed top velocity versus WJ shear traction','stock pressure-fixed natural-velocity outlet versus complete WJ stress traction','stock projected analytic initial U versus retained .04 WJ velocity','stock contains no k/epsilon/AWF equations']))
    return dest

def setup(root):
    from wj_h16_observer import compatibility
    for name in ('runtime-build-v1.json','runtime-dependencies-v1.json','compiled-source-objects.json'):
        old=OLD/'science'/name;put_same_or_new(root/'science'/name,load(old))
    runtime=load(root/'science/runtime-build-v1.json')
    for e in runtime['executables']:verify(e['ELF']);verify(e['wrapper'])
    for v in load(root/'science/runtime-dependencies-v1.json')['libraries']:verify(v['actual'])
    commit=load(OLD/'science/M/monitor/step-16.qualification-commit.json')
    anchor={k:commit[k]for k in ('checkpoint','callback','report')};anchor['qualified_commit']=identity(OLD/'science/M/monitor/step-16.qualification-commit.json')
    anchor['completion']=identity(OLD/'science/M/monitor/step-16.completion.json')
    for v in anchor.values():verify(v)
    selected=load(OLD/'inputs/selected.json')
    inputs=dict(anchor=anchor,own_cp0=selected['cp0']['M'],original=selected['original'],support=selected['common_support'])
    put_same_or_new(root/'inputs/baseline.json',inputs)
    oldcompat=load(OLD/'science/target-cache-compatibility.json')
    compatibility(root,verify(oldcompat['historical_manifest']),root/'science/cache-compatible-v2.json')
    original=verify(selected['cp0']['M']['xml']);stock_case(original,root/'science/S0-v2',8)
    folder=root/'science/W0';folder.mkdir();x=ET.parse(original).getroot()
    for k,v in dict(NumSteps=8,IO_CheckSteps=1,IO_InfoSteps=1).items():set_parameter(x,k,v)
    for k,v in dict(RANSWallRestart=str(verify(anchor['checkpoint']).parent),RANSWallCheckpoint=str(folder/'CP'),RANSWallRetention='BoundedV1',
        RANSWallMonitoringPolicy='CheckpointAuditV1',RANSWallMonitorIntent=str(folder/'intent.json'),RANSWallMonitorFault='Off',
        RANSWallReport=str(folder/'native.json'),RANSWallPimpleSnapshots='Off',RANSWallPimpleFlowVectors='Off',RANSWallPimpleArithmetic='Off',RANSWallPimpleVerify='Off',RANSWallImmutableRowReuse='Off').items():set_info(x,k,v)
    write(x,folder/'case.xml')
    atomic(root/'science/setup-v2.json',dict(status='PASS',inputs=identity(root/'inputs/baseline.json'),case_S0=identity(root/'science/S0-v2/case.xml'),
        case_W0=identity(folder/'case.xml'),native=runtime,compatibility=identity(root/'science/cache-compatible-v2.json'),
        source_files=[identity(p)for p in sorted(HERE.glob('*h25*.py'))],native_changes=0))

def executable(root,name):return verify(next(e['wrapper']for e in load(root/'science/runtime-build-v1.json')['executables']if Path(e['wrapper']['path']).name==name))

def stock(root):
    folder=root/'science/S0-v2';started=time.monotonic()
    with(folder/'native.log').open('x')as f:code=subprocess.run([str(executable(root,'IncNavierStokesSolver')),str(folder/'case.xml')],stdout=f,stderr=subprocess.STDOUT).returncode
    text=(folder/'native.log').read_text();samples=[dict(step=int(s),time=float(t),seconds=float(w))for s,t,w in re.findall(r'Steps:\s*(\d+)\s+Time:\s*([\d.eE+-]+)\s+CPU Time:\s*([\d.eE+-]+)s',text)]
    atomic(folder/'result.json',dict(status='PASS'if code==0 and len(samples)==8 else'FAIL',exit=code,seconds=time.monotonic()-started,samples=samples,producer=identity(folder/'native.log'),
        SVV=False,meaning='cost reference, different equations/pressure/boundaries disclosed'))
    assert code==0 and len(samples)==8,'stock launch/step failure: retained log'

def baseline(root):
    folder=root/'science/W0';d=folder/'monitor';d.mkdir();ip=folder/'intent.json';inputs=load(root/'inputs/baseline.json')
    v=dict(schema='WJMonitorIntentV1',policy='CheckpointAuditV1',directory=str(d),nonce=uuid.uuid4().hex,
        scope_sha256=identity(root/'science/setup-v2.json')['sha256'],source_sha256=identity(root/'science/runtime-build-v1.json')['sha256'],observer_sha256=identity(HERE/'wj_h25_observer.py')['sha256'],inventory_rule_sha256=hashlib.sha256(RULE.encode()).hexdigest(),
        original_qh0_checkpoint=inputs['own_cp0']['checkpoint'],original_qh0_state_sha256=inputs['original']['coupled_state_sha256'],
        spatial_lineage_anchor=inputs['original']['checkpoint'],restored_checkpoint=inputs['anchor']['checkpoint'],restart_completion=inputs['anchor']['completion'],
        completed_anchor=inputs['anchor'],maximum_absolute_step=24,restored_absolute_step=16,recovery_stride=1,scientific_stride=8,dt_bits=bits(.0025),
        ack_timeout_seconds=120.,substep_evidence_steps=[],flow_capture='Off',snapshot_capture='Off',preserve_qualification_artifacts=True,
        cache_compatibility=identity(root/'science/cache-compatible-v2.json'),common_support=inputs['support'])
    obs=native=None;code=None;error=None;began=time.monotonic()
    with(folder/'observer.log').open('x')as ol,(folder/'native.log').open('x')as nl:
        try:
            obs=subprocess.Popen([sys.executable,'-B',str(HERE/'wj_h25_observer.py'),'--root',str(root),'--intent',str(ip)],stdout=ol,stderr=subprocess.STDOUT)
            v['observer_pid']=obs.pid;atomic(ip,v)
            while not(d/'observer-started.json').exists():assert obs.poll()is None and time.monotonic()-began<30;time.sleep(.025)
            native=subprocess.Popen([str(executable(root,'IncNavierStokesSolver')),str(folder/'case.xml')],stdout=nl,stderr=subprocess.STDOUT)
            atomic(folder/'native-launch.json',dict(pid=native.pid,xml=identity(folder/'case.xml'),intent=identity(ip),additional=8))
            while native.poll()is None:assert obs.poll()is None,'observer exited before native';time.sleep(.025)
            code=native.returncode
        except Exception as e:error=repr(e)
        finally:
            if native and native.poll()is None:native.terminate();native.wait(timeout=5)
            atomic(d/'observer-stop.json',dict(native_exit=code,reason='coordinator after native exit/failure'))
            if obs:
                try:obs.wait(timeout=5)
                except subprocess.TimeoutExpired:obs.terminate();obs.wait(timeout=5)
    p=folder/'case.evolution.json';r=native_load(p)if p.exists()else None
    actual=int(r['implicit_callbacks'])if r else 8
    atomic(Path(os.environ['WJ_H25_OPERATION'])/'attempt-counts.json',dict(actual_attempts=actual,known=r is not None))
    atomic(folder/'outcome.json',dict(status='PASS'if code==0 and error is None else'FAIL',code=code,error=error,seconds=time.monotonic()-began,producer=identity(p)if r else None))
    assert code==0 and error is None,(code,error,r.get('failure')if r else None)
    assert r['initialization_work']==[0,0]and r['initializer_calls']==0 and int(r['steps'])==24
    assert int(load(d/'latest-fully-verified-checkpoint.json')['binding']['absolute_step'])==24

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('setup','stock','baseline'),required=True)
    a=p.parse_args();globals()[a.mode](a.root)
