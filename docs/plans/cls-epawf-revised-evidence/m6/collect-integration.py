#!/usr/bin/env python3
"""Close M6 only from passed, revision-bound Linux evidence and actual outputs."""
import argparse,hashlib,itertools,json,os,shutil,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
import xml.etree.ElementTree as ET
os.environ.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
sys.dont_write_bytecode=True
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
for name in ('qualification','regression','kernel','quadrature','fingerprint'):p.add_argument('--'+name,required=True)
p.add_argument('--output',type=Path,required=True);a=p.parse_args();root=a.root;source=root/'source'
sys.path.insert(0,str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
from validate_cls_epawf_m4 import Space
from validate_cls_epawf_m5_geometry import geometry
def read(path):return json.loads(path.read_text())
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
qdir=root/a.qualification;rdir=root/a.regression;kdir=root/a.kernel;quad=root/a.quadrature
q=read(qdir/'qualification/validation.json');qr=read(qdir/'receipt.json');rr=read(rdir/'receipt.json');kr=read(kdir/'receipt.json');qc=read(quad/'receipt.json');fp=read(root/a.fingerprint)
revision=qr['source_commit'];assert all(r['source_commit']==revision for r in (rr,kr,qc,fp))
assert q['status']=='passed' and qr['status']==rr['status']==qc['status']=='passed' and kr['exit_code']==0 and fp['success']
assert all(not r.get('source_end_status') for r in (qr,rr,kr,qc))
assert git(source,'rev-parse','HEAD')==revision and not git(source,'status','--porcelain')
assert rr['main_commit']==rr['main_end_commit']=='f4b70baec10016dd098005bc5a5507b87e4a57bf' and not rr['main_end_status']
assert git(root.parent.parent/'nektar-v2026','rev-parse','HEAD')==rr['main_commit']
assert not git(root.parent.parent/'nektar-v2026','status','--porcelain')
originals={'library/SolverUtils/Turbulence/RANSModels.cpp':'86a7154ba3cb27e4d3319aa42d156a228fd1b180f35f9b217159822ec9d66b05',
 'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp':'5baf29819fed6b81ad1b042bc09258a3e0779d15de13076b688bf32f80ddcf7e',
 'docs/plans/cls-epawf-revised-contract-v1.json':'4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a',
 'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex':'d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f'}
assert all(sha(source/name)==value for name,value in originals.items())
registered=json.loads(subprocess.check_output(['ctest','--test-dir',str(root/'build'),'--show-only=json-v1'],text=True))
names={test['name'] for test in registered['tests']};assert len(names)==30
coverage={}
for path,kind in [(rdir/'ctest.xml','registered CTest regression'),(kdir/'ctest.xml','registered M6 kernel CTest')]:
 cases=ET.parse(path).getroot().findall('.//testcase')
 for case in cases:
  assert not any(case.findall(tag) for tag in ('failure','error','skipped'))
  assert case.get('name') not in coverage
  coverage[case.get('name')]={'status':'passed','execution':kind,'source_commit':revision,'receipt':str(path),'seconds':case.get('time')}
assert len(coverage)==29
coverage['RANSCLSEPAWFM6Oracle']={'status':'passed','execution':'direct execution of the registered qualification entry point, --mode all',
 'source_commit':revision,'receipt':str(qdir/'qualification/validation.json'),'commands':len(q['commands'])}
assert set(coverage)==names
kernel=read(root/'build/cls-epawf-m6-kernel/validation.json');assert kernel['success']
assert len(q['sections']['operators'])==3 and len(q['sections']['temporal'])==3 and len(q['sections']['rejections'])==11
assert q['sections']['standard_control']['stage']['maximum_scaled_residual']<1e-9
# The omitted chi derivative must fail the original 1e-10 trace criterion,
# well above observed FD comparison noise; the looser all-row FD tolerance
# alone is not used to establish that this negative control is detectable.
omissions=[]
for op in q['sections']['operators']:
 for row in op['directional'].get('blocks',[]):
  block=row['block'];error=q['errors']['implicit_trace_action_block_'+str(block)]
  omitted=row['trace_omitted_chi_derivative_maximum']
  assert error<1e-10 and omitted>1e-10 and omitted>100*error
  omissions.append({'block':block,'error':error,'omitted':omitted,'ratio_to_noise':omitted/error})
stages=[];outputs=[]
for path in sorted((qdir/'qualification').glob('*/case.evolution.json')):
 data=read(path)
 if data['stage_audits']:
  assert data['eta']==0 and data['upstream_transition']
  assert max(row[4] for row in data['stage_audits'])<1e-9
  assert np.max(np.abs(np.sum(data['regional_stage_loads'],axis=0)-data['final_stage_loads']),initial=0)<1e-10
  assert abs(data['state'][-1])<1e-9 and abs(data['energy_budget'][7])<1e-9
  stages.append({'case':path.parent.name,'ranks':data['ranks'],'order':data['order'],'steps':data['steps'],
    'unknowns':len(data['state']),'maximum_scaled_stage_norm':max(row[4] for row in data['stage_audits']),
    'gauge_multiplier':data['state'][-1],'energy_identity_error':data['energy_budget'][7]})
 if not data['mapped']:continue
 files=sorted(path.parent.glob('case-*.vtu'));chart=list(path.parent.glob('case-*.chart.fld'));d=data['dimension'];H=data['parameters'][1]
 assert len(files)==len(chart)==2*d+4
 definition=ET.parse(path.parent/'case.xml').getroot();fields=[Space(r,d) for r in data['fields']]
 velocity=['u','v'] if d==2 else ['u','v','w'];ce=ve=0.;count=0
 for item in files:
  mesh=ET.parse(item).getroot().find('UnstructuredGrid/Piece');points=np.fromstring(mesh.find('Points/DataArray').text,sep=' ').reshape(-1,3)
  values=np.fromstring(mesh.find('PointData/DataArray').text,sep=' ');assert len(points)==len(values)==int(mesh.get('NumberOfPoints'))
  assert np.all(np.isfinite(points)) and np.all(np.isfinite(values))
  conn=np.fromstring(mesh.find('Cells/DataArray[@Name="connectivity"]').text,sep=' ',dtype=int)
  assert np.all((conn>=0)&(conn<len(points)))
  name=item.name[5:-4]
  if name in velocity:block=velocity.index(name);part='outer';variable=name
  elif name.endswith('_w') and name[:-2] in velocity:block=velocity.index(name[:-2]);part='inner';variable=name[:-2]
  else:block={'k':d,'eps_tilde':d+1,'pi_w':d+2,'p':d+3}[name];part='all';variable='p' if name in ('p','pi_w') else name
  expansion=next(e for e in definition.find('EXPANSIONS') if variable in e.get('FIELDS').split(','))
  nq=[int(n) for n in expansion.get('NUMPOINTS').split(',')]
  nodes=[np.r_[-1.,np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),1.] for n in nq]
  expected_points=[];expected_values=[];field=fields[block]
  for elem in field.elements:
   if (part=='inner' and elem['upper'][2]>H+1e-11) or (part=='outer' and elem['lower'][2]<H-1e-11):continue
   chart_points=[]
   for index in itertools.product(*(range(n) for n in nq[::-1])):
    point=np.zeros(3)
    for axis,j,node in zip(field.axes,index[::-1],nodes):point[axis]=elem['lower'][axis]+(node[j]+1)*(elem['upper'][axis]-elem['lower'][axis])/2
    chart_points.append(point)
   chart_points=np.asarray(chart_points);expected_points.extend(geometry(data['map_parameters'],chart_points)['point'])
   ids,phi,_=field.rows(elem,chart_points);expected_values.extend(phi@np.asarray(data['state'])[data['offsets'][block]+ids])
  assert np.asarray(expected_points).shape==points.shape
  c=float(np.max(np.abs(points-expected_points),initial=0));v=float(np.max(np.abs(values-expected_values),initial=0))
  assert c<1e-10 and v<1e-10;ce=max(ce,c);ve=max(ve,v);count+=len(points)
 outputs.append({'case':path.parent.name,'vtu_files':len(files),'chart_files':len(chart),'points':count,'coordinate_error':ce,'field_error':ve})
contract=read(source/'docs/plans/cls-epawf-revised-contract-v1.json')
owners=[entry for entry in contract['equation_owners'] if entry['delivery_milestone']=='M6'] if 'equation_owners' in contract else []
if not owners:
 owners=[entry for value in contract.values() if isinstance(value,list) for entry in value if isinstance(entry,dict) and entry.get('delivery_milestone')=='M6']
assert {entry['id'] for entry in owners}=={'scalar_forcing_response','transition_and_adapted_baseline'}
index=[]
for directory in [qdir,rdir,kdir,quad]:
 for path in sorted(directory.rglob('*')):
  if path.is_file():index.append({'path':str(path),'bytes':path.stat().st_size,'sha256':sha(path)})
receipt={'schema':1,'status':'verified_before_orientation_guard','milestone':'M6','source_commit':revision,'collected_utc':datetime.now(timezone.utc).isoformat(),
 'registered_test_coverage':coverage,'total_registered_tests':len(names),'main_checks':q['checks'],'kernel_checks':kernel['checks'],
 'additional_independent_weak_checks':sum(op['additional_weak_checks'] for op in q['sections']['operators']),
 'equation_delivery':owners,'frozen_originals':originals,'native_stage_reports':stages,'physical_output_audits':outputs,
 'omission_controls':omissions,'quadrature':qc,'fingerprint':str(root/a.fingerprint),'fingerprint_sha256':sha(root/a.fingerprint),
 'preserved_linux_main_commit':rr['main_commit'],'physical_failure_gates':'original CLS stability and general inner-k positivity remain failed; no robust terrain accuracy claim'}
a.output.mkdir(parents=True,exist_ok=False)
copies={'verification.json':qdir/'qualification/validation.json','qualification-receipt.json':qdir/'receipt.json',
 'regression.xml':rdir/'ctest.xml','regression-receipt.json':rdir/'receipt.json','kernel.xml':kdir/'ctest.xml',
 'kernel-receipt.json':kdir/'receipt.json','kernel-verification.json':root/'build/cls-epawf-m6-kernel/validation.json',
 'quadrature-receipt.json':quad/'receipt.json','qualification-fingerprint.json':root/a.fingerprint}
for name,path in copies.items():shutil.copyfile(path,a.output/name)
(a.output/'raw-artifact-index.json').write_text(json.dumps(index,indent=2)+'\n')
(a.output/'integration.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'status':receipt['status'],'registered_tests':len(names),'stages':len(stages),'physical_output_sets':len(outputs),'raw_artifacts':len(index)}))
