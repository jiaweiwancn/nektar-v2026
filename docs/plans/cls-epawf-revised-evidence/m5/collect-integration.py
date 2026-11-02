#!/usr/bin/env python3
"""Close M5 only from successful source-bound Linux qualification/regressions."""
import argparse,hashlib,json,subprocess,os,sys,itertools
os.environ.update({name:'1' for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
from datetime import datetime,timezone
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--qualification',required=True);p.add_argument('--regression',required=True);p.add_argument('--note',required=True);p.add_argument('--replay',required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
root=a.root;source=root/'source';qdir=root/a.qualification;rdir=root/a.regression
sys.path.insert(0,str(source/'solvers/IncNavierStokesSolver/Tests/RANS'))
from validate_cls_epawf_m4 import Space
from validate_cls_epawf_m5_geometry import geometry
q=json.loads((qdir/'oracle/validation.json').read_text());qr=json.loads((qdir/'receipt.json').read_text());rr=json.loads((rdir/'receipt.json').read_text())
assert q['verification_status']=='passed' and q['milestone_closure']
assert qr['status']=='passed'
fix_dir=root/a.replay;fix=json.loads((fix_dir/'receipt.json').read_text());replay=json.loads((fix_dir/'replay/validation.json').read_text())
assert fix['status']=='passed' and replay['success']
assert q['source_commit']==qr['source_commit']
note_dir=root/a.note;nr=json.loads((note_dir/'receipt.json').read_text());note=json.loads((note_dir/'note/validation.json').read_text())
assert nr['status']=='passed' and note['success'] and nr['source_commit']==rr['source_commit']
qfp=json.loads((qdir/'qualification-fingerprint.json').read_text());rfp=json.loads((fix_dir/'replay-fingerprint.json').read_text())
assert qfp['source_commit']==q['source_commit'] and rfp['source_commit']==fix['source_commit']
assert rfp['identical_equations_and_solver_runtime']
changed=rfp['probe_session_filename_only_changes']
assert all(rfp[k].get(name)==value or name in changed[k] for k in ('source_inputs','runtime_artifacts') for name,value in qfp[k].items())
assert not qr['source_end_status'] and not rr['source_end_status'] and not rr['main_end_status']
assert rr['main_commit']==rr['main_end_commit']=='f4b70baec10016dd098005bc5a5507b87e4a57bf'
tests=ET.parse(rdir/'ctest.xml').getroot();cases=tests.findall('.//testcase')
assert len(cases)==26 and all(not case.findall('error') and not case.findall('skipped') for case in cases)
failed=[c for c in cases if c.findall('failure')]
assert len(failed)==1 and failed[0].get('name')=='RANSEquilibriumEvolutionNative'
repeated=sorted((fix_dir/'replay').glob('legacy-repeat-*.xml'))
assert len(repeated)==2
for item in repeated:
 checks=ET.parse(item).getroot().findall('.//testcase')
 assert len(checks)==1 and checks[0].get('name')=='RANSEquilibriumEvolutionNative' and not checks[0].findall('failure')
assert any(c.get('name')=='RANSCLSEPAWFM4Oracle' for c in cases)
assert all(c.get('name')!='RANSCLSEPAWFM5Oracle' for c in cases)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def git(*args):return subprocess.check_output(['git','-c','gc.auto=0','-C',str(source),*args],text=True).strip()
assert git('rev-parse','HEAD')==fix['source_commit'] and not git('status','--porcelain')
a.output.mkdir(exist_ok=False,parents=True)
# Preserve the complete command/results reports and exact runner receipts.
for src,name in [(qdir/'oracle/validation.json','verification.json'),(qdir/'receipt.json','qualification-receipt.json'),(rdir/'receipt.json','regression-receipt.json'),(rdir/'ctest.xml','original-regression.xml'),(note_dir/'note/validation.json','note-verification.json'),(note_dir/'receipt.json','note-receipt.json'),(qdir/'qualification-fingerprint.json','qualification-fingerprint.json'),(fix_dir/'replay-fingerprint.json','replay-fingerprint.json'),(fix_dir/'replay/validation.json','session-replay-verification.json'),(fix_dir/'receipt.json','session-replay-receipt.json')]:
 (a.output/name).write_bytes(src.read_bytes())
# Compose coverage transparently from 25 unchanged tests and the corrected
# repeated probe. The original failed run and both passing replay XMLs remain.
combined=ET.Element('testsuite',{'name':'M5 composed regression coverage','tests':'26','failures':'0','errors':'0'})
for case in cases:
 if case.findall('failure'):continue
 props=ET.SubElement(case,'properties');ET.SubElement(props,'property',{'name':'source_commit','value':rr['source_commit']});combined.append(case)
case=ET.parse(repeated[-1]).getroot().findall('.//testcase')[0];props=ET.SubElement(case,'properties');ET.SubElement(props,'property',{'name':'source_commit','value':fix['source_commit']});combined.append(case)
ET.ElementTree(combined).write(a.output/'regression.xml',encoding='utf-8',xml_declaration=True)
for i,item in enumerate(repeated,1):(a.output/f'legacy-repeat-{i}.xml').write_bytes(item.read_bytes())
receipts=[]
for path in sorted(root.glob('m5-attempt-*/receipt.json')):
 data=json.loads(path.read_text());receipts.append(dict(path=str(path),sha256=sha(path),receipt=data))
# Actual final-stage loads, physical outputs and max constraints are collected
# from solver reports, rather than inferred from successful process exits.
metric_audits=[]
for path in sorted((qdir/'oracle').glob('operator-*/probe.json')):
 data=json.loads(path.read_text());levels=[]
 for level in data['weak_metric']:
  levels.append(dict(order=level['order'],maximum_all_mode_weak_metric=max(f['maximum'] for f in level['fields']),
                     maximum_shared_lateral_load=max(f['seam_maximum'] for f in level['fields']),
                     maximum_constant_flux_boundary=max(max(map(abs,f['constant_flux_boundary'])) for f in level['fields'])))
 assert all(item['maximum_all_mode_weak_metric']<1e-10 and item['maximum_shared_lateral_load']<1e-10 and item['maximum_constant_flux_boundary']<1e-10 for item in levels if item['order']>=12)
 metric_audits.append(dict(case=path.parent.name,levels=levels))
stages=[];outputs=[]
for path in sorted((qdir/'oracle').glob('*/case.evolution.json')):
 data=json.loads(path.read_text())
 if not data['mapped']:continue
 if data['stage_audits']:
  assert max(record[4] for record in data['stage_audits'])<=1e-9
  if data['eta']==0:
   R=np.asarray(data['final_stage_loads']);regions=np.asarray(data['regional_stage_loads'])
   assert np.max(np.abs(regions.sum(axis=0)-R),initial=0)<1e-10
   assert abs(data['state'][-1])<1e-9
  stages.append(dict(case=path.parent.name,unknowns=len(data['state']),ranks=data['ranks'],chi=data['chi'],eta=data['eta'],steps=data['steps'],order=data['order'],
                     maximum_stage_norm=max(r[4] for r in data['stage_audits']),gauge_multiplier=data['state'][-1],boundary_surface_rms_correction=data['boundary_flux_correction'],energy_budget=data['energy_budget']))
 files=sorted(path.parent.glob('case-*.vtu'));chart=list(path.parent.glob('case-*.chart.fld'))
 assert len(files)==2*data['dimension']+4 and len(chart)==len(files)
 count=0;coordinate_error=0.;value_error=0.
 definition=ET.parse(path.parent/'case.xml').getroot();dimension=data['dimension'];H=data['parameters'][1]
 velocity_names=['u','v'] if dimension==2 else ['u','v','w']
 native_spaces=[Space(record,dimension) for record in data['fields']]
 for item in files:
  mesh=ET.parse(item).getroot().find('UnstructuredGrid/Piece');assert mesh is not None
  points=np.fromstring(mesh.find('Points/DataArray').text,sep=' ').reshape(-1,3)
  values=np.fromstring(mesh.find('PointData/DataArray').text,sep=' ')
  assert len(points)==len(values)==int(mesh.get('NumberOfPoints')) and np.all(np.isfinite(points)) and np.all(np.isfinite(values))
  conn=mesh.find('Cells/DataArray[@Name="connectivity"]');ids=np.fromstring(conn.text,sep=' ',dtype=int)
  assert np.all((ids>=0)&(ids<len(points)))
  name=item.name[len('case-'):-len('.vtu')]
  if name in velocity_names:block=velocity_names.index(name);part='outer';variable=name
  elif name.endswith('_w') and name[:-2] in velocity_names:block=velocity_names.index(name[:-2]);part='inner';variable=name[:-2]
  else:
   block={'k':dimension,'eps_tilde':dimension+1,'pi_w':dimension+2,'p':dimension+3}[name];part='all';variable='p' if name in ('p','pi_w') else name
  expansion=next(e for e in definition.find('EXPANSIONS') if variable in e.get('FIELDS').split(','))
  nq=[int(n) for n in expansion.get('NUMPOINTS').split(',')]
  nodes=[np.r_[-1.,np.polynomial.legendre.Legendre.basis(n-1).deriv().roots(),1.] for n in nq]
  expected_points=[];expected_values=[];field=native_spaces[block]
  for elem in field.elements:
   if part=='inner' and elem['upper'][2]>H+1e-11 or part=='outer' and elem['lower'][2]<H-1e-11:continue
   chart_points=[]
   for index in itertools.product(*(range(n) for n in nq[::-1])):
    node_chart=np.zeros(3)
    for axis,j,node in zip(field.axes,index[::-1],nodes):node_chart[axis]=elem['lower'][axis]+(node[j]+1)*(elem['upper'][axis]-elem['lower'][axis])/2
    chart_points.append(node_chart)
   chart_points=np.asarray(chart_points);expected_points.extend(geometry(data['map_parameters'],chart_points)['point'])
   local,phi,_=field.rows(elem,chart_points);expected_values.extend(phi@np.asarray(data['state'])[data['offsets'][block]+local])
  assert np.asarray(expected_points).shape==points.shape
  ce=float(np.max(np.abs(points-np.asarray(expected_points)),initial=0.));ve=float(np.max(np.abs(values-np.asarray(expected_values)),initial=0.))
  assert ce<1e-10 and ve<1e-10,(item,ce,ve)
  coordinate_error=max(coordinate_error,ce);value_error=max(value_error,ve);count+=len(points)
 outputs.append(dict(case=path.parent.name,physical_vtu_files=len(files),explicit_chart_field_files=len(chart),finite_physical_points=count,maximum_cartesian_coordinate_error=coordinate_error,maximum_independent_field_value_error=value_error))
# Read the frozen delivery owner keys without changing their planned baseline.
contract=json.loads((source/'docs/plans/cls-epawf-revised-contract-v1.json').read_text())
def keys(obj):
 if isinstance(obj,dict):
  for key,value in obj.items():
   if key.startswith('M5.'):yield key
   yield from keys(value)
 elif isinstance(obj,list):
  for value in obj:yield from keys(value)
ids=sorted(set(keys(contract)));assert len(ids)==10
sections={'normal_geometry':'geometry; weak_operators geometry audits', 'curved_continuity':'all pressure residual blocks; actual stage constraints; boundary projection',
 'curved_momentum_and_roughness_force':'all Cartesian momentum blocks; force/work and energy budgets',
 'tangential_and_normal_momentum':'full traction oracle; circular normal acceleration and manufactured pressure refinement',
 'shared_momentum_and_pressure':'regional coupled-stage trace loads; common physical pressure gauge',
 'curved_scalar_transport':'independent k/outer-epsilon loads; physical h/p refinement',
 'shared_k_flux':'independent interface k fluxes; free shared coupled-stage reactions',
 'dissipation_constraint_reactions':'dependent-field analytic gradient, inner diagnostic residual and outer flux reaction',
 'discrete_metric_identity':'all-mode actual native weak Piola 8/12/16; lateral and constant-flux loads',
 'curved_manufactured_examples':'circular/variable-H/full-3D/C1 MMS; partition/MPI/restart; flat/plane native IMEX repeats'}
delivery={key:dict(status='closed',evidence=sections[key.split('.',1)[1]],qualification=str(qdir/'oracle/validation.json')) for key in ids}
index=[]
for directory in [*sorted(root.glob('m5-attempt-*')),root/'build']:
 for item in sorted(directory.rglob('*')):
  if not item.is_file() or item.is_symlink():continue
  if directory.name=='build' and not ('cls-epawf-m4-' in str(item) or item.name in ('IncNavierStokesSolver','RANSCurvedWallProbe','RANSCurvedGeometryProbe','RANSCurvedNoteProbe')):continue
  index.append(dict(path=str(item),bytes=item.stat().st_size,sha256=sha(item)))
(a.output/'raw-artifact-index.json').write_text(json.dumps(index,indent=2)+'\n')
result=dict(schema=1,milestone='M5',status='closed',qualified_source_commit=q['source_commit'],supplementary_and_regression_source_commit=rr['source_commit'],session_replay_source_commit=fix['source_commit'],identical_equations_and_solver_runtime=True,probe_session_filename_only_changes=changed,added_inputs=rfp.get('added_inputs',[]),connection='tailscale:100.105.247.110',
            qualification=str(qdir/'oracle/validation.json'),regression_receipt=str(rdir/'receipt.json'),
            registered_test_coverage=dict(total=28,direct_registered_M5_oracle=1,direct_registered_note_oracle=1,other_passing_CTests=25,corrected_CTest_passes=1,corrected_CTest_repeated_executions=2,duplicate_M5_execution=False),
            stage_audits=stages,physical_output_audits=outputs,actual_native_metric_audits=metric_audits,equation_delivery=delivery,attempts=receipts,
            raw_artifact_index=dict(path=str(a.output/'raw-artifact-index.json'),sha256=sha(a.output/'raw-artifact-index.json'),files=len(index)),
            preserved_main=dict(branch='codex/hit-inflow-generator',commit=rr['main_end_commit'],clean=True),
            frozen_sources=q['frozen_sources'],compiler_workers=16,maximum_sequential_MPI_ranks=4,thread_environment=q['threads'],
            physical_gates=dict(original_CLS_stability='failed as specified',general_inner_k_positivity='failed as specified'),
            finished_utc=datetime.now(timezone.utc).isoformat())
(a.output/'integration.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(status='closed',registered_tests=28,actual_stage_reports=len(stages),equation_delivery=len(delivery),artifacts=len(index))))
