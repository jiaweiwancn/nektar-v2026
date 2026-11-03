#!/usr/bin/env python3
"""Bind M7 verification and failed physical gates to immutable Linux evidence."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import xml.etree.ElementTree as ET

p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();root=a.root.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
assert platform.system()=='Linux'
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(name):return json.loads((root/name).read_text())
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
checks=0
def require(value,message):
 global checks
 checks+=1
 if not value:raise RuntimeError(message)
receipt=dict(schema=1,milestone='M7',status='failed',started_utc=datetime.now(timezone.utc).isoformat(),
             collector_sha256=sha(__file__),connection='tailscale:100.105.247.110')
try:
 native_parent=read('m7-native-final-001/receipt.json');native=read('m7-native-final-001/qualification/receipt.json')
 legacy=read('m7-regression-001/receipt.json');batch=read('m7-final-audits-001/receipt.json')
 replay=read('m7-replay-003/receipt.json');transverse=read('m7-strict-transverse-002/validation.json')
 dense=read('m7-dense-guard-003/receipt.json');active=read('m7-energy-preflight-003/receipt.json')
 flow=read('m7-flow-replay-001/validation.json');flow_original=read('m7-flow-diagnosis-001/validation.json')
 for name,obj in [('native parent',native_parent),('native',native),('legacy',legacy),('audit batch',batch),('replay',replay),('dense',dense),('active energy',active)]:
  require(obj.get('status')=='passed',name+' not passed')
 require(native_parent.get('finished_utc') and legacy.get('finished_utc') and batch.get('finished_utc'),'unfinished runner receipt')
 require(transverse['success'] and transverse['physical_stability_gate']=='failed_as_specified','missing original-law growth witness')
 require(flow['replay_status']=='passed' and flow_original['status']=='diagnosis_collected' and not flow_original['collection_errors'],'incomplete stopped-flow diagnosis')
 require(not flow['passage_qualified'] and not flow_original['passage_qualified'],'unqualified stopped trajectory promoted')
 require(dense['result']['sole_failed_admissibility_set']=='reconstruction_dense','focused reconstruction guard missing')
 tree=ET.parse(root/'m7-regression-001/ctest.xml')
 tests=tree.findall('.//testcase');require(len(tests)==31,'legacy CTest count changed')
 require(all(t.find('failure') is None and t.find('error') is None and t.find('skipped') is None for t in tests),'legacy failure/skip')
 registry=json.loads(subprocess.check_output(['ctest','--test-dir',str(root/'build'),'--show-only=json-v1'],text=True))
 names={x['name'] for x in registry['tests']};old_names={x.attrib['name'] for x in tests}
 new_names={'RANSCLSEPAWFM7NativeOracle','RANSCLSEPAWFM7TransverseOracle'}
 require(len(names)==33 and old_names|new_names==names and not old_names&new_names,'registered test coverage mismatch')
 covariance=[]
 for record in flow['cases']:
  require(record['observation_verified'] and record['status']=='verified_stop','unverified/still-running flow observation')
  require(record['exit_code']!=0 and not record['trajectory_completed'],'flow stop misclassified')
  covariance.append(dict(name=record['name'],last_time=record['observation']['last_accepted_time'],
      minimum_outer_covariance=record['independent_last_accepted_energy']['direct']['outer_covariance_minimum_on_declared_Gauss_rule'],
      min_k=record['accepted_prefix']['min_k'],classification=record['observation']['classification']))
 successes=[];rejections=[];source=root/'m7-audit-source'
 for path in sorted((root/'m7-native-final-001/qualification').glob('*/case.evolution.json')):
  raw=json.loads(path.read_text())
  if raw.get('failure'):
   require(not raw.get('accepted_states'),'invalid initial-state control accepted');rejections.append(path.parent.name);continue
  require(raw.get('diagnostics_schema')==1,'missing final accepted diagnostics')
  require(len(raw['accepted_states'])==len(raw['stage_audits'])+1,'accepted stage coverage')
  require(all(audit['admissible']==[1] for audit in raw['state_audits']),'invalid accepted state')
  require(all(row[4]<=1e-9 for row in raw['stage_audits']),'unconverged accepted native stage')
  successes.append(dict(name=path.parent.name,states=len(raw['accepted_states']),
      stage_residual=max([row[4] for row in raw['stage_audits']],default=0.),report_sha256=sha(path)))
 require(len(successes)==25 and len(rejections)==3,'final native case coverage changed')
 require(sum(x['states'] for x in successes)==135,'native state audit count changed')
 pin=read('m7-attempt-004/runtime-fingerprint.json')
 qualified=read('m7-qualified-sync/qualified-runtime-fingerprint.json')
 require(qualified['success'] and qualified['source_commit']==native_parent['source_commit'],'qualified runtime snapshot missing')
 for name in ('runtime_artifacts','runtime_bindings','runtime_environment','python_packages','tool_images'):
  require(qualified[name]==pin[name],'campaign-to-qualified runtime context changed: '+name)
 for key,expected in pin['runtime_artifacts'].items():require(sha(pin['runtime_resolved_paths'][key])==expected,'runtime changed: '+key)
 require(not git(source,'status','--porcelain'),'audit source dirty')
 require(not git(root/'source','status','--porcelain'),'operator source dirty')
 main=Path('/media/jay/新加卷/nektar-v2026');require(not git(main,'status','--porcelain'),'user Linux checkout dirty')
 require(git(main,'rev-parse','HEAD')==legacy['main_commit'],'user Linux revision changed')
 # Frozen originals are enforced independently of optional replay-driver changes.
 frozen={
 'docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex':'d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f',
 'docs/plans/cls-epawf-revised-contract-v1.json':'4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a',
 'library/SolverUtils/Turbulence/RANSModels.cpp':'86a7154ba3cb27e4d3319aa42d156a228fd1b180f35f9b217159822ec9d66b05',
 'library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp':'5baf29819fed6b81ad1b042bc09258a3e0779d15de13076b688bf32f80ddcf7e'}
 for path,expected in frozen.items():require(sha(source/path)==expected,'frozen input changed: '+path)
 coverage=[dict(name=name,source_commit=legacy['source_commit'],evidence='m7-regression-001/ctest.xml') for name in sorted(old_names)]
 coverage += [dict(name='RANSCLSEPAWFM7NativeOracle',source_commit=native_parent['source_commit'],compiled_physics_commit=native_parent['compiled_physics_commit'],evidence='m7-native-final-001/qualification/receipt.json'),
              dict(name='RANSCLSEPAWFM7TransverseOracle',source_commit=transverse['source_commit'],compiled_physics_commit=pin['source_commit'],evidence='m7-strict-transverse-002/validation.json')]
 equation_owners={
  'M7.matching_height_dependence':dict(status='verified_reference_only',evidence='m7-native-final-001/qualification/receipt.json',key='results.matching_height_study'),
  'M7.full_energy_accounting':dict(status='verified_accounting_identity_not_stability_bound',evidence=['m7-native-final-001/qualification/receipt.json','m7-energy-preflight-003/receipt.json','m7-flow-replay-001/validation.json'])}
 # Index full raw results outside the repository, including failed attempts.
 artifacts=[]
 for directory in sorted(root.glob('m7-*')):
  if not directory.is_dir() or directory.name=='m7-audit-source' or directory==out:continue
  for path in sorted(directory.rglob('*')):
   if path.is_file():artifacts.append(dict(path=str(path),bytes=path.stat().st_size,sha256=sha(path)))
 (out/'artifact-index.json').write_text(json.dumps(artifacts,indent=2)+'\n')
 receipt.update(status='verified',verification_scope='faithful-model implementation and diagnosis; robust model qualification failed',
  source_commits=dict(compiled_physics=pin['source_commit'],native=native_parent['source_commit'],audits=git(source,'rev-parse','HEAD'),operator_checkout=git(root/'source','rev-parse','HEAD'),user_linux=git(main,'rev-parse','HEAD')),
  registered_test_coverage=coverage,registered_test_count=len(coverage),native_checks=native['checks'],native_cases=successes,
  initial_rejections=rejections,audited_native_states=sum(x['states'] for x in successes),
  maximum_stage_residual=max(x['stage_residual'] for x in successes),transverse_checks=transverse['checks'],
  independent_replay_checks=replay['checks'],unchanged_runtime_images=len(pin['runtime_artifacts']),identical_campaign_runtime_context=True,frozen_inputs=frozen,
  equation_owners=equation_owners,stopped_flow_observations=covariance,
  physical_gates=dict(original_CLS_stability='failed_as_specified',inner_k_positivity='not_generally_invariant_M2',unforced_transition_passage='not_qualified',robust_M8_variant='not_adopted'),
  worker_allocation=dict(compiler_jobs=16,ctest_reserved_ranks=12,maximum_overlapping_numerical_ranks=15,numerical_threads=1,hard_limit=72),
  artifact_index=dict(path=str(out/'artifact-index.json'),sha256=sha(out/'artifact-index.json'),files=len(artifacts)))
except Exception as error:receipt['failure']=str(error)
receipt.update(checks=checks,finished_utc=datetime.now(timezone.utc).isoformat())
(out/'integration-audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt.get(k) for k in ('status','checks','registered_test_count','audited_native_states','failure')}),flush=True)
raise SystemExit(0 if receipt['status']=='verified' else 1)
