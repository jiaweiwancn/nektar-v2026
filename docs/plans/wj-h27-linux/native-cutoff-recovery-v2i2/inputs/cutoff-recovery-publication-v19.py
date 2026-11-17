from pathlib import Path
import json,shutil,tarfile
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/cutoff-recovery-v19';out=r/'evidence/native-cutoff-recovery-v2i2';out.mkdir(parents=True);raw=[]
def read(p):return json.loads(p.read_text(),parse_int=lambda token:-0.0 if token=='-0'else int(token))
assert read(o/'independent.json')['foundation_gate_pass']
def copy(p):
 target=out/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(identity(p))
drop={'coefficients','elements','input_history','current_advection','previous_advection','free_velocity_ids','k_free_ids','epsilon_free_ids','essential_flow_lift','essential_k_lift','essential_epsilon_lift','mms_bulk_body_rhs','mms_sampled_body_rhs'}
def compact(v,parent=''):
 if isinstance(v,dict):return{k:compact(value,k)for k,value in v.items()if not(v.get('schema')=='WJH27NativeCoupledStateV2'and k in drop)and not(parent in('profiles','owner_interfaces')and k=='rows')}
 if isinstance(v,list):return[compact(value,parent)for value in v]
 return v
for p in sorted(o.rglob('*')):
 if not p.is_file():continue
 raw.append(identity(p))
 if p.suffix in('.bin','.fld'):continue
 if p.suffix in('.xml','.opt','.jsonl')or(p.suffix=='.json'and p.stat().st_size<1024*1024):copy(p)
 elif p.suffix=='.json'and'.science-'not in p.name:
  value=compact(read(p));value['preserved_raw_artifact']=identity(p);value['publication_scope']='Compact publication; full signed native and independent owner/history/profile data retained on Linux.'
  target=out/'compact'/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);atomic(target,value)
for p in sorted((r/'inputs').iterdir()):
 if p.is_file()and('cutoff-recovery'in p.name or p.name=='retained_recovery_check-v19.py'):copy(p)
for p in (r/'tools-cutoff-recovery-v19/overlay.json',r/'tools-cutoff-recovery-v19/h27/retained_recovery_check.py',r/'tools-cutoff-recovery-v19/admission.md',r/'launches/cutoff-recovery-v19.log'):copy(p)
account=read(r/'accounting.json');labels=[key for key in account['entries']if key.startswith('cutoff-recovery-')]
for label in labels:
 for name in('request.json','receipt.json','command.log'):
  p=r/'operations'/label/name
  if p.exists():copy(p)
shutil.copy2(r/'accounting.json',out/'accounting-open-snapshot.json')
atomic(out/'assessment.json',dict(schema='WJH27RetainedRecoveryCheckpointV1',status='PASS_RETAINED_RECOVERY_CONTINUATION_ONLY',source_identity=read(o/'prelaunch.json')['source_identity'],
 full_H27_qualification=False,TASK3_READY=False,parent=identity(r/'evidence/native-cutoff-panels-v2i/inventory.json'),independent=identity(o/'independent.json'),
 account_closed=account['closed'],aggregate_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in account['entries'].values()),
 expected_negative_operation='cutoff-recovery-prior-source-native-v19 is a charged expected scientific-signature rejection, not a passed native import',
 commands_and_inputs='Exact argv, prelaunch input SHA and both scientific/retained history bindings are preserved in operation receipts and the independent report.'))
atomic(out/'inventory.json',dict(schema='WJH27RetainedRecoveryEvidenceV1',TASK3_READY=False,account_closed=False,
 raw_artifacts=list({v['path']:v for v in raw}.values()),phase_labels=labels,files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]))
archive=r/'evidence/native-cutoff-recovery-v2i2.tar.gz'
with tarfile.open(archive,'w:gz')as tar:tar.add(out,arcname=out.name)
print(json.dumps(dict(archive=identity(archive),status='PASS_RETAINED_RECOVERY_CONTINUATION_ONLY')))
