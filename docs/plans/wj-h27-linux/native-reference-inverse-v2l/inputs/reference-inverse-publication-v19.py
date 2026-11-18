from pathlib import Path
import json,shutil,tarfile,subprocess
from resources import identity,atomic
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/reference-inverse-v19';t=r/'reference-inverse-project-v19';b=r/'build-reference-inverse-v19';out=r/'evidence/native-reference-inverse-v2l';out.mkdir();raw=[]
def read(p):return json.loads(p.read_text())
def copy(p):
 target=out/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(identity(p))
j=read(o/'independent.json');a=read(o/'reuse-audit.json');assert j['inverse_guard_gate_pass']is True and a['status']=='PASS_REUSED_NUMERICAL_FOUNDATIONS_ONLY'
project=read(t/'manifest.json')
for f in project['libraries']:assert identity(Path(f['path']))==f
for p in sorted(o.rglob('*')):
 if p.is_file():copy(p)
for p in sorted(t.rglob('*')):
 if p.is_file():copy(p)
for p in sorted((r/'inputs').glob('*')):
 if p.is_file()and(p.name.startswith('reference-inverse-')or p.name=='verification-reuse-audit-v19.py'):copy(p)
for n in ('compile_commands.json','CMakeFiles/WJH27ReferenceInverseProbe.dir/link.txt','CMakeFiles/WJH27ReferenceInverseProbe.dir/flags.make'):
 copy(b/n)
for p in(r/'science/runtime-v19.json',r/'launches/reference-inverse-v19.log'):copy(p)
account=read(r/'accounting.json');labels=[k for k in account['entries']if k.startswith('reference-inverse-')];entries={k:account['entries'][k]for k in labels}
for label in labels:
 for n in('request.json','receipt.json','command.log'):
  p=r/'operations'/label/n
  if p.exists():copy(p)
shutil.copy2(r/'accounting.json',out/'accounting-open-snapshot.json')
old=Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime');probe=identity(b/'WJH27ReferenceInverseProbe');loader=identity(old/'glibc/ld-linux-x86-64.so.2')
atomic(out/'probe-runtime.json',dict(probe=probe,loader=loader,library_path=str(old/'glibc')+':'+str(old/'geometry')+':/home/jay/.local/opt/nektar-v2026-mpi/lib64',static_libraries=project['libraries'],source=project['source'],selected_solver_runtime=identity(r/'science/runtime-v19.json'),selected_application_rebuilt=False))
def size(p):return sum(v.stat().st_size for v in p.rglob('*')if v.is_file())
source_bytes=sum(size(p)for p in r.iterdir()if p.is_dir()and p.name.startswith(('source-v','build-v','runtime-v')))+size(t)+size(b);science_bytes=size(r/'science');free=shutil.disk_usage(r).free
assert source_bytes<8*1024**3 and science_bytes<8*1024**3 and free>15*1024**3
resources=dict(source_build_runtime_and_thin_probe_bytes=source_bytes,science_bytes=science_bytes,free_bytes=free,thin_project_and_build_bytes=size(t)+size(b),account_closed=False,campaign_ceiling=43200,evidence_protected=3600,aggregate_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in account['entries'].values()),phase_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in entries.values()),phase_max_host_workers=max(v.get('peak_host_workers',0)for v in entries.values()),phase_max_tree_rss_bytes=max(v.get('peak_tree_rss_bytes',0)for v in entries.values()),ranks=1,all_six_thread_limits_one=True,compiler_parallelism=4)
assessment=dict(schema='WJH27ReferenceInverseCheckpointV1',status='PASS_FOCUSED_INVERSE_AND_REUSE_FOUNDATIONS_ONLY',full_H27_qualification=False,TASK3_READY=False,active_model='FiniteNuIsotropicWallFullWJCollarV2',v4_implemented=False,parent='d205335c80458009387430d2fa54f3c3e39ef560',source_identity=a['source_identity'],native=identity(o/'native.json'),probe=probe,independent=identity(o/'independent.json'),reuse_audit=identity(o/'reuse-audit.json'),maximum_scaled_errors=j['maximum_scaled_errors'],roundtrip_queries=129,independent_radius_queries=129,invalid_radius_controls=7,synthetic_mutations=4,default_identity_rejections=4,constructor_guard_rejections=4,legacy_foundations=a['legacy_reports'],resources=resources,scope='New source19 inverse guard checks plus source-bound reuse of unchanged mixed/gauge/Q3 and constitutive/admissibility tests. No full coupled/physical/operational/cost closure or production recipe selection.')
atomic(out/'assessment.json',assessment)
primary='/media/jay/新加卷1/nektar-v2026';head=subprocess.check_output(['git','-C',primary,'rev-parse','HEAD'],text=True).strip();status=subprocess.check_output(['git','-C',primary,'status','--porcelain'],text=True);assert head=='f4b70baec10016dd098005bc5a5507b87e4a57bf'and status==' M docs/developer-guide/src/library\n';atomic(out/'preserved-checkout.json',dict(path=primary,head=head,status=status,preserved=True))
atomic(out/'inventory.json',dict(schema='WJH27ReferenceInverseEvidenceV1',TASK3_READY=False,account_closed=False,parent=identity(r/'evidence/native-cutoff-owner-transition-v2k/inventory.json'),raw_artifacts=raw+[probe,loader]+project['libraries'],phase_labels=labels,files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]))
archive=r/'evidence/native-reference-inverse-v2l.tar.gz'
with tarfile.open(archive,'w:gz')as f:f.add(out,arcname=out.name)
print(json.dumps(dict(archive=identity(archive),status=assessment['status'],resources=resources)))
