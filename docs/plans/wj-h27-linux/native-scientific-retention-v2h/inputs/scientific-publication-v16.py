from pathlib import Path
import json,shutil,tarfile,hashlib,tempfile,sys,math,subprocess,re
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/scientific-retention-v16';t=r/'tools-scientific-retention-v16b';out=r/'evidence/native-scientific-retention-v2h';raw=[]
def read(p):return json.loads(p.read_text(),parse_int=lambda token:-0.0 if token=='-0'else int(token))
assert read(o/'post-retention-check-v2.json')['status']=='PASS_BOUNDED_SCIENTIFIC_RETENTION_FOUNDATION_ONLY'
assert read(o/'retention-assessment-v2.json')['status']=='PASS_SCIENTIFIC_REDUCTION_AND_ROTATION_ONLY'
out.mkdir(parents=True)
def copy(p,target=None):
 relative=p.relative_to(r);target=out/relative if target is None else out/target
 target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(identity(p))
drop={'coefficients','elements','input_history','current_advection','previous_advection','free_velocity_ids','k_free_ids','epsilon_free_ids',
      'essential_flow_lift','essential_k_lift','essential_epsilon_lift','mms_bulk_body_rhs','mms_sampled_body_rhs'}
def compact(v,parent=''):
 if isinstance(v,dict):
  native=v.get('schema')=='WJH27NativeCoupledStateV2'
  return{k:compact(value,k)for k,value in v.items()if not(native and k in drop)and not(parent=='profiles'and k=='rows')}
 if isinstance(v,list):return[compact(value,parent)for value in v]
 return v
for p in sorted(o.rglob('*')):
 if not p.is_file():continue
 raw.append(identity(p))
 if p.suffix in ['.npz','.bin','.fld','.gz']:continue
 if p.suffix in ['.xml','.opt','.jsonl']or(p.suffix=='.json'and p.stat().st_size<1024*1024):copy(p)
 elif p.suffix=='.json'and'.science-'not in p.name:
  value=compact(read(p));value['preserved_raw_artifact']=identity(p)
  value['publication_scope']='Compact evidence; complete raw signed owner/coefficient/history/profile payloads remain on Linux.'
  target=out/'compact'/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);atomic(target,value)
# Independently bind all four recovery writes to the archived accepted science.
sys.path.insert(0,str(t/'h27'));import physical_cache_history_check as history
receipt=read(o/'protected-retention-v2/retention-receipt.json');archive=Path(receipt['archive']['path']);assert identity(archive)==receipt['archive']
recovery=[]
with tarfile.open(archive,'r:gz')as tar,tempfile.TemporaryDirectory(prefix='recovery-publication-',dir=o)as temp:
 for step,epoch in [(3,1),(6,2),(9,3),(12,3)]:
  rel='epoch'+str(epoch)+'.recovery-'+str(step)+'.bin';original=o/rel
  if original.exists():path=original
  else:
   member=tar.extractfile(rel);assert member is not None;path=Path(temp)/rel
   with path.open('xb')as f:f.write(member.read())
  sci='epoch'+str(epoch)+'.science-'+str(step)+'.json';record=json.loads(tar.extractfile(sci).read(),parse_int=lambda token:-0.0 if token=='-0'else int(token))
  states,binary=history.binary_history(path,record,.005,__import__('time').monotonic()+15)
  assert binary['step']==step and all(binary['hashes'][role]==record[key]for role,key in zip(history.HISTORY_NAMES,history.HASH_NAMES))
  expected=next(item for item in receipt['retained']+receipt['removed']if item['path']==str(original))
  assert identity(path)['sha256']==expected['sha256']
  recovery.append(dict(step=step,original=expected,archived_or_retained_input='retained'if original.exists()else'protected_archive',
    four_complete_history_hashes=binary['hashes'],accepted_scientific_state_sha256=record['state_sha256'],status='MET'))
atomic(out/'four-recovery-write-audit.json',dict(schema='WJH27FourScheduledRecoveryAuditV1',status='PASS_FOUR_SCHEDULED_WRITES',archive=identity(archive),writes=recovery,
    reader=identity(t/'h27/physical_cache_history_check.py'),scope='All four exact recovery binaries independently bind current/previous/explicit histories to their accepted scientific frames; latest9/12 separately cold-imported after actual rotation.'))
for folder in [t,r/'tools-scientific-retention-v16']:
 for p in sorted(folder.rglob('*')):
  if p.is_file()and p.suffix in ['.py','.md','.json','.csv']:copy(p)
for p in sorted((r/'inputs').glob('*')):
 if p.is_file()and('scientific' in p.name or p.name.startswith('v16b-')):copy(p)
for name in ['scientific-data-v16.log','scientific-retention-v16b.log']:copy(r/'launches'/name)
for name in ['source-v15.json','science/runtime-v15.json']:copy(r/name)
account=read(r/'accounting.json');labels=[key for key in account['entries']if key.startswith('scientific-')]
for label in labels:
 for name in ['request.json','receipt.json','command.log']:
  p=r/'operations'/label/name
  if p.exists():copy(p)
copy(r/'accounting.json','accounting-open-snapshot.json')
cases={}
for name in ['control','epoch1','epoch2','epoch3','reload9','reload12']:
 p=read(o/(name+'.result.json'));s=p['state'];logs=[readline for line in(o/(name+'.steps.jsonl')).read_text().splitlines()for readline in[json.loads(line)]]
 cases[name]={key:p[key]for key in ['construction_seconds','initialisation_seconds','complete_advancing_seconds','final_export_seconds','complete_step_seconds','scientific_frame_seconds','scientific_index_seconds','recovery_seconds']}
 cases[name].update(source_identity=s['source_identity'],signature=s['signature'],end_step=s['step'],end_time=s['time'],dt=s['dt'],
   newton_iterations=[row['nonlinear_iterations']for row in logs],krylov_iterations=[row['coupled_newton']['krylov_iterations']for row in logs],
   whole_operation_seconds=account['entries']['scientific-'+name+'-native-v16'if not name.startswith('reload')else'scientific-retained-'+name+'-native-v16b']['seconds'])
v=cases['control']['complete_step_seconds'];means=[math.fsum(v[i:i+3])/3 for i in [0,3,6,9]];N=math.ceil(505.02466521278404/.005)
conditional=dict(dt=.005,steps=N,initial_block_mean=means[0],warm_block_means=means[1:],largest_warm_mean=max(means[1:]),
    advancement_only_extrapolation_seconds=N*max(means[1:]),one_point_two_five_advancement_only_seconds=1.25*N*max(means[1:]),
    four_trajectory_advancement_only_seconds=4*N*max(means[1:]),
    scope='Unqualified short700-cell startup work-class extrapolation only. Not a complete per-trajectory/campaign forecast; selected h/p/q/dt, long-state costs, monitoring, statistics, restart and fixed costs remain unqualified.',production_forecast_qualified=False)
rt=read(r/'science/runtime-v15.json');elf=rt['executables'][0]['executable'];assert identity(Path(elf['path']))==elf
listing=subprocess.check_output([rt['loader']['path'],'--library-path',rt['libraries'],'--list',elf['path']],text=True)
paths=set(re.findall(r'(?:=>\s+)?(/[^\s()]+)',listing))
atomic(out/'runtime-dependencies.json',dict(binary=identity(Path(elf['path'])),loader=identity(Path(rt['loader']['path'])),listing=listing,
    dependencies=[identity(Path(p))for p in sorted(paths)if Path(p).is_file()]))
primary=Path('/media/jay/新加卷1/nektar-v2026');status=subprocess.check_output(['git','-C',str(primary),'status','--short'],text=True)
head=subprocess.check_output(['git','-C',str(primary),'rev-parse','HEAD'],text=True).strip();assert head=='f4b70baec10016dd098005bc5a5507b87e4a57bf' and status==' M docs/developer-guide/src/library\n'
source_dirs=[p for p in r.iterdir()if p.is_dir()and p.name.startswith(('source-v','build-v','runtime-v'))]
def bytes_in(folder):return sum(p.stat().st_size for p in folder.rglob('*')if p.is_file())
source_bytes=sum(bytes_in(p)for p in source_dirs);science_bytes=bytes_in(r/'science');free=shutil.disk_usage(r).free
assert source_bytes<8*1024**3 and science_bytes<8*1024**3 and free>15*1024**3
entries=[account['entries'][key]for key in labels if account['entries'][key]['status']!='RESERVED']
atomic(out/'resource-preservation.json',dict(primary=dict(path=str(primary),head=head,status=status,preserved=True),
    source_build_runtime_bytes=source_bytes,science_bytes=science_bytes,free_bytes=free,max_scientific_host_workers=max(v.get('peak_host_workers',0)for v in entries),
    max_scientific_tree_rss_bytes=max(v.get('peak_tree_rss_bytes',0)for v in entries),threads_all_six_one=True,rank_count=1,worker_ceiling=72,
    prior_evidence=identity(r/'evidence/native-current-data-offset-v2g2/inventory.json'),account_closed=False,
    aggregate_closed_prior_and_current_seconds=sum(v.get('seconds',0)for v in account['entries'].values()),campaign_ceiling=43200,evidence_protected=3600))
atomic(out/'assessment.json',dict(schema='WJH27ScientificRetentionCheckpointV1',status='PASS_BOUNDED_SCIENTIFIC_RETENTION_FOUNDATION_ONLY',
    active_model='FiniteNuIsotropicWallFullWJCollarV2',unapproved_v4_implemented=False,no_new_cpp_build=True,source_identity='47fb3597c4e56348361cfab6eadf9a34ba8d89027e57ed01153ca4d9c5ab10d1',
    reference_identity=history.po.TABLE_SHA,cases=cases,conditional_cost=conditional,
    reduction=identity(o/'retention-assessment-v2.json'),cold_imports=identity(o/'post-retention-check-v2.json'),
    all_four_recovery_write_audit=identity(out/'four-recovery-write-audit.json'),full_H27_qualification=False,H27_D_pass=False,TASK3_READY=False,
    reason='Only short persisted moments/owned rotation qualified; same-inlet physical reference, selected numerical recipe, required-duration pilots, monitoring and complete cost remain open.'))
atomic(out/'inventory.json',dict(scope='Source15 executable with V2 persisted statistics, independent moments, actual owned retention and cold imports. No H27-C/D or Task3 closure.',
    TASK3_READY=False,account_closed=False,raw_artifacts=list({v['path']:v for v in raw}.values()),files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]))
archive=r/'evidence/native-scientific-retention-v2h.tar.gz'
with tarfile.open(archive,'w:gz')as tar:tar.add(out,arcname='native-scientific-retention-v2h')
print(json.dumps(dict(archive=identity(archive),status='PASS_BOUNDED_SCIENTIFIC_RETENTION_FOUNDATION_ONLY',conditional_cost=conditional)))
