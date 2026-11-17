from pathlib import Path
import os,sys,json,shutil,tarfile,subprocess,re
from resources import atomic,identity
r=Path('/home/jay/build/wj-h27-coupling-r1-20261117');o=r/'science/cutoff-panels-v19';out=r/'evidence/native-cutoff-panels-v2i';raw=[]
def read(p):return json.loads(p.read_text(),parse_int=lambda token:-0.0 if token=='-0'else int(token))
assert read(o/'assessment.json')['status']=='PASS_CUTOFF_INTEGRATION_FOUNDATION_ONLY'
out.mkdir(parents=True)
def copy(p,target=None):
 target=out/p.relative_to(r)if target is None else out/target;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target);raw.append(identity(p))
drop={'coefficients','elements','input_history','current_advection','previous_advection','free_velocity_ids','k_free_ids','epsilon_free_ids','essential_flow_lift','essential_k_lift','essential_epsilon_lift','mms_bulk_body_rhs','mms_sampled_body_rhs'}
def compact(value,parent=''):
 if isinstance(value,dict):
  native=value.get('schema')=='WJH27NativeCoupledStateV2'
  return{k:compact(v,k)for k,v in value.items()if not(native and k in drop)and k not in ('best_unknowns','best_raw_residual')and not(parent in ('profiles','owner_interfaces')and k=='rows')}
 if isinstance(value,list):return[compact(v,parent)for v in value]
 return value
folders=[r/'science/physical-reference-r2-v17',r/'science/htrace-regression-v18',r/'science/radial-quadrature-v18',o,r/'trace-geometry-v17',r/'unit-cutoff-panels-v19']
for folder in folders:
 for p in sorted(folder.rglob('*')):
  if not p.is_file():continue
  raw.append(identity(p))
  if p.suffix in ('.bin','.npz','.fld','.elf','.a','.o','.gz')or'/.cache/'in str(p)or'matplotlib-cache'in str(p):continue
  if p.suffix in ('.xml','.opt','.py','.cpp','.h','.png','.svg','.jsonl','.md')or(p.suffix=='.json'and p.stat().st_size<1024*1024):copy(p)
  elif p.suffix=='.json'and'.science-'not in p.name:
   value=compact(read(p));value['preserved_raw_artifact']=identity(p);value['publication_scope']='Compact evidence; raw signed fields, histories and per-owner profile/interface rows remain at the exact Linux artifact identity.'
   target=out/'compact'/p.relative_to(r);target.parent.mkdir(parents=True,exist_ok=True);atomic(target,value)
for p in sorted((r/'science').glob('trace-geometry-*.json')):copy(p)
copy(r/'science/cutoff-panels-unit-v19.json')
for v in (18,19):
 for p in (r/('source-v'+str(v)+'.json'),r/'science'/('runtime-v'+str(v)+'.json')):copy(p)
 source=r/('source-v'+str(v))
 for p in sorted((source/'solvers/IncNavierStokesSolver/EquationSystems').glob('WJNative*')):
  if p.is_file():copy(p)
 for p in sorted((source/'solvers/IncNavierStokesSolver/Tests/RANS/h27').iterdir()):
  if p.is_file()and p.suffix in ('.py','.cpp','.cmake','.txt'):copy(p)
for p in sorted((r/'inputs').iterdir()):
 if p.is_file()and p.suffix in ('.py','.json','.md')and any(key in p.name for key in ('htrace','quadrature','trace-geometry','cutoff-panels','materialize-v18','materialize-v19')):copy(p)
for p in sorted((r/'launches').iterdir()):
 if p.is_file()and any(key in p.name for key in ('htrace','quadrature','trace-geometry','cutoff-panels','radial-reference')):copy(p)
account=read(r/'accounting.json');prior=read(r/'evidence/native-scientific-retention-v2h/accounting-open-snapshot.json')
labels=[key for key in account['entries']if key not in prior['entries']]
for label in labels:
 for name in ('request.json','receipt.json','command.log'):
  p=r/'operations'/label/name
  if p.exists():copy(p)
copy(r/'accounting.json','accounting-open-snapshot.json');copy(o/'assessment.json','assessment.json')
dependencies={}
for v in (18,19):
 rt=read(r/'science'/('runtime-v'+str(v)+'.json'))
 for executable in rt['executables']:
  elf=executable['executable'];assert identity(Path(elf['path']))==elf
  listing=subprocess.check_output([rt['loader']['path'],'--library-path',rt['libraries'],'--list',elf['path']],text=True)
  paths=set(re.findall(r'(?:=>\s+)?(/[^\s()]+)',listing));dependencies[Path(elf['path']).name+'-v'+str(v)]=dict(executable=elf,loader=identity(Path(rt['loader']['path'])),listing=listing,dependencies=[identity(Path(p))for p in sorted(paths)if Path(p).is_file()])
atomic(out/'runtime-dependencies.json',dependencies)
atomic(out/'inventory.json',dict(schema='WJH27CutoffPanelEvidenceInventoryV1',scope='Approved-v2 H-trace and moving-cutoff quadrature repair foundations, including all preserved failed/conditional-not-run outcomes; no H27-C/D or Task3 closure',
 TASK3_READY=False,account_closed=False,parent_inventory=identity(r/'evidence/native-scientific-retention-v2h/inventory.json'),
 phase_labels=labels,aggregate_completed_seconds_open_snapshot=sum(v.get('seconds',0)for v in account['entries'].values()),
 raw_artifacts=list({v['path']:v for v in raw}.values()),files=[identity(p)for p in sorted(out.rglob('*'))if p.is_file()]))
archive=r/'evidence/native-cutoff-panels-v2i.tar.gz'
with tarfile.open(archive,'w:gz')as tar:tar.add(out,arcname=out.name)
print(json.dumps(dict(archive=identity(archive),status=read(o/'assessment.json')['status'],TASK3_READY=False)))
