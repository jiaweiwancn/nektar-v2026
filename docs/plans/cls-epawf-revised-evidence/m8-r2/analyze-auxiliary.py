"""Read-only, immutable-prefix analysis of one completed native R2 solve record."""
import argparse,hashlib,json,os,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
for name in ('source','preflight','preflight-session','stage-session','trace','runtime-pin','output'):
 p.add_argument('--'+name,required=True,type=Path)
p.add_argument('--expected-commit',required=True)
p.add_argument('--record-index',type=int,default=0)
p.add_argument('--matrix-root',type=Path)
a=p.parse_args();a.source=a.source.resolve();a.output.mkdir(parents=True,exist_ok=False)
keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
def sha(path):
 digest=hashlib.sha256()
 with Path(path).open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
 return digest.hexdigest()
def canonical(value):return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def git(source,*args):return subprocess.check_output(['git','-c','gc.auto=0','-c','pack.threads=1','-c','checkout.workers=1','-C',str(source),*args],text=True).strip()
def save(path,value):
 temporary=path.with_suffix(path.suffix+'.tmp');temporary.write_text(json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n');temporary.replace(path)
def prefix_hash(path,count):
 digest=hashlib.sha256();remaining=count
 with Path(path).open('rb') as stream:
  while remaining:
   block=stream.read(min(1024*1024,remaining))
   if not block:raise ValueError('Native trace prefix was truncated')
   digest.update(block);remaining-=len(block)
 return digest.hexdigest()
helper=a.source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_stage_analysis.py'
inputs=(a.preflight,a.preflight_session,a.stage_session,a.runtime_pin,helper,Path(__file__))
receipt=dict(schema=1,success=False,started_utc=datetime.now(timezone.utc).isoformat(),native_calls=0,
 source_commit=git(a.source,'rev-parse','HEAD'),source_clean_start=not bool(git(a.source,'status','--porcelain')),
 numerical_threads={key:os.environ.get(key) for key in keys},cpu_affinity=sorted(os.sched_getaffinity(0)),
 command_argv=sys.argv,wrapper_sha256=sha(__file__),record_index=a.record_index,
 inputs={str(path.resolve()):dict(bytes=path.stat().st_size,sha256=sha(path)) for path in inputs},failures=[],
 scope='One completed native nonlinear record only; no final case protocol, physical advancement or target-status claim')
save(a.output/'receipt.json',receipt);started=time.monotonic();result=None
try:
 if receipt['source_commit']!=a.expected_commit or not receipt['source_clean_start']:raise ValueError('Analysis driver source is not the clean requested pin')
 if receipt['cpu_affinity']!=[48] or any(receipt['numerical_threads'][key]!='1' for key in keys):raise ValueError('Wrong CPU/thread allocation')
 if a.record_index<0:raise ValueError('Negative native record index')
 digest=hashlib.sha256();count=0;offset=0;line=None
 with a.trace.open('rb') as stream:
  for index in range(a.record_index+1):
   current=stream.readline()
   if not current or not current.endswith(b'\n'):raise ValueError('Requested native JSONL record is not fully published; no analysis attempted')
   offset=count;count+=len(current);digest.update(current);line=current
 snapshot=a.output/('native-record-'+str(a.record_index)+'.raw.jsonl');snapshot.write_bytes(line)
 receipt['active_trace_prefix']=dict(path=str(a.trace.resolve()),prefix_bytes=count,prefix_sha256=digest.hexdigest(),
  original_trace_bytes=a.trace.stat().st_size,selected_offset=offset,selected_bytes=len(line),selected_sha256=sha(snapshot),
  snapshot=str(snapshot),later_append_permitted=True)
 receipt['inputs'][str(snapshot.resolve())]=dict(bytes=snapshot.stat().st_size,sha256=sha(snapshot));save(a.output/'receipt.json',receipt)
 pin=json.loads(a.runtime_pin.read_text())
 if not pin['success'] or not pin['clean']:raise ValueError('Producing native runtime pin failed')
 producing=Path(pin['source_root'])
 if git(producing,'rev-parse','HEAD')!=pin['source_commit'] or git(producing,'status','--porcelain'):raise ValueError('Producing source no longer matches the runtime pin')
 solver_keys=[key for key in pin['runtime_artifacts'] if Path(key).name=='IncNavierStokesSolver']
 if len(solver_keys)!=1:raise ValueError('No unique producing solver image')
 solver=Path(pin['runtime_resolved_paths'][solver_keys[0]])
 if sha(solver)!=pin['runtime_artifacts'][solver_keys[0]]:raise ValueError('Producing native solver image changed')
 receipt['producing_source_commit']=pin['source_commit'];receipt['producing_solver']=dict(path=str(solver),sha256=sha(solver))
 receipt['producing_runtime_sections']={key:canonical(pin[key]) for key in ('source_inputs','runtime_artifacts','runtime_bindings','runtime_environment','python_packages','tool_images')}
 sys.dont_write_bytecode=True;sys.path.insert(0,str(helper.parent))
 import cls_epawf_m8_r2_stage_analysis as analysis
 imported={str(Path(module.__file__).resolve()):sha(module.__file__) for module in list(sys.modules.values()) if getattr(module,'__file__',None) and str(Path(module.__file__).resolve()).startswith(str(a.source))}
 receipt['imported_source_files_start']=imported
 raw=json.loads(a.preflight.read_text());record=json.loads(line)
 result=analysis.verify_auxiliary(raw,record,a.preflight_session,a.stage_session,a.source,a.matrix_root or a.trace.parent)
 result['input_receipt']=str(a.output/'receipt.json');save(a.output/'validation.json',result)
 receipt['success']=bool(result['analysis_success']);receipt['selected_alpha']=record['alpha']
 for path,item in result['dense_matrix_files'].items():receipt['inputs'][path]=item
 for path,expected in imported.items():
  if sha(path)!=expected:raise ValueError('Imported driver changed: '+path)
except Exception as error:
 receipt['success']=False;receipt['failures'].append(dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc()))
for path,expected in receipt['inputs'].items():
 if sha(path)!=expected['sha256'] or Path(path).stat().st_size!=expected['bytes']:
  receipt['success']=False;receipt['failures'].append(dict(type='ProvenanceError',message='Retained input changed: '+path))
if 'active_trace_prefix' in receipt:
 prefix=receipt['active_trace_prefix']
 try:
  prefix['prefix_sha256_end']=prefix_hash(a.trace,prefix['prefix_bytes']);prefix['trace_bytes_end']=a.trace.stat().st_size
  if prefix['prefix_sha256_end']!=prefix['prefix_sha256']:raise ValueError('Analyzed native prefix changed while a later solve ran')
 except Exception as error:receipt['success']=False;receipt['failures'].append(dict(type='ProvenanceError',message=str(error)))
receipt.update(source_commit_end=git(a.source,'rev-parse','HEAD'),source_clean_end=not bool(git(a.source,'status','--porcelain')),seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat())
if receipt['source_commit_end']!=receipt['source_commit'] or not receipt['source_clean_end']:
 receipt['success']=False;receipt['failures'].append(dict(type='ProvenanceError',message='Analysis source changed'))
save(a.output/'receipt.json',receipt)
print(json.dumps({key:receipt.get(key) for key in ('success','record_index','selected_alpha','seconds','failures')}),flush=True)
raise SystemExit(0 if receipt['success'] else 1)
