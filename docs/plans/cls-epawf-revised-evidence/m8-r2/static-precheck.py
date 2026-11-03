"""Read-only static adapter and exact affine weak-continuity value precheck."""
import argparse,hashlib,json,os,subprocess,sys,time,traceback
from datetime import datetime,timezone
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--source',required=True,type=Path);p.add_argument('--expected-commit',required=True)
p.add_argument('--case',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
a=p.parse_args();source=a.source.resolve();a.output.mkdir(exist_ok=False)
raw_path,session=a.case/'case.evolution.json',a.case/'case.xml'
helper=source/'solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m8_r2_stage_analysis.py'
keys=('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')
sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()
git=lambda *args:subprocess.check_output(['git','-c','gc.auto=0','-c','pack.threads=1','-C',str(source),*args],text=True).strip()
started=time.monotonic();receipt=dict(schema=1,success=False,started_utc=datetime.now(timezone.utc).isoformat(),
 source_commit=git('rev-parse','HEAD'),source_clean_start=not bool(git('status','--porcelain')),
 role='Read-only constructor, actual static values and separately exact affine weak-continuity values',
 native_calls=0,weak_residual_calls=0,cpu_affinity=sorted(os.sched_getaffinity(0)),numerical_threads={k:os.environ.get(k) for k in keys},
 command_argv=sys.argv,wrapper_sha256=sha(__file__),inputs={str(path):dict(bytes=path.stat().st_size,sha256=sha(path)) for path in (raw_path,session,helper)},failures=[])
def save():
 temporary=a.output/'receipt.json.tmp';temporary.write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n');temporary.replace(a.output/'receipt.json')
save()
try:
 assert receipt['source_commit']==a.expected_commit and receipt['source_clean_start']
 assert receipt['cpu_affinity']==[48] and all(receipt['numerical_threads'][k]=='1' for k in keys)
 sys.dont_write_bytecode=True;sys.path.insert(0,str(helper.parent))
 import cls_epawf_m8_r2_stage_analysis as analysis
 raw=json.loads(raw_path.read_text());analyzer=analysis.Analyzer(raw,session,source)
 receipt.update(constructor_passed=True,unknowns=analyzer.size,static_coverage=analyzer.static_coverage,
  static_group_count=len(analyzer.static),matching_epsilon_rows=len(analyzer.trace_rows),external_fixed_rows=len(analyzer.external),
  native_source_anchors=len(analyzer.source_anchors),native_trace_anchors=len(analyzer.trace_anchors))
 names={group['name'] for group in analyzer.static}
 observations=[row for row in raw['initial_scalar_domain']['observations'] if row['inventory'] in names]
 effective=analyzer.external_lift(raw['initial']);analyzer.observations(observations,effective,effective)
 receipt['actual_static_field_values_checked']=len(observations)
 receipt['independent_affine_continuity_values']=analyzer.affine_continuity_values(raw['fixed_reference'],raw['raw_reference_loads'])
 receipt['independent_affine_continuity_assembly']=analyzer.continuity_metadata
 assert receipt['independent_affine_continuity_values']['fixed_1e_minus7_gate_met']
 receipt.update(checks=analyzer.audit.checks,errors=analyzer.audit.errors,
  scope='No full weak reassembly or native call; initial fixed-reference weak continuity/gauge values only. No current-Jacobian or stage convergence claim.')
 receipt['imported_source_files']={str(Path(module.__file__).resolve()):sha(module.__file__) for module in list(sys.modules.values()) if getattr(module,'__file__',None) and str(Path(module.__file__).resolve()).startswith(str(source))}
 receipt['success']=True
except Exception as e:receipt['failures'].append(dict(type=type(e).__name__,message=str(e),traceback=traceback.format_exc()))
for path,item in receipt['inputs'].items():
 if sha(path)!=item['sha256'] or Path(path).stat().st_size!=item['bytes']:receipt['success']=False;receipt['failures'].append(dict(type='ProvenanceError',message='Input changed: '+path))
receipt.update(source_commit_end=git('rev-parse','HEAD'),source_clean_end=not bool(git('status','--porcelain')),seconds=time.monotonic()-started,finished_utc=datetime.now(timezone.utc).isoformat())
if receipt['source_commit_end']!=receipt['source_commit'] or not receipt['source_clean_end']:receipt['success']=False;receipt['failures'].append(dict(type='ProvenanceError',message='Source changed'))
save();print(json.dumps({key:receipt.get(key) for key in ('success','constructor_passed','unknowns','checks','seconds','static_coverage','errors','independent_affine_continuity_values','failures')}),flush=True)
raise SystemExit(0 if receipt['success'] else 1)
