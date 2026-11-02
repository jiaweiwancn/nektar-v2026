import argparse, hashlib, json, os, platform, re, subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('attempt');args=p.parse_args()
root=Path('/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised')
main=Path('/media/jay/新加卷/nektar-v2026');source=root/'source';build=root/'build';out=root/args.attempt

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def git(path,*args):return subprocess.check_output(['git','-C',str(path),*args],text=True).strip()
def read(path):return json.loads(Path(path).read_text())

direct=read(out/'oracle/validation.json');receipt=read(out/'receipt.json')
assert receipt['status']=='passed' and len(receipt['commands']) in (3,4) and all(c['exit_code']==0 for c in receipt['commands'])
assert direct['verification_status']=='passed' and direct['milestone_closure']
assert git(source,'rev-parse','HEAD')==direct['source_commit']==receipt['source_commit']
assert not git(source,'status','--porcelain') and not git(main,'status','--porcelain')
assert '100% tests passed, 0 tests failed out of 22' in (out/'suite.log').read_text()

def current_report(pattern):
 found=[(f,read(f)) for f in build.glob(pattern) if read(f)['source_commit']==direct['source_commit']]
 assert len(found)==1,(pattern,len(found))
 f,r=found[0];assert r['verification_status']=='passed' and r['source_status']==[]
 return f,r
ctest_path,ctest=current_report('cls-epawf-m3-*/validation.json')
assert ctest['milestone_closure'] and ctest['checks']==direct['checks']
native_sha=sha(out/'oracle/patch.json')
previous_direct=root/'m3-attempt-016/oracle'
prior=read(previous_direct/'validation.json')
assert prior['verification_status']=='passed' and prior['milestone_closure']
assert sha(previous_direct/'patch.json')==native_sha
assert all(c.get('binary_sha256')==sha(build/'IncNavierStokesSolver') for c in prior['commands'][1:] if c['command'][0].endswith('/IncNavierStokesSolver'))
assert all(v==prior['source_sha256'][f] for f,v in direct['source_sha256'].items()
           if not f.endswith('/validate_cls_epawf_m3.py'))

assert native_sha==ctest['artifacts_sha256']['patch.json']
comparison={'numeric_values':0,'max_absolute_difference':0.,'max_scaled_difference':0.,'worst_path':''}
def compare(a,b,path=''):
 if isinstance(a,(int,float)) and not isinstance(a,bool) and isinstance(b,(int,float)) and not isinstance(b,bool):
  error=abs(a-b);scaled=error/max(1.,abs(a),abs(b));comparison['numeric_values']+=1
  comparison['max_absolute_difference']=max(comparison['max_absolute_difference'],error)
  if scaled>comparison['max_scaled_difference']:comparison.update(max_scaled_difference=scaled,worst_path=path)
  assert scaled<1e-12,(path,a,b)
 elif isinstance(a,dict):
  assert isinstance(b,dict) and a.keys()==b.keys()
  for key in a:compare(a[key],b[key],path+'/'+key)
 elif isinstance(a,list):
  assert isinstance(b,list) and len(a)==len(b)
  for key,(u,v) in enumerate(zip(a,b)):compare(u,v,path+'/'+str(key))
 else:assert a==b,(path,a,b)
compare(read(out/'oracle/patch.json'),read(build/'rans-wall-native/patch.json'))
comparison['byte_identical']=native_sha==sha(build/'rans-wall-native/patch.json')
comparison['standalone_native_sha256']=sha(build/'rans-wall-native/patch.json')

old_m1=read(source/'docs/plans/cls-epawf-revised-evidence/m1/verification.json')
_,m1=current_report('cls-epawf-m1-*/validation.json')
unchanged_m1={n:m1['artifacts_sha256'][n]==old_m1['artifacts_sha256'][n] for n in ('kernels.json','constitutive.json')}
assert all(unchanged_m1.values())
m2_path,m2=current_report('cls-epawf-m2-*/validation.json')
old_m2=read(source/'docs/plans/cls-epawf-revised-evidence/m2/verification.json')
assert m2['artifacts_sha256']['profiles.json']==old_m2['artifacts_sha256']['profiles.json']
names=['IncNavierStokesSolver','RANSKernelProbe','RANSWallLayerProbe','RANSWallGeometryProbe',
       'RANSObservationKernelProbe','RANSStressTangentAudit','RANSConstitutiveProbe',
       'RANSEquilibriumProfilesProbe','RANSEquilibriumPatchProbe','RANSFieldAudit','RANSWallFieldAudit','RANSObservationAudit']
libraries={};ldd_hashes={}
for name in names:
 text=subprocess.check_output(['ldd',str(build/name)],text=True);assert 'not found' not in text
 f=out/(name+'.ldd');f.write_text(text);ldd_hashes[name]=sha(f)
 for line in text.splitlines():
  for token in line.split():
   if token.startswith('/'):
    path=Path(token).resolve();assert path.is_file();libraries[str(path)]=sha(path)
cache=Path('/home/jay/build/nektar-v2026-mpi/CMakeCache.txt').read_text()
configuration={key:re.search(r'^'+key+r':[^=]+=(.*)$',cache,re.M).group(1)
 for key in ('NEKTAR_USE_MPI','CMAKE_BUILD_TYPE','CMAKE_CXX_COMPILER','CMAKE_CXX_FLAGS','CMAKE_CXX_FLAGS_RELEASE')}
assert configuration['NEKTAR_USE_MPI']=='ON'
assert str((build/'geometry-lib/libSpatialDomains.so.5.9.0').resolve()) in libraries
previous=read(source/'docs/plans/cls-epawf-revised-evidence/m2/integration.json')['dynamic_libraries_sha256']
assert libraries==previous
integration={'schema':1,'status':'passed','source_commit':direct['source_commit'],
 'source_version':(source/'VERSION').read_text().strip(),'installation_configuration':configuration,
 'machine':platform.machine(),'current_suite_passed':22,'m2_suite_passed':19,'m1_suite_passed':17,'exact_base_suite_passed':16,
 'registered_oracle_and_native_reports_numerically_equivalent':True,'native_report_comparison':comparison,'native_artifact_sha256':native_sha,
 'unchanged_native_results_from_prior_direct_qualification':True,'prior_direct_source_commit':prior['source_commit'],
 'ctest_verification_path':str(ctest_path),'ctest_verification_sha256':sha(ctest_path),
 'ctest_verification_checks':ctest['checks'],'ctest_m2_verification_path':str(m2_path),
 'unchanged_M1_native_reports':unchanged_m1,'unchanged_M2_native_profiles':True,
 'dynamic_libraries_unchanged_from_M2':True,
 'unchanged_main_checkout':{'revision':git(main,'rev-parse','HEAD'),'status':[]},'isolated_source_status':[],
 'binary_sha256':{n:sha(build/n) for n in names},'static_kernel_archive_sha256':sha(build/'libRANSKernels.a'),
 'dynamic_libraries_sha256':libraries,'ldd_log_sha256':ldd_hashes,
 'runtime_LD_LIBRARY_PATH':os.environ.get('LD_LIBRARY_PATH',''),
 'probe_compile_flags':(build/'CMakeFiles/RANSEquilibriumPatchProbe.dir/flags.make').read_text(),
 'm3_solver_sources_sha256':{str(f.relative_to(source)):sha(f) for f in (source/'solvers/IncNavierStokesSolver/EquationSystems').glob('RANSWall*.cpp')},
 'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
 'cmake':subprocess.check_output(['cmake','--version'],text=True).splitlines()[0],
 'mpi':subprocess.check_output(['mpirun','--version'],text=True).splitlines()[0],
 'source_manifest_verified':all(sha(source/f)==h for f,h in direct['source_sha256'].items()),
 'test_inventory':subprocess.check_output(['ctest','--test-dir',str(build),'-N'],text=True).splitlines(),
 'collector_sha256':sha(__file__)}
assert integration['source_manifest_verified']
raw_index={str(f):sha(f) for directory in {(out/'oracle').resolve(),ctest_path.parent.resolve(),build/'rans-wall-native',build/'rans-wall-stage'}
           for f in directory.rglob('*') if f.is_file()}
(out/'artifact-index.json').write_text(json.dumps(raw_index,indent=2)+'\n')
integration['raw_artifact_index']={'path':str(out/'artifact-index.json'),'sha256':sha(out/'artifact-index.json'),'files':len(raw_index)}
integration['direct_session_sha256']={str(f.relative_to(out/'oracle')):sha(f) for f in (out/'oracle').rglob('*.xml')}
integration['solver_compile_flags']=(build/'CMakeFiles/IncNavierStokesSolver.dir/flags.make').read_text()
(out/'integration.json').write_text(json.dumps(integration,indent=2)+'\n')
attempts={'schema':1,'current_attempt':args.attempt,'attempts':[]}
for directory in sorted(root.glob('m3-attempt-*')):
 f=directory/'receipt.json'
 if not f.exists():continue
 r=read(f)
 item={'path':str(directory),'receipt_sha256':sha(f),'source_commit':r['source_commit'],
       'status':r['status'],'commands':r['commands']}
 v=directory/'oracle/validation.json'
 if v.exists():
  validation=read(v);item.update(oracle_status=validation['verification_status'],oracle_error=validation.get('error'))
 attempts['attempts'].append(item)
attempts['exploratory_diagnostics']=[{'path':str(d),'files_sha256':{str(f.relative_to(d)):sha(f) for f in d.rglob('*') if f.is_file()}}
                                    for d in sorted(root.glob('m3-diagnosis-*'))]
(out/'attempts.json').write_text(json.dumps(attempts,indent=2)+'\n')
print(json.dumps({'status':'passed','source_commit':direct['source_commit'],'binaries':len(names),'libraries':len(libraries),
 'checks':direct['checks'],'source_pins':len(direct['source_sha256']),
 'command_seconds':{str(i):c['seconds'] for i,c in enumerate(receipt['commands'])}},indent=2))
