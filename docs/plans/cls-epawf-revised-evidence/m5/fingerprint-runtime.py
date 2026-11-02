#!/usr/bin/env python3
"""Record exact M5 qualification inputs and linked Linux runtime artifacts."""
import argparse,hashlib,json,os,subprocess
from pathlib import Path
os.environ.update({k:'1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')})
p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--compare',type=Path);p.add_argument('--allow-replayed-probes',action='store_true');a=p.parse_args()
root=a.root;source=root/'source';build=root/'build'
def git(*args):return subprocess.check_output(['git','-c','gc.auto=0','-C',str(source),*args],text=True).strip()
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
assert not git('status','--porcelain')
paths=git('ls-tree','-r','--name-only','HEAD','library/SolverUtils/Turbulence','library/UnitTests/SolverUtils','solvers/IncNavierStokesSolver').splitlines()
inputs={name:digest(source/name) for name in paths if Path(name).suffix in ('.cpp','.h','.py','.xml','.tst')}
programs=['IncNavierStokesSolver','RANSCurvedGeometryProbe','RANSCurvedWallProbe','RANSEquilibriumEvolutionProbe','RANSEquilibriumProfilesProbe']
artifacts={str(build/name):digest(build/name) for name in programs}
env=dict(os.environ);env['LD_LIBRARY_PATH']=str(build/'geometry-lib')+':'+env.get('LD_LIBRARY_PATH','')
linked=subprocess.check_output(['ldd',str(build/'IncNavierStokesSolver')],text=True,env=env)
for line in linked.splitlines():
 if '=>' in line:
  name=line.split('=>',1)[1].strip().split(' ',1)[0]
  if name.startswith('/') and Path(name).is_file():artifacts[name]=digest(Path(name))
result=dict(schema=1,source_commit=git('rev-parse','HEAD'),clean=True,source_inputs=inputs,runtime_artifacts=artifacts,ldd=linked)
if a.compare:
 previous=json.loads(a.compare.read_text());changed={kind:[name for name,value in previous[kind].items() if result[kind].get(name)!=value] for kind in ('source_inputs','runtime_artifacts')}
 allowed_sources={'library/UnitTests/SolverUtils/RANSEquilibriumEvolutionProbe.cpp','library/UnitTests/SolverUtils/RANSCurvedWallProbe.cpp'}
 allowed_binaries={str(build/'RANSEquilibriumEvolutionProbe'),str(build/'RANSCurvedWallProbe')}
 if any(changed.values()):
  assert a.allow_replayed_probes
  assert set(changed['source_inputs'])<=allowed_sources and set(changed['runtime_artifacts'])<=allowed_binaries,changed
  for name in changed['source_inputs']:
   original=subprocess.check_output(['git','-C',str(source),'show',previous['source_commit']+':'+name],text=True)
   updated=(source/name).read_text()
   updated='\n'.join(line for line in updated.splitlines() if not line.strip().startswith('// SessionReader prepends optimiser sidecars') and not line.strip().startswith('// the original XML as the clone'))+'\n'
   assert updated==original.replace('auto filenames=world->GetFilenames();','std::vector<std::string> filenames{argv[1]};'),name
 result['identical_qualification_inputs_and_runtime']=not any(changed.values())
 result['identical_equations_and_solver_runtime']=True
 result['probe_session_filename_only_changes']=changed
 result['compared_source_commit']=previous['source_commit'];result['added_inputs']=sorted(set(inputs)-set(previous['source_inputs']))

a.output.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(source_commit=result['source_commit'],inputs=len(inputs),runtime_files=len(artifacts),identical=result.get('identical_qualification_inputs_and_runtime'))))
