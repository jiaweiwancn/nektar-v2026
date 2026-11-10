#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Content-level mathematical and history provenance for valid evidence reuse."""
import argparse,ast,hashlib,json,subprocess
from pathlib import Path
from rans_qualification import write_json

def digest(b):return hashlib.sha256(b).hexdigest()
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);a=p.parse_args();r=a.root.resolve();s=r/'source';baseline='3b0961e12142d7eca56d75962feebfffa34630cc'
 git=lambda *args:subprocess.check_output(['git','-C',str(s),'-c','core.filemode=false',*args])
 files=[str(p.relative_to(s)) for p in (s/'library/SolverUtils/Turbulence').glob('*') if p.suffix in ('.cpp','.h')]
 base='solvers/IncNavierStokesSolver/EquationSystems/'
 files += [base+n for n in ('RANSWallEvolutionPimple.cpp','RANSWallEvolutionPicard.cpp','RANSWallEvolutionEvaluation.cpp','RANSWallEvolutionNonpressure.cpp','RANSWallEvolutionContact.cpp','RANSWallEvolutionPatch.cpp','RANSWallGlobalField.cpp','RANSWallTimeIntegration.cpp','RANSHistory.cpp','RANSWallPhysicsSignature.h','RANSWallSolveDiagnostics.h')]
 checks=[]
 for name in sorted(files):
  old=git('show',baseline+':'+name);new=(s/name).read_bytes();assert old==new,name;checks.append(dict(file=name,sha256=digest(new),unchanged=True))
 oracles={}
 for name in ('rans_preservation_oracle.py','validate_a93_a95_e94e95.py','a93_a95_e94e95_note_wall_oracle.py'):
  path='solvers/IncNavierStokesSolver/Tests/RANS/'+name;old=git('show',baseline+':'+path);new=(s/path).read_bytes();assert old==new,path;oracles[name]=dict(file_sha256=digest(new),semantic_AST_sha256=digest(ast.dump(ast.parse(new),include_attributes=False).encode()))
 specs={}
 for path in (s/'docs/literature').glob('*A93*A95*.tex'):specs[path.name]=digest(path.read_bytes())
 output=dict(status='pass',baseline=baseline,source_revision=git('rev-parse','HEAD').decode().strip(),governing_files=checks,independent_oracles=oracles,specifications=specs,
  reuse_scope='unchanged constitutive/wall kernels and analytical identities; unchanged evolution/history/cache/rollback implementation. New model-aware term export and diagnostic eligibility receive fresh targeted controls. No claim that new observation code is identical.',
  historical_epsilon_orders='eight observations retain original unresolved classifications')
 path=r/'mathematical-source-reuse.json';assert not path.exists();write_json(path,output);print(len(checks));return 0
if __name__=='__main__':raise SystemExit(main())
