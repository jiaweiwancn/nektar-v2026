#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze an isolated Linux build without changing any installation/runtime."""
import argparse,os,platform,shutil,subprocess,sys
from pathlib import Path
from rans_qualification import THREADS,write_json
from run_cls_epawf_picard import sha
from run_cls_epawf_pimple import identity
from continue_cls_epawf_pimple import freeze

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--build',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();s=a.source.resolve();b=a.build.resolve();out=a.output.resolve();assert platform.system()=='Linux';os.environ.update(THREADS)
 git=lambda *args:subprocess.check_output(['git','-c','core.filemode=false','-C',str(s),*args],text=True).strip()
 assert not git('status','--porcelain'),'source contents must be clean before freezing';revision=git('rev-parse','HEAD');solver,pinned=freeze(b/'IncNavierStokesSolver',out)
 probes={}
 for name in ('RANSPreservationProbe','RANSNewModelsProbe','RANSNativeDiagnosticChecks'):
  shutil.copy2(b/name,out/name);probes[name]=identity(out/name)
 files={}
 for directory in ('library/SolverUtils/Turbulence','library/UnitTests/SolverUtils','solvers/IncNavierStokesSolver'):
  for f in (s/directory).rglob('*'):
   if f.is_file() and (f.suffix in ('.cpp','.h') or f.name=='CMakeLists.txt'):files[str(f.relative_to(s))]=sha(f)
 record=dict(revision=revision,compiled_cpp_revision=revision,source_status='',source_files=files,solver=identity(solver),probes=probes,ranks=1,threads=THREADS,compiler=subprocess.check_output(['c++','--version'],text=True),os=subprocess.check_output(['uname','-a'],text=True),cmake_cache_sha256=sha(b/'CMakeCache.txt'),build_flags=(b/'CMakeFiles/IncNavierStokesSolver.dir/flags.make').read_text())
 write_json(out/'source-build.json',record);print(revision,pinned['sha256']);return 0
if __name__=='__main__':sys.exit(main())
