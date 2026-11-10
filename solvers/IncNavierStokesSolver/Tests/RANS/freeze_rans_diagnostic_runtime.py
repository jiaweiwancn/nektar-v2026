#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Freeze a new probe runtime; preserve the trajectory and prior probe runtimes."""
import argparse,os,shutil,subprocess
from pathlib import Path
from rans_request_manifest import sha,runtime_identity
from rans_qualification import THREADS,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output-name',default='runtime-resolution-v3');p.add_argument('--parent-name',default='runtime-resolution-v2');a=p.parse_args();r=a.root.resolve();s=r/'source';out=r/a.output_name;assert not out.exists();shutil.copytree(r/a.parent_name,out)
    for name in ('source-build.json','lineage.json','actual-loaded-libraries.json'):
        if (out/name).exists():(out/name).rename(out/('inherited-'+a.parent_name+'-'+name))
    shutil.copy2(r/'build/RANSPreservationProbe',out/'RANSPreservationProbe');os.environ.update(THREADS,LD_LIBRARY_PATH=str(out/'lib'))
    names=[]
    for folder in ('library/SolverUtils/Turbulence','solvers/IncNavierStokesSolver/EquationSystems'):
        names.extend(p for p in (s/folder).iterdir() if p.suffix in ('.h','.cpp') and ('RANS' in p.name or folder.startswith('library')))
    names.append(s/'library/UnitTests/SolverUtils/RANSPreservationProbe.cpp')
    revision=subprocess.check_output(['git','-C',str(s),'rev-parse','HEAD'],text=True).strip()
    write_json(out/'source-build.json',dict(schema='ProbeOnlyRuntimeV1',source_revision=revision,source_file_sha256={str(p.relative_to(s)):sha(p) for p in sorted(names)},cmake_cache_sha256=sha(r/'build/CMakeCache.txt'),compiler_flags='-O3 -DNDEBUG -fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS',build_command=['cmake','--build',str(r/'build'),'--parallel','2','--target','RANSPreservationProbe'],runtime=runtime_identity(out/'RANSPreservationProbe'),trajectory_solver_sha256=sha(r/'runtime/IncNavierStokesSolver'),copied_solver_sha256=sha(out/'IncNavierStokesSolver'),scope='only this new directory receives a rebuilt read-only diagnostic probe; all evolution executables and existing runtimes remain unchanged',threads=THREADS))
    print(revision,sha(out/'RANSPreservationProbe'));return 0


if __name__=='__main__':raise SystemExit(main())
