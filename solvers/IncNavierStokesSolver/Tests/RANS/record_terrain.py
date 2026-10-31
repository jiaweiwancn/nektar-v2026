#!/usr/bin/env python3
"""Collect compact Linux provenance and passed terrain acceptance receipts."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--suite',type=Path,action='append',required=True,help='Passed evidence.json; may be repeated')
    parser.add_argument('--ctest',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Collect authoritative build/solver evidence on Linux')
    source = Path(__file__).resolve().parents[4]
    build = args.build.resolve()
    suites = {}
    for path in args.suite:
        value = json.loads(path.read_text())
        if value['status'] != 'PASS' or not value['cases'] or not all(v['passed'] for v in value['cases'].values()):
            raise RuntimeError('An acceptance suite is incomplete: '+str(path))
        suites[value['suite']] = value
    if args.ctest and '100% tests passed, 0 tests failed' not in args.ctest.read_text():
        raise RuntimeError('CTest did not pass')
    files = set()
    for pattern in ('solvers/IncNavierStokesSolver/EquationSystems/RANS*.cpp',
                    'solvers/IncNavierStokesSolver/EquationSystems/RANS*.h',
                    'library/SolverUtils/Turbulence/RANS*.*',
                    'library/SolverUtils/Turbulence/AnalyticalWallFunction.*',
                    'library/UnitTests/SolverUtils/RANS*Probe.cpp',
                    'library/UnitTests/SolverUtils/RANSWallFieldAudit.cpp'):
        files.update(source.glob(pattern))
    files.update(source/p for p in ('library/SpatialDomains/MeshGraphIOXml.cpp',
        'library/SolverUtils/CMakeLists.txt','solvers/IncNavierStokesSolver/CMakeLists.txt',
        'solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt',
        'solvers/IncNavierStokesSolver/Tests/RANS/generate_terrain.py',
        'solvers/IncNavierStokesSolver/Tests/RANS/terrain_initial.py',
        'solvers/IncNavierStokesSolver/Tests/RANS/terrain_continuity.py',
        'solvers/IncNavierStokesSolver/Tests/RANS/validate_terrain.py',
        'solvers/IncNavierStokesSolver/Tests/RANS/terrain_campaigns.py',
        'solvers/IncNavierStokesSolver/Tests/RANS/record_terrain.py'))
    sources = {str(p.relative_to(source)):sha(p) for p in sorted(files) if p.is_file()}
    binaries = {str(p.relative_to(build)):sha(p) for p in (
        build/'IncNavierStokesSolver',build/'RANSWallFieldAudit',
        build/'geometry-lib/libSpatialDomains.so')}
    for suite in suites.values():
        for name,digest in suite['binaries'].items():
            path = 'geometry-lib/'+name if name == 'libSpatialDomains.so' else name
            if binaries[path] != digest: raise RuntimeError('Binary changed after acceptance: '+name)
    def command(args):
        return subprocess.check_output(args,text=True).strip()
    geometry_resolution = [line.strip() for line in command(['ldd',str(build/'IncNavierStokesSolver')]).splitlines()
                           if 'libSpatialDomains.so' in line]
    if len(geometry_resolution) != 1 or str(build/'geometry-lib') not in geometry_resolution[0]:
        raise RuntimeError('Select the isolated source geometry library through LD_LIBRARY_PATH')
    result = dict(status='PASS',timestamp_utc=datetime.now(timezone.utc).isoformat(),
        source_revision=command(['git','-C',str(source),'rev-parse','HEAD']),
        working_tree=command(['git','-C',str(source),'status','--short']),
        source_sha256=sources,binary_sha256=binaries,
        environment=dict(platform=platform.platform(),python=sys.version,
            geometry_library_resolution=geometry_resolution,
            compiler=command(['c++','--version']).splitlines()[0],
            cmake=command(['cmake','--version']).splitlines()[0],
            mpi=command(['mpirun','--version']).splitlines()[0],
            threads={key:os.environ.get(key) for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}),
        suites=suites)
    if args.ctest: result['ctest'] = dict(path=str(args.ctest),sha256=sha(args.ctest),text=args.ctest.read_text())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
