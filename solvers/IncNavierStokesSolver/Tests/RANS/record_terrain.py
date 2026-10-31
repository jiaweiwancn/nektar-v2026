#!/usr/bin/env python3
"""Collect compact Linux provenance and passed terrain acceptance receipts."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build',type=Path,required=True)
    parser.add_argument('--suite',type=Path,action='append',required=True,help='Passed evidence.json; may be repeated')
    parser.add_argument('--ctest',type=Path)
    parser.add_argument('--legacy',type=Path,help='Passed standard-AWF checks report, including independent kernel/CLS short regressions')
    parser.add_argument('--workflow',type=Path,help='Public example generation/preparation receipt')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Collect authoritative build/solver evidence on Linux')
    source = Path(__file__).resolve().parents[4]
    build = args.build.resolve()
    suites = {}
    completion = {}
    for path in args.suite:
        value = json.loads(path.read_text())
        if value['status'] != 'PASS' or not value['cases'] or not all(v['passed'] for v in value['cases'].values()):
            raise RuntimeError('An acceptance suite is incomplete: '+str(path))
        suites[value['suite']] = value
        if value['suite'] in ('pilot','pilot-2d','preservation'):
            for name,record in value['cases'].items():
                directory = path.parent/name
                root = ET.parse(directory/'neutral.xml').getroot()
                params = {p.text.split('=',1)[0].strip():float(p.text.split('=',1)[1])
                          for p in root.findall('CONDITIONS/PARAMETERS/P')
                          if p.text.split('=',1)[0].strip() in ('NumSteps','TimeStep')}
                expected = int(params['NumSteps'])
                with (directory/'neutral.rans.csv').open() as stream:
                    stages = list(csv.DictReader(stream))
                end = expected*params['TimeStep']
                if len(stages) != expected or not math.isclose(float(stages[-1]['time']),end,rel_tol=1e-11,abs_tol=1e-14):
                    raise RuntimeError('Incomplete physical time advancement: '+name)
                count,min_k,min_eps = 0,math.inf,math.inf
                with (directory/'field.csv').open() as stream:
                    for row in csv.DictReader(stream):
                        values = {key:float(value) for key,value in row.items()}
                        if not all(math.isfinite(value) for value in values.values()):
                            raise RuntimeError('Nonfinite saved physical field: '+name)
                        if not math.isclose(values['time'],end,rel_tol=1e-11,abs_tol=1e-14):
                            raise RuntimeError('Saved field has the wrong physical time: '+name)
                        count += 1
                        min_k = min(min_k,values['k'])
                        min_eps = min(min_eps,values['eps'] if 'eps' in values else values['eps_tilde'])
                if not count or min_k <= 0 or min_eps <= 0:
                    raise RuntimeError('Empty/nonpositive saved physical field: '+name)
                completion[value['suite']+'/'+name] = dict(steps=expected,dt=params['TimeStep'],final_time=end,
                    saved_physical_points=count,minimum_saved_k=min_k,minimum_saved_epsilon=min_eps,
                    stage_sha256=sha(directory/'neutral.rans.csv'),field_csv_sha256=sha(directory/'field.csv'))
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
        'solvers/IncNavierStokesSolver/Tests/RANS/record_terrain.py',
        'solvers/IncNavierStokesSolver/Tests/RANS/plot_terrain.py',
        'solvers/IncNavierStokesSolver/Examples/TerrainHill/generate.py'))
    binaries = {str(p.relative_to(build)):sha(p) for p in (
        build/'IncNavierStokesSolver',build/'RANSWallFieldAudit',
        build/'geometry-lib/libSpatialDomains.so')}
    legacy = None
    if args.legacy:
        legacy = {}
        for label,path,status in (
                ('standard_awf_checks',args.legacy,'passed'),
                ('independent_kernels',args.legacy.parent/'legacy-kernels/validation.json','passed'),
                ('cls_short_regression',args.legacy.parent/'legacy-cls/validation.json','pilot')):
            report = json.loads(path.read_text())
            if report['status'] != status:
                raise RuntimeError('Legacy regression did not complete: '+str(path))
            for name,digest in report.get('binaries',report.get('binary_sha256',{})).items():
                if sha(build/name) != digest: raise RuntimeError('Legacy binary changed: '+name)
                binaries[name] = digest
            for name,digest in report['source_sha256'].items():
                if sha(source/name) != digest: raise RuntimeError('Legacy source changed: '+name)
                files.add(source/name)
            legacy[label] = dict(path=str(path),report_sha256=sha(path),report=report)
    sources = {str(p.relative_to(source)):sha(p) for p in sorted(files) if p.is_file()}
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
        completed_physical_times=completion,
        environment=dict(platform=platform.platform(),python=sys.version,
            geometry_library_resolution=geometry_resolution,
            compiler=command(['c++','--version']).splitlines()[0],
            cmake=command(['cmake','--version']).splitlines()[0],
            mpi=command(['mpirun','--version']).splitlines()[0],
            threads={key:os.environ.get(key) for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS')}),
        suites=suites)
    if legacy: result['legacy_regression'] = legacy
    if args.workflow:
        workflow = json.loads(args.workflow.read_text())
        if workflow['status'] != 'PASS' or not workflow['cases']:
            raise RuntimeError('Public workflow did not pass')
        for case in workflow['cases']:
            initial = case['initial_audit']
            if sha(case['command'][2]) != case['fixture_sha256']:
                raise RuntimeError('Public prepared session changed')
            if sha(initial['command'][-1]) != initial['field_sha256']:
                raise RuntimeError('Public prepared field changed')
        result['public_workflow'] = dict(path=str(args.workflow),report_sha256=sha(args.workflow),report=workflow)
    if args.ctest: result['ctest'] = dict(path=str(args.ctest),sha256=sha(args.ctest),text=args.ctest.read_text())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
