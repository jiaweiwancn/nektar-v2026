#!/usr/bin/env python3
"""Collect completed Linux CV2 gates into a small, reproducible receipt."""
import argparse
from datetime import datetime, timezone
import json
import platform
from pathlib import Path
import subprocess
import sys

sys.dont_write_bytecode = True
from cao_tamura_experiment import ROOT, DATA, THREADS, frozen_data, sha, write_json


def command(args):
    return subprocess.check_output(args, text=True).strip()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--build',type=Path,required=True)
    for name in ('probes','controls','generator','workflow','utilities','preservation','legacy'):
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--ctest-log',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if not sys.platform.startswith('linux'):
        raise RuntimeError('Linux evidence is authoritative')
    frozen_data(); suites={}
    for name in ('probes','controls','generator','workflow','utilities','preservation','legacy'):
        directory=getattr(a,name).resolve()
        path=directory/('validation-checks.json' if name=='legacy' else 'evidence.json')
        report=json.loads(path.read_text())
        if report.get('status','').lower() not in ('pass','passed'):
            raise RuntimeError('Incomplete suite: '+str(path))
        suites[name]=dict(path=str(path),sha256=sha(path),report=report)
    ctest=a.ctest_log.read_text()
    if '100% tests passed, 0 tests failed out of 11' not in ctest:
        raise RuntimeError('All eleven CTests must pass')
    build=a.build.resolve()
    binaries={name:sha(build/name) for name in ('IncNavierStokesSolver','RANSWallFieldAudit',
              'RANSObservationAudit','RANSObservationKernelProbe')}
    for name in ('controls','workflow','probes'):
        report=suites[name]['report']
        used=('IncNavierStokesSolver','RANSWallFieldAudit') if name=='controls' else tuple(binaries)[:3]
        for binary in used:
            if report['binaries'][binary]!=binaries[binary]:
                raise RuntimeError('Changed tested executable: '+name+'/'+binary)
    if suites['generator']['report']['binaries']['RANSWallFieldAudit']!=binaries['RANSWallFieldAudit']:
        raise RuntimeError('Changed native preparer')
    geometry=sha(build/'geometry-lib/libSpatialDomains.so')
    if any(suites[n]['report']['geometry_sha256']!=geometry for n in ('controls','workflow','probes','generator')):
        raise RuntimeError('Changed tested geometry library')
    kernel=build/'rans-observation-kernels.json';kernels=json.loads(kernel.read_text())
    if kernels['status']!='PASS' or len(kernels['reference_cases'])!=32:
        raise RuntimeError('Missing reference/reconstruction kernel checks')
    if kernels.get('auxiliary_H_shear_flux_max_scaled_error',1.)>5e-10:
        raise RuntimeError('Missing/failed auxiliary interface shear check')
    legacy_children={}
    for name,status in (('legacy-kernels','passed'),('legacy-cls','pilot')):
        path=a.legacy/name/'validation.json';report=json.loads(path.read_text())
        if report['status']!=status:
            raise RuntimeError('Incomplete independent legacy regression: '+name)
        legacy_children[name]=dict(path=str(path.resolve()),sha256=sha(path),status=status)
    files=command(['git','-C',str(ROOT),'ls-files','--modified','--others','--exclude-standard']).splitlines()
    files=[s for s in files if not s.startswith('docs/literature/') and
           (ROOT/s).is_file() and Path(s).suffix in ('.h','.cpp','.py','.txt','.tex','.json','.md')]
    # Full receipts remain in the run root. Compact evidence retains all source
    # and report hashes, the actual cases, comparisons and unsupported coverage.
    compact={}
    for name,entry in suites.items():
        report=entry['report']
        selected={k:v for k,v in report.items() if k in ('results','cases','comparisons','spatial_refinement','status','suite')}
        if name=='legacy':
            selected={k:report[k] for k in ('status','suite','comparisons','oracle_maximum_scaled_errors')}
            selected['case_ids']=list(report['cases'])
        compact[name]=dict(path=entry['path'],sha256=entry['sha256'],**selected)
    result=dict(status='PASS',milestone='CV2',timestamp_utc=datetime.now(timezone.utc).isoformat(),
                source_parent=command(['git','-C',str(ROOT),'rev-parse','HEAD']),
                source_sha256={s:sha(ROOT/s) for s in files},
                binaries=binaries,geometry_library_sha256=geometry,
                environment=dict(platform=platform.platform(),python=sys.version,
                    compiler=command(['c++','--version']).splitlines()[0],
                    mpi=command(['mpirun','--version']).splitlines()[0],threads=THREADS,
                    project_worker_cap=16,compiler_workers=12,maximum_MPI_ranks=4,ctest_concurrency=2),
                cv1_snapshot_sha256=sha(DATA/'snapshot.json'),
                ctest=dict(log=str(a.ctest_log.resolve()),sha256=sha(a.ctest_log),passed=11,total=11),
                kernels=kernels,suites=compact,independent_legacy=legacy_children,
                validation_commands={n:['taskset','-c',suites[n]['report']['cpus'],'python3',
                    str(ROOT/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cao_tamura.py'),
                    '--build',str(build),'--output',str(getattr(a,n).resolve()),'--cpus',suites[n]['report']['cpus'],'--suite',n]
                    for n in ('controls','workflow','probes','generator','utilities')},
                collector_command=sys.argv,
                interpretation=dict(experimental_calibration_complete=False,stationarity_established=False,
                    displaced_reference='Canonical equations retained; finite-nu pure log and physical-distance displaced CLS are not exact equilibria.',
                    unsupported=['below-hm k/epsilon/full Reynolds tensor','log-wall auxiliary layer','exact C1-foot auxiliary layer',
                                 'body-forced or dealiased auxiliary forcing','unknown experimental delta75 physical heights',
                                 'unscaled dimensional amplitudes','scalar/ratio endpoints awaiting additional matched sampling'],
                    continuity='Integrated strip continuity retained, including nonzero transpiration; C1 maximum convergence is not claimed.'))
    supplemental=build.parent/'kernel-final.log'
    if supplemental.exists():
        result['ctest']['supplemental_kernel']=dict(path=str(supplemental),sha256=sha(supplemental),
            command=['ctest','--test-dir',str(build),'-R','^RANSObservationKernels$','--output-on-failure'])
    write_json(a.output,result)
    print('PASS: CV2 receipt',a.output)


if __name__=='__main__':
    main()
