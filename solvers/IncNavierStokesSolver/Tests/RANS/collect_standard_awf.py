#!/usr/bin/env python3
"""Collect the Linux standard-AWF acceptance, source and retained-attempt evidence."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

sys.dont_write_bytecode = True


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(value, message):
    if not value:
        raise RuntimeError(message)


def command(argv):
    return subprocess.check_output(list(map(str, argv)), text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(sys.platform.startswith('linux'), 'Collect authoritative evidence on Linux')
    source = Path(__file__).resolve().parents[4]
    run, output = args.run.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    checks = json.loads((run/'checks-final/validation-checks.json').read_text())
    study = json.loads((run/'study/validation-study.json').read_text())
    preparation = json.loads((run/'public-preparation/validation.json').read_text())
    require(all(r['status'] == 'passed' for r in (checks, study, preparation)), 'Incomplete acceptance')
    require(len(checks['cases']) == 20 and len(study['cases']) == 12 and len(preparation['cases']) == 14,
            'Incomplete case matrix')
    require(checks['binaries'] == study['binaries'], 'Study and checks used different executables')
    require(all(digest(run/'build'/name) == sha for name,sha in checks['binaries'].items()),
            'Executable changed after acceptance')
    for report in (checks, study):
        require(all((source/name).is_file() and digest(source/name) == sha
                    for name, sha in report['source_sha256'].items()), 'Tested source differs from delivery')
    required = [f'five_d{dim}_{mode}' for dim in (2, 3) for mode in ('fixed', 'switched')]+['richards_norris']
    require(all(study['cases'][n]['engineering_target_met'] and study['cases'][n]['spec']['flow_times'] >= 5
                for n in required), 'Reference preservation gate failed')
    for filename, src in [('validation-checks.json', run/'checks-final/validation-checks.json'),
                          ('validation-study.json', run/'study/validation-study.json'),
                          ('preparation.json', run/'public-preparation/validation.json')]:
        if src.resolve() != (output/filename).resolve():
            shutil.copyfile(src, output/filename)
    for path in sorted((run/'figures').glob('*.png')):
        shutil.copyfile(path, output/path.name)
    source_paths = set()
    for directory in (source/'library/SolverUtils/Turbulence', source/'solvers/IncNavierStokesSolver'):
        source_paths.update(p for p in directory.rglob('*') if p.is_file() and p.suffix in ('.cpp', '.h'))
    inputs = ['README.md', 'docs/plans/standard-kepsilon-equilibrium-awf.md',
              'docs/plans/standard-kepsilon-awf-evidence/README.md', 'docs/plans/neutral-abl-rans-awf-proposal.md',
              'docs/user-guide/solvers/incompressible-ns.tex', 'docs/developer-guide/solvers/incompressible-solver.tex',
              'library/UnitTests/SolverUtils/RANSKernelProbe.cpp', 'library/UnitTests/SolverUtils/RANSFieldAudit.cpp',
              'library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp', 'solvers/IncNavierStokesSolver/CMakeLists.txt',
              'solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt',
              'solvers/IncNavierStokesSolver/Tests/RANS/README.md',
              'solvers/IncNavierStokesSolver/Tests/RANS/StandardAWFEquilibrium.xml',
              'solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md',
              'solvers/IncNavierStokesSolver/Examples/NeutralABL/generate.py',
              'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py']
    source_paths.update(source/p for p in inputs)
    source_paths.update((source/'solvers/IncNavierStokesSolver/Tests/RANS').glob('*.py'))
    # Retain unsuccessful numerical attempts and the controls that explained
    # comparison assumptions. These paths are optional in a fresh reproduction.
    retained = {}
    directories = [(run/'checks/original', 'small_step_linear_1e13_failure'),
                   (run/'checks-final/switched_off-small-step-failure', 'strong_departure_small_step_failure'),
                   (run/'checks-final/small_P_fallback_failure', 'small_P_fallback_failure'),
                   (run/'checks-final/partial_switch_inlet', 'locally_mixed_inlet_control'),
                   (run/'checks-final/partial_switch_original_inlet', 'locally_mixed_original_control'),
                   (run/'checks-final/periodic_budget_default_tolerance', 'budget_at_profile_tolerance')]
    for parent in ('iteration-pilots', 'precision-pilots'):
        if (run/parent).exists():
            directories.extend((p, parent+'/'+p.name) for p in sorted((run/parent).iterdir()) if p.is_dir())
    for directory, label in directories:
        if not directory.exists():
            continue
        log = directory/'solver.log'
        stages = list(directory.glob('*.rans.csv'))
        accepted = len(list(csv.DictReader(stages[0].open()))) if stages else 0
        retained[label] = dict(path=str(directory), completed_steps=accepted,
            iteration_limit=log.exists() and 'coupled implicit stage failed to converge' in log.read_text(),
            sha256={str(p.relative_to(directory)): digest(p) for p in sorted(directory.iterdir())
                    if p.is_file() and p.suffix in ('.xml', '.csv', '.log', '.json')})
    budgets = {}
    for name in ('periodic_budget', 'periodic_budget_default_tolerance'):
        path = run/'checks-final'/name/(('periodic_budget' if name.endswith('tolerance') else name)+'.awf.csv')
        if path.exists():
            data = list(csv.DictReader(path.open()))
            budgets[name] = dict(maximum_absolute_momentum=max(abs(float(r['momentum_budget_x'])) for r in data),
                                 maximum_absolute_k=max(abs(float(r['k_budget'])) for r in data),
                                 available=all(float(r['budget_available']) == 1 for r in data), sha256=digest(path))
    manuals = {}
    for name in ('user-guide', 'developer-guide'):
        pdf = run/'manuals'/name/(name+'.pdf')
        require(pdf.is_file(), 'Compile '+name)
        info = command(['pdfinfo', pdf])
        pages = next(int(line.split(':')[1]) for line in info.splitlines() if line.startswith('Pages:'))
        manuals[name] = dict(pdf=str(pdf), pages=pages, sha256=digest(pdf),
            source_sha256=digest(source/('docs/'+name+'/solvers/'+('incompressible-ns.tex' if name == 'user-guide' else 'incompressible-solver.tex'))))
    linked = {}
    for p in sorted((run/'manuals/qa').glob('*.png')):
        linked[str(p.relative_to(run))] = digest(p)
    # These are observations, not mutations of the primary Linux checkout.
    primary = Path('/media/jay/新加卷/nektar-v2026')
    ctest = run/'checks-final/ctest.log'
    require('100% tests passed' in ctest.read_text() and 'out of 9' in ctest.read_text(), 'Missing CTest pass')
    report = dict(status='passed', utc=datetime.now(timezone.utc).isoformat(),
        meaning='Five-flow degree-four preservation, 2D convergence and explicit workflow gates passed; retained limits remain separate',
        source_parent=command(['git', '-C', source, 'rev-parse', 'HEAD']),
        source_sha256={str(p.relative_to(source)): digest(p) for p in sorted(source_paths)},
        binary_sha256=checks['binaries'], platform=platform.platform(),
        python=platform.python_version(), packages={p:version(p) for p in ('numpy','scipy','sympy','matplotlib')},
        compiler=command(['g++', '--version']).splitlines()[0], cmake=command(['cmake', '--version']).splitlines()[0],
        mpi=command(['mpirun', '--version']).splitlines()[0],
        allocation=dict(project_cpu_limit=16, build_jobs=8, study_jobs=2, maximum_study_ranks=8,
                        maximum_check_ranks=4, numerical_threads_per_rank=1, affinity=list(range(112,128))),
        ctest=dict(tests=9, log_sha256=digest(ctest)),
        inputs_sha256={p.name:digest(p) for p in sorted(output.iterdir()) if p.name != 'delivery.json' and p.is_file()},
        retained_attempts=retained, conservation=budgets, manuals=manuals, rendered_manual_pages_sha256=linked,
        primary_linux=dict(head=command(['git','-C',primary,'rev-parse','HEAD']),
                           tracked_status=command(['git','-C',primary,'status','--porcelain','--untracked-files=no'])),
        remaining_limits=['Six-layer pointwise shear-stress target failure',
                          'Recorded small-step and stiff-P coupled-iteration failures',
                          'No long-time standard-AWF 3D polynomial convergence claim',
                          'Existing CLS degree-five 3D growing mode remains unresolved',
                          'Flat neutral HighRe MatrixFree benchmark; no terrain/thermal validation'])
    (output/'delivery.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='passed', source_files=len(source_paths), numerical_sources=sum(p.suffix in ('.cpp','.h') for p in source_paths),
                         retained_attempts=len(retained), conservation=budgets), indent=2))


if __name__ == '__main__':
    main()
