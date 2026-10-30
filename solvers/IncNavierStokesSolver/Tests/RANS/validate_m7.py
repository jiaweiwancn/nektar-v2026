#!/usr/bin/env python3
"""Linux delivery acceptance of the documented neutral-ABL public workflows."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote

sys.dont_write_bytecode = True
THREADS = {key: '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(sys.platform.startswith('linux'), 'Run delivery acceptance on Linux')
    source = Path(__file__).resolve().parents[4]
    build, output = args.build.resolve(), args.output.resolve()
    require(not output.is_relative_to(source), 'Keep generated delivery data outside the source checkout')
    output.mkdir(parents=True, exist_ok=True)
    os.sched_setaffinity(0, sorted(os.sched_getaffinity(0))[-16:])
    solver, audit = build/'IncNavierStokesSolver', build/'RANSFieldAudit'
    require(solver.is_file() and audit.is_file(), 'Build the isolated solver and audit first')
    example = source/'solvers/IncNavierStokesSolver/Examples/NeutralABL'
    commands, cases, comparisons = [], {}, {}

    def execute(command, directory, label, timeout=1200, rejection=None):
        command = list(map(str, command))
        commands.append(dict(command=command, cwd=str(directory), timeout=timeout))
        print('Running '+label+' in '+str(directory), flush=True)
        result = subprocess.run(command, cwd=directory, env=os.environ.copy(), text=True,
                                capture_output=True, timeout=timeout)
        raw = result.stdout+result.stderr
        (directory/(label+'.log')).write_text(raw)
        require((result.returncode != 0 and rejection in raw) if rejection else result.returncode == 0,
                'Inspect '+str(directory/(label+'.log')))
        return raw

    def audit_value(command, directory, label):
        raw = execute(command, directory, label)
        return json.loads(raw[raw.index('{"time"'):])

    def solve(name, session, steps, ranks=1):
        directory = session.parent
        command = (['mpirun', '--bind-to', 'none', '--map-by', 'slot', '-np', str(ranks)] if ranks > 1 else [])
        execute(command+[solver, session.name], directory, 'solver')
        field = directory/(session.stem+'.fld')
        value = audit_value([audit, session.name, field.name, field.name, 'profiles.csv'], directory, 'audit')
        rows = [{key: float(v) for key, v in row.items()}
                for row in csv.DictReader((directory/(session.stem+'.rans.csv')).open())]
        require(len(rows) == steps, 'Incomplete steps: '+name)
        require(all(math.isfinite(v) for row in rows for v in row.values()), 'Nonfinite diagnostics: '+name)
        require(all(row['min_k'] > 0 and row['min_eps'] > 0 and row['guard_activations'] == 0 for row in rows),
                'Unexpected state guard: '+name)
        require(all(math.isfinite(v) for field_value in value['fields'].values() for v in field_value.values()),
                'Nonfinite field audit: '+name)
        require(value.get('history_directory') and Path(value['history_directory']).is_dir(), 'Missing history: '+name)
        artifacts = {p.name: digest(p) for p in directory.iterdir()
                     if p.is_file() and p.suffix in ('.xml', '.log', '.csv', '.fld', '.chk')}
        if field.is_dir():
            artifacts.update({str(p.relative_to(directory)): digest(p) for p in field.rglob('*') if p.is_file()})
        cases[name] = dict(session=str(session), field=str(field), ranks=ranks, steps=steps, audit=value,
            maximum_iterations=max(row['iterations'] for row in rows),
            maximum_weak_divergence=max(row['weak_divergence_l2'] for row in rows), artifacts_sha256=artifacts)
        return cases[name]

    def prepare(name, model='cls_switched', dimension=2, steps=4, dt=.5, ranks=1, perturbation=0):
        directory = output/name
        directory.mkdir(exist_ok=True)
        session = directory/'neutral.xml'
        execute([sys.executable, example/'generate.py', session, '--model', model, '--dimension', dimension,
                 '--steps', steps, '--dt', dt, '--checkpoint', 2, '--perturbation', perturbation,
                 '--linear-tolerance', '1e-15'], directory, 'prepare')
        return solve(name, session, steps, ranks)

    def resume(name, prefix, steps, ranks=4, field=None, mode='History', history=None):
        directory = output/name
        directory.mkdir(exist_ok=True)
        session = directory/'resume.xml'
        command = [sys.executable, example/'restart.py', prefix['session'], field or prefix['field'],
                   session, '--steps', steps, '--mode', mode]
        if history is not None:
            command += ['--history-path', history]
        execute(command, directory, 'prepare-restart')
        return solve(name, session, steps, ranks)

    def compare(name, first, second, tolerance=1e-7):
        value = audit_value([audit, first['session'], first['field'], second['field']], output, name)
        largest = max(v['difference_linf'] for v in value['fields'].values())
        require(largest < tolerance, 'Continuation differs from uninterrupted fields: '+name)
        comparisons[name] = dict(maximum_absolute_linf=largest, tolerance=tolerance,
                                 fields=value['fields'], time=value['time'])

    documents = [source/p for p in ('README.md',
        'solvers/IncNavierStokesSolver/Examples/README.md',
        'solvers/IncNavierStokesSolver/Examples/NeutralABL/README.md',
        'solvers/IncNavierStokesSolver/Tests/RANS/README.md',
        'docs/user-guide/solvers/incompressible-ns.tex',
        'docs/developer-guide/solvers/incompressible-solver.tex')]
    links = []
    for document in documents:
        if document.suffix != '.md':
            continue
        for target in re.findall(r'\[[^\]]*\]\(([^\s]+)\)', document.read_text()):
            if ':' in target or target.startswith('#'):
                continue
            path = (document.parent/unquote(target.split('#', 1)[0])).resolve()
            require(path.exists(), 'Broken documentation link: '+str(document)+' -> '+target)
            links.append(dict(document=str(document.relative_to(source)), target=target))

    baseline = json.loads((source/'docs/plans/neutral-abl-rans-awf-milestone-6-figures/evidence.json').read_text())
    numerical = {name: sha for name, sha in baseline['source_sha256'].items() if Path(name).suffix in ('.cpp', '.h')}
    require(all(digest(source/name) == sha for name, sha in numerical.items()), 'Numerical source changed since M6')
    ctest = execute(['ctest', '--test-dir', build, '--parallel', '1', '--output-on-failure'], output, 'ctest')
    require('100% tests passed' in ctest and 'out of 8' in ctest, 'Expected eight smoke/kernel tests')

    for model in ('standard', 'cls_original', 'cls_fixed', 'cls_switched'):
        prepare(model, model=model)
    saved_field = Path(cases['standard']['field'])
    saved_hash = digest(saved_field)
    execute([sys.executable, example/'generate.py', saved_field], output,
            'reject-generator-field-output', rejection='.xml output path')
    execute([sys.executable, example/'restart.py', cases['standard']['session'], saved_field,
             saved_field, '--steps', 1], output, 'reject-restart-field-output', rejection='.xml output path')
    require(digest(saved_field) == saved_hash, 'A preparation command modified the saved field')
    prepare('standard3d', model='standard', dimension=3, steps=2, dt=.05, ranks=4)
    prepare('cls3d', dimension=3, steps=2, dt=.05, ranks=4, perturbation=.01)
    for model in ('standard', 'cls_switched'):
        full = prepare('full_'+model, model=model, steps=8)
        continuation = resume('resume_'+model, cases[model], 4)
        compare('compare_'+model, full, continuation)
    checkpoint = Path(cases['cls_switched']['session']).parent/'neutral_1.chk'
    require(checkpoint.exists(), 'Missing documented first checkpoint')
    continuation = resume('checkpoint_cls', cases['cls_switched'], 6, field=checkpoint)
    compare('compare_checkpoint', cases['full_cls_switched'], continuation)
    moved = output/'relocated-history'
    require(not moved.exists(), 'Use a fresh output directory for the delivery check')
    shutil.copytree(cases['cls_switched']['audit']['history_directory'], moved)
    continuation = resume('relocated_cls', cases['cls_switched'], 4, history=moved)
    compare('compare_relocated', cases['full_cls_switched'], continuation)
    resume('cold_standard', cases['standard'], 1, ranks=1, mode='Cold')
    require('RANS cold restart' in (output/'cold_standard/solver.log').read_text(), 'Cold mode did not report startup')
    full = prepare('full_cls3d', dimension=3, steps=4, dt=.05, ranks=4, perturbation=.01)
    continuation = resume('resume_cls3d', cases['cls3d'], 2, ranks=2)
    compare('compare_3d_repartition', full, continuation, 1e-6)

    pollution = [str(p.relative_to(source)) for p in (source/'solvers/IncNavierStokesSolver/Tests/RANS').iterdir()
                 if p.suffix in ('.fld', '.chk', '.csv', '.his', '.opt') or '.rans-' in p.name]
    require(not pollution, 'Smoke tests wrote generated data into the source checkout: '+str(pollution))
    files = documents+[example/'generate.py', example/'restart.py', Path(__file__).resolve(),
        source/'solvers/IncNavierStokesSolver/CMakeLists.txt',
        source/'solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt']
    report = dict(schema=1, status='passed', milestone=7, utc=datetime.now(timezone.utc).isoformat(),
        tested_parent=subprocess.check_output(['git', '-C', source, 'rev-parse', 'HEAD'], text=True).strip(),
        source_sha256={str(p.relative_to(source)): digest(p) for p in files},
        unchanged_numerical_source_sha256=numerical,
        m6_evidence_sha256=digest(source/'docs/plans/neutral-abl-rans-awf-milestone-6-figures/evidence.json'),
        binary_sha256={p.name: digest(p) for p in (solver, audit)},
        platform=platform.platform(), python=platform.python_version(), threads=THREADS,
        maximum_project_workers=4, maximum_mpi_ranks=4, concurrent_cases=1,
        cpu_affinity=sorted(os.sched_getaffinity(0)), ctest_count=8, ctest_log_sha256=digest(output/'ctest.log'),
        documentation_links=links, source_pollution=pollution, cases=cases, comparisons=comparisons, commands=commands,
        protected_field_output=dict(rejected_commands=2,field_sha256=saved_hash),
        interpretation='Short delivery checks; five-flow accuracy and higher-degree limits remain the M6 findings')
    (output/'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status='passed', solver_cases=len(cases), ctests=8,
        maximum_restart_difference=max(v['maximum_absolute_linf'] for v in comparisons.values())), indent=2), flush=True)


if __name__ == '__main__':
    main()
