#!/usr/bin/env python3
"""Supervise one declared R3 benchmark or exact target; no auxiliary sequence.

Execution belongs to the authoritative Linux host and an explicit allocation.
No default time budget is supplied. Target execution requires the actual full-J
benchmark receipt and a frozen budget note, as well as current action gates.
"""
import argparse
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


def pin(path):
    path = Path(path).resolve()
    return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))


def load(path):
    return json.loads(Path(path).read_text())


def emit(message):
    if emit.disabled:
        return
    try:
        print(message, flush=True)
    except OSError as error:
        if not isinstance(error, BrokenPipeError) and error.errno != errno.EPIPE:
            raise
        emit.disabled = True
        try:
            descriptor = os.open(os.devnull, os.O_WRONLY)
            try: os.dup2(descriptor, sys.stdout.fileno())
            finally: os.close(descriptor)
        except (OSError, ValueError, AttributeError):
            sys.stdout = None
emit.disabled = False


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('epoch', 'session', 'manifest', 'output'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--action-gate', type=Path)
    p.add_argument('--phase', choices=('benchmark', 'target', 'r4-local', 'r5-local'), required=True)
    p.add_argument('--expected-source', required=True)
    p.add_argument('--seconds', type=int, required=True)
    p.add_argument('--cpu', type=int, default=0); p.add_argument('--metadata-cpu', type=int, default=50)
    p.add_argument('--cost-receipt', type=Path); p.add_argument('--budget-note', type=Path)
    a = p.parse_args()
    require(platform.system() == 'Linux' and a.seconds > 0, 'Authoritative Linux and a positive explicit budget are required')
    epoch, source, build = a.epoch.resolve(), a.epoch.resolve()/'source', a.epoch.resolve()/'build'
    session, output = a.session.resolve(), a.output.resolve()
    require(not output.exists() and not session.is_relative_to(source) and not output.is_relative_to(source),
        'Use fresh external outputs; never write into the producing checkout')
    output.mkdir(parents=True)
    lock = (epoch/'validation.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_SH)
    def git(*words):
        return subprocess.check_output(['git', '-C', str(source), *words], text=True).strip()
    require(git('rev-parse', 'HEAD') == a.expected_source and not git('status', '--porcelain'), 'Wrong/unclean producing checkout')
    root = ET.parse(session).getroot()
    info = {row.get('PROPERTY', '').upper():row.get('VALUE') for row in root.findall('.//SOLVERINFO/I')}
    r4 = info.get('RANSWALLGLOBALIZATION') == 'CoupledR4'
    r5 = info.get('RANSWALLGLOBALIZATION') == 'JointRestorationR5'
    expected_stage = {'benchmark': 'JacobianBenchmarkR3', 'target': 'ExactTargetR3',
                      'r4-local': 'DirectionExperimentR4', 'r5-local': 'DirectionExperimentR5'}[a.phase]
    require(info.get('RANSWALLSTAGETEST') == expected_stage and
        info.get('RANSWALLNONLINEARPOLICY') == 'ScalarFeasibleR3' and info.get('RANSWALLACTIONPOLICY') == 'VerifiedR3Action',
        'Changed explicit R3 phase/policy')
    require(a.phase != 'r4-local' or r4, 'Local R4 requires the opt-in direction')
    require(a.phase != 'r5-local' or r5, 'Local R5 requires explicit restoration')
    progress = Path(info.get('RANSWALLR3PROGRESSFILE', 'case.r3-progress.jsonl'))
    progress = progress if progress.is_absolute() else session.parent/progress
    report = session.with_suffix('.evolution.json')
    require(not progress.exists() and not report.exists(), 'Native outputs already exist; no silent retry/reuse')
    gate = load(a.action_gate) if a.action_gate else {}
    built = load(epoch/'runtime-fingerprint.json')
    executable = pin(build/'IncNavierStokesSolver')
    require(a.phase in ('r4-local', 'r5-local') or (gate.get('schema') == 1 and gate.get('source_commit') == a.expected_source and
        gate.get('executable_sha256') == executable['sha256'] and gate.get('absolute_action_tolerance') == 1e-7 and
        gate.get('required_state_action_gates_met') is True and gate.get('active_dependencies_covered') is True),
        'Current independent absolute action/active-dependency gates are not complete')
    manifest = load(a.manifest)
    require(manifest['current_target_xml'] == pin(session) and manifest['current_executable'] == executable,
        'Manifest/current native input or executable changed')
    if a.phase == 'target':
        require(a.cost_receipt is not None and a.budget_note is not None, 'Target budget requires the measured complete-J receipt and frozen note')
        cost = load(a.cost_receipt)
        require(cost.get('phase') == ('r5-local' if r5 else 'r4-local' if r4 else 'benchmark') and cost.get('status') == 'completed' and
            cost.get('source_commit') == a.expected_source and cost.get('executable') == executable and
            cost.get('runtime_unchanged') is True and (r4 or r5 or cost.get('durable_collection', {}).get('complete_J_benchmark') is True),
            'Target cost basis is not the completed current-source full-J benchmark')
        if r4:
            require(a.seconds == 600 and gate.get('local_evidence') == pin(a.cost_receipt) and
                gate.get('useful_feasible_native_descent') is True and
                gate.get('base_state') == manifest['state'],
                'R4 pilot requires the measured useful descent at the exact saved guess and a 600s cap')
        if r5:
            require(a.seconds <= 900 and gate.get('local_evidence') == pin(a.cost_receipt) and
                gate.get('useful_joint_native_descent') is True and gate.get('base_state') == manifest['state'] and
                0 < int(info['RANSWALLR5UPDATES']) == gate.get('maximum_updates') <= 32,
                'R5 continuation requires joint descent at its exact saved base and a bounded measured budget')
        require(bool(a.budget_note.read_text().strip()), 'Empty target budget rationale')
    environment = dict(os.environ, **{key:'1' for key in THREADS}, PYTHONDONTWRITEBYTECODE='1',
        TMPDIR=str(epoch/'tmp'), LD_LIBRARY_PATH=str(build/'geometry-lib'))
    receipt = dict(schema=1, phase=a.phase, status='running', source_commit=a.expected_source,
        source_clean_start=True, started_utc=datetime.now(timezone.utc).isoformat(), runner=pin(__file__),
        session=pin(session), manifest=pin(a.manifest), action_gate=pin(a.action_gate) if a.action_gate else None, executable=executable,
        numerical_threads={key:'1' for key in THREADS}, ranks=1, native_cpu=a.cpu, metadata_cpu=a.metadata_cpu,
        maximum_seconds=a.seconds, cost_receipt=pin(a.cost_receipt) if a.cost_receipt else None,
        budget_note=pin(a.budget_note) if a.budget_note else None, commands=[], failures=[])
    def save():
        path = output/'receipt.json'; temporary = output/'receipt.json.tmp'
        with temporary.open('w') as stream:
            json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
        temporary.replace(path)
    def metadata(name, argv):
        log = output/(name+'.log'); start = time.monotonic()
        command = dict(name=name, argv=['taskset', '-c', str(a.metadata_cpu), *map(str, argv)])
        receipt['commands'].append(command); save()
        with log.open('w') as stream:
            completed = subprocess.run(command['argv'], cwd=output, env=environment, stdout=stream,
                stderr=subprocess.STDOUT, timeout=240)
        command.update(exit_code=completed.returncode, seconds=time.monotonic()-start, log=pin(log)); save()
        require(completed.returncode == 0, 'Metadata command failed: '+name)
    fingerprint = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
    process = None; before_complete = False; native_started = None
    def stop(child):
        try: os.killpg(child.pid, signal.SIGTERM)
        except ProcessLookupError: pass
        try: child.wait(timeout=10)
        except subprocess.TimeoutExpired:
            try: os.killpg(child.pid, signal.SIGKILL)
            except ProcessLookupError: pass
            child.wait()
    save()
    try:
        metadata('fingerprint-before', ['python3', '-B', fingerprint, '--root', epoch,
            '--output', output/'runtime-before.json', '--compare', epoch/'runtime-fingerprint.json'])
        before_complete = True
        current = load(output/'runtime-before.json')
        require(current['source_commit'] == a.expected_source and current['runtime_artifacts'] == built['runtime_artifacts'], 'Producing runtime changed')
        argv = ['taskset', '-c', str(a.cpu), executable['path'], str(session)]
        native = dict(name='native', argv=argv, cwd=str(session.parent), timeout_seconds=a.seconds)
        receipt['commands'].append(native); save(); native_started = time.monotonic()
        log_path = output/'native.log'
        with log_path.open('w') as stream:
            process = subprocess.Popen(argv, cwd=session.parent, env=environment, stdout=stream,
                stderr=subprocess.STDOUT, start_new_session=True)
            native['pid'] = process.pid; save(); last = native_started
            while process.poll() is None:
                now = time.monotonic()
                if now-native_started >= a.seconds:
                    stop(process)
                    native['exit_code'] = 124; break
                if now-last >= 30:
                    emit('R3 '+a.phase+' alive '+str(round(now-native_started, 1))+'s'); last = now
                time.sleep(.2)
            else:
                native['exit_code'] = process.returncode
        native.update(seconds=time.monotonic()-native_started, log=pin(log_path)); save()
        receipt['status'] = 'native_returned' if native['exit_code'] == 0 else 'native_stopped'
    except Exception as error:
        receipt['status'] = 'failed'; receipt['failures'].append(dict(type=type(error).__name__, message=str(error)))
    finally:
        # Even a non-stdout supervisor exception must not leave a native child
        # running while a supposed post-completion fingerprint is captured.
        if process is not None and process.poll() is None:
            stop(process)
            native.update(exit_code=process.returncode, seconds=time.monotonic()-native_started,
                interrupted_by_supervisor_failure=True, log=pin(log_path)); save()
        if before_complete:
            try:
                metadata('fingerprint-after', ['python3', '-B', fingerprint, '--root', epoch,
                    '--output', output/'runtime-after.json', '--compare', output/'runtime-before.json'])
                receipt['runtime_unchanged'] = True
            except Exception as error:
                receipt['runtime_unchanged'] = False; receipt['failures'].append(dict(type=type(error).__name__, message=str(error)))
        if progress.exists():
            try:
                spec = importlib.util.spec_from_file_location('r3_durable_collection', Path(__file__).with_name('collect-progress.py'))
                collector = importlib.util.module_from_spec(spec); spec.loader.exec_module(collector)
                collection = collector.collect(progress)
                collection_path = output/'durable-collection.json'
                with collection_path.open('x') as stream:
                    json.dump(collection, stream, indent=2, allow_nan=False); stream.write('\n')
                receipt['durable_collection'] = dict(file=pin(collection_path),
                    complete_J_benchmark=collection['complete_J_benchmark'], terminal_record_available=collection['terminal_record_available'],
                    completed_iterate_count=collection['completed_iterate_count'], native_converged_observation=collection['native_converged_observation'])
            except Exception as error:
                receipt['failures'].append(dict(type=type(error).__name__, message=str(error), phase='durable_collection'))
        if report.exists(): receipt['native_report'] = pin(report)
        receipt.update(source_commit_end=git('rev-parse', 'HEAD'), source_clean_end=not bool(git('status', '--porcelain')),
            finished_utc=datetime.now(timezone.utc).isoformat())
        if receipt['status'] == 'native_returned' and receipt.get('runtime_unchanged') and not receipt['failures']:
            receipt['status'] = 'completed'
        save(); emit(json.dumps(dict(status=receipt['status'], durable=receipt.get('durable_collection'))))
    return 0 if receipt['status'] == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
