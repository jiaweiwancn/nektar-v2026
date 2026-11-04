#!/usr/bin/env python3
"""Run one immutable R3 fixed-state probe on allocated Linux CPU 0.

No nonlinear stage, integrator history or trajectory is advanced. XML/packet
copies, images, complete runtime fingerprints and all produced artifacts are
bound before/after. A stopped probe remains a stopped observation, not a pass.
"""
import argparse
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time

THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
    'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
RUNTIME_KEYS = ('source_inputs', 'runtime_artifacts', 'runtime_bindings',
    'runtime_environment', 'python_packages', 'tool_images')
PROGRAMS = {'action':'RANSNonpressureActionProbe', 'resolution':'RANSWallResolutionProbe'}


def require(value, message):
    if not value:
        raise ValueError(message)


def pin(path):
    path = Path(path).resolve(); value = hashlib.sha256()
    with path.open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            value.update(part)
    return dict(path=str(path), bytes=path.stat().st_size, sha256=value.hexdigest())


def load(path):
    return json.loads(Path(path).read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


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
            fd = os.open(os.devnull, os.O_WRONLY)
            try: os.dup2(fd, sys.stdout.fileno())
            finally: os.close(fd)
        except (OSError, ValueError, AttributeError):
            sys.stdout = None
emit.disabled = False


def stop(process):
    try: os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError: pass
    try: process.wait(timeout=10)
    except subprocess.TimeoutExpired:
        try: os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError: pass
        process.wait()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--epoch', '--source-epoch', dest='epoch', type=Path, required=True)
    p.add_argument('--xml', type=Path, required=True)
    p.add_argument('--packet', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--kind', choices=tuple(PROGRAMS), required=True)
    p.add_argument('--seconds', type=int, required=True)
    p.add_argument('--expected-source', required=True)
    p.add_argument('--direction', type=Path, help='Optional retained action witness binary; action kind only')
    p.add_argument('--metadata-cpu', type=int, default=50)
    a = p.parse_args()
    require(platform.system() == 'Linux' and a.seconds > 0 and a.metadata_cpu >= 0,
        'Authoritative Linux, an explicit positive timeout and one allocated metadata CPU are required')
    require(a.kind == 'action' or a.direction is None, 'Resolution probe does not accept a direction argument')
    epoch, output = a.epoch.resolve(), a.output.resolve()
    source, build = epoch/'source', epoch/'build'
    require(not output.exists() and output.is_relative_to(epoch.parent) and
        not output.is_relative_to(source) and not output.is_relative_to(build),
        'Use a fresh external per-case output below the run root')
    for key in ('OMPI_COMM_WORLD_SIZE', 'PMI_SIZE', 'PMIX_SIZE'):
        require(key not in os.environ or int(os.environ[key]) == 1, 'Fixed-state probe must not inherit multiple MPI ranks')
    lock = (epoch/'validation.lock').open('a'); fcntl.flock(lock, fcntl.LOCK_SH)
    git = lambda *words: subprocess.check_output(['git', '-C', str(source), *words], text=True).strip()
    require(git('rev-parse', 'HEAD') == a.expected_source and not git('status', '--porcelain'),
        'Producing checkout is not the requested clean source')
    baseline_path = epoch/'runtime-fingerprint.json'; baseline = load(baseline_path)
    require(baseline.get('success') is True and baseline.get('clean') is True and
        baseline.get('source_commit') == a.expected_source, 'No successful matching producing runtime fingerprint')
    image_key = 'build/'+PROGRAMS[a.kind]
    require(image_key in baseline['runtime_artifacts'] and image_key in baseline['runtime_resolved_paths'],
        'Probe image is missing from the frozen producing runtime')
    executable = pin(baseline['runtime_resolved_paths'][image_key])
    require(executable['sha256'] == baseline['runtime_artifacts'][image_key] and
        Path(executable['path']) == (build/PROGRAMS[a.kind]).resolve(), 'Producing probe image/path changed')
    output.mkdir(parents=True); case = output/'case'; case.mkdir()
    xml, packet = case/'case.xml', case/'input.json'
    original_xml, original_packet = pin(a.xml), pin(a.packet)
    shutil.copy2(a.xml, xml); shutil.copy2(a.packet, packet)
    require(pin(xml)['sha256'] == original_xml['sha256'] and pin(packet)['sha256'] == original_packet['sha256'],
        'Copied XML/packet bytes differ')
    report = case/'report.json'; matrix = case/'report.json.J.bin'
    inputs = dict(original_xml=original_xml, original_packet=original_packet, xml=pin(xml), packet=pin(packet),
        executable=executable, baseline_runtime=pin(baseline_path))
    if a.direction is not None:
        inputs['direction'] = pin(a.direction)
    environment = dict(os.environ, **{key:'1' for key in THREADS}, PYTHONDONTWRITEBYTECODE='1',
        TMPDIR=str(epoch/'tmp'), LD_LIBRARY_PATH=str(build/'geometry-lib'))
    receipt = dict(schema=1, status='preparing', scope='One fixed-state probe; no solver/history/trajectory advancement',
        kind=a.kind, source_commit=a.expected_source, source_clean_start=True,
        started_utc=datetime.now(timezone.utc).isoformat(), runner=pin(__file__), inputs_before=inputs,
        ranks=1, numerical_threads={key:'1' for key in THREADS}, native_cpu=0, metadata_cpu=a.metadata_cpu,
        maximum_concurrent_workers=1, timeout_seconds=a.seconds, child_environment={key:environment[key] for key in
            (*THREADS, 'LD_LIBRARY_PATH', 'TMPDIR', 'PYTHONDONTWRITEBYTECODE')}, commands=[], failures=[],
        artifacts_before=dict(report=dict(path=str(report), exists=False), matrix=dict(path=str(matrix), exists=False)),
        history_advanced=False, trajectory_calls=0, nonlinear_stage_calls=0)
    def save():
        path = output/'receipt.json'; temporary = output/'receipt.json.tmp'
        with temporary.open('w') as stream:
            json.dump(receipt, stream, indent=2, allow_nan=False); stream.write('\n'); stream.flush(); os.fsync(stream.fileno())
        temporary.replace(path)
    def fingerprint(name, comparison):
        result = output/(name+'.json'); log = output/(name+'.log')
        tool = source/'docs/plans/cls-epawf-revised-evidence/m7/fingerprint-runtime.py'
        argv = ['taskset', '-c', str(a.metadata_cpu), 'python3', '-B', str(tool), '--root', str(epoch),
            '--output', str(result), '--compare', str(comparison)]
        row = dict(name=name, argv=argv, cwd=str(output), timeout_seconds=240)
        receipt['commands'].append(row); save(); began = time.monotonic()
        try:
            with log.open('w') as stream:
                completed = subprocess.run(argv, cwd=output, env=environment, stdout=stream,
                    stderr=subprocess.STDOUT, timeout=240)
            row.update(exit_code=completed.returncode, seconds=time.monotonic()-began, log=pin(log)); save()
            require(completed.returncode == 0, 'Runtime fingerprint failed: '+name)
        except subprocess.TimeoutExpired:
            row.update(exit_code=124, seconds=time.monotonic()-began, log=pin(log)); save(); raise
        current = load(result)
        require(current.get('success') is True and current.get('clean') is True and
            current.get('source_commit') == a.expected_source, 'Fingerprint source/runtime scope changed')
        return current
    process = None; before = None; native = None; native_started = None; log = output/'native.log'
    save()
    try:
        before = fingerprint('runtime-before', baseline_path)
        receipt['runtime_before'] = pin(output/'runtime-before.json')
        receipt['runtime_sections_before'] = {key:canonical(before[key]) for key in RUNTIME_KEYS}
        argv = ['taskset', '-c', '0', executable['path'], str(xml), str(packet), str(report)]
        if a.direction is not None: argv.append(inputs['direction']['path'])
        native = dict(name='fixed-state-'+a.kind, argv=argv, cwd=str(case), timeout_seconds=a.seconds)
        receipt['commands'].append(native); receipt['status'] = 'running'; save()
        native_started = time.monotonic()
        with log.open('w') as stream:
            process = subprocess.Popen(argv, cwd=case, env=environment, stdout=stream,
                stderr=subprocess.STDOUT, start_new_session=True)
            native['pid'] = process.pid; save(); last = native_started
            while process.poll() is None:
                now = time.monotonic()
                if now-native_started >= a.seconds:
                    native['timeout_requested'] = True; stop(process); native['exit_code'] = 124; break
                if now-last >= 30:
                    emit('R3 fixed-state '+a.kind+' alive '+str(round(now-native_started, 1))+'s'); last = now
                time.sleep(.2)
            else: native['exit_code'] = process.returncode
        native.update(process_returncode=process.returncode, seconds=time.monotonic()-native_started, log=pin(log)); save()
        receipt['status'] = 'native_returned' if native['exit_code'] == 0 else 'native_stopped'
    except Exception as error:
        receipt['status'] = 'failed'; receipt['failures'].append(dict(type=type(error).__name__, message=str(error)))
    finally:
        if process is not None and process.poll() is None:
            stop(process)
            native.update(exit_code=process.returncode, process_returncode=process.returncode,
                seconds=time.monotonic()-native_started, interrupted_by_supervisor_failure=True, log=pin(log)); save()
        if before is not None:
            try:
                after = fingerprint('runtime-after', output/'runtime-before.json')
                receipt['runtime_after'] = pin(output/'runtime-after.json')
                receipt['runtime_sections_after'] = {key:canonical(after[key]) for key in RUNTIME_KEYS}
                receipt['runtime_unchanged'] = all(before[key] == after[key] for key in RUNTIME_KEYS)
            except Exception as error:
                receipt['runtime_unchanged'] = False
                receipt['failures'].append(dict(type=type(error).__name__, message=str(error), phase='runtime-after'))
        try:
            receipt['inputs_after'] = {name:pin(item['path']) for name,item in inputs.items()}
            require(receipt['inputs_after'] == inputs, 'Immutable XML/packet/direction/program/runtime input changed')
            receipt['input_hashes_unchanged'] = True
        except Exception as error:
            receipt['input_hashes_unchanged'] = False
            receipt['failures'].append(dict(type=type(error).__name__, message=str(error), phase='input_pins'))
        artifacts = {}
        for name,path in (('report',report), ('matrix',matrix)):
            artifacts[name] = pin(path) if path.is_file() else dict(path=str(path), available=False)
        receipt['artifacts_after'] = artifacts
        if report.is_file():
            receipt['report'] = artifacts['report']
            try:
                raw = load(report)
                require(raw.get('schema') == 1 and raw.get('dimension') == 2 and raw.get('mapped') is False,
                    'Unexpected fixed-state report schema/scope')
                if a.kind == 'action':
                    require(matrix.is_file() and Path(raw['matrix_file']).resolve() == matrix and
                        raw['matrix_size'] == len(raw['state']) and matrix.stat().st_size == 8*len(raw['state'])**2 and
                        raw.get('ranks') == 1 and raw.get('action_policy') == 'VerifiedR3Action',
                        'Incomplete/changed fixed-state action report or matrix')
                    receipt['matrix'] = artifacts['matrix']
            except Exception as error:
                receipt['failures'].append(dict(type=type(error).__name__, message=str(error), phase='report_binding'))
        elif receipt['status'] == 'native_returned':
            receipt['failures'].append(dict(type='MissingReport', message='Native process returned without a complete report'))
        receipt.update(source_commit_end=git('rev-parse', 'HEAD'), source_clean_end=not bool(git('status', '--porcelain')),
            finished_utc=datetime.now(timezone.utc).isoformat())
        if (receipt['status'] == 'native_returned' and receipt.get('runtime_unchanged') is True and
            receipt.get('input_hashes_unchanged') is True and receipt['source_clean_end'] and
            receipt['source_commit_end'] == a.expected_source and not receipt['failures']):
            receipt['status'] = 'completed'
        receipt['success'] = receipt['status'] == 'completed'
        save(); emit(json.dumps(dict(status=receipt['status'], output=str(output), failures=receipt['failures'])))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
