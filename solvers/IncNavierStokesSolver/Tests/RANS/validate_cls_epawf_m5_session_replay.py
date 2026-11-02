#!/usr/bin/env python3
"""Targeted cached-session replay after the native probe filename fix.

Reuses the qualified XML and mathematical results. No nonlinear solve or full
M4/M5 qualification is repeated. Optimizer sidecars are recorded and retained.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
import validate_cls_epawf_m5 as physics
import validate_cls_epawf_m5_geometry as metrics
from cls_epawf_m5_boundary_oracle import boundary_projection

TEST = 'RANSEquilibriumEvolutionNative'
FAMILIES = ('plane', 'circle', 'smooth', 'variable', 'surface3d', 'feet')


def digest(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            result.update(block)
    return result.hexdigest()


class Audit:
    def __init__(self):
        self.checks, self.errors, self.numeric_values = 0, {}, 0

    def require(self, good, message):
        self.checks += 1
        if not good:
            raise RuntimeError(message)

    def near(self, actual, expected, name, scale=1., tolerance=1e-10):
        a, b = np.asarray(actual, float), np.asarray(expected, float)
        self.require(a.shape == b.shape and np.all(np.isfinite(a)) and np.all(np.isfinite(b)), 'invalid '+name)
        error = float(np.max(np.abs(a-b)/np.maximum(scale, np.abs(b)), initial=0.))
        self.numeric_values += int(a.size)
        self.errors[name] = max(error, self.errors.get(name, 0.))
        self.require(error <= tolerance, f'{name}: {error} > {tolerance}')

    def compare(self, actual, expected, name):
        """All JSON data, with exact integer topology and finite float checks."""
        structural = name.rsplit('.', 1)[-1] in {
            'ids', 'sign', 'offsets', 'sizes', 'cell_ownership', 'schema', 'dimension',
            'ranks', 'owner', 'geometry', 'block', 'size', 'degree', 'order',
            'quadrature_order', 'quadrature_orders', 'total_points'}
        if isinstance(expected, dict):
            self.require(isinstance(actual, dict) and actual.keys() == expected.keys(), 'different keys: '+name)
            for key in expected:
                self.compare(actual[key], expected[key], name+'.'+key)
        elif isinstance(expected, list):
            self.require(isinstance(actual, list) and len(actual) == len(expected), 'different list shape: '+name)
            try:
                a, b = np.asarray(actual), np.asarray(expected)
                numeric = a.dtype.kind in 'biuf' and b.dtype.kind in 'biuf' and a.shape == b.shape
            except ValueError:
                numeric = False
            if numeric:
                if structural:
                    self.require(np.array_equal(a, b), 'different native integer structure: '+name)
                else:
                    self.near(a, b, name)
            else:
                for index, (a, b) in enumerate(zip(actual, expected)):
                    self.compare(a, b, f'{name}[{index}]')
        elif isinstance(expected, (bool, str)) or expected is None:
            self.require(type(actual) is type(expected) and actual == expected, 'different exact metadata: '+name)
        elif structural:
            self.require(actual == expected, 'different native integer metadata: '+name)
        else:
            self.near(actual, expected, name)


def optimizer_files(directory):
    directory = Path(directory)
    files = set(directory.rglob('*.opt')) if directory.exists() else set()
    if (directory/'.opt').is_file():
        files.add(directory/'.opt')
    return [dict(path=str(path.relative_to(directory)), bytes=path.stat().st_size,
                 mtime_ns=path.stat().st_mtime_ns, sha256=digest(path))
            for path in sorted(files) if path.is_file()]


class Replay:
    def __init__(self, build, output, audit, report):
        self.build, self.output, self.audit, self.report = build, output, audit, report
        self.environment = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
        self.environment['LD_LIBRARY_PATH'] = str(build/'geometry-lib')+':'+self.environment.get('LD_LIBRARY_PATH', '')
        self.logs = output/'logs'
        self.logs.mkdir()
        report['environment'] = dict(THREADS, LD_LIBRARY_PATH=self.environment['LD_LIBRARY_PATH'],
                                     TMPDIR=self.environment.get('TMPDIR'), PYTHONDONTWRITEBYTECODE='1')

    def save(self):
        self.report.update(checks=self.audit.checks, maximum_scaled_errors=self.audit.errors,
                           numeric_values_compared=self.audit.numeric_values)
        (self.output/'validation.json').write_text(json.dumps(self.report, indent=2, allow_nan=False)+'\n')

    def run(self, name, command, cwd, timeout=1800):
        command = list(map(str, command)); log = self.logs/(name+'.log'); start = time.monotonic()
        code, failure = 127, None
        try:
            with log.open('w') as stream:
                result = subprocess.run(command, cwd=cwd, env=self.environment, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=timeout, check=False)
                code = result.returncode
        except subprocess.TimeoutExpired:
            code, failure = 124, 'command timed out'
        except OSError as error:
            failure = str(error)
        if failure:
            with log.open('a') as stream:
                stream.write('\n'+failure+'\n')
        record = dict(name=name, argv=command, cwd=str(cwd), exit_code=code, ranks=1, threads=THREADS,
                      seconds=time.monotonic()-start, log=str(log), log_sha256=digest(log))
        self.report['commands'].append(record); self.save()
        print(f'{name}: exit {code} ({record["seconds"]:.2f}s)', flush=True)
        return code, log.read_text()

    def failure(self, section, error):
        self.report['failures'].append(dict(section=section, error=str(error)))
        self.save()


def git(replay, source, name, *arguments):
    code, text = replay.run(name, ['git', '-C', source, *arguments], source, timeout=60)
    replay.audit.require(code == 0, 'git metadata command failed: '+name)
    return text.strip()


def junit(audit, filename):
    audit.require(filename.is_file(), 'missing CTest JUnit report')
    cases = ET.parse(filename).getroot().findall('.//testcase')
    audit.require(len(cases) == 1 and cases[0].get('name') == TEST, 'CTest did not execute exactly the requested test')
    audit.require(not any(cases[0].find(tag) is not None for tag in ('failure', 'error', 'skipped')), 'CTest JUnit has a failed/skipped test')
    return dict(path=str(filename), sha256=digest(filename), tests=1)


def ctest_replay(replay):
    audit, build = replay.audit, replay.build
    code, text = replay.run('ctest-discovery', ['ctest', '--test-dir', build, '--show-only=json-v1', '-R', '^'+TEST+'$'], build)
    audit.require(code == 0, 'CTest discovery failed')
    tests = json.loads(text)['tests']
    audit.require(len(tests) == 1 and tests[0]['name'] == TEST, 'requested CTest was not found uniquely')
    properties = {item['name']: item['value'] for item in tests[0].get('properties', [])}
    directory = Path(properties.get('WORKING_DIRECTORY', build/'rans-wall-evolution-native')).resolve()
    audit.require(directory.is_dir(), 'existing cached CTest working directory is missing')
    records = []
    for repeat in (1, 2):
        receipt = dict(repeat=repeat, cwd=str(directory), optimizer_before=optimizer_files(directory))
        filename = replay.output/f'legacy-repeat-{repeat}.xml'
        code, _ = replay.run(f'ctest-repeat-{repeat}', ['ctest', '--test-dir', build, '--parallel', '1',
                             '--output-on-failure', '-R', '^'+TEST+'$', '--output-junit', filename], build)
        receipt.update(exit_code=code, optimizer_after=optimizer_files(directory))
        try:
            audit.require(code == 0, f'cached CTest repeat {repeat} failed')
            receipt['junit'] = junit(audit, filename)
            before = {item['path'] for item in receipt['optimizer_before']}
            after = {item['path'] for item in receipt['optimizer_after']}
            audit.require(before <= after, 'cached CTest removed an existing optimizer sidecar')
        except Exception as error:
            replay.failure(f'ctest-repeat-{repeat}', error)
            receipt['error'] = str(error)
        records.append(receipt)
    replay.report['ctest_repeats'] = records


def m4_spatial(replay, source):
    directory = replay.output/'m4-spatial'
    code, _ = replay.run('m4-spatial', [sys.executable, source/'solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m4.py',
                         '--section', 'spatial', '--build', replay.build, '--output', directory], replay.output, timeout=1800)
    replay.audit.require(code == 0, 'targeted M4 spatial qualification failed')
    path = directory/'validation.json'; raw = json.loads(path.read_text())
    replay.audit.require(raw['verification_status'] == 'passed' and raw['section'] == 'spatial', 'M4 spatial receipt did not pass')
    replay.report['m4_spatial'] = dict(path=str(path), sha256=digest(path), checks=raw['checks'],
                                     maximum_scaled_errors=raw['maximum_scaled_errors'], commands=raw['commands'])


def operator_replay(replay, qualification, family):
    audit = replay.audit; original = qualification/('operator-'+family); directory = replay.output/('operator-'+family)
    directory.mkdir(); xml = directory/'case.xml'; shutil.copy2(original/'case.xml', xml)
    original_path = original/'probe.json'; baseline = json.loads(original_path.read_text())
    record = dict(family=family, original_report=str(original_path), original_report_sha256=digest(original_path),
                  original_xml_sha256=digest(original/'case.xml'), copied_xml_sha256=digest(xml), runs=[])
    audit.require(record['original_xml_sha256'] == record['copied_xml_sha256'], 'qualified XML changed during copy')
    reports = []
    for repeat in (1, 2):
        path = directory/f'report{repeat}.json'
        item = dict(repeat=repeat, optimizer_before=optimizer_files(directory))
        code, _ = replay.run(f'operator-{family}-{repeat}', [replay.build/'RANSCurvedWallProbe', xml, path], directory)
        item.update(exit_code=code, optimizer_after=optimizer_files(directory))
        if path.is_file():
            item.update(report=str(path), report_sha256=digest(path))
        record['runs'].append(item)
        if code == 0 and path.is_file():
            reports.append(json.loads(path.read_text()))
    replay.report['operators'].append(record); replay.save()
    audit.require(len(reports) == 2 and all(item['exit_code'] == 0 for item in record['runs']), 'both repeated operator commands must pass: '+family)
    for item in record['runs']:
        before = {cache['path'] for cache in item['optimizer_before']}
        after = {cache['path'] for cache in item['optimizer_after']}
        audit.require(before <= after, 'operator replay removed an existing optimizer cache: '+family)
    record['optimizer_cache_exercised'] = any(Path(item['path']).name == 'case.opt' for item in record['runs'][0]['optimizer_after'])
    audit.require(record['optimizer_cache_exercised'], 'first operator run did not leave case.opt for cached replay: '+family)
    audit.compare(reports[1], reports[0], family+'_repeat2_vs_repeat1')
    audit.compare(reports[0], baseline, family+'_repeat1_vs_qualified')
    audit.compare(reports[1], baseline, family+'_repeat2_vs_qualified')
    final = reports[1]
    start = time.monotonic(); initial_checks = physics.checks
    record['oracles'] = dict(residual=physics.residual_oracle(final), interface=physics.interface_oracle(final),
                             boundary=boundary_projection(final, physics.spaces(final), metrics.geometry, physics.near, physics.require),
                             weak_geometry=metrics.validate_wall(final))
    record.update(oracle_seconds=time.monotonic()-start, physics_checks=physics.checks-initial_checks,
                  physics_maximum_scaled_errors=dict(physics.errors), mathematical_replay_equal=True)
    replay.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--qualification', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); source = Path(__file__).resolve().parents[4]
    build, qualification, output = args.build.resolve(), args.qualification.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    audit = Audit(); report = dict(schema=1, purpose='M5 targeted cached-session replay', started_utc=datetime.now(timezone.utc).isoformat(),
                                  build=str(build), qualification=str(qualification), commands=[], failures=[], operators=[], success=False, verification_status='running',
                                  scope='probe XML filename fix only; cached CTest twice, M4 spatial only, six curved operator replays; no nonlinear solve')
    replay = Replay(build, output, audit, report); revision = None
    try:
        audit.require(platform.system() == 'Linux', 'native replay requires the authoritative Linux environment')
        revision = git(replay, source, 'source-start-commit', 'rev-parse', 'HEAD')
        status = git(replay, source, 'source-start-status', 'status', '--porcelain')
        report.update(source_commit=revision, source_start_status=status)
        audit.require(not status, 'targeted replay requires a clean source worktree')
        original = qualification/'validation.json'; previous = json.loads(original.read_text())
        audit.require(previous['verification_status'] == 'passed', 'original main qualification did not pass')
        report['original_qualification'] = dict(path=str(original), sha256=digest(original), source_commit=previous.get('source_commit'))
        report['binary_sha256'] = {name: digest(build/name) for name in ('IncNavierStokesSolver', 'RANSEquilibriumEvolutionProbe', 'RANSCurvedWallProbe', 'RANSEquilibriumTimeProbe')}
        for section, action in [('cached-ctest', lambda: ctest_replay(replay)), ('m4-spatial', lambda: m4_spatial(replay, source))]:
            try:
                action()
            except Exception as error:
                replay.failure(section, error)
        for family in FAMILIES:
            try:
                operator_replay(replay, qualification, family)
            except Exception as error:
                replay.failure('operator-'+family, error)
    except Exception as error:
        replay.failure('setup', error)
    finally:
        try:
            report['source_end_commit'] = git(replay, source, 'source-end-commit', 'rev-parse', 'HEAD')
            report['source_end_status'] = git(replay, source, 'source-end-status', 'status', '--porcelain')
            audit.require(revision == report['source_end_commit'] and not report['source_end_status'], 'source changed or was dirty during replay')
            if 'binary_sha256' in report:
                current = {name: digest(build/name) for name in report['binary_sha256']}
                audit.require(current == report['binary_sha256'], 'a tested executable changed during replay')
                report['binary_end_sha256'] = current
        except Exception as error:
            replay.failure('source-or-binary-end-state', error)
        report.update(finished_utc=datetime.now(timezone.utc).isoformat(), success=not report['failures'], verification_status='failed' if report['failures'] else 'passed',
                      physics_checks=physics.checks, physics_maximum_scaled_errors=physics.errors)
        replay.save()
    print(json.dumps(dict(status=report['verification_status'], checks=audit.checks, failures=report['failures'])), flush=True)
    return 1 if report['failures'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
