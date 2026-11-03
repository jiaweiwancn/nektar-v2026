#!/usr/bin/env python3
"""Targeted physical-orientation qualification and unchanged M6 operator replay.

The incoming vector is reconstructed independently in the physical plane
frame. This tests a constructor admission rule; no nonlinear evolution is
repeated. Optional replays compare the complete native mathematical reports
against the prior passed M6 qualification, including all weak loads/actions.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import xml.etree.ElementTree as ET

sys.dont_write_bytecode = True
THREADS = {name: '1' for name in
           ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'BLIS_NUM_THREADS')}
os.environ.update(THREADS)
import numpy as np
from cls_epawf_m6_cases import transition_session
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m3_cases import set_info, set_parameter, write

EXACT = {'ids', 'sign', 'offsets', 'sizes', 'cell_ownership', 'schema', 'dimension',
         'ranks', 'owner', 'geometry', 'block', 'size', 'degree', 'order',
         'quadrature_order', 'quadrature_orders', 'total_points'}
ADDED_ORIENTATION_METADATA = {'constructor_only', 'incoming_flow_cartesian',
                              'actual_wall_coordinate_rate', 'reference_departure_maximum'}


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

    def near(self, actual, expected, name, tolerance=1e-10):
        actual, expected = np.asarray(actual, float), np.asarray(expected, float)
        self.require(actual.shape == expected.shape and np.all(np.isfinite(actual)) and
                     np.all(np.isfinite(expected)), 'Invalid numerical data: '+name)
        value = float(np.max(np.abs(actual-expected)/np.maximum(1., np.abs(expected)), initial=0.))
        self.numeric_values += int(actual.size)
        self.errors[name] = max(value, self.errors.get(name, 0.))
        self.require(value <= tolerance, f'{name}: {value} > {tolerance}')

    def compare(self, actual, expected, name='operator', root=True):
        """Retain every original JSON field; topology is exact, floats tight."""
        key = name.rsplit('.', 1)[-1]
        if isinstance(expected, dict):
            self.require(isinstance(actual, dict), 'Different JSON type: '+name)
            added, removed = set(actual)-set(expected), set(expected)-set(actual)
            self.require(not removed and (not added or root and added <= ADDED_ORIENTATION_METADATA),
                         'Different mathematical JSON keys: '+name+' added='+str(sorted(added))+' removed='+str(sorted(removed)))
            for field in expected:
                self.compare(actual[field], expected[field], name+'.'+field, False)
        elif isinstance(expected, list):
            self.require(isinstance(actual, list) and len(actual) == len(expected), 'Different JSON list shape: '+name)
            try:
                a, b = np.asarray(actual), np.asarray(expected)
                numeric = a.shape == b.shape and a.dtype.kind in 'biuf' and b.dtype.kind in 'biuf'
            except ValueError:
                numeric = False
            if numeric:
                if key in EXACT:
                    self.require(np.array_equal(a, b), 'Different exact native structure: '+name)
                else:
                    self.near(a, b, 'replay_'+name)
            else:
                for index, (a, b) in enumerate(zip(actual, expected)):
                    self.compare(a, b, f'{name}[{index}]', False)
        elif isinstance(expected, (bool, str)) or expected is None:
            self.require(type(actual) is type(expected) and actual == expected, 'Different exact metadata: '+name)
        elif key in EXACT:
            self.require(actual == expected, 'Different exact native integer: '+name)
        else:
            self.near(actual, expected, 'replay_'+name)


class Runner:
    def __init__(self, build, output, audit, receipt):
        self.build, self.output, self.audit, self.receipt = build, output, audit, receipt
        self.environment = dict(os.environ, **THREADS, PYTHONDONTWRITEBYTECODE='1')
        self.environment['LD_LIBRARY_PATH'] = str(build/'geometry-lib')+':'+self.environment.get('LD_LIBRARY_PATH', '')
        self.logs = output/'logs'
        self.logs.mkdir()
        receipt['environment'] = dict(THREADS, LD_LIBRARY_PATH=self.environment['LD_LIBRARY_PATH'],
                                      PYTHONDONTWRITEBYTECODE='1', TMPDIR=self.environment.get('TMPDIR'),
                                      cpu_affinity=sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None)

    def save(self):
        self.receipt.update(checks=self.audit.checks, maximum_scaled_errors=self.audit.errors,
                            numeric_values_compared=self.audit.numeric_values)
        (self.output/'validation.json').write_text(json.dumps(self.receipt, indent=2, allow_nan=False)+'\n')

    def command(self, name, argv, cwd, timeout=180):
        argv, log, start = list(map(str, argv)), self.logs/(name+'.log'), time.monotonic()
        code, failure = 127, None
        try:
            with log.open('w') as stream:
                result = subprocess.run(argv, cwd=cwd, env=self.environment, stdout=stream,
                                        stderr=subprocess.STDOUT, timeout=timeout, check=False)
                code = result.returncode
        except subprocess.TimeoutExpired:
            code, failure = 124, 'Command timed out'
        except OSError as error:
            failure = str(error)
        if failure:
            with log.open('a') as stream:
                stream.write('\n'+failure+'\n')
        self.receipt['commands'].append(dict(name=name, argv=argv, cwd=str(cwd), exit_code=code,
                                             ranks=1, threads=THREADS, seconds=time.monotonic()-start,
                                             log=str(log), log_sha256=digest(log)))
        self.save()
        print(f'{name}: exit {code}', flush=True)
        return code, log.read_text()

    def git(self, source, name, *arguments):
        code, text = self.command(name, ['git', '-c', 'gc.auto=0', '-C', source, *arguments], source, 60)
        self.audit.require(code == 0, 'Git metadata command failed: '+name)
        return text.strip()


def frame(sx, sy, dx, dy):
    normal = np.array([-sx, -sy, 1.])/np.sqrt(1+sx*sx+sy*sy)
    tangent1 = np.array([1., 0., sx])/np.sqrt(1+sx*sx)
    tangent2 = np.cross(normal, tangent1)
    flow = dx*tangent1+dy*tangent2
    return flow, float(np.sqrt(1+sx*sx)*flow[0]), normal


def constructor_session(template, directory, plane, dx, dy):
    xml = transition_session(template, directory/'case.xml', family='flat', dimension=3,
                             degree=3, steps=0, probe=True, quadrature=8)
    root = ET.parse(xml).getroot()
    sx, sy = (1., 1.) if plane else (0., 0.)
    if plane:
        physical_ell = np.array([-4., -3., -2., -1., 1., 2.])
        tensor_mesh(root, [0., 2., 4.], 3, 3, 1., 1., xknots=(physical_ell/np.sqrt(1+sx*sx)).tolist())
        set_info(root, 'RANSWallMap', 'Plane')
    for name, value in {'RANSWallMapSlopeX': sx, 'RANSWallMapSlopeY': sy,
                        'RANSWallDirectionX': dx, 'RANSWallDirectionY': dy,
                        'RANSWallProbeConstructorOnly': 1}.items():
        set_parameter(root, name, repr(value))
    return write(root, xml), sx, sy


def constructor_checks(run, template):
    root = math.sqrt(.99)
    cases = (
        ('tilted-positive-parameter-negative-physical', True, .1, root, False),
        ('tilted-negative-parameter-positive-physical', True, -.1, -root, True),
        ('tilted-q1-tangent', True, 1., 0., True),
        ('tilted-positive-parameter-positive-physical', True, .1, -root, True),
        ('flat-positive', False, 1., 0., True),
        ('flat-negative', False, -1., 0., False),
    )
    audit, results = run.audit, run.receipt['constructors']
    for name, plane, dx, dy, accept in cases:
        directory = run.output/name
        directory.mkdir()
        xml, sx, sy = constructor_session(template, directory, plane, dx, dy)
        expected_flow, expected_rate, normal = frame(sx, sy, dx, dy)
        audit.near(expected_flow@normal, 0., 'independent_frame_tangent', 1e-14)
        audit.near(expected_flow@expected_flow, 1., 'independent_frame_unit', 1e-14)
        audit.require((expected_rate > 0) == accept, 'Independent orientation fixture has wrong sign')
        report = directory/'probe.json'
        result = dict(name=name, plane=plane, slopes=[sx, sy], parameter_direction=[dx, dy],
                      expected_accept=accept, expected_incoming_flow_cartesian=expected_flow.tolist(),
                      expected_actual_wall_coordinate_rate=expected_rate, xml=str(xml), xml_sha256=digest(xml))
        results.append(result)
        code, text = run.command(name, [run.build/'RANSCurvedWallProbe', xml, report], directory)
        result['exit_code'] = code
        if accept:
            audit.require(code == 0 and report.is_file(), 'Valid physical incoming orientation rejected: '+name)
            raw = json.loads(report.read_text())
            audit.require(raw.get('schema') == 1 and raw.get('constructor_only') is True and
                          raw.get('dimension') == 3 and raw.get('ranks') == 1, 'Wrong constructor-only probe output')
            audit.near(raw['incoming_flow_cartesian'], expected_flow, 'native_cartesian_incoming_'+name, 1e-12)
            audit.near(raw['actual_wall_coordinate_rate'], expected_rate, 'native_physical_orientation_'+name, 1e-12)
            departure = float(raw['reference_departure_maximum'])
            audit.require(np.isfinite(departure) and 0 <= departure <= 1e-24,
                          'Accepted orientation does not preserve the actual projected reference sensors')
            result.update(native_report=str(report), native_report_sha256=digest(report),
                          actual_incoming_flow_cartesian=raw['incoming_flow_cartesian'],
                          actual_wall_coordinate_rate=raw['actual_wall_coordinate_rate'],
                          reference_departure_maximum=departure)
        else:
            audit.require(code not in (0, 124, 127), 'Invalid physical orientation was accepted or did not execute: '+name)
            lower = text.lower()
            audit.require('coordinate' in lower and any(word in lower for word in ('incoming', 'direction', 'orientation', 'flow')),
                          'Negative orientation failed for a different initialization reason: '+name)
            audit.require(not report.exists(), 'Rejected constructor emitted a success report')
            result['rejection_log'] = str(run.logs/(name+'.log'))
        run.save()
    audit.require(sum(case['expected_accept'] for case in results) == 4 and
                  sum(not case['expected_accept'] for case in results) == 2, 'Missing accepted/rejected orientation coverage')


def operator_replay(run, template, qualification):
    audit = run.audit
    original_validation = qualification/'validation.json'
    audit.require(original_validation.is_file(), 'Missing original M6 qualification receipt')
    original = json.loads(original_validation.read_text())
    audit.require(original.get('status', original.get('verification_status')) == 'passed' and
                  original.get('success', True), 'Original M6 qualification did not pass')
    run.receipt['original_qualification'] = dict(path=str(original_validation), sha256=digest(original_validation),
                                               source_commit=original.get('source_commit', original.get('tested_source_commit')))
    for family, dimension in (('flat', 2), ('feet', 2), ('flat', 3)):
        name = f'operator-{family}-d{dimension}'
        original_report = qualification/name/'probe.json'
        audit.require(original_report.is_file(), 'Missing original native operator '+name)
        expected = json.loads(original_report.read_text())
        directory = run.output/name
        directory.mkdir()
        xml = transition_session(template, directory/'case.xml', family=family, dimension=dimension,
                                 degree=3, steps=0, probe=True, quadrature=8)
        report = directory/'probe.json'
        record = dict(name=name, original_report=str(original_report), original_report_sha256=digest(original_report),
                      original_xml_sha256=digest(qualification/name/'case.xml'),
                      replay_xml=str(xml), replay_xml_sha256=digest(xml))
        run.receipt['operator_replays'].append(record)
        code, _ = run.command(name, [run.build/'RANSCurvedWallProbe', xml, report], directory, 1800)
        audit.require(code == 0 and report.is_file(), 'Native operator replay failed: '+name)
        actual = json.loads(report.read_text())
        before = audit.numeric_values
        audit.compare(actual, expected, name)
        added = sorted(set(actual)-set(expected))
        if 'constructor_only' in added:
            audit.require(actual['constructor_only'] is False, 'Full replay ran only its constructor')
        if 'incoming_flow_cartesian' in added:
            audit.near(actual['incoming_flow_cartesian'], [1., 0., 0.], 'replay_added_flow_metadata', 1e-12)
        if 'actual_wall_coordinate_rate' in added:
            audit.near(actual['actual_wall_coordinate_rate'], 1., 'replay_added_orientation_metadata', 1e-12)
        if 'reference_departure_maximum' in added:
            value = float(actual['reference_departure_maximum'])
            audit.require(np.isfinite(value) and 0 <= value <= 1e-24, 'Invalid added reference metadata')
        record.update(replay_report=str(report), replay_report_sha256=digest(report),
                      compared_numeric_values=audit.numeric_values-before, added_orientation_metadata=added,
                      original_mathematical_fields_unchanged=True)
        run.save()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--qualification', type=Path)
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Native orientation/replay qualification belongs on Linux')
    output, build = args.output.resolve(), args.build.resolve()
    output.mkdir(parents=True, exist_ok=False)
    audit = Audit()
    receipt = dict(schema=1, milestone='M6', success=False, verification_status='failed',
                   started_utc=datetime.now(timezone.utc).isoformat(), commands=[], failures=[],
                   constructors=[], operator_replays=[], build=str(build), output=str(output))
    run = Runner(build, output, audit, receipt)
    source = Path(__file__).resolve().parents[4]
    template = Path(__file__).with_name('EquilibriumWallFlat.xml')
    try:
        receipt['source_commit'] = run.git(source, 'git-start-revision', 'rev-parse', 'HEAD')
        status = run.git(source, 'git-start-status', 'status', '--porcelain=v1', '--untracked-files=all')
        audit.require(not status, 'Source checkout is not clean at orientation start')
        receipt['source_clean_start'] = True
        program = build/'RANSCurvedWallProbe'
        audit.require(program.is_file(), 'Missing native curved-wall probe')
        receipt['program'] = dict(path=str(program), sha256=digest(program))
        constructor_checks(run, template)
        if args.qualification:
            operator_replay(run, template, args.qualification.resolve())
        receipt['success'] = True
    except Exception as error:
        receipt['failures'].append(str(error))
    try:
        revision = run.git(source, 'git-end-revision', 'rev-parse', 'HEAD')
        status = run.git(source, 'git-end-status', 'status', '--porcelain=v1', '--untracked-files=all')
        receipt.update(source_commit_end=revision, source_clean_end=not status)
        audit.require(not status and revision == receipt.get('source_commit'), 'Source changed during orientation/replay checks')
        if 'program' in receipt:
            end = digest(Path(receipt['program']['path']))
            receipt['program']['sha256_end'] = end
            audit.require(end == receipt['program']['sha256'], 'Native probe binary changed during orientation/replay checks')
    except Exception as error:
        receipt['success'] = False
        receipt['failures'].append(str(error))
    receipt.update(verification_status='passed' if receipt['success'] else 'failed',
                   finished_utc=datetime.now(timezone.utc).isoformat())
    run.save()
    print(json.dumps(dict(success=receipt['success'], checks=audit.checks, constructors=len(receipt['constructors']),
                          operator_replays=len(receipt['operator_replays']), failures=receipt['failures'])))
    return 0 if receipt['success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
