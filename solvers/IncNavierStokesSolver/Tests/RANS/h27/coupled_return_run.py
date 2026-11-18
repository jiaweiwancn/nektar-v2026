#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Run the prospectively bounded source19 cache-return verification on Linux."""
import argparse
import hashlib
import json
import shutil
import struct
import subprocess
import sys
from pathlib import Path

SOURCE = '56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
ELF = '37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b'
TABLE = 'b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d'
CASES = ('warm-a-first', 'warm-b', 'warm-a-return', 'retry-a', 'cold-a', 'cold-b', 'failed-a')


def read(path):
    return json.loads(path.read_text())


def history_parts(path):
    """Keep the binary header and four indivisible state/action byte ranges."""
    content = path.read_bytes()
    payload, checksum = content[:-64], content[-64:]
    if hashlib.sha256(payload).hexdigest().encode('ascii') != checksum:
        raise ValueError(f'checkpoint checksum differs: {path}')
    offset = payload.index(b'\n', payload.index(b'\n')+1)+1
    offset += struct.calcsize('<Qdd')
    header = payload[:offset]
    states = []
    for _ in range(4):
        begin = offset
        for _ in range(3):
            length, = struct.unpack_from('<Q', payload, offset)
            offset += 8*(1+length)
        states.append(payload[begin:offset])
    if offset != len(payload):
        raise ValueError(f'checkpoint layout differs: {path}')
    return header, states


def fixtures(root, package, name):
    from resources import atomic, identity

    output = root/'science'/name
    project = root/('project-'+name)
    output.mkdir()
    shutil.copytree(package/'coupled_return', project/'coupled_return')
    for file in ('coupled_return_probe.cpp', 'admission.md'):
        shutil.copy2(package/file, project/file)
    owner = root/'science/cutoff-owner-transition-v19'
    first = owner/'prefix.final.bin'
    second = owner/'cached.recovery-18.bin'
    header, parts = history_parts(first)
    _, later = history_parts(second)
    if parts[2] == later[2] or len(parts[2]) != len(later[2]):
        raise ValueError('the chosen stale action has no comparable nonzero change')
    parts[2] = later[2]
    payload = header+b''.join(parts)
    stale = output/'stale-A.bin'
    with stale.open('xb') as stream:
        stream.write(payload+hashlib.sha256(payload).hexdigest().encode('ascii'))
    atomic(output/'inputs.json', dict(schema='WJH27CoupledReturnInputsV1',
                                     A=identity(first), B=identity(second),
                                     stale_A=identity(stale)))
    for suffix in ('', '.opt'):
        shutil.copy2(owner/('prefix.xml'+suffix), output/('case.xml'+suffix))
    atomic(project/'manifest.json', dict(
        schema='WJH27CoupledReturnProjectV1', source=identity(root/'source-v19.json'),
        selected_runtime=identity(root/'science/runtime-v19.json'),
        libraries=[identity(root/'build-v19'/path) for path in
                   ('libWJNativeMixedFlow.a', 'h26-kernel/libWJNativeFlow.a')],
        project=[identity(path) for path in sorted(project.rglob('*')) if path.is_file()],
        package_files=[identity(path) for path in sorted(package.rglob('*')) if path.is_file()],
        helper_overlay=identity(root/'tools-cutoff-owner-transition-v19/transition-overlay.json'),
        case=identity(output/'case.xml'), options=identity(output/'case.xml.opt'),
        inputs=identity(output/'inputs.json')))
    print('COUPLED_RETURN_FIXTURES_READY', flush=True)


def run(root, package, name):
    from resources import atomic, identity

    account = read(root/'accounting.json')
    if account['closed'] or any(v['status'] == 'RESERVED' for v in account['entries'].values()):
        raise ValueError('account is closed or an existing operation requires reconciliation')
    used = sum(v.get('seconds', 0.) for v in account['entries'].values())
    if used+1120+account['evidence_protected'] > account['ceiling']:
        raise ValueError('the whole admitted return matrix no longer fits')
    for path, expected in ((root/'science/reference-tables-v1.json', TABLE),
                           (root/'runtime-v19/IncNavierStokesSolver.elf', ELF)):
        if identity(path)['sha256'] != expected:
            raise ValueError(f'selected input differs: {path}')
    if read(root/'source-v19.json')['expected_source_identity'] != SOURCE:
        raise ValueError('wrong immutable numerical source identity')
    runner = [sys.executable, '-B', str(root/'tools/run.py'), '--root', str(root)]

    def operation(label, phase, seconds, *command):
        subprocess.run(runner+['--label', name+'-'+label, '--phase', phase,
                              '--seconds', str(seconds), '--', *map(str, command)], check=True)

    output = root/'science'/name
    project = root/('project-'+name)
    build = root/('build-'+name)
    operation('fixtures', 'analysis', 30, sys.executable, '-B', __file__,
              '--root', root, '--name', name, '--fixtures')
    operation('configure', 'build', 40, 'cmake', '-S', project/'coupled_return', '-B', build,
              '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_CXX_FLAGS=-ffp-contract=off -fno-fast-math',
              '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake',
              '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON',
              '-DH27_PINNED_SOURCE='+str(root/'source-v19'),
              '-DH27_PINNED_BUILD='+str(root/'build-v19'))
    operation('build', 'build', 180, 'cmake', '--build', build, '--parallel', '4')
    runtime = Path('/home/jay/build/wj-h22-early-phase-replay-20261115/runtime')
    loader = runtime/'glibc/ld-linux-x86-64.so.2'
    libraries = ':'.join(map(str, (runtime/'glibc', runtime/'geometry',
                                  '/home/jay/.local/opt/nektar-v2026-mpi/lib64')))
    probe = build/'WJH27CoupledReturnProbe.elf'
    shutil.copy2(build/'WJH27CoupledReturnProbe', probe)
    runtime_manifest = output/'probe-runtime.json'
    atomic(runtime_manifest, dict(probe=identity(probe), loader=identity(loader),
                                  library_path=libraries, project=identity(project/'manifest.json')))
    native = output/'native'
    operation('native', 'science', 180, loader, '--library-path', libraries,
              probe, output/'case.xml', output/'inputs.json', native)

    inputs = read(output/'inputs.json')
    table = root/'science/reference-tables-v1.json'
    checker = package/'coupled_return_check.py'
    helpers = root/'tools-cutoff-owner-transition-v19/h27'
    check = [sys.executable, '-B', '-c',
             'import runpy,sys;sys.path.insert(0,sys.argv.pop(1));'
             'sys.argv[0]=sys.argv.pop(1);runpy.run_path(sys.argv[0],run_name="__main__")',
             helpers, checker]
    manifest = dict(schema='WJH27CoupledReturnManifestV1', root=str(output), source_identity=SOURCE,
                    inputs={key: inputs[key] for key in ('A', 'B', 'stale_A')},
                    runtime=identity(runtime_manifest),
                    crossing_geometry=identity(root/'science/mms-requalification-v19/crossing/geometry.json'),
                    stale_rejection=identity(native/'stale-rejection.json'),
                    runs={name: identity(native/(name+'.result.json')) for name in CASES})
    original = output/'manifest-native.json'
    atomic(original, manifest)
    for case in CASES[:-1]:
        operation(case+'-current', 'analysis', 80, *check, '--manifest', original,
                  '--table', table, '--case', case, '--out', output/(case+'-current.json'),
                  '--seconds', '75')
    manifest['initial_manifest'] = identity(original)
    manifest['current_equations'] = {case: identity(output/(case+'-current.json')) for case in CASES[:-1]}
    final = output/'manifest.json'
    atomic(final, manifest)
    operation('independent', 'analysis', 120, *check, '--manifest', final,
              '--table', table, '--out', output/'independent.json', '--seconds', '115')
    print('SOURCE19_COUPLED_RETURN_QUALIFICATION_FINISHED', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--name', default='coupled-return-v19r1')
    parser.add_argument('--fixtures', action='store_true')
    args = parser.parse_args()
    if not sys.platform.startswith('linux') or args.root != Path('/home/jay/build/wj-h27-coupling-r1-20261117'):
        parser.error('requires authoritative Linux and the existing H27 root')
    if not args.name.startswith('coupled-return-v19') or not args.name.replace('-', '').isalnum():
        parser.error('requires a create-only coupled-return label')
    sys.path.insert(0, str(args.root/'tools'))
    package = Path(__file__).resolve().parent
    action = fixtures if args.fixtures else run
    action(args.root, package, args.name)


if __name__ == '__main__':
    main()
