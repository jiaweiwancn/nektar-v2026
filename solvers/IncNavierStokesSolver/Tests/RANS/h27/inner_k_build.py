#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Materialize the isolated inner-k diagnostic from immutable source19."""

import argparse
import json
from pathlib import Path
import shutil
import sys

from resources import atomic, identity


ROOT = Path('/home/jay/build/wj-h27-coupling-r1-20261117')
SOURCE19 = '56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf'
ELF19 = '37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b'
NATIVE_SUFFIXES = {'.cpp', '.hpp', '.h', '.cmake'}
NATIVE = Path('solvers/IncNavierStokesSolver/EquationSystems')
H27 = Path('solvers/IncNavierStokesSolver/Tests/RANS/h27')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--version',type=int,default=20)
    parser.add_argument('--diagnostic-main',type=Path)
    args = parser.parse_args()
    assert args.version>=20
    assert sys.platform.startswith('linux') and ROOT.is_dir()
    base, source = ROOT/'source-v19', ROOT/f'source-v{args.version}'
    base_index = ROOT/'source-v19.json'
    records = json.loads(base_index.read_text())
    assert records['expected_source_identity'] == SOURCE19
    assert identity(ROOT/'runtime-v19/IncNavierStokesSolver.elf')['sha256'] == ELF19
    expected = {x['path']: x['sha256'] for x in records['files']}
    count = 0
    for path in base.rglob('*'):
        if path.is_file() and (path.suffix in NATIVE_SUFFIXES or path.name == 'CMakeLists.txt'):
            assert identity(path)['sha256'] == expected[str(path)], path
            count += 1
    assert count == 1853, count
    package = ROOT/'inputs/inner-k-f2/overlay'
    allowed = {NATIVE/(f+'.'+suffix) for f in ('WJNativeReference','WJNativeCoupling') for suffix in ('cpp','h')}
    allowed.update(H27/name for name in ('inner_k_fixture.cpp','CMakeLists.txt','build.py','physical_oracle.py',
                                          'inner_k_mechanism.py','inner_k_columns.py'))
    assert {p.relative_to(package) for p in package.rglob('*') if p.is_file()} == allowed
    shutil.copytree(base, source, ignore=shutil.ignore_patterns('__pycache__'))
    overlay = []
    for rel in sorted(allowed):
        out = source/rel
        incoming = args.diagnostic_main if rel==H27/'inner_k_fixture.cpp' and args.diagnostic_main else package/rel
        shutil.copyfile(incoming, out)
        overlay.append(dict(relative_path=str(rel), **identity(out)))
    atomic(ROOT/f'source-v{args.version}.json', dict(schema='WJInnerKDiagnosticSourceV1', source19=identity(base_index),
        exact_source19_files=count, source19_identity=SOURCE19, verified_source19_elf=ELF19,
        overlay=overlay, files=[identity(p) for p in sorted(source.rglob('*')) if p.is_file()],
        scope='Isolated steady diagnostic source. The deployed normal application retains source19.'))
    print(f'INNER_K_SOURCE{args.version}_READY',count,flush=True)


if __name__=='__main__':
    main()
