#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reject a default-shell dependency seal when CTest selected another library."""
import argparse
import hashlib
import os
from pathlib import Path
import subprocess


def selected(probe, env):
    lines = subprocess.check_output(['ldd', str(probe)], env=env, text=True).splitlines()
    line = next(s for s in lines if 'libSpatialDomains.so' in s)
    assert 'not found' not in line
    return Path(line.split('=>', 1)[1].split()[0]).resolve()


def main():
    p = argparse.ArgumentParser(); p.add_argument('--probe', type=Path, required=True)
    p.add_argument('--library', type=Path, required=True); a = p.parse_args()
    actual = selected(a.probe, os.environ.copy())
    assert actual == a.library.resolve(), 'native loader ignored source-matched SpatialDomains'
    print('Actual loaded library:', actual)
    print('SHA256:', hashlib.sha256(actual.read_bytes()).hexdigest())
    # This is a loader inspection, with no governing/initialization operation.
    assert a.library.resolve().parent in [Path(s).resolve() for s in os.environ['LD_LIBRARY_PATH'].split(':') if s]
    print('PASS: selected runtime directory and ELF dependency agree; zero governing solves')


if __name__ == '__main__': main()
