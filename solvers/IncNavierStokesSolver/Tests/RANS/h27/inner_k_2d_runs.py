#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Admit and run bounded serial Linux tranches for the inner-k 2D fixture."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess
import sys

from resources import identity, size, workers


OPENING = 19772.376549492823
V2_ELF = '37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b'
ROOT = Path('/home/jay/build/wj-h27-coupling-r1-20261117')
SCIENCE = ROOT/'science'
OUT = SCIENCE/'inner-k-f3'
CP2 = SCIENCE/'cp2-refinement-v19'
SOURCE = ROOT/'source-v21/solvers/IncNavierStokesSolver/Tests/RANS/h27'
INPUT = ROOT/'inputs/inner-k-f3'
APP = ROOT/'runtime-v21/WJH27InnerKFixture'
BOUNDARY = SCIENCE/'reference-tables-v1.json'
LOADER = ('import runpy,sys;sys.path[:0]=[sys.argv.pop(1),sys.argv.pop(1)];'
          'sys.argv[0]=sys.argv.pop(1);runpy.run_path(sys.argv[0],run_name="__main__")')
PY = [sys.executable,'-B','-c',LOADER,INPUT,SOURCE]


def read(path):
    return json.loads(path.read_text())


def operation(label, phase, cap, *argv):
    args = [sys.executable,'-B',ROOT/'tools/run.py','--root',ROOT,'--label','inner-k-f3-'+label,
            '--phase',phase,'--seconds',cap,'--',*argv]
    subprocess.run(list(map(str,args)),check=True)


def admission(caps, memory_guess, gate):
    a = read(ROOT/'accounting.json')
    charges = sum(e['seconds'] for e in a['entries'].values())
    package = sum(e['seconds'] for k,e in a['entries'].items() if k.startswith('inner-k-'))
    assert (ROOT/'evidence/inner-k-f2/settled.json').is_file() and read(gate)['F2']=='PASS'
    runtime = read(SCIENCE/'runtime-v21.json')
    assert identity(APP)==runtime['executables'][0]['wrapper']
    assert identity(Path(str(APP)+'.elf'))==runtime['executables'][0]['executable']
    assert identity(ROOT/'runtime-v19/IncNavierStokesSolver.elf')['sha256']==V2_ELF
    assert not a['closed'] and not any(e['status']=='RESERVED' for e in a['entries'].values())
    assert abs(charges-OPENING-package)<1e-6
    assert package+sum(caps.values())<=6000 and charges+sum(caps.values())<=a['ceiling']-a['evidence_protected']
    native_bytes = sum(size(p) for p in ROOT.iterdir() if p.is_dir() and p.name.startswith(('source','build','runtime')))
    science_bytes = size(SCIENCE)
    free = shutil.disk_usage(ROOT).free
    host_workers = workers()['workers']
    assert memory_guess<8*1024**3 and native_bytes<=8*1024**3 and science_bytes+.9*1024**3<=8*1024**3
    assert free>=15.9*1024**3 and host_workers+12<=72
    return dict(aggregate_seconds=charges, fixture_seconds=package, caps=caps, maximum_forecast_science_bytes=memory_guess,
                source_build_runtime_bytes=native_bytes, science_bytes=science_bytes, free_bytes=free,
                observed_host_workers=host_workers, f2_gate=identity(gate), runtime=identity(SCIENCE/'runtime-v21.json'))


def prepare_base(gate, forecast):
    from inner_k_fixture_sessions import SOURCE19, axes, edit_session
    from inner_k_fixture_oracle import atomic_out
    from physical_oracle import TABLE_SHA

    OUT.mkdir()
    base = CP2/'base.xml'
    old = CP2/'base.result.json'
    old_audit = CP2/'base-current.json'
    evidence = read(old_audit)
    assert evidence['native']['sha256']==identity(old)['sha256'] and evidence['core_steady_equations_status']=='MET'
    x, z = axes(base)
    assert len(x)==29 and len(z)==101
    session = edit_session(base,OUT/'base.xml',
        dict(WJInitialFields=old,WJInitialSourceIdentity=SOURCE19,WJFixtureInitialReferenceIdentity=TABLE_SHA,
             WJOutputPrefix=OUT/'base'),
        dict(WJSteadySeconds=810))
    atomic_out(OUT/'base-inputs.json',dict(schema='WJInnerKMatchedBoundaryInputsV1',session=session,
        table=identity(BOUNDARY), v2_session=identity(base),v2_state=identity(old),v2_audit=identity(old_audit),
        cells=[len(x)-1,len(z)-1],x=x,physical_z=z, admission=forecast,
        held_fixed=['physical geometry','essential lifts','v2 traction and scalar exterior loads','steady and quadrature settings'],
        same_boundary_reference=TABLE_SHA,f2_gate=identity(gate)))


def audit(label, cap, native, out, gate, order, *extra):
    operation(label,'analysis',cap,*PY,INPUT/'inner_k_2d.py','--table',BOUNDARY,'--gate',gate,
              '--seconds',cap-10,'--order',order,'--native',native,'--out',out,*extra)


def base(args):
    caps = dict(inputs=35,solve=860,current=235,matched=180,publication=60)
    forecast = admission(caps,4500000000,args.gate)
    operation('base-inputs','analysis',caps['inputs'],*PY,__file__,'prepare-base',args.gate,json.dumps(forecast))
    operation('base-steady','science',caps['solve'],APP,'steady-fixture',OUT/'base.xml',OUT/'base')
    native = OUT/'base/steady.json'
    audit('base-current',caps['current'],native,OUT/'base-current.json',args.gate,32,'--localize')
    audit('base-matched-v2',caps['matched'],native,OUT/'base-matched-v2.json',args.gate,20,
          '--other',CP2/'base.result.json','--control','baseline')
    print('INNER_K_F3_MATCHED_CASE_FINISHED',flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=('base','prepare-base'))
    parser.add_argument('gate',type=Path)
    parser.add_argument('forecast',nargs='?')
    args = parser.parse_args()
    if args.stage=='prepare-base':
        prepare_base(args.gate,json.loads(args.forecast))
    else:
        base(args)


if __name__=='__main__':
    main()
