#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Passive source/compiler/test inventory and explicit unlaunched-phase record."""
import argparse
import json
import platform
from pathlib import Path
import re
import subprocess
import sys
from wj_immutable import read,file_identity,file_hash,write_once,artifacts
from run_cls_epawf_picard import require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--compact',type=Path,required=True);a=p.parse_args()
    r=a.root.resolve();out=a.compact.resolve();sources={}
    for folder in [*sorted(r.glob('source-v*')),*sorted(r.glob('publication-source-v*'))]:
        git=lambda *args:subprocess.check_output(['git','-c','core.filemode=false','-C',str(folder),*args],text=True).strip()
        sources[str(folder)]=dict(revision=git('rev-parse','HEAD'),tree=git('rev-parse','HEAD^{tree}'),status=git('status','--porcelain'),
            python_sources=artifacts((folder/'solvers/IncNavierStokesSolver/Tests/RANS').glob('*.py')))
    inventory_text=(r/'operations/final-test-inventory/command.log').read_text()
    inventory,_=json.JSONDecoder().raw_decode(inventory_text[inventory_text.index('{'):])
    require(len(inventory['tests'])==63,'actual registered inventory differs')
    testlog=r/'checks-final/Testing/Temporary/LastTest.log';text=testlog.read_text()
    match=re.search(r'Ran (\d+) tests',text);require(match and int(match[1])==11 and '\nOK\n' in text,'final registered fixtures did not execute all11')
    refusal=(r/'operations/expected-joint-admission-refusal/command.log').read_text()
    require('joint previous barrier/admission did not pass' in refusal,'refusal reason differs')
    readbacks={};attempts={}
    for policy in ('EP','Budget'):
        s=read(r/'study-v2'/policy/'summary.json');readbacks[policy]=[];attempts[policy]=0
        for name,e in s['execution'].items():
            d=read(e['report']['path']);require(e['actual_attempts']==len(d['pimple_attempts'])==0,'resource refusal launched native evolution')
            readbacks[policy].append(dict(case=name,report=e['report'],steps=d['steps'],time=d['time'],dt=d['dt'],
                initialization_lift=d['initialization_lift_calls_iterations'],implicit_callbacks=d['implicit_callback_count'],
                checkpoint=file_identity(Path(d['checkpoint'])/'State.bin')))
        require(len(readbacks[policy])==3,'readback coverage differs')
        require(not any((r/'study-v2'/policy/name).exists() for name in ('coarse-64','fine-96','fine-128','baseline-64')),
                'unexpected denied-cadence output exists')
    spec={name:file_identity(r/'source-v8/docs/literature'/name) for name in ('wj_earsm_neutral_abl_ep_awf_note.tex','wj_earsm_budget_cutoff_awf_appendix.tex')}
    require([v['sha256'] for v in spec.values()]==['1959e96ea0b61f9deb9ff9ad1f5b26488a932a8451c2de6d3b6f7b7c9dbd9124','0ab52e0d813b8fb344973f75ca6d8da97ff389ef9f93476500ca9d1df1fbc47c'],'specifications differ')
    value=dict(schema='WJ7SourceAndTestLineageV1',status='PASS',sources=sources,specifications=spec,
        compiler=subprocess.check_output(['c++','--version'],text=True),platform=platform.platform(),python=sys.version,
        installed_environment=file_identity('/home/jay/.local/opt/nektar-v2026-mpi/env.sh'),
        registered_tests=63,selected_focused_ctests=1,executed_focused_fixtures=11,reused_prior_selected_ctests=38,
        readbacks=readbacks,new_governing_attempts=attempts,joint_refusal_verified=True,
        unlaunched=dict(coarse=[64,96,128],fine=[128,192,256],baseline=[64,96,128],new_arithmetic_per_policy=9),
        expected_setup_duplicate_refusal=dict(scope='tool-observed, before child launch; no numerical work; conservative source/I/O allowance covers it',
            original_preserved_source=str(r/'source-v6'),original_revision=sources[str(r/'source-v6')]['revision'],requested_revision='f33231f'),
        logs=artifacts([testlog,r/'operations/final-test-inventory/command.log',r/'operations/expected-joint-admission-refusal/command.log']),code=file_identity(__file__))
    write_once(out/'source-test-lineage.json',value);write_once(out/'test-inventory.json',inventory)
    with (out/'final-ctest-LastTest.log').open('xb') as stream:stream.write(testlog.read_bytes())
    print('PASS:63 registered,1 selected,11 fixtures,6 zero-solve imports,0 continuation attempts')


if __name__=='__main__':main()
