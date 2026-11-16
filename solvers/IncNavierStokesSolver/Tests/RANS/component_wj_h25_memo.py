#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact current-domain scalar memo, including ordered adverse-case inventories."""
import argparse,os,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from cls_epawf_m3_cases import set_info,write


def run(root):
    old=Path('/home/jay/build/wj-h16-repaired-monitor-20261114/science/small/A')
    pins=[dict(path=str(old/'case.xml'),bytes=60217,sha256='decfaa0d7028a423e810247899b01c023e54660520c0d0ed4c71edc7d908e252'),
        dict(path=str(old/'view/request.json'),bytes=1070,sha256='077cde9cc818eb442d389e279edad49fb69e891f349b3ac3826a971a1d37c43d')]
    for pin in pins:verify(pin)
    base=root/'resumption-r1';d=base/'science/memo-component-v1';d.mkdir()
    x=ET.parse(old/'case.xml').getroot()
    set_info(x,'RANSWallNativeLookup','TensorExactV1');set_info(x,'RANSWallScalarEvaluation','Off');write(x,d/'case.xml')
    request=load(old/'view/request.json');request.update(operation='scalar_value_checks',memo=True)
    atomic(d/'request.json',request)
    runtime_path=base/'science/runtime-build-memo-v1.json';runtime=load(runtime_path)
    probe=verify(next(e['wrapper']for e in runtime['executables']if Path(e['wrapper']['path']).name=='RANSWallResolutionProbe'))
    atomic(d/'source-bindings.json',dict(inputs=pins,runtime=identity(runtime_path),
        distinction='new exact scalar values for repeated identical points in one const lifted-state domain; full ordered inventory/failure hashes compared'))
    subprocess.run([str(probe),str(d/'case.xml'),str(d/'request.json'),str(d/'native.json')],check=True)
    r=load(d/'native.json');assert r['status']=='PASS' and r['tested_policy']=='MemoExactV1'
    assert len(r['actual_scalar_checks'])==16 and r['evaluation_work']['scalar_point_memo_hits']>0
    print('memo ordered inventory/failure parity PASS',r['evaluation_work']['scalar_point_memo_hits'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
