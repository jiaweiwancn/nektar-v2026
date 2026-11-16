#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only native H25 value/inventory checks on the sealed small H16 state."""
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
    for p in pins:verify(p)
    d=root/'science/scalar-value-component-v1';d.mkdir();x=ET.parse(old/'case.xml').getroot()
    set_info(x,'RANSWallNativeLookup','TensorExactV1');set_info(x,'RANSWallScalarEvaluation','ValueExactV1');write(x,d/'case.xml')
    request=load(old/'view/request.json');request.update(operation='scalar_value_checks')
    atomic(d/'request.json',request);atomic(d/'source-bindings.json',dict(inputs=pins,source='H16 sealed remote-heavy-manifest; completed producer/history retained',
        current_runtime=identity(root/'science/runtime-build-scalar-v1.json')))
    subprocess.run([str(root/'runtime-scalar-v1/RANSWallResolutionProbe'),str(d/'case.xml'),str(d/'request.json'),str(d/'result.json')],check=True)
    assert load(d/'result.json')['status']=='PASS'

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
