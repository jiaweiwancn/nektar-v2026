#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only native sufficient-positivity prototype on a qualified saved state."""
import argparse,os,subprocess
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import identity,verify,atomic
from cls_epawf_m3_cases import set_info,write


def run(root):
    base=root/'resumption-r1';d=base/'science/positive-prototype-v1';d.mkdir()
    old=root/'science/W1-lookup-v1/monitor';q=load(old/'step-20.qualification-commit.json')
    request=dict(operation='h25_positive_checks',inner_stress_model='Linear',import_role='QualifiedV1',
        native_report=q['report']['path'],native_report_sha256=q['report']['sha256'],
        native_checkpoint=str(verify(q['checkpoint']).parent),checkpoint_sha256=q['checkpoint']['sha256'],
        completed_sidecar=q['callback'],qualified_cursor=identity(old/'step-20.qualification-commit.json'))
    x=ET.parse(verify(load(root/'inputs/baseline.json')['own_cp0']['xml'])).getroot()
    for k,v in dict(RANSWallNativeLookup='TensorExactV1',RANSWallScalarEvaluation='Off',RANSWallAssemblyScatter='Off',
        RANSWallRetention='Off',RANSWallMonitoringPolicy='Off').items():set_info(x,k,v)
    write(x,d/'case.xml');atomic(d/'request.json',request)
    runtime=load(base/'science/runtime-build-scatter-v1.json')
    probe=verify(next(e['wrapper']for e in runtime['executables']if Path(e['wrapper']['path']).name=='RANSWallResolutionProbe'))
    subprocess.run([str(probe),str(d/'case.xml'),str(d/'request.json'),str(d/'native.json')],check=True)
    r=load(d/'native.json');assert r['status']=='PASS_READ_ONLY_PROTOTYPE'
    atomic(d/'decision.json',dict(status='READ_ONLY_PROTOTYPE_NOT_PRODUCTION_GUARD_REPLACEMENT',native=identity(d/'native.json'),
        original_scientific_checks='all retained; zero guards skipped',
        reason='typed certified lower bounds and actually evaluated minima are distinct; acceptance/correction diagnostics still require exact observations',
        limitation='moving-point cost coverage uses measured native row/contraction work on saved evolving states; does not establish all-in production benefit'))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
