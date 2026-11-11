#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Capture actual ELF mappings during an immutable, zero-step physical restart."""
import argparse
import os
from pathlib import Path
import subprocess
import time
import numpy as np
from compare_cls_epawf_pimple import same_history
from run_wj_h2 import load, save, sha
from wj_h2_cases import case_xml


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--solver',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    producer=a.root/'study-v1/base-pilot-on/case.evolution.json';reference=load(producer);cp=Path(reference['checkpoint'])
    immutable={str(v):sha(v) for v in (producer,cp/'State.bin')}
    xml=case_xml(a.output/'readback','physical',0,.02,reuse='ExactStepV1',restart=cp)
    mapped=set();start=time.monotonic()
    with (a.output/'native.log').open('w') as f:
        child=subprocess.Popen([str(a.solver),str(xml)],cwd=xml.parent,stdout=f,stderr=subprocess.STDOUT)
        while child.poll() is None:
            try:
                for line in Path(f'/proc/{child.pid}/maps').read_text().splitlines():
                    fields=line.split(maxsplit=5)
                    if len(fields)==6 and fields[5].startswith('/'):
                        path=Path(fields[5]);mapped.add(str(path.resolve()))
            except FileNotFoundError:pass
            time.sleep(.05)
        assert child.wait()==0,'native physical frozen readback failed'
    after=load(xml.with_suffix('.evolution.json'));same_history(reference,after)
    assert not after['pimple_attempts'] and after['initialization_lift_calls_iterations']==[0,0]
    assert np.asarray(reference['fixed_reference'],dtype='<f8').tobytes()==np.asarray(after['fixed_reference'],dtype='<f8').tobytes()
    assert all(sha(v)==h for v,h in immutable.items()),'changed original producer bytes'
    identities=[dict(path=v,bytes=Path(v).stat().st_size,sha256=sha(v)) for v in sorted(mapped) if Path(v).is_file()]
    spatial=[v for v in identities if Path(v['path']).name=='libSpatialDomains.so.5.9.0']
    assert len(spatial)==1 and spatial[0]['sha256']=='0c23e0de180a233bd88f70b265ae9865aeac456e6482f8ff0fa1b4be84d69dab'
    save(a.output/'summary.json',dict(status='PASS',method='actual /proc/main-native-pid/maps throughout frozen restart',
        solver=str(a.solver),solver_sha256=sha(a.solver),loaded_files=identities,LD_LIBRARY_PATH=os.environ.get('LD_LIBRARY_PATH'),
        elapsed_seconds=time.monotonic()-start,governing_attempts=0,initialization_Lift=0,
        primary_gauge_history_time_fixed_reference_exact=True,unchanged_original=immutable,
        preceding_callback_report=dict(path=str(producer),sha256=sha(producer)),
        startup_report=dict(path=str(xml.with_suffix('.evolution.json')),sha256=sha(xml.with_suffix('.evolution.json'))),
        scope='startup alpha=0 audit does not replace the original completed callback'))
    print('PASS: actual loaded source-built library and immutable zero-step physical restart')


if __name__=='__main__':main()
