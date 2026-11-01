#!/usr/bin/env python3
"""Unscored full-history continuation of a retained, valid R1 checkpoint."""
import argparse
import json
import os
from pathlib import Path
import shutil
import sys
import xml.etree.ElementTree as ET

sys.dont_write_bytecode=True
from cao_tamura_experiment import ROOT,THREADS,allocation,environment,execute,run,sha,write_json
from validate_cao_tamura_r1_map import audit
from generate_m4 import parameter


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('parent','checkpoint','build','output'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--steps',type=int,required=True);p.add_argument('--cpus',required=True)
    p.add_argument('--ranks',type=int,default=4);p.add_argument('--timeout',type=int,default=900)
    p.add_argument('--debug',action='store_true',help='Log nonlinear and independent action residuals')
    a=p.parse_args()
    if not sys.platform.startswith('linux'):raise ValueError('Linux is authoritative')
    allocation(a.ranks,a.cpus);os.environ.update(THREADS)
    os.sched_setaffinity(0,set(map(int,a.cpus.split(','))))
    build=environment(a.build);out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    shutil.copy2(__file__,out/'producer.py')
    original=json.loads((a.parent/'run.json').read_text())
    if original['case_sha256']!=sha(a.parent/'case.xml') or original['solver_sha256']!=sha(build/'IncNavierStokesSolver'):
        raise ValueError('Parent source/binary binding differs')
    result=dict(status='RUNNING',scope='full-history numerical-conditioning control; not developed or calibrated',
        parent=str(a.parent),parent_run_sha256=sha(a.parent/'run.json'),
        checkpoint=str(a.checkpoint),producer_sha256=sha(out/'producer.py'),
        checkpoint_files={str(f.relative_to(a.checkpoint)):sha(f) for f in a.checkpoint.rglob('*') if f.is_file()}
            if a.checkpoint.is_dir() else {a.checkpoint.name:sha(a.checkpoint)})
    write_json(out/'evidence.json',result)
    try:
        restart=ROOT/'solvers/IncNavierStokesSolver/Examples/NeutralABL/restart.py'
        execute([sys.executable,restart,a.parent/'case.xml',a.checkpoint,out/'case.xml','--steps',a.steps],out,'prepare')
        if a.debug:
            tree=ET.parse(out/'case.xml');parameter(tree.getroot(),'RANSDebug',1)
            tree.write(out/'case.xml',encoding='utf-8',xml_declaration=True)
        manifest=json.loads((a.parent/'manifest.json').read_text())
        if a.debug:manifest.setdefault('numerics',{})['RANSDebug']=1
        manifest['case_sha256']=sha(out/'case.xml');manifest['config']['steps']=a.steps
        write_json(out/'manifest.json',manifest)
        result['run']=run(out/'case.xml',build,a.ranks,a.cpus,a.timeout)
        result.update(status='PASS',audits=audit(out))
    except Exception as error:
        result.update(status='FAIL',error=repr(error));write_json(out/'evidence.json',result);raise
    write_json(out/'evidence.json',result);print('full-history conditioning control PASS',flush=True)


if __name__=='__main__':main()
