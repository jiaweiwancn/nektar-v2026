#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Focused native exact-reuse contracts on the fixed accepted/active fixtures."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import os
import xml.etree.ElementTree as ET
from diagnose_cls_accepted_covariance import checkpoint,common_packet,run_probe
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import THREADS,set_info,write,sha,require


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--probe',type=Path,required=True)
    p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--interface-xml',type=Path,required=True);a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(parents=True,exist_ok=False)
    def case(name):
        f=a.fixtures/name;cp=f/'start';q=checkpoint(cp)
        packet=dict(common_packet(cp),operation='evaluation_reuse_checks',state=q['coupled'],history=q['coupled'],
            time=q['time'],alpha=0.,inner_stress_model='Linear' if name=='BudgetH2' else 'SaturatedR2')
        d,rc=run_probe(a.probe,f/'production.xml',packet,a.output/name)
        require(d['passed'] and d['governing_solves']==d['physical_commits']==0,'reuse contract failed')
        return name,dict(checks=d,receipt=rc)
    out=dict(status='running',probe_sha256=sha(a.probe),cases={});atomic(a.output/'summary.json',out)
    with ThreadPoolExecutor(max_workers=2) as pool:
        for name,result in pool.map(case,['BudgetH2','UpstreamOne']):out['cases'][name]=result;atomic(a.output/'summary.json',out)
    tree=ET.parse(a.interface_xml).getroot();set_info(tree,'RANSWallEvaluationReuse','ExactStepV1');write(tree,a.output/'interface.xml')
    d,rc=run_probe(a.probe,a.output/'interface.xml',dict(operation='interface_provenance',inner_stress_model='SaturatedR2'),a.output/'interface-cached')
    require(d['passed'] and d['tangential_radial_neighbors'] and d['incompatible_side_rejections']==4,'cached interface provenance failed')
    out.update(status='passed',interface=dict(checks=d,receipt=rc));atomic(a.output/'summary.json',out)

if __name__=='__main__':main()
