#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Zero-step paired represented/analytical flat residuals; same fixed mesh."""
import argparse,os,resource
from pathlib import Path
from cao_tamura_h2d_state import load
from cao_tamura_budget_h2g import case
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import THREADS,sha
from continue_cls_epawf_pimple import atomic

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for n in ('probe','preflight','h1','output'):p.add_argument('--'+n,type=Path,required=True)
 a=p.parse_args();os.environ.update(THREADS);pre=load(a.preflight);a.output.mkdir(exist_ok=False);out=dict(version='H2gPairedFlatSpatialV1',results={},governing_solves=0)
 for model in ('cls','standard'):
  xml=case(a.output/model,a.h1,pre,flat=True,model=model);d,r=run_probe(a.probe,xml,dict(operation='h2d_flat',inner_stress_model='Linear',expected_model='CLSHighRe' if model=='cls' else 'StandardKEpsilon',xml_sha256=sha(xml)),a.output/(model+'-probe'));out['results'][model]=dict(receipt=r,diagnostic=load(Path(r['output'])),model=load(Path(r['output']+'.model.json')),peak_RSS_KiB_cumulative=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss);atomic(a.output/'diagnosis.json',out)
if __name__=='__main__':main()
