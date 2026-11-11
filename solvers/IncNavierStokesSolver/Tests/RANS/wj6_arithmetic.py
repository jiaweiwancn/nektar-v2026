#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One bounded offline arithmetic/SVD process; zero governing solves."""
import argparse
from pathlib import Path
from wj3_measure import arithmetic
from wj_immutable import read,write_once,artifacts,file_hash
from run_cls_epawf_picard import require,state_hash


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--report',type=Path,required=True)
    p.add_argument('--measurement',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    d=read(a.report);m=read(a.measurement)['measurements']
    require(state_hash(d['state'])==m['state_sha256'] and d['time']==m['time'],'arithmetic observation differs from native producer')
    inputs=artifacts([a.report,a.measurement]);value=arithmetic(d,m)
    value.update(input_lineage=inputs,code=artifacts([Path(__file__),Path(__file__).with_name('wj3_measure.py')]),governing_solves=0)
    require(all(file_hash(n)==i['sha256'] for n,i in inputs.items()),'arithmetic mutated inputs')
    write_once(a.output,value);print(a.output)


if __name__=='__main__':main()
