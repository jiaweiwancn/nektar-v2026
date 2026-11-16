#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent q8/12 initial-data screen, never a changed-signature restart."""
import argparse,os
from pathlib import Path
import xml.etree.ElementTree as ET
from run_wj_h4 import THREADS,load
os.environ.update(THREADS)
from wj_h15_observer import verify,atomic,identity
from cls_epawf_m3_cases import set_parameter,write
from seed_wj_h25 import run


def initialize(root):
    base=root/'resumption-r1/science';input_dir=base/'q812-input-v1';input_dir.mkdir()
    source=load(root/'inputs/baseline.json')['own_cp0']['xml'];x=ET.parse(verify(source)).getroot()
    set_parameter(x,'RANSWallTangentialOrder',8);set_parameter(x,'RANSWallNormalOrder',12)
    write(x,input_dir/'case.xml')
    atomic(input_dir/'intent.json',dict(source=source,target=identity(input_dir/'case.xml'),
        native_representation_points=[8,8],governing_quadrature=[8,12],degree=3,
        distinction='quadrature-only target initial construction; no source import under target signature; no evolving accuracy admission'))
    run(root,'seed700-q812-r1-v1',.0025,False,input_dir/'case.xml',base/'runtime-build-scatter-v1.json')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);initialize(p.parse_args().root)
