#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Same-size geometry identity and strict incompatible restart check (zero steps)."""
import argparse,os,shutil
from pathlib import Path
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import CampaignLock,atomic
from diagnose_cls_accepted_covariance import checkpoint
from cls_awf_native_metrics import flat_partition
from run_cls_awf_late_ivp import native
from run_cls_epawf_picard import THREADS,load,require,sha,set_info,write


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('solver','control','baseline','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output);os.environ.update(THREADS)
    xml=a.control/'initial/case.xml';base=a.baseline/'case.xml';x,z=flat_partition(xml);xb,zb=flat_partition(base)
    require(x==xb and z==[.00125 if v==.0025 else v for v in zb],'control changed more than first knot')
    left,right=[checkpoint(q) for q in (a.control/'initial/checkpoint',a.baseline/'checkpoint')]
    checks=dict(same_coupled_size=len(left['coupled'])==len(right['coupled'])==2861,
        distinct_physical_signature=left['signature']!=right['signature'],both_native_zero_time=left['steps']==right['steps']==0 and left['time']==right['time']==0)
    folder=a.output/'incompatible-same-size';folder.mkdir();shutil.copy2(a.control/'initial/definition.json',folder/'definition.json')
    root=ET.parse(xml).getroot();set_info(root,'RANSWallRestart',str(a.baseline/'checkpoint'));set_info(root,'RANSWallCheckpoint',str(folder/'checkpoint'))
    target=folder/'case.xml';write(root,target);_,rc=native(a.solver,target,1800,lock,False)
    text=(folder/'solver.log').read_text();checks['ordinary_same_size_incompatible_restart_rejected']='checkpoint signature mismatch' in text
    require(all(checks.values()),'first-cell native identity failed: '+str(checks))
    atomic(a.output/'summary.json',dict(status='passed',checks=checks,rejection_receipt=rc,baseline_checkpoint_sha256=sha(a.baseline/'checkpoint/State.bin'),control_checkpoint_sha256=sha(a.control/'initial/checkpoint/State.bin')))


if __name__=='__main__':main()
