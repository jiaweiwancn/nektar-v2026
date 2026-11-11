#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Fresh, zero-physical-step same-size first-knot representation control."""
import argparse
import os
from pathlib import Path
import platform
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import CampaignLock,atomic
from run_cls_accepted_awf import case_xml
from run_cls_awf_late_ivp import native
from run_cls_epawf_picard import THREADS,set_info,write,sha,require,load


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('solver','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();require(platform.system()=='Linux','native Linux initialization required')
    a.output.mkdir(parents=True,exist_ok=False);lock=CampaignLock(a.output);os.environ.update(THREADS)
    xml=case_xml('Accepted-AWF-05G',a.output/'initial',0,companion='regularized',
        companion_reason='same-size first-cell representation; faithful covariance blocker retained',
        inner_refinement='cutoff-v1',first_cell='half-first-v1')
    tree=ET.parse(xml).getroot();set_info(tree,'RANSWallEvaluationReuse','ExactStepV1');write(tree,xml)
    definition=load(xml.parent/'definition.json');definition['xml_sha256']=sha(xml);atomic(xml.parent/'definition.json',definition)
    q,rc=native(a.solver,xml,1800,lock)
    physical=q.get('accepted_awf_final') or q['accepted_awf_initial']
    require(len(q['state'])==2861 and q['steps']==0 and q['implicit_callback_count']==0 and not q['pimple_attempts'],
        'zero-time control changed native dimensions/schedule')
    atomic(a.output/'summary.json',dict(status='passed',unknowns=len(q['state']),
        native_dimensions=physical['native_degrees_of_freedom'],physical_steps=0,implicit_callbacks=0,
        solver_sha256=sha(a.solver),receipt=rc,initialization_work=q['initialization_lift_calls_iterations'],
        checkpoint=q['checkpoint'],checkpoint_sha256=sha(Path(q['checkpoint'])/'State.bin')))


if __name__=='__main__':main()
