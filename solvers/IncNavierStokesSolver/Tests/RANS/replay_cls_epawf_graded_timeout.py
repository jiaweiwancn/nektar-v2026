#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact U05-graded attempt-65 projection replay; no governing solves or Lift."""
import argparse
import math
import os
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET
from continue_cls_epawf_pimple import header
from run_cls_epawf_picard import load,run,require,write_json,sha,state_hash,set_info,set_parameter,write,THREADS
from run_cls_epawf_pimple import identity

PACKET='54740a789b17757d5bc74d42d47018b1a7267056d93fda0293ec9c65215b5222'
CHECKPOINT='4a953092ce4e85c6a6fec9c172ce5eee83bf820c8b0600762637afd106a9abe0'
RAW='17a7a306373898d129f06e7c7a0815ce99608d708fb02d3f07f2d78f983ae659'
HISTORY='b7dab6b2f21a94089bae47f9ebda24b91ec90813f8e0a32ccea7384c0f0d055b'


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--finite-solver',choices=('WorkingSet','RetainedPool'),default='RetainedPool')
    p.add_argument('--seconds',type=float,default=900.);p.add_argument('--full-trace',action='store_true')
    p.add_argument('--total',type=int,default=20000);p.add_argument('--per-add',type=int,default=2000)
    p.add_argument('--separations',type=int,default=2000)
    p.add_argument('--prepare-only',action='store_true');a=p.parse_args();os.environ.update(THREADS)
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    require(sha(a.archive/'attempt.json')==PACKET and sha(a.archive/'checkpoint/State.bin')==CHECKPOINT,'saved failure identities differ')
    old=load(a.archive/'manifest.json')
    for name,digest in old['files'].items():require(sha(a.archive/name)==digest,'immutable artifact changed: '+name)
    q=load(a.archive/'attempt.json');require(state_hash(q['raw'])==RAW and state_hash(q['history'])==HISTORY,'raw candidate or actual callback history differs')
    cp=header(a.archive/'checkpoint',True);require(cp['steps']==64 and state_hash(cp['coupled_state'])==state_hash(q['committed']),'accepted checkpoint differs from attempt base')
    shutil.copytree(a.archive/'checkpoint',a.output/'checkpoint64')
    tree=ET.parse(a.archive/'case.xml').getroot()
    set_info(tree,'RANSWallPimpleKBoundSolver',a.finite_solver)
    set_info(tree,'RANSWallPimpleKBoundTrace','Full' if a.full_trace else 'Summary')
    set_info(tree,'RANSWallPimpleKBoundTraceFile',str(a.output/'trace.jsonl'))
    set_info(tree,'RANSWallPimpleKBoundDumpFile',str(a.output/'pool.json'))
    for name,value in {'Seconds':a.seconds,'Total':a.total,'Add':a.per_add,'Separations':a.separations}.items():
        set_parameter(tree,'RANSWallPimpleKBoundMax'+name,value)
    write(tree,a.output/'case.xml')
    request=load(a.archive/'replay-input.json');request['native_checkpoint']=str(a.output/'checkpoint64')
    request['witnesses']=[q['k_bound'][name]['coordinate'] for name in ('raw_minimum','corrected_minimum')]
    write_json(a.output/'input.json',request)
    manifest=dict(status='prepared',archive=str(a.archive.resolve()),packet_sha256=PACKET,checkpoint_sha256=CHECKPOINT,
        raw_sha256=RAW,history_sha256=HISTORY,transport_A_sha256=q['matrix_sha256'],
        projection_metric='reconstructed native PicardKMass Mff on saved free rows, not transport A_k',
        finite_solver=a.finite_solver,controls=dict(per_add=a.per_add,total=a.total,separations=a.separations,seconds=a.seconds),
        xml_sha256=sha(a.output/'case.xml'),input_sha256=sha(a.output/'input.json'),governing_solves=0,ranks=1,threads=THREADS)
    write_json(a.output/'manifest.json',manifest)
    if a.prepare_only:return
    manifest['probe']=identity(a.probe.resolve())
    receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],a.output,'probe.log',math.ceil(a.seconds+600))
    manifest['receipt']=receipt;write_json(a.output/'manifest.json',manifest)
    require(receipt['exit_code']==0,'native projection probe failed: inspect probe.log')
    r=load(a.output/'result.json');b=r['k_bound'];require(state_hash(r['raw'])==RAW and r['history_sha256']==HISTORY,'replay input changed')
    require(r['initialization_lift_calls_iterations']==[0.,0.] and not b['projection_only'] and r['pre_k_inventory']['admissible'],'unexpected entry waiver or h0 construction')
    require(r['non_k_change']==0 and r['essential_k_change']==0 and r['shared_cg_jump']<1e-12,'protected field/CG invariance failed')
    lo,hi=map(int,r['offsets'][2:4]);require(state_hash(r['raw'][:lo]+r['raw'][hi:])==state_hash(r['corrected'][:lo]+r['corrected'][hi:]),'non-k bits changed')
    if b['success']:
        require(b['primal_floor_defect']==0 and b['ground_floor_defect']==0 and b['current_inventory_passed'] and b['ground_passed'],'physical floor/inventory check failed')
        require(r['scalar_domain']['admissible'] and not r['complete_audit_failure'] and r['complete_hard_admissible'],'complete hard audit failed')
        require(b['stationarity_mass_dual']<1e-10 and b['dual_feasibility_defect']<1e-12,'finite-QP KKT check failed')
        if a.finite_solver=='RetainedPool':require(b['finite_pool_complete'],'finite constraint pool not solved')
    manifest.update(status='repaired_and_hard_audited' if b['success'] else 'projection_rejected',bound=b,
        result_sha256=sha(a.output/'result.json'),corrected_sha256=state_hash(r['corrected']),
        hard_admissible=r.get('complete_hard_admissible',False),audit_failure=r['complete_audit_failure'],
        replay_context=r['replay_context'],non_k_change=r['non_k_change'],essential_k_change=r['essential_k_change'],cg_jump=r['shared_cg_jump'])
    write_json(a.output/'manifest.json',manifest);print(a.output/'manifest.json')


if __name__=='__main__':main()
