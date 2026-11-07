#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Exact native attempt-73 projection replay; zero governing transport solves."""
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,run,require,write_json,state_hash,sha,set_info,set_parameter,write,THREADS,ROOT
from run_cls_epawf_pimple import identity

PACKET_SHA='17180bc0162008914722ee85c6b529c69baed2c6708d2ca0c12733c8586bf78d'
RAW_SHA='ec7eb83a29bcfd691eaa0e207eb084cfc5fdd65c824b56adaf5b3e9f3339f367'
CP_SHA='0ebcf803b486d63c42d6ed93db08efd4208d87b98e90dd918e139a77f8c26510'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--original',type=Path,required=True)
    p.add_argument('--uniform',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--total',type=int,default=2000);p.add_argument('--per-add',type=int,default=2000)
    p.add_argument('--separations',type=int,default=2000);p.add_argument('--seconds',type=float,default=900.)
    p.add_argument('--prepare-only',action='store_true');a=p.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False)
    segment=a.original.resolve()/'campaign/trajectory/steps-000065-000080'
    packet=segment/'pimple-states/k-t0.045625000000000027-31c63a118607.json'
    checkpoint=a.uniform.resolve()/'upstream-terminal/checkpoint-last-accepted'
    source=a.uniform.resolve()/'sources-common32/upstream-h0/result.json'
    require(sha(packet)==PACKET_SHA and sha(checkpoint/'State.bin')==CP_SHA,'reviewed artifact identity mismatch')
    q=load(packet);require(state_hash(q['raw'])==RAW_SHA,'reviewed raw vector mismatch')
    require(sha(Path(q['matrix_file']))==q['matrix_sha256'],'frozen transport A identity differs')
    report=load(segment/'case.evolution.json');boundary=load(source)['fixed_external_coefficients']
    shutil.copy2(packet,a.output/'archived-packet.json');shutil.copytree(checkpoint,a.output/'checkpoint72')
    shutil.copy2(Path(q['matrix_file']),a.output/'transport_A.bin')
    write_json(a.output/'archived-native-boundary.json',dict(source=str(source),source_sha256=sha(source),fixed_external_coefficients=boundary))
    xml=a.output/'case.xml';tree=ET.parse(segment/'case.xml').getroot()
    for key,value in {'RANSWallPimpleKBoundMaxAdd':a.per_add,'RANSWallPimpleKBoundMaxTotal':a.total,
        'RANSWallPimpleKBoundMaxSeparations':a.separations,'RANSWallPimpleKBoundMaxSeconds':a.seconds}.items():set_parameter(tree,key,value)
    for key,value in {'RANSWallPimpleKBoundTrace':'Full','RANSWallPimpleKBoundTraceFile':str(a.output/'trace.jsonl'),
        'RANSWallPimpleKBoundDumpFile':str(a.output/'projection-pool.json')}.items():set_info(tree,key,value)
    write(tree,xml)
    request=dict(operation='k_bound',state=q['raw'],pre_k=q['pre_k'],history=q['history'],committed=q['committed'],
        time=q['time'],alpha=q['alpha'],transport_rows=q['transport_rows'],raw_sha256=RAW_SHA,
        native_checkpoint=str(a.output/'checkpoint72'),checkpoint_sha256=CP_SHA,native_signature=report['signature'],
        fixed_external_coefficients=boundary,projection_only=False,require_history_h0=False,inner_stress_model='SaturatedR2')
    write_json(a.output/'input.json',request)
    manifest=dict(schema=1,status='prepared',source_revision=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        source_status=subprocess.check_output(['git','-C',str(ROOT),'status','--porcelain'],text=True).strip(),
        original_packet=dict(path=str(packet),sha256=PACKET_SHA),checkpoint=dict(path=str(checkpoint),sha256=CP_SHA),
        raw_sha256=RAW_SHA,pre_k_sha256=state_hash(q['pre_k']),history_sha256=state_hash(q['history']),
        original_corrected_sha256=state_hash(q['corrected']),original_k_bound=q['k_bound'],
        transport_A=dict(path=str(a.output/'transport_A.bin'),sha256=q['matrix_sha256'],definition='frozen transport A_k; never used as projection metric'),
        xml_sha256=sha(xml),input_sha256=sha(a.output/'input.json'),ranks=1,threads=THREADS,
        controls=dict(total=a.total,per_add=a.per_add,separations=a.separations,seconds=a.seconds),governing_solves=0)
    if not a.prepare_only:
        manifest['probe']=identity(a.probe.resolve())
        receipt=run([str(a.probe.resolve()),'case.xml','input.json','result.json'],a.output,'probe.log',math.ceil(a.seconds+300))
        manifest['execution']=receipt;write_json(a.output/'manifest.partial.json',manifest)
        require(receipt['exit_code']==0,'projection diagnostic infrastructure failed; inspect probe.log')
        r=load(a.output/'result.json');b=r['k_bound']
        require(r['initialization_lift_calls_iterations']==[0.,0.] and not b['projection_only'],'native replay bypassed entry or rebuilt h0')
        require(r['history_sha256']==state_hash(q['history']) and state_hash(r['raw'])==RAW_SHA,'saved history/raw changed')
        require(r['non_k_change']==0 and r['essential_k_change']==0 and r['shared_cg_jump']<1e-12,'invariance or shared CG violation')
        require(r['pre_k_inventory']['admissible'],'unexpected native pre-k incompatibility')
        if a.total==2000 and a.separations==2000:
            require(not b['success'] and b['rejected_total_limit_checks']==1 and b['active_set_iterations']==2000,'old total limit did not reproduce')
            require(b['legacy_guard_counter']==2001 and b['kkt_solves']==2000,'old work differs')
            require(state_hash(r['corrected'])==state_hash(q['corrected']),'old arithmetic/cut sequence changed retained repaired state')
        if b['success']:
            require(r['scalar_domain']['admissible'] and r['complete_audit_failure']=='','repair failed full hard acceptance audit')
            require(b['ground_floor_defect']<1e-12 and b['primal_floor_defect']<1e-12,'repair lacks physical floor feasibility')
            require(b['stationarity_mass_dual']<1e-10 and b['dual_feasibility_defect']<1e-12,'repair optimality check failed')
        manifest.update(status='projection_repaired_and_audited' if b['success'] else 'projection_rejected',
            replay_context=r['replay_context'],bound=b,result_sha256=sha(a.output/'result.json'),
            corrected_sha256=state_hash(r['corrected']),non_k_change=r['non_k_change'],essential_k_change=r['essential_k_change'],
            shared_cg_jump=r['shared_cg_jump'],complete_audit_failure=r['complete_audit_failure'])
        require(all(identity(a.probe.resolve())[k]==manifest['probe'][k] for k in ('sha256','libraries')),'probe runtime changed')
    write_json(a.output/'manifest.json',manifest);print(a.output/'manifest.json')


if __name__=='__main__':main()
