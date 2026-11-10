#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent physical scalar flux budget; finite defects are not row identities."""
import argparse,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import exact_vector_hash,write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity
from rans_duration_contract import require
from rans_preservation_oracle import Reference
from qualify_rans_preservation import load,replay_request,Q
from diagnose_rans_balanced_spatial import run_probe
from run_cls_epawf_picard import set_parameter,write


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--label',required=True);p.add_argument('--time',type=float,required=True);p.add_argument('--output-label',required=True);p.add_argument('--quadrature-index',type=int,choices=(1,2),default=2);p.add_argument('--face-order',type=int,choices=(16,24,32),default=16);p.add_argument('--runtime-name',default='runtime-budget-v5');p.add_argument('--resume',action='store_true');a=p.parse_args();r=a.root.resolve();runtime=r/a.runtime_name;os.environ.update(LD_LIBRARY_PATH=str(runtime/'lib'),TMPDIR=str(r/'tmp'));folder=r/'campaigns'/a.label;req=load(folder/'request-manifest.json');endpoint=folder/('T-'+format(a.time,'.4f'));require(RequestManifest(folder,req['request'],True).completed(endpoint),'budget_endpoint_not_sealed')
    native=load(endpoint/'case.evolution.json');states=native['accepted_states'];require(len(states)>=2,'budget_needs_consecutive_states');previous=states[-2];require(exact_vector_hash(previous)==native['pimple_attempts'][-1]['incoming_sha256'],'budget_stage_reconstruction_mismatch')
    dt=req['request']['dt'];require(abs(native['audit_times'][-1]-native['audit_times'][-2]-dt)<1e-13,'budget_history_time_mismatch')
    out=r/'physical-budgets'/a.output_label
    intent=dict(schema='PhysicalScalarBudgetRequestV1',configuration={k:req['request'][k] for k in ('model','policy','wall','descriptor')},time=native['time'],dt=dt,orders=Q[a.quadrature_index],face_order=a.face_order,trajectory_request_sha256=req['sha256'],source_report_sha256=sha(endpoint/'case.evolution.json'),source_xml_sha256=sha(endpoint/'case.xml'),state_sha256=exact_vector_hash(native['state']),previous_state_sha256=exact_vector_hash(previous),runtime=runtime_identity(runtime/'RANSPreservationProbe'))
    seal=RequestManifest(out,intent,a.resume)
    if seal.completed(out):return 0 if load(out/'results.json')['status']=='pass' else 1
    shutil.copy2(endpoint/'case.xml',out/'case.xml');xml=ET.parse(out/'case.xml').getroot()
    for function in xml.findall('./CONDITIONS/FUNCTION'):
        if function.get('NAME') in ('ScalarSourceInner','ScalarSourceOuter'):
            for entry in function:
                require(entry.tag=='E' and float(entry.get('VALUE','nan'))==0.,'physical_budget_requires_zero_external_scalar_sources')
    for name,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),Q[a.quadrature_index],strict=True):set_parameter(xml,name,value)
    write(xml,out/'case.xml');seed=load(folder/'initial-packet.json');request=replay_request(seed,[],state=native['state'],time=native['time']);request.update(fixed_boundary=native['initialization_receipt']['live_fixed_boundary'],export_physical_scalar_budget=True,physical_face_order=a.face_order,previous_state=previous);packet,elapsed=run_probe(runtime,out,request);require(exact_vector_hash(packet['observation']['state'])==exact_vector_hash(native['state']),'budget_modified_current_state')
    b=packet.get('physical_scalar_budget');require(b is not None and b['schema']=='PhysicalFrozenStageScalarBudgetV1','budget_probe_unavailable');require(b['face_order']==a.face_order,'physical_face_order_not_consumed');require(b['wall_local_reaction_rate']==0,'prescribed_sink_local_rate_not_zero');ref=Reference(req['request']['model'],req['request']['policy'],req['request']['wall']);scaled={}
    for name,values in b['equations'].items():
        scale=(ref.A/ref.H if name=='outer_epsilon' else ref.K)*ref.us/ref.H*values['volume'];scaled[name]={key:val/scale for key,val in values.items() if key!='volume'}
        x=values['instantaneous_physical_defect']-values['frozen_stage_physical_defect']-values['frozen_operator_splitting'];require(abs(x)<=1e-10,'physical_budget_partition_mismatch')
    result=dict(status='pass',execution='completed',configuration={k:req['request'][k] for k in ('model','policy','wall','descriptor')},time=native['time'],dt=dt,physical_budget=b,scaled=scaled,source_report_sha256=sha(endpoint/'case.evolution.json'),state_sha256=exact_vector_hash(native['state']),previous_state_sha256=exact_vector_hash(previous),history=load(endpoint/'metrics.json')['history_identity'],runtime=runtime_identity(runtime/'RANSPreservationProbe'),orders=Q[a.quadrature_index],packet_sha256=sha(out/'packet.json.gz'),seconds=elapsed,accuracy_classification='physical flux/volume reconstruction completed; finite budget defects require separate quadrature/space/time interpretation and are not certified conservation',scope='no prescribed-row reaction used; no inner epsilon equation; actual dependent boundary storage included; no source/field adjustment')
    write_json(out/'results.json',result);seal.seal(out,['results.json','case.xml','request.json','packet.json.gz','probe.log']);print('physical budget reconstructed; defect accuracy assessed separately');return 0


if __name__=='__main__':raise SystemExit(main())
