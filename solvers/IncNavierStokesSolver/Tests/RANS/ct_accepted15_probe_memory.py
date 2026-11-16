#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One fixed diagnostic working-set remedy; no physical history is advanced."""
import argparse,gc,os,resource,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted15_exec import native_run
from ct_accepted15_equivalence import numerical
from ct_accepted15_pair import numerical_receipts
from ct_accepted_pilot import budget_classification
from ct_accepted_open import require
from ct_accepted_feasibility import write_json
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,sha,runtime_identity

def main():
    p=argparse.ArgumentParser()
    for key in ('root','runtime','pair','output','contract'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--sparse-only',action='store_true',help='Explicit partial evidence on the two sparse controls; never qualifies the complete working-set policy')
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';os.environ.pop('NEKTAR_RANS_TIMING_FILE',None)
    pair=n.read(a.pair/'results.json');require(pair['qualification']['status']=='pass','fresh candidate pair required')
    rt=runtime_identity(a.runtime/'RANSPreservationProbe');original=n.read(a.pair/'request-manifest.json')['request'];require(rt==original['probe'],'foreign diagnostic runtime')
    cases=[('4000-q64',a.root/'timestep-v1/4000-q64','frozen','sparse'),('8000-q64',a.root/'timestep-v1/8000-q64','frozen','sparse')]+[(label,a.pair/'0.000125'/label,'budget','scalar') for label in ('q64f16','q64f32')]
    if a.sparse_only:cases=cases[:2]
    inputs={label:dict(request_sha256=sha(folder/(tag+'-request.json')),configuration_sha256=sha(folder/'case.xml'),oracle_packet_sha256=sha(folder/(tag+'.json'))) for label,folder,tag,kind in cases}
    m=RequestManifest(a.output,dict(schema='CT15ProbeWorkingSet65536V1',source=sha(Path(__file__)),contract=sha(a.contract),pair=sha(a.pair/'results.json'),probe=rt,maximum_immutable_rows=65536,inputs=inputs,partial_sparse_only=a.sparse_only,physical_advance=False))
    result=dict(schema='CT15ProbeWorkingSet65536V1',status='incomplete',cases={},maximum_rows=65536,physical_advance=False)
    for label,folder,tag,kind in cases:
        dest=a.output/label;dest.mkdir();root=ET.parse(folder/'case.xml').getroot();set_parameter(root,'RANSWallImmutableMaximumRows',65536);write(root,dest/'case.xml');shutil.copy2(folder/(tag+'-request.json'),dest/'request.json')
        os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(dest/'sparse-receipts.jsonl')
        e=native_run([a.runtime/'RANSPreservationProbe',dest/'case.xml',dest/'request.json',dest/'packet.json'],dest,'probe.log',900,1200)
        os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None)
        old=n.read(folder/(tag+'.json'));new=n.read(dest/'packet.json') if (dest/'packet.json').is_file() else {}
        require(e['exit_code']==0 and not new.get('failure'),'bounded working-set probe failed:'+label)
        require(old['observation']['state']==new['observation']['state'] and old['observation']['identity']==new['observation']['identity'],'work capacity changed physical state/identity')
        key='sparse_backend' if kind=='sparse' else 'physical_scalar_budget'
        require(numerical(old[key])==numerical(new[key]),'work capacity changed numerical arrays/decisions:'+label)
        row=dict(status='pass',execution=e,oracle_packet_sha256=inputs[label]['oracle_packet_sha256'],candidate_packet_sha256=sha(dest/'packet.json'),candidate_configuration_sha256=sha(dest/'case.xml'),all_numerical_arrays_equal=True,maximum_scaled_action_or_state_difference=0.)
        if kind=='scalar':
            before=budget_classification(old,pair['case']);after=budget_classification(new,pair['case']);require(before==after and after['reconstruction_status']=='pass','physical scalar classification changed');row['classification']=after
        else:row['numerical_receipts']=numerical_receipts(dest/'sparse-receipts.jsonl')
        result['cases'][label]=row;write_json(a.output/'progress.json',result);del old,new;gc.collect();print('CT15 probe working set',label,'pass',e['process_cpu_seconds'],e['peak_address_space_KiB'],flush=True)
    result['maximum_measured_AS_bytes']=max(v['execution']['peak_address_space_KiB']*1024 for v in result['cases'].values());result['AS_with_headroom_bytes']=1.25*result['maximum_measured_AS_bytes'];result['memory_admission']=('pass' if result['AS_with_headroom_bytes']<=16*1024**3 else 'fail') if not a.sparse_only else 'unqualified_partial_scope';result['status']='partial' if a.sparse_only else 'pass';result['scope']='diagnostic working-set equivalence and capacity only; complete campaign resource admission remains separate'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
