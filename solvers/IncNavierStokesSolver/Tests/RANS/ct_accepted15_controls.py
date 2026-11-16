#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT15 immutable-policy controls and clock-off same-dt training workloads."""
import argparse,copy,json,os,resource,shutil,time,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_controls import enable,compare_packets,maximum
from ct_accepted9_pressure_run import sealed
from ct_accepted14_component import ordering as baseline_ordering
from ct_accepted14_pair import numerical_receipts
from ct_accepted13_orderings import native_run
from ct_accepted13_memory import equal_primitives
from ct_accepted_open import probe,require
from ct_accepted_pilot import replay
from compare_cls_epawf_pimple import same_history
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha
from ct_accepted_feasibility import write_json

POLICY='ImmutableTensorV1'
COMPONENT='CT15-ImmutableTensorV1'
def ordering(folder):
    baseline_ordering(folder);policy(folder)

def policy(folder,value=POLICY):
    root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallImmutableEvaluation',value);write(root,folder/'case.xml')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct14','runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--phase',choices=['small','training','baseline','diagnostic'],required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';os.environ.pop('NEKTAR_RANS_TIMING_FILE',None)
    ct9=a.ct14.parent/'ct-accepted-9-20261114'
    m=RequestManifest(a.output,dict(schema='CT15ControlsV1',phase=a.phase,source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),baseline=str(a.ct14),policy=POLICY,detailed_timing=False))
    out=dict(schema='CT15ControlsV1',status='incomplete',phase=a.phase,cases={})
    def evolve(label,c,steps,parent=None,template=None,value=POLICY,reject=0,dt=1e-6):
        folder=a.output/label;n.make(folder,a.runtime,c,dt,steps,parent,template);enable(folder);ordering(folder);policy(folder,value)
        root=ET.parse(folder/'case.xml').getroot();set_parameter(root,'RANSWallPimpleHardAuditRejectStep',reject);write(root,folder/'case.xml')
        os.environ['NEKTAR_RANS_SPARSE_RECEIPTS']=str(folder/'sparse-receipts.jsonl')
        if a.phase=='diagnostic':os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl')
        e=native_run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log',1200,1500)
        report=folder/'case.evolution.json';q=n.read(report) if report.exists() else dict(failure=(folder/'solver.log').read_text()[-3000:])
        os.environ.pop('NEKTAR_RANS_SPARSE_RECEIPTS',None)
        return folder,q,e
    if a.phase=='small':
        old=sealed(a.ct14/'component-v1');require(old['status']=='pass','baseline component not sealed')
        for label in ('d0-p3','d0-p4','displaced-p3','displaced-p4'):
            src=a.ct14/'component-v1'/(label+'-VectorCOLAMDV1');oracle=n.read(src/'operator.json');c=n.read(ct9/'controls-v4'/(label+'-sparse')/'construction.json');req=n.read(src/'operator-request.json')
            dest=a.output/label;dest.mkdir();shutil.copy2(src/'case.xml',dest/'case.xml');policy(dest)
            q,e=probe(a.runtime,dest,req,'operator');require(e['exit_code']==0,'immutable operator failed:'+str(q.get('failure')))
            require(COMPONENT in q['observation']['identity'],'missing new execution identity')
            result=compare_packets(oracle,q,c)
            out['cases'][label]=dict(status='pass',checks=result,execution=e,oracle_sha256=sha(src/'operator.json'),packet_sha256=sha(dest/'operator.json'))
            if label=='displaced-p3':
                request=replay(q,q['observation']['state'],0.,[],alpha=0.,history=q['observation']['state'])
                request['export_immutable_checks']=True
                os.environ['NEKTAR_RANS_VERIFY_OWNER_INDEX']='1'
                checks,execution=probe(a.runtime,dest,request,'dependencies')
                os.environ.pop('NEKTAR_RANS_VERIFY_OWNER_INDEX',None)
                require(execution['exit_code']==0 and checks['passed'] and checks['immutable_only'],'immutable dependency/owner controls')
                out['dependencies']=dict(status='pass',checks=checks,execution=execution,packet_sha256=sha(dest/'dependencies.json'))
            write_json(a.output/'progress.json',out);print('CT15 small',label,'pass',flush=True)
        c=n.read(ct9/'controls-v4/one/construction.json');template=ct9/'controls-v4/one/case.xml'
        full,fr,fe=evolve('full',c,4,template=template)
        two,tr,te=evolve('two',c,2,template=template)
        cont,cr,ce=evolve('continued',c,2,Path(tr['checkpoint']),template)
        rollback,rr,re=evolve('rollback',c,4,template=template,reject=3)
        oracle,orr,oe=evolve('oracle',c,4,template=template,value='Off')
        require(all(e['exit_code']==0 for e in (fe,te,ce,oe)) and all(not q['failure'] for q in (fr,tr,cr,orr)),'small native update failure')
        require(same_history(fr,cr),'genuine candidate restart histories differ')
        require(re['exit_code']!=0 and 'deliberate' in rr['failure'] and same_history(tr,rr),'full rejected-step rollback differs')
        defect=maximum(np.asarray(fr['state'])-orr['state']);require(defect<=1e-8,'small complete update equivalence')
        negatives={}
        for label,parent,value in [('old-to-new',Path(orr['checkpoint']),POLICY),('new-to-old',Path(tr['checkpoint']),'Off'),('unknown',None,'UnknownV1')]:
            folder,q,e=evolve(label,c,0,parent,template,value)
            message=q.get('failure','')+(folder/'solver.log').read_text()[-3000:]
            require(e['exit_code']!=0 and ('signature' if parent else 'immutable') in message.lower(),'wrong identity rejection:'+label)
            negatives[label]=dict(status='pass',execution=e,failure=message)
        out['histories']=dict(status='pass',full_continuation_equal=True,full_rollback_equal=True,oracle_maximum_state_difference=defect,negatives=negatives,executions=[fe,te,ce,re,oe],reports={str(f/'case.evolution.json'):sha(f/'case.evolution.json') for f in (full,two,cont,rollback,oracle)})
    else:
        base=sealed(a.ct14/'static-v1');c=base['case'];out['case']=c
        for steps in ((0,1,2) if a.phase=='training' else (1,)):
            dest,q,e=evolve('steps-'+str(steps),c,steps,dt=.00025,value='Off' if a.phase=='baseline' else POLICY)
            require(e['exit_code']==0 and not q['failure'] and q['steps']==steps,'selected training trajectory failure')
            require(exact_vector_hash(q['accepted_states'][0])=='dd0af02695e58934d6469580aec40249700b404b7f079c0abe8db10c99d3c77a','actual selected initialization changed')
            receipts=numerical_receipts(dest/'sparse-receipts.jsonl')
            out['cases'][str(steps)]=dict(status='pass',steps=steps,execution=e,state_sha256=exact_vector_hash(q['state']),report_sha256=sha(dest/'case.evolution.json'),receipts=receipts,output_bytes=sum(f.stat().st_size for f in dest.rglob('*') if f.is_file()))
            out['cases'][str(steps)]['detailed_timing']=False
            out['cases'][str(steps)]['broad_diagnostic_timing']=a.phase=='diagnostic'
            if steps:
                old=n.read(a.ct14/'pair-v2/0.00025/case.evolution.json');truth=old['accepted_states'][steps]
                difference=maximum(np.asarray(q['state'])-truth);require(difference<=1e-8,'selected update differs from CT14 oracle')
                out['cases'][str(steps)]['maximum_absolute_coefficient_difference']=difference
            write_json(a.output/'progress.json',out);print('CT15 timing-off training',steps,e['process_cpu_seconds'],flush=True)
    out['status']='pass';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
