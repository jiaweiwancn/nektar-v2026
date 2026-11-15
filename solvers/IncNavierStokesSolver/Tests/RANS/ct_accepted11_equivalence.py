#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Native streamed/serialized observer equivalence; no governing evolution."""
import argparse,copy,json,math,os,resource,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted_pilot import replay,budget_classification
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity
from run_cls_epawf_picard import set_parameter,write

POLICY='NativeOwnerVolumeStreamV1'
CASES={'d0-p3','d0-p4','displaced-p3','displaced-p4','baseline-T008','candidate-T020'}

def difference(a,b):
    """All primitive/vector entries, with no silent truncation or mean shift."""
    top=0.;absolute=0.;entries=0
    def walk(x,y,path):
        nonlocal top,absolute,entries
        if isinstance(x,dict):
            require(isinstance(y,dict) and set(x)==set(y),'object coverage:'+path)
            for k in x:walk(x[k],y[k],path+'/'+str(k))
        elif isinstance(x,list):
            require(isinstance(y,list) and len(x)==len(y),'array length:'+path)
            for i in range(len(x)):walk(x[i],y[i],path+'/'+str(i))
        elif isinstance(x,(float,int)) and not isinstance(x,bool):
            require(isinstance(y,(float,int)) and math.isfinite(x) and math.isfinite(y),'nonfinite:'+path)
            d=abs(x-y);absolute=max(absolute,d);top=max(top,d/max(1.,abs(x),abs(y)));entries+=1
        else:require(x==y,'metadata:'+path)
    walk(a,b,'')
    require(top<=2e-10,'observer primitive/vector disagreement:'+str(top))
    return dict(maximum_scaled_difference=top,maximum_absolute_difference=absolute,compared_numeric_entries=entries)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct9','ct10','runtime','old-runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    controls=sealed(a.ct9/'controls-v4');base=sealed(a.ct9/'equivalence-v2');fine=sealed(a.ct10/'candidate-T020-v2')
    fixtures={}
    for name in sorted(CASES-{'baseline-T008','candidate-T020'}):
        f=a.ct9/'controls-v4'/(name+'-sparse');packet=n.read(f/'operator.json');require(sha(f/'operator.json')==controls['checks'][name]['artifacts'][1]['sha256'],'changed small fixture')
        c=n.read(f/'construction.json');state=packet['observation']['state'];pars={k.strip():v.strip() for x in ET.parse(f/'case.xml').findall('CONDITIONS/PARAMETERS/P') for k,v in [x.text.split('=',1)]};dt=float(pars['TimeStep'])
        # Static frozen-field fixtures do not pretend to be accepted histories.
        req=replay(packet,state,dt,[],export_physical_scalar_budget=True,export_frozen_scalar_reconciliation=True,previous_state=state,stage_history=state,actual_flow_stage=state,actual_pre_k_stage=state)
        fixtures[name]=dict(directory=f,request=req,construction=c,scope='static frozen-field observer equivalence; no physical budget-accuracy claim')
    for name,row in [('baseline-T008',base['members']['0.001']),('candidate-T020',fine['members']['0.0005'])]:
        f=Path(row['directory']);r=n.read(f/'case.evolution.json');require(sha(f/'case.evolution.json')==row['native_report_sha256'],'changed sealed strip report')
        for label,_,_ in n.CONTROLS:require(sha(f/label/'budget.json')==row['budgets']['budgets'][label]['packet_sha256'],'changed sealed physical budget')
        fixtures[name]=dict(directory=f,construction=n.read(f/'construction.json'),scope='actual consecutive accepted states and frozen PIMPLE stage snapshots',report_sha256=sha(f/'case.evolution.json'))
    manifest=RequestManifest(a.output,dict(schema='CT11ObserverEquivalenceV1',contract=sha(a.contract),source=sha(Path(__file__)),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),old_runtime=runtime_identity(a.old_runtime/'RANSPreservationProbe'),policy=POLICY,fixtures={k:dict(directory=str(v['directory']),configuration_sha256=sha(v['directory']/'case.xml'),scope=v['scope']) for k,v in fixtures.items()},controls=n.CONTROLS))
    out=dict(schema='CT11ObserverEquivalenceV1',status='incomplete',policy=POLICY,checks={})
    for name,f in fixtures.items():
        rows={};dest=a.output/name;dest.mkdir();c=f['construction']
        for label,orders,face in n.CONTROLS:
            folder=dest/label;folder.mkdir();root=ET.parse(f['directory']/'case.xml').getroot()
            for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),orders,strict=True):set_parameter(root,k,v)
            write(root,folder/'case.xml')
            if 'request' in f:
                req=copy.deepcopy(f['request']);req['physical_face_order']=face
                old,oe=probe(a.old_runtime,folder,req,'old');require(oe['exit_code']==0,'sealed old observer failed:'+name+':'+label)
                old_packet_sha=sha(folder/'old.json');old_execution=oe
            else:
                req=n.read(f['directory']/label/'budget-request.json');old=n.read(f['directory']/label/'budget.json');old_packet_sha=sha(f['directory']/label/'budget.json');old_execution=None
            req['export_basis']=False;req['export_faces_only']=False;req['physical_budget_volume_policy']=POLICY
            q,e=probe(a.runtime,folder,req,'stream');require(e['exit_code']==0,'stream observer failed:'+str(q.get('failure')))
            require(q['observation']['state']==req['state'] and q['observation']['identity']==old['observation']['identity'] and q['reference']==old['reference'] and q['fixed_boundary']==old['fixed_boundary'],'observer changed state/reference/boundary/identity')
            b=copy.deepcopy(q['physical_scalar_budget']);meta=b.pop('volume_observer');legacy=copy.deepcopy(old['physical_scalar_budget']);legacy.pop('volume_observer',None);d=difference(legacy,b)
            require(meta['policy']==POLICY and meta['materialized_global_volume_rows'] is False and meta['inner_points']>0 and meta['outer_points']>0,'missing streamed volume identity')
            receipt=q['initialization_receipt'];require(receipt['open_explicit_implicit_agreement'],'observer live boundary mismatch')
            classification=budget_classification(q,c);require(set(classification['equations'])=={'inner_k','outer_k','outer_epsilon','whole_k'} and set(classification['reconstruction'])=={'k','outer_epsilon'} and classification['reconstruction_status']=='pass','all-row reconstruction failed')
            provenance=None
            if label=='native':
                req['physical_budget_volume_policy']='SerializedVolumeV1';z,ze=probe(a.runtime,folder,req,'serialized');require(ze['exit_code']==0,'rebuilt legacy observer failed')
                before=copy.deepcopy(z['physical_scalar_budget']);zm=before.pop('volume_observer');oracle=difference(legacy,before);require(zm['native_rows_provenance_fnv1a64']==meta['native_rows_provenance_fnv1a64'] and zm['inner_points']==meta['inner_points'] and zm['outer_points']==meta['outer_points'] and zm['physical_weights']==meta['physical_weights'],'native row/owner/weight/order provenance changed');provenance=dict(status='pass',legacy=zm,stream=meta,old_binary=oracle,execution=ze,packet_sha256=sha(folder/'serialized.json'))
            row=dict(status='pass',all_primitive_difference=d,volume_observer=meta,native_provenance=provenance,classification=classification,execution=e,old_execution=old_execution,old_packet_sha256=old_packet_sha,packet_sha256=sha(folder/'stream.json'),request_sha256=sha(folder/'stream-request.json'),state_sha256=exact_vector_hash(req['state']),previous_state_sha256=exact_vector_hash(req['previous_state']),flow_stage_sha256=exact_vector_hash(req['actual_flow_stage']),pre_k_stage_sha256=exact_vector_hash(req['actual_pre_k_stage']),configuration_sha256=sha(folder/'case.xml'),scope=f['scope'],raw_state_unchanged=True,free_prescribed={k:{'free':len(v['free_rows']),'prescribed':len(v['prescribed_or_dependent_rows'])} for k,v in b['reconciliation']['assembled'].items()},dependent_epsilon_storage_maximum=max(abs(x) for x in b['reconciliation']['regions']['outer_epsilon']['backward_storage']))
            rows[label]=row;out['checks'][name]=rows;write_json(a.output/'progress.json',out);print('CT11 observer',name,label,'pass',flush=True)
        require(set(rows)=={x[0] for x in n.CONTROLS},'observer control coverage')
    require(set(out['checks'])==CASES,'fixture coverage');require(any(q['dependent_epsilon_storage_maximum']>0 for name,rows in out['checks'].items() if name.startswith(('baseline','candidate')) for q in rows.values()),'dependent epsilon storage test has no signal')
    out['status']='pass';write_json(a.output/'results.json',out);manifest.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
