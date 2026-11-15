#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Locate/retest the exact evolved aggregate q64 allocation, separate from scalars."""
import argparse,os,resource,shutil,copy,gc,math
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def equal_primitives(a,b):
    if isinstance(a,dict):
        require(isinstance(b,dict) and a.keys()==b.keys(),'primitive dictionary coverage changed')
        return max((equal_primitives(a[k],b[k]) for k in a),default=0.)
    if isinstance(a,list):
        require(isinstance(b,list) and len(a)==len(b),'primitive list coverage changed')
        return max((equal_primitives(x,y) for x,y in zip(a,b,strict=True)),default=0.)
    if type(a) in (int,float) and type(b) in (int,float):
        require(math.isfinite(a) and math.isfinite(b),'nonfinite observer primitive');return abs(a-b)
    require(a==b,'non-numeric primitive changed');return 0.

def small_controls(a):
    from ct_accepted7_top import make as make_zero,request as zero_request
    from ct_accepted9_controls import enable
    from ct_accepted8_reference import column
    oldrt=a.ct11/'build-v2/runtime';case=sealed(a.ct11/'target-v1')['case'];out={}
    for label in ('d0','displaced'):
        dest=a.output/('small-'+label)
        if label=='d0':make_zero(dest,a.runtime,degree=3);enable(dest,True);req=zero_request({})
        else:
            c=copy.deepcopy(case);c.update(Z=.5,x=[0.,.125,.25],degree=3,pressure_degree=1)
            col=column(c['K0'],c['H'],c['roughness_input'],c['nu'],c['displacement'])
            c['knots']=sorted(set([0.,c['H']/4,c['H']/2,col.cutoff,c['roughness_input'],c['H'],.3,.5]))
            n.make(dest,a.runtime,c,.001,0);enable(dest);req=n.request(c,True,False)
        os.environ.pop('NEKTAR_RANS_PROBE_OBSERVER',None)
        old,oe=probe(oldrt,dest,req,'legacy');require(oe['exit_code']==0,'supported legacy small fixture failed')
        os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1'
        new,ne=probe(a.runtime,dest,req,'receipt');require(ne['exit_code']==0,'supported native-receipt fixture failed')
        checks={k:equal_primitives(old[k],new[k]) for k in ('observation','initialization_receipt','test_basis','native_constitutive','open_boundary','open_frozen_systems')}
        require(max(checks.values())<=2e-10 and old['physical_admissibility']==new['physical_admissibility'],'small old/new receipt/operator mismatch')
        out[label]=dict(status='pass',defects=checks,legacy_execution=oe,new_execution=ne,legacy_packet_sha256=sha(dest/'legacy.json'),new_packet_sha256=sha(dest/'receipt.json'),configuration_sha256=sha(dest/'case.xml'))
    os.environ['NEKTAR_RANS_PROBE_OBSERVER']='UnknownLifecycleV0'
    bad,e=probe(a.runtime,dest,req,'wrong-policy');require(e['exit_code']!=0 and 'unknown probe observer lifecycle' in str(bad.get('failure')),'unknown observer policy was not rejected')
    out['unknown_policy']=dict(status='pass',failure=bad.get('failure'),execution=e)
    os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';return out

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct12','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--observer-policy',choices=('LegacyTreeV1','NativeAuditReceiptV1'),default='LegacyTreeV1')
    p.add_argument('--ct10',type=Path);p.add_argument('--ct11',type=Path)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_PHASES']='1'
    prior=sealed(a.ct12/'pressure-v5');m=RequestManifest(a.output,dict(schema='CT13EvolvedAggregateMemoryV2',source=sha(Path(__file__)),contract=sha(a.contract),prior=sha(a.ct12/'pressure-v5/results.json'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),observer_policy=a.observer_policy,address_space_GiB=16))
    if a.observer_policy=='NativeAuditReceiptV1':
        require(a.ct10 and a.ct11,'receipt qualification requires native history/observation bindings');small=small_controls(a)
    else:os.environ.pop('NEKTAR_RANS_PROBE_OBSERVER',None);small={}
    controls={}
    for order in (48,64):
        src=a.ct12/'pressure-v4'/('held-q'+str(order));dest=a.output/('q'+str(order));dest.mkdir();shutil.copy2(src/'case.xml',dest/'case.xml')
        old=next(q for q in prior['held_quadrature'] if q['orders'][0]==order)
        require(sha(src/'loads.json')==old['packet_sha256'],'foreign sealed quadrature packet')
        req=n.read(src/'loads-request.json');require(req.get('export_basis',False)==False and req['ct10_export_points']==False,'unexpected original full-export request')
        packet,e=probe(a.runtime,dest,req,'loads');log=dest/'loads.log'
        phases=[x for x in log.read_text().splitlines() if 'RANS probe phase' in x] if log.exists() else []
        row=dict(execution=e,request_sha256=sha(src/'loads-request.json'),configuration_sha256=sha(src/'case.xml'),packet_sha256=sha(dest/'loads.json'),phase_trace=phases,status='capacity_pass' if e['exit_code']==0 else 'blocked_capacity',failure=packet.get('failure'),state_sha256=exact_vector_hash(req['state']),scope='same actual evolved state/configuration/request; no governing steps')
        if order==48:
            require(e['exit_code']==0,'supported q48 old/new control failed')
            oldpacket=n.read(src/'loads.json');diffs={}
            for part in ('actual','analytic'):
                left=packet['ct10_momentum'][part];right=oldpacket['ct10_momentum'][part];require(left.keys()==right.keys(),'momentum component coverage changed')
                for k in left:diffs[part+'/'+k]=float(np.max(abs(np.asarray(left[k])-right[k]),initial=0.))
            require(max(diffs.values())<=2e-10 and packet['observation']['state']==oldpacket['observation']['state'],'supported old/new aggregate rows differ')
            require(packet['observation']['audit_alpha_zero']==oldpacket['observation']['audit_alpha_zero'],'mandatory native diagnostics changed')
            observation_defect=equal_primitives(oldpacket['observation'],packet['observation']);require(observation_defect<=2e-10,'observation rows/weights/owners/derivatives changed')
            row.update(all_row_defects=diffs,all_observation_primitive_defect=observation_defect,mandatory_diagnostics_exact=True,effective_input_policy=packet['ct10_momentum'].get('input_policy'))
        elif e['exit_code']!=0:
            require(packet.get('failure')=='std::bad_alloc','unexpected exact q64 failure')
        if order==64 and e['exit_code']==0 and a.observer_policy=='NativeAuditReceiptV1':
            from ct_accepted12_pressure import Frozen
            target=sealed(a.ct11/'target-v1');member=target['members']['0.0005'];report=n.read(Path(member['directory'])/'case.evolution.json')
            require(sha(Path(member['directory'])/'case.evolution.json')==member['native_report_sha256'],'changed actual native report')
            att=report['pimple_attempts'][-1];require(sha(Path(att['flow_snapshot']))==att['flow_snapshot_sha256'],'changed actual frozen history')
            initial=n.read(a.ct10/'target-preflight-v1/target/readiness.json');weights=n.inventory(target['case'])[2]
            fr=Frozen(n.read(att['flow_snapshot']),initial,target['case'],weights,len(weights))
            meta=prior['stage_loads']['0.0005/accepted']['metadata'];require(sha(Path(meta['path']))==meta['sha256'],'foreign native comparison load')
            base=n.read(Path(meta['path']))['ct10_momentum'];mom=packet['ct10_momentum']
            require(mom['maximum_reconstruction_error']<=2e-10,'q64 native row reconstruction')
            change=np.asarray(mom['total_actual'])-base['total_actual'];v=fr.image(change,'CT13-q64');tight=fr.image(change,'CT13-q64-tight',2e-14)
            err=fr.measure(v);repeat=fr.measure(v-tight);require(err['maximum_scaled_value']<1e-8 and repeat['maximum_scaled_value']<1e-8,'q64 pressure-image uncertainty')
            np.save(a.output/'q64-pressure-image.npy',v)
            row.update(pressure_image=err,tighter_repeat=repeat,pressure_image_sha256=sha(a.output/'q64-pressure-image.npy'),maximum_native_load_change=float(np.max(abs(change),initial=0.)),all_row_reconstruction=mom['maximum_reconstruction_error'],actual_frozen_snapshot_sha256=att['flow_snapshot_sha256'],solves=fr.receipts)
            del fr;gc.collect()
        controls[str(order)]=row;write_json(a.output/'progress.json',dict(controls=controls));print('CT13 memory',order,row['status'],phases[-1:] if phases else [],flush=True)
        if order==48:del oldpacket
        del packet;gc.collect()
    result=dict(schema='CT13EvolvedAggregateMemoryV2',status=controls['64']['status'],observer_policy=a.observer_policy,small_controls=small,controls=controls,historical_CT12='q64 failure unchanged; four passing scalar observer controls are separate',physical_continuation=False)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
