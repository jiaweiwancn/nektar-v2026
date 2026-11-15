#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Qualify the distinct CT14 predictor component, complete updates and histories."""
import argparse,copy,json,os,resource,shutil,xml.etree.ElementTree as ET
from functools import lru_cache
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_controls import enable,compare_packets,maximum
from ct_accepted9_pressure_run import sealed
from ct_accepted13_memory import equal_primitives
from ct_accepted_open import probe,run,require
from ct_accepted5_traction import profiles,difference
from compare_cls_epawf_pimple import same_history
from ct_accepted_feasibility import write_json
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha

COMPONENT='CT14-ScotchVectorND0V1'
def ordering(folder,value='ScotchVectorND0V1'):
    root=ET.parse(folder/'case.xml').getroot();set_info(root,'RANSWallPredictorOrdering',value);write(root,folder/'case.xml')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct9','ct10','ct11','runtime','capacity','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_OBSERVER']='NativeAuditReceiptV1';os.environ.pop('NEKTAR_RANS_TIMING_FILE',None);n.reference=lru_cache(None)(n.reference)
    cap=sealed(a.capacity);require(cap['status']=='pass','frozen capacity gate');controls=sealed(a.ct9/'controls-v4');oldpair=sealed(a.ct11/'target-v1')
    m=RequestManifest(a.output,dict(schema='CT14ComponentQualificationV1',source=sha(Path(__file__)),capacity=sha(a.capacity/'results.json'),contract=sha(a.contract),old_controls=sha(a.ct9/'controls-v4/results.json'),old_original_pair=sha(a.ct11/'target-v1/results.json'),runtime=runtime_identity(a.runtime/'IncNavierStokesSolver'),probe=runtime_identity(a.runtime/'RANSPreservationProbe'),component=COMPONENT))
    result=dict(schema='CT14ComponentQualificationV1',status='incomplete',small={},histories={},original={})
    for label in ('d0-p3','d0-p4','displaced-p3','displaced-p4'):
        src=a.ct9/'controls-v4'/(label+'-sparse');old=n.read(src/'operator.json');c=n.read(src/'construction.json');req=n.read(src/'operator-request.json')
        require(sha(src/'operator.json')==controls['checks'][label]['artifacts'][1]['sha256'],'old complete operator witness changed')
        packets=[];executions=[];artifacts=[]
        for policy in ('COLAMD','ScotchVectorND0V1'):
            dest=a.output/(label+'-'+policy);dest.mkdir();shutil.copy2(src/'case.xml',dest/'case.xml');ordering(dest,policy)
            q,e=probe(a.runtime,dest,req,'operator');require(e['exit_code']==0,'component small operation:'+str(q.get('failure')))
            require((COMPONENT in q['observation']['identity'])==(policy=='ScotchVectorND0V1'),'component identity')
            packets.append(q);executions.append(e);artifacts.append(dict(file=str(dest/'operator.json'),sha256=sha(dest/'operator.json')))
        unchanged={key:equal_primitives(old[key],packets[0][key]) for key in ('open_frozen_systems','open_boundary','test_basis')}
        require(max(unchanged.values())<=2e-10,'rebuilt canonical operators differ from qualified original')
        numerical=compare_packets(packets[0],packets[1],c)
        result['small'][label]=dict(status='pass',unchanged_old_operators=unchanged,complete_independent_actions=numerical,executions=executions,artifacts=artifacts)
        write_json(a.output/'progress.json',result);print('CT14 small component',label,'pass',flush=True)
    # Genuine same-component histories use new checkpoints created by this runtime.
    c=n.read(a.ct9/'controls-v4/one/construction.json');template=a.ct9/'controls-v4/one/case.xml'
    def evolve(label,steps=2,parent=None,policy='ScotchVectorND0V1',reject=0,foreign_fit=False):
        dest=a.output/label;n.make(dest,a.runtime,c,1e-6,steps,parent,template);enable(dest);ordering(dest,policy)
        root=ET.parse(dest/'case.xml').getroot();set_parameter(root,'RANSWallPimpleHardAuditRejectStep',reject)
        if foreign_fit:set_info(root,'RANSOpenFitPolicy','ForeignFitV1')
        write(root,dest/'case.xml');e=run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log');path=dest/'case.evolution.json'
        q=n.read(path) if path.exists() else dict(failure=(dest/'solver.log').read_text()[-3000:]);return dest,q,e
    full,fr,fe=evolve('full');one,orr,oe=evolve('one',1)
    require(fe['exit_code']==oe['exit_code']==0 and not fr['failure'] and not orr['failure'],'small native updates')
    cont,cr,ce=evolve('continued',1,Path(orr['checkpoint']));rollback,rr,re=evolve('rollback',reject=2)
    require(ce['exit_code']==0 and same_history(fr,cr),'genuine restart/full history equality')
    require(re['exit_code']!=0 and 'deliberate' in rr['failure'] and same_history(orr,rr),'full rejected-step rollback')
    default,dr,de=evolve('default-full',policy='COLAMD');require(de['exit_code']==0 and not dr['failure'],'default update control')
    require(maximum(np.asarray(fr['state'])-dr['state'])<=1e-8,'small component complete state update')
    negatives={}
    previous=n.read(a.ct9/'controls-v4/one/case.evolution.json')
    for label,parent,policy,foreign in [('old-component',Path(previous['checkpoint']),'ScotchVectorND0V1',False),('new-to-old',Path(orr['checkpoint']),'COLAMD',False),('foreign-fit',Path(orr['checkpoint']),'ScotchVectorND0V1',True),('unknown-order',None,'UnknownOrderingV1',False)]:
        dest,q,e=evolve(label,0,parent,policy,foreign_fit=foreign);message=q.get('failure','')+(dest/'solver.log').read_text()[-3000:]
        require(e['exit_code']!=0 and ('ordering' if label=='unknown-order' else 'signature') in message.lower(),'incompatible history/policy rejection:'+label)
        negatives[label]=dict(status='pass',failure=message,execution=e)
    result['histories']=dict(status='pass',uninterrupted_state=exact_vector_hash(fr['state']),continued_state=exact_vector_hash(cr['state']),rollback_state=exact_vector_hash(rr['state']),rollback_failure=rr['failure'],default_state_maximum_difference=maximum(np.asarray(fr['state'])-dr['state']),negatives=negatives,reports={str(f/'case.evolution.json'):sha(f/'case.evolution.json') for f in (full,one,cont,rollback,default)},executions=[fe,oe,ce,re,de])
    write_json(a.output/'progress.json',result)
    # Original full-domain control is solvable by both components. Reuse the
    # exact old trajectory and canonical observation rows; never relabel it.
    oldrow=oldpair['members']['0.001'];oldfolder=Path(oldrow['directory']);old=n.read(oldfolder/'case.evolution.json');require(sha(oldfolder/'case.evolution.json')==oldrow['native_report_sha256'],'old original trajectory changed')
    c=oldpair['case'];dest=a.output/'original-full';n.make(dest,a.runtime,c,.001,2);enable(dest);ordering(dest)
    e=run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log');new=n.read(dest/'case.evolution.json');require(e['exit_code']==0 and not new['failure'] and new['time']==.002,'original trajectory equivalence control')
    require(new['accepted_states'][0]==old['accepted_states'][0] and COMPONENT in new['signature'],'original initialization/component')
    ini=n.read(a.ct10/'target-preflight-v1/target/readiness.json');weights=n.inventory(c)[2];count=len(weights);common=[]
    for t,x,y in zip(old['audit_times'],old['accepted_states'],new['accepted_states'],strict=True):
        fields=n.norms(n.observations(ini,x,count),n.observations(ini,y,count),weights,ini['observation']['grid_regions'][:count],c);tr=difference(profiles(ini,x,c),profiles(ini,y,c),c)
        field_max=max(v['scaled'] for v in fields.values());traction_max=max(tr['profiles'][k]['normalized_weighted_RMS'] for k in ('ground','interface'))
        require(max(field_max,traction_max)<=1e-8,'original component field/absolute-pressure/traction update equivalence')
        common.append(dict(time=t,fields=fields,traction=tr,maximum_scaled_field=field_max,maximum_scaled_traction=traction_max))
    result['original']=dict(status='pass',common_times=common,execution=e,old_report_sha256=oldrow['native_report_sha256'],new_report_sha256=sha(dest/'case.evolution.json'),initial_sha256=exact_vector_hash(new['accepted_states'][0]),old_physical_failure=oldpair['qualification'])
    result['status']='pass';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT14 component pass',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
