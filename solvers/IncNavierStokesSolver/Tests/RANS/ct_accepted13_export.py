#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Actual-history pre-factor exports, canonical snapshot agreement, one mesh."""
import argparse,copy,json,os,resource,shutil
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted_pilot import replay
from ct_accepted_open import probe,run,require
from ct_accepted10_diagnose import matrix
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def maximum(x):return float(np.max(np.abs(x),initial=0.))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct10','ct11','ct12','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_PHASES']='1'
    old=sealed(a.ct11/'target-v1');ini=n.read(a.ct10/'target-preflight-v1/target/readiness.json')
    require(old['qualification']['reasons']==['temporal:Pi'],'changed historical target failure')
    m=RequestManifest(a.output,dict(schema='CT13PreFactorExportQualificationV1',source=sha(Path(__file__)),contract=sha(a.contract),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),target=sha(a.ct11/'target-v1/results.json'),static=sha(a.ct12/'static-v1/results.json')))
    rows={};systems={}
    algebra=run([a.runtime/'CT9SparseAlgebraChecks'],a.output,'algebra.log');require(algebra['exit_code']==0,'unchanged sparse algebra controls')
    for label,dt,which in [('original-initial','0.001',0),('original-evolved','0.0005',-1)]:
        directory=Path(old['members'][dt]['directory']);report=n.read(directory/'case.evolution.json')
        paths=[Path(v['flow_snapshot']) for v in report['pimple_attempts'] if v.get('flow_snapshot')]
        s=n.read(paths[which]);dest=a.output/label;dest.mkdir();shutil.copy2(directory/'case.xml',dest/'case.xml')
        req=replay(ini,s['input'],s['time'],[],export_pre_factor_frozen=True,actual_history=s['history'],sparse_alpha=s['alpha'],frozen_request_kind='EvolvedStageV1',pre_factor_prefix=str(dest/'frozen'))
        packet,e=probe(a.runtime,dest,req,'prefactor');require(e['exit_code']==0,'native export:'+str(packet.get('failure')))
        f=packet['pre_factor_frozen'];require(f['rows']==s['rows'] and f['effective_input']==s['input'] and f['history']==s['history'],'actual stage/history/mask mismatch')
        defects={};exact={}
        for name,left,right in [('system',f['frozen_flow_matrix'],s['frozen_flow_matrix']),('mass',f['free_velocity_mass'],s['free_velocity_mass'])]+[(name,f[name],s['pressure_operator'][name]) for name in ('G','D','C')]:
            delta=matrix(left)-matrix(right);defects[name]=maximum(delta.data);exact[name]=left['sha256']==right['sha256']
        defects['residual']=maximum(np.asarray(f['frozen_flow_residual'])-s['frozen_flow_residual'])
        require(max(defects.values())<=2e-10,'pre-factor all-row/snapshot disagreement')
        require(f['factorizations']==0 and f['updates']==0,'pre-factor performed a solve/update')
        rows[label]=dict(status='pass',execution=e,source_snapshot=str(paths[which]),source_snapshot_sha256=sha(paths[which]),packet_sha256=sha(dest/'prefactor.json'),defects=defects,exact_CSR_bytes=exact,history_sha256=f['history_sha256'],state_sha256=f['input_sha256'],actual_history_differs_from_state=f['history']!=s['input'],alpha=f['alpha'],time=f['time'])
        systems[label]=dict(packet=str(dest/'prefactor.json'),packet_sha256=sha(dest/'prefactor.json'),A=f['A'],rhs=[-x for x in f['frozen_flow_residual'][:f['velocity_size']]])
        write_json(a.output/'progress.json',dict(controls=rows))
    # The already selected fresh zero-step candidate is not an evolved history.
    candidate=n.read(a.ct12/'static-v1/target/case.evolution.json');require(candidate['time']==0 and len(candidate['state'])==62533,'foreign selected candidate')
    for alpha in (1000.,2000.):
        label='candidate-alpha'+str(int(alpha));dest=a.output/label;dest.mkdir();shutil.copy2(a.ct12/'static-v1/target/case.xml',dest/'case.xml')
        req=dict(operation='functions',grid=[],export_pre_factor_frozen=True,actual_history=candidate['state'],sparse_alpha=alpha,frozen_request_kind='InitializedFrozenV1',pre_factor_prefix=str(dest/'frozen'))
        packet,e=probe(a.runtime,dest,req,'prefactor');require(e['exit_code']==0,'candidate no-solve export:'+str(packet.get('failure')))
        f=packet['pre_factor_frozen'];require(f['input']==candidate['state'] and f['history']==candidate['state'],'candidate fresh projection changed')
        require(f['factorizations']==0 and f['updates']==0,'candidate diagnostic performed solve/update')
        rows[label]=dict(status='pass',execution=e,packet_sha256=sha(dest/'prefactor.json'),actual_initialized_state_sha256=exact_vector_hash(f['input']),initialized_frozen_fixture=True,alpha=alpha,slots=len(f['input']),velocity_rows=f['velocity_size'],matrix_nonzeros=f['frozen_flow_matrix']['nnz'])
        systems[label]=dict(packet=str(dest/'prefactor.json'),packet_sha256=sha(dest/'prefactor.json'),A=f['A'],rhs=[-x for x in f['frozen_flow_residual'][:f['velocity_size']]])
        write_json(a.output/'progress.json',dict(controls=rows))
    negatives={};dest=a.output/'original-evolved';base=n.read(dest/'prefactor-request.json')
    for label,change,reason in [('missing_history',lambda q:q.pop('actual_history'),'actual history'),('wrong_kind',lambda q:q.update(frozen_request_kind='Unspecified'),'history role'),('short_history',lambda q:q.update(actual_history=[0.]),'pre-factor frozen request')]:
        q=copy.deepcopy(base);change(q);packet,e=probe(a.runtime,dest,q,label)
        require(e['exit_code']!=0 and reason in str(packet.get('failure')),'wrong negative rejection:'+label)
        negatives[label]=dict(status='pass',intended_reason=reason,failure=packet.get('failure'),execution=e)
    result=dict(schema='CT13PreFactorExportQualificationV1',status='pass',controls=rows,systems=systems,negatives=negatives,algebra=algebra,scope='no factorization,update,new history or physical continuation; original stage histories and initialized candidate fixtures explicitly separate')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT13 export pass',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
