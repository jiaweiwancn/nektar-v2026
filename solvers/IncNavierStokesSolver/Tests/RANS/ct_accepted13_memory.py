#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Locate/retest the exact evolved aggregate q64 allocation, separate from scalars."""
import argparse,os,resource,shutil
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('ct12','runtime','contract','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');os.environ['NEKTAR_RANS_PROBE_PHASES']='1'
    prior=sealed(a.ct12/'pressure-v5');m=RequestManifest(a.output,dict(schema='CT13EvolvedAggregateMemoryV1',source=sha(Path(__file__)),contract=sha(a.contract),prior=sha(a.ct12/'pressure-v5/results.json'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),address_space_GiB=16))
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
            row.update(all_row_defects=diffs,mandatory_diagnostics_exact=True,effective_input_policy=packet['ct10_momentum'].get('input_policy'))
        elif e['exit_code']!=0:
            require(packet.get('failure')=='std::bad_alloc','unexpected exact q64 failure')
        controls[str(order)]=row;write_json(a.output/'progress.json',dict(controls=controls));print('CT13 memory',order,row['status'],phases[-1:] if phases else [],flush=True)
    result=dict(schema='CT13EvolvedAggregateMemoryV1',status=controls['64']['status'],controls=controls,historical_CT12='q64 failure unchanged; four passing scalar observer controls are separate',physical_continuation=False)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0
if __name__=='__main__':raise SystemExit(main())
