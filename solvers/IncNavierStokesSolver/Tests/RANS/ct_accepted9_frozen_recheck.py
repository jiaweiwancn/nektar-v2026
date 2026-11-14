#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Same coefficient state, new numerical-component signature, read-only only."""
import argparse,os,resource
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import ct_accepted8_native as n
from ct_accepted_open import probe,require
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','parent','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');resource.setrlimit(resource.RLIMIT_CORE,(0,0));resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));native=n.read(a.parent/'case.evolution.json');c=n.read(a.parent/'construction.json');request=n.read(a.parent/'frozen-request.json')
    m=RequestManifest(a.output,dict(schema='CT9FrozenNumericalComponentRecheckV1',parent_report=sha(a.parent/'case.evolution.json'),parent_configuration=sha(a.parent/'case.xml'),state_sha256=exact_vector_hash(native['state']),old_signature=native['signature'],runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),request=request,source=sha(Path(__file__)),scope='fresh provider evaluation must reproduce the old coefficients exactly; no old checkpoint/signature is edited, resumed or supplied as new-component history'))
    folder=a.output/'work';n.make(folder,a.runtime,c,.001,0,template=a.parent/'case.xml');os.environ['NEKTAR_RANS_TIMING_FILE']=str(folder/'timing.jsonl');q,e=probe(a.runtime,folder,request,'frozen')
    result=dict(status='fail',execution=e,packet_sha256=sha(folder/'frozen.json'),old_signature=native['signature'],failure=q.get('failure'))
    if e['exit_code']==0:
        require(q['observation']['state']==native['state'],'recheck changed the represented physical state');result.update(status='pass',same_coefficient_bytes=exact_vector_hash(q['observation']['state'])==exact_vector_hash(native['state']),new_identity=q['observation']['identity'],sparse_backend=q['sparse_backend'])
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT9 frozen recheck',result['status'],result.get('failure'),flush=True);return e['exit_code'] or (0 if result['status']=='pass' else 1)
if __name__=='__main__':raise SystemExit(main())
