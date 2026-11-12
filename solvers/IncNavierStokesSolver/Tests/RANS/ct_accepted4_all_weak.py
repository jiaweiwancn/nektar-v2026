#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only independent continuous weak balances for all nine CT4 references."""
import argparse,json,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted4_analysis import weak
from ct_accepted4_jobs import CASES,label
from ct_accepted4_reference import write_json
from ct_accepted_pilot import replay
from ct_accepted_open import probe,require
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','initial-root','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    cases=[x for x in CASES if x[2]==3]
    m=RequestManifest(a.output,dict(schema='CT4AllReferenceWeakV1',runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),
        sources={f:sha(Path(__file__).with_name(f)) for f in ('ct_accepted4_all_weak.py','ct_accepted4_analysis.py','ct_accepted4_reference.py')},
        initial_receipts={label(*x):sha(a.initial_root/label(*x)/'results.json') for x in cases}))
    result=dict(status='fail',cases={},scope='independent continuous equations on actual native test rows, not the residual of projected/evolved coefficients; no field/source changes')
    for model,h,degree in cases:
        key=label(model,h,degree);folder=a.initial_root/key/'0.002';initial=json.loads((folder/'initial.json').read_text())
        c=json.loads((folder/'construction.json').read_text());dest=a.output/key;dest.mkdir();shutil.copy2(folder/'case.xml',dest/'case.xml')
        root=ET.parse(dest/'case.xml').getroot()
        for k,v in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),(64,16,16,16,16),strict=True):set_parameter(root,k,v)
        write(root,dest/'case.xml')
        packet,e=probe(a.runtime,dest,replay(initial,initial['observation']['state'],0.,[],export_basis=True,export_resolution=True),'weak')
        require(e['exit_code']==0,'native_weak_packet_failed:'+key)
        r=weak(packet,c);r.update(execution=e,packet_sha256=sha(dest/'weak.json'),initial_packet_sha256=sha(folder/'initial.json'))
        result['cases'][key]=r;write_json(a.output/'progress.json',result)
    result['status']='pass' if len(result['cases'])==9 and all(q['status']=='pass' for q in result['cases'].values()) else 'fail'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT4 all-reference native weak',result['status']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
