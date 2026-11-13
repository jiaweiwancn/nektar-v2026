#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Frozen actual RN coefficients: explicitly selected external-face control."""
import argparse,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
from ct_accepted4_analysis import weak
from ct_accepted_pilot import replay
from ct_accepted_open import probe,require
from run_cls_epawf_picard import set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','input','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--face-orders',type=int,nargs=2,default=[24,32])
    a=p.parse_args();labels=['A95-p3','A93-p3','E94E95-p3','A95-p4']
    require(a.face_orders in ([24,32],[64,96]),'outside the bounded CT7 face controls')
    m=RequestManifest(a.output,dict(schema='CT7RNFrozenFaceQuadratureV1',runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),
        source=sha(Path(__file__)),packets={n:sha(a.input/n/'native.json') for n in labels},face_orders=a.face_orders,limit=2e-10))
    rows={}
    for name in labels:
        old=a.input/name;ini=json.loads((old/'native.json').read_text());c=json.loads((old/'construction.json').read_text());rows[name]={}
        for order in a.face_orders:
            folder=a.output/(name+'-f'+str(order));folder.mkdir();root=ET.parse(old/'case.xml').getroot();set_parameter(root,'RANSWallInterfaceOrder',order);write(root,folder/'case.xml')
            packet,e=probe(a.runtime,folder,replay(ini,ini['observation']['state'],0.,[],export_basis=True,export_resolution=True),'frozen')
            require(e['exit_code']==0,'frozen face quadrature native rejection')
            require(packet['observation']['state']==ini['observation']['state'] and packet['reference']==ini['reference'],'changed frozen actual coefficients/reference')
            keys=('explicit_boundary_sha256','implicit_boundary_sha256','fixed_reference_sha256')
            require(all(packet['initialization_receipt'][k]==ini['initialization_receipt'][k] for k in keys),'changed consumed frozen boundary/reference')
            q=weak(packet,c);maximum=max(r['scaled_free_maximum'] for r in q['equations'].values())
            rows[name][str(order)]=dict(analytical_weak=q,maximum_scaled_free_load=maximum,limit=2e-10,status='pass' if maximum<=2e-10 else 'fail',
                frozen_state_and_boundary_identical=True,packet_sha256=sha(folder/'frozen.json'),execution=e)
            write_json(a.output/'progress.json',rows);print('CT7 frozen top face',name,order,maximum,flush=True)
    result=dict(schema='CT7RNFrozenFaceQuadratureV1',cases=rows,status='pass' if all(r['status']=='pass' for v in rows.values() for r in v.values()) else 'fail',physical_admission=False)
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if result['status']=='pass' else 1

if __name__=='__main__':raise SystemExit(main())
