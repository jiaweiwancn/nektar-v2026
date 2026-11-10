#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Initialization-first native spatial diagnosis on explicitly declared spaces."""
import argparse,gzip,hashlib,json,os,subprocess,sys,time
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
from rans_qualification import THREADS,write_json,exact_vector_hash
from rans_preservation_oracle import Reference,physical_grid,representation,weak_balance
from rans_request_manifest import runtime_identity
from qualify_rans_preservation import load,Q,nested,replay_request
from run_a93_a95_e94e95_matrix import case_xml
from rans_spatial_resolution import analyze,common_test

def run_probe(runtime,folder,request):
    write_json(folder/'request.json',request);t=time.monotonic()
    with (folder/'probe.log').open('w') as f:c=subprocess.run([str(runtime/'RANSPreservationProbe'),'case.xml','request.json','packet.json'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,timeout=1200)
    if not (folder/'packet.json').exists():raise RuntimeError('missing native packet; exit '+str(c.returncode))
    data=load(folder/'packet.json')
    if c.returncode or data.get('status')!='pass':raise RuntimeError('native diagnostic failed '+str({k:v for k,v in data.items() if k not in ('observation','resolution','test_basis')}))
    elapsed=time.monotonic()-t
    with (folder/'packet.json').open('rb') as src,gzip.open(folder/'packet.json.gz','wb',compresslevel=3) as dst:__import__('shutil').copyfileobj(src,dst)
    (folder/'packet.json').unlink()
    return data,elapsed

def specs(ref):
    base=ref.knots();noaux=[x for x in base if x!=ref.kernel.ze]
    pairs=list(zip(base[:-1],base[1:]));interface=[(.02,.04),(.04,.1),(.1,.25),(.25,.5),(.5,1.),(1.,2.),(2.,2.25),(2.5,3.)]
    split=lambda intervals:sorted(set(base+[(a+b)/2 for a,b in intervals]))
    return dict(BaseP5=dict(degree=5,knots=base),Interface8P5=dict(degree=5,knots=split(interface)),Nearwall8P5=dict(degree=5,knots=split(pairs[:8])),H1P5=dict(degree=5,knots=nested(base,1)),H1P4=dict(degree=4,knots=nested(base,1)),NoAuxP5=dict(degree=5,knots=noaux))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--model',default='E94E95',choices=('E94E95','A95','A93'));p.add_argument('--policy',default='Natural',choices=('Natural','Prescribed'));p.add_argument('--wall',required=True,choices=('EP','Accepted'));p.add_argument('--label',required=True);p.add_argument('--names',nargs='+');p.add_argument('--specifications',type=Path);p.add_argument('--quadrature-index',type=int,default=1,choices=(0,1,2));a=p.parse_args();r=a.root.resolve();out=r/'diagnosis'/a.label;out.mkdir(parents=True,exist_ok=False)
    os.environ.update(THREADS,LD_LIBRARY_PATH=str(r/'runtime/lib'),TMPDIR=str(r/'tmp'))
    ref=Reference(a.model,a.policy,a.wall);choices=load(a.specifications) if a.specifications else specs(ref);names=a.names or list(choices)
    active=[x for x in ref.knots() if x!=ref.kernel.ze];grid,weights=physical_grid(sorted(set(nested(ref.knots(),1)+[x for s in choices.values() for x in s['knots']])),10,8)
    intent=dict(schema='BalancedSpatialDiagnosisV1',model=a.model,policy=a.policy,wall=a.wall,spaces={n:choices[n] for n in names},common_test=dict(degree=3,pressure_degree=1,knots=active,scope='fixed native p3 tests; omit unused original epsilon breakpoint in this common space only'),quadrature=Q[a.quadrature_index],runtime={n:runtime_identity(r/'runtime'/n) for n in ('RANSPreservationProbe','RANSNewModelsProbe')},evaluation=dict(grid=grid,weights=weights),alpha=0,time=0)
    write_json(out/'immutable-request.json',intent)
    common=out/'common';case_xml(common,a.model,a.policy,a.wall,steps=0,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=active,degree=3,quadrature=Q[a.quadrature_index],boundary_policy='FrozenReferenceV1')
    cp,_=run_probe(r/'runtime',common,dict(operation='neutral',grid=[],export_basis=True,export_resolution=True));checks=[]
    for name in names:
        try:
            spec=choices[name];folder=out/name;case_xml(folder,a.model,a.policy,a.wall,steps=0,descriptor_probe=r/'runtime/RANSNewModelsProbe',radial_knots=spec['knots'],degree=spec['degree'],quadrature=Q[a.quadrature_index],boundary_policy='FrozenReferenceV1')
            packet,elapsed=run_probe(r/'runtime',folder,dict(operation='neutral',grid=grid,export_basis=True,export_resolution=True));o=packet['observation']
            cross=out/(name+'-common-cross');cross.mkdir();__import__('shutil').copy2(common/'case.xml',cross/'case.xml');request=replay_request(cp,[q['coordinate'] for q in o['volume']]);request['export_test_basis']=True
            crosspacket,_=run_probe(r/'runtime',cross,request)
            record=dict(name=name,status='pass',space=dict(spec,pressure_degree=spec['degree']-2,dof=sum(o['sizes']),block_sizes=o['sizes']),evolution_eligible=sum(o['sizes'])<=3000,
                exact_source_runtime=intent['runtime'],initial_state_sha256=exact_vector_hash(o['state']),reference_sha256=exact_vector_hash(packet['reference']),initialization=packet['initialization_receipt'],
                initialization_and_full_export_seconds=elapsed,actual_quadrature=o['actual_quadrature'],representation=representation(ref,o,weights),analytical_weak_balance=weak_balance(ref,o),
                spatial=analyze(packet,ref),common_test=common_test(cp,packet,crosspacket,ref),native_alpha_zero=o['audit_alpha_zero'],source_packet_sha256=hashlib.sha256((folder/'packet.json.gz').read_bytes()).hexdigest(),source=str(folder))
            if record['analytical_weak_balance']['status']!='pass':record['status']='fail'
            # Large coefficient/boundary arrays remain in the sealed native packet.
            record['initialization']={k:v for k,v in record['initialization'].items() if k not in ('effective_initial','live_fixed_boundary')}
            write_json(folder/'analysis.json',record);checks.append(record)
        except Exception as e:checks.append(dict(name=name,status='fail',failure=repr(e)))
        write_json(out/'results.json',dict(status='running',checks=checks));print(name,checks[-1]['status'],flush=True)
    status='pass' if len(checks)==len(names) and all(c['status']=='pass' for c in checks) else 'fail';write_json(out/'results.json',dict(status=status,checks=checks,scope='initialization-only screen, no evolved spatial qualification'));return 0 if status=='pass' else 1
if __name__=='__main__':sys.exit(main())
