#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT6 zero-step native initialization, projection and operator qualification."""
import argparse,json,os,subprocess,traceback
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
from ct_accepted6_cases import make,request
from ct_accepted6_metrics import metrics
from ct_accepted6_operator import reconstruct,continuous_weak,geometry_qualification
from ct_accepted_open import run,probe,require
from ct_accepted_pilot import artifact_hash
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('runtime','contract','mathematics','output'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    mathematical={}
    for name in ('continuous','kernels'):
        path=a.mathematics/name/'results.json';receipt=json.loads(path.read_text())
        require(receipt['status']=='pass','required_independent_'+name+'_gate')
        mathematical[name]=dict(path=str(path),sha256=sha(path))
    manifest=RequestManifest(a.output,dict(schema='CT6NativeInitializationV1',contract_sha256=sha(a.contract),
        runtime={name:runtime_identity(a.runtime/name) for name in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
        required_mathematics=mathematical,source_sha256={n:sha(Path(__file__).with_name(n)) for n in
            ('ct_accepted6_initial.py','ct_accepted6_metrics.py','ct_accepted6_operator.py','ct_accepted6_cases.py')},
        cases=[['A95',3],['A95',4],['A93',4],['E94E95',4]],time=0.,dt=.002,zero_step_history='fresh native IMEX1',
        counts=dict(P3=1273,P4=2283),pressure='absolute, no gauge fitting'))
    result=dict(schema='CT6NativeInitializationV1',cases={})
    for model,degree in [('A95',3),('A95',4),('A93',4),('E94E95',4)]:
        label=model+'-p'+str(degree);folder=a.output/label;row=dict(status='incomplete');result['cases'][label]=row
        try:
            make(folder,a.runtime,model=model,degree=degree)
            execution=run([a.runtime/'IncNavierStokesSolver','case.xml'],folder,'solver.log');row['execution']=execution
            require((folder/'case.evolution.json').exists(),'native initialization produced no report')
            native=json.loads((folder/'case.evolution.json').read_text());row['native_report_sha256']=sha(folder/'case.evolution.json')
            require(execution['exit_code']==0 and 'state' in native and not native.get('failure'),'native zero-step rejection: '+str(native.get('failure',native)))
            row.update(stored_slots=len(native['state']),state_sha256=exact_vector_hash(native['state']),
                checkpoint=native['checkpoint'],history_sha256=artifact_hash(native['checkpoint']),
                raw_native_alpha_zero=native['raw_reference_norms'],initialization_receipt=native['initialization_receipt'])
            require(row['stored_slots']==(1273 if degree==3 else 2283),'actual native count differs from contracted count')
            req=request(False);req['export_tensor_validation']=False
            packet,e=probe(a.runtime,folder,req,'initial');row['probe_execution']=e
            require(e['exit_code']==0,'native initial probe rejected: '+str(packet.get('failure')))
            require(packet['observation']['state']==native['state'],'driver and probe effective initialized coefficients differ')
            row.update(packet_sha256=sha(folder/'initial.json'),packet_path=str(folder/'initial.json'),
                physical_metrics=metrics(packet,model),all_row=reconstruct(packet,model),
                continuous_weak=continuous_weak(packet,model),geometry=geometry_qualification(packet),
                raw_native_audit=packet['observation']['audit_alpha_zero'],
                actual_quadrature=packet['observation']['actual_quadrature'])
            row['status']='pass' if all(row[k]['status']=='pass' for k in ('physical_metrics','all_row','continuous_weak','geometry')) else 'fail'
        except Exception as error:
            row.update(status='fail',failure=repr(error));traceback.print_exc()
        write_json(a.output/'progress.json',result);print(label,row['status'],flush=True)
    result['status']='pass' if all(row['status']=='pass' for row in result['cases'].values()) else 'fail'
    write_json(a.output/'results.json',result)
    manifest.seal(a.output,['results.json']+[str(p.relative_to(a.output)) for p in a.output.rglob('case.xml')])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
