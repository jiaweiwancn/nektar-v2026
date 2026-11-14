#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Recompute saved CT9 operator controls after a collector key correction.

No native operation or trajectory is rerun. The original failed receipt stays
sealed; complete native packets and all history identities are bound first.
"""
import argparse,os
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted9_controls import compare_packets
from ct_accepted9_pressure_run import sealed
from ct_accepted7_top import classify
import ct_accepted8_native as n
from ct_accepted_open import require
from compare_cls_epawf_pimple import same_history
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','ct8','source','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');old=sealed(a.source);request=n.read(a.source/'request-manifest.json')['request'];require(request['runtime']==runtime_identity(a.runtime/'IncNavierStokesSolver') and request['probe']==runtime_identity(a.runtime/'RANSPreservationProbe'),'changed source native runtime')
    files={str(f.relative_to(a.source)):sha(f) for f in a.source.rglob('*') if f.is_file() and (f.name in ('operator.json','operator-request.json','operator-execution.json','case.xml','case.evolution.json') or f.suffix=='.json' and 'checkpoint' in str(f))}
    m=RequestManifest(a.output,dict(schema='CT9BackendControlsReplayV2',source_receipt=sha(a.source/'results.json'),source=sha(Path(__file__)),control_implementation=sha(Path(__file__).with_name('ct_accepted9_controls.py')),runtime=request['runtime'],native_artifacts=files,correction='advection_load is the actual native key; previous KeyError retained, no changed operator, thresholds, fields or histories'))
    result=dict(schema='CT9BackendControlsV2',checks={},status='incomplete')
    for zero,degree in ((True,3),(True,4),(False,3),(False,4)):
        label=('d0' if zero else 'displaced')+'-p'+str(degree);packets=[];artifacts=[]
        try:
            for sparse in (False,True):
                folder=a.source/(label+('-sparse' if sparse else '-dense'));q=n.read(folder/'operator.json');c=n.read(folder/'construction.json');execution=n.read(folder/'operator-execution.json');require(execution['exit_code']==0,'source native operator failed')
                require(execution['actual_mapped_files'][str(a.runtime/'RANSPreservationProbe')]==sha(a.runtime/'RANSPreservationProbe'),'changed actually consumed native executable');require(('AcceptedSerialSparseV1' in q['observation']['identity'])==sparse,'foreign backend binding')
                independent=classify(q,c,**({} if zero else dict(column=n.col(c),reference_function=n.reference)));packets.append(q);artifacts.append(dict(path=str(folder/'operator.json'),sha256=sha(folder/'operator.json'),independent=independent))
            result['checks'][label]=dict(compare_packets(*packets,c),artifacts=artifacts)
        except Exception as error:result['checks'][label]=dict(status='fail',failure=repr(error))
        write_json(a.output/'progress.json',result);print('CT9 replay',label,result['checks'][label]['status'],flush=True)
    history=old['checks']['restart_rollback_backend_identity'];require(history['status']=='pass','historical native history control missing')
    for folder,digest in history['native_reports'].items():require(sha(Path(folder)/'case.evolution.json')==digest,'changed native history witness')
    full,continued,one,rejected,dense=[n.read(a.source/f/'case.evolution.json') for f in ('full','continued','one','rollback','dense-full')]
    require(same_history(full,continued) and same_history(one,rejected) and same_history(dense,n.read(a.ct8/'controls-v2/full/case.evolution.json')),'actual history/rollback/default-dense equality failed')
    result['checks']['restart_rollback_backend_identity']=history
    result['status']='pass' if len(result['checks'])==5 and all(v['status']=='pass' for v in result['checks'].values()) else 'fail';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
