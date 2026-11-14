#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Recompute saved CT9 operator controls after a collector key correction.

No existing native operation or trajectory is rerun. A new, explicitly shifted
absolute-pressure control is added. The original failed receipt stays sealed;
complete native packets and all history identities are bound first.
"""
import argparse,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS
os.environ.update(THREADS)
from ct_accepted9_controls import compare_packets
from ct_accepted9_controls import enable,maximum
from ct_accepted9_pressure_run import sealed
from ct_accepted7_top import classify
import ct_accepted8_native as n
from ct_accepted_open import require,run
from run_cls_epawf_picard import set_parameter,write
import numpy as np
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
    try:
        c=request['case'];c.update(degree=3,pressure_degree=1);dest=a.output/'pressure-offset';n.make(dest,a.runtime,c,1e-6,2);enable(dest);root=ET.parse(dest/'case.xml').getroot();set_parameter(root,'RANSOpenPressureLevel',.125);write(root,dest/'case.xml')
        execution=run([a.runtime/'IncNavierStokesSolver','case.xml'],dest,'solver.log');shifted=n.read(dest/'case.evolution.json');require(execution['exit_code']==0 and not shifted['failure'],'shifted sparse pressure-level control rejected')
        packet=n.read(a.source/'displaced-p3-sparse/operator.json');off=packet['observation']['offsets'];expected=np.zeros(len(full['state']))
        for q in packet['observation']['pressure_mass']:
            row=q['pressure'];expected[np.asarray(row['ids'],int)+off[4+q['region']]]=.125*np.asarray(row['constant'])
        error=maximum(np.asarray(shifted['state'])-full['state']-expected);require(error<2e-10,'absolute traction/pressure offset changed sparse physical fields')
        result['checks']['absolute_pressure_level']=dict(status='pass',shift=.125,maximum_coefficient_defect=error,native_report_sha256=sha(dest/'case.evolution.json'),execution=execution,scope='one new two-tiny-step pressure-level control; prior native histories are reused unchanged')
    except Exception as error:result['checks']['absolute_pressure_level']=dict(status='fail',failure=repr(error))
    result['status']='pass' if len(result['checks'])==6 and all(v['status']=='pass' for v in result['checks'].values()) else 'fail';write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);return 0 if result['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
