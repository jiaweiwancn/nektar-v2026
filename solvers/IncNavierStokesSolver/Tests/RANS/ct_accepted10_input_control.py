#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""V2 observer: unchanged initial loads and explicit governing-residual lift."""
import argparse,os,shutil
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted_open import probe,require
from ct_accepted_pilot import replay
from ct_accepted_feasibility import write_json
from rans_request_manifest import RequestManifest,sha,runtime_identity


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','ct9','runtime','old-runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');pre=sealed(a.root/'static-v1');end=sealed(a.root/'candidate-T020-v2')
    m=RequestManifest(a.output,dict(schema='CT10ActualInputControlV1',source=sha(Path(__file__)),new_runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),old_runtime=runtime_identity(a.old_runtime/'RANSPreservationProbe'),preflight=sha(a.root/'static-v1/results.json'),continuation=sha(a.root/'candidate-T020-v2/results.json'),scope='observer only; no governing trajectory or history modified'))
    out=dict(schema='CT10ActualInputControlV1',status='pass',checks={})
    for label,initial_path,folder,state,time in [
        ('baseline-initial',a.ct9/'equivalence-v2/zero-0.001/initial.json',a.ct9/'equivalence-v2/0.001',None,0.),
        ('candidate-initial',Path(pre['initial']['candidate-0.001']['directory'])/'initial.json',Path(end['members']['0.0005']['directory']),None,0.),
        ('candidate-evolved',Path(pre['initial']['candidate-0.001']['directory'])/'initial.json',Path(end['members']['0.0005']['directory']),'evolved',.02)]:
        ini=n.read(initial_path);native=n.read(folder/'case.evolution.json');s=native['state'] if state else ini['observation']['state'];dest=a.output/label;dest.mkdir();shutil.copy2(folder/'case.xml',dest/'case.xml')
        req=replay(ini,s,time,[],export_ct10_momentum=True,ct10_export_points=False);q,e=probe(a.runtime,dest,req,'actual');require(e['exit_code']==0 and q['observation']['state']==s,'new observer changed/rejected actual state')
        req=replay(ini,s,time,[],export_resolution=True,resolution_localize=False);old,oe=probe(a.old_runtime,dest,req,'original');require(oe['exit_code']==0 and old['observation']['state']==s,'original residual replay failed')
        mom=q['ct10_momentum'];require(mom['schema']=='CT10SignedMomentumV2','wrong observer input schema');v=np.array(old['resolution']['native']['spatial_unmasked'])[:mom['momentum_size']];diff=float(max(abs(np.asarray(mom['total_actual'])-v)))
        require(diff<=2e-10,'effective input residual disagreement')
        row=dict(status='pass',time=time,state_sha256=exact_vector_hash(s),maximum_original_residual_difference=diff,governing_matching_enabled=mom['governing_matching_enabled'],maximum_input_lift_difference=mom['maximum_input_lift_difference'],lifted_rows=mom['input_lift_rows'],lifted_delta=mom['input_lift_delta'],new_packet_sha256=sha(dest/'actual.json'),old_packet_sha256=sha(dest/'original.json'),execution=e,original_execution=oe,scope='governing-residual effective input; raw accepted coefficients and their traction/epsilon-lag diagnostics are unchanged')
        if label.endswith('initial'):
            old_label='baseline-attribution-v3' if label.startswith('baseline') else 'candidate-attribution-v1';receipt=sealed(a.root/old_label);captured=n.read(Path(receipt['packet_path']))['ct10_momentum'];maximum=max(float(max(abs(np.asarray(mom[k][name])-np.asarray(captured[k][name])))) for k in ('actual','analytic') for name in mom[k]);require(maximum<=2e-10 and mom['maximum_input_lift_difference']==0.,'sealed raw/effective initial attribution changed');row['maximum_sealed_initial_component_difference']=maximum
        out['checks'][label]=row;write_json(a.output/'progress.json',out)
    write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT10 actual-input observer controls pass',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
