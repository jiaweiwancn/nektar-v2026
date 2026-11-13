#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Discriminating CT6 read-only mutations and unsupported-route controls."""
import argparse,copy,json,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json
os.environ.update(THREADS)
import numpy as np
from ct_accepted6_cases import make
from ct_accepted6_operator import reconstruct
from ct_accepted_open import probe,require
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','operator','targeted','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();m=RequestManifest(a.output,dict(schema='CT6DiscriminatingControlsV1',
        runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),operator_sha256=sha(a.operator/'results.json'),
        targeted_sha256=sha(a.targeted/'results.json'),source_sha256=sha(Path(__file__))))
    summary=json.loads((a.operator/'results.json').read_text());checks={}
    path=a.operator/'p4'/'native.json'
    # The original packet location is explicit in its immutable result.
    if not path.exists():
        paths=list(a.operator.rglob('native.json'));require(paths,'missing sealed native operator packet')
        path=next((q for q in paths if '4' in str(q.parent.name)),paths[0])
    original=json.loads(path.read_text())
    base=reconstruct(original,'A95');require(base['status']=='pass','unmutated native all-row control failed')
    for name in ('omit_metric','wrong_outlet_orientation'):
        packet=copy.deepcopy(original)
        if name=='omit_metric':
            for row in packet['observation']['volume']:row['weight']/=row['volume_metric']
        else:
            for row in packet['observation']['faces']:
                if row['coordinate'][0]==1.2:row['area']=[-x for x in row['area']]
        value=reconstruct(packet,'A95')
        checks[name]=dict(status='pass' if value['status']=='fail' else 'fail',intended_outcome='independent all-row reconstruction must fail',observed=value)
    targeted=json.loads((a.targeted/'results.json').read_text())
    for model in ('A95','A93','E94E95'):
        row=targeted['checks'][model+'-static']['packet']
        for kind,key in [('omit_tilt','omitted_tilt_signal'),('fixed_chart_derivative','fixed_chart_derivative_signal')]:
            v=row[key];checks[model+'-'+kind]=dict(status='pass' if v>1e-8 else 'fail',primitive_signal=v,minimum_discriminating_signal=1e-8)
    # Curvature and normal motion contributions in the actual exact Cartesian
    # gradient must be observable, rather than a parallel-shear substitute.
    signal=0.
    for row in original['observation']['faces']:
        if not row['interface']:continue
        data=row['CT6_continuous'];T=np.asarray(data['T']);G=np.asarray(data['gradients'])[:2]
        signal=max(signal,abs(float(T[2,2]+data['values'][4])))
    checks['omit_normal_motion_traction']=dict(status='pass' if signal>1e-8 else 'fail',normal_viscous_or_nonlinear_stress_signal=signal,minimum_discriminating_signal=1e-8)
    variants=[('ordinary_curved_production',{'RANSWallCurvedCapability':'Off','RANSWallManufactured':'Off',
                  'RANSWallBoundaryPolicy':'Closed','RANSWallReference':'Neutral'}, {},'curved'),
              ('curved_EP',{'RANSWallFormulation':'EPAWF'}, {},'Accepted'),
              ('variable_height_evolution',{}, {'RANSWallMapHeightAmplitude':.03},'height'),
              ('unqualified_cache',{'RANSWallEvaluationReuse':'ExactStepV1'}, {},'caching')]
    for name,infos,params,word in variants:
        folder=a.output/name;make(folder,a.runtime,model='A95',degree=3,small=True)
        root=ET.parse(folder/'case.xml').getroot()
        for k,v in infos.items():set_info(root,k,v)
        for k,v in params.items():set_parameter(root,k,v)
        write(root,folder/'case.xml')
        packet,e=probe(a.runtime,folder,dict(operation='functions'),'negative')
        reason=packet.get('failure','')
        checks[name]=dict(status='pass' if e['exit_code']!=0 and word.lower() in reason.lower() else 'fail',
            intended_reason=word,observed_reason=reason,execution=e,packet_sha256=sha(folder/'negative.json'))
        write_json(a.output/'progress.json',checks)
    result=dict(schema='CT6DiscriminatingControlsV1',status='pass' if all(x['status']=='pass' for x in checks.values()) else 'fail',
                checks=checks,scope='mutations affect diagnostic copies only; no modified load is passed to governing evolution')
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
