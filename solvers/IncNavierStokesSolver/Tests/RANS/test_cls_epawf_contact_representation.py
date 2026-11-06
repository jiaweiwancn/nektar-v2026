#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""One fixed contact state: source quadrature and a single near-wall mesh split.

No physical trajectory or checkpoint is created. Mesh refinement restricts
the same represented polynomial and essential traces into the new space.
"""
import argparse
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cls_epawf_picard import load,run,write_json,write,set_parameter,require,sha
from cls_epawf_m4_cases import tensor_mesh
from cls_epawf_m8_r2_cases import RADIAL,X


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--prior',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output=a.output.resolve();a.output.mkdir(parents=True,exist_ok=False);a.prior=a.prior.resolve();a.probe=a.probe.resolve()
    xml=a.prior/'fixed-horizon/bound-0/case.xml';report=load(xml.with_suffix('.evolution.json'))
    row=next(v for v in report['pimple_attempts'] if abs(v['time']-.0025)<1e-13);packet=load(Path(row['k_snapshot']))
    state=packet['corrected'];summary=dict(status='running',state_snapshot=str(Path(row['k_snapshot'])),cases={})
    for name in ('baseline','normal32','near_split'):
        out=a.output/name;out.mkdir();root=ET.parse(xml).getroot()
        if name=='normal32':set_parameter(root,'RANSWallNormalOrder',32)
        if name=='near_split':tensor_mesh(root,[0.,.005,*RADIAL['outer6'][1:]],3,2,6.,1.,points=[8,8],xknots=X)
        write(root,out/'case.xml')
        request=dict(operation='contact',state=state,history=state,time=.0025,alpha=0.,require_history_h0=False,
            inner_stress_model='SaturatedR2',h=[1e-6,1e-8,1e-9],witnesses=[packet['k_bound']['corrected_ground']['point']])
        if name=='near_split':request['transfer_source_xml']=str(xml)
        write_json(out/'input.json',request)
        receipt=run([str(a.probe),'case.xml','input.json','result.json'],out,'probe.log',1800)
        require(receipt['exit_code']==0,'representation replay process failed: '+name)
        result=load(out/'result.json')
        summary['cases'][name]=dict(receipt=receipt,xml_sha256=sha(out/'case.xml'),result_sha256=sha(out/'result.json'),
            result={k:v for k,v in result.items() if k not in ('scan',)})
        if 'scan' in result:
            summary['cases'][name]['scan']=[dict(h=v['h'],points=v['points'],rate_error=v['raw_minus_full_spatial_rate_mass_norm'],
                bound_rate_mass_norm=v['bound_rate_mass_norm'],bound_success=v['bound']['success']) for v in result['scan']]
        write_json(a.output/'summary.partial.json',summary)
    summary['status']='recorded';write_json(a.output/'summary.json',summary);print(a.output/'summary.json')


if __name__=='__main__':main()
