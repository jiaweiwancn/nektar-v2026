#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Small deterministic production-audit fixtures. No time advancement."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from rans_qualification import THREADS,dependencies,write_json
dependencies('numpy','scipy')
from run_a93_a95_e94e95_matrix import case_xml
from cls_epawf_m4_cases import tensor_mesh
from run_cls_epawf_picard import set_info,set_parameter,write,sha

def fixture(folder,label,descriptor):
    wall='Accepted' if label=='accepted_invalid_cutoff' else 'EP'
    model='A95' if label in ('inactive_invalid_constitutive','near_boundary_A95_forward') else 'E94E95'
    xml=case_xml(folder,model,'Natural',wall,steps=0,descriptor_probe=descriptor)
    root=ET.parse(xml).getroot();tensor_mesh(root,[0.,1.,2.,3.,4.],3,2,3.,1.,points=[8,8],xknots=[-4.,-3.,-2.25,-2.,-1.])
    set_info(root,'RANSWallReference','Functions')
    set_info(root,'RANSWallTransition','Fixed' if wall=='Accepted' else 'UniformEPFlat' if label=='active_invalid_cutoff' else 'Upstream')
    set_parameter(root,'RANSWallRoughness',0.)
    set_parameter(root,'RANSWallChi',1. if label in ('outside_spatial_support','active_invalid_cutoff') else 0.)
    K=1e-12;cl=2.55;nu=1.5e-5;epsilon=K*K/(2*cl*cl*nu)
    for function in list(root.findall('./CONDITIONS/FUNCTION')):root.find('CONDITIONS').remove(function)
    for name in ('InitialConditionsInner','InitialConditionsOuter','BoundaryValuesInner','BoundaryValuesOuter'):
        fn=ET.SubElement(root.find('CONDITIONS'),'FUNCTION',{'NAME':name})
        for field,value in [('u','x*y/2' if label=='inactive_invalid_constitutive' else '0'),('v','0'),
                            ('k',str(K)),('eps_tilde',str(-epsilon if label=='inactive_invalid_scalar' else epsilon)),('p','0')]:
            ET.SubElement(fn,'E',{'VAR':field,'VALUE':value})
    if label=='near_boundary_A95_forward':
        H=2.;a=6.5/(math.sqrt(18)-3)*(1-1e-13);K=1e-8;sink=K**1.5/(cl*H)
        tensor_mesh(root,[0.,1.,2.,3.,4.],4,2,3.,1.,points=[8,8],xknots=[-4.,-3.,-2.25,-2.,-1.])
        r='((y-2)/2)';shape='(3*'+r+'^2-2*'+r+'^3)'
        for fn in root.findall('./CONDITIONS/FUNCTION'):
            inner=fn.get('NAME').endswith('Inner')
            for e in fn:
                field=e.get('VAR')
                value={'u':'0' if inner else str(-a)+'*x*'+shape,
                    'v':'0' if inner else str(2*a)+'*('+r+'^3-0.5*'+r+'^4)',
                    'k':str(K),'eps_tilde':str(sink) if inner else str(sink)+'+'+str(K-sink)+'*'+shape,'p':'0'}[field]
                e.set('VALUE',value)
    write(root,xml);return xml

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--solver',type=Path,required=True);p.add_argument('--descriptor-probe',type=Path,required=True)
    p.add_argument('--output',type=Path);p.add_argument('--before',action='store_true');a=p.parse_args()
    temporary=None
    if a.output is None:temporary=tempfile.TemporaryDirectory(prefix='native-rans-audits-');a.output=Path(temporary.name)
    else:a.output.mkdir(parents=True,exist_ok=False)
    result=dict(schema=1,checks=[],status='running',scope='initialization and actual native audit only; no steps')
    for label in ('inactive_chi_zero','outside_spatial_support','active_invalid_cutoff','accepted_invalid_cutoff','inactive_invalid_scalar','inactive_invalid_constitutive','near_boundary_A95_forward'):
        folder=a.output/label;fixture(folder,label,a.descriptor_probe)
        with (folder/'solver.log').open('w') as f:
            completed=subprocess.run([str(a.solver),'case.xml'],cwd=folder,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,**THREADS),timeout=180)
        log=(folder/'solver.log').read_text();report=folder/'case.evolution.json';r=json.loads(report.read_text()) if report.exists() else {}
        baseline=label in ('inactive_chi_zero','outside_spatial_support')
        diagnostic=label=='near_boundary_A95_forward'
        expected_reason=('covariance_or_budget' if a.before and baseline else
            'radicand' if label in ('inactive_invalid_constitutive','near_boundary_A95_forward') else
            'epsilon' if label=='inactive_invalid_scalar' else 'inadmissib')
        if baseline and not a.before:
            success=completed.returncode==0 and not r.get('failure') and r.get('steps')==0
            if success:
                audit=r['state_audits'][0];success=audit['wall_cutoff_activity'][0]==0 and audit['wall_cutoff_required_budget'][0]==1 and audit['accepted_awf_scalar_domain'][7]<1
        elif diagnostic and not a.before:
            audit=(r.get('state_audits') or [{}])[0]
            success=completed.returncode==0 and not r.get('failure') and audit.get('accepted_awf_hard_admissible')==[1] and audit.get('outer_tangent_diagnostic_available')==[0] and audit.get('native_A95_radicand_range',[-1])[0]>0
        else:
            success=completed.returncode!=0 and expected_reason.lower() in log.lower()
            if a.before and label in ('inactive_invalid_scalar','inactive_invalid_constitutive'):
                success=completed.returncode!=0 and 'inadmissible evolution volume/sampler' in log
                expected_reason='legacy collective physical-state rejection (specific cause masked)'
        result['checks'].append(dict(name=label,status='pass' if success else 'fail',execution_exit=completed.returncode,
            intended_outcome='before_defect_rejection' if a.before and baseline else 'baseline_defined' if baseline else 'forward_physical_state_valid' if diagnostic and not a.before else 'expected_rejection',
            expected_reason=None if (baseline or diagnostic) and not a.before else expected_reason,xml_sha256=sha(folder/'case.xml'),log=str(folder/'solver.log'),failure=r.get('failure')))
        write_json(a.output/'audits.json',result)
    result['status']='pass' if len(result['checks'])==7 and all(x['status']=='pass' for x in result['checks']) else 'fail';write_json(a.output/'audits.json',result)
    print(json.dumps(result));return 0 if result['status']=='pass' else 1
if __name__=='__main__':sys.exit(main())
