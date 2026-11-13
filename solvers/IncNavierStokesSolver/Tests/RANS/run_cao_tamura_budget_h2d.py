#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H2d bounded native setup/replay and affine spatial diagnostics. Linux only."""
import argparse,json,os,platform,time,resource
from pathlib import Path
import xml.etree.ElementTree as ET
from run_cao_tamura_budget_h2c import case as h2c_case
from cls_epawf_m4_cases import tensor_mesh
from run_cls_epawf_picard import THREADS,load,sha,set_info,write,require
from diagnose_cls_accepted_covariance import run_probe
from continue_cls_epawf_pimple import freeze,atomic

def case(folder,h1,pre,policy='PhysicalHeightInterfaceV2',flat=False,z=None,degree=3):
    xml=h2c_case(folder,h1,pre,.0065,'affine' if flat else 'hill',fixture=flat,seed=policy)
    t=ET.parse(xml).getroot()
    if flat and z is not None:tensor_mesh(t,z,degree,2,.8,1.,points=[16,16],xknots=[-.2,.2,.6])
    set_info(t,'RANSWallCompositeSeed',policy);write(t,xml)
    d=load(folder/'definition.json');d.update(version='CaoTamuraBudgetH2dV1',initializer_version=policy,physics='unchanged H2c primary calibration',xml_sha256=sha(xml))
    if z is not None:
        nx=2;ni=sum(v<=.0065 for v in z[1:]);no=len(z)-1-ni;pp=degree-2
        sizes=[(degree*nx+1)*(degree*(ni+no)+1)]*3+[(degree*nx+1)*(degree*no+1),(pp*nx+1)*(pp*ni+1),(pp*nx+1)*(pp*no+1),1]
        require(sum(sizes)<=3000,'H2d diagnostic exceeds retained backend size limit')
        d.update(sizes=sizes,slots=sum(sizes),mesh=dict(x=[-.2,.2,.6],z=z,degree=degree,regional_pressure_degree=pp,points=[16,16]),physical_steps=0)
    atomic(folder/'definition.json',d);atomic(folder/'h2d-definition.json',d)
    return xml

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('probe','h1','preflight','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--phase',choices=['replay','setup','flat','candidate'],required=True);p.add_argument('--candidate',type=Path);p.add_argument('--combined',action='store_true')
    a=p.parse_args();require(platform.system()=='Linux','native Linux only');os.environ.update(THREADS)
    a.output.mkdir(parents=True,exist_ok=False);probe,runtime=freeze(a.probe.resolve(),a.output/'runtime')
    pre=load(a.preflight);out=dict(status='running',phase=a.phase,runtime=runtime,preflight_sha256=sha(a.preflight),results={},physical_steps=0)
    atomic(a.output/'progress.json',out)
    jobs=[('65-hill',None,3)]
    if a.phase=='replay':jobs=[(k,None,3) for k in ('TangentialNeutralV2','PhysicalHeightInterfaceV2')]
    if a.phase=='flat':
        # Same interval and physical data, two tangential panels; no evolution.
        h=pre['heights'][-1];rapid=pre['nu']/(pre['cMu_wall']*pre['cLength_wall']*h['K0']**.5)-h['viscosity_intercept']
        z=[0,rapid/8,.00005,.0002,h['cutoff'],h['h_s'],.0065,.012,.025,.05,.1,.2,.4,.8]
        refined=sorted(z+[rapid/16,.0001,.0004,.0015,.004,.009,.018,.035,.075,.15,.3,.6])
        jobs=[('flat-combined-confirmation',sorted(z+[rapid/16,.0001,.0004,.0015,.004,.009,.018,.035]),4)] if a.combined else [('flat-baseline',z,3),('flat-radial-h',refined,3),('flat-p4-diagnostic',z,4)]
    try:
        for name,z,degree in jobs:
            policy=name if a.phase=='replay' else 'PhysicalHeightInterfaceV2'
            xml=case(a.output/name,a.h1,pre,policy,a.phase=='flat',z,degree)
            packet=dict(operation='h2c_seed' if a.phase=='replay' else 'h2d_'+a.phase,inner_stress_model='Linear',xml_sha256=sha(xml))
            if a.phase=='candidate':
                c=load(a.candidate);packet.update(state=c['state'],state_sha256=c['state_sha256'],base_sha256=c['base_sha256'])
            result,receipt=run_probe(probe,xml,packet,a.output/(name+'-probe'))
            receipt['child_peak_RSS_KiB_cumulative_max']=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
            out['results'][name]=dict(receipt=receipt,result_file=str(a.output/(name+'-probe/result.json')))
            # run_probe output convention is recorded explicitly by receipt.
            out['results'][name]['result_file']=receipt['output']
            atomic(a.output/'progress.json',out)
        out['status']='completed'
    except Exception as e:out.update(status='diagnostic_failed',failure=str(e));raise
    finally:atomic(a.output/'progress.json',out)
if __name__=='__main__':main()
