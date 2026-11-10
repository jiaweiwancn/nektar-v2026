#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Production forward probe: declared sides override the interface tolerance."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from rans_qualification import THREADS, write_json
os.environ.update(THREADS)
from rans_request_manifest import sha, runtime_identity
from rans_preservation_oracle import Reference
from run_a93_a95_e94e95_matrix import case_xml


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--probe',type=Path,required=True);p.add_argument('--descriptor-probe',type=Path,required=True)
    p.add_argument('--legacy-probe',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();out=a.output or Path(tempfile.mkdtemp(prefix='rans-forward-sides-',dir=os.environ.get('TMPDIR')))
    if a.output:out.mkdir(parents=True,exist_ok=False)
    print('output='+str(out),flush=True)
    rows=[]
    def run(folder,probe,request):
        write_json(folder/'request.json',request)
        with (folder/'probe.log').open('x') as f:
            code=subprocess.run([str(probe),'case.xml','request.json','packet.json'],cwd=folder,
                stdout=f,stderr=subprocess.STDOUT,timeout=14400).returncode
        data=json.loads((folder/'packet.json').read_text()) if (folder/'packet.json').exists() else None
        return code,data
    for model in ('A93','A95','E94E95'):
        for wall in ('EP','Accepted'):
            folder=out/(model+'-'+wall);ref=Reference(model,'Natural',wall)
            case_xml(folder,model,'Natural',wall,steps=0,degree=3,radial_knots=ref.knots(),
                     descriptor_probe=a.descriptor_probe,boundary_policy='FrozenReferenceV1')
            near=[[0.,0.,ref.H+sign*5e-13] for sign in (-1,1) for side in (-1,1)]
            request=dict(operation='neutral',grid=near,grid_sides=[-1,1,-1,1],
                         export_constitutive=True,export_test_basis=True)
            code,data=run(folder,a.probe,request)
            assert code==0 and data['status']=='pass' and data['performed_time_steps']==0
            obs=data['observation'];samples=data['native_constitutive']['samples']
            assert obs['grid_regions']==[0,1,0,1] and obs['grid_sides']==request['grid_sides']
            assert [s['sample_index'] for s in samples]==[1,3]
            for j,row in enumerate(data['test_basis']):
                geometry=obs['grid_geometry_by_field'][j]
                assert row['velocity_k']['geometry']==geometry[0]==geometry[2]
                assert row['pressure']['geometry']==geometry[4]
                assert (row['epsilon'] is None) == (j in (0,2))
                if row['epsilon'] is not None:assert row['epsilon']['geometry']==geometry[3]==geometry[0]
            trace=max(max(abs(x-y) for x,y in zip(obs['grid_values'][0][:3],v[:3],strict=True)) for v in obs['grid_values'][1:])
            assert trace<1e-10
            rows.append(dict(model=model,wall=wall,status='pass',grid_regions=obs['grid_regions'],
                tensor_sample_indices=[s['sample_index'] for s in samples],maximum_common_velocity_k_trace_difference=trace,
                one_sided_velocity_gradients=[obs['grid_gradients'][j][0][2] for j in (0,1)],
                packet_sha256=sha(folder/'packet.json'),configuration_sha256=sha(folder/'case.xml')))
            if a.legacy_probe:
                ordinary=dict(request,grid=[[0.,0.,1.5],[0.,0.,ref.H],[0.,0.,ref.H],[0.,0.,2.5]],grid_sides=[0,-1,1,0])
                packets=[]
                for name,probe in (('legacy',a.legacy_probe),('corrected',a.probe)):
                    target=folder/name;target.mkdir();shutil.copy2(folder/'case.xml',target/'case.xml')
                    code,body=run(target,probe,ordinary);assert code==0 and body['status']=='pass';packets.append(body)
                fields=('state','grid_values','grid_gradients','grid_traction','grid_regions','grid_sides','grid_geometry_by_field')
                equal=all(packets[0]['observation'][key]==packets[1]['observation'][key] for key in fields)
                equal &= packets[0]['native_constitutive']==packets[1]['native_constitutive']
                assert equal;rows[-1]['ordinary_samples_bitwise_equal']=True
                if model=='A95' and wall=='Accepted':
                    target=folder/'legacy-near';target.mkdir();shutil.copy2(folder/'case.xml',target/'case.xml')
                    code,body=run(target,a.legacy_probe,request)
                    assert code!=0 and 'physical sample not found on requested radial side' in (target/'probe.log').read_text()
                    rows[-1]['before_reproduction']=dict(exit_code=code,reason='physical sample not found on requested radial side',
                        ordinary_physical_state_valid=True,scope='diagnostic selector failure, not a constitutive-domain failure')
    result=dict(status='pass',checks=rows,probe=runtime_identity(a.probe),legacy_probe=runtime_identity(a.legacy_probe) if a.legacy_probe else None,
                scope='six natural model/wall initialization audits only; no timestep or governing solve, exact/interior samples unchanged')
    write_json(out/'results.json',result);print(str(out));return 0


if __name__=='__main__':raise SystemExit(main())
