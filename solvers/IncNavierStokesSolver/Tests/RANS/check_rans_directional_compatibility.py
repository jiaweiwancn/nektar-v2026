#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Pinned old/new isotropic rows, mass, residual, signatures and field values."""
import argparse,os
from pathlib import Path
os.environ.update(dict.fromkeys(('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','NUMEXPR_NUM_THREADS'),'1'))
import numpy as np
from run_a93_a95_e94e95_matrix import case_xml
from qualify_rans_preservation import Q
from diagnose_rans_balanced_spatial import run_probe
from rans_preservation_oracle import Reference,physical_grid
from rans_duration_contract import require
from rans_request_manifest import RequestManifest,runtime_identity,sha
from rans_qualification import write_json


def equal(a,b,path=''):
    if isinstance(a,dict):
        require(isinstance(b,dict) and set(a)==set(b),'isotropic_keys:'+path)
        return max((equal(a[k],b[k],path+'/'+str(k)) for k in a),default=0.)
    if isinstance(a,list):
        require(isinstance(b,list) and len(a)==len(b),'isotropic_length:'+path)
        return max((equal(x,y,path+'/'+str(i)) for i,(x,y) in enumerate(zip(a,b,strict=True))),default=0.)
    if isinstance(a,(int,float)):
        require(isinstance(b,(int,float)) and np.isfinite(a) and np.isfinite(b),'isotropic_nonfinite:'+path)
        require(a==b,'isotropic_arithmetic_changed:'+path+':'+str(a-b))
        return abs(a-b)
    require(a==b,'isotropic_identity_changed:'+path);return 0.


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--old-runtime',type=Path,required=True);p.add_argument('--runtime',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    manifest=RequestManifest(a.output,dict(schema='DirectionalIsotropicCompatibilityV1',degrees=[3,5],
        old_runtime=runtime_identity(a.old_runtime/'RANSPreservationProbe'),runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),
        criterion='exact parsed numeric equality and unchanged native identity; all array lengths and keys checked'))
    ref=Reference('A95','Natural','Accepted');grid,_=physical_grid(ref.knots(),4,4)
    grid += [[-.37,0.,ref.H],[-.37,0.,ref.H]];sides=[0]*(len(grid)-2)+[-1,1];rows=[]
    try:
        for degree in (3,5):
            packets=[]
            for label,runtime in (('old',a.old_runtime),('new',a.runtime)):
                folder=a.output/(str(degree)+'-'+label)
                case_xml(folder,'A95','Natural','Accepted',steps=0,degree=degree,quadrature=Q[1],
                    descriptor_probe=runtime/'RANSNewModelsProbe',radial_knots=ref.knots(),boundary_policy='FrozenReferenceV1')
                os.environ['LD_LIBRARY_PATH']=str(runtime/'lib')
                packet,_=run_probe(runtime,folder,dict(operation='neutral',grid=grid,grid_sides=sides,export_basis=True,export_resolution=True))
                packets.append(packet)
            old,new=packets;checks={}
            for key in ('reference','fixed_boundary','initialization_receipt'):
                checks[key]=equal(old[key],new[key],key)
            for key in ('identity','state','offsets','sizes','profiles','volume','faces','pressure_mass','audit_alpha_zero',
                        'grid','grid_regions','grid_sides','grid_geometry_by_field','grid_traction','grid_values','grid_gradients'):
                checks[key]=equal(old['observation'][key],new['observation'][key],key)
            checks['resolution']=equal(old['resolution'],new['resolution'],'resolution')
            rows.append(dict(degree=degree,status='pass',exact_maximum_differences=checks,
                packet_sha256={label:sha(a.output/(str(degree)+'-'+label)/'packet.json.gz') for label in ('old','new')}))
        result=dict(status='pass',checks=rows,scope='frozen isotropic arithmetic and identity only; historical trajectories remain bound to their pinned binary')
    except Exception as e:result=dict(status='fail',failure=repr(e),checks=rows)
    write_json(a.output/'results.json',result)
    manifest.seal(a.output,[p for p in a.output.rglob('*') if p.is_file() and p.name not in ('request-manifest.json','completed-artifacts.json')])
    print(result['status']);return 0 if result['status']=='pass' else 1


if __name__=='__main__':raise SystemExit(main())
