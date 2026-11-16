# SPDX-License-Identifier: MIT
"""Common full-height neutral cases; model calibration is never refitted."""
from pathlib import Path
import hashlib,json,math
import xml.etree.ElementTree as ET
from cao_tamura_budget_h2g import case as historical_case, canonical_bytes
from cao_tamura_h2d_state import load
from cls_epawf_m4_cases import tensor_mesh
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha,set_info,set_parameter,write


def radial_mesh(pre, level):
    base=pre['common_mesh']['z'][:]
    if level==0:return base
    first=base[1]
    z=sorted(set(base+[first/2,.0001,.0004,.0015,.004,.009,.018,.035,.075,.15,.3,.6]))
    if level==1:return z
    for _ in range(min(level-1,1)):
        added=[]
        for lo,hi in zip(z,z[1:]):
            if hi<=.0002 or (lo>0 and hi/lo>(1.4 if hi<=.05 else 1.8)):
                added.append(math.sqrt(lo*hi) if lo else hi/2)
        z=sorted(set(z+added))
    if level==3:
        # The measured level2 U residual has reached an outer-panel floor;
        # further submicron inner splitting is not a useful next experiment.
        z=sorted(set(z+[(lo+hi)/2 for lo,hi in zip(z,z[1:]) if lo>=.05]))
    if level>=4:
        # Signed globally assembled mass-dual projections locate the remaining
        # U defect around H. Retain level2, without the ineffective high-outer
        # splits, and bisect the last inner / first two outer panels.
        z=sorted(set(z+[(lo+hi)/2 for lo,hi in zip(z,z[1:]) if lo>=.005 and hi<=.012]))
    if level==5:
        # A predeclared finer confirmation of the same interface allocation,
        # plus the largest remaining inner derivative-error panels.
        z=sorted(set(z+[(lo+hi)/2 for lo,hi in zip(z,z[1:])
                        if (lo>=.005 and hi<=.012) or (lo>=.0002 and hi<=.0015)]))
    return z


def case(folder,h1,pre,model,level=0,order=16,steps=0,dt=.001,restart=None,pout=0.,projection='RadialL2V1'):
    require(model in ('CLS','standard') and 0<=level<=5 and order in (16,24,32),'WP2 declared case space')
    mesh=dict(pre['common_mesh']);mesh.update(z=radial_mesh(pre,level),degree=3,regional_pressure_degree=1,points=[order,order])
    xml=historical_case(folder,h1,pre,flat=True,model=model,shared_mesh=mesh)
    t=ET.parse(xml).getroot();tensor_mesh(t,mesh['z'],3,2,.8,1.,points=[order,order],xknots=mesh['x'])
    manifest=dict(pre['model_manifest_standard' if model=='standard' else 'model_manifest_CLS'])
    manifest.update(boundary_policy='NeutralEquilibriumOpenV1',upper_boundary='W essential; driving full traction and prescribed epsilon conormal flux; no inlet tail')
    mf=folder/'model-manifest.json';mf.write_bytes(canonical_bytes(manifest)+b'\n')
    require(projection in ('RadialL2V1','RadialEnergyV1'),'unknown representation policy')
    for k,v in dict(RANSWallBoundaryPolicy='NeutralEquilibriumOpenV1',RANSWallGeometry='LegacyV1',RANSWallMap='Flat',RANSWallNeutralProjection=projection,
                    RANSWallModelManifestFile=str(mf),RANSWallModelManifestSHA256=sha(mf)).items():set_info(t,k,v)
    set_info(t,'RANSWallReportingZ',','.join(format(x,'.17g') for x in radial_mesh(pre,5)))
    if level>0:
        # Explicit numerical policy, supported by measured 5339/8002/9882-slot
        # full-inventory diagnostics (peak 3.72GiB). First-step cost must still
        # qualify affordable evolution; default historical cases keep 3000.
        set_info(t,'RANSWallResourcePolicy','WP2MeasuredDenseV1')
        set_parameter(t,'RANSWallDenseSizeLimit',12000)
    for k,v in dict(NumSteps=steps,TimeStep=dt,RANSWallEta=0,RANSWallOutletTotalPressure=pout,
                    RANSWallNormalOrder=order,RANSWallTangentialOrder=order).items():set_parameter(t,k,v)
    if restart:set_info(t,'RANSWallRestart',str(restart))
    write(t,xml)
    nx=len(mesh['x'])-1;ni=sum(z<=pre['H'] for z in mesh['z'][1:]);no=len(mesh['z'])-1-ni
    sizes=[(3*nx+1)*(3*(ni+no)+1)]*3+[(3*nx+1)*(3*no+1),(nx+1)*(ni+1),(nx+1)*(no+1),1]
    atomic(folder/'wp2-definition.json',dict(version='WP2CommonNeutralCaseV1',model=model,level=level,
        mesh=mesh,estimated_sizes=sizes,estimated_slots=sum(sizes),native_dimensions_required=True,
        common_mesh_sha256=hashlib.sha256(canonical_bytes(mesh)).hexdigest(),model_manifest_sha256=sha(mf),
        physical_calibration='unchanged H2g full-precision model-specific roots and deployed CLS roughness',
        boundary=manifest['boundary_policy'],Pout=pout,dt=dt,steps=steps,eta=0,
        inlet='same immutable analytical log/AWF composite over 0<=z<=.8; no tail; native trace converges',
        pressure='pi=p/rho+2k/3; streamwise outer total-normal-pressure Pout; no pressure Dirichlet or closed gauge',
        projected_initial_fields='U,nonconstant inner k,complete outer epsilon and model-correct regional pi; W=0',
        initial_projection=projection,old_hill_certificates_inapplicable=True,xml_sha256=sha(xml)))
    return xml
