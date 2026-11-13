# SPDX-License-Identifier: MIT
"""H2g common mesh and explicit model/calibration cases; no hidden recalibration."""
from pathlib import Path
import copy,hashlib,json
import xml.etree.ElementTree as ET
from cao_tamura_h2d_state import load
from run_cao_tamura_budget_h2d import case as previous_case
from cls_epawf_m4_cases import tensor_mesh
from run_cls_epawf_picard import set_info,set_parameter,write,sha,require
from continue_cls_epawf_pimple import atomic

def canonical_bytes(d):return json.dumps(d,sort_keys=True,separators=(',',':')).encode()

def case(folder,h1,pre,flat=False,model='standard',manifest=None,shared_mesh=None):
    # Pass the retained CLS calibration solely to obtain the historical geometry.
    # Model/reference/scalar data are subsequently changed together below.
    baseline=pre['cls_preflight'];xml=previous_case(folder,h1,baseline,flat=flat)
    t=ET.parse(xml).getroot();mesh=pre['common_mesh'] if shared_mesh is None else shared_mesh
    tensor_mesh(t,mesh['z'],3,2,.8,1.,points=[16,16],xknots=mesh['x'])
    standard=model=='standard';cal=dict(pre['models'][1 if standard else 0]);
    if not standard:cal['h_s']=.0022472751006336437
    spec=pre['model_manifest_standard'] if standard else pre['model_manifest_CLS']
    manifest_file=folder/'model-manifest.json';manifest_file.write_bytes(canonical_bytes(spec)+b'\n');manifest_hash=sha(manifest_file)
    if standard:
        for k,v in dict(TurbulenceModel='StandardKEpsilon' if standard else 'CLSHighRe',RANSWallModelCalibration='StandardNeutralCalibratedV1' if standard else '',RANSWallFormulation='StandardCalibratedBudgetAWF' if standard else 'CLSCalibratedBudgetAWF',RANSWallVariant='Original' if standard else 'Regularized',RANSWallDestruction='Original',RANSWallModelManifestSHA256=manifest_hash,RANSWallModelManifestFile=str(manifest_file)).items():set_info(t,k,v)
        for k,v in dict(RANSWallRoughness=cal['h_s'],RANSWallCubicBeta=1 if standard else .4,RANSWallRegularizationScale=0 if standard else 1,RANSWallRegularizationGamma=0,Qstar=cal['q'],Kappa=cal['kappa'],K0=cal['K0'],RANSCmu=.09,RANSCepsilon1=1.44,RANSCepsilon2=1.92,RANSSigmaK=1.,RANSSigmaEpsilon=cal['sigmaEpsilon'],RANSYapCoefficient=0 if standard else .83).items():set_parameter(t,k,v)
    write(t,xml)
    d=load(folder/'definition.json');d.update(version='CaoTamuraBudgetH2gV1',model_specification=spec,model_manifest_sha256=manifest_hash,h_s=cal['h_s'],common_mesh_sha256=hashlib.sha256(canonical_bytes(mesh)).hexdigest(),mesh=mesh,xml_sha256=sha(xml),slots=2804,initializer_strain_weight='common kappa*u_tau*physical incoming height; independent of physical Gamma')
    atomic(folder/'definition.json',d)
    return xml

def prepare(calibration,cls_preflight,mesh_definition,output):
    cal=load(calibration);old=load(cls_preflight);mesh=load(mesh_definition)['mesh']
    # Same actual knot arrays, degree and production quadrature, independent of
    # standard roughness/cutoff. Both models retain their own integration panels.
    mesh=dict(x=mesh['x'],z=mesh['z'],degree=3,regional_pressure_degree=1,points=[16,16])
    specs=[]
    for i,m in enumerate(cal['models']):
        specs.append(dict(version='StandardNeutralCalibratedV1' if i else 'DeployedCLSV1',model=m['name'],coefficients=dict(Cmu=.09,Cepsilon1=1.44,Cepsilon2=1.92,sigmaK=1.,sigmaEpsilon=m['sigmaEpsilon'],YapCoefficient=0 if i else .83),q=m['q'],kappa=m['kappa'],K0=m['K0'],wall=dict(cMu=m['cMu_wall'],cLength=m['cLength_wall'],h_s=m['h_s'],cutoff_policy='native local physical budget root; never prescribed flat cutoff'),common_fit=old['fit'],coordinates='unshifted physical incoming height',H=.0065,u_tau=cal['u_tau'],nu=cal['nu'],z0_eff=cal['z0_eff'],outer_epsilon='u_tau^3/(kappa*z)',upper_boundary='declared compatibility fixture or RoughFlatOpenV1 facility slip roof; no constant tail',physical_terms=dict(quadratic=not bool(i),cubic_beta=0 if i else .4,Gamma=0 if i else 1,gamma0=0,Yap=not bool(i))))
    d=dict(cal,version='H2gComparisonCalibrationV1',cls_preflight=old,common_mesh=mesh,model_manifest_CLS=specs[0],model_manifest_standard=specs[1],fit_sha256=hashlib.sha256(json.dumps(old['fit'],sort_keys=True).encode()).hexdigest())
    # Adapt the already validated native/GL20 certification helpers only through
    # the reference values they use for normalization, never a model dispatch.
    d.update(kappa_CLS=cal['kappa'],heights=[dict(H=.0065,K0=cal['models'][1]['K0'],E_C_H=cal['models'][1]['E_C_H'])])
    atomic(output,d);return d
