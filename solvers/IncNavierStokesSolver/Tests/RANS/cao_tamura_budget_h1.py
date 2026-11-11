#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Accepted budget-AWF H1 adapter. Frozen flat data only; no EP route imports.

Dimensional convention U(Hh)=1 m/s, Hh=.04 m, nu=U(Hh)*Hh/10500.
This is equivalent to the frozen Hh/U(Hh) convention, not a measured speed.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

from run_cls_epawf_picard import ROOT,require,sha,set_info,set_parameter,write,write_json
from run_cls_accepted_awf import case_xml
from cls_epawf_m4_cases import tensor_mesh

DATA=ROOT/'docs/plans/cao-tamura-validation-data'
VERSION='CaoTamuraBudgetH1V1'
HH=.04
UH=1.
UH_OVER_UINF=.6867027875750893  # recomputed from sealed CSV precision at Hh
UINF=UH/UH_OVER_UINF
USTAR=.051*UINF
NU=UH*HH/10500
KAPPA=.440553905186659
Q=3.190675081624288


def frozen():
    snapshot=json.loads((DATA/'snapshot.json').read_text())
    for name,digest in snapshot['files'].items():require(sha(DATA/name)==digest,'frozen data differs: '+name)
    rows=list(csv.DictReader((DATA/'profiles.csv').open()))
    series=json.loads((DATA/'series.json').read_text())
    scalars=list(csv.DictReader((DATA/'reported_scalars.csv').open()))
    require(len(rows)==1151 and len(series)==98 and len(scalars)==39,'frozen inventory differs')
    rough=sorted((HH*float(r['y_over_H']),float(r['value'])) for r in rows if r['series_id']=='f3a_R' and float(r['weight'])>0)
    require(abs(linear(rough,HH)-UH_OVER_UINF)<2e-15,'frozen Hh velocity scale changed')
    return snapshot,rows


def pchip(x,y):
    """Fritsch-Carlson monotone slopes; exact C1 Hermite derivative is native."""
    h=[b-a for a,b in zip(x,x[1:])];d=[(b-a)/v for a,b,v in zip(y,y[1:],h)]
    m=[d[0]]
    for i in range(1,len(x)-1):
        if d[i-1]*d[i]<=0:m.append(0.)
        else:
            w1=2*h[i]+h[i-1];w2=h[i]+2*h[i-1]
            m.append((w1+w2)/(w1/d[i-1]+w2/d[i]))
    m.append(d[-1]);return m


def linear(points,z):
    if z<=points[0][0]:return points[0][1]
    for a,b in zip(points,points[1:]):
        if z<=b[0]:return a[1]+(b[1]-a[1])*(z-a[0])/(b[0]-a[0])
    return points[-1][1]


def make_profile(output):
    output=Path(output);output.mkdir(parents=True,exist_ok=False);snapshot,rows=frozen()
    permitted={}
    for name in ('f3a_R','f3b_R','f3c_R'):
        permitted[name]=[(HH*float(q['y_over_H']),float(q['value']),float(q['weight']),q['quality'])
                         for q in rows if q['series_id']==name and float(q['weight'])>0]
    u=[(0.,0.)]+[(z,v*UINF) for z,v,w,q in permitted['f3a_R']]
    # Preserve the frozen Hh/U(Hh) scale explicitly; this inserted anchor is
    # the allowed Fig3a linear interpolation, not an additional observation.
    u=sorted(u+[(HH,UH)])
    # Recorded last measured mean, including its 0.0104% extraction-level
    # departure from Uinf, is retained; not fitted to any hill observation.
    u += [(.4,u[-1][1]),(.8,u[-1][1])]
    x=[q[0] for q in u];y=[q[1] for q in u];m=pchip(x,y);m[-1]=m[-2]=0.
    rms=[(z,v*USTAR) for z,v,w,q in permitted['f3c_R']]
    # Deliberate one-member completion: Ryy=.5 Rxx, Rzz=.25 Rxx.
    # These missing components are assumptions, not Fig.3d fitting/isotropy.
    # A .3% Uinf isotropic weak-freestream level is reached by z=.4.
    kfree=1.5*(.003*UINF)**2
    kt=[(z,.875*v*v) for z,v in rms]+[(.4,kfree),(.8,kfree)]
    table=[]
    for z,U,dU in zip(x,y,m):
        k=linear(kt,z)
        ell=KAPPA*max(z,.0002)
        ell=min(ell,.025)  # prescribed .025m outer integral-length completion
        epsilon=(Q**-2)**.75*k**1.5/ell
        table.append([z,U,dU,k,epsilon])
    profile=output/'finite-inflow.txt'
    profile.write_text('FiniteInflowV1 SI:z,U,dU_dz,k,epsilon\n'+''.join(' '.join(format(v,'.17g') for v in q)+'\n' for q in table))
    manifest=dict(version=VERSION,data_version=snapshot['version'],snapshot_sha256=sha(DATA/'snapshot.json'),
        data_files=snapshot['files'],source_pdf_sha256=snapshot['source_sha256'],
        dataset_counts=dict(series=98,points=1151,scalar_facts=39),
        application='rough inflow, rough hill; matching rough-flat control',
        dimensional_convention=dict(Hh_m=HH,U_at_Hh_m_per_s=UH,Uinf_m_per_s=UINF,nu_m2_per_s=NU,
            Re_H=10500,U_at_Hh_over_Uinf=UH_OVER_UINF,u_star_over_Uinf=.051,u_star_m_per_s=USTAR,
            status='explicit equivalent scale; absolute Uinf and nu not separately measured',
            roughness_Re_derived=.051*UINF*.0002/NU,roughness_Re_reported=3.97,
            inconsistency='retain approximately 1.8% mismatch; no simultaneous calibration'),
        geometry=dict(Hh_m=HH,L_m=.1,cube_b_m=.005,h_eff_m=None,aerodynamic_z0_m=.0002,
            displacement_d0_m=.0035,nominal_delta_m=.25,roof_m=.8,planned_hill_streamwise_extent_m=[-1.,1.2],
            hill='Hh*cos(pi*x/(2L))^2 for abs(x)<L; zero outside',
            feet_are_explicit_interfaces=True,wall_normal_Hw_candidates_m=[.01,.015,.02]),
        model=dict(formulation='CLSCalibratedBudgetAWF',outer='Regularized',inner_stress='Linear',
            cMu_wall=Q**-2,cLength_wall=KAPPA*Q**1.5,beta=.4,Gamma=1.,gamma0=0.,
            source='P_actual-E_C(K); zero redistribution; no E_C/k reaction',
            endpoint='z_K=z_m=Hw; K=transported k at Gamma_H',
            cutoff_reference='counterfactual flat parallel zero-force p_A0(K), not finite-reference production'),
        inlet=dict(permitted_series=list(permitted),reserved_flat_check=['f3d_R'],
            held_out='all hill series, crest speedup at vertical .01m, ground-extrapolated x_R/Hh=6.5',
            mean='visible Fig3a point interpolation; fixed height weights retained for error reporting',
            shear='Fig3b only used for same-profile tensor discrepancy, not forced into CLS',
            k_completion='k=.5*(Rxx+Ryy+Rzz); Rxx=Fig3c sigma_u^2, Ryy=.5Rxx, Rzz=.25Rxx below measured edge; linear-height k completion',
            epsilon_completion='UnshiftedCLSLengthV1: knot values Cstar^.75*k^1.5/ell, ell=min(kappa*max(z,z0),.025m), then positive loglinear interpolation; unmeasured assumption. Native trace is separately constrained to exact E_C(Hw); report interpolation/projection mismatch',
            freestream='last measured U constant by z=.4; positive .3% Uinf weak turbulence; fixed .025m length',
            low_height='U Hermite to zero at floor; k held at lowest allowed measurement below coverage; epsilon length regularized below z0. d0 converts measured ordinates only, never shifts CLS/Yap distance',
            interpolation='C1 Hermite U/exact derivative, positive loglinear k/epsilon in native provider',
            profile_file=str(profile.resolve()),profile_sha256=sha(profile),knot_columns=['z','U','dU_dz','k','epsilon'],knots=table,
            native_inner_completion='AWF no-slip U anchored at finite U(Hw), nonconstant budget-reference k anchored at completed K(Hw), zero ground k flux; not the raw RMS completion below Hw. The transported EPS is only outer and is constrained natively at Hw',
            reference_frame=dict(origin=[-.2,0,0],flow=[1,0,0],normal=[0,0,1],height='global physical z above upstream floor'),
            fitted_values_modified=False,weighted_mean_error_at_measured_knots=0.),
        data_policy=dict(point_visibility='retain original; zero-weight points not fitted',weights='frozen normalized per-series Voronoi weights',
            uncertainty='extraction envelope is not a statistical confidence interval; retain all frozen ambiguity/coverage rules'),
        boundaries=dict(compatibility_implemented='native external policy fixes both x boundaries and top; ground velocity no-slip, k natural, regional pressures/common gauge',
            precursor_planned='inlet profile, developed/weak outlet, top treatment and flux/gauge qualification remain H2'),
        old_generator_reuse='data definitions and affine tensor topology only; no EP/collar/P-N/fallback/old wall coefficients')
    queries=[]
    for series,quantity,scale in [('f3a_R','U',UINF),('f3b_R','minusRxz',UINF**2),('f3c_R','sigmaU',USTAR)]:
        for r in rows:
            if r['series_id']==series and float(r['weight'])>0:
                queries.append(f"{HH*float(r['y_over_H']):.17g} {quantity} {scale*float(r['value']):.17g} {float(r['weight']):.17g} {scale*float(r['value_bound']):.17g}\n")
    (output/'permitted-queries.txt').write_text(''.join(queries))
    manifest['inlet']['query_file']=str((output/'permitted-queries.txt').resolve())
    manifest['inlet']['query_sha256']=sha(output/'permitted-queries.txt')
    write_json(output/'manifest.json',manifest);return manifest


def compatibility_case(folder,profile,preflight,steps=8,dt=1e-5,restart=None,reject=0):
    """Small physical rough-flat pilot with both streamwise sides prescribed."""
    profile=Path(profile).resolve();preflight=json.loads(Path(preflight).read_text())
    H=preflight['selected']['Hw'];h=preflight['selected']['h_eff']
    folder=Path(folder)
    xml=case_xml('Accepted-AWF-05G',folder,steps,dt=dt,restart=restart,companion='regularized',
        companion_reason='finite-inflow compatibility only, not neutral/experimental precursor',checkpoint_every=1)
    root=ET.parse(xml).getroot()
    z=sorted(set([0.,.00025,.001,.003,H,.03,.06,.1,.2,.25,.3,.4,.5,.8]))
    x=[-.2,-.1,0.,.1,.2]
    tensor_mesh(root,z,3,2,.4,1.,points=[8,8],xknots=x)
    for key,value in dict(RANSWallReference='FiniteInflowV1',RANSWallRegularizedReference='FiniteInflowV1',
        RANSWallReferenceFile=str(profile),RANSWallReferenceSHA256=sha(profile),RANSWallComparison='Off',
        RANSWallEvaluationReuse='ExactStepV1').items():set_info(root,key,value)
    for key,value in dict(RANSWallHeight=H,RANSWallSampleHeight=H,RANSWallTailStart=H,
        RANSWallRoughness=h,Kinvis=NU,RANSWallFrictionVelocity=USTAR,RANSWallReferenceOriginX=-.2,
        RANSWallPimpleHardAuditRejectStep=reject,K0=Q*USTAR*USTAR,
        CaoTamuraHillHeight=HH,CaoTamuraHalfFootprint=.1,CaoTamuraCubeHeight=.005,
        CaoTamuraAerodynamicZ0=.0002,CaoTamuraDisplacementD0=.0035,
        CaoTamuraNominalDelta=.25,CaoTamuraUAtHillHeight=UH,CaoTamuraUInf=UINF,
        RANSWallMapAmplitude=HH,RANSWallMapLengthX=.1).items():set_parameter(root,key,value)
    # All old source/perturbation helpers must be identically zero: no inherited
    # cancellation, terrain/EP policy or old neutral helper boundary data.
    for fun in root.findall('./CONDITIONS/FUNCTION'):
        if fun.get('NAME','').startswith(('BodyForce','ScalarSource','InitialPerturbation')):
            for e in fun.findall('E'):e.set('VALUE','0')
    write(root,xml)
    ni=sum(v<=H for v in z[1:]);no=len(z)-1-ni;nx=len(x)-1
    sizes=[(3*nx+1)*(3*(ni+no)+1)]*3+[(3*nx+1)*(3*no+1),(nx+1)*(ni+1),(nx+1)*(no+1),1]
    require(sum(sizes)<=3000,'compatibility mesh exceeds retained dense backend')
    definition=dict(version=VERSION,case='rough-flat-finite-profile-compatibility-tail-v1',physical_steps=steps,dt=dt,
        endpoint_time=steps*dt,scaling='SI U(Hh)=1, Hh=.04, Re_H=10500',Hw=H,z_K=H,z_m=H,h_eff=h,
        mesh=dict(x=x,z=z,velocity_scalar_degree=3,pressure_degree=1,points=[8,8],
            revision='completion-tail-v1',added_outer_knots=[.25,.4],reason='preserved coarse T0 epsilon undershoot near .304; resolve steep measured-to-freestream completion, no field clipping'),
        native_size_estimate=sum(sizes),native_field_sizes=sizes,
        profile_sha256=sha(profile),fresh_native_initialization=restart is None,
        boundary_scope='both streamwise boundaries and top fixed; compatibility test, not developed experimental precursor',
        physical_map='affine flat; hill production remains guarded',initialization='native finite profile projection, AWF inner U/k completion, dependent epsilon trace Lift once during setup',
        restart=str(restart) if restart else None,xml_sha256=sha(xml))
    (folder/'definition.json').rename(folder/'unexecuted-template-definition.json')
    write_json(folder/'definition.json',definition)
    write_json(folder/'h1-definition.json',definition);return xml,definition


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();make_profile(a.output);print(a.output/'manifest.json')


if __name__=='__main__':main()
