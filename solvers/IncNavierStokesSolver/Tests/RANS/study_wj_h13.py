#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded offline H13 assessment. Native executables are never invoked."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time
import xml.etree.ElementTree as ET
import numpy as np
from scipy.integrate import quad
from wj_h13_cases import primary,prepare,save,boundary
from wj_h13_domain import H,D0,TOP,X,identity,topology,geometry_checks,mapped,dof_counts,preflight
from wj_h6_projection import expression
from wj_h5_neutral import KAPPA,CD,NE,observation_regions
from wj_h2_profiles import constitutive

def table(p,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('x',newline='')as f:w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)

def function_values(xml):
    return {(f.get('NAME'),e.get('VAR')):e.get('VALUE')for f in ET.parse(xml).findall('./CONDITIONS/FUNCTION')for e in f.findall('E')}

def selected_parent(source,new_xml,z):
    root=Path('/home/jay/build/wj-h10-flat-inlet-20261113');local=Path(source)/'docs/plans/wj-h10-linux/publication-v1'
    for name in ('manifest.json','inputs.json','seal.json'):
        if identity(local/name)['sha256']!=identity(root/'publication-v1'/name)['sha256']:raise ValueError('H10 source manifest differs')
    inputs=json.loads((local/'inputs.json').read_text());entry=next(v for v in inputs['cases']if v['branch']=='coarse')['cp0']
    if identity(entry['path'])!=entry:raise ValueError('selected retained H10 XML bytes differ')
    manifest=json.loads((local/'manifest.json').read_text());assert entry in manifest['heavy']
    old=ET.parse(entry['path']).getroot();oldz=sorted(set(float(v.text.split()[1])for v in old.findall('./GEOMETRY/VERTEX/V')))
    if np.asarray(oldz,dtype='<f8').tobytes()!=np.asarray(z,dtype='<f8').tobytes():raise ValueError('H6 native normal knot bytes differ')
    if function_values(entry['path'])!=function_values(new_xml):raise ValueError('frozen physical Function bytes differ')
    return dict(status='PASS',parent_XML=entry,normal_chart_binary64_sha256=hashlib.sha256(np.asarray(z,dtype='<f8').tobytes()).hexdigest(),
        physical_Function_values_byte_identical=True,parent_mesh='Nx1 original; reuse input/reference evidence only; no native launch',
        new_mesh='Nx14 fresh signature/history required; no old checkpoint conversion',selected_manifest=identity(local/'manifest.json'))

def reference_checks(c,xml):
    b=boundary(c);B=b['top']['B'];rT=c.top-c.d0
    functions={key:expression(v)for key,v in function_values(xml).items()};errorU=0.;errork=0.;rows=[]
    for r in sorted(set(np.linspace(0,c.H,101).tolist()+[c.w['rc'],c.h,c.H])):
        errorU=max(errorU,abs(float(functions['InitialConditionsInner','u'](0,c.d0+r))-c.U_inner(r)))
        errork=max(errork,abs(float(functions['InitialConditionsInner','k'](0,c.d0+r))-c.k_inner(r)))
    if max(errorU,errork)>2e-11:raise ValueError('reference expression approximation failed')
    p_integral=quad(c.production,0,c.H,epsabs=2e-14,epsrel=2e-12,points=[c.w['rc'],c.h])[0]
    e_integral=quad(c.destruction,0,c.H,epsabs=2e-14,epsrel=2e-12,points=[c.w['rc']])[0]
    eps_integral=quad(lambda r:B/r**2,c.H,rT,epsabs=2e-15,epsrel=2e-12)[0]
    checks=dict(expression_U_max_abs=errorU,expression_k_max_abs=errork,budget_production_integral=p_integral,budget_destruction_integral=e_integral,
        budget_integral_defect=e_integral-p_integral,Q_endpoints=[c.w['Q'](0),c.w['Q'](c.H)],
        epsilon_B=B,epsilon_reaction_coefficient=b['top']['epsilon_equation_coefficient'],epsilon_coefficient_defect=B-b['top']['epsilon_equation_coefficient'],
        epsilon_gradient_flux_integral=eps_integral,epsilon_boundary_difference=b['top']['integrated_gradient_flux'],
        epsilon_integral_defect=eps_integral-b['top']['integrated_gradient_flux'],top_outward_epsilon_flux=b['top']['epsilon_outward_diffusive_flux'],
        interface_U_difference=c.U_inner(c.H)-c.U_outer(c.H),interface_K_difference=c.k_inner(c.H)-c.K,
        total_inner_traction=c.D(c.H)*c.u_tau**2/c.D(c.H),outer_shear=c.u_tau**2,finite_nu_gradient_ratio=(c.u_tau**2/c.D(c.H))/(c.u_tau/KAPPA/c.H),
        dependent_epsilon=c.A/c.H,wall_E_H=c.destruction(c.H),pi_jump=-c.K/NE,
        full_height_top_U_over_Uinf=c.U_outer(rT),molecular_ratio_at_H=c.nu/(KAPPA*c.u_tau*c.H),
        no_native_projection_or_discrete_spatial_residual_evaluated=True)
    if abs(checks['budget_integral_defect'])>2e-12*c.A or abs(checks['epsilon_integral_defect'])>2e-14 or abs(checks['epsilon_coefficient_defect'])>2e-18:
        raise ValueError('continuum scalar budget identity failed')
    # Full3D independent Python WJ algebra, not a native probe.
    for r in (c.H,.25,1.,6.25,rT):
        G=np.zeros((3,3));G[0,2]=c.u_tau/KAPPA/r;R,nt,N=constitutive(c.K,c.A/r,G)
        if max(float(np.max(abs(R-c.Rc))),abs(nt-KAPPA*c.u_tau*r),abs(N-NE))>2e-12:raise ValueError('full3D neutral WJ reference mismatch')
        rows.append(dict(r=r,physical_z=c.d0+r,N=float(N),nu_eff=float(nt),nu_d=CD*c.K*c.K/(c.A/r),raw_production=float(-np.sum(R*G)),epsilon=c.A/r,
            covariance=np.linalg.eigvalsh(R).tolist()))
    checks.update(status='PASS',outer_reference_rows=rows,qualification='continuum/reference expression checks only; no native equilibrium preservation result')
    return checks

def profiles(c,radii_knots):
    # Preserve the generator's original flat-r values, especially exact H;
    # subtracting d0 from a rounded physical sum can lose its interface role.
    rphys=list(radii_knots)
    observed=[v['r_mm']/40 for v in observation_regions(c)]
    radii=sorted(set([0.,c.H,c.w['rc'],c.h]+rphys+observed+[c.top-c.d0]))
    rows=[]
    for r in radii:
        if not 0<=r<=c.top-c.d0:continue
        for region in (('inner','outer')if r==c.H else('inner',)if r<c.H else('outer',)):
            if region=='inner':
                U=c.U_inner(r);k=c.k_inner(r);ep=None;G=c.u_tau**2/c.D(r);R=c.covariance_inner(r);nu=c.nt(r);pi=-c.K/NE
                prod=c.production(r);destruction=c.destruction(r);total_shear=c.D(r)*G;molecular_shear=c.nu*G
            else:
                U=c.U_outer(r);k=c.K;ep=c.A/r;G=c.u_tau/KAPPA/r;R=c.Rc;nu=KAPPA*c.u_tau*r;pi=0.;prod=ep;destruction=ep;total_shear=c.u_tau**2;molecular_shear=0.
            row=dict(physical_z_code=c.d0+r,physical_z_mm=40*(c.d0+r),effective_r_code=r,effective_r_mm=40*r,region=region,
                U_over_Uinf=U,W=0.,k=k,epsilon_outer=ep,pi=pi,mechanical_pressure=pi-2*k/3,Gxz=G,nu_t=nu,
                total_shear=total_shear,molecular_shear=molecular_shear,production=prod,destruction=destruction,
                k_flux=-c.w['Q'](r)if region=='inner'else 0.,covariance_min=float(np.linalg.eigvalsh(R).min()))
            for i,v in enumerate(R.ravel()):row['R'+str(i//3)+str(i%3)]=v
            rows.append(row)
    return rows

def assess(root):
    start=time.monotonic();source=root/'source';out=root/'science';out.mkdir(exist_ok=False)
    manifest=prepare(out/'paired-case',source);c,_=primary(source);xml=manifest['flat_session']['path'];t=topology(xml)
    save(out/'parent-input-equivalence.json',selected_parent(source,xml,manifest['normal_chart_knots']))
    ref=reference_checks(c,xml);save(out/'reference-boundary-checks.json',ref)
    geometry={name:geometry_checks(t,hill)for name,hill in (('flat',False),('terrain',True))};save(out/'paired-geometry-checks.json',geometry)
    table(out/'reference-profiles.csv',profiles(c,manifest['inner_normal_knots']+manifest['outer_physical_flat_r'][1:]))
    table(out/'counts-storage.csv',[dict(case='target14',**v)for v in manifest['storage']['entries']])
    # Check compatibility/blocked dimensions on the actual target, not by launch.
    gate=preflight(xml,manifest,source);assert gate['status']=='OBSTRUCTED'and not gate['native_launch_allowed']
    result=dict(paired_domain_construction='PASS',continuum_boundary_compatibility=ref['status'],geometry_certificate='PASS',target_capacity='OBSTRUCTED',
        fresh_native_initialization='NOT_RUN',target_domain_evolution='NOT_RUN',duration_preservation='NOT_RUN',experimental_comparison='NOT_RUN',
        all_native_invocations=0,native_builds=0,global_dense_matrix_allocations=0,target_count=manifest['counts']['state'],
        shared_topology_sha256=t['sha256'],source=identity(__file__),elapsed_science_seconds=time.monotonic()-start)
    save(out/'component-gates.json',result)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(10,3.7),layout='constrained')
    s=np.linspace(X[0],X[-1],601)
    for hill,name in ((False,'flat'),(True,'terrain')):
        ground=np.array([mapped(float(x),0.,hill)['physical']for x in s]);interface=np.array([mapped(float(x),H,hill)['physical']for x in s])
        ax[0].plot(40*ground[:,0],40*ground[:,2],label=name+' effective wall');ax[0].plot(40*interface[:,0],40*interface[:,2],ls='--',label=name+' matching H')
    for x in X:ax[0].axvline(40*x,color='.85',lw=.5)
    ax[0].set(xlabel='physical x (mm)',ylabel='original-floor physical z (mm)',ylim=(0,60),title='Shared chart columns; distinct physical maps (wall zoom)');ax[0].legend(fontsize=7)
    rr=np.geomspace(c.H,c.top-c.d0,200);ax[1].plot([c.U_outer(r)for r in rr],40*(c.d0+rr),label='selected full-height neutral U')
    ri=np.linspace(0,c.H,80);ax[1].plot([c.U_inner(r)for r in ri],40*(c.d0+ri),label='finite-nu Budget wall U')
    ax[1].axvline(1.,color='.5',ls=':',label='U_infinity normalization');ax[1].set(xlabel='U/U_infinity',ylabel='physical z (mm)',ylim=(0,800),title='Idealized reference: top speed is not forced to 1');ax[1].legend(fontsize=7)
    fig.savefig(out/'paired-domain-reference.png',dpi=160);plt.close(fig)
    print('H13 offline assessment PASS; target capacity OBSTRUCTED, all native stages NOT_RUN',json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();assess(a.root)
