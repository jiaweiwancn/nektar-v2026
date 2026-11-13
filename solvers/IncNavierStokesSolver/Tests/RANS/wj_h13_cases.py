#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Reusable full-height neutral paired-domain preparation. Never executes native code."""
import argparse
import csv
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.integrate import quad
from wj_h5_neutral import NeutralReference,NE,QE,CD,KAPPA,SCALE_MM,observation_regions
from wj_h5_cases import case_xml,mesh_knots
from cls_epawf_m3_cases import set_info,write
from wj_h13_domain import X,H,D0,TOP,L,identity,digest,topology,mapped,dof_counts,storage,source_capacity,preflight

def save(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')

def model_metadata():
    return dict(outer='WJHighRe',tensor='full3D Original',gradient_policy='DevStrainRawWorkV1',production='-R:Graw',
        inner='WJCalibratedBudgetAWF:Linear',c1=1.8,c2=5/9,Ce1=1.44,Ce2=1.82,sigma_k=1.,sigma_epsilon=1.30,
        Ne=NE,qe=QE,Cstar=CD,Cd=CD,kappa=KAPPA,outer_molecular_viscosity=0.,wall_molecular_viscosity=.051*.005/3.97,
        cutoff='metric counterfactual production at current represented K; actual P_A-E_c drives transported k_w',
        dependent_epsilon='A(K_current)/H, outer only',chi='no physical chi',source_caps=False)

def primary(source):
    record=Path(source)/'docs/plans/wj-h5-linux/calibration.json';data=json.loads(record.read_text())
    v=next(v for v in data['candidates']if v['H_w_mm']==6.)
    c=NeutralReference(v['H_w'],v['z0_e'])
    for k,a in (('H_w',c.H),('z0_e',c.z0_e),('h_s',c.h),('u_tau',c.u_tau),('nu',c.nu),('d0',c.d0),('top',c.top),('K',c.K),('A',c.A)):
        if float(a).hex()!=float(v[k]).hex():raise ValueError('retained reference bytes differ: '+k)
    return c,dict(calibration=identity(record),frozen_primary=v,historical_heights=[{k:r[k]for k in ('H_w_mm','h_s_mm','h_s','physical_z_H_mm')}for r in data['candidates']],
                  lower_mean_fit=data['fits']['primary'],no_new_fit=True)

def boundary(c):
    rT=c.top-c.d0;B=CD*c.K*c.K/1.30
    return dict(pressure=dict(variable='pi=pbar+2k/3',inner=-c.K/NE,outer=0.,gauge='absent; complete outlet traction fixes pressure level',role='initial/reference seed only'),
        inlet=dict(U='u_tau/kappa*log((physical_z-d0)/z0_e) above H; finite-nu J0 below',W=0.,k='nonconstant k_w below H; K above',epsilon='outer only A/r; H belongs to dependent current-K trace',no_pressure_Dirichlet=True,ground_corner_precedence='no slip'),
        top=dict(selector='ReferenceTractionV1',physical_z=c.top,physical_r=rT,normal_velocity=0.,tangent_traction=[c.u_tau**2,0.,0.],normal_traction='current impermeability reaction; not prescribed',
            k_outward_diffusive_flux=0.,epsilon_outward_diffusive_flux=B/rT,epsilon_signed_gradient_flux=-B/rT,
            unspecialized_coefficient='B=Cd*K^2/sigma_epsilon=(Cd/Cstar)*u_tau^4/sigma_epsilon',B=B,
            epsilon_equation_coefficient=(1.82-1.44)*c.A*c.A/c.K,epsilon_diffusion_derivative='B/r^2',integrated_gradient_flux=B*(1/c.H-1/rT),
            off_reference_law='fixed physical outward reference flux; not fixed gradient, dynamic Robin or Dirichlet'),
        outlet=dict(normal=[1,0,0],inner=[c.K/NE,0.,c.u_tau**2],outer=[-c.K/NE,0.,c.u_tau**2],scalar_outward_diffusive_flux=[0.,0.],backflow='native explicit stop; no clamp/extrapolation'),
        ground=dict(velocity=[0.,0.,0.],k_flux=0.,inner_epsilon_PDE=False),
        interface=dict(U='shared',k='shared endpoint K',epsilon='dependent outer value A(K_current)/H',normal_stress='conservative weak coupling; no extra pressure jump equation',individual_R_continuity=False))

def clocks(c,x):
    rT=c.top-c.d0;length=x[-1]-x[0]
    # Independent physical cross-section integration; no native projection/mass.
    cuts=sorted(set([0.,c.H]+[v for v in (c.w['rc'],c.h)if 0<v<c.H]))
    ui=sum(quad(c.U_inner,a,b,epsabs=2e-13,epsrel=2e-13)[0]for a,b in zip(cuts[:-1],cuts[1:]))
    uo=c.u_tau/KAPPA*((rT*math.log(rT/c.z0_e)-rT)-(c.H*math.log(c.H/c.z0_e)-c.H))
    bulk=(ui+uo)/rT;Tb=length/bulk;slow=length/c.U_outer(c.H)
    inner_diff=c.H*c.H/min(c.D(0),c.D(c.H));top_turb=c.K/(c.A/rT)
    horizon=max(3*Tb,2*slow,inner_diff,top_turb)
    return dict(Lx=length,cross_section_physical_height=rT,volume_per_unit_span=length*rT,reference_flow_rate_per_unit_span=ui+uo,U_bulk=bulk,T_bulk=Tb,
        T_slowest_outer_transit=slow,T_top_transit=length/c.U_outer(rT),T_inner_diffusion_conservative=inner_diff,
        T_interface_k_over_epsilon=c.K*c.H/c.A,T_top_k_over_epsilon=top_turb,
        proposed_preservation_horizon=horizon,horizon_rule='max(3*T_bulk,2*T_slowest_outer,T_inner_diffusion,T_top_k/epsilon)',
        terminal_window=[horizon-slow,horizon],horizon_native_admission='NOT_RUN; future measured cost and frozen application budgets required',
        temporal_steps_at_metadata_dt_0p01=math.ceil(horizon/.01),runtime_estimate=None,
        note='slowest OUTER transit excludes no-slip infinite advection time; clocks and finite monitoring do not prove stationarity')

def protocol(c,t):
    x=t['xknots'];stations=sorted(set(x+[-2.5,0.,2.5,1.25,3.75,5.,6.25,7.5]+[3.75+.25*i for i in range(18)]));time=clocks(c,x)
    return dict(schema='WJH13TargetPreservationProtocolV1',status='PREPARED_NOT_ADMITTED',fresh_mesh_signature_required=True,old_checkpoint_conversion=False,
        initial_states=dict(q0='analytical coupled neutral reference',qh0='actual accepted state after projection, external/dependent lifting and FlatInletTraceV1',
            internal_reference='separate native pre-initializer frozen reference/defect; not assumed equal to qh0',qh0_identity=None),
        stages=dict(A=dict(kind='future fresh target CP0 plus short admission pilot',maximum_pilot_steps=4,metadata_dt=.01,admission='independent future contract, capacity/runtime and measured cost first',
            CP0_checks=['native IDs/signs and trace compatibility','all regional Q1 continuity moments <=2e-15','physical div <=1e-12*max(1,||G||)','unchanged scalars/root/covariance/principal/cutoff/outlet/ground predicates',
                        'half-neutral quality m/nu_eff>=.002068763923454, distinct from hard physical sign','actual K, cutoff and epsilon_H=A(K)/H','fixed reference/defect and actual qh0 distinct']),
            B=dict(kind='target-domain duration preservation',clocks=time,application_error_budgets=None,comparison_floors=None,
                   before_native_launch='freeze budgets/floors/horizon/monitor cadence and finite resources; no result-dependent tuning',
                   terminal_monitoring='sample station and inventory drift over the complete declared terminal window; report trend, not a stationarity proof',
                   insufficient_duration='PARTIAL_DURATION/UNRESOLVED'),C=dict(kind='separately authorized targeted temporal/spatial sensitivity',automatic=False)),
        errors=dict(E_repr='||qh0-q0||/S',E_drift='||qh(t)-qh0||/S',E_ref='||qh(t)-q0||/S',E_fetch='||qh(x,t)-qh(x_in,t)||/S',station_to_analytical='also compare every station with q0',
            temporal_difference='separate diagnostic, not preservation or continuum error'),
        norms=dict(volume='physical regional Gauss8 L2, faces and wall/cutoff/roughness breakpoints split; Jacobian once',finite_maximum='all normal faces with both owners, native points and deterministic inward offsets; freeze exact inventory before native run',
            scales=dict(U=1.,alternative_velocity=c.u_tau,k=c.K,R=c.u_tau**2,pressure=c.u_tau**2,epsilon_fixed=c.A/c.H),epsilon_additional='outer-only height-resolved relative error; do not hide top errors with A/H scale',zero_components='absolute fixed scales, no division by no-slip U=0',pressure='traction-fixed convention; no independent regauge'),
        stations=[dict(x_code=s,x_mm=40*s,role='inlet'if s==x[0]else'outlet'if s==x[-1]else'upstream_foot'if s==-2.5 else'downstream_foot'if s==2.5 else'virtual_crest'if s==0 else'profile/wake diagnostic',
                       mesh_vertex=s in x,physical_height='original floor z; flat r=z-d0; retain both face/H owners')for s in stations],
        required_diagnostics=['U,W,k,outer epsilon,pi,all Rij','outer -Rxz/u_tau^2, inner total molecular+turbulent shear, wall traction','P/epsilon where meaningful, actual versus counterfactual production/destruction',
            'current K/cutoff/dependent trace','signed storage/advection/stress/roughness/diffusion/source/exterior/interface budgets','local test moments, regional inventories, weak continuity and strong divergence separately','limiter/bound/fallback and raw/corrected differences','continuum, expression, projection, quadrature, assembled spatial residual and genuine completed-stage defects separately'],
        no_equilibrium_repairs=['no eta anchor','no EP redistribution','no reference reset','no extra forcing','no retuning'],
        experimental_roles='Fig3 at flat virtual crest x=0, not moved to inlet=-10; selected neutral inlet is an explicit model choice, not measured finite-depth reconstruction')

def next_backend(source,count):
    source=Path(source);eq=Path('solvers/IncNavierStokesSolver/EquationSystems')
    needles=[(eq/'RANSWallGlobalField.cpp','RANSWallGlobalField::WallMatrix RANSWallGlobalField::Mass()'),
        (eq/'RANSWallGlobalField.cpp','RANSWallGlobalField::Vector RANSWallGlobalField::Project('),
        (eq/'RANSWallEvolutionPatch.cpp','RANSWallEvolutionPatch::WallMatrix RANSWallEvolutionPatch::PhysicalMass('),
        (eq/'RANSWallEvolutionPatch.cpp','RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::PhysicalProject('),
        (eq/'RANSWallEvolutionPatch.cpp','RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Norms('),
        (eq/'RANSWallEvolutionPicard.cpp','system.matrix.assign(n, Vector(n))'),
        (eq/'RANSWallEvolutionPimple.cpp','WallMatrix a(nv,Vector(nv))'),
        (eq/'RANSWallEvolutionPimple.cpp','const Factor momentumFactor(a)'),
        (eq/'RANSWallEvolutionPimple.cpp','WallMatrix bg(np,Vector(nv)),schur(np,Vector(np))'),
        (eq/'RANSWallEvolutionPimple.cpp','const auto delta=Factor(block.matrix).Apply(rhs)'),
        (eq/'RANSWallEvolutionKBound.cpp','std::unique_ptr<RANSWallKBound::Projection>'),
        (eq/'RANSWallEvolutionDiagnostics.cpp','{"pressure_mass",J::Rows(m_mass[b])}'),
        (Path('library/LibUtilities/LinearAlgebra/NekLinSysIterGMRES.h'),'class NekLinSysIterGMRES')]
    locations=[]
    for path,needle in needles:
        lines=(source/path).read_text().splitlines();hits=[i+1 for i,s in enumerate(lines)if needle in s]
        if not hits:raise ValueError('backend source call graph moved: '+str(path)+' '+needle)
        locations.append(dict(relative_path=str(path),lines=hits,needle=needle,identity=identity(source/path)))
    return dict(schema='WJH13MappedBackendPrerequisiteV1',status='IMPLEMENTATION_REQUIRED',target=count,
        proposed_opt_in='MappedBlockV1 (proposal only; keep legacy dense default and guards)',source_call_graph=locations,
        smallest_complete_extension=[
            'element/face/ray-local sparse or matrix-free mass and frozen Picard actions using existing native IDs/signs and full mapped rows; consistent mass in startup projection as well as governing steps',
            'replace global dense PhysicalMass/Project and residual mass-dual factors with qualified consistent-mass applications/solves; preserve physical quadrature, constraints and reference freeze order',
            'predictor solve with alpha*M+L; nonsymmetric Schur action D*(alpha*M_free)^-1*G-C via actual complete pressure action and free-row restrictions; no dense BG column corpus',
            'k then epsilon transport, retaining dependent epsilon identity rows AND their transport/mass columns, current-K profiles/cutoff and refreshed boundary reactions',
            'operator-backed bound mass inverse/retained constraint pool, diagnostics/strict reader and streamed sparse rows; do not leave dense startup/observation hidden behind a new governing backend selector'],
        existing_infrastructure='LibUtilities NekLinSysIterGMRES/NekSysOperators and native AssemblyMapCG; requires a custom mapped-row adapter, not only a session GlobalSysSoln change',
        forbidden_shortcuts=['mass lumping','full momentum inverse in place of mass','Gp=-D^T','omitted roughness/exterior/interface pressure terms','assumed SPD Schur','physics/space/source changes','normal coarsening or guard-only increase'],
        qualification_sequence=['small existing dense-qualified mapped cases, including H10 two-column synthetic/native fixture: mass/flow/pressure/scalar actions and solves versus dense with identical constraints',
            'preserve linear residual tolerances, native hard/quality checks, current-K trace, bound adverse witnesses and one-sweep counters; freeze scaled comparison tolerances before tests',
            'restart/history/initializer invocation and rollback equivalence; pressure level and signed exterior/interface identities; input/generation cache invalidation',
            'only then actual 24688-state target fresh CP0 with finite memory/output/worker budget and all source-aware dimensions/native row maps verified',
            'retained short target pilot measures startup, operator, solve, predicate and observation cost; admit full-duration preservation only from those measurements'],
        measured_target_cost=None,historical_timing='H2 small-case four sweeps Off10.2503s/ExactStep3.55065s, domain.compute11.1269s/2.14574s, dense factorization .01698/.01680s; timers overlap; not target performance',
        implementation_in_H13=False)

def prepare(output,source):
    output=Path(output);output.mkdir(parents=True,exist_ok=False);source=Path(source);c,cal=primary(source)
    flat=case_xml(output/'flat',c,bisect_once=True,xknots=X);tree=ET.parse(flat).getroot()
    set_info(tree,'RANSWallInitialVelocity','FlatInletTraceV1');set_info(tree,'RANSWallDenseEvolutionCapacity','3000')
    tree.insert(0,ET.Comment('NOT_RUN: 24688-state mapped target. The external preflight MUST refuse before native allocation; NumSteps=0 is not a capacity guard.'))
    write(tree,flat);t=topology(flat);count=dof_counts(t);inner,outer,z=mesh_knots(c,bisect_once=True)
    if (count['Nx'],count['Ni'],count['No'],count['state'])!=(14,16,34,24688):raise ValueError('wrong frozen target')
    boundaries=boundary(c);policy=source_capacity(source)
    partners={}
    for name,hill in (('flat',False),('terrain',True)):
        vertices=[]
        for key,p in t['vertices'].items():
            region='inner'if p[1]<=c.H else'outer';v=mapped(p[0],p[1],hill,region=region)
            vertices.append(dict(id=int(key),chart=[p[0],0.,p[1]],physical=v['physical']))
        value=dict(schema='WJH13PairedGeometryV1',partner=name,shared_topology_sha256=t['sha256'],map='Plane'if not hill else'CosineHill',amplitude=0. if not hill else 1.,
            continuation='CompositeTunnelV1',H=c.H,d0=c.d0,top=c.top,hill_half_width=L,slopeX=0.,vertices=vertices,
            native_route_status='NOT_RUN_CAPACITY'if not hill else'GEOMETRY_ARTIFACT_ONLY: curved physical OpenTunnel/initialization unqualified')
        p=output/(name+'-geometry.json');save(p,value);partners[name]=identity(p)
    # Geometry-only terrain XML contains no CONDITIONS/solver or runnable route.
    terrain=ET.Element('NEKTAR');terrain.append(ET.fromstring(ET.tostring(tree.find('GEOMETRY'))));write(terrain,output/'terrain-chart-topology.xml')
    manifest=dict(schema='WJH13FullHeightPairedDomainV1',objective='user-selected full-height neutral inlet; supersedes finite-depth/ambient designs',
        primary_reference=c.summary(),model=model_metadata(),reference_scope='flat analytical reference/calibration; terrain artifact is geometry only and has no initialized field or imported flat cutoff',
        calibration_lineage=cal,shared_topology=t,inner_normal_knots=inner,outer_physical_flat_r=outer,
        normal_chart_knots=z,physical_flat_knots=[c.d0+r for r in inner]+[c.d0+r for r in outer[1:]],
        geometry_partners=partners,boundary_data=boundaries,counts=count,storage=storage(count),capacity_policy=policy,flat_session=identity(flat),
        quadrature=dict(native_representation_GLL=[8,8],canonical_mass_projection_Gauss=8,external_lift_GLL=8,physical_mass_projection_Gauss=16,operator_tangential=12,operator_normal=16,expression_degree=16),
        experimental=dict(Fig3_source_plane='flat virtual crest x=0; data not moved upstream to selected inlet=-10',physical_7mm_effective_r_mm=3.5,
            physical_7mm_region='inner AWF',primary_inverse_lower_H_mm=4.413515542,all_7mm_outside_upper_H_mm=3.5,intersection='EMPTY',no_physical_height_optimality_claim=True),
        native_execution=dict(status='NOT_RUN',all_invocations=0,reason='capacity prerequisite; hill physical OpenTunnel separately unqualified'))
    save(output/'paired-domain.json',manifest);save(output/'preservation-protocol.json',protocol(c,t));save(output/'mapped-backend-prerequisite.json',next_backend(source,count))
    observations=observation_regions(c)
    save(output/'observations.json',dict(flat_reference_plane_x=0.,inlet_x=X[0],rows=observations,
        terrain_anchors=[dict(name='crest_speedup',x=0.,physical_z=1.+.25,vertical_clearance=.25,region_at_crest='outer; effective r=.1625>H',not_fit_target=True),
            dict(name='reattachment',x=6.5,role='wall traction/dividing streamline, not an invented velocity height'),
            dict(name='downstream_7mm_U',x_range=[3.75,8.],x_interval=.25,physical_z=.175,effective_r=.0875,region='inner AWF',epsilon='no inner transported field',provenance='Cao–Tamura printed13:150–320mm,10mm interval,7mm above floor'),
            dict(name='Fig12_downstream_7mm_rms',physical_z=.175,effective_r=.0875,region='inner AWF',digitized_values=None,provenance='printed15 Fig12, height stated printed13; no new digitization')],
        height_rule='z is original-floor physical height; y_prime=z-f(s) is vertical clearance, not wall-normal r; future hill uses certified inverse and owner enumeration'))
    guard=preflight(flat,manifest,source);save(output/'native-preflight.json',guard)
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path);p.add_argument('--source',type=Path,default=Path(__file__).resolve().parents[4]);p.add_argument('--preflight',type=Path);a=p.parse_args()
    if a.preflight:
        m=json.loads(a.preflight.read_text());v=preflight(m['flat_session']['path'],m,a.source);print(json.dumps(v,indent=2));raise SystemExit(2 if not v['native_launch_allowed']else 0)
    if a.output is None:p.error('--output or --preflight required')
    m=prepare(a.output,a.source);print('PREPARED',m['counts']['state'],'native NOT_RUN: capacity OBSTRUCTED')
