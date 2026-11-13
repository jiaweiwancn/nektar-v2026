#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Offline paired chart topology, composite maps and integer capacity accounting.

No native executable, projection, global mass/momentum matrix or solve is used.
Topological keys are independent canonical IDs, not exported native numbering.
"""
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

X=(-10.,-7.5,-5.,-2.5,-1.25,0.,1.25,2.5,4.,5.5,6.5,7.5,10.,12.5,15.)
H=.15;D0=.0875;TOP=20.;L=2.5
GLL5=(-1.,-math.sqrt(3/7),0.,math.sqrt(3/7),1.)

def digest(v):return hashlib.sha256(json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
def identity(p):
    p=Path(p);return dict(path=str(p),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())

def topology(xml,H=H):
    r=ET.parse(xml).getroot()if not isinstance(xml,ET.Element)else xml
    g=r.find('GEOMETRY');vertices={int(v.get('ID')):list(map(float,v.text.split()))for v in g.findall('./VERTEX/V')}
    edges={int(e.get('ID')):list(map(int,e.text.split()))for e in g.findall('./EDGE/E')}
    seen={}
    for key,pair in edges.items():
        if len(pair)!=2 or pair[0]==pair[1]or any(p not in vertices for p in pair):raise ValueError('invalid edge')
        tag=tuple(sorted(pair))
        if tag in seen:raise ValueError('duplicate rather than shared edge')
        seen[tag]=key
    elements=[];owners={e:[]for e in edges}
    for el in g.findall('./ELEMENT/Q'):
        eid=int(el.get('ID'));ee=list(map(int,el.text.split()));vv=set(v for e in ee for v in edges[e]);p=[vertices[v]for v in vv]
        if len(ee)!=4 or len(vv)!=4:raise ValueError('nonquad topology')
        xs=sorted(set(v[0]for v in p));zs=sorted(set(v[1]for v in p))
        if len(xs)!=2 or len(zs)!=2 or any(v[2]!=0 for v in p):raise ValueError('non-affine planar chart')
        if zs[0]<H<zs[1]:raise ValueError('element crosses matching interface')
        corner=[next(v for v in vv if vertices[v][:2]==[x,z])for x,z in ((xs[0],zs[0]),(xs[1],zs[0]),(xs[1],zs[1]),(xs[0],zs[1]))]
        expected=[seen[tuple(sorted((a,b)))]for a,b in zip(corner,corner[1:]+corner[:1])]
        if expected!=ee:raise ValueError('noncanonical quad edge order')
        region='inner'if zs[1]<=H else 'outer'
        elements.append(dict(id=eid,edges=ee,vertices=corner,lower=[xs[0],0.,zs[0]],upper=[xs[1],0.,zs[1]],region=region))
        for e in ee:owners[e].append(eid)
    if sorted(v['id']for v in elements)!=list(range(len(elements))):raise ValueError('noncontiguous element IDs')
    x=sorted(set(v[0]for v in vertices.values()));z=sorted(set(v[1]for v in vertices.values()))
    if len(elements)!=(len(x)-1)*(len(z)-1)or z.count(H)!=1:raise ValueError('incomplete tensor partition/interface')
    roles={k:[]for k in ('inlet','outlet','ground','top','interface','interior')}
    for eid,(a,b)in edges.items():
        p,q=vertices[a],vertices[b];role='interface'if p[1]==q[1]==H else 'ground'if p[1]==q[1]==z[0]else 'top'if p[1]==q[1]==z[-1]else 'inlet'if p[0]==q[0]==x[0]else 'outlet'if p[0]==q[0]==x[-1]else 'interior'
        if len(owners[eid])!=(1 if role in ('inlet','outlet','ground','top')else 2):raise ValueError('unpaired/orphan face')
        roles[role].append(dict(edge=eid,vertices=[a,b],owners=owners[eid]))
    result=dict(vertices={str(k):v for k,v in vertices.items()},edges={str(k):v for k,v in edges.items()},elements=elements,roles=roles,xknots=x,chart_normal_knots=z,
                spaces=dict(differential_degree=3,pressure_degree=1,basis='Modified_A',representation_GLL_points=[8,8]))
    result['sha256']=digest(result);return result

def graph(s,hill,side=0):
    if not hill:return 0.,0.,0.
    k=math.pi/L;at_foot=abs(abs(s)-L)<=32*np.finfo(float).eps*max(1.,L,abs(s))
    if at_foot:return 0.,0.,.5*k*k if(s<0 and side>0)or(s>0 and side<0)else 0.
    if abs(s)>L:return 0.,0.,0.
    return .5*(1+math.cos(k*s)),-.5*k*math.sin(k*s),-.5*k*k*math.cos(k*s)

def mapped(s,rho,hill,*,side=0,region='inner',H=H,d0=D0,top=TOP):
    if not 0<=rho<=top or region not in ('inner','outer'):raise ValueError('outside declared map')
    f,p,pp=graph(s,hill,side);g=math.hypot(1.,p);curvature=-pp/g**3;n=(-p/g,1/g)
    reach=d0+(rho if region=='inner'else H);stretch=1+reach*curvature
    px=s+reach*n[0];pz=f+reach*n[1];xs=stretch;zs=p*stretch;xr,zr=n
    eta=0.
    if region=='outer':
        if rho<H:raise ValueError('outer coordinate below H')
        eta=(rho-H)/(top-H);xr=(s-px)/(top-H);zr=(top-pz)/(top-H)
        px=(1-eta)*px+eta*s;pz=(1-eta)*pz+eta*top;xs=(1-eta)*xs+eta;zs*=(1-eta)
    elif rho>H:raise ValueError('inner rays cannot continue to tunnel ceiling')
    det=xs*zr-xr*zs
    if det<=0 or stretch<=0:raise ValueError('invalid composite determinant/stretch')
    # 2x2 x,z physical Jacobian; a passive unit span has been suppressed.
    return dict(physical=[px,0.,pz],jacobian_xz=[[xs,xr],[zs,zr]],determinant=det,
        cofactor_s=[zr,0.,-xr],cofactor_r=[-zs,0.,xs],normal=[n[0],0.,n[1]],
        kappa_phys=curvature,kappa_eff=curvature/(1+d0*curvature),base_area=g*(1+d0*curvature),
        ray_J=1+rho*curvature/(1+d0*curvature)if region=='inner'else None,
        inner_normal_distance=rho if region=='inner'else None,outer_eta=eta if region=='outer'else None)

def certificate(hill,H=H,d0=D0,top=TOP):
    bound=math.pi**2/(2*L*L)if hill else 0.;slope=math.pi/(2*L)if hill else 0.;amp=1. if hill else 0.
    reach=d0+H;a=1-reach*bound;b=top-amp-reach*math.hypot(1,slope);c=top-amp-reach
    if reach*bound>=.9 or b<=0:raise ValueError('native minimum-stretch/full-height certificate fails')
    return dict(status='PASS',physical_curvature_bound=bound,normal_reach=reach,reach_times_curvature=reach*bound,
        minimum_normal_stretch_lower_bound=a,minimum_required_stretch=.1,outer_ceiling_bracket_lower_bound=b,
        outer_determinant_lower_bound=min(a*b,c)/(top-H),
        global_injectivity='x_s>0; det>0 gives z increasing at fixed physical x; straight flat lateral boundaries give connected eta intervals',
        region_nonoverlap='outer segments stay above the offset-interface graph by eta*(top-f-reach*sqrt(1+M^2))>0; normal tube stays below',
        whole_height=True,normal_ray_scope='only d0+H; outer eta is an independent continuation')

def geometry_checks(t,hill):
    proof=certificate(hill);minimum=(math.inf,None);samples=0;cofactor_error=0.;coordinate_error=0.;face_error=0.;element_metrics=[]
    for el in t['elements']:
        lo,hi=el['lower'],el['upper'];local_min=math.inf
        for a in GLL5:
            s=lo[0]if a==-1 else hi[0]if a==1 else (lo[0]+hi[0])/2+(hi[0]-lo[0])/2*a
            foot=0
            if hill and abs(abs(s)-L)<1e-14:foot=-1 if(lo[0]+hi[0])/2<s else 1
            for b in GLL5:
                rho=lo[2]if b==-1 else hi[2]if b==1 else (lo[2]+hi[2])/2+(hi[2]-lo[2])/2*b
                v=mapped(s,rho,hill,side=foot,region=el['region']);samples+=1;local_min=min(local_min,v['determinant'])
                if v['determinant']<minimum[0]:minimum=(v['determinant'],dict(element=el['id'],chart=[s,0.,rho],foot_side=foot,region=el['region']))
                if rho==TOP:coordinate_error=max(coordinate_error,abs(v['physical'][2]-TOP),abs(v['physical'][0]-s))
                if s in (t['xknots'][0],t['xknots'][-1]):coordinate_error=max(coordinate_error,abs(v['physical'][0]-s))
        element_metrics.append(dict(element=el['id'],region=el['region'],owned_samples=25,minimum_sampled_determinant=local_min))
    # Every shared chart face uses each legitimate owner and opposite orientations.
    by_id={v['id']:v for v in t['elements']}
    vertices={int(k):v for k,v in t['vertices'].items()}
    for role in ('interface','interior'):
        for edge in t['roles'][role]:
            p,q=[vertices[i]for i in edge['vertices']]
            axis='s'if p[0]==q[0]else 'r'
            for u in GLL5:
                s=p[0]if u==-1 else q[0]if u==1 else (p[0]+q[0])/2+(q[0]-p[0])/2*u
                rho=p[1]if u==-1 else q[1]if u==1 else (p[1]+q[1])/2+(q[1]-p[1])/2*u
                rows=[]
                for eid in edge['owners']:
                    el=by_id[eid];foot=0
                    if hill and abs(abs(s)-L)<1e-14:foot=-1 if(el['lower'][0]+el['upper'][0])/2<s else 1
                    v=mapped(s,rho,hill,side=foot,region=el['region']);position=s if axis=='s'else rho
                    sign=-1 if position==el['lower'][0 if axis=='s'else 2]else 1
                    rows.append((v,sign*np.array(v['cofactor_'+axis])))
                face_error=max(face_error,float(max(abs(np.array(rows[0][0]['physical'])-rows[1][0]['physical']))))
                cofactor_error=max(cofactor_error,float(max(abs(rows[0][1]+rows[1][1]))))
    if max(face_error,cofactor_error,coordinate_error)>2e-12:raise ValueError('shared physical face/cofactor/top mismatch')
    feet=[]
    if hill:
        for s in (-L,L):
            feet.append(dict(s=s,owners=[dict(side=side,**mapped(s,H,True,side=side,region='inner'))for side in (-1,1)]))
    return dict(status='PASS',certificate=proof,owned_sample_count=samples,minimum_sampled_determinant=minimum[0],witness=minimum[1],
                physical_shared_face_error=face_error,opposite_cofactor_error=cofactor_error,top_lateral_coordinate_error=coordinate_error,
                foot_limits=feet,element_metrics=element_metrics)

def cg_keys(nx,nz,p,start=0):
    """Enumerate entities without matrices or a native AssemblyMapCG allocation."""
    out=set()
    for j in range(start,start+nz):
        for i in range(nx):
            for r in range(p+1):
                for x in range(p+1):
                    if x<2 and r<2:key=('v',i+x,j+r)
                    elif r<2:key=('h',i,j+r,x)
                    elif x<2:key=('v_edge',i+x,j,r)
                    else:key=('cell',i,j,x,r)
                    out.add(key)
    return out

def dof_counts(t):
    nx=len(t['xknots'])-1;z=t['chart_normal_knots'];ni=z.index(H);no=len(z)-1-ni;n=ni+no
    all3=cg_keys(nx,n,3);inner3=cg_keys(nx,ni,3);outer3=cg_keys(nx,no,3,ni)
    pinner=cg_keys(nx,ni,1);pouter=cg_keys(nx,no,1,ni)
    # Intersecting pressure coordinates belong to different regional fields.
    fields=dict(U=len(all3),W=len(all3),k=len(all3),epsilon=len(outer3),pi_inner=len(pinner),pi_outer=len(pouter),gauge=0)
    line=3*nx+1;inlet=3*n+1;eps_inlet=3*no+1
    essential=dict(U=inlet+line-1,W=inlet+2*line-2,k=inlet,epsilon_excluding_dependent_corner=eps_inlet-1)
    trace=line
    independent=dict(U=fields['U']-essential['U'],W=fields['W']-essential['W'],k=fields['k']-essential['k'],
        epsilon=fields['epsilon']-essential['epsilon_excluding_dependent_corner']-trace,pi_inner=fields['pi_inner'],pi_outer=fields['pi_outer'])
    velocity=independent['U']+independent['W'];pressure=fields['pi_inner']+fields['pi_outer']
    # Picard Epsilon retains dependent trace rows as identity equations and
    # keeps their columns in transport, rather than removing those unknowns.
    solves=dict(momentum_predictor=velocity,pressure=pressure,assembled_flow=velocity+pressure,k=independent['k'],
        epsilon=independent['epsilon']+trace,epsilon_independent_transport=independent['epsilon'],epsilon_identity_trace_rows=trace)
    total=sum(fields.values());formula=3*(3*nx+1)*(3*n+1)+(3*nx+1)*(3*no+1)+(nx+1)*(n+2)
    if total!=formula or len(inner3&outer3)!=trace:raise ValueError('topological count/shared trace formula disagreement')
    return dict(Nx=nx,Ni=ni,No=no,fields=fields,state=total,state_formula=formula,shared_U_k_interface_DOFs=trace,
        regional_velocity=dict(inner=len(inner3),outer=len(outer3)),regional_pressure_interface_vertices_each=nx+1,
        essential=essential,dependent_epsilon_DOFs=trace,independent=independent,free_solve_blocks=solves,
        independent_state=sum(independent.values()),count_kind='source-aware canonical topology; native numbering/sign verification deferred to fresh target CP0')

def storage(count):
    f=count['fields'];b=count['free_solve_blocks'];rows=[]
    def add(name,n,m=None,width=8,copies=1,scope='double payload lower bound only'):
        m=n if m is None else m;rows.append(dict(name=name,rows=n,columns=m,bytes_per_entry=width,copies=copies,bytes=n*m*width*copies,scope=scope))
    for k in ('U','W','k','epsilon','pi_inner','pi_outer'):add('stored_mass_'+k,f[k])
    add('GlobalField_Mass_long_double_temporary',f['U'],width=16,scope='current x86 long-double accumulation, plus double result; report actual ABI')
    add('Picard_flow',b['assembled_flow']);add('Pimple_A',b['momentum_predictor']);add('Pimple_M',b['momentum_predictor'])
    add('Pimple_A_LU',b['momentum_predictor']);add('Pimple_M_LU',b['momentum_predictor'])
    for name in ('G','D','BGp'):add(name,b['momentum_predictor'],b['pressure'])
    add('Schur',b['pressure']);add('scalar_k',b['k']);add('scalar_epsilon_with_trace',b['epsilon'])
    add('legacy_full_state_square_if_requested',count['state'],scope='sizing only, not asserted to be the selected PIMPLE allocation')
    return dict(entries=rows,stored_field_mass_payload_bytes=sum(v['bytes']for v in rows if v['name'].startswith('stored_mass_')),
        simultaneous_Pimple_payload_lower_bound=sum(v['bytes']for v in rows if v['name']in ('Picard_flow','Pimple_A','Pimple_M','Pimple_A_LU','Pimple_M_LU','G','D','BGp','Schur')),
        measured_target_performance=None,note='integer sizing only; no global matrix allocated; excludes row containers, factor workspaces, diagnostics, copies and bound inventories')

def source_capacity(source):
    source=Path(source);base=source/'solvers/IncNavierStokesSolver/EquationSystems';p=(base/'RANSWallEvolutionPatch.cpp').read_text();h=(base/'RANSWallEvolutionPatch.h').read_text();m=(source/'library/SolverUtils/Turbulence/RANSWallNormalMap.cpp').read_text()
    block=p[p.index('if(m_denseCapacity>3000 ||'):p.index('Require(!ColumnAverageK()',p.index('if(m_denseCapacity>3000 ||'))]
    if '!Mapped()'not in block or 'requested<=4331.'not in p or 'm_denseCapacity=3000'not in h or 'return !compositeTunnel && PlanarConstantHeight()'not in m:
        raise ValueError('capacity capability changed: explicit review/admission package required')
    return dict(default_capacity=3000,opt_in_maximum=4331,opt_in_requires_unmapped=True,CompositeTunnelV1_is_mapped=True,
        source=[identity(base/'RANSWallEvolutionPatch.cpp'),identity(base/'RANSWallEvolutionPatch.h'),identity(source/'library/SolverUtils/Turbulence/RANSWallNormalMap.cpp')])

def preflight(xml,manifest,source):
    expected=manifest['flat_session'];actual=identity(xml)
    if actual['sha256']!=expected['sha256']or actual['bytes']!=expected['bytes']:raise ValueError('session identity differs')
    t=topology(xml)
    if t['sha256']!=manifest['shared_topology']['sha256']:raise ValueError('target topology differs')
    policy=source_capacity(source)
    def signature(v):return ({k:x for k,x in v.items()if k!='source'},[(Path(x['path']).name,x['bytes'],x['sha256'])for x in v['source']])
    if signature(policy)!=signature(manifest['capacity_policy']):raise ValueError('capacity source identity changed; requalification required')
    r=ET.parse(xml).getroot();params={v.text.split('=')[0].strip():v.text.split('=')[1].strip()for v in r.findall('./CONDITIONS/PARAMETERS/P')}
    info={v.get('PROPERTY'):v.get('VALUE')for v in r.findall('./CONDITIONS/SOLVERINFO/I')}
    if params['NumSteps']!='0'or info['RANSWallMap']!='Plane'or info['RANSWallInitialVelocity']!='FlatInletTraceV1':raise ValueError('H13 offline scope/initializer differs')
    count=dof_counts(t);blocked=count['state']>policy['default_capacity']
    return dict(status='OBSTRUCTED'if blocked else 'NOT_RUN',native_launch_allowed=False,state=count['state'],effective_capacity=policy['default_capacity'],
        reason='mapped target exceeds legacy dense route; 4331 opt-in excludes CompositeTunnelV1'if blocked else 'H13 is offline-only even if dimensions fit',
        NumSteps_zero_does_not_bypass_capacity=True,all_native_invocations=0)
