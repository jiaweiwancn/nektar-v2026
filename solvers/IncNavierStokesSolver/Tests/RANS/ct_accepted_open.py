#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""CT-Accepted-1 case definitions and independent conservative outlet checks."""
import argparse,json,math,os,resource,subprocess,time,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.linalg import block_diag,cholesky,solve_triangular,svdvals
from run_a93_a95_e94e95_matrix import case_xml
from run_cls_epawf_picard import set_info,set_parameter,write
from cls_epawf_m4_cases import tensor_mesh
from ct_accepted_feasibility import Column,parameters
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,runtime_identity,sha


def require(ok,reason):
    if not ok:raise ValueError(reason)


def functions(root,name,values):
    c=root.find('CONDITIONS')
    for f in list(c.findall('FUNCTION')):
        if f.get('NAME')==name:c.remove(f)
    f=ET.SubElement(c,'FUNCTION',{'NAME':name})
    for key,value in values.items():ET.SubElement(f,'E',{'VAR':key,'VALUE':str(value)})


def fixture(folder,runtime,model='A95',steps=0,dt=1e-7,pressure=0.,traction_shift=0.,flux_shift=0.,cache='Off',backflow=False,restart=None,reject=0):
    path=case_xml(folder,model,'Natural','Accepted',steps=steps,dt=dt,descriptor_probe=runtime/'RANSNewModelsProbe',
        radial_knots=[0.,2.,4.],degree=3,xknots=[0.,.5,1.],cache=cache,restart=restart,reject=reject,
        quadrature=[32,8,10,10,10])
    root=ET.parse(path).getroot();p=parameters(model);K=1.;H=2.;nu=1.5e-5
    c=Column(model,K,H,0.,nu);eps=c.A/H
    for key,value in dict(RANSWallBoundaryPolicy='AcceptedOpenFlatV1',RANSOpenTop='DirichletCompatibleV1',
        RANSOpenBackflow='RejectV1',RANSWallReference='Functions',RANSWallPimpleSnapshots='On').items():set_info(root,key,value)
    for key,value in dict(RANSWallRoughness=0.,K0=K,RANSWallFrictionVelocity=math.sqrt(K/p['q']),RANSOpenVelocityScale=1.).items():set_parameter(root,key,value)
    u=('-' if backflow else '')+'.02*y*(8-y)'
    values=dict(u=u,v='0',k='1',eps_tilde=repr(eps),p=repr(pressure))
    for suffix in ('Inner','Outer'):
        functions(root,'InitialConditions'+suffix,values)
        functions(root,'BoundaryValues'+suffix,values)
        functions(root,'OutletTotalTraction'+suffix,dict(u=repr(-pressure-traction_shift),v='.01'))
        functions(root,'OutletDiffusiveFlux'+suffix,dict(k=repr(.003+flux_shift),eps_tilde=repr(.004+2*flux_shift)))
    write(root,path);return path


def run(command,folder,log):
    start=time.monotonic();before=resource.getrusage(resource.RUSAGE_CHILDREN)
    with (folder/log).open('x') as f:r=subprocess.run([str(x) for x in command],cwd=folder,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,**THREADS))
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    return dict(command=[str(x) for x in command],exit_code=r.returncode,wall_seconds=time.monotonic()-start,
        process_cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,maximum_child_RSS_KiB=after.ru_maxrss)


def probe(runtime,folder,request,name='sample'):
    write_json(folder/(name+'-request.json'),request)
    execution=run([runtime/'RANSPreservationProbe',folder/'case.xml',folder/(name+'-request.json'),folder/(name+'.json')],folder,name+'.log')
    packet=json.loads((folder/(name+'.json')).read_text());write_json(folder/(name+'-execution.json'),execution)
    return packet,execution


def sampling_request():
    return dict(operation='functions',grid=[[.25,0,.3],[.75,0,2.],[.75,0,2.],[.75,0,3.]],grid_sides=[0,-1,1,0],
        export_basis=True,export_constitutive=True,export_open_verification=True,export_mixed_space=True,
        export_tensor_validation=True)


def row_value(row,state,offset):
    ids=np.asarray(row['ids'],int)+offset;phi=np.asarray(row['phi']);d=np.asarray(row['derivative'])
    return phi@state[ids],d@state[ids]


def independent_operator(packet,model='A95'):
    """Published tensors and conservative weak form on exported *native* rows.

    This fixture has zero roughness so no roughness-volume pressure forcing.
    Geometry/basis is reused; stress, sources, signs and tensor assembly are not.
    """
    o=packet['observation'];s=np.array(o['state']);off=np.array(o['offsets'],int);n=len(s)
    masses=[np.array(m) for m in packet['open_frozen_systems']['masses']]
    M=block_diag(*masses,np.zeros((1,1)));A=M.copy();r=np.zeros(n)
    params=parameters(model);c=Column(model,1.,2.,0.,1.5e-5)
    for point in o['volume']:
        region=int(point['region']);v=point['velocity_k'];pr=point['pressure'];ids=np.array(v['ids'],int)
        phi=np.array(v['phi']);d=np.array(v['derivative'])[[0,2]];w=point['weight']
        u=np.array([row_value(v,s,off[b])[0] for b in (0,1)]);G=np.zeros((3,3))
        for b,axis in enumerate((0,2)):G[axis]=row_value(v,s,off[b])[1]
        k,_=row_value(v,s,off[2]);pi,_=row_value(pr,s,off[4+region]);nt=c.D(point['coordinate'][2])-c.nu
        if region:
            eps,_=row_value(point['epsilon'],s,off[3]);law=oracle(model,k,eps,G)
            nt=law['nut'];T=-law['R']+np.eye(3)*(2*k/3-pi);visc=nt
            rates={2:(eps-min(0.,law['production']))/k,3:(1.9 if model=='E94E95' else 1.92)*eps/k}
            sources={2:max(0.,law['production']),3:law['eps_production']}
            diffusion={2:nt,3:nt/params['sigma_epsilon']}
        else:
            T=(c.nu+nt)*(G+G.T)-np.eye(3)*pi;visc=c.nu+nt
            P=nt*np.sum((G+G.T)*G)
            rates={2:0.};sources={2:P-c.A/max(point['coordinate'][2],c.cutoff)};diffusion={2:c.nu+nt}
        adv=d.T@u
        for a,axis in enumerate((0,2)):
            ii=ids+off[a];r[ii]+=w*(d.T@T[axis,[0,2]]-adv*u[a])
            for b,other in enumerate((0,2)):
                jj=ids+off[b];block=visc*np.outer(d[b],d[a])
                if a==b:block+=visc*d.T@d-np.outer(adv,phi)
                A[np.ix_(ii,jj)]+=w*block
            pp=np.array(pr['ids'],int)+off[4+region];pphi=np.array(pr['phi'])
            A[np.ix_(ii,pp)]-=w*np.outer(d[a],pphi);A[np.ix_(pp,ii)]+=w*np.outer(pphi,d[a])
        pp=np.array(pr['ids'],int)+off[4+region];r[pp]+=w*np.array(pr['phi'])*np.trace(G)
        for b in ([2,3] if region else [2]):
            q=point['epsilon'] if b==3 else v;ix=np.array(q['ids'],int)+off[b];ph=np.array(q['phi']);der=np.array(q['derivative'])[[0,2]]
            f,grad=row_value(q,s,off[b]);ad=der.T@u
            r[ix]+=w*(ph*(rates[b]*f-sources[b])-ad*f+diffusion[b]*(der.T@grad[[0,2]]))
            A[np.ix_(ix,ix)]+=w*(rates[b]*np.outer(ph,ph)-np.outer(ad,ph)+diffusion[b]*der.T@der)
    boundary=dict(advection=np.zeros(n),traction=np.zeros(n),diffusion=np.zeros(n))
    for face in packet['open_boundary']['faces']:
        v=face['velocity_k_row'];u=np.array(face['velocity']);area=np.array(face['area'])[[0,2]];un=u@area;w=face['weight'];length=np.linalg.norm(area)
        for b in ([0,1,2,3] if face['region'] else [0,1,2]):
            q=face['epsilon_row'] if b==3 else v;ids=np.array(q['ids'],int)+off[b];ph=np.array(q['phi'])
            f=face['velocity'][b] if b<2 else face['k'] if b==2 else face['epsilon']
            boundary['advection'][ids]+=w*ph*un*f;A[np.ix_(ids,ids)]+=w*un*np.outer(ph,ph)
            if b<2:boundary['traction'][ids]-=w*ph*length*face['traction'][b]
            else:boundary['diffusion'][ids]-=w*ph*length*face['k_diffusive_datum' if b==2 else 'epsilon_diffusive_datum']
    for x in boundary.values():r+=x
    return r,A,boundary


def classify(packet):
    require(packet['status']=='pass','native_actual_state_rejected');o=packet['observation'];f=packet['open_frozen_systems']
    r,A,loads=independent_operator(packet);native=np.array(f['raw_residual']);scale=max(1.,np.max(abs(native)))
    defects=dict(raw_weak=float(np.max(abs(r-native)))/scale)
    require(defects['raw_weak']<=2e-10,'independent_conservative_weak_load')
    for name,load in loads.items():
        defect=float(np.max(abs(load-np.array(packet['open_boundary'][name+'_load' if name!='diffusion' else 'diffusive_load']))))
        defects[name+'_outlet']=defect;require(defect<=2e-10,'outlet_'+name+'_sign')
    off=np.array(o['offsets'],int);s=np.array(o['state']);external=dict(f['external_constraints'])
    Hrows=set()
    for face in o['faces']:
        if face['interface'] and face['region']==1:Hrows.update(off[3]+int(i) for i in face['epsilon']['ids'] if abs(face['epsilon']['phi'][face['epsilon']['ids'].index(i)])>1e-14)
    # The dependent-epsilon rows are replaced by frozen trace equations, not
    # an additional inner equation or two-sided epsilon-flux condition.
    for block in f['blocks']:
        ids=np.array(block['rows'],int);a=A[np.ix_(ids,ids)].copy()
        if block['block']==3:
            for j,row in enumerate(ids):
                if row in Hrows:a[j]=0.;a[j,j]=1.
        d=float(np.max(abs(a-np.array(block['matrix']))));defects['matrix_'+str(block['block'])]=d
        require(d<=2e-10,'independent_frozen_matrix_'+str(block['block']))
        d=float(np.max(abs(np.array(block['residual'])-np.array(f['constrained_residual'])[ids])))
        defects['frozen_residual_'+str(block['block'])]=d;require(d<=2e-10,'frozen_original_residual_consistency')
    mixed=packet['mixed_space'];rows=np.array(mixed['rows'],int);nv=int(mixed['velocity_size']);a=np.array(mixed['matrix']);np_=len(rows)-nv
    masses=[np.array(m) for m in f['masses']];mass=block_diag(*masses,np.ones((1,1)))[np.ix_(rows,rows)]
    L=cholesky(mass,lower=True);white=solve_triangular(L,a,lower=True);white=solve_triangular(L,white.T,lower=True).T
    G=a[:nv,nv:];D=a[nv:,:nv];C=a[nv:,nv:];schur=D@np.linalg.solve(mass[:nv,:nv],G)-C
    P=cholesky(mass[nv:,nv:],lower=True);sw=solve_triangular(P,schur,lower=True);sw=solve_triangular(P,sw.T,lower=True).T
    spectra={k:svdvals(v) for k,v in dict(mixed=white,pressure_action=G,divergence=D,pimple_pressure=sw).items()}
    ranks={k:int(np.count_nonzero(v>1e-11*v[0])) for k,v in spectra.items()}
    require(ranks==dict(mixed=len(rows),pressure_action=np_,divergence=np_,pimple_pressure=np_),'open_pressure_extra_null_mode')
    constant=np.array(mixed['pressure_constant_with_zero_gauge']);level=float(np.max(abs(a@constant)))
    require(level>1e-5,'absolute_traction_did_not_determine_pressure_level')
    require(off[6] not in rows and not any(off[4]<=int(i)<off[6] for i in external),'physical_pressure_constrained_or_unused_gauge_in_solve')
    return dict(status='pass',unknowns=len(s),physical_unknowns=len(s)-1,defects=defects,ranks=ranks,
        singular_values={k:v.tolist() for k,v in spectra.items()},constant_pressure_action=level,
        rank_relative_threshold=1e-11,pressure='absolute traction; all physical tests free; unused zero slot excluded',
        scope='actual conservative production residual and frozen matrices, finite two-column fixture; no general inf-sup theorem')


def verify(runtime,output):
    manifest=RequestManifest(output,dict(schema='CTAcceptedOpenVerificationV1',runtime={n:runtime_identity(runtime/n) for n in ('IncNavierStokesSolver','RANSPreservationProbe','RANSNewModelsProbe')},
        fixture='two columns, p3/p1, H=2 Z=4, smooth mathematical source with positive actual constant k',
        limits=dict(weak=2e-10,matrix=2e-10,pressure_offset=2e-10,rank_relative=1e-11)))
    result=dict(status='fail',checks={},executions=[])
    try:
        packets={}
        for name,options in [('base',{}),('data_signs',dict(traction_shift=.125,flux_shift=.007)),('pressure_offset',dict(pressure=.125)),('backflow',dict(backflow=True)),('unsupported_cache',dict(cache='ExactStepV1'))]:
            folder=output/name;fixture(folder,runtime,**options);p,e=probe(runtime,folder,sampling_request());packets[name]=p;result['executions'].append(e)
            if name in ('backflow','unsupported_cache'):
                reason='outlet_backflow_unsupported' if name=='backflow' else 'evaluation caching is not qualified'
                require(e['exit_code']!=0 and reason in json.dumps(p),'negative_control_wrong_reason:'+name)
                result['checks'][name]=dict(status='pass',expected_rejection=reason);continue
            result['checks'][name]=classify(p)
        a=packets['base'];b=packets['data_signs'];c=packets['pressure_offset']
        require(a['observation']['state']==b['observation']['state'],'prescribed_data_changed_initial_coefficients')
        for name,p in [('data_signs',b),('pressure_offset',c)]:
            free=sorted(set(range(len(a['observation']['state'])))-{int(i) for i,v in a['open_frozen_systems']['external_constraints']})
            delta=np.array(p['open_frozen_systems']['raw_residual'])-a['open_frozen_systems']['raw_residual']
            if name=='data_signs':
                expected=np.array(p['open_boundary']['traction_load'])-a['open_boundary']['traction_load']+np.array(p['open_boundary']['diffusive_load'])-a['open_boundary']['diffusive_load']
                defect=float(np.max(abs(delta-expected)));require(defect<=2e-10,'native_outlet_data_signs_not_consumed')
            else:
                defect=float(np.max(abs(delta[free])));require(defect<=2e-10,'pressure_offset_velocity_or_scalar_weak_action')
            result['checks'][name+'_discriminating']=dict(status='pass',maximum_defect=defect)
        result['status']='pass'
    except Exception as error:result['failure']=repr(error)
    write_json(output/'results.json',result)
    if result['status']=='pass':manifest.seal(output,['results.json']+[str(p.relative_to(output)) for p in output.rglob('case.xml')]+[str(p.relative_to(output)) for p in output.rglob('sample.json')])
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--runtime',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib');r=verify(a.runtime,a.output);print(json.dumps(dict(status=r['status'],failure=r.get('failure'))));return 0 if r['status']=='pass' else 1
if __name__=='__main__':raise SystemExit(main())
