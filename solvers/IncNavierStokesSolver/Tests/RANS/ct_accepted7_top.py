#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""RN2011 physical top: independent all-row/frozen/rank and live history tests.

The d=0 controls are explicitly nonexperimental. No displaced physical case
is substituted by these controls or admitted by their successful execution.
"""
import argparse,copy,json,math,os,xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS,write_json,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from scipy.linalg import block_diag
from ct_accepted_feasibility import Column,parameters
from ct_accepted4_run import make as legacy_make
from ct_accepted4_analysis import weak
from ct_accepted3_pressure import pressure_action
from ct_accepted_open import probe,run,require,row_value,independent_operator
from ct_accepted_pilot import replay,budget_classification
from ct_accepted_budget_gate import physical_budgets
from run_cls_epawf_picard import set_info,set_parameter,write
from rans_request_manifest import RequestManifest,runtime_identity,sha
from validate_a93_a95_e94e95 import oracle


def make(folder,rt,model='A95',degree=3,steps=0,dt=1e-6,restart=None,level=0.,old_top=False,interface_order=32):
    p=parameters(model);us=.4;H=.5;nu=1.5e-5;hr=.4;col=Column(model,p['q']*us*us,H,hr,nu,order=96)
    ref=col.result();G=np.zeros((3,3));G[0,2]=us/(p['kappa']*H);law=oracle(model,col.K,col.A/H,G)
    ref['outer_equilibrium']=dict(outlet_traction=[float(law['Q'][2,2]-law['Q'][0,0]),us*us])
    c=dict(model=model,H=H,H_mm=20.,Z=1.,nu=nu,ustar=us,B=col.t0*col.J0-us/p['kappa']*math.log(H),K0=col.K,
        roughness_input=hr,reference_column=ref,x=[0.,.5,1.],knots=[0.,.0625,.5,1.],degree=degree,pressure_degree=degree-2,
        provider='NeutralCompositeV1',fixture='CT7RNSyntheticFixtureV1',scope='synthetic d0 RN boundary verification, not experimental inlet admission')
    xml=legacy_make(folder,rt,c,dt,steps,parent=restart);root=ET.parse(xml).getroot()
    if not old_top:
        set_info(root,'RANSWallBoundaryPolicy','AcceptedRichardsNorrisFlatV1')
        set_info(root,'RANSOpenTop','RichardsNorris2011V1')
    if level:set_parameter(root,'RANSOpenPressureLevel',level)
    set_parameter(root,'RANSWallMassOrder',12);set_parameter(root,'RANSWallVolumeOrder',12);set_parameter(root,'RANSWallInterfaceOrder',interface_order)
    c['integration']=dict(normal=32,tangential=8,mass=12,volume=12,interface=interface_order)
    write(root,xml);write_json(folder/'construction.json',dict(c,configuration_sha256=sha(xml),dt=dt,steps=steps))
    return c


def request(c):
    return dict(operation='functions',grid=[[.25,0,.2],[.75,0,.5],[.75,0,.5],[.75,0,1.]],grid_sides=[0,-1,1,0],
        export_basis=True,export_resolution=True,export_constitutive=True,export_tensor_validation=True,
        export_open_verification=True,export_mixed_space=True)


def classify(packet,c):
    require(packet['status']=='pass','RN actual initialized state rejected')
    receipt=packet['initialization_receipt']
    require(receipt['open_explicit_implicit_agreement'] and receipt['explicit_boundary_sha256']==receipt['implicit_boundary_sha256'],
        'RN explicit/implicit live essential boundary mismatch')
    o=packet['observation'];s=np.asarray(o['state']);off=o['offsets'];f=packet['open_frozen_systems'];top=packet['open_boundary']['top_faces']
    require(bool(top),'no consumed RN top rows')
    col=Column(c['model'],c['K0'],c['H'],c['roughness_input'],c['nu'],order=96)
    r,A,loads=independent_operator(packet,c['model'],column=col)
    # Tangential roughness convection and grad(K_H) vanish for these radial
    # control states. Check that premise rather than suppressing actual data.
    radial_defect=0.
    for face in o['faces']:
        if not face['interface']:continue
        for b in (0,1,2):radial_defect=max(radial_defect,abs(row_value(face['velocity_k'],s,off[b])[1][0]))
    require(radial_defect<2e-12,'radial roughness oracle used on tangentially varying data')
    # Pressure derivatives at the actual rough-volume x and outer H owner,
    # independently reconstructed from P1/P2 polynomials, not native dphi.
    for j in range(off[4],off[6]):
        direction=np.zeros(len(s));direction[j]=1.
        action,rough,_=pressure_action(o,direction,c['roughness_input'],c['pressure_degree'])
        A[:off[2],j]=action[:off[2]]
        r+=s[j]*rough
    native=np.asarray(f['raw_residual']);defects=dict(all_row=float(np.max(abs(r-native)))/max(1.,float(np.max(abs(native)))))
    for name,load in loads.items():
        key='diffusive_load' if name=='diffusion' else name+'_load'
        defects['boundary_'+name]=float(np.max(abs(load-np.asarray(packet['open_boundary'][key]))))
    fixed={int(i) for i,_ in packet['resolution']['constraints']}
    for b in f['blocks']:
        ids=np.array(b['rows'],int);expected=A[np.ix_(ids,ids)].copy()
        if b['block']==3:
            for j,row in enumerate(ids):
                if row in fixed:expected[j]=0.;expected[j,j]=1.
        defects['frozen_'+str(b['block'])]=float(np.max(abs(expected-np.asarray(b['matrix']))))
        defects['residual_'+str(b['block'])]=float(np.max(abs(np.asarray(b['residual'])-np.asarray(f['constrained_residual'])[ids])))
    require(max(defects.values())<=2e-10,'RN independent all-row/frozen operator mismatch:'+str(defects))
    external={int(i) for i,_ in f['external_constraints']};mode_ids={b:set() for b in (0,1,2,3)};e_load=np.zeros(len(s));tx_load=np.zeros(len(s));length=0.
    for face in top:
        w=face['weight']*float(np.linalg.norm(face['area']));length+=w
        for b in (0,1,2,3):
            q=face['epsilon_row'] if b==3 else face['velocity_k_row']
            active=np.abs(q['phi'])>1e-13;ids=np.asarray(q['ids'],int)[active]+off[b];mode_ids[b].update(ids.tolist())
        v=face['velocity_k_row'];e=face['epsilon_row']
        tx_load[np.asarray(v['ids'],int)+off[0]]-=w*np.asarray(v['phi'])*face['tangential_total_traction']
        e_load[np.asarray(e['ids'],int)+off[3]]-=w*np.asarray(e['phi'])*face['epsilon_diffusive_datum']
        require(face['epsilon_diffusive_datum']<0 and face['k_diffusive_datum']==0 and not face['independent_normal_traction'],'RN flux/normal-traction contract')
    free_top={str(b):sorted(mode_ids[b]-external) for b in mode_ids}
    require(all(free_top[str(b)] for b in (0,2,3)) and not free_top['1'],'top essential/free modes wrong')
    q=weak(packet,c);maximum=max(v['scaled_free_maximum'] for v in q['equations'].values())
    require(maximum<=2e-10,'continuous RN reference balance on native tests:'+str(q))
    # Deleting either physical datum leaves a resolved load in actual free
    # top modes. All trace modes, including zero integral odd modes, are kept.
    free=sorted(set(range(len(s)))-fixed)
    eps_signal=float(np.max(abs(e_load[free])));shear_signal=float(np.max(abs(2*tx_load[free])))
    require(eps_signal>1e-8 and shear_signal>1e-8,'nondiscriminating top negative controls')
    omitted=copy.deepcopy(packet)
    for f in omitted['open_boundary']['top_faces']:f['epsilon_diffusive_datum']=0.
    zero_flux_balance=weak(omitted,c)
    require(zero_flux_balance['equations']['outer_epsilon']['status']=='fail',
        'zero epsilon top flux did not break the continuous reference identity')
    mixed=packet['mixed_space'];rows=np.asarray(mixed['rows'],int);nv=int(mixed['velocity_size']);mat=np.asarray(mixed['matrix'])
    mass=block_diag(*[np.asarray(m) for m in f['masses']],np.ones((1,1)))[np.ix_(rows,rows)]
    schur=mat[nv:,:nv]@np.linalg.solve(mass[:nv,:nv],mat[:nv,nv:])
    row_scale=np.maximum(np.linalg.norm(schur,axis=1),np.finfo(float).tiny);normalized=schur/row_scale[:,None]
    normalized/=np.maximum(np.linalg.norm(normalized,axis=0),np.finfo(float).tiny)[None,:]
    singular=np.linalg.svd(normalized,compute_uv=False);ratio=float(singular[-1]/singular[0])
    action=float(np.linalg.norm(mat@np.asarray(mixed['pressure_constant_with_zero_gauge'])))
    require(ratio>1e-10 and action>1e-8 and not any(off[4]<=i<off[6] for i in external),'RN absolute-pressure rank/constraints')
    require(off[6] not in rows,'unused gauge included in physical pressure solve')
    for field in packet['tensor_validation']['fields']:
        require(max(field['defects_value_dx_dz_dxz_dzz_mass_polynomial_faces'])<=2e-10,'native all-mode basis/face checks')
    return dict(status='pass',actual_stored_slots=len(s),defects=defects,continuous_native_test_maximum=maximum,
        top_trace_mode_ids={str(b):sorted(v) for b,v in mode_ids.items()},free_top_mode_ids=free_top,
        top_physical_length=length,prescribed_epsilon_outward_conservative_flux=col.um**4/(parameters(c['model'])['sigma_epsilon']*c['Z']),
        zero_epsilon_flux_negative_signal=eps_signal,reversed_shear_negative_signal=shear_signal,
        zero_epsilon_flux_continuous_negative=zero_flux_balance['equations']['outer_epsilon'],
        initialization_receipt=receipt,
        pressure_minimum_singular_ratio=ratio,constant_pressure_action=action,additional_pressure_null_modes=int(np.sum(singular<=singular[0]*1e-10)),
        scope='complete actual governing operator/frozen matrices and finite mixed-pressure fixture; no formal inf-sup or physical preservation claim')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('runtime','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--interface-order',type=int,choices=(12,24,32),default=32)
    a=p.parse_args();m=RequestManifest(a.output,dict(schema='CT7RichardsNorrisVerificationV1',runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source=sha(Path(__file__)),interface_order=a.interface_order,controls='all3 P3 and A95 P4; same synthetic d0 compatible reference; absolute traction'))
    result=dict(schema='CT7RichardsNorrisVerificationV1',checks={},execution='incomplete',physical_experimental_admission=False)
    try:
        for model,degree in [('A95',3),('A93',3),('E94E95',3),('A95',4)]:
            label=model+'-p'+str(degree);folder=a.output/label;c=make(folder,a.runtime,model,degree,interface_order=a.interface_order)
            try:
                packet,e=probe(a.runtime,folder,request(c),'native');require(e['exit_code']==0,'RN native construction:'+str(packet.get('failure')))
                result['checks'][label]=dict(classify(packet,c),execution=e,packet_sha256=sha(folder/'native.json'),construction_sha256=sha(folder/'construction.json'))
            except Exception as error:result['checks'][label]=dict(status='fail',failure=repr(error),directory=str(folder))
            write_json(a.output/'progress.json',result);print('CT7 RN top',label,result['checks'][label]['status'],flush=True)
        require(result['checks']['A95-p3']['status']=='pass','base A95 operator unavailable for dependent offset/history controls')
        base=a.output/'A95-p3';initial=json.loads((base/'native.json').read_text());off=initial['observation']['offsets']
        folder=a.output/'pressure-level';c=make(folder,a.runtime,level=.125,interface_order=a.interface_order)
        shifted,e=probe(a.runtime,folder,request(c),'native');require(e['exit_code']==0,'shifted pressure-level native rejection')
        delta=np.asarray(shifted['observation']['state'])-initial['observation']['state'];expected=np.zeros(len(delta))
        for q in initial['observation']['pressure_mass']:
            row=q['pressure'];expected[np.asarray(row['ids'],int)+off[4+q['region']]]=.125*np.asarray(row['constant'])
        defect=float(np.max(abs(delta-expected)));free=sorted(set(range(len(delta)))-{int(i) for i,_ in initial['resolution']['constraints']})
        change=np.asarray(shifted['open_frozen_systems']['raw_residual'])-initial['open_frozen_systems']['raw_residual']
        require(max(defect,float(np.max(abs(change[free]))))<=2e-10,'absolute pressure/traction shift changed physical fields/free weak action')
        result['checks']['pressure_offset']=dict(status='pass',coefficient_defect=defect,free_residual_defect=float(np.max(abs(change[free]))),execution=e)
        # A top K mode not touching inlet or interface changes actual nu_t,
        # while the prescribed epsilon physical datum must remain immutable.
        state=np.asarray(initial['observation']['state']).copy();fixed={int(i) for i,_ in initial['resolution']['constraints']}
        ids=sorted({off[2]+i for f in initial['open_boundary']['top_faces'] for i,v in zip(f['velocity_k_row']['ids'],f['velocity_k_row']['phi'],strict=True) if abs(v)>1e-12}-fixed)
        state[ids[-1]]+=1e-3*c['K0'];req=replay(initial,state.tolist(),0.,[],export_basis=True)
        changed,e=probe(a.runtime,base,req,'changed-top-k');require(e['exit_code']==0,'admissible top coefficient control rejected')
        old=[f['epsilon_diffusive_datum'] for f in initial['open_boundary']['top_faces']];new=[f['epsilon_diffusive_datum'] for f in changed['open_boundary']['top_faces']]
        require(old==new,'prescribed epsilon flux changed with evolving state')
        wrong=[]
        for face in changed['open_boundary']['top_faces']:
            v=face['velocity_k_row'];er=face['epsilon_row'];G=np.zeros((3,3));s=np.asarray(changed['observation']['state'])
            for b,axis in enumerate((0,2)):G[axis]=row_value(v,s,off[b])[1]
            k=row_value(v,s,off[2])[0];eps=row_value(er,s,off[3])[0];law=oracle('A95',k,eps,G)
            A=c['ustar']**3/parameters('A95')['kappa']
            wrong.append(-A/c['Z']**2*law['nut']/parameters('A95')['sigma_epsilon'])
        signal=max(abs(x-y) for x,y in zip(wrong,new,strict=True))
        require(signal>1e-8,'evolving-diffusivity negative control is not discriminating')
        result['checks']['state_independent_epsilon_flux']=dict(status='pass',changed_global_K_row=ids[-1],actual_state_sha256=exact_vector_hash(changed['observation']['state']),datum=old[0],wrong_evolving_diffusivity_signal=signal,execution=e)
        for name,modify,word in [('displaced_original_wall',lambda root:set_parameter(root,'RANSOpenDisplacement',.0875),'displacement'),
             ('unsupported_cache',lambda root:set_info(root,'RANSWallEvaluationReuse','ExactStepV1'),'caching')]:
            folder=a.output/name;c=make(folder,a.runtime,interface_order=a.interface_order);root=ET.parse(folder/'case.xml').getroot();modify(root);write(root,folder/'case.xml')
            packet,e=probe(a.runtime,folder,dict(operation='functions'),'negative');require(e['exit_code']!=0 and word in packet.get('failure','').lower(),'wrong expected rejection:'+name)
            result['checks'][name]=dict(status='pass',reason=packet['failure'],execution=e)
        result['execution']='complete';result['status']='pass' if all(r['status']=='pass' for r in result['checks'].values()) else 'fail'
    except Exception as error:result['failure']=repr(error);result['status']='fail'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json'])
    return 0 if result['status']=='pass' else 1

if __name__=='__main__':raise SystemExit(main())
