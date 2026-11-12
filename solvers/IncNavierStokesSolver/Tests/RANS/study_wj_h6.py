#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H6 retained-state replay, one fixed mesh, and finite common-point gates."""
import argparse
from pathlib import Path
import subprocess
import numpy as np
from numpy.testing import assert_allclose
from run_wj_h6 import load,save,sha
from study_wj_h5 import configurations,identity,verify,assess
from wj_h5_cases import case_xml,mesh_knots
from wj_h5_neutral import KAPPA,QE,mean_fit,measurements,SCALE_MM
from wj_h6_projection import reconstruct,Functions,gll
from wj_h4_points import owners,evaluate,physical,mapping,vector_hash
from wj_h3_analysis import principal
from validate_wj1 import Probe


PARENT='wj-h5-equilibrium-inlet-20261112'
def configuration(height):return next(c for c in configurations()if abs(c.H*SCALE_MM-height)<1e-12)
def record(root,height,coarse=False):
    folder=(root.parent/PARENT if coarse else root)/'cases'/f'H{height:g}mm'
    path=folder/'case.evolution.json';return load(path),path,folder/'case.xml'


def import_h5(root):
    verify(root)
    source=Path(__file__).resolve().parents[4];manifest=load(source/'docs/plans/wj-h5-linux/manifest.json')
    expected={v['path']:v for v in manifest['heavy']};records=[]
    for height in (5.,5.5,6.):
        r,p,xml=record(root,height,True)
        for path in (p,xml):assert identity(path)==expected[str(path)]
        assert r['steps']==0 and r['time']==0 and not r['pimple_attempts'] and r['implicit_callback_count']==0
        records.append(dict(height_mm=height,producer=identity(p),xml=identity(xml),failure=r['failure'],
            role='uncommitted startup candidate'if r['failure']else 'accepted CP0',state_sha256=vector_hash(r['state'])))
    fit=mean_fit();bounds=[]
    for row in fit['points']:
        low=(row['r_mm']-row['reading_height_mm'])*np.exp(-KAPPA*(row['value']+row['reading_value'])/fit['u_tau'])
        high=(row['r_mm']+row['reading_height_mm'])*np.exp(-KAPPA*(row['value']-row['reading_value'])/fit['u_tau'])
        bounds.append(dict(point=row['point'],compatible_intercept_mm=[float(low),float(high)]))
    save(root/'imports.json',dict(parent_manifest=identity(source/'docs/plans/wj-h5-linux/manifest.json'),records=records,
        refitted_estimator_reading_bounds_mm=fit['reading_bounds_mm']['combined'],individual_point_bounds=bounds,
        point_compatible_intersection_mm=[max(v['compatible_intercept_mm'][0]for v in bounds),min(v['compatible_intercept_mm'][1]for v in bounds)],
        scope='two different deterministic reading constructions, neither a confidence interval; fit/window unchanged'))


def generate(root):
    intent=[]
    for c in configurations():
        height=c.H*SCALE_MM;old,old_report,old_xml=record(root,height,True)
        path=case_xml(root/'cases'/f'H{height:g}mm',c,bisect_once=True)
        # The original XML, not regenerated expected values, binds the reference.
        import xml.etree.ElementTree as ET
        funcs=lambda p:[ET.tostring(v)for v in ET.parse(p).findall('./CONDITIONS/FUNCTION')]
        assert funcs(path)==funcs(old_xml),'reference Functions must remain byte-identical'
        design=load(path.with_name('case-design.json'));assert design['predicted_state_size']==2328
        intent.append(dict(height_mm=height,xml=identity(path),design=identity(path.with_name('case-design.json')),
            parent_xml=identity(old_xml),expected_dimensions=[16,34,2328],governing_timesteps=0))
    save(root/'run-intent.json',intent)
    print('Exactly three fresh zero-step one-bisection inputs frozen',flush=True)


def packets(r):return r['wj_audit_packets']+[v['wj_packet']for v in r.get('physical_range_events',[])if v.get('wj_packet')]


def inventory(reports,c,dense=False):
    points={}
    def add(p,tag):
        p=list(map(float,p))
        if not(0<=p[0]<=8 and c.d0<=p[2]<=c.top):return
        key=tuple(v.hex()for v in p);points.setdefault(key,dict(physical=p,tags=[]))['tags'].append(tag)
    for x in (0.,4.,8.):
        for z in (c.d0,c.d0+c.H,7/SCALE_MM):add([x,0,z],'ground/H/Fig12')
        for v in measurements():add([x,0,v['z_mm']/SCALE_MM],'Fig3')
    for r in reports:
        for packet in packets(r):
            for name,w in packet['principal_witnesses'].items():
                add(w['position'],'native_'+name)
                for axis in (0,2):
                    for sign in (-1,1):
                        q=w['position'].copy();q[axis]+=sign*1e-6;add(q,'witness_physical_offset')
        for row in r['open_boundary']['interface']:add(row['physical'],'native_interface')
        if dense:
            for row in r['open_boundary']['boundary_rows']:add(row['physical'],'native_exterior')
            for e in r['fields'][0]['mesh']:
                z=np.r_[-1.,np.polynomial.legendre.leggauss(12)[0],gll(8)[0],1.]
                for rr in (e['lower'][2]+e['upper'][2])/2+(e['upper'][2]-e['lower'][2])/2*z:
                    for x in (0.,4.,8.):add(physical(r,[x,0,rr]),'face/Gauss12/LGL8')
    from wj_h4_points import digest
    v=dict(points=list(points.values()),dense=dense,scope='finite physical union, all owners; exported native witnesses plus declared quadrature/face samples')
    v['sha256']=digest(v);return v


def replay(root,height,coarse=False,version='v2'):
    if not version.replace('-','').isalnum():raise ValueError('simple immutable replay version required')
    c=configuration(height);r,rp,xml=record(root,height,coarse);out=root/('replays-'+version)/('coarse'if coarse else 'refined')/f'H{height:g}mm';out.mkdir(parents=True)
    before=identity(rp);stages,meta=reconstruct(r,xml)
    functions=Functions(xml,r);g=np.polynomial.legendre.leggauss(16)[0];xml_errors=dict(U=0.,k=0.,epsilon=0.)
    for e in r['fields'][0]['mesh']:
        rr=(e['lower'][2]+e['upper'][2])/2+(e['upper'][2]-e['lower'][2])/2*g;x=np.full_like(rr,4.)
        zz=np.array([physical(r,[4.,0,v])[2]-c.d0 for v in rr]);inside=e['upper'][2]<=c.H
        expected_u=np.array([c.U_inner(v)if inside else c.U_outer(v)for v in zz])
        expected_k=np.array([c.k_inner(v)if inside else c.K for v in zz])
        xml_errors['U']=max(xml_errors['U'],float(np.max(abs(functions.block(0,x,rr)-expected_u))))
        xml_errors['k']=max(xml_errors['k'],float(np.max(abs(functions.block(2,x,rr)-expected_k))))
        if not inside:xml_errors['epsilon']=max(xml_errors['epsilon'],float(np.max(abs(functions.block(3,x,rr)-c.A/zz))))
    assert max(xml_errors.values())<2e-11,'XML approximation differs from the frozen continuum'
    save(out/'reconstruction.json',dict(meta,producer=before,xml=identity(xml),XML_continuum_max_errors=xml_errors))
    np.savez_compressed(out/'reconstructed-stages.npz',**stages)
    if meta['status']!='PASS':print('UNRESOLVED final coefficient replay; attribution withheld',flush=True);return
    inv=inventory([r],c);rows=[];grad_errors=[];field_errors=[];witness=[]
    for point in inv['points']:
        p=point['physical']
        for owner in owners(r,p):
            native=evaluate(r,r['state'],p,owner);final=evaluate(r,stages['dependent_epsilon_lift'],p,owner)
            grad_errors.append(float(np.max(abs(np.asarray(native['G'])-final['G']))))
            field_errors.extend([float(np.max(abs(np.asarray(native['U'])-final['U']))),
                abs(native['k']-final['k'])/c.K,abs(native['pi']-final['pi'])/c.u_tau**2,
                float(np.max(abs(np.asarray(native['R'])-final['R'])))/c.K])
            if native['epsilon'] is not None:field_errors.append(abs(native['epsilon']-final['epsilon'])/c.epsilon(c.H))
            if native['principal']is not None:
                assert max(abs(native['principal'][key]-final['principal'][key])for key in ('in_plane','transverse'))<2e-11
            for name,s in dict(stages,native=np.asarray(r['state'])).items():
                v=evaluate(r,s,p,owner);rows.append(dict(stage=name,**v))
    assert max(grad_errors)<2e-11 and max(field_errors)<2e-11,'final field/gradient replay differs; no attribution allowed'
    kernel=Probe(str(root.parent/'wj-h2-20261111/build-v1/WJKernelProbe'))
    try:
        for packet in packets(r):
            for name,w in packet['principal_witnesses'].items():
                p=w['position'];o=next(o for o in owners(r,p)if o['geometry']==w['geometry'])
                v=evaluate(r,r['state'],p,o);G=np.asarray(v['G']);rphys=p[2]-c.d0;Gc=np.zeros((3,3));Gc[0,2]=c.u_tau/(KAPPA*rphys)
                H=np.asarray(w['H_row_major']).reshape(3,3);val=kernel.closure(v['k'],v['epsilon'],G);action=kernel.closure(v['k'],v['epsilon'],G,(0,0,H))
                compiled=-np.sum(H*action['R'])/val['scalar'][1];assert_allclose(compiled,w['minimum_over_nu_eff'],rtol=2e-9,atol=2e-11)
                vs=evaluate(r,stages['physical_L2'],p,o)
                witness.append(dict(kind=name,native=w,compiled_replay=float(compiled),
                    stages={label:next(vv for vv in rows if vv['stage']==label and vv['owner']['geometry']==o['geometry']and vv['physical']==list(map(float,p)))for label in (*stages,'native')},
                    diagnostic_hybrids=dict(exact_G_native_scalars=principal(v['k'],v['epsilon'],Gc),native_G_exact_scalars=principal(c.K,c.epsilon(rphys),G),exact=principal(c.K,c.epsilon(rphys),Gc)),
                    shear_error=v['G'][0][2]/Gc[0,2]-1,epsilon_error=v['epsilon']/c.epsilon(rphys)-1,
                    outlet_Gxx_identity=dict(actual=v['G'][0][0],predicted=(vs['U'][0]-c.U_outer(rphys))/8,
                        scope='witness on outlet/normal face; no state replacement')))
    finally:kernel.close()
    assert identity(rp)==before
    save(out/'stage-rows.json',rows);save(out/'witnesses.json',witness)
    save(out/'result.json',dict(status='PASS',reconstruction=identity(out/'reconstruction.json'),stages=identity(out/'reconstructed-stages.npz'),
        inventory=inv,field_replay_max=max(field_errors),physical_G_replay_max=max(grad_errors),witnesses=witness,
        stage_label='passive reconstruction only',replay_version=version,producer=before))
    print(f'PASS: {height:g}mm '+('H5'if coarse else 'H6')+' final coefficients/fields/gradients replay and compiled witnesses',flush=True)


def weak_status(executed,value):
    if not executed:return 'NOT_RUN'
    if value not in ('PASS','UNRESOLVED','OBSTRUCTED','FAIL'):raise ValueError('explicit executed weak status required')
    return value


def compare(root,height):
    c=configuration(height);a,ap,ax=record(root,height,True);b,bp,bx=record(root,height)
    inv=inventory([a,b],c,True);out=root/'comparisons'/f'H{height:g}mm';out.mkdir(parents=True)
    save(out/'inventory.json',inv);summary={};allrows={}
    for label,r in (('H5',a),('H6',b)):
        errors=dict(U=0.,G=0.,Gxz=0.,k=0.,epsilon=0.,covariance=0.,shear_parameter=0.)
        rows=[];min_in=float('inf');min_trans=float('inf');min_cov=float('inf')
        for point in inv['points']:
            p=point['physical'];distance=p[2]-c.d0
            for o in owners(r,p):
                if o['region']!='outer':continue
                v=evaluate(r,r['state'],p,o);g=c.u_tau/(KAPPA*distance);ep=c.epsilon(distance)
                exactG=np.zeros((3,3));exactG[0,2]=g
                local=dict(U=abs(v['U'][0]-c.U_outer(distance)),G=float(np.max(abs(np.asarray(v['G'])-exactG))),
                    Gxz=abs(v['G'][0][2]/g-1),k=abs(v['k']/c.K-1),epsilon=abs(v['epsilon']/ep-1),
                    covariance=float(np.max(abs(np.asarray(v['R'])-c.Rc)))/c.K,
                    shear_parameter=abs(v['k']/v['epsilon']*v['G'][0][2]/QE-1))
                for key in errors:errors[key]=max(errors[key],local[key])
                min_in=min(min_in,v['principal']['in_plane']);min_trans=min(min_trans,v['principal']['transverse']);min_cov=min(min_cov,min(v['covariance']))
                rows.append(dict(**v,reference_errors=local))
        exact=principal(c.K,c.epsilon(c.H),np.array([[0,0,c.u_tau/KAPPA/c.H],[0,0,0],[0,0,0.]]))['in_plane']
        summary[label]=dict(errors=errors,min_in_plane=min_in,min_transverse=min_trans,min_covariance=min_cov,
            exact_min=exact,retained_fraction=min_in/exact,outer_owner_rows=len(rows))
        save(out/(label+'-rows.json'),rows);allrows[label]=identity(out/(label+'-rows.json'))
    gates=dict(principal_half_margin=summary['H6']['retained_fraction']>=.5,
        transverse=summary['H6']['min_transverse']>0,covariance=summary['H6']['min_covariance']>0,
        gradient_reduction=summary['H6']['errors']['Gxz']<summary['H5']['errors']['Gxz'],
        epsilon_reduction=summary['H6']['errors']['epsilon']<summary['H5']['errors']['epsilon'])
    save(out/'result.json',dict(status='PASS'if all(gates.values())else 'OBSTRUCTED',gates=gates,inventory=identity(out/'inventory.json'),
        summaries=summary,rows=allrows,producer_H5=identity(ap),producer_H6=identity(bp),
        scope='same physical point union/all legitimate owners; one bisection, no order or global-stability claim'))
    print(height,summary,gates,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--mode',choices=('import','generate','replay','observe','compare'),required=True)
    p.add_argument('--height-mm',type=float);p.add_argument('--coarse',action='store_true');p.add_argument('--observer',type=Path);p.add_argument('--replay-version',default='v2')
    a=p.parse_args()
    if a.mode=='import':import_h5(a.root)
    elif a.mode=='generate':generate(a.root)
    elif a.mode=='replay':replay(a.root,a.height_mm,a.coarse,a.replay_version)
    elif a.mode=='observe':assess(a.root,a.height_mm,a.observer)
    else:compare(a.root,a.height_mm)
