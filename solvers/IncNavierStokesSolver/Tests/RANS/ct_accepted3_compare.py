#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Common physical six-point observations; bounded contrasts, no order claim."""
import argparse, json, os, hashlib, shutil, xml.etree.ElementTree as ET
from pathlib import Path
from rans_qualification import THREADS, write_json, exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted3_metrics import grid, inventory, point
from ct_accepted3_run import CASES
from ct_accepted3_pressure import projection, projected_normal_load
from ct_accepted_open import probe, require, row_value
from ct_accepted_pilot import replay, compare
from ct_accepted2_audit import values
from ct_accepted2_diagnosis import WITNESS, function_identity, analytic
from ct_accepted2_storage import canonical_location
from rans_request_manifest import RequestManifest, sha, canonical, runtime_identity


def reports(root,case):
    files=[]
    for f in (root/case).rglob('case.evolution.json'):
        if f.parent.name not in ('0.002','0.001'):continue
        files.append(f)
    return sorted(files)


def read_space(root,case,ct2=False):
    if ct2:
        zeros={dt:root/'companion-initial-v1/work'/('Zminus-'+dt) for dt in ('0.002','0.001')}
        files=[root/group/'work'/('Zminus-'+dt)/'case.evolution.json' for group in ('companion-v1','continue02-v1','continue2-v1') for dt in zeros]
    else:
        candidates=list((root/case).rglob('initial.json'))
        zeros={dt:next(f.parent for f in candidates if f.parent.name==dt) for dt in ('0.002','0.001')}
        files=reports(root,case)
    c=json.loads((zeros['0.002']/'construction.json').read_text());states={};sources={};initials={}
    for dt,folder in zeros.items():
        packet=json.loads((folder/'initial.json').read_text());initials[dt]=packet
        states[dt]={0.:packet['observation']['state']};sources[dt]={0.:dict(file=str(folder/'initial.json'),sha256=sha(folder/'initial.json'))}
    for f in files:
        dt=f.parent.name.removeprefix('Zminus-');r=json.loads(f.read_text());digest=sha(f)
        require(abs(float(dt)-r['dt'])<1e-15,'foreign_report_timestep')
        require(len(r['audit_times'])==len(r['accepted_states']),'accepted_state_inventory_length')
        for t,state in zip(r['audit_times'],r['accepted_states'],strict=True):
            t=round(t,12)
            if t in states[dt]:require(states[dt][t]==state,'duplicate_restart_boundary_state_differs')
            states[dt][t]=state;sources[dt][t]=dict(file=str(f),sha256=digest)
    require(exact_vector_hash(states['0.002'][0.])==exact_vector_hash(states['0.001'][0.]),'same_space_comparison_initial_mismatch')
    return dict(construction=c,zeros=zeros,initials=initials,states=states,sources=sources)


def normal_load(initial,state,c):
    o=initial['observation'];off=o['offsets'];s=np.array(state);terms=np.zeros(4);elements={}
    L=c['x'][-1]-c['x'][0]
    for q in o['volume']:
        if not q['region']:continue
        v=point(q,s,off,c['model'])
        if v['status']!='admissible':return dict(status='undefined_original_constitutive_domain',witness=v)
        x,_,z=q['coordinate'];factor=np.exp(-((z-6.25)/1.5)**2)
        phi=factor*x*(L-x);dx=factor*(L-2*x);dz=-2*(z-6.25)/1.5**2*phi
        stress=-np.array(v['R'])+np.eye(3)*2*v['k']/3
        u=row_value(q['velocity_k'],s,off[0])[0];w=row_value(q['velocity_k'],s,off[1])[0]
        p=-v['Pi']*dz;visc=stress[2,0]*dx+stress[2,2]*dz;adv=-w*(u*dx+w*dz)
        nonlinear=-np.array(v['Q'])[2,0]*dx-np.array(v['Q'])[2,2]*dz
        load=q['weight']*np.array([p,visc,adv,nonlinear]);terms+=load
        gid=str(q['generating_geometry']);elements.setdefault(gid,np.zeros(4));elements[gid]+=load
    return dict(status='pass',pressure=float(terms[0]),full_constitutive_stress=float(terms[1]),conservative_advection=float(terms[2]),
        nonlinear_subset=float(terms[3]),sum=float(terms[:3].sum()),elements={k:v.tolist() for k,v in elements.items()},
        fixed_test='x*(L-x)*exp(-((z-6.25)/1.5)^2)',scope='outer volume load against unchanged physical test; interface boundary term not included; not a steady residual or qualification by cancellation')


def sample_summary(basis,state,c,n,weights):
    sampled=values(basis,state,len(basis['test_basis']));o=sampled['observation'];s=np.array(state)
    fixed=inventory(basis['test_basis'],o['grid'],s,basis['observation']['offsets'],c['model'],o['grid_sides'])
    g=np.array(o['grid_gradients'])[:n];w=np.array(weights);v=np.array(o['grid_values'])[:n];div=g[:,0,0]+g[:,1,2]
    outer=np.array(o['grid_regions'][:n])==1;iso=-2*v[outer,2]*div[outer]/3
    named=[]
    for i,p in enumerate(o['grid']):
        if p==WITNESS or p==[0.,0.,6.109375] or (p[0] in (0.,.25) and p[2]==c['H']):
            q=basis['test_basis'][i];entry=dict(coordinate=p,side=o['grid_sides'][i],geometry=q['velocity_k']['geometry'])
            if q['epsilon'] is not None:entry['invariants']=point(q,s,basis['observation']['offsets'],c['model'])
            else:entry['values']=o['grid_values'][i];entry['gradients']=o['grid_gradients'][i]
            named.append(entry)
    return sampled,dict(domain=fixed,named_points=named,divergence_weighted_RMS=float(np.sqrt(w@(div*div)/w.sum())),
        divergence_maximum=float(np.max(abs(div))),isotropic_work_RMS=float(np.sqrt(w[outer]@(iso*iso)/w[outer].sum())),
        isotropic_work_maximum=float(np.max(abs(iso))))


def radicands(sampled):
    o=sampled['observation'];mask=np.array(o['grid_regions'])==1
    gradients=np.array(o['grid_gradients'])[mask];v=np.array(o['grid_values'])[mask]
    G=np.zeros((len(v),3,3));G[:,0,:]=gradients[:,0,:];G[:,2,:]=gradients[:,1,:]
    div=np.trace(G,axis1=1,axis2=2);S=(G+G.transpose(0,2,1))/2-div[:,None,None]*np.eye(3)/3
    O=(G-G.transpose(0,2,1))/2;sn=np.sqrt(np.sum(S*S,axis=(1,2)));on=np.sqrt(np.sum(O*O,axis=(1,2)))
    eig=np.linalg.eigvalsh(S);tau=v[:,2]/v[:,3];C=np.empty(len(v));regular=sn!=0
    C[regular]=1/(6.5+3*eig[regular,-1]/sn[regular]*np.hypot(sn[regular],on[regular])*tau[regular])
    C[~regular]=1/6.5
    require(not np.any((~regular)&(on!=0)),'undefined_pure_rotation_in_readonly_rad_comparison')
    return 1-9*C*C*sn*sn*tau*tau, np.flatnonzero(mask)


def spatial(a,b,ba,bb,ca,cb,points,weights,n,t):
    av,sa=sample_summary(ba,a,ca,n,weights);bv,sb=sample_summary(bb,b,cb,n,weights)
    local={}
    if ca['model']=='A95' and (sa['domain']['domain_failure_count'] or sb['domain']['domain_failure_count']):
        local=dict(status='failed_actual_domain_on_common_grid',reason='original scalar/constitutive evaluation undefined',radicand_comparison='unavailable')
    elif ca['model']=='A95':
        ar,ix=radicands(av);br,bx=radicands(bv);require(np.array_equal(ix,bx),'local_side_inventory_mismatch')
        j=int(np.argmax(abs(ar-br)));difference=float(abs(ar[j]-br[j]));index=int(ix[j]);margin=float(min(ar.min(),br.min()))
        local=dict(maximum_absolute_radicand_difference=difference,minimum_margin_both_spaces=margin,
            difference_over_margin=difference/margin if margin>0 else None,coordinate=points[index],side=ba['observation']['grid_sides'][index],
            status='failed_actual_domain_on_common_grid' if margin<0 else 'below_evaluation_floor' if difference<=1e-10 else 'bounded_on_declared_grid' if margin>0 and difference<=.25*margin else 'unresolved_spatial_domain_sensitivity',
            declared_limit_fraction=.25,evaluation_absolute_floor=1e-10)
    for q in (av,bv):
        q['observation']['time']=t
        for key in ('grid','grid_sides','grid_regions','grid_values','grid_gradients'):q['observation'][key]=q['observation'][key][:n]
    fields=compare(av,bv,weights,ca)
    # Global field differences and actual local constitutive sensitivity are
    # separate: no positive global RMS can excuse a negative local radicand.
    return dict(time=t,fields=fields,local=local,coarse=sa,fine=sb,
                classification='bounded_global_contrast' if fields['status']=='pass' else 'failed_global_difference_limit',
                convergence_order='not_claimed')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('root','ct2','output'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();runtime=a.ct2/'runtime';os.environ['LD_LIBRARY_PATH']=str(runtime/'lib')
    manifest=RequestManifest(a.output,dict(schema='CT3CommonPhysicalComparisonsV1',runtime=runtime_identity(runtime/'RANSPreservationProbe'),
        source_sha256=sha(Path(__file__)),contract_sha256=sha(a.root/'source/docs/plans/cao-tamura-accepted-3-evidence/contract.json'),
        primary_receipts={str(f.relative_to(a.root)):sha(f) for f in a.root.glob('**/results.json') if '/launch-' not in str(f)},
        CT2_runtime_bridge_sha256=sha(a.ct2/'relocation-runtime-bridge.json'),quadrature=6,pressure='absolute, no shift'))
    spaces={'A95-Zminus-P3-CT2':read_space(a.ct2,'',True)}
    for case in CASES:
        if list((a.root/case).rglob('initial.json')):spaces[case]=read_space(a.root,case)
    result=dict(schema='CT3BoundedSpatialEvidenceV1',diagnostic_execution='pass',spaces={},contrasts={},reuse={})
    expected=json.loads((a.ct2/'runtime-identity.json').read_text())['runtimes']
    require(canonical_location(runtime_identity(runtime/'IncNavierStokesSolver'))==canonical_location(expected['IncNavierStokesSolver']), 'CT2_runtime_reuse_identity_mismatch')
    result['reuse']=dict(status='pass',runtime=runtime_identity(runtime/'IncNavierStokesSolver'),
        source='sealed CT2 Zminus, .008/.02/.2; old histories evaluated read-only; no evolution rerun',
        bridge_sha256=sha(a.ct2/'relocation-runtime-bridge.json'))
    bases={};samples={}
    for name,space in spaces.items():
        c=space['construction'];points,sides,weights,n=grid(c)
        folder=a.output/name;folder.mkdir();shutil.copy2(space['zeros']['0.002']/'case.xml',folder/'case.xml')
        ini=space['initials']['0.002'];req=replay(ini,space['states']['0.002'][0.],0.,points,export_basis=True,export_test_basis=True)
        req['grid_sides']=sides
        packet,e=probe(runtime,folder,req,'common-basis');require(e['exit_code']==0,'common_physical_basis_export_failed')
        bases[name]=packet;samples[name]={};row=dict(construction=c,common_packet_sha256=sha(folder/'common-basis.json'),execution=e,
            initial_native_hashes={dt:exact_vector_hash(sp[0.]) for dt,sp in space['states'].items()},states={},temporal=[],projection={})
        pr=projection(ET.parse(folder/'case.xml').getroot(),c,c['pressure_degree'],packet['observation'],packet['test_basis'][-1])
        lo,hi=[projected_normal_load(ET.parse(folder/'case.xml').getroot(),c,pr,k) for k in (48,96)]
        row['projection']=dict(pr,independent_normal_load48=lo,independent_normal_load96=hi,quadrature_difference=abs(lo-hi))
        common=set(space['states']['0.002'])&set(space['states']['0.001'])
        selected=sorted(t for t in common if t in (0.,.002,.004,.006,.008,.02,.2))
        for dt in ('0.002','0.001'):
            samples[name][dt]={}
            for t in selected:
                state=space['states'][dt][t];v,summary=sample_summary(packet,state,c,n,weights)
                samples[name][dt][t]=v
                row['states'][dt+'@'+str(t)]=dict(time=t,state_sha256=exact_vector_hash(state),source=space['sources'][dt][t],
                    observations=summary,normal_load=normal_load(packet,state,c))
        for t in selected:
            ps=[]
            for dt in ('0.002','0.001'):
                q=samples[name][dt][t];o=q['observation'];small={k:v[:n] if isinstance(v,list) else v for k,v in o.items()};small['time']=t;ps.append(dict(observation=small))
            comparison=compare(*ps,weights,c)
            if c['model']=='A95':
                valid=all(not row['states'][dt+'@'+str(t)]['observations']['domain']['domain_failure_count'] for dt in ('0.002','0.001'))
                if valid:
                    ar,ix=radicands(samples[name]['0.002'][t]);br,bx=radicands(samples[name]['0.001'][t])
                    require(np.array_equal(ix,bx),'temporal_local_side_inventory_mismatch')
                    j=int(np.argmax(abs(ar-br)));comparison['maximum_radicand_dt_difference']=float(abs(ar[j]-br[j]));comparison['radicand_dt_location']=points[int(ix[j])]
                else:comparison['radicand_comparison']='unavailable_original_domain'
            row['temporal'].append(comparison)
        bad=[dict(sample=k,time=v['time'],domain=v['observations']['domain']) for k,v in row['states'].items()
             if v['observations']['domain']['domain_failure_count'] or v['observations']['domain']['minima'].get('covariance',{}).get('value',-1)<-1e-10*c['K0']]
        row['common_grid_domain']=dict(status='pass' if not bad else 'fail',first_failure=bad[0] if bad else None,
            scope='additional common physical inventory; historical CT2 receipt unchanged')
        result['spaces'][name]=row;write_json(a.output/'progress.json',result)
    for label,coarse,fine in [('streamwise_P3','A95-Zminus-P3-CT2','A95-XZminus-P3'),('degree_same_partition','A95-XZminus-P3','A95-XZminus-P4')]:
        if coarse not in spaces or fine not in spaces:continue
        ca,cb=spaces[coarse]['construction'],spaces[fine]['construction'];points,sides,weights,n=grid(ca)
        require(bases[coarse]['observation']['grid']==bases[fine]['observation']['grid'] and bases[coarse]['observation']['grid_sides']==bases[fine]['observation']['grid_sides'],'cross_space_physical_samples_differ')
        common=set(spaces[coarse]['states']['0.002'])&set(spaces[fine]['states']['0.002'])&set(spaces[coarse]['states']['0.001'])&set(spaces[fine]['states']['0.001'])
        times=sorted(t for t in common if t in (0.,.002,.004,.006,.008,.02,.2));rows={}
        for dt in ('0.002','0.001'):
            rows[dt]=[spatial(spaces[coarse]['states'][dt][t],spaces[fine]['states'][dt][t],bases[coarse],bases[fine],ca,cb,points,weights,n,t) for t in times]
        for q in rows['0.001']:
            t=q['time'];ta=next(x for x in result['spaces'][coarse]['temporal'] if x['time']==t);tb=next(x for x in result['spaces'][fine]['temporal'] if x['time']==t)
            uncertainty={}
            for field,d in q['fields']['fields'].items():
                dtmax=max(x['fields'][field]['scaled_weighted_RMS'] for x in (ta,tb));limit=max(1e-7,.1*d['scaled_weighted_RMS'])
                uncertainty[field]=dict(maximum_dt_difference=dtmax,limit=limit,status='pass' if dtmax<=limit else 'unresolved_temporal_uncertainty')
            localdt=max(ta.get('maximum_radicand_dt_difference',0.),tb.get('maximum_radicand_dt_difference',0.))
            if 'maximum_absolute_radicand_difference' in q['local']:
                q['local']['maximum_dt_radicand_difference']=localdt
                q['local']['dt_uncertainty_limit']=max(1e-10,.1*q['local']['maximum_absolute_radicand_difference'])
                if localdt>q['local']['dt_uncertainty_limit'] and q['local']['status']!='failed_actual_domain_on_common_grid':q['local']['status']='unresolved_temporal_uncertainty'
            q['field_spatial_interpretation_uncertainty']=uncertainty
        evolved=[q for q in rows['0.001'] if q['time']>0]
        result['contrasts'][label]=dict(coarse=coarse,fine=fine,common_times=times,comparisons=rows,
            status='bounded_on_declared_grid' if evolved and all(q['fields']['status']=='pass' and q['local'].get('status')=='bounded_on_declared_grid'
                and all(x['status']=='pass' for x in q['field_spatial_interpretation_uncertainty'].values()) for q in evolved) else 'unresolved_spatial_sensitivity',
            initial_projection='t=0 reported separately; later differences include evolution; no formal order inferred')
        write_json(a.output/'progress.json',result)
    result['evolved_spatial_status']='bounded_contrasts_only' if result['contrasts'] and all(x['status']=='bounded_on_declared_grid' for x in result['contrasts'].values()) else 'unresolved'
    write_json(a.output/'results.json',result);manifest.seal(a.output,['results.json']);print(result['evolved_spatial_status']);return 0
if __name__=='__main__':raise SystemExit(main())
