#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only CT4 analytical/native-test weak checks and bounded degree contrast."""
import argparse,json,os,shutil,xml.etree.ElementTree as ET
from pathlib import Path
from functools import lru_cache
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
from ct_accepted4_run import expected
from ct_accepted4_reference import Column,reference,write_json,parameters
from ct_accepted_open import probe,require
from ct_accepted2_audit import values
from ct_accepted_pilot import replay,compare
from run_cls_epawf_picard import write,set_parameter
from validate_a93_a95_e94e95 import oracle
from rans_request_manifest import RequestManifest,sha,runtime_identity


def weak(packet,c,wrong_pressure=False,constant_inner_k=False):
    o=packet['observation'];off=o['offsets'];r=np.zeros(len(o['state']));col=Column(c['model'],c['K0'],c['H'],c['roughness_input'],c['nu'],order=96)
    at=lru_cache(None)(lambda z,outer:reference(col,z,outer));p=parameters(c['model'])
    def add(q,b,load):r[np.array(q['ids'],int)+off[b]]+=load
    for q in o['volume']:
        z=q['coordinate'][2];region=bool(q['region']);v=at(z,region);u=np.array([v[0],0,0]);G=np.zeros((3,3));G[0,2]=v[4]
        test=q['velocity_k'];d=np.array(test['derivative']);phi=np.array(test['phi']);w=q['weight'];pi=0. if wrong_pressure else v[3]
        if region:
            law=oracle(c['model'],v[1],v[2],G);T=-law['R']+np.eye(3)*(2*v[1]/3-pi)
            source=law['production']-v[2];diff=law['nut'];eg=q['epsilon'];ed=np.array(eg['derivative']);ephi=np.array(eg['phi'])
            add(eg,3,w*((law['nut']/p['sigma_epsilon'])*v[6]*ed[2]-law['eps_source']*ephi-v[2]*(u@ed)))
        else:
            T=v[8]*(G+G.T)-np.eye(3)*pi;source=float(col.production(z)-v[2]);diff=v[8]
        for b,axis in enumerate((0,2)):add(test,b,w*(T[axis]@d-u[axis]*(u@d)))
        k=col.K if not region and constant_inner_k else v[1];grad=0. if constant_inner_k and not region else v[5]
        add(test,2,w*(diff*grad*d[2]-source*phi-k*(u@d)))
    for f in o['faces']:
        if f['interface'] or f['coordinate'][0]!=c['x'][-1]:continue
        region=bool(f['region']);v=at(f['coordinate'][2],region);q=f['velocity_k'];w=f['weight'];phi=np.array(q['phi']);U=v[0]
        tr=c['reference_column']['outer_equilibrium']['outlet_traction'];tx=tr[0] if region else 0.
        add(q,0,w*phi*(U*U-tx));add(q,1,-w*phi*tr[1]);add(q,2,w*phi*U*(col.K if constant_inner_k and not region else v[1]))
        if region:add(f['epsilon'],3,w*np.array(f['epsilon']['phi'])*U*v[2])
    fixed={int(i) for i,_ in packet['resolution']['constraints']};out={}
    scales=[c['ustar']**2*c['Z']]*2+[c['K0']*c['ustar']*c['Z'],col.A/c['H']*c['ustar']*c['Z'],c['ustar']*c['Z'],c['ustar']*c['Z']]
    for b,name in enumerate(('momentum_x','momentum_z','k_shared_once','outer_epsilon','inner_continuity','outer_continuity')):
        ids=[i for i in range(off[b],off[b+1]) if i not in fixed]
        maximum=float(np.max(abs(r[ids]))) if ids else 0.
        out[name]=dict(free_maximum=maximum,scale=scales[b],scaled_free_maximum=maximum/scales[b],limit=2e-9,
            status='pass' if maximum/scales[b]<=2e-9 else 'fail',free_rows=len(ids),prescribed_rows=len(range(off[b],off[b+1]))-len(ids))
    return dict(status='pass' if all(q['status']=='pass' for q in out.values()) else 'fail',equations=out,
        prescribed_row_reaction_maximum=float(np.max(abs(r[list(fixed)]))),
        scope='independent analytical values/gradients and model sources on native integration/test rows; original conservative outlet terms; interface shared rows counted once; constrained row reactions distinct; not residual cancellation or actual-state equilibrium')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','runtime','initial-root','evolution-root','output'):p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--require-evolved',action='store_true')
    a=p.parse_args();os.environ['LD_LIBRARY_PATH']=str(a.runtime/'lib')
    m=RequestManifest(a.output,dict(schema='CT4ReadOnlySpatialV1',runtime=runtime_identity(a.runtime/'RANSPreservationProbe'),source_sha256=sha(Path(__file__)),
        initial_receipts={str(f):sha(f) for f in a.initial_root.glob('*/results.json')},evolved_receipts={str(f):sha(f) for f in a.evolution_root.glob('*/*/results.json')}))
    result=dict(diagnostic_execution='pass',analytical_weak={},degree_comparison={},actual_alpha_zero={},stopping='no physical continuation authorized by this diagnostic')
    spaces={}
    for degree in (3,4):
        case='A95-6.0-P'+str(degree);folder=a.initial_root/case/'0.002';ini=json.loads((folder/'initial.json').read_text());c=json.loads((folder/'construction.json').read_text())
        q=a.output/('weak-P'+str(degree));q.mkdir();shutil.copy2(folder/'case.xml',q/'case.xml')
        for key,value in zip(('RANSWallNormalOrder','RANSWallTangentialOrder','RANSWallMassOrder','RANSWallVolumeOrder','RANSWallInterfaceOrder'),(64,16,16,16,16),strict=True):
            root=ET.parse(q/'case.xml').getroot();set_parameter(root,key,value);write(root,q/'case.xml')
        packet,e=probe(a.runtime,q,replay(ini,ini['observation']['state'],0.,[],export_basis=True,export_resolution=True),'weak')
        require(e['exit_code']==0,'reference_frozen_weak_probe_failed')
        out=weak(packet,c);out.update(execution=e,raw_packet_sha256=sha(q/'weak.json'),native_actual_alpha_zero_norms=packet['resolution']['norms'])
        result['analytical_weak'][case]=out
        result['analytical_weak'][case]['missing_pressure_offset_negative']=weak(packet,c,wrong_pressure=True)
        result['analytical_weak'][case]['constant_inner_k_negative']=weak(packet,c,constant_inner_k=True)
        require(out['missing_pressure_offset_negative']['status']=='fail' and out['constant_inner_k_negative']['status']=='fail','analytical_native_test_negative_not_discriminating')
        states={};receipts={}
        for dt in ('0.002','0.001'):
            z=json.loads((a.initial_root/case/dt/'case.evolution.json').read_text());states[dt]={0.:z['state']}
            f=a.evolution_root/case/'T0.008/results.json'
            if f.exists():
                receipt=json.loads(f.read_text());member=receipt['members'][dt];d=Path(member['directory']);r=json.loads((d/'case.evolution.json').read_text())
                require(sha(d/'case.evolution.json')==member['native_report_sha256'],'changed_evolved_report')
                for t,s in zip(r['audit_times'],r['accepted_states'],strict=True):states[dt][round(t,12)]=s
                receipts[dt]=member
                if abs(r['time']-.008)<1e-13:
                    dest=a.output/(case+'-evolved-'+dt);dest.mkdir();shutil.copy2(d/'case.xml',dest/'case.xml')
                    pp,ee=probe(a.runtime,dest,replay(ini,r['state'],r['time'],[],export_resolution=True),'actual-alpha-zero')
                    result['actual_alpha_zero'][case+'/'+dt]=dict(execution=ee,norms=pp.get('resolution',{}).get('norms'),packet_sha256=sha(dest/'actual-alpha-zero.json'),time=r['time'],scope='unforced alpha=0 residual of actual evolved coefficients, distinct from stage mass/history residual')
        spaces[degree]=dict(initial=ini,c=c,states=states,receipt=receipts)
    aa,bb=spaces[3],spaces[4];n=json.loads((a.initial_root/'A95-6.0-P3/results.json').read_text())['global_count'];weights=json.loads((a.initial_root/'A95-6.0-P3/results.json').read_text())['weights']
    require(aa['initial']['observation']['grid']==bb['initial']['observation']['grid'] and aa['initial']['observation']['grid_sides']==bb['initial']['observation']['grid_sides'],'degree_comparison_grid_or_side_mismatch')
    common=set.intersection(*(set(s['states'][dt]) for s in spaces.values() for dt in ('0.002','0.001')))
    for t in sorted(x for x in common if x in (0.,.002,.004,.006,.008)):
        per_dt={}
        for dt in ('0.002','0.001'):
            packets=[]
            for s in (aa,bb):
                v=values(s['initial'],s['states'][dt][t],n);v['observation']['time']=t;packets.append(v)
            per_dt[dt]=compare(*packets,weights,aa['c'])
        uncertainty={}
        for field,v in per_dt['0.001']['fields'].items():
            diffs=[]
            for s in (aa,bb):
                p1=values(s['initial'],s['states']['0.002'][t],n);p2=values(s['initial'],s['states']['0.001'][t],n)
                p1['observation']['time']=p2['observation']['time']=t
                diffs.append(compare(p1,p2,weights,s['c'])['fields'][field]['scaled_weighted_RMS'])
            limit=max(1e-7,.1*v['scaled_weighted_RMS']);uncertainty[field]=dict(maximum_dt_difference=max(diffs),limit=limit,status='pass' if max(diffs)<=limit else 'unresolved_temporal_uncertainty')
        result['degree_comparison'][str(t)]=dict(per_dt=per_dt,interpretation_uncertainty=uncertainty,initial_projection=(t==0),
            scope='identical continuous inputs and partition, compatible P3/P1 vs P4/P2; physical weighted RMS and absolute pressure, no convergence order')
    result['analytical_weak_status']='pass' if all(q['status']=='pass' for q in result['analytical_weak'].values()) else 'fail'
    end=result['degree_comparison'].get('0.008')
    result['spatial_sensitivity']='bounded_selected_contrast' if end and all(q['status']=='pass' for q in end['interpretation_uncertainty'].values()) and all(q['status']=='pass' for q in end['per_dt'].values()) else 'unresolved_or_incomplete'
    write_json(a.output/'results.json',result);m.seal(a.output,['results.json']);print('CT4 read-only spatial diagnostics',result['spatial_sensitivity'])
    return 0 if result['analytical_weak_status']=='pass' and (not a.require_evolved or end is not None) else 1
if __name__=='__main__':raise SystemExit(main())
