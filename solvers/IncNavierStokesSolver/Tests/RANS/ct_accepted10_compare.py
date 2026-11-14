#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Common-coordinate CT10 projection, drift, pressure and spatial comparisons."""
import argparse,os
from functools import lru_cache
from pathlib import Path
from rans_qualification import THREADS,exact_vector_hash
os.environ.update(THREADS)
import numpy as np
import ct_accepted8_native as n
from ct_accepted9_pressure_run import sealed
from ct_accepted9_gate import pair
from ct_accepted2_audit import values
from ct_accepted3_metrics import inventory
from ct_accepted5_traction import profiles,difference,diffusivity
from ct_accepted_feasibility import write_json,require
from ct_accepted_open import row_value
from ct_accepted3_metrics import point
from rans_request_manifest import RequestManifest,sha


def regional_pressure(a,b,w,regions,scale):
    rows={};d=a-b
    for label,mask in [('whole',np.ones(len(w),bool)),('inner',regions==0),('outer',regions==1)]:
        x=d[mask];wt=w[mask];mean=float(wt@x/wt.sum());rows[label]=dict(signed_mean=mean,weighted_RMS=float(np.sqrt(wt@(x*x)/wt.sum())),spatial_variation_RMS=float(np.sqrt(wt@((x-mean)**2)/wt.sum())),scale=scale)
    return rows


def gradients(p,s,c,weights,count):
    v=values(p,s,count)['observation'];g=np.asarray(v['grid_gradients']);w=np.asarray(weights);r=np.asarray(v['grid_regions']);expected=np.zeros_like(g);col=n.col(c)
    at=lru_cache(None)(lambda z,outer:n.reference(col,z,bool(outer)))
    for i,(q,outer) in enumerate(zip(v['grid'],r,strict=True)):
        f=at(q[2],outer);expected[i,0,2]=f[4];expected[i,2,2]=f[5]
        if outer:expected[i,3,2]=f[6]
    result={};scales=[c['ustar']/c['H'],c['ustar']/c['H'],c['K0']/c['H'],col.A/(c['H']-c['displacement'])/c['H'],c['ustar']**2/c['H']]
    for j,name in enumerate(('U','w','k','outer_epsilon','Pi')):
        mask=r==1 if j==3 else np.ones(count,bool);ds=g[mask,j]-expected[mask,j];ww=w[mask];square=np.sum(ds*ds,axis=1);rms=float(np.sqrt(ww@square/ww.sum()));ids=np.flatnonzero(mask);i=int(ids[np.argmax(square)])
        result[name]=dict(weighted_vector_RMS=rms,scale=scales[j],scaled=rms/scales[j],maximum_vector_norm=float(np.sqrt(max(square))),witness=dict(coordinate=v['grid'][i],side=v['grid_sides'][i],geometry=p['test_basis'][i]['velocity_k']['geometry']))
    div=g[:,0,0]+g[:,1,2];result['divergence']=dict(weighted_RMS=float(np.sqrt(w@(div*div)/w.sum())),maximum=float(max(abs(div))),signed_mean=float(w@div/w.sum()))
    return result


def outlet(p,s,c):
    s=np.asarray(s);off=p['observation']['offsets'];col=n.col(c);ref=n.reference(col,c['H'],True);G=np.zeros((3,3));G[0,2]=ref[4]
    from validate_a93_a95_e94e95 import oracle
    Q=oracle('A95',col.K,ref[2],G)['Q'];expected=[np.array([0.,col.t0]),np.array([Q[2,2]-Q[0,0],col.t0])];square=den=0.;witness=None
    for f in p['observation']['faces']:
        if f['interface'] or f['coordinate'][0]!=c['x'][-1] or f['area'][0]<=0:continue
        q=f['velocity_k'];grad=np.zeros((3,3))
        for b,axis in enumerate((0,2)):grad[axis]=row_value(q,s,off[b])[1]
        if f['region']:
            law=point(f,s,off,'A95');require('total_stress' in law,'undefined outlet constitutive evaluation');T=np.array(law['total_stress'])
        else:
            x=f['coordinate'][0];Hrow=next(p['test_basis'][i]['velocity_k'] for i,(v,side) in enumerate(zip(p['observation']['grid'],p['observation']['grid_sides'],strict=True)) if v==[x,0.,c['H']] and side==1);K=row_value(Hrow,s,off[2])[0];Pi=row_value(f['pressure'],s,off[4])[0];T=diffusivity(K,f['coordinate'][2],c)*(grad+grad.T)-np.eye(3)*Pi
        d=T[[0,2],0]-expected[f['region']];value=float(d@d);square+=f['weight']*value;den+=f['weight']
        if witness is None or value>witness['square']:witness=dict(square=value,coordinate=f['coordinate'],geometry=q['geometry'],region=f['region'],actual=T[[0,2],0].tolist(),prescribed=expected[f['region']].tolist())
    return dict(weighted_RMS=float(np.sqrt(square/den)),normalized_weighted_RMS=float(np.sqrt(square/den))/c['ustar']**2,witness=witness,scope='actual reconstructed natural-face stress versus datum, not a prescribed-row reaction or pointwise outlet admission requirement')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('preflight','baseline','candidate','output'):p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();pre=sealed(a.preflight);base=sealed(a.baseline);fine=sealed(a.candidate);w=np.asarray(pre['weights']);count=pre['global_count'];c=pre['cases']['baseline'];n.reference=lru_cache(None)(n.reference)
    m=RequestManifest(a.output,dict(schema='CT10CommonSpaceComparisonV1',preflight=sha(a.preflight/'results.json'),baseline=sha(a.baseline/'results.json'),candidate=sha(a.candidate/'results.json'),source=sha(Path(__file__)),rule='same physical grid,weights and explicit sides; separate actual initial states; absolute pressure; no two-space order inference'))
    out=dict(schema='CT10CommonSpaceComparisonV1',baseline_qualification=pair(base),candidate_qualification=pair(fine),initial_gradients={},times={},spatial_convergence='two representations, no order claim',physical_target_admission=False)
    ip={};reports={}
    for label,q in [('baseline',base),('candidate',fine)]:
        folder=Path(pre['initial']['baseline' if label=='baseline' else 'candidate-0.001']['directory']);ip[label]=n.read(folder/'initial.json');out['initial_gradients'][label]=gradients(ip[label],ip[label]['observation']['state'],c,w,count)
        reports[label]={}
        for dt,row in q['members'].items():
            file=Path(row['directory'])/'case.evolution.json';require(sha(file)==row['native_report_sha256'],'changed compared native state');r=n.read(file);reports[label][dt]={round(t,12):s for t,s in zip(r['audit_times'],r['accepted_states'],strict=True)}
    for key in ('grid','grid_sides','grid_regions'):require(ip['baseline']['observation'][key]==ip['candidate']['observation'][key],'cross-space physical sample identity mismatch:'+key)
    regions=np.asarray(ip['baseline']['observation']['grid_regions'][:count]);all_times=sorted(set.intersection(*(set(v) for q in reports.values() for v in q.values())))
    require(all_times and all_times[0]==0. and all_times[-1]==.008,'missing complete equal-time comparison')
    for t in all_times:
        rows={};temporal={}
        for dt in ('0.001','0.0005'):
            x,y=[n.observations(ip[k],reports[k][dt][t],count) for k in ('baseline','candidate')]
            rows[dt]=dict(fields=n.norms(x,y,w,regions,c),pressure=regional_pressure(x[:,4],y[:,4],w,regions,c['ustar']**2),traction=difference(profiles(ip['baseline'],reports['baseline'][dt][t],c),profiles(ip['candidate'],reports['candidate'][dt][t],c),c),states={k:exact_vector_hash(reports[k][dt][t]) for k in reports})
        for label in reports:
            pp=[n.observations(ip[label],reports[label][dt][t],count)[:,4] for dt in ('0.001','0.0005')];temporal[label]=regional_pressure(pp[0],pp[1],w,regions,c['ustar']**2)
        rows['temporal_pressure']=temporal
        if t in (0.,.001,.002,.004,.008):
            rows['state_checks']={}
            for label in reports:
                state=np.asarray(reports[label]['0.0005'][t]);ini=ip[label]
                rows['state_checks'][label]=dict(gradients=gradients(ini,state,c,w,count),constitutive=inventory(ini['test_basis'],ini['observation']['grid'],state,ini['observation']['offsets'],'A95',ini['observation']['grid_sides']),outlet=outlet(ini,state,c))
        out['times'][str(t)]=rows
    out['status']='comparison_completed';write_json(a.output/'results.json',out);m.seal(a.output,['results.json']);print('CT10 common-space comparison complete',flush=True);return 0
if __name__=='__main__':raise SystemExit(main())
