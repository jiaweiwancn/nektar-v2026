#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only H2f native tensor replay and constraint-aware accessibility reduction."""
import argparse,os,time
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from cao_tamura_initializer_context import NativeContext
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,require
from validate_cls_epawf_m1 import indexed
from search_cao_tamura_budget_h2d import tensors

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('root','preflight','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();r=a.root;pre=load(a.preflight);K0=pre['heights'][-1]['K0'];c=NativeContext(r/'enriched-context/context.json');T,Ts,rank=c.transforms();base=np.array(c.array('base')[0]);basis=np.array(c.array('basis'));nv=c.counts['velocity'];Z=np.vstack([T.T@basis[:nv],Ts.T@basis[nv:]]);raw=c.block('velocity_basis');inner=c.block('inner').astype(bool);cov=c.block('covariance').astype(bool);dep=c.block('dependent_epsilon').astype(bool);coords=c.block('chart');physical=c.block('physical');stages={}
    def independent(w,nt=None):
        g=np.array(w['gradient']);k=w['k'];S=g+g.T;terms=dict(isotropic=2*k/3*np.eye(3))
        if nt is not None:terms['linear']=-nt*S
        else:
            e=w['epsilon'];model=indexed(g,k,e,'physical');z=w['physical'][2];sr=np.zeros((3,3));sr[0,2]=sr[2,0]=pre['u_tau']/(pre['kappa_CLS']*z);terms.update(linear=-float(model['nut'].real)*S,quadratic=model['quadratic'].real,cubic=.4*model['cubic'].real,regularization=-pre['kappa_CLS']*pre['u_tau']*z*(S-sr))
        R=sum(terms.values());evec,V=np.linalg.eigh(R);v=V[:,0];error=float(np.max(abs(R-np.array(w['tensor']))));require(error<1e-11,'independent tensor mismatch')
        return dict(eigenvalues=evec.tolist(),critical_eigenvector=v.tolist(),signed_critical_contributions={n:float(v@A@v) for n,A in terms.items()},tensor_max_error=error,gradient=g.tolist(),physical=w['physical'],chart=w['coordinate'],geometry=w['geometry'],element=w['element'],k=k,epsilon=w['epsilon'],production=-float(np.sum(R*g)))
    for stage,folder,cert in [('A','stage-A-final','certify-A'),('B','stage-B','certify-B')]:
        candidate=load(r/folder/'best-candidate.json');progress=load(r/folder/'progress.json');receipt=load(r/cert/'certification.json');native=load(Path(receipt['probe_output_file']));x=np.array(candidate['conditioned_velocity_coefficients']);sc=np.array(candidate['conditioned_scalar_coefficients']);amplitude=np.r_[x,sc];state=base.copy();free=np.array(c.meta['rows']['velocity']+(c.meta['rows']['k']+c.meta['rows']['epsilon'] if stage=='B' else []),int);state[free]+=(amplitude@Z[:len(amplitude)])[free];require(state_hash(state)==candidate['state_sha256'],'terminal state reconstruction differs')
        d=np.einsum('nbc,b->nc',raw,T@x);g=c.gradient()
        for col,(i,j) in enumerate(((0,0),(0,2),(2,0),(2,2))):g[:,i,j]+=d[:,col]
        k=c.block('k').copy();e=c.block('epsilon').copy()
        if len(sc):coeff=Ts@sc;k+=c.block('k_basis')@coeff[:c.counts['k']];e+=c.block('epsilon_basis')@coeff[c.counts['k']:]
        R=tensors(g,k,e,inner,c.block('inner_nu_t'),physical,pre['u_tau'],pre['kappa_CLS']);ev=np.linalg.eigvalsh(R)[:,0];mins=[float(ev[inner&cov].min()),float(ev[~inner&cov].min())];require(np.max(abs(np.array(mins)-receipt['audit']['accepted_awf_covariance']))<1e-10,'batch/native inventories disagree');witnesses={}
        for region in ('inner','outer'):
            w=next(w for w in native['seed_target_witnesses'] if w['inner']==(region=='inner') and w['is_audit_minimizer']);nw=dict(w,gradient=w['represented_gradient'],tensor=w['native_tensor'],physical=w['physical_coordinate']);witnesses[region]=independent(nw,receipt['audit']['accepted_awf_inner_minimum_point'][5] if region=='inner' else None)
        witnesses['independent_GL20']=independent(receipt['extra_sampling']['witness'])
        gradient_metrics={}
        for label,mask in [('inner',inner),('outer',~inner),('roof',coords[:,2]==.8)]:
            norm=np.linalg.norm(g[mask],axis=(1,2));idx=np.flatnonzero(mask)[int(np.argmax(norm))];grad_delta=g[mask]-c.gradient()[mask];gradient_metrics[label]=dict(max_gradient_Frobenius=float(norm.max()),gradient_witness_chart=coords[idx].tolist(),gradient_witness_physical=physical[idx].tolist(),max_gradient_change=float(np.max(np.linalg.norm(grad_delta,axis=(1,2)))),max_absolute_divergence=float(np.max(abs(g[mask,0,0]+g[mask,2,2]))),max_component_velocity_increment=float(np.max(abs(d[mask,4:]))))
        stages[stage]=dict(candidate_sha256=candidate['state_sha256'],termination=progress['termination'],major_iterations=len(progress['history']),accepted_states=progress['accepted_states'],minimum_accepted_inner_margin=progress['minimum_accepted_inner_margin'],native_covariance=receipt['audit']['accepted_awf_covariance'],extra=receipt['extra_sampling'],gates=receipt['gates'],contract=receipt['contract'],witnesses=witnesses,gradient_metrics=gradient_metrics,search_seconds=progress['search_seconds'],peak_RSS_KiB=progress['peak_RSS_KiB'],full_tensor_scans=progress['full_tensor_scans'],complete_inner_scans=progress['complete_inner_scans'],conservative_full_scan_work=progress['conservative_full_scan_work'],active_derivative_batches=progress['active_derivative_batches'],max_inner_cuts=max(x['inner_cuts'] for x in progress['history']),max_response_caps=max(x['active_response_caps'] for x in progress['history']),max_scalar_caps=max(x['active_scalar_caps'] for x in progress['history']),mass_distance=candidate['velocity_mass_distance'],max_component_response=candidate['maximum_sampled_component_increment'],K_envelope=candidate['K_envelope'],epsilon_envelope=candidate['epsilon_envelope'],terminal_coefficients_reconstructed_bitwise=True,native_certification_seconds=receipt['receipt']['seconds'])
    # Constraint-aware local response diagnostics. These do not advance a field
    # or consume additional optimizer trials; each tests one retained direction.
    n=load(r/'native-combined/native-data.json');nr=n['context'];folder=Path(n['array_directory']);get=lambda name:np.fromfile(folder/nr['files'][name]['file'],dtype='<f8').reshape(nr['files'][name]['shape']);M=get('mass_free');Q=get('continuity_Q');global_r=get('global_responses');gg=get('gradients');free=np.array(nr['free_rows'],int);orth=Z[:nv,free];den=np.einsum('ni,ij,nj->n',global_r,M,global_r);eta=np.sum((orth@M@global_r.T)**2,axis=0)/den;old_eta=np.sum((orth[:32]@M@global_r.T)**2,axis=0)/den;native_eta=np.array([w['eta_squared'] for w in nr['witnesses']]);require(np.max(abs(old_eta-native_eta))<1e-9,'independent accessibility disagrees')
    qres=Q@global_r.T;kkt=M@global_r.T-gg.T;mu=Q@(-kkt);projected=kkt+Q.T@mu;kkt_error=float(np.max(abs(projected))/max(1,float(np.max(abs(gg)))))
    oldA=load('/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2e-20261113/stage-A-final/best-candidate.json');x=np.zeros(nv);x[:32]=oldA['conditioned_velocity_coefficients'];response=np.einsum('nbc,b->nc',raw,T@x)[:,4:];rawD=np.array(c.array('divergence'));Dt=np.zeros((nv+1,nv+1));Dt[0,0]=1;Dt[1:,1:]=T;D=Dt.T@rawD@Dt;z=np.r_[1,x];Dlimit=np.sqrt(D[0,0])*(1+1e-6)+1e-9;decisions=load(r/'enrichment/enrichment.json')['decisions'];new=0;directions=[];member=c.array('membership');geom=np.full(c.N,-1,int);geom[member[:,1]]=member[:,3]
    for entry in decisions:
        if not entry['retained']:continue
        j=32+new;new+=1;w=nr['witnesses'][entry['witness']];direction=raw[:,j,4:];positive=direction>0;negative=direction<0;cap=.04
        if np.any(positive):cap=min(cap,float(np.min((.5-response)[positive]/direction[positive])))
        if np.any(negative):cap=min(cap,float(np.min((-.5-response)[negative]/direction[negative])))
        v=np.zeros(nv+1);v[j+1]=1;aa=float(v@D@v);bb=float(2*v@D@z);cc=float(z@D@z-Dlimit**2);root=(-bb+np.sqrt(max(0,bb*bb-4*aa*cc)))/(2*aa) if aa>1e-24 else np.inf;mass=np.sqrt(max(0,.12**2-oldA['velocity_mass_distance']**2));limit=max(0,min(cap,root,mass));limiter='strong_divergence' if root<=min(cap,mass) else ('response_or_trust' if cap<=mass else 'mass')
        match=np.flatnonzero(np.all(coords==np.array(w['coordinate']),axis=1)&(geom==w['geometry'])&(~dep)&cov);require(len(match)>0,'witness lost native geometry/context');idx=int(match[0]);g=np.array(w['gradient']);dg=raw[idx,j,:4]
        for col,(i,jj) in enumerate(((0,0),(0,2),(2,0),(2,2))):g[i,jj]+=limit*dg[col]
        t=tensors(g[None],np.array([w['k']]),np.array([w['epsilon']]),np.array([False]),np.array([0.]),np.array([w['physical']]),pre['u_tau'],pre['kappa_CLS'])[0];ev=np.array(w['eigenvector']);slope=float(gg[entry['witness']]@Z[j,free]);directions.append(dict(source_witness=entry['witness'],basis_index=j,old_eta_squared=float(old_eta[entry['witness']]),enriched_eta_squared=float(eta[entry['witness']]),positive_step=limit,limiting_constraint=limiter,predicted_rayleigh_gain_over_K0=limit*slope,actual_local_rayleigh_gain_over_K0=float((ev@t@ev-w['critical_Rayleigh'])/K0),inner_gradient_increment_exactly_zero=True))
    out=dict(version='H2fCompactNativeV1',physical_hill_attempts=0,physical_hill_steps=0,bounded_searches=stages,accessibility=dict(old_eta_squared=old_eta.tolist(),enriched_eta_squared=eta.tolist(),relative_projected_KKT_error=kkt_error,one_direction_constraints= directions,native_response_witnesses=nr['witnesses']),schema_zero_replay=load(r/'zero-replay.json'),enrichment=load(r/'enrichment/enrichment.json'),sparse_mass=load(r/('sparse-mass-final' if (r/'sparse-mass-final/sparse-mass.json').exists() else 'sparse-mass')/'sparse-mass.json'),symbolic_mass=load(r/'symbolic-mass.json'),search_classification='bounded negative result; full native freedom not exhausted; no global incompatibility certificate; no hill startup certificate issued')
    atomic(a.output,out)
if __name__=='__main__':main()
