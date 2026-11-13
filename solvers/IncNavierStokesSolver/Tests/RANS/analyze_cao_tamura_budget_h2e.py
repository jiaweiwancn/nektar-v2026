#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Read-only compact reduction and independent tensor replay for H2e."""
import argparse,json,math,os,time
from pathlib import Path
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from scipy.optimize import minimize,brentq
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,require
from validate_cls_epawf_m1 import indexed
from search_cao_tamura_budget_h2d import tensors
from search_cao_tamura_budget_h2e import condition

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--preflight',type=Path,required=True);p.add_argument('--old-native',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root;pre=load(a.preflight);doc=load(r/'context-v1/context.json');f=doc['context'];arrays=Path(doc['array_directory'])
    def arr(name):d=f['files'][name];require(sha(arrays/d['file'])==d['sha256'],'context changed');return np.memmap(arrays/d['file'],mode='r',dtype=d['dtype'],shape=tuple(d['shape']))
    T,vc=condition(np.array(arr('mass_gram')));G=np.array(arr('scalar_gram'));Tk,kc=condition(G[:8,:8]);Te,ec=condition(G[8:,8:]);Ts=np.zeros((16,16));Ts[:8,:8]=Tk;Ts[8:,8:]=Te
    obs=arr('observations');Z=np.array(arr('basis'));raw=obs[:,18:210].reshape(len(obs),32,6);base=np.array(arr('base')[0]);cov=obs[:,16].astype(bool);inner=obs[:,15].astype(bool);nt=obs[:,14];phys=obs[:,3:6];K0=pre['heights'][-1]['K0'];stages={}
    def replay(w):
        g=np.array(w['gradient']);k=w['k'];e=w['epsilon'];require(not w['inner'],'outer replay expected');t=indexed(g,k,e,'physical');S=g+g.T;Om=g-g.T;sr=np.zeros((3,3));sr[0,2]=sr[2,0]=pre['u_tau']/(pre['kappa_CLS']*w['physical'][2]);gamma=pre['kappa_CLS']*pre['u_tau']*w['physical'][2]
        terms=dict(isotropic=2*k/3*np.eye(3),linear=-float(t['nut'].real)*S,quadratic=t['quadratic'].real,cubic=.4*t['cubic'].real,regularization=-gamma*(S-sr));R=sum(terms.values());ev,V=np.linalg.eigh(R);v=V[:,0];native=np.array(w['tensor']);err=float(np.max(abs(R-native)));require(err<1e-11,'independent tensor differs')
        qs=k/e*np.linalg.norm(S)/math.sqrt(2);qo=k/e*np.linalg.norm(Om)/math.sqrt(2)
        return dict(eigenvalues=ev.tolist(),critical_eigenvector=v.tolist(),critical_contributions={n:float(v@A@v) for n,A in terms.items()},tensor_error=err,native_minimum_error=abs(float(ev[0])-w['eigenvalues'][0]),branch='strain' if qs>=qo else 'rotation',q_strain=qs,q_rotation=qo,C_mu=float(t['cmu'].real) if 'cmu' in t else None,nu_t=float(t['nut'].real),gamma=gamma,production=-float(np.sum(R*g)),Original_same_state_minimum=float(np.linalg.eigvalsh(terms['isotropic']+terms['linear']+terms['quadratic']+t['cubic'].real)[0]))
    for stage,folder,cert in [('A','stage-A-final','certify-A'),('B','stage-B','certify-B')]:
        c=load(r/folder/'best-candidate.json');s=load(r/folder/'progress.json');n=load(r/cert/'certification.json');native=load(Path(n['receipt']['output']));x=np.array(c['conditioned_velocity_coefficients']);ZH=np.vstack([T.T@Z[:32],Ts.T@Z[32:]]);sg=np.zeros(16);amplitudes=x
        if stage=='B':sg=np.array(c['conditioned_scalar_coefficients'])@Ts.T;amplitudes=np.r_[x,c['conditioned_scalar_coefficients']]
        delta=amplitudes@ZH[:len(amplitudes)]
        free=np.array(f['rows']['velocity']+(f['rows']['k']+f['rows']['epsilon'] if stage=='B' else []),int);actual=base.copy();actual[free]+=delta[free];require(state_hash(actual)==c['state_sha256'],'exact state reconstruction differs')
        d=np.einsum('nbc,b->nc',raw,x@T.T);g=np.zeros((len(obs),3,3));g[:,0,0]=obs[:,6]+d[:,0];g[:,0,2]=obs[:,7]+d[:,1];g[:,2,0]=obs[:,8]+d[:,2];g[:,2,2]=obs[:,9]+d[:,3];k=np.array(obs[:,12]);e=np.array(obs[:,13]);k+=obs[:,210:218]@sg[:8];e+=obs[:,218:]@sg[8:]
        began=time.monotonic();R=tensors(g,k,e,inner,nt,phys,pre['u_tau'],pre['kappa_CLS']);val=np.linalg.eigvalsh(R[cov])[:,0];scan=time.monotonic()-began;require(abs(float(val.min())/K0-c['min_lambda_over_K0'])<1e-9,'saved terminal scan omitted a changed state')
        roof=obs[:,2]==.8;far=phys[:,2]>=.4;valid=~inner;kr=k[valid]/obs[valid,12];er=e[valid]/obs[valid,13];lim=np.flatnonzero(cov)[int(np.argmin(val))]
        # A local envelope diagnostic keeps the actual limiting velocity fixed.
        # It is neither a global field search nor an incompatibility proof.
        gl=g[lim:lim+1];kl=obs[lim,12];el=obs[lim,13]
        def margin(ke):return float(np.linalg.eigvalsh(tensors(gl,np.array([kl*ke[0]]),np.array([el*ke[1]]),np.array([False]),np.array([0.]),phys[lim:lim+1],pre['u_tau'],pre['kappa_CLS']))[0,0])
        grid=[(margin((u,v)),u,v) for u in np.linspace(.5,2,41) for v in np.linspace(.5,2,41)];best=max(grid);local=minimize(lambda ke:-margin(ke),[best[1],best[2]],bounds=[(.5,2),(.5,2)],method='L-BFGS-B')
        # Directional finite-difference checks retain the actual max branch.
        direction=np.array([[.2,0,-.3],[0,0,0],[.1,0,-.2]]);eps=2e-6
        def local_tensor(step):return tensors(gl+step*direction,np.array([k[lim]]),np.array([e[lim]]),np.array([False]),np.array([0.]),phys[lim:lim+1],pre['u_tau'],pre['kappa_CLS'])[0]
        deriv=(local_tensor(eps)-local_tensor(-eps))/(2*eps);half=(local_tensor(eps/2)-local_tensor(-eps/2))/eps
        extra=n['extra_sampling'];qgrad=np.r_[1.,x];roofQ=np.array(arr('roof'))
        stages[stage]=dict(candidate_sha256=c['state_sha256'],status=s['status'],termination=s['termination'],iterations=len(s['history']),search_seconds=s['search_seconds'],full_scans=s['full_scans'],active_derivative_batches=s['active_derivative_batches'],peak_RSS_KiB=s['peak_RSS_KiB'],gates=n['gates'],audit=n['audit'],extra_sampling=extra,contract=n['contract'],mass_distance=c['velocity_mass_distance'],component_response_cap=c['maximum_sampled_component_increment'],K_envelope=[float(kr.min()),float(kr.max())],epsilon_envelope=[float(er.min()),float(er.max())],roof_sampled_U=[float((obs[roof,10]+d[roof,4]).min()),float((obs[roof,10]+d[roof,4]).max())],roof_maximum_delta_U=float(np.max(abs(d[roof,4]))),far_field_maximum_component_increment=float(np.max(abs(d[far,4:]))),roof_streamwise_gradient_integral=float(qgrad@roofQ@qgrad),base_roof_streamwise_gradient_integral=float(roofQ[0,0]),terminal_saved_state_reconstructed_bitwise=True,postprocessing_full_scans=1,postprocessing_scan_seconds=scan,independent_extra_tensor_replay=replay(extra['witness']),local_scalar_envelope=dict(minimum_on_best_grid=best[0],optimized_minimum=margin(local.x),k_ratio=float(local.x[0]),epsilon_ratio=float(local.x[1]),scope='fixed one local velocity gradient; finite grid plus local bounded optimizer; no field/global infeasibility certificate'),directional_check=dict(central_vs_half_max_error=float(np.max(abs(deriv-half))),gradient_direction=direction.tolist(),step=eps))
    old=next(w for w in load(a.old_native)['best_native']['seed_target_witnesses'] if not w['inner'] and w['is_audit_minimizer']);gg=np.array(old['represented_gradient']);kk=old['k'];ee=old['epsilon'];pp=np.array([old['physical_coordinate']]);fun=lambda scale:float(np.linalg.eigvalsh(tensors(gg[None],np.array([kk*scale]),np.array([ee*scale]),np.array([False]),np.array([0.]),pp,pre['u_tau'],pre['kappa_CLS']))[0,0]);scale=brentq(fun,1,10)
    out=dict(version='H2eCompactEvidenceV1',physical_hill_steps=0,physical_hill_attempts=0,context=dict(base_sha256=f['base_sha256'],seconds=f['seconds'],inclusive_seconds=doc['receipt']['seconds'],peak_RSS_KiB=doc['peak_RSS_KiB'],binary_bytes=sum((arrays/d['file']).stat().st_size for d in f['files'].values()),unique_observations=f['files']['observations']['shape'][0],memberships=f['files']['membership']['shape'][0],physical_volume=f['volume'],base_divergence_L2=f['base_divergence_L2'],base_maximum_quadrature_divergence=f['base_maximum_quadrature_divergence'],homogeneous_continuity_max=f['homogeneous_continuity_max'],velocity_rank=vc,k_rank=kc,epsilon_rank=ec),stages=stages,old_H2d_roof_equal_scalar_scale_diagnostic=scale,classification='native and extra covariance rejection; bounded nonconvex search exhausted; no hill solution nonexistence claim',native_equations_unchanged=True,flat_spatial_gate='unclosed; retained 4399-5339-slot plan only')
    atomic(a.output,out)
if __name__=='__main__':main()
