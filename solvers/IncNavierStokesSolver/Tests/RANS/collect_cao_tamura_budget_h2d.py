#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H2d executed evidence, independent tensors/derivatives and plots."""
import argparse,json,hashlib,math,os
from pathlib import Path
from datetime import datetime,timezone,timedelta
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',BLIS_NUM_THREADS='1')
import numpy as np
from validate_cls_epawf_m1 import indexed
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,load,require

def small(v):
    if isinstance(v,dict):return {k:small(x) for k,x in v.items() if k not in ('state','base','basis','observations','physical_samples','gradient_basis','velocity_basis')}
    if isinstance(v,list):return [small(x) for x in v]
    return v

def selected(g,k,e,z,us,kap,branch='physical'):
    v=indexed(g,k,e,branch);S=g+g.T;sr=np.zeros((3,3));sr[0,2]=sr[2,0]=us/(kap*z)
    return 2*k*np.eye(3)/3-v['nut']*S+v['quadratic']+.4*v['cubic']-kap*us*z*(S-sr)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--revision',required=True);a=p.parse_args();r=a.root;a.output.mkdir(exist_ok=False)
    def read(n):return load(r/n)
    def save(n,d):atomic(a.output/(n+'.json'),d)
    def result(root,name):
        p=read(root+'/progress.json');require(p['status']=='completed','unfinished '+root);return load(Path(p['results'][name]['result_file']))
    final=read('final-checks-v3/progress.json');require(final['status']=='passed','final native contract tests incomplete')
    pre=load(r.parent/'cls-cao-tamura-budget-h2c-20261112/preflight-final.json');us=pre['u_tau'];kap=pre['kappa_CLS'];H=.0065;K=pre['heights'][-1]['K0']
    replay={name:result('replay-final-v1',name) for name in ('TangentialNeutralV2','PhysicalHeightInterfaceV2')}
    best=result('best-final-v1','65-hill');search=read('search-v1/progress.json');candidate=read('search-v1/best-candidate.json')
    require(best['state_sha256']==candidate['state_sha256'],'saved optimizer/native candidate differs')
    require(best['audit']['accepted_awf_hard_admissible']==[0] and search['status']=='bounded_search_rejected','unexpected changed hill classification')
    scalar_start=int(best['seed_constitutive']['offsets'][2]);require(state_hash(best['state'][scalar_start:])==state_hash(replay['PhysicalHeightInterfaceV2']['state'][scalar_start:]),'scalar/pressure bits changed during velocity search')
    save('native',dict(capture_Asia_Shanghai=datetime.now(timezone(timedelta(hours=8))).isoformat(),status='initialization_rejected',physical_hill_steps=0,physical_hill_attempts=0,
        conditional_pair='unexecuted: independent complete covariance gate failed',best_candidate_sha256=best['state_sha256'],best_native=small(best),
        setup_state_sha256=replay['PhysicalHeightInterfaceV2']['state_sha256'],all_nonvelocity_preserved=True,
        nonvelocity_verification='native candidate reader compares every non-velocity and prescribed coefficient exactly to the reconstructed, hash-validated base',
        conclusion='bounded frozen-scalar velocity space/response limits did not find covariance feasibility; no global infeasibility or ellipticity inference'))
    save('curvature',dict(replay=small(replay),reference=small(read('independent-replay-v2.json')),old_defect='sorted kappa1 was zero on negative shoulders; signed stream curvature is t_stream^T B t_stream',
        historical_H2c_roots_unchanged=True,initializer_versions=['TangentialNeutralV2','PhysicalHeightInterfaceV2']))
    save('search',small(search))
    tensor_checks=[];derivatives=[]
    for label,w in [('minimum_inner',best['seed_constitutive']['minimum_inner']),('minimum_outer',best['seed_constitutive']['minimum_outer'])]:
        g=np.array(w['velocity_gradient']);k=w['k'];e=w['epsilon'];z=w['physical_coordinate'][2]
        if w['inner']:R=2*k*np.eye(3)/3-w['selected']['eddy_viscosity']*(g+g.T)
        else:R=selected(g,k,e,z,us,kap).real
        err=float(np.max(abs(R-np.array(w['selected']['tensor']))));require(err<1e-11,'independent tensor mismatch')
        tensor_checks.append(dict(label=label,geometry=w['geometry'],element=w['element'],chart=w['coordinate'],physical=w['physical_coordinate'],tensor_error=err,eigenvalues=np.linalg.eigvalsh(R).tolist()))
    w=best['extra_sampling']['witness'];R=selected(np.array(w['gradient']),w['k'],w['epsilon'],w['physical'][2],us,kap).real
    require(np.max(abs(R-np.array(w['tensor'])))<1e-11,'independent extra-grid tensor mismatch')
    tensor_checks.append(dict(label='extra_GL12',tensor_error=float(np.max(abs(R-np.array(w['tensor'])))),eigenvalues=np.linalg.eigvalsh(R).tolist()))
    w=best['seed_constitutive']['minimum_outer'];k=w['k'];e=w['epsilon'];z=w['physical_coordinate'][2]
    D=np.array([[.1,0,.5],[0,0,0],[-.3,0,-.1]]);D/=np.linalg.norm(D)
    neutral=np.zeros((3,3));neutral[0,2]=us/(kap*z)
    for label,g in [('active',np.array(w['velocity_gradient'])),('branch_equal_neutral',neutral),('near_branch',neutral+1e-10*D)]:
        levels=[]
        for sign in (1,-1):
            delta=sign*D;trial=g+1e-6*delta;s=trial+trial.T;o=trial-trial.T;branch='strain' if np.sum(s*s)>=np.sum(o*o) else 'rotation'
            derivative=selected(g.astype(complex)+1j*1e-30*delta,k,e,z,us,kap,branch).imag/1e-30
            errors=[]
            for h in (2e-5,1e-5,5e-6):
                fd=(selected(g+h*delta,k,e,z,us,kap).real-selected(g,k,e,z,us,kap).real)/h
                errors.append(float(np.max(abs(fd-derivative))))
            require(errors[-1]<1e-5,'branch-aware directional derivative unresolved')
            levels.append(dict(sign=sign,selected_branch=branch,one_sided_FD_errors=errors))
        derivatives.append(dict(case=label,levels=levels,semantics='one-sided directional q=max branches; no differentiability assertion at equality'))
    save('verification',dict(final_native=final,initial_import_suite=read('import-validation-v1/progress.json'),tensor_replays=tensor_checks,directional_checks=derivatives,kernel_log=(r/'kernel-checks.log').read_text(),
        no_accepted_hill_certificate=True,certified_flat_no_Lift=True,ordinary_changed_dt_and_geometry_rejected=True))
    flats={n:result('flat-flux-final-v1',n) for n in ('flat-baseline','flat-radial-h','flat-p4-diagnostic')}
    flats['flat-combined-confirmation']=result('flat-flux-combined-v1','flat-combined-confirmation')
    volume=[.8*.8]*3+[.8*(.8-H)];rates=[us*us/H,us*us/H,K*us/H,(us**3/(kap*H))*us/H]
    rows=[]
    for name,d in flats.items():
        row=dict(case=name,size=d['size'],selected_scaled=d['selected_spatial_norms'][:4],physical_free_mass_dual=[v*s*math.sqrt(V) for v,s,V in zip(d['selected_spatial_norms'],rates,volume)],reference_errors=d['physical_reference']['flat_reference_errors'],radial_sum_error=d['radial_global_free_sum_error'])
        require(row['radial_sum_error']<1e-10,'native panel/global shared-row sum does not close');rows.append(row)
    save('flat',dict(cases=small(flats),comparison=rows,equilibrium_target=1e-3,all_blocks_pass=False,
        statement='Two tangential panels are a diagnostic reduction. p4 uses supported tensor p4/P2; no isolated radial-p order claim. Only physical representation errors and shared-row residual contractions are compared.',
        analytic_Yap='nonzero exact analytical epsilon_yap contribution retained; its cancellation with other physical terms is part of exact equilibrium, not an error',
        transfer=dict(original_pilot_slots=2804,full_h_allocation_six_panel_slots=5339,combined_cubic_allocation_six_panel_slots=4399,limit=3000,conclusion='full useful allocation does not fit original pilot limit; no guard increase or evolving new mesh')))
    save('source-revisions',dict(base='c18e9ab7e6e22a4e268b04d07a0a05ec25083532',revision=a.revision,source_export=str(r/'source'),source_hashes=read('source-final-hashes.json'),
        final_native_runtime=final['runtime'],final_probe_runtime=final['probe'],data_volume_UUID='5AB8A669B8A642FD',run_root=str(r),
        calibration_preflight_sha256=sha(r.parent/'cls-cao-tamura-budget-h2c-20261112/preflight-final.json'),resource_contract='build parallel 4, pilots rank 1/all computational threads 1; independent external build 2 preserved; aggregate below 72'))
    save('cost',dict(search_usage_wall_user_system_RSS_KiB=(r/'search-v1-usage.txt').read_text(),search_tensor_scan_seconds=search['seconds_per_full_tensor_scan'],native_setup_receipt=read('setup-v1/65-hill-probe/receipt.json'),
        best_audit=read('best-final-v1/progress.json')['results']['65-hill']['receipt'],flat=read('flat-flux-final-v1/progress.json'),final_native=final.get('native'),hill_step_cost='unmeasured: no admissible seed'))
    # Reproducible scientific figures, not an interactive surrogate.
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    fig,ax=plt.subplots(1,2,figsize=(10,3.4));labels=['U','W','k','epsilon'];xx=np.arange(4)
    for row in rows:ax[0].semilogy(xx,row['selected_scaled'],'o-',label=row['case'].replace('flat-',''))
    ax[0].axhline(1e-3,color='k',linestyle='--',linewidth=.8);ax[0].set_xticks(xx,labels);ax[0].set_ylabel('native free-row mass-dual residual');ax[0].legend(fontsize=7)
    h=search['history'];ax[1].plot([v['seconds'] for v in h],[v['minimum_lambda_over_K0'] for v in h]);ax[1].axhline(0,color='k',linewidth=.8);ax[1].set_xlabel('search wall time (s)');ax[1].set_ylabel('minimum covariance eigenvalue / K0');fig.tight_layout();fig.savefig(a.output/'qualification.png',dpi=180);plt.close(fig)
    files=[]
    for path in sorted(r.rglob('*')):
        if not path.is_file() or path.is_symlink() or a.output in path.parents or path.relative_to(r).parts[0] in ('source','build') or path.suffix=='.tar':continue
        files.append(dict(path=str(path.relative_to(r)),bytes=path.stat().st_size,sha256=sha(path)))
    save('artifact-index',dict(root=str(r),files=files,total_bytes=sum(x['bytes'] for x in files)))
if __name__=='__main__':main()
