#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Collect compact evidence from completed H2i jobs; never starts optimization."""
import argparse,os,time,hashlib
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_h2i_convex import Context
from cao_tamura_h2h_scalar import Budget
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,require,THREADS
from run_cao_tamura_budget_h2h import runtime_identity

def file(p):return dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size)

def root_limit(Q,limit,x,d):
    y=np.r_[1.,x];z=np.r_[0.,d];A=float(z@Q@z);B=float(2*y@Q@z);C=min(0.,float(y@Q@y-limit))
    if A<1e-30:return -C/B if B>0 else float('inf')
    return max(0.,float((-B+np.sqrt(max(0.,B*B-4*A*C)))/(2*A)))

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--root',type=Path,required=True);ap.add_argument('--parent',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);ap.add_argument('--revision',required=True);a=ap.parse_args();os.environ.update(THREADS);start=time.monotonic();r=a.root;dest=a.output;dest.mkdir(exist_ok=False)
    ctx=Context(r/'joint-context/context.json');d=load(r/'independent-reduction.json');search=load(r/'joint-search/result.json');theta=search['retention']['theta'];stamp=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
    provenance=dict(capture_Asia_Shanghai=stamp,reviewed_parent='409f7d98b21ddee999dfe67b1600810b7091b769',implementation_revision=a.revision,root=str(r))
    out=dict(context_sha256=sha(ctx.path),native_response_verification=ctx.meta['verification'],directions=ctx.meta['direction_counts'],conditioning=ctx.doc['conditioning'],observation_blocks=ctx.blocks,native_memberships=ctx.meta['native_memberships'],all_memberships=ctx.meta['all_memberships'],unique_evaluations=ctx.N,families=ctx.meta['families'],arrays=ctx.meta['files'],basis_sha256=ctx.doc['basis_sha256'],algebra={k:v for k,v in ctx.meta['algebra'].items() if k!='free_velocity_rows'},identity={k:ctx.doc[k] for k in ('geometry_identity','inventory_identity','reference_sha256','defect_sha256','patch_signature')},model=ctx.meta['model_specification'])
    atomic(dest/'context.json',out)
    out.update(native_coefficient_offsets=ctx.meta['native_coefficient_offsets'],native_coefficient_sizes=ctx.meta['native_coefficient_sizes'],permitted_rows=ctx.meta['rows'])
    atomic(dest/'context.json',out)
    phases=[]
    for p in search['phases']:
        v=p['result'];phases.append(dict(anchor=p['anchor'],buffered=p['buffered'],termination=v['termination'],delta=v['delta'],epsilon_band=v['epsilon_band'],best_summary=v['best_summary'],rounds=len(v['rounds']),cut_count=v['cut_count'],dual_certificate=v['proof']['receipt']['certificate']))
    atomic(dest/'optimizer.json',dict(provenance,classification=search['classification'],budget=search['budget'],retention=search['retention'],anchors=load(r/'joint-search/anchors.json'),bounds=search['bounds'],phases=phases,anchor_span_roundoff_max=search['anchor_span_roundoff_max'],peak_RSS_KiB=search['peak_RSS_KiB']))
    replays=[]
    members=ctx.arrays['members'];points=ctx.arrays['points']
    for replay in d['replays']:
        v=dict(replay);v['limiting_constraints']=[]
        for m in replay['limiting_constraints']:
            m=dict(m)
            if m['kind']=='epsilon_bound':
                nout=np.flatnonzero(~ctx.block('inner').astype(bool));j=m['hard_row'];m['point']=int(nout[j%len(nout)]);m['epsilon_bound']='upper' if j<len(nout) else 'lower'
            if 'point' in m:
                i=m['point'];m.update(chart=ctx.block('chart')[i].tolist(),physical=ctx.block('physical')[i].tolist(),memberships=members[members[:,1]==i].tolist())
            v['limiting_constraints'].append(m)
        replays.append(v)
    atomic(dest/'active-constraints.json',dict(provenance,replays=replays,independent_reduction=file(r/'independent-reduction.json'),replay_scans=d['cumulative_replay_scans'],replay_seconds=d['cumulative_replay_seconds']))
    native={};access={};extra_scans=0
    for anchor in ('A','B'):
        v=d['native'][anchor];receipt=load(r/('native-'+anchor)/'receipt.json');raw=load(Path(receipt['output']));native[anchor]=dict(v,reference_sha256=raw['reference_sha256'],defect_sha256=raw['defect_sha256'],outer_initializer_contract=raw['outer_initializer_contract'],native_budget=raw['native_budget'],physical_sample_artifact=dict(file=receipt['output'],field='physical_samples',sha256=receipt['output_sha256']),receipt_file=file(r/('native-'+anchor)/'receipt.json'))
        P=ctx.problem(anchor,theta);candidate=load(r/'joint-search'/(anchor+'-actual')/'best-candidate.json');x=np.array(candidate['joint_coefficients']);s=P.scan(x,1e-4*P.K0,Budget(seconds=60));extra_scans+=1
        rr=load(r/('accessibility-'+anchor+'-receipt.json'));result=load(Path(rr['output']))['context'];G=P.mass[1:,1:];wits=[]
        for w in result['witnesses']:
            w=dict(w);overlap=np.array(w['span_derivative']);response=np.linalg.solve(G,overlap);norm=np.sqrt(max(0.,response@G@response));direction=response/norm if norm>1e-20 else response
            limits=dict(mass=root_limit(P.mass,.12**2,x[:41],direction))
            for order,Q,limit in zip((16,20,24),P.div,P.div_limits):limits['divergence_GL'+str(order)]=root_limit(Q,limit,x[:41],direction)
            du=np.einsum('njc,j->nc',P.VB[:,:,4:],direction,optimize=False);cap=float('inf')
            for sign in (1.,-1.):
                h=sign*du;mask=h>1e-25
                if np.any(mask):cap=min(cap,float(np.min((.5-sign*s['velocity_increment'][mask])/h[mask])))
            limits['sampled_cap']=max(0.,cap)
            w.update(mass_unit_span_self_Rayleigh_derivative=float(overlap@direction),single_direction_cap_steps={k:(v if np.isfinite(v) else None) for k,v in limits.items()},single_direction_available_step=min(limits.values()),scope='constrained-Riesz accessibility and one fixed-epsilon span ascent ray; not joint feasibility or an enrichment experiment')
            wits.append(w)
        access[anchor]=dict(result,witnesses=wits,receipt=rr)
    atomic(dest/'native.json',dict(provenance,diagnostic_states=native,classification='neither diagnostic state is a startup seed',repeatable_checks=load(r/'candidate-checks.json'),physical_attempts=0,accepted_steps=0,fresh_GL28=0,certificates=0))
    atomic(dest/'accessibility.json',dict(provenance,results=access,additional_complete_scans=extra_scans))
    flat={}
    oldnative=load(r/'historical-evidence/h2g-native.json')
    for m,f in d['flat'].items():
        flat[m]={k:f[k] for k in ('base_sha256','state_sha256','H1_length','setup','scalar_domain','fields_changed')}
        require(f['base_sha256']==oldnative['models'][m]['cold']['flat']['state_sha256'],'flat immutable baseline differs from H2g')
        for kind in ('baseline','projected'):
            z=f[kind];flat[m][kind]={k:z[k] for k in ('size','selected_spatial_norms','analytical_free_load_norms','radial_errors','error_columns','one_sided_derivative_jumps','radial_global_free_sum_error','global_components','analytical_component_norms')};flat[m][kind]['physical_reference']={k:z['physical_reference'][k] for k in ('flat_reference_errors','wall_reference_scales','outer_native_scalar_minima','stress_work_relative_defect')};flat[m][kind]['signed_radial_projections']={k:v['signed_dual_projection_on_total'] for k,v in z['radial_components'].items()}
        flat[m]['baseline_hash_reproduces_H2g']=True
    atomic(dest/'flat.json',dict(provenance,models=flat,adopted=False,physical_steps=0,interpretation='same-space M+H_A^2 stiffness projection of free U/k; incompatible induced continuity and remaining U/k residual preclude automatic adoption'))
    frozen=load(r/'probe-runtime/identity.json');same=sorted(frozen['original']['libraries'].values())==sorted(frozen['frozen']['libraries'].values()) and frozen['original']['sha256']==frozen['frozen']['sha256'];require(same,'flat ambient closure differs from frozen closure')
    current=runtime_identity(r/'probe-runtime/RANSWallResolutionProbe');require(current['sha256']==frozen['frozen']['sha256'] and current['libraries']==frozen['frozen']['libraries'],'validation runtime changed')
    solver=runtime_identity(a.parent/'runtime-v2/IncNavierStokesSolver')
    test=load(r/'unit-qualified.json');require(test['passed'],'focused tests failed')
    atomic(dest/'verification.json',dict(provenance,focused=test,native_affine=ctx.meta['verification'],dual_reconstruction_error=max(v['cut_reconstruction_error'] for v in replays),independent_tensor_error=max(q[k]['tensor_error'] for n in d['native'].values() for q in n['seen'].values() for k in ('inner_covariance_witness_independent','outer_covariance_witness_independent')),flat_baseline_hashes_match_H2g=True,solver_changed=False,completed_H2g_native_restart_import_rollback_evidence_reused=True))
    source=r/'source';names=['library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp','library/UnitTests/SolverUtils/RANSWallJointFeasibilityContext.h','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionJointInitializer.cpp','solvers/IncNavierStokesSolver/Tests/RANS/CMakeLists.txt']
    names += ['solvers/IncNavierStokesSolver/Tests/RANS/'+s for s in ('cao_tamura_h2i_convex.py','run_cao_tamura_budget_h2i.py','analyse_cao_tamura_budget_h2i.py','collect_cao_tamura_budget_h2i.py','test_cao_tamura_budget_h2i.py','run_cao_tamura_budget_h2f_startup.py')]
    atomic(dest/'source-revisions.json',dict(provenance,source_hashes={n:sha(source/n) for n in names},selected_probe=current,selected_solver=solver,flat_actual_probe=frozen['original'],flat_matches_frozen_closure=same,storage_UUID='5AB8A669B8A642FD',threads=THREADS,build_workers=4,native_ranks=1,production_guard=3000))
    times={}
    for p in r.glob('*usage.txt'):
        try:w,u,s,mem=map(float,p.read_text().split());times[p.name]=dict(wall_seconds=w,user_seconds=u,system_seconds=s,peak_RSS_KiB=mem)
        except ValueError:pass
    total_scans=search['budget']['complete_scans']+d['cumulative_replay_scans']+extra_scans;require(total_scans<=120,'aggregate complete-scan budget')
    atomic(dest/'cost.json',dict(provenance,measurements=times,accessibility_seconds={a:load(r/('accessibility-'+a+'-receipt.json'))['seconds'] for a in ('A','B')},LP_rounds=search['budget']['rounds'],complete_optimizer_scans=search['budget']['complete_scans'],independent_replay_scans=d['cumulative_replay_scans'],accessibility_reduction_scans=extra_scans,aggregate_complete_scans=total_scans,optimizer_seconds=search['budget']['seconds'],replay_seconds=d['cumulative_replay_seconds'],collection_seconds=time.monotonic()-start))
    calibration=load(a.parent/'comparison.json')
    terminals={k:dict(state_sha256=v['candidate_sha256'],native_covariance=v['native_covariance'],seen_covariance={o:[s['inner_minimum_covariance'],s['outer_minimum_covariance']] for o,s in v['seen'].items()}) for k,v in d['native'].items()}
    atomic(dest/'comparison-manifest.json',dict(provenance,version='H2iContinuingTwoModelComparisonV1',status='open; initializer sampled obstruction plus unresolved spatial/development/temporal/experimental gates',calibration_manifest=file(a.parent/'comparison.json'),unchanged_full_precision_calibration=calibration,CLS='deployed model and retained H2f/H2g results unchanged; no qualifying hill seed',standard_cold_sha256=ctx.meta['cold_sha256'],standard_B_sha256=ctx.meta['B_sha256'],anchors=load(r/'joint-search/anchors.json'),standard_terminal_diagnostics=terminals,physical_hill_attempts=0,physical_hill_accepted_steps=0,next_stages=['native velocity/interface representation and bounded feasible initialization','common flat spatial qualification','bounded developed flow','temporal/spatial sensitivity','held-out experimental assessment']))
    gates=dict(passed=['calibration/model preservation','full 41+8 native response and continuity contract','actual-limit finite sampled dual obstruction independently replayed','repeatable seen-grid witnesses and scalar/profile checks','fail-fast invocation sequencing','one paired zero-step flat H1 diagnostic'],failed=['hill full covariance','GL24 inner startup margin','flat U/k 1e-3 spatial target'],unexecuted=['fresh GL28: seen/native prerequisite failed','V2 issuance/import: no qualified state','physical coarse8/fine16 and 4+4 restart/curved rollback: no qualified state'])
    atomic(dest/'final-status.json',dict(provenance,classification='verified precisely scoped sampled obstruction for both declared fixed-k anchors',gates=gates,repeatable_native_checks=2,fresh_GL28_audits=0,physical_hill_attempts=0,physical_hill_accepted_steps=0,long_continuations=0))
    entries=[]
    for folder in ('joint-context','joint-search','native-A','native-B','accessibility-A','accessibility-B','flat-H1'):
        for p in sorted((r/folder).rglob('*')):
            if p.is_file():entries.append(file(p))
    for p in sorted(r.iterdir()):
        if p.is_file() and p.suffix in ('.json','.txt','.log') and p.name not in ('collection-usage.txt','collection.log'):entries.append(file(p))
    atomic(dest/'artifact-index.json',dict(provenance,files=entries,runtimes=[file(r/'probe-runtime/identity.json')],source_hashes=file(dest/'source-revisions.json'),heavy_data_in_Git=False))

if __name__=='__main__':main()
