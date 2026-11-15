#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Adaptive numerical supports under one frozen H2k physical/aggregate policy."""
import argparse, os, time, resource
from pathlib import Path
import numpy as np
from scipy import linalg
from cao_tamura_h2d_state import load
from cao_tamura_h2i_convex import Context, phase, summary
from cao_tamura_h2h_scalar import Budget, BudgetExhausted
from cao_tamura_initializer_context import NativeContext
from run_cao_tamura_budget_h2i import form_anchor, make_candidate, B_HASH, COLD_HASH, now
from run_cao_tamura_budget_h2j import arrays, native
from run_cao_tamura_budget_h2h import runtime_identity
from cao_tamura_budget_h2g import case
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha, state_hash, require, THREADS
from analyse_cao_tamura_budget_h2i import replay
from cao_tamura_h2k_lp import equilibrated_lp
from cao_tamura_h2k_exact_caps import exact_caps_phase


def charge(a, event, **work):
    path=a.root/'aggregate.json'
    origin_file=a.root/'origin.json'
    if not origin_file.exists():atomic(origin_file,dict(started_epoch=time.time(),started_Asia_Shanghai=now(),wall_limit_seconds=28800))
    origin=load(origin_file)
    d=load(path) if path.exists() else dict(version='H2kAggregateWorkV1',created=now(),totals={},events=[])
    for key,value in work.items():d['totals'][key]=d['totals'].get(key,0)+value
    d['events'].append(dict(event=event,time=now(),work=work))
    atomic(path,d)
    limits=dict(attempts=96,construction_seconds=2400,rounds=240,scans=720,search_seconds=3600,epochs=4,native_checks=8,fresh_audits=1)
    for key,limit in limits.items():require(d['totals'].get(key,0)<=limit,'aggregate H2k budget '+key)
    require(time.time()-origin['started_epoch']<origin['wall_limit_seconds'],'H2k total wall budget including build/export')
    return d


def old(a):
    f=load(a.previous/'prepare/frozen.json');ctx=Context(a.retained_context if a.retained_context else Path(f['context']))
    B=np.array(load(a.parent/'standard-B/best-candidate.json')['state'])
    cold=np.array(NativeContext(a.parent/'standard-context/context.json').array('base')[0])
    require(state_hash(B)==B_HASH and state_hash(cold)==COLD_HASH,'immutable state binding')
    nv=ctx.meta['direction_counts']['velocity']
    Z=np.array(np.memmap(ctx.doc['basis_file'],dtype='<f8',mode='r',shape=(nv+16,len(B))))
    require(sha(Path(ctx.doc['basis_file']))==ctx.doc['basis_sha256'],'prior basis hash')
    return ctx,B,cold,Z


def patch(a):
    a.output.mkdir(exist_ok=False);ctx,B,cold,Z=old(a)
    anchor=form_anchor(ctx,B,cold,ctx.doc['saved_retention']['theta'],'A')
    # Select actual dual-supported locations/vectors, not geometry-number guesses.
    if a.proof_source:
        proof=[load(a.proof_source/(label+'-actual')/'independent-replay.json') for label in ('A','B')]
        geometry=np.empty(ctx.N,int);members=ctx.arrays['members'];geometry[members[:,1]]=members[:,5]
        for p in proof:
            for w in p['limiting_constraints']:
                if 'point' in w:w.update(chart=ctx.block('chart')[w['point']].tolist(),geometry=int(geometry[w['point']]))
    else:proof=[p for p in load(a.previous/'compact/obstruction.json')['replays'] if p['anchor']=='A' and not p['buffered']]
    groups=[]
    for pr in proof:
        group=[];geometries=set()
        for w in sorted(pr['limiting_constraints'],key=lambda r:abs(r['dual_multiplier']),reverse=True):
            if w['kind']!='PSD' or w['inner'] or w['geometry'] in geometries:continue
            geometries.add(w['geometry']);group.append(w)
        groups.append(group)
    selected=[];seen=set()
    for j in range(max(map(len,groups))):
        for group in groups:
            if j>=len(group):continue
            w=group[j];key=(w['geometry'],tuple(w['chart']),tuple(np.round(w['vector'],6)))
            if key in seen:continue
            seen.add(key);selected.append(dict(chart=w['chart'],geometry=w['geometry'],eigenvector=w['vector']))
            if len(selected)==a.target_count:break
        if len(selected)==a.target_count:break
    geometries=[int(x) for x in a.geometries.split(',')]
    xml=case(a.output/'case',a.h1,load(a.preflight),flat=False);folder=a.output/'arrays';folder.mkdir()
    d,rr=native(a,xml,dict(operation='h2k_patch',state=anchor.tolist(),state_sha256=state_hash(anchor),base_sha256=COLD_HASH,patch_level=1,patch_geometries=geometries,witnesses=selected,context_directory=str(folder)),a.output/'native')
    atomic(a.output/'patch.json',dict(context=d['context'],receipt=rr,requested_geometries=geometries,targets=selected,source_context_sha256=sha(ctx.path)))


def construct(a):
    begun=time.monotonic();a.output.mkdir(exist_ok=False);ctx,B,cold,Z=old(a);nv=len(Z)-16
    doc=load(a.patch/'patch.json');meta=doc['context'];P=arrays(meta,a.patch/'arrays')
    rows=np.array(meta['coefficient_rows'],int);free=np.array(ctx.meta['algebra']['free_velocity_rows'],int)
    index={int(i):j for j,i in enumerate(free)};ix=[index[int(i)] for i in rows]
    M=np.array(P['mass']);W=linalg.solve_triangular(np.linalg.cholesky(M).T,np.eye(len(rows)),lower=False)
    Q=np.array(P['Q'])@W;qn=np.linalg.norm(Q,axis=1);Q=Q[qn>1e-12]/qn[qn>1e-12,None]
    _,sq,Vq=linalg.svd(Q,full_matrices=True);rq=int(np.count_nonzero(sq>sq[0]*1e-10));Nq=W@Vq[rq:].T
    actions=np.vstack([P['divergence_action_GL'+str(o)] for o in (16,20,24)])
    C=actions@Nq;_,sd,V=linalg.svd(C,full_matrices=False)
    rd=int(np.count_nonzero(sd>max(1e-12,sd[0]*1e-10))) if len(sd) else 0
    Ns=Nq@V[rd:].T
    F=Z[:nv,free];fullM=ctx.algebra['mass_free'];G=F@fullM@F.T;added=[];history=[]
    massCross=F@fullM[:,ix]
    # Restrict within the divergence-null space before removing old overlap.
    # Subtracting a weak-only old direction afterwards would destroy control.
    if Ns.shape[1]:
        H=massCross@Ns;_,ss,Vs=linalg.svd(H,full_matrices=True);rh=int(np.count_nonzero(ss>max(1e-12,ss[0]*1e-10))) if len(ss) else 0;Ns=Ns@Vs[rh:].T
    candidates=[]
    for j,g in enumerate(P['Rayleigh_gradients']):
        if Ns.shape[1]:candidates.append((j,'exact_strong_null',Ns@(Ns.T@g)))
    if a.compensated:
        # The weak-null responses retain their actual divergence. They require
        # explicit finite compensation/cap checks; a penalty is not a gate.
        tangent=np.vstack([P['base_divergence_GL'+str(o)][0]@P['divergence_action_GL'+str(o)]@Nq for o in (16,20,24)]+list(massCross@Nq))
        norms=np.linalg.norm(tangent,axis=1);tangent=tangent[norms>1e-12]/norms[norms>1e-12,None]
        _,st,Vt=linalg.svd(tangent,full_matrices=True);rt=int(np.count_nonzero(st>st[0]*1e-10));Nt=Nq@Vt[rt:].T
        # These constraints cancel the linear divergence change. The positive
        # direction-square term remains and must pass the finite native gate.
        Dt=(actions@Nt).T@(actions@Nt)
        for eta in map(float,a.weights.split(',')):
            A=np.eye(Nt.shape[1])+eta*Dt
            for j,g in enumerate(P['Rayleigh_gradients']):candidates.append((j,'weak_null_tangent_compensated_weight_'+str(eta),Nt@linalg.solve(A,Nt.T@g,assume_a='pos')))
        if a.divergence_compensator:
            _,so,Vo=linalg.svd(massCross@Nq,full_matrices=True);ro=int(np.count_nonzero(so>max(1e-12,so[0]*1e-10)));N0=Nq@Vo[ro:].T
            D0=actions@N0;h=sum(P['divergence_action_GL'+str(o)].T@P['base_divergence_GL'+str(o)][0] for o in (16,20,24))
            candidates.append((-1,'explicit_divergence_reducing_compensator',-N0@linalg.solve(np.eye(N0.shape[1])+D0.T@D0,N0.T@h,assume_a='pos')))
    for target,method,u in candidates:
        if nv+len(added)>=128:break
        charge(a,'constrained response',attempts=1)
        norm=float(np.sqrt(max(0,u@M@u)));entry=dict(target=target,method=method,unscaled_mass_norm=norm)
        if norm<1e-10:entry['retained']=False;history.append(entry);continue
        u=u/norm;v=np.zeros(len(B));v[rows]=u;current=np.vstack((Z[:nv],added)) if added else Z[:nv];Fc=current[:,free]
        cross=Fc@fullM[:,ix]@u;GG=Fc@fullM@Fc.T;aug=np.block([[GG,cross[:,None]],[cross[None,:],np.array([[1.]])]])
        ev=np.linalg.eigvalsh(aug);entry.update(Gram_extremes=[float(ev[0]),float(ev[-1])],Q_max=float(np.max(abs(P['Q']@u))),divergence_norms=[float(np.linalg.norm(P['divergence_action_GL'+str(o)]@u)) for o in (16,20,24)],base_divergence_cross=[float(P['base_divergence_GL'+str(o)][0]@(P['divergence_action_GL'+str(o)]@u)) for o in (16,20,24)])
        if ev[0]>1e-10*ev[-1] and entry['Q_max']<1e-9:
            added.append(v);entry.update(retained=True,index=nv+len(added)-1,Rayleigh_derivatives=(P['Rayleigh_gradients']@u).tolist())
        else:entry['retained']=False
        history.append(entry)
    if not added:
        atomic(a.output/'construction.json',dict(version='H2kNoUsefulNativeResponseV1',patch_coefficients=len(rows),weak_constraint_rank=rq,weak_nullity=Nq.shape[1],strong_action_rank=rd,strong_nullity=V.shape[0]-rd,strong_null_old_orthogonal_rank=Ns.shape[1],directions=history,numerical_seconds=time.monotonic()-begun))
        charge(a,'no useful response',construction_seconds=time.monotonic()-begun)
        return
    Znew=np.vstack((Z[:nv],added,Z[nv:]));file=a.output/'basis.f64';Znew.astype('<f8').tofile(file)
    report=dict(version='H2kNativeConstructionV1',old_velocity_count=nv,new_velocity_count=len(Znew)-16,added_count=len(added),patch_file=str(a.patch/'patch.json'),patch_geometries=meta['geometries'],patch_coefficients=len(rows),weak_constraint_rank=rq,weak_nullity=Nq.shape[1],strong_action_rank=rd,strong_nullity=V.shape[0]-rd,strong_null_old_orthogonal_rank=Ns.shape[1],strong_singular_values=sd.tolist(),directions=history,basis_file=str(file),basis_sha256=sha(file),old_basis_sha256=ctx.doc['basis_sha256'],old_columns_bits_preserved=Znew[:nv].tobytes()==Z[:nv].tobytes(),scalar_basis_bits_preserved=Znew[-16:].tobytes()==Z[-16:].tobytes(),numerical_seconds=time.monotonic()-begun)
    atomic(a.output/'construction.json',report);charge(a,'response construction complete',construction_seconds=report['numerical_seconds'])


def context(a):
    a.output.mkdir(exist_ok=False);ctx,B,cold,Z=old(a);c=load(a.construction/'construction.json');meta=load(Path(c['patch_file']))['context'];file=Path(c['basis_file']);nv=c['new_velocity_count'];basis=np.memmap(file,dtype='<f8',mode='r',shape=(nv+16,len(B)))
    # A target with nonzero outer derivative; choose a bounded finite amplitude.
    selected=next(d for d in c['directions'] if d['retained'] and d['target']>=0 and not meta['witnesses'][d['target']]['inner'] and abs(d['Rayleigh_derivatives'][d['target']])>1e-6)
    j=selected['index'];derivative=selected['Rayleigh_derivatives'][selected['target']];amplitude=min(.0001,4e-4*ctx.meta['model_specification']['K0']/abs(derivative))
    plus=np.zeros(nv);plus[j]=amplitude*np.sign(derivative)
    compensators=[x for x in c['directions'] if x.get('retained') and x['target']==-1]
    if compensators:plus[compensators[-1]['index']]=min(1e-9,.01*amplitude)
    trials=np.array([np.zeros(nv),plus,-plus]);folder=a.output/'arrays';folder.mkdir();xml=case(a.output/'case',a.h1,load(a.preflight),flat=False)
    d,rr=native(a,xml,dict(operation='h2k_context',state=B.tolist(),state_sha256=B_HASH,base_sha256=COLD_HASH,basis_file=str(file),basis_sha256=sha(file),velocity_direction_count=nv,B_scalar_coefficients=load(a.parent/'standard-B/best-candidate.json')['conditioned_scalar_coefficients'],velocity_verification=trials.tolist(),context_directory=str(folder)),a.output/'native')
    d.update(array_directory=str(folder),basis_file=str(file),basis_sha256=sha(file),receipt=rr,runtime=runtime_identity(a.probe),conditioning=dict(version='RawSupportedFullGramV1',velocity=f'{len(Z)-16} prior columns unchanged; native constrained additions and full physical Gram',rank=nv),saved_retention=ctx.doc['saved_retention'],construction_file=str(a.construction/'construction.json'))
    cf=a.output/'context.json';atomic(cf,d);new=Context(cf)
    for key in ('geometry_identity','inventory_identity','reference_sha256','defect_sha256','patch_signature'):require(ctx.doc[key]==new.doc[key],'physical identity '+key)
    for name in ctx.blocks:
        oldb=ctx.block(name);newb=new.block(name)
        if name=='velocity_basis':newb=newb[:,:len(Z)-16]
        require(oldb.tobytes()==newb.tobytes(),'original response block changed '+name)
    for name in ('members','member_weights'):require(ctx.arrays[name].tobytes()==new.arrays[name].tobytes(),'membership changed')
    P=new.problem('A',d['saved_retention']['theta']);budget=Budget(0,12,120);base=P.scan(np.zeros(P.dim),1e-4*P.K0,budget);w=meta['witnesses'][selected['target']];ids=np.flatnonzero(np.all(new.block('chart')==w['chart'],axis=1)&(P.geometry==w['geometry'])&P.cov);require(len(ids)>0,'target not in declared inventory');i=int(ids[0]);v=w['eigenvector'];row,_=P.rayleigh_cut(i,v,1e-4*P.K0);results=[]
    for sign in (1.,-1.):
        x=np.r_[sign*plus,np.zeros(8)];s=P.scan(x,1e-4*P.K0,budget);g=s['gradient'][i]-base['gradient'][i];nt=P.nt[i] if P.inner[i] else .09*P.k[i]**2/base['epsilon'][i];gain=-2*nt*(v[0]**2*g[0]+v[0]*v[2]*(g[1]+g[2])+v[2]**2*g[3]);results.append(dict(sign=sign,amplitude=amplitude,gain=float(gain),predicted=-float(row[:nv]@(sign*plus))*P.K0*(P.ec[i]/base['epsilon'][i] if not P.inner[i] else 1),summary=summary(s)))
    useful=any(z['gain']>=1e-4*P.K0 and z['summary']['hard_violation']<=1e-9 for z in results)
    verification=d['context']['verification'];mass_error=max(abs(results[t]['summary']['mass_squared']-verification['meaningful_trial_mass_squared'][t+1]) for t in range(2));div_error=max(abs(results[t]['summary']['strong_divergence_squared'][0]-verification['meaningful_trial_divergence_squared'][t+1]) for t in range(2))
    require(mass_error<1e-12 and div_error<1e-12,'finite native physical quadratic actions differ')
    atomic(a.output/'finite-control.json',dict(passed=useful,selected_direction=j,target=w,trials=results,base=summary(base),native_mass_squared_error=mass_error,native_GL16_divergence_squared_error=div_error,compensator_coefficients={str(x['index']):float(plus[x['index']]) for x in compensators},old_response_blocks_bitwise=True,old_memberships_bitwise=True,context_sha256=sha(cf)))
    charge(a,'context finite control',scans=budget.scans,construction_seconds=budget.report()['seconds'])
    require(useful,'no useful finite controlled response; do not optimize unchanged/inadequate control')
    atomic(a.output/'frozen.json',dict(context=str(cf),context_sha256=sha(cf),basis_sha256=sha(file),velocity_count=nv,finite_control_passed=True))


def search(a):
    a.output.mkdir(exist_ok=False);ctx=Context(a.context);_,B,cold,_=old(a);nv=ctx.meta['direction_counts']['velocity'];Z=np.memmap(ctx.doc['basis_file'],dtype='<f8',mode='r',shape=(nv+16,len(B)));theta=ctx.doc['saved_retention']['theta']
    if a.resume_frozen:
        prior=load(a.resume_frozen/'result.json');require(prior['context_sha256']==sha(ctx.path) and prior['classification']=='inconclusive','only identical inconclusive frozen context can resume numerical organization')
    ledger=charge(a,'equivalent LP numerical repair in frozen epoch' if a.resume_frozen else 'frozen search epoch',epochs=0 if a.resume_frozen else 1);t=ledger['totals'];budget=Budget(240-t.get('rounds',0),720-t.get('scans',0),3600-t.get('search_seconds',0));out=dict(version='H2kJointSearchV1',velocity_count=nv,epsilon_count=8,retention=ctx.doc['saved_retention'],phases=[],classification='inconclusive',context_sha256=sha(ctx.path),created_Asia_Shanghai=now(),resumed_frozen=str(a.resume_frozen) if a.resume_frozen else None);count=0
    try:
        for label in ('A','B'):
            require(label=='A' or out.get('A_actual_obstruction'),'B requires independently replayed actual A obstruction')
            P=ctx.problem(label,theta);anchor=form_anchor(ctx,B,cold,theta,label)
            for buffered in (True,False):
                tag=label+('-buffered' if buffered else '-actual');folder=a.output/tag;folder.mkdir()
                def save(event):
                    nonlocal count
                    count+=1;atomic(folder/('event-%04d.json'%count),event)
                    if event['event']=='best_hard_diagnostic':atomic(folder/'best-candidate.json',make_candidate(ctx,anchor,B,Z,event['coefficients'],label,theta))
                initial=load(a.initial_candidate)['joint_coefficients'] if a.initial_candidate else None
                r=(exact_caps_phase if a.exact_caps else phase)(P,(2 if buffered else 1)*1e-4*P.K0,[.501,1.999] if buffered else [.5,2],budget,save,initial=initial,lp_solver=equilibrated_lp if a.equilibrated else None);ph=dict(anchor=label,buffered=buffered,result=r);out['phases'].append(ph)
                if r['best_coefficients'] is not None:
                    candidate=make_candidate(ctx,anchor,B,Z,r['best_coefficients'],label,theta);candidate.update(optimizer_summary=r['best_summary'],phase=tag);atomic(a.output/'candidate.json',candidate)
                if r['termination']=='seen_feasible':out.update(classification='seen_feasible_candidate',selected_anchor=label);break
                if r['termination']!='verified_relaxation_obstruction':out.update(classification='inconclusive',termination=r['termination']);break
                initial=P.scan(np.zeros(P.dim),r['delta'],budget);proof=replay(P,ph,initial);atomic(folder/'independent-replay.json',proof)
                if not buffered:out[label+'_actual_obstruction']=proof['verified']
                atomic(a.output/'result.json',out)
            if out['classification']=='seen_feasible_candidate' or not out.get(label+'_actual_obstruction'):break
            if label=='B':out['classification']='verified_both_anchors_obstruction'
    except BudgetExhausted as e:out.update(termination=str(e),classification='inconclusive')
    finally:
        out.update(budget=budget.report(),peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,physical_attempts=0);atomic(a.output/'result.json',out);charge(a,'search and independent replay',rounds=budget.rounds,scans=budget.scans,search_seconds=budget.report()['seconds'])


def reserve(a):
    """One optional stricter construction solve, never a physical obstruction."""
    a.output.mkdir(exist_ok=False);ctx=Context(a.context);c=load(a.candidate);_,B,cold,_=old(a);P=ctx.problem(c['anchor'],c['theta']);P.div_limits=np.square([.06030809942199109-1e-9,.06030809942199109+1e-8-1e-9,.06030809942199109+1e-8-1e-9])
    nv=P.nv;Z=np.memmap(ctx.doc['basis_file'],dtype='<f8',mode='r',shape=(nv+16,len(B)));t=charge(a,'startup reserve construction')['totals'];budget=Budget(240-t.get('rounds',0),720-t.get('scans',0),3600-t.get('search_seconds',0));count=0
    def save(event):
        nonlocal count
        count+=1;atomic(a.output/('event-%04d.json'%count),event)
    r=(exact_caps_phase if a.exact_caps else phase)(P,20e-4*P.K0,[.501,1.999],budget,save,initial=c['joint_coefficients'],lp_solver=equilibrated_lp);chosen=c
    if r['best_coefficients'] is not None and min(r['best_summary']['inner_min'],r['best_summary']['outer_min'])>min(c['optimizer_summary']['inner_min'],c['optimizer_summary']['outer_min']):
        anchor=form_anchor(ctx,B,cold,c['theta'],c['anchor']);chosen=make_candidate(ctx,anchor,B,Z,r['best_coefficients'],c['anchor'],c['theta']);chosen.update(optimizer_summary=r['best_summary'],phase='startup_reserve_20delta')
    atomic(a.output/'candidate.json',chosen);atomic(a.output/'result.json',dict(version='H2kStartupReserveV1',scope='optional stricter numerical construction only; acceptance remains delta and original divergence limits; no physical obstruction inference',result=r,selected_state_sha256=chosen['state_sha256'],changed=chosen['state_sha256']!=c['state_sha256'],search_divergence_limits=np.sqrt(P.div_limits).tolist()))
    charge(a,'startup reserve work',rounds=budget.rounds,scans=budget.scans,search_seconds=budget.report()['seconds'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('patch','construct','context','search','reserve'))
    for name in ('root','parent','previous','probe','output','patch','construction','context','candidate'):p.add_argument('--'+name,type=Path)
    p.add_argument('--geometries');p.add_argument('--compensated',action='store_true');p.add_argument('--divergence-compensator',action='store_true');p.add_argument('--weights',default='.01,1');p.add_argument('--target-count',type=int,default=24);p.add_argument('--retained-context',type=Path);p.add_argument('--proof-source',type=Path);p.add_argument('--resume-frozen',type=Path);p.add_argument('--initial-candidate',type=Path);p.add_argument('--equilibrated',action='store_true');p.add_argument('--exact-caps',action='store_true');a=p.parse_args();require(1<=a.target_count<=32,'explicit target count');a.preflight=a.parent/'comparison.json';a.h1=Path('/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111/compatibility-resolved');os.environ.update(THREADS);globals()[a.phase](a)


if __name__=='__main__':main()
