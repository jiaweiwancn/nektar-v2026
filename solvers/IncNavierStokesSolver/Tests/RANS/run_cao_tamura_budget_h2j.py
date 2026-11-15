#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded native support/divergence enrichment and positively tested V2 handoff."""
import argparse,os,time,resource
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy import linalg
from cao_tamura_h2d_state import load
from cao_tamura_h2i_convex import Context
from cao_tamura_h2h_scalar import Budget
from cao_tamura_initializer_context import NativeContext
from run_cao_tamura_budget_h2i import form_anchor,B_HASH,COLD_HASH
from cao_tamura_budget_h2g import case
from diagnose_cls_accepted_covariance import run_probe
from run_cao_tamura_budget_h2h import runtime_identity
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import sha,state_hash,require,THREADS
from cao_tamura_startup_receipt import complete_receipt,validate_receipt
from run_cao_tamura_budget_h2e import certificate_v2
from run_cao_tamura_budget_h2f_startup import StartupRun

def defaults(a):
    a.preflight=a.parent/'comparison.json'
    a.h1=Path('/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h1-20261111/compatibility-resolved')

def native(a,xml,packet,folder):
    packet.update(xml_sha256=sha(xml),expected_model='StandardKEpsilon',inner_stress_model='Linear')
    _,rr=run_probe(a.probe,xml,packet,folder);return load(Path(rr['output'])),rr

def arrays(meta,root):
    out={}
    for name,e in meta['files'].items():
        p=root/e['file'];require(sha(p)==e['sha256'],'patch array hash');out[name]=np.memmap(p,dtype=e['dtype'],mode='r',shape=tuple(e['shape']))
    return out

def prepare(a):
    start=time.monotonic();a.output.mkdir(exist_ok=False);old=Context(a.previous/'joint-context/context.json');B=np.array(load(a.parent/'standard-B/best-candidate.json')['state']);cold=np.array(NativeContext(a.parent/'standard-context/context.json').array('base')[0]);require(state_hash(B)==B_HASH and state_hash(cold)==COLD_HASH,'H2j immutable state hashes')
    saved=load(a.previous/'joint-search/result.json')['retention'];anchor=form_anchor(old,B,cold,saved['theta'],'A');require(state_hash(anchor)==load(a.previous/'joint-search/anchors.json')['anchors']['A']['state_sha256'],'H2i exact anchor changed')
    Zold=np.array(np.memmap(old.doc['basis_file'],dtype='<f8',mode='r',shape=(57,len(B))));require(sha(Path(old.doc['basis_file']))==old.doc['basis_sha256'],'old conditioned basis changed')
    fullM=old.algebra['mass_free'];free=old.meta['algebra']['free_velocity_rows'];ix={int(j):i for i,j in enumerate(free)};oldFree=Zold[:41,free];attempts=0;numerical=0.;history=[]
    for level in (1,2):
        folder=a.output/('level'+str(level));folder.mkdir();datafolder=folder/'arrays';datafolder.mkdir();xml=case(folder/'case',a.h1,load(a.preflight),flat=False)
        d,rr=native(a,xml,dict(operation='h2j_patch',state=anchor.tolist(),state_sha256=state_hash(anchor),base_sha256=COLD_HASH,patch_level=level,context_directory=str(datafolder)),folder/'native');meta=d['context'];P=arrays(meta,datafolder);rows=np.array(meta['coefficient_rows'],int);indices=np.array([ix[int(i)] for i in rows]);tic=time.monotonic();M=np.array(P['mass']);L=np.linalg.cholesky(M);W=linalg.solve_triangular(L.T,np.eye(len(rows)),lower=False)
        pieces=[np.array(P['Q'])]+[np.array(P['divergence_action_GL'+str(o)]) for o in (16,20,24)];C=np.vstack(pieces)@W;norm=np.linalg.norm(C,axis=1);C=C[norm>1e-12]/norm[norm>1e-12,None];_,sv,V=linalg.svd(C,full_matrices=False);rank=int(np.count_nonzero(sv>sv[0]*1e-10));N=W@V[rank:].T;new=[];construction=[]
        for j,g in enumerate(P['Rayleigh_gradients']):
            attempts+=1;require(attempts<=12,'construction attempts exhausted');u=N@(N.T@g);size=float(np.sqrt(max(0.,u@M@u)))
            if size<1e-12:construction.append(dict(target=j,reason='no strong-divergence-null response'));continue
            u/=size;v=np.zeros(len(B));v[rows]=u;current=np.vstack((Zold[:41],new)) if len(new) else Zold[:41];F=current[:,free];mf=fullM[:,indices]@u;cross=F@mf;G=F@fullM@F.T;aug=np.block([[G,cross[:,None]],[cross[None,:],np.ones((1,1))]]);eigen=np.linalg.eigvalsh(aug);rank_ok=eigen[0]>1e-10*eigen[-1]
            entry=dict(target=j,mass_norm=1.,Rayleigh_derivative=float(g@u),rank_ok=bool(rank_ok),Gram_min_max=[float(eigen[0]),float(eigen[-1])],Q_max=float(np.max(abs(P['Q']@u))),strong_direction_norms={str(o):float(np.linalg.norm(P['divergence_action_GL'+str(o)]@u)) for o in (16,20,24)})
            require(entry['Q_max']<1e-9 and max(entry['strong_direction_norms'].values())<1e-8,'native response not divergence controlled')
            if rank_ok and len(new)<6:new.append(v);entry['added_index']=40+len(new)
            construction.append(entry)
        numerical+=time.monotonic()-tic;require(numerical<=600,'numerical construction time exhausted');summary=dict(level=level,patch=meta,rank=rank,nullity=N.shape[1],singular_values=sv.tolist(),attempts=attempts,numerical_seconds=numerical,responses=construction,old_columns_unchanged=True,conditioning='qualified raw mass-normalized responses retained; full Gram used, no subtraction');atomic(folder/'construction.json',summary);history.append(summary)
        if not new:continue
        # Choose the H-plus response when present; no added quota or scalar enrichment.
        chosen=next((e for e in construction if e.get('target')==1 and 'added_index' in e),next(e for e in construction if 'added_index' in e));nv=41+len(new);Z=np.vstack((Zold[:41],new,Zold[41:]));basis=folder/'basis.f64';Z.astype('<f8').tofile(basis);alpha=min(1e-4,max(1e-6,4e-4*old.meta['model_specification']['K0']/abs(chosen['Rayleigh_derivative'])));positive=np.zeros(nv);positive[chosen['added_index']]=alpha;trials=np.array([np.zeros(nv),positive,-positive]);contextdir=folder/'context-arrays';contextdir.mkdir()
        d2,rr2=native(a,xml,dict(operation='h2j_context',state=B.tolist(),state_sha256=B_HASH,base_sha256=COLD_HASH,basis_file=str(basis),basis_sha256=sha(basis),velocity_direction_count=nv,B_scalar_coefficients=load(a.parent/'standard-B/best-candidate.json')['conditioned_scalar_coefficients'],velocity_verification=trials.tolist(),context_directory=str(contextdir)),folder/'context-native')
        d2.update(array_directory=str(contextdir),basis_file=str(basis),basis_sha256=sha(basis),receipt=rr2,runtime=runtime_identity(a.probe),conditioning=dict(velocity='old 41 unchanged plus native supported strong-divergence-null raw responses',rank=nv),saved_retention=saved,construction_file=str(folder/'construction.json'));contextfile=folder/'context.json';atomic(contextfile,d2);ctx=Context(contextfile)
        # Equality is checked on all old response columns and memberships, not just a tuple.
        for key in ('geometry_identity','inventory_identity','reference_sha256','defect_sha256','patch_signature'):require(old.doc[key]==ctx.doc[key],'old physical identity changed '+key)
        require(old.block('velocity_basis').tobytes()==ctx.block('velocity_basis')[:,:41].tobytes(),'old velocity responses changed')
        for key in ('chart','physical','gradient_B','velocity_B','velocity_cold','k_B','k_cold','epsilon_B','epsilon_cold','k_basis','epsilon_basis','inner_nu','dependent','inner','covariance','weight'):
            require(old.block(key).tobytes()==ctx.block(key).tobytes(),'old response bits changed '+key)
        require(old.arrays['members'].tobytes()==ctx.arrays['members'].tobytes(),'inventory memberships changed');require(old.arrays['member_weights'].tobytes()==ctx.arrays['member_weights'].tobytes(),'inventory weights changed')
        tic=time.monotonic();problem=ctx.problem('A',saved['theta']);local=[];budget=Budget(0,8,120);base=problem.scan(np.zeros(nv+8),1e-4*problem.K0,budget)
        for sign in (1.,-1.):
            x=np.r_[sign*positive,np.zeros(8)];scan=problem.scan(x,1e-4*problem.K0,budget);w=meta['witnesses'][chosen['target']];mask=np.all(ctx.block('chart')==w['chart'],axis=1)&(problem.geometry==w['geometry'])&problem.cov;ids=np.flatnonzero(mask);i=int(ids[0]);v=np.array(w['eigenvector']);g=scan['gradient'][i];gb=base['gradient'][i];nu=problem.nt[i] if problem.inner[i] else .09*problem.k[i]**2/scan['epsilon'][i];gain=-2*nu*(v[0]**2*(g[0]-gb[0])+v[0]*v[2]*(g[1]+g[2]-gb[1]-gb[2])+v[2]**2*(g[3]-gb[3]));pred=chosen['Rayleigh_derivative']*sign*alpha
            local.append(dict(sign=sign,amplitude=alpha,Rayleigh_gain=float(gain),predicted_gain=float(pred),covariance_before=float(base['lambda_covariance'][i]),covariance_after=float(scan['lambda_covariance'][i]),mass_L2=float(np.sqrt(scan['mass_squared'])),divergence_L2=np.sqrt(scan['strong_divergence_squared']).tolist(),sampled_cap=scan['velocity_max'],hard_violation=scan['hard_violation']))
        numerical+=time.monotonic()-tic;useful=any(x['Rayleigh_gain']>=1e-4*problem.K0 and x['hard_violation']<=1e-9 for x in local);require(numerical<=600,'construction time exhausted');control=dict(useful=useful,trials=local,old_response_bits_unchanged=True,construction_scans=budget.scans,aggregate_numerical_seconds=numerical);atomic(folder/'finite-control.json',control)
        if useful:
            frozen=dict(version='H2jFrozenAugmentedBasisV1',context=str(contextfile),context_sha256=sha(contextfile),basis_file=str(basis),basis_sha256=sha(basis),velocity_count=nv,added_count=len(new),selected_level=level,attempts=attempts,numerical_seconds=numerical,total_setup_seconds=time.monotonic()-start,finite_control=control,source_context_sha256=sha(old.path),scalar_basis_unchanged=True);atomic(a.output/'frozen.json',frozen);return
    atomic(a.output/'frozen.json',dict(version='H2jNoUsefulEnrichmentV1',useful=False,attempts=attempts,numerical_seconds=numerical,history=history))

def handoff(a):
    a.output.mkdir(exist_ok=False);old=a.parent/'standard-flat-final';candidate=load(old/'changed-candidate.json');xml=old/'candidate/case.xml';d,rr=native(a,xml,dict(operation='h2e_candidate',state=candidate['state'],state_sha256=candidate['state_sha256'],base_sha256=candidate['base_sha256'],coupled=True),a.output/'native')
    model=load(Path(rr['output']+'.model.json'));solver=runtime_identity(a.solver);probe=runtime_identity(a.probe);receipt=complete_receipt(Path(rr['output']),model,candidate['divergence_limit'],probe,solver);validate_receipt(receipt);receiptfile=a.output/'issuance-receipt.json';atomic(receiptfile,receipt)
    src=Path(__file__).resolve().parents[4];names=['solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp','library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp','solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_startup_receipt.py'];prov=dict(K0=model['K0'],divergence_limit=candidate['divergence_limit'],basis_sha256=sha(old/'changed-candidate.json'),context_sha256=sha(old/'changed-candidate.json'),native_receipt_file=str(receiptfile),source_revision='H2j native closure and exact source hashes',source_hashes={n:sha(src/n) for n in names});cert=certificate_v2(d,receipt,d['fixed_reference'],solver,probe,prov);cert['model_specification']=model;certpath=a.output/'certificate.json';atomic(certpath,cert)
    builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=True)
    run=StartupRun(a.solver,certpath,a.preflight,a.h1,a.output/'import',flat=True,builder=builder);h0=run.checked_native('zero',2e-6);require(state_hash(h0['state'])==candidate['state_sha256'] and state_hash(h0['state'])!=state_hash(h0['fixed_reference']),'changed candidate/reference import distinction');require(h0['initialization_lift_calls_iterations']==[0,0],'V2 import used Lift/projection');require(state_hash(h0['fixed_reference'])==cert['fixed_reference_sha256'] and state_hash(h0['fixed_defect'])==cert['fixed_defect_sha256'],'import reference/defect differs');atomic(a.output/'ledger.json',dict(case='independent affine standard fixture',passed=True,state_sha256=candidate['state_sha256'],base_sha256=candidate['base_sha256'],hill_GL28_consumed=False,steps=0,complete_receipt=str(receiptfile),native_import=run.out))

def flat(a):
    a.output.mkdir(exist_ok=False);results={}
    for model in ('cls','standard'):
        xml=case(a.output/model,a.h1,load(a.preflight),flat=True,model=model);packet=dict(operation='h2j_flat',xml_sha256=sha(xml),inner_stress_model='Linear',expected_model='CLSHighRe' if model=='cls' else 'StandardKEpsilon');_,rr=run_probe(a.probe,xml,packet,a.output/(model+'-probe'));results[model]=dict(receipt=rr,output=rr['output'])
    atomic(a.output/'results.json',results)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('prepare','handoff','flat'))
    for n in ('parent','previous','probe','solver','output'):p.add_argument('--'+n,type=Path)
    a=p.parse_args();defaults(a);os.environ.update(THREADS);os.environ['OMP_THREAD_LIMIT']='1';globals()[a.phase](a)
if __name__=='__main__':main()
