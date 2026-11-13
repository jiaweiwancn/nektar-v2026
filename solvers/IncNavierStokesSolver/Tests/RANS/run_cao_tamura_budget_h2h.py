#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Bounded fixed-B scalar feasibility, one fresh audit, conditional V2 startup."""
import argparse,os,time,resource,struct
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
from cao_tamura_h2d_state import load
from cao_tamura_initializer_context import NativeContext
from cao_tamura_budget_h2g import case
from cao_tamura_h2h_scalar import Context,Problem,Budget,solve,BLOCKS
from continue_cls_epawf_pimple import atomic,identity
from diagnose_cls_accepted_covariance import run_probe
from run_cls_epawf_picard import sha,state_hash,require,THREADS
from run_cao_tamura_budget_h2e import certificate_v2
from run_cao_tamura_budget_h2f_startup import StartupRun

B_HASH='c83da50b65d2e9980325d365779a9d53394fbf10e6431d1b87f79728366ae13d'
COLD_HASH='038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c'
def now():return datetime.now(ZoneInfo('Asia/Shanghai')).isoformat()
def bytes_equal(a,b):return np.asarray(a,dtype='<f8').tobytes()==np.asarray(b,dtype='<f8').tobytes()
def require_fresh(receipt,state,base,model):
 require(receipt and receipt.get('quadrature_order')==24,'H2h certificate requires fresh GL24 receipt')
 require(receipt.get('state_sha256')==state and receipt.get('cold_base_sha256')==base,'H2h GL24 candidate/base mismatch')
 require(receipt.get('model_specification')==model,'H2h GL24 model mismatch')
 require(receipt.get('minimum_covariance',-float('inf'))>=1e-4*model['K0']-1e-12,'H2h GL24 positive covariance margin required')
 require(receipt.get('scalar_positive_finite') and receipt.get('scalar_envelopes_passed') and receipt.get('integrals_complete'),'H2h GL24 scalar gate')

def prepare(a):
 start=time.monotonic();old=NativeContext(a.parent/'standard-context/context.json');B=load(a.parent/'standard-B/best-candidate.json');require(state_hash(B['state'])==B_HASH and old.meta['base_sha256']==COLD_HASH and B['base_sha256']==COLD_HASH,'H2g B/cold identity')
 _,Ts,conditioning=old.transforms();nv=old.counts['velocity'];raw=np.array(old.array('basis')[nv:]);require(raw.shape==(16,2804),'saved scalar basis size');Z=Ts.T@raw
 a.output.mkdir(exist_ok=False);file=a.output/'scalar-basis.f64';Z.astype('<f8').tofile(file);np.asarray(Ts,dtype='<f8').tofile(a.output/'scalar-transform.f64')
 cold=np.array(old.array('base')[0]);krows=old.meta['rows']['k'];erows=old.meta['rows']['epsilon'];allowed=set(krows+erows)
 for i in range(2804):
  if i not in allowed:require(np.all(Z[:,i]==0),'saved scalar basis changes protected row')
 zero=np.array(B['state'],dtype='<f8');require(bytes_equal(zero,B['state']),'zero-increment B copy changed bits')
 xml=case(a.output/'case',a.h1,load(a.preflight),flat=False);arrays=a.output/'arrays';arrays.mkdir()
 packet=dict(operation='h2h_scalar_context',inner_stress_model='Linear',xml_sha256=sha(xml),expected_model='StandardKEpsilon',state=B['state'],state_sha256=B_HASH,base_sha256=COLD_HASH,scalar_basis_file=str(file),scalar_basis_sha256=sha(file),B_scalar_coefficients=B['conditioned_scalar_coefficients'],context_directory=str(arrays))
 d,r=run_probe(a.probe,xml,packet,a.output/'native-export');d=load(Path(r['output']));d.update(array_directory=str(arrays),basis_file=str(file),basis_sha256=sha(file),conditioning=conditioning,transform_file=str(a.output/'scalar-transform.f64'),transform_sha256=sha(a.output/'scalar-transform.f64'),parent_context=str(old.path),parent_context_sha256=sha(old.path),parent_B_file=str(a.parent/'standard-B/best-candidate.json'),parent_B_file_sha256=sha(a.parent/'standard-B/best-candidate.json'),receipt=r,probe_runtime=identity(a.probe),observation_blocks={k:dict(offset=v[0],width=v[1]) for k,v in BLOCKS.items()},created_Asia_Shanghai=now(),setup_wall_seconds=time.monotonic()-start,peak_RSS_KiB_cumulative=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
 require(d['context']['state_sha256']==B_HASH and d['context']['base_sha256']==COLD_HASH,'native context parent mismatch');atomic(a.output/'context.json',d)

def search(a):
 setup=time.monotonic();ctx=Context(a.context);B=load(a.parent/'standard-B/best-candidate.json');pre=load(a.preflight);K0=pre['models'][1]['K0'];P=ctx.problem();Z=np.memmap(ctx.doc['basis_file'],dtype='<f8',mode='r',shape=(16,2804));require(sha(Path(ctx.doc['basis_file']))==ctx.doc['basis_sha256'],'saved scalar basis changed');rows=ctx.meta['scalar_rows'];allowed=np.array(sorted(rows['k']+rows['epsilon']),int)
 a.output.mkdir(exist_ok=False);budget=Budget();setup_seconds=time.monotonic()-setup
 # Implied coefficient box and native centring bounds quantify the floating
 # affine interpretation, without reconstructing velocity or pressure.
 box,box_info=P.coefficient_box();outer=ctx.block('inner')==0;err=ctx.block('basis_error')[outer];bound=np.full(16,box);errors=ctx.block('k_error')[outer]+err[:,:8]@bound[:8];eerrors=ctx.block('epsilon_error')[outer]+err[:,8:]@bound[8:]
 coefficient_ulp=64*np.finfo(float).eps*(abs(np.array(B['state']))+np.abs(Z).T@bound)
 amplification=ctx.block('row_amplification')[outer];kcoeff=float(np.max(coefficient_ulp[rows['k']]));ecoeff=float(np.max(coefficient_ulp[rows['epsilon']]))
 errors+=amplification[:,0]*kcoeff;eerrors+=amplification[:,1]*ecoeff
 rounding=float(np.max((np.maximum(P.b,0)*errors+(2/3)*eerrors)/P.ec))
 def save(out):atomic(a.output/'progress.json',out)
 out,coeff=solve(P,K0,budget,save,rounding_allowance=rounding);out.update(created_Asia_Shanghai=now(),context_file=str(a.context),context_sha256=sha(a.context),B_sha256=B_HASH,base_sha256=COLD_HASH,K0=K0,setup_seconds=setup_seconds,rounding_allowance_normalized_zero_margin=rounding,implied_coefficient_box=box_info)
 if coeff is not None:
  coeff=np.asarray(coeff);state=np.array(B['state'],dtype='<f8').copy()
  if np.any(coeff):state[allowed]+=coeff@Z[:,allowed]
  protected=np.ones(2804,bool);protected[allowed]=False;require(bytes_equal(state[protected],np.asarray(B['state'])[protected]),'fixed B bits changed');require(bytes_equal(state[:1520],B['state'][:1520]) and bytes_equal(state[2698:],B['state'][2698:]),'B velocity/pressure changed')
  candidate=dict(version='H2hFixedBScalarCandidateV1',stage='B',state=state.tolist(),state_sha256=state_hash(state),base_sha256=COLD_HASH,B_sha256=B_HASH,scalar_increments=coeff.tolist(),context_sha256=sha(a.context),basis_sha256=ctx.doc['basis_sha256'],divergence_limit=.06030809942199109,velocity_and_pressure_B_bitwise=True,classification=out['classification'])
  atomic(a.output/'candidate.json',candidate);sc=P.scan(coeff,1e-4*K0,budget);j=sc['worst'];record=int(P.indices[j]);out['terminal']=dict(candidate_file=str(a.output/'candidate.json'),candidate_sha256=candidate['state_sha256'],minimum_covariance=sc['minimum'],normalized_covariance=sc['minimum']/K0,k_ratio=sc['k_ratio'],epsilon_ratio=sc['epsilon_ratio'],worst_record=record,physical=ctx.block('physical')[record].tolist(),gradient=ctx.block('gradient')[record].reshape(3,3).tolist(),k=float(sc['k'][j]),epsilon=float(sc['e'][j]),cold_k=float(P.kc[j]),cold_epsilon=float(P.ec[j]),strain_max_eigenvalue=float(P.s[j]))
 out['budget']=budget.report();out['peak_RSS_KiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss;atomic(a.output/'result.json',out)

def certify(a):
 c=load(a.candidate);B=load(a.parent/'standard-B/best-candidate.json');ctx=Context(a.context);pre=load(a.preflight);K0=pre['models'][1]['K0'];a.output.mkdir(exist_ok=False);ledger=a.output.parent/'native-candidate-ledger.json';entries=load(ledger) if ledger.exists() else []
 require(len(entries)<3,'H2h native candidate-check budget');entries.append(dict(candidate=c['state_sha256'],started_Asia_Shanghai=now(),fresh_order=24));atomic(ledger,entries)
 xml=case(a.output/'case',a.h1,pre,flat=False);p=dict(operation='h2h_candidate',inner_stress_model='Linear',xml_sha256=sha(xml),expected_model='StandardKEpsilon',state=c['state'],state_sha256=c['state_sha256'],base_sha256=COLD_HASH,coupled=True,B_state=B['state'],B_sha256=B_HASH)
 started=time.monotonic();d,r=run_probe(a.probe,xml,p,a.output/'native');d=load(Path(r['output']));model=load(Path(r['output']+'.model.json'));au=d['audit'];ex=d['extra_sampling'];fresh=d['fresh_sampling'];require(d['state_sha256']==c['state_sha256'],'candidate changed in audit');limit=c['divergence_limit'];delta=1e-4*K0
 gates=dict(native_hard=au['accepted_awf_hard_admissible']==[1] and au['admissible']==[1],native_positive_margin=min(au['accepted_awf_covariance'])>=delta-1e-12,optimizer_interior_margin=min(au['accepted_awf_covariance'])>=delta-1e-12,budget=au['accepted_awf_scalar_domain'][8]<=1e-10,algebraic=max(au['coupled_norms'][4:])<1e-9,GL20_margin=ex['minimum_covariance']>=delta-1e-12,GL20_scalars=ex['scalar_positive_finite'] and ex['scalar_envelopes_passed'],fresh_GL24_margin=fresh['minimum_covariance']>=delta-1e-12,fresh_GL24_scalars=fresh['scalar_positive_finite'] and fresh['scalar_envelopes_passed'],strong_divergence=max(ex['strong_divergence_L2'],fresh['strong_divergence_L2'])<=limit+1e-8,velocity=max(ex['velocity_change_L2'],fresh['velocity_change_L2'])<=.12*(1+1e-10) and max(ex['maximum_sampled_velocity_increment'],fresh['maximum_sampled_velocity_increment'])<=.5+1e-10,B_velocity_pressure_exact=c['velocity_and_pressure_B_bitwise'] and d['H2h_fixed_B_contract']['non_scalar_bits_exact'],profile_invariant=d['H2h_fixed_B_contract']['profile_invariant'],fresh_identity=fresh.get('state_sha256')==c['state_sha256'] and fresh.get('cold_base_sha256')==COLD_HASH and fresh.get('quadrature_order')==24 and fresh.get('model_specification')==model)
 receipt=dict(version='H2hNativeFreshAuditV1',passed=all(gates.values()),gates=gates,K0=K0,state_sha256=c['state_sha256'],base_sha256=COLD_HASH,B_sha256=B_HASH,probe_output_file=r['output'],probe_output_sha256=sha(Path(r['output'])),model_specification=model,audit=au,extra_sampling=ex,fresh_sampling=fresh,probe_runtime=identity(a.probe),solver_runtime=identity(a.solver),context_sha256=sha(a.context),basis_sha256=ctx.doc['basis_sha256'],geometry_identity=d['geometry_identity'],inventory_identity=d['inventory_identity'],seconds=time.monotonic()-started,created_Asia_Shanghai=now(),governing_solves=0,physical_commits=0)
 atomic(a.output/'native-receipt.json',receipt);entries[-1].update(passed=receipt['passed'],receipt=str(a.output/'native-receipt.json'));atomic(ledger,entries)
 if not receipt['passed']:return
 require_fresh(fresh,c['state_sha256'],COLD_HASH,model)
 source=Path(__file__).resolve().parents[4];names=['library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp','library/UnitTests/SolverUtils/RANSWallScalarFeasibilityContext.h','solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_h2h_scalar.py'];prov=dict(K0=K0,divergence_limit=limit,basis_sha256=ctx.doc['basis_sha256'],context_sha256=sha(a.context),native_receipt_file=str(a.output/'native-receipt.json'),source_revision='H2h exact source manifest',source_hashes={n:sha(source/n) for n in names});cert=certificate_v2(d,receipt,d['fixed_reference'],identity(a.solver),identity(a.probe),prov);cert['model_specification']=model;cert['H2h_fresh_audit']=dict(receipt=str(a.output/'native-receipt.json'),sha256=sha(a.output/'native-receipt.json'),candidate=c['state_sha256'],required_order=24);atomic(a.output/'certificate.json',cert)

def pilot(a):
 cert=load(a.certificate);binding=cert.get('H2h_fresh_audit');require(binding,'H2h runner needs bound fresh-audit prerequisite');path=Path(binding['receipt']);require(sha(path)==binding['sha256'],'H2h fresh receipt hash');r=load(path);require(r['passed'] and all(r['gates'].values()),'H2h native/fresh gates failed');require_fresh(r['fresh_sampling'],cert['state_sha256'],COLD_HASH,cert['model_specification']);require(r['solver_runtime']==identity(a.solver),'H2h solver runtime changed')
 run=StartupRun(a.solver,a.certificate,a.preflight,a.h1,a.output,builder=lambda folder,h1,pre,flat=False:case(folder,h1,pre,flat=flat));run.pair()

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('phase',choices=('prepare','search','certify','pilot'))
 for name in ('parent','probe','solver','h1','preflight','context','candidate','certificate'):p.add_argument('--'+name,type=Path)
 p.add_argument('--output',required=True,type=Path);a=p.parse_args();os.environ.update(THREADS);globals()[a.phase](a)
if __name__=='__main__':main()
