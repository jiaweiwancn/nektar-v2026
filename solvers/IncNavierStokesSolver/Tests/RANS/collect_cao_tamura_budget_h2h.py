#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Compact H2h evidence; large immutable fields remain in Linux storage."""
import argparse,hashlib,json,os,subprocess,sys,ast
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic,identity
from run_cls_epawf_picard import sha,state_hash,require

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();r=a.root;old=r.parent/'cls-cao-tamura-budget-h2g-20261113';src=Path(__file__).resolve().parents[4];a.output.mkdir(exist_ok=False);date=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat();result=load(r/'scalar-search/result.json');analysis=load(r/'active-constraint-analysis.json');receipt=load(r/'terminal-audit/native-receipt.json');native=load(Path(receipt['probe_output_file']));ctx=load(r/'scalar-context/context.json');B=load(old/'standard-B/best-candidate.json');candidate=load(r/'scalar-search/candidate.json');cal=load(old/'comparison.json');previous=load(old/'evidence-final/native.json');cls_B_record=load(src/'docs/plans/cls-cao-tamura-budget-h2f-evidence/native.json')['bounded_searches']['B'];probe=r/'probe-runtime/RANSWallResolutionProbe';solver=old/'runtime-v2/IncNavierStokesSolver'
 require(not (r/'terminal-audit/certificate.json').exists() and not receipt['passed'],'collector expects the bounded negative H2h outcome')
 # Independent actual tensors, not only a formula-level optimizer check.
 tensor_checks=[]
 for name in ('extra_sampling','fresh_sampling'):
  w=receipt[name]['witness'];g=np.array(w['gradient']);k=w['k'];e=w['epsilon'];R=(2/3)*k*np.eye(3)-.09*k*k/e*(g+g.T);ev,V=np.linalg.eigh(R);err=float(np.max(abs(R-np.array(w['tensor']))));require(err<1e-12 and abs(ev[0]-w['eigenvalues'][0])<1e-12,'independent terminal tensor differs');tensor_checks.append(dict(set=name,tensor_error=err,eigenvalue=float(ev[0]),eigenvalue_error=float(abs(ev[0]-w['eigenvalues'][0])),critical_eigenvector=V[:,0].tolist(),eddy_viscosity=.09*k*k/e,production=float(-np.sum(R*g)),isotropic=(2/3)*k,physical=w['physical']))
 n=np.array(candidate['state'],dtype='<f8');b=np.array(B['state'],dtype='<f8');rows=ctx['context']['scalar_rows'];allowed=set(rows['k']+rows['epsilon']);frozen=[i for i in range(len(n)) if i not in allowed];require(n[frozen].tobytes()==b[frozen].tobytes(),'candidate changed protected B bits')
 fields={}
 for name,lo,hi in [('U',0,760),('W',760,1520),('all_pressure',2698,2804)]:fields[name]=dict(B_sha256=state_hash(b[lo:hi]),candidate_sha256=state_hash(n[lo:hi]),bitwise_equal=b[lo:hi].tobytes()==n[lo:hi].tobytes())
 atomic(a.output/'policy.json',load(src/'docs/plans/cls-cao-tamura-budget-h2h-evidence/policy.json'))
 atomic(a.output/'optimizer.json',result);atomic(a.output/'active-constraints.json',analysis)
 context={k:v for k,v in ctx['context'].items() if k!='scalar_rows'};context['permitted_row_counts']={k:len(v) for k,v in rows.items()};context['permitted_rows_sha256']=hashlib.sha256(json.dumps(rows,separators=(',',':')).encode()).hexdigest();context.update(array_directory=ctx['array_directory'],context_sha256=sha(r/'scalar-context/context.json'),conditioning=ctx['conditioning'],basis_sha256=ctx['basis_sha256'],transform_sha256=ctx['transform_sha256'],setup_wall_seconds=ctx['setup_wall_seconds'],observation_blocks=ctx['observation_blocks'],raw_context_bytes=sum((Path(ctx['array_directory'])/v['file']).stat().st_size for v in ctx['context']['files'].values()),context_usage=(r/'context-usage.txt').read_text().strip())
 atomic(a.output/'context.json',context)
 atomic(a.output/'native.json',dict(version='H2hNativeOutcomeV1',capture_Asia_Shanghai=date,classification=result['classification'],physical_hill_attempts=0,physical_hill_accepted_steps=0,governing_solves=0,certificate_issued=False,pilot_reason='failed complete covariance and fresh GL24 margin/envelope gates',receipt=receipt,exact_B_velocity_pressure=fields,complete_fixed_B_contract=native['H2h_fixed_B_contract'],reference_sha256=native['reference_sha256'],defect_sha256=native['defect_sha256'],independent_tensor_checks=tensor_checks,normalized_uncancelled_alpha0_coupled_norms=native['audit']['coupled_norms'],audit_usage=(r/'audit-usage.txt').read_text().strip(),search_usage=(r/'search-usage.txt').read_text().strip()))
 tests=load(r/('scalar-unit-tests-publish.json' if (r/'scalar-unit-tests-publish.json').exists() else 'scalar-unit-tests-final.json'));tests['native_affine_regression']=load(r/'native-regression.json');tests['H2g_model_AWF_regression']=load(r/'kernel-regression.json');tests['native_context_verification']=ctx['context']['verification'];atomic(a.output/'verification.json',tests)
 if (r/'io-cost.json').exists():atomic(a.output/'io-cost.json',load(r/'io-cost.json'))
 # Full-precision H2g manifest is referenced by bytes, never silently rounded
 # or repurposed as a new calibration.
 atomic(a.output/'comparison-manifest.json',dict(version='H2hContinuingTwoModelComparisonV1',capture_Asia_Shanghai=date,comparison_status='open: common flat spatial, developed-flow, temporal/spatial sensitivity and held-out experimental gates remain',calibration_manifest=dict(path=str(old/'comparison.json'),sha256=sha(old/'comparison.json')),unchanged_full_precision_calibration=cal,CLS_prior=previous['models']['cls'],CLS_retained_B={k:cls_B_record[k] for k in ('candidate_sha256','native_covariance','extra','gates','mass_distance','max_component_response','K_envelope','epsilon_envelope')},standard_prior=previous['models']['standard'],standard_B_sha256=B['state_sha256'],standard_cold_sha256=candidate['base_sha256'],H2h_terminal_candidate_sha256=candidate['state_sha256'],H2h_outer_comparison=dict(B=analysis['B'],terminal=analysis['terminal']),H2h_native_inner_outer=native['audit']['accepted_awf_covariance'],H2h_GL20=receipt['extra_sampling'],H2h_fresh_GL24=receipt['fresh_sampling'],H2h_scalar_envelopes='relative to immutable standard cold; no resetting at B',H2h_exact_B_velocity_pressure=fields,physical_hill_attempts=0,physical_hill_steps=0,physical_step_cost='no H2h hill step qualified; no hill-step measurement',historical_flat_step_cost='H2g: about 4 seconds per continued native step; not a matched optimizer benchmark',flat_U_k_spatial_gate='open: 1e-3 target not met',normal_stress_qualification='standard neutral normal stresses are isotropic; no experimental RMS calibration',hill_observations='held out; no recalibration'))
 source_names=['library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp','library/UnitTests/SolverUtils/RANSWallScalarFeasibilityContext.h','solvers/IncNavierStokesSolver/EquationSystems/RANSWallModelSpec.h','solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp']+['solvers/IncNavierStokesSolver/Tests/RANS/'+n for n in ('cao_tamura_h2h_scalar.py','run_cao_tamura_budget_h2h.py','test_cao_tamura_h2h_scalar.py','analyze_cao_tamura_budget_h2h.py','collect_cao_tamura_budget_h2h.py')]
 def node(path,name):
  tree=ast.parse(path.read_text());return ast.dump(next(n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef)) and n.name==name),include_attributes=False)
 oldmodule=r/'LP-used-source/cao_tamura_h2h_scalar.py';newmodule=src/'solvers/IncNavierStokesSolver/Tests/RANS/cao_tamura_h2h_scalar.py';require(all(node(oldmodule,k)==node(newmodule,k) for k in ('Problem','solve','dual_check','Budget')),'optimization arithmetic changed after LP')
 try:cpu=subprocess.check_output(['c++','--version'],text=True).splitlines()[0]
 except Exception:cpu='unavailable'
 import scipy,scipy.optimize._linprog_highs as highs,inspect
 highs_source=inspect.getsource(highs);atomic(a.output/'source-revisions.json',dict(version='H2hSourceRuntimeV1',capture_Asia_Shanghai=date,reviewed_parent='f30c309d6944dc7a117f14d6d7e2b929bbcf0d54',implementation_parent='07b680966481d331d1d1d727d118e425829b3f51',source_hashes={v:sha(src/v) for v in source_names},source_verification=load(r/'source-verification.json'),LP_used_sources={p.name:sha(p) for p in (r/'LP-used-source').iterdir()},optimization_arithmetic_unchanged_after_context_validation_hardening=True,solver=identity(solver),probe=identity(probe),compiler=cpu,flags='Release; -fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS',SciPy=scipy.__version__,numpy=np.__version__,HiGHS_wrapper_sha256=hashlib.sha256(highs_source.encode()).hexdigest(),HiGHS_serial_dual_strategy='simplex_strategy' in highs_source and 'HIGHS_SIMPLEX_STRATEGY_DUAL' in highs_source,resources=dict(build_workers=4,native_ranks=1,computational_threads=1,OMP_THREAD_LIMIT=1,total_limit=72),storage_UUID='5AB8A669B8A642FD',storage_mount='/media/jay/新加卷1',parent_verification=load(r/'parent-verification.json'),native_candidate_checks=load(r/'native-candidate-ledger.json')))
 # Index large outputs without copying the isolated source/build tree into Git.
 files=[]
 for p in sorted(r.rglob('*')):
  if not p.is_file() or p.is_symlink():continue
  relative=p.relative_to(r)
  if relative.parts[0] in ('source','build','LP-used-source') or relative.parts[0].startswith('evidence') or relative.name.startswith('h2h-') or relative.name.endswith('.tgz') or relative.name=='source.tar':continue
  files.append(dict(path=str(relative),bytes=p.stat().st_size,sha256=sha(p)))
 atomic(a.output/'artifact-index.json',dict(version='H2hArtifactIndexV1',root=str(r),capture_Asia_Shanghai=date,files=files,source_tree='verified by source-revisions/source-verification; build isolated; historical runtimes unchanged',reproduction='see H2h report commands and immutable context/basis/LP-used-source hashes'))
 print('H2h compact evidence complete',flush=True)
if __name__=='__main__':main()
