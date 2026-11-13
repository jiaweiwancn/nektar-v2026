#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Complete native/GL20 H2f seed gates, then exact-runtime V2 issuance if passing."""
import argparse,os
from pathlib import Path
from cao_tamura_h2d_state import load
from run_cao_tamura_budget_h2e import certify,certificate_v2
from continue_cls_epawf_pimple import atomic,identity
from run_cls_epawf_picard import sha,require,THREADS

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('context','candidate','preflight','probe','output'):p.add_argument('--'+n,type=Path,required=True)
    p.add_argument('--solver',type=Path);a=p.parse_args();os.environ.update(THREADS);certify(a);path=a.output/'certification.json';r=load(path);K0=load(a.preflight)['heights'][-1]['K0'];r['K0']=K0
    r['gates']['native_inner_search_margin']=r['audit']['accepted_awf_covariance'][0]>=1e-4*K0-1e-12;r['gates']['native_positive_startup_margin']=min(r['audit']['accepted_awf_covariance'])>=1e-4*K0-1e-12;r['gates']['extra_positive_startup_margin']=r['extra_sampling']['minimum_covariance']>=1e-4*K0-1e-12
    r['passed']=all(r['gates'].values());atomic(path,r)
    if a.solver and r['passed']:
        solver=a.solver.resolve();os.environ['LD_LIBRARY_PATH']=str(solver.parent/'lib');image=identity(solver);native=load(Path(r['probe_output_file']));base=load(a.context);from cao_tamura_initializer_context import NativeContext
        ctx=NativeContext(a.context);source=Path(__file__).resolve().parents[4];names=['solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionComposite.cpp','library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp'];prov=dict(K0=K0,divergence_limit=load(a.candidate)['divergence_limit'],basis_sha256=ctx.meta['files']['basis']['sha256'],context_sha256=sha(a.context),native_receipt_file=str(path),source_revision='H2f exact source hashes',source_hashes={n:sha(source/n) for n in names});cert=certificate_v2(native,r,ctx.array('base')[0].tolist(),image,r['runtime'],prov);atomic(a.output/'startup-certificate.json',cert)
if __name__=='__main__':main()
