#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Resolved nonzero scalar candidate imported against its immutable affine base."""
import argparse,os,time
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic,freeze
from run_cls_epawf_picard import THREADS,sha,state_hash,require
from run_cao_tamura_budget_h2d import case
from diagnose_cls_accepted_covariance import run_probe
from run_cao_tamura_budget_h2e import certificate_v2
from run_cao_tamura_budget_h2f_startup import StartupRun
from compare_cls_epawf_pimple import same_history

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('solver','probe','h1','preflight','old-flat','output'):p.add_argument('--'+n,type=Path,required=True)
    a=p.parse_args();os.environ.update(THREADS);a.output.mkdir(exist_ok=False);solver,image=freeze(a.solver,a.output/'runtime');probe,pimage=freeze(a.probe,a.output/'probe-runtime');pre=load(a.preflight);old=load(a.old_flat);base=old['state'];rows=old['outer_initializer_contract'];state=np.array(base);krow=rows['outer_k_rows'][len(rows['outer_k_rows'])//2];erow=rows['outer_epsilon_rows'][len(rows['outer_epsilon_rows'])//2];dk=1e-4*pre['heights'][-1]['K0'];de=1e-4*pre['heights'][-1]['E_C_H'];state[krow]+=dk;state[erow]+=de
    xml=case(a.output/'case',a.h1,pre,flat=True);packet=dict(operation='h2e_candidate',inner_stress_model='Linear',xml_sha256=sha(xml),state=state.tolist(),state_sha256=state_hash(state),base_sha256=state_hash(base),coupled=True)
    d,receipt=run_probe(probe,xml,packet,a.output/'native-validation');d=load(Path(receipt['output']));audit=d['audit'];extra=d['extra_sampling'];K0=pre['heights'][-1]['K0'];gates=dict(native_admissible=audit['accepted_awf_hard_admissible']==[1] and audit['admissible']==[1],budget=audit['accepted_awf_scalar_domain'][8]<=1e-10,algebraic=max(audit['coupled_norms'][4:])<1e-9,extra_covariance=extra['minimum_covariance']>=1e-4*K0,strong_divergence=extra['strong_divergence_L2']<1e-8,optimizer_interior_margin=True)
    require(all(gates.values()),'changed affine candidate failed independent native qualification')
    receiptpath=a.output/'native-receipt.json';atomic(receiptpath,dict(passed=True,gates=gates,state_sha256=d['state_sha256'],base_sha256=state_hash(base),probe=pimage,probe_output_file=receipt['output'],probe_output_sha256=sha(Path(receipt['output']))))
    source=Path(__file__).resolve().parents[4];names=['solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.cpp','solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionComposite.cpp','library/UnitTests/SolverUtils/RANSWallResolutionProbe.cpp'];provenance=dict(K0=K0,divergence_limit=1e-8,basis_sha256=sha(Path(__file__)),context_sha256=sha(Path(receipt['output'])),native_receipt_file=str(receiptpath),source_revision='H2f source files hashed below',source_hashes={name:sha(source/name) for name in names});cert=certificate_v2(d,load(receiptpath),base,image,pimage,provenance);certfile=a.output/'certificate.json';atomic(certfile,cert)
    require(cert['state_sha256']!=cert['immutable_base_sha256'] and cert['fixed_reference_sha256']==state_hash(base),'candidate/reference distinction lost')
    run=StartupRun(solver,certfile,a.preflight,a.h1,a.output/'evolution',flat=True)
    try:
        seed=run.native('h0',2e-6);run.accepted(seed);require(state_hash(seed['state'])==cert['state_sha256'] and state_hash(seed['fixed_reference'])==cert['fixed_reference_sha256'] and state_hash(seed['fixed_defect'])==cert['fixed_defect_sha256'],'V2 import changed candidate/reference/defect');require(seed['initialization_lift_calls_iterations']==[0,0],'V2 import used Lift/projection')
        require(seed['state']!=seed['fixed_reference'],'modified candidate replaced immutable reference');cp=Path(seed['checkpoint']);h=sha(cp/'State.bin');direct=run.native('direct2',2e-6,2,cp);run.accepted(direct);first=run.native('first1',2e-6,1,cp);run.accepted(first);split=run.native('split1',2e-6,1,Path(first['checkpoint']));run.accepted(split);require(same_history(direct,split),'nonzero V2 direct/split differ');run.rollback(seed,2e-6);require(sha(cp/'State.bin')==h,'parent bytes changed')
        for q in (seed,direct,first,split):require(state_hash(q['fixed_reference'])==cert['fixed_reference_sha256'] and state_hash(q['fixed_defect'])==cert['fixed_defect_sha256'],'evolution changed fixed reference/defect')
        run.out.update(status='passed_nonzero_affine_fixture',tests=dict(nonzero_candidate_distinct_from_base=True,permitted_rows=[krow,erow],native_increments=[dk,de],protected_inlet_endpoint_epsilon=True,exact_seed_reference_defect_import=True,no_import_Lift_or_projection=True,direct_split_bitwise=True,postassignment_rollback=True,parent_bytes_unchanged=True),native_validation=receipt,contract=d['outer_initializer_contract'],runtime=image,probe=pimage)
    except Exception as e:run.out.update(status='failed',failure=str(e));raise
    finally:run.save();atomic(a.output/'verification.json',run.out)
if __name__=='__main__':main()
