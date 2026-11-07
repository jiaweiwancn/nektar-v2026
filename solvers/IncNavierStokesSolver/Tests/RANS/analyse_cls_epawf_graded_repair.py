#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent mass/QP/frozen-transport accounting on exact attempt-65 replays."""
import argparse
import datetime
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve,cholesky
from run_cls_epawf_picard import load,write_json,require,sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    p.add_argument('--old-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    archive=a.old_root/'terminal-U05-graded';raw=load(archive/'attempt.json');rows=np.array(raw['transport_rows'],dtype=int)
    n=len(rows);A=np.fromfile(archive/'transport_A.bin',dtype='<f8').reshape((n,n))
    require(sha(archive/'transport_A.bin')==raw['matrix_sha256'],'archived A_k changed')
    reference_mass=load(a.old_root/'discrete-budget/U05-graded/mass-and-bands.json')
    require(np.array_equal(rows,np.array(reference_mass['free_k_rows'],dtype=int)),'physical mass free rows differ')
    M=np.array(reference_mass['free_k_mass']);factor=cho_factor(M,lower=True);L=cholesky(M,lower=True)
    massnorm=lambda v:float(np.sqrt(max(0.,v@M@v)))
    dualnorm=lambda v:float(np.sqrt(max(0.,v@cho_solve(factor,v))))
    lo=int(reference_mass['offsets'][2]);ell=np.array(reference_mass['k_integral_load'])[rows-lo]
    center=np.array(raw['raw']);pre=np.array(raw['pre_k']);r=np.array(raw['transport_residual']);step=np.array(raw['raw_increment'])
    rawResidual=r+A@step;require(max(abs(rawResidual))<1e-12,'raw transport system no longer accurate')
    out=dict(captured_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),ranks=1,threads=1,
        governing_solves=0,metric='actual consistent PicardKMass Mff, not frozen A_k',replays={},cross_pools={})
    solutions={};pools={}
    for name in ('replay-original','replay-pool-v1'):
        folder=a.root/name;q=load(folder/'result.json');b=q['k_bound'];pool=load(folder/'pool.json')
        actualM=np.fromfile(folder/'pool.json.Mff.bin',dtype='<f8').reshape((n,n))
        require(np.array_equal(M,actualM),'replay and independent native mass differ')
        delta=(np.array(q['corrected'])-center)[rows];y=np.array(pool['whitened_primal']);active=np.array(pool['solved_active'],dtype=int)
        C=np.fromfile(folder/'pool.json.Cwhite.bin',dtype='<f8').reshape((-1,n));rhs=np.array([c['rhs'] for c in pool['pool']])
        lam=np.array(pool['solved_lambda']);action=np.array(np.sum(C[active].astype(np.longdouble)*lam.astype(np.longdouble)[:,None],axis=0),dtype=float)
        objective=.5*massnorm(delta)**2;dual=float(lam@rhs[active]-.5*(action@action))
        stationarity=float(np.linalg.norm(y-action));violation=float(max(0.,max(rhs-C@y)))
        added=A@delta;storage=raw['alpha']*(M@delta);corrected=r+A@(step+delta)
        require(max(abs(corrected-added-rawResidual))<1e-12,'corrected frozen residual identity failed')
        require(abs(float(ell@delta)-b['integral_total'])<1e-13,'physical correction integral differs from native accounting')
        require(abs(objective-b['primal_objective'])<1e-16,'original-raw physical mass objective differs')
        out['replays'][name]=dict(success=b['success'],failure=b['failure'],bound=b,
            result_sha256=sha(folder/'result.json'),pool_sha256=sha(folder/'pool.json'),mass_sha256=sha(folder/'pool.json.Mff.bin'),
            physical_objective=objective,finite_pool_dual_objective=dual,physical_primal_dual_gap=objective-dual,
            independent_stationarity=stationarity,independent_pool_normalized_violation=violation,
            independent_dual_violation=float(max(0.,max(-lam))),
            independent_complementarity=float(max(abs(lam*(C[active]@y-rhs[active])))),
            raw_linear_residual_max=float(max(abs(rawResidual))),corrected_frozen_residual_max=float(max(abs(corrected))),
            added_frozen_residual_max=float(max(abs(added))),storage_reaction_max=float(max(abs(storage))),
            added_frozen_residual_mass_dual=dualnorm(added),storage_reaction_mass_dual=dualnorm(storage),
            raw_step_mass_norm=massnorm(step),correction_mass_norm=massnorm(delta),accepted_step_mass_norm=massnorm(step+delta),
            correction_integral=float(ell@delta),increment_identity_defect=float(max(abs((step+delta)-step-delta))),
            hard_audit=q.get('complete_hard_admissible'),kkt_describes_full_pool=bool(b.get('finite_pool_complete',False)))
        solutions[name]=delta;pools[name]=(C,rhs,y)
    for target,(C,rhs,y) in pools.items():
        out['cross_pools'][target]={name:float(max(0.,max(rhs-C@(L.T@delta)))) for name,delta in solutions.items()}
    out['between_replay_mass_difference']=massnorm(solutions['replay-pool-v1']-solutions['replay-original'])
    out['interpretation']='Each replay retains every cut it discovers. Resolving the finite QP changes subsequent moving-cut discovery; the failed old pool is not native checkpoint history.'
    write_json(a.output,out);print(a.output)


if __name__=='__main__':main()
