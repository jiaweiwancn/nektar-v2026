#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Matched-T=.3 pressure-space forensic control; no acceptance inferred."""
import argparse,os,time
from pathlib import Path
import numpy as np
from cao_tamura_h2d_state import load
from continue_cls_epawf_pimple import atomic
from run_cls_epawf_picard import require,sha,state_hash,THREADS
from analyze_wp2_recovery_mesh_control import same_physical_problem
from wp2_all_state import Rows,US
from wp2_common_quadrature import resample
from wp2_independent_reference import evaluate
from audit_wp2_dense_profiles import gates
from run_cao_tamura_budget_wp2 import charge


class PressureMassProjection:
    """Observer-only physical pressure mass; production dense backend unchanged."""
    def __init__(self,rows):
        from scipy.sparse import coo_matrix
        from scipy.sparse.linalg import splu
        ids=np.asarray(rows.a['p_ids']);phi=np.asarray(rows.a['p_phi']);width=ids.shape[1]
        active=np.ones_like(ids,dtype=bool) if 'p_counts' not in rows.a else np.arange(width)[None,:]<rows.a['p_counts']
        point=np.broadcast_to(np.arange(len(ids))[:,None],ids.shape)
        self.offset=rows.offset[4];self.count=rows.offset[6]-self.offset
        self.S=coo_matrix((phi[active],(point[active],ids[active]-self.offset)),shape=(len(ids),self.count)).tocsr()
        self.w=rows.weight;self.M=(self.S.T@self.S.multiply(self.w[:,None])).tocsc();self.factor=splu(self.M)
        constant=np.zeros(self.count);known=np.zeros(self.count,dtype=bool)
        for j,c in zip((ids[active]-self.offset).flat,np.asarray(rows.a['p_constant'])[active].flat):
            if known[j]:require(constant[j]==c,'inconsistent signed pressure CG sharing')
            constant[j]=c;known[j]=True
        require(np.all(known) and max(abs(self.S@constant-1))<2e-12,'pressure mass constant reproduction')
        self.constant_volume_error=abs(float(constant@(self.M@constant))-sum(self.w))
        require(self.constant_volume_error<2e-12,'pressure mass physical volume differs')
    def evaluate(self,state,tuple_data):
        pressure=np.asarray(state)[self.offset:self.offset+self.count]
        action_error=float(max(abs(self.S@pressure-tuple_data['values'][:,4])))
        require(action_error<2e-11,'signed pressure mass action differs from native-row evaluation')
        div=tuple_data['G'][:,0,0]+tuple_data['G'][:,2,2]
        rhs=self.S.T@(self.w*div);c=self.factor.solve(rhs);projected=self.S@c;orth=div-projected
        residual=float(np.linalg.norm(self.M@c-rhs)/max(np.linalg.norm(rhs),1e-300))
        norms=[float(np.dot(self.w,v*v)) for v in [div,projected,orth]]
        cross=float(np.dot(self.w*projected,orth));require(residual<1e-9 and abs(norms[0]-norms[1]-norms[2])<1e-10*max(norms[0],1e-20),'pressure divergence projection identity failed')
        return dict(strong_L2_squared=norms[0],actual_C0_pressure_projection_L2_squared=norms[1],orthogonal_L2_squared=norms[2],
            projected_fraction=norms[1]/norms[0] if norms[0]>1e-25 else None,orthogonality_cross=cross,mass_solve_relative_residual=residual,
            physical_mass_nnz=int(self.M.nnz),pressure_DOFs=self.count,constant_volume_error=self.constant_volume_error,pressure_action_error=action_error,
            scope='Actual conforming regional pressure test space, all signed CG rows; exact flat polynomial integration GL8; observer mass solve, no governing solve or production backend change. Tiny projected fractions are numerical-floor diagnostics.')


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--project-divergence',action='store_true');a=p.parse_args()
    os.environ.update(THREADS);started=time.monotonic();result={};bindings={}
    rows=[Rows(a.root/'geometry-rows-selected/result.json'),Rows(a.root/'pressure-row-verification/mixed_p2/result.json')]
    projections=[PressureMassProjection(r) for r in rows] if a.project_divergence else None
    for model in ['CLS','standard']:
        old=load(a.root/('preservation-recovery-'+model)/'progress.json')
        obs=next(o for o in old['observations'] if abs(o['time']-.3)<2e-14)
        reports=[Path(obs['report']),a.root/('pressure-recovery-'+model)/'step30/case.evolution.json']
        data=[load(f) for f in reports];same_physical_problem(*data)
        current=load(a.root/('pressure-recovery-'+model)/'progress.json')
        require(old['solver_sha256']==current['solver_sha256'],'pressure comparison changes frozen native runtime')
        require(data[0]['dt']==data[1]['dt']==.01 and abs(data[0]['time']-data[1]['time'])<2e-14,'pressure comparison times differ')
        offset=int(data[1]['offsets'][4]);require(offset==int(data[0]['offsets'][4]),'changed velocity/scalar spaces')
        require(state_hash(data[0]['accepted_states'][0][:offset])==state_hash(data[1]['accepted_states'][0][:offset]),
                'initial nonpressure coefficients differ')
        metrics=[]
        for row,d,f in zip(rows,data,reports):
            v=row.tuple(d['state'],d);actual=resample(dict(d,wp2_common_samples=dict(d['wp2_common_samples'],values=v['values']),physical_samples=dict(wall_K=v['endpoint_K'])),36)
            reference=evaluate(d,actual['physical']);check=gates(actual,reference)
            div=actual['G'][:,0,0]+actual['G'][:,2,2];w=actual['weights'];D=float(np.sqrt(np.dot(w,div*div)/sum(w))/(US/.0065))
            R=actual['R'];ev=np.linalg.eigvalsh(R)[:,0]
            metrics.append(dict(slots=row.size,actual_time=d['time'],strong_divergence_RMS_over_ustar_H=D,
                covariance_minimum=float(min(ev)),profile_gates=check,
                interface_normal_full_kflux_native=d['state_audits'][-1]['interface'][2:6],
                native_free_block_norms=d['wp2_per_state_monitoring'][-1]['physical_checks']['unanchored_free_block_norms']))
            if projections:metrics[-1]['actual_pressure_space_divergence_projection']=projections[len(metrics)-1].evaluate(d['state'],v)
            bindings[str(f)]=sha(f)
        result[model]=dict(P1=metrics[0],inner_P2=metrics[1],strong_divergence_ratio_P2_over_P1=metrics[1]['strong_divergence_RMS_over_ustar_H']/metrics[0]['strong_divergence_RMS_over_ustar_H'],
            initial_nonpressure_coefficients_bitwise_equal=True,physics_and_runtime_matched=True)
    charge(a.root,'matched pressure-space recovery forensic at T=.3',readback_comparisons=1)
    atomic(a.output,dict(version='WP2MatchedPressureRecoveryV2' if a.project_divergence else 'WP2MatchedPressureRecoveryV1',models=result,bindings=bindings,source_sha256=sha(Path(__file__)),
        row_catalogues=[dict(path=str(r.path),sha256=sha(r.path)) for r in rows],seconds=time.monotonic()-started,governing_solves=0,
        scope='Same model, runtime, dt, geometry, velocity/scalar spaces and initial nonpressure bits; only inner pressure representation differs. Actual T=.3, not a recovery or continuum-stability certificate.'))


if __name__=='__main__':main()
