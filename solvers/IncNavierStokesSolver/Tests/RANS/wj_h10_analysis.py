#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""H10 state invariants and saved-system matvecs; never solve a flow system."""
import numpy as np
from wj_h4_points import vector_hash
from wj_h9_local import same,external_ids,phase_bindings,MomentAudit

MATRIX=(('coarse',.01),('fine',.005),('quarter',.0025))
QUALITY=.002068763923454

def remaining_matrix(book,resume=False):
    if not resume:return MATRIX
    assert not book.get('h10_halt')and not book.get('coverage_halt')
    assert (book['cp0_reservations'],book['governing_attempts'],book['strict_readbacks'])==(1,1,1)
    return MATRIX[1:]


def preservation(r):
    d=r['initial_velocity_diagnostic'];a=np.asarray(d['pre_state']);b=np.asarray(d['post_state']);off=r['offsets']
    assert r['initial_velocity_mode']=='FlatInletTraceV1'and r['initial_velocity_calls']==1
    assert vector_hash(a)==d['pre_sha256']and vector_hash(b)==d['post_sha256']==vector_hash(r['state'])
    assert same(a[off[2]:],b[off[2]:])and same(a,r['fixed_reference'])
    assert vector_hash(r['fixed_reference'])==d['fixed_reference_sha256']and vector_hash(r['fixed_defect'])==d['fixed_defect_sha256']
    essential=sorted(i for i in external_ids(r)if i<off[2]);assert same(a[essential],b[essential])
    return dict(status='PASS',pre_sha256=vector_hash(a),post_sha256=vector_hash(b),essential_velocity_count=len(essential),essential_bits=True,nonvelocity_bits=True,reference_defect_bits=True)


def gates(summary,moments,*,startup):
    if summary['status']!='PASS':return summary['status']
    if summary['min_in_plane']<QUALITY:return 'QUALITY_LOST'
    if not all(np.isfinite(x)for x in (summary['min_in_plane'],summary['divergence_scaled_max'],summary['min_covariance'],summary['min_transverse'])):return 'UNRESOLVED'
    if set(moments)!={'4','5'} or any(len(moments[k]['moments'])!=n or not np.all(np.isfinite(moments[k]['moments']))for k,n in [('4',34),('5',70)]):return 'UNRESOLVED'
    if max(abs(x)for v in moments.values()for x in v['moments'])>2e-15:return 'OBSTRUCTED'
    if startup and summary['divergence_scaled_max']>1e-12:return 'OBSTRUCTED'
    return 'PASS'


def matrix(path,descriptor):
    from study_wj_h9 import check
    check(dict(path=str(path),sha256=descriptor['sha256'],bytes=8*int(np.prod(descriptor['shape']))))
    return np.fromfile(path,dtype='<f8').reshape(*map(int,descriptor['shape']))


def residual_check(actual,expected,scale):
    actual=np.asarray(actual);expected=np.asarray(expected);scale=np.asarray(scale)
    error=abs(actual-expected);limit=64*np.finfo(float).eps*scale+2e-14
    assert np.all(np.isfinite(error))and np.all(error<=limit)
    return dict(maximum_error=float(np.max(error)),maximum_limit=float(np.max(limit)),status='PASS')


def arithmetic(r,snapshot,momentum,pressure,mass,BG,Schur,Au,mom,step=1):
    # Defaults preserve H10. Restarted producers require an explicit absolute step.
    first=r['steps']+1-len(r['accepted_states']);local_attempt=step-first-1
    assert 0<=local_attempt<len(r['pimple_attempts']),'selected step has no completed attempt in this producer'
    states,binding=phase_bindings(r,snapshot,step);assert binding['P_C_bitwise']
    vu=np.array(momentum['rows'],int);pp=np.array(pressure['rows'],int);predictor=np.asarray(momentum['solution']);dpi=np.asarray(pressure['solution']);rhs=np.asarray(pressure['rhs'])
    assert np.all(vu<r['offsets'][2])and np.all(pp>=r['offsets'][4])and len(pp)==104
    correction=dpi@BG;du=states['F'][vu]-states['C'][vu];actual_dpi=states['F'][pp]-states['C'][pp]
    checks=dict(velocity_update=residual_check(du,predictor-correction,abs(predictor)+abs(correction)),pressure_update=residual_check(actual_dpi,dpi,abs(dpi)+abs(states['C'][pp])),
        predictor_system=residual_check(Au@predictor-np.asarray(momentum['rhs']),momentum['linear_residual'],abs(Au)@abs(predictor)+np.abs(momentum['rhs'])),
        pressure_system=residual_check(Schur@dpi-rhs,pressure['linear_residual'],abs(Schur)@abs(dpi)+abs(rhs)))
    def direction(value):
        x=np.zeros(len(states['C']));x[vu]=value;return x
    def ordered(packet):
        v={r['offsets'][int(b)]+i:x for b,p in packet.items()for i,x in enumerate(p['moments'])}
        return np.array([v[int(i)]for i in pp])
    rc=mom.moments(states['C']);rf=mom.moments(states['F']);dp=ordered(mom.moments(direction(predictor),False));dc=ordered(mom.moments(direction(correction),False));dd=ordered(mom.moments(direction(du),False));c=ordered(rc);f=ordered(rf)
    checks.update(pressure_rhs=residual_check(rhs,dp+c,abs(dp)+abs(c)),continuity_linearity=residual_check(f-c,dd,abs(f)+abs(c)+abs(dd)),
        full_Schur_action=residual_check(Schur@dpi,dc,abs(Schur)@abs(dpi)+abs(dc)))
    assert max(abs(c))<=2e-15 and max(abs(f))<=2e-15
    assert r['initialization_lift_calls_iterations']==[0,0]and r['initial_velocity_calls']==0
    a=r['pimple_attempts'][local_attempt]
    for name in ('momentum_predictor_solves','pressure_system_solves','k_transport_solves','epsilon_transport_solves','bound_calls'):assert a[name]==1
    assert a['bound_activations']==0
    return dict(status='PASS',binding=binding,checks=checks,velocity_rows=vu.tolist(),pressure_rows=pp.tolist(),predictor=predictor.tolist(),BGp_delta_pi=correction.tolist(),pressure_correction_velocity=(-correction).tolist(),complete_Gp_delta_pi=((mass@correction)/r['dt']).tolist(),correction_sign='delta_u = predictor + pressure_correction_velocity; BGp_delta_pi is subtracted',total_velocity_increment=du.tolist(),
        pressure_increment=dpi.tolist(),pressure_rhs=rhs.tolist(),D_predictor=dp.tolist(),r_p_C=c.tolist(),r_p_F=f.tolist(),D_delta_u=dd.tolist(),
        predictor_per_dt=(predictor/r['dt']).tolist(),consistent_mass_correction_action_max=float(np.max(abs(mass@correction))),consistent_mass_norms={name:float(np.sqrt(max(0.,v@(mass@v))))for name,v in [('predictor',predictor),('pressure_correction',correction),('total',du),('predictor_per_dt',predictor/r['dt'])]},
        scope='native saved vectors/maps and matvecs, complete BGp; Cpp=0 from source; no new governing factorization/solve')
