#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent reference/actions and steady-state checks for the inner-k fixture."""

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import physical_oracle as po
from inner_k_columns import cutoff, historical
from inner_k_mechanism import TOL, read
import numpy as np


POLICY = 'AuxiliaryHistoricalInnerKFixtureV1'
POLICY_ID = hashlib.sha256(f'WJH27|{POLICY}|retained-v2-boundaries|{po.TABLE_SHA}'.encode()).hexdigest()


class FixtureReference(po.Reference):
    def __init__(self, table):
        super().__init__(table)
        self.identity = POLICY_ID
        self.inner_k_transport = POLICY

    def cutoff(self, K):
        K = np.asarray(K, dtype=float)
        po.require(np.all(np.isfinite(K)&(K >= self.lo)&(K <= self.hi)), 'historical fixture K left its interval')
        a = self.kappa*np.sqrt(K/self.qe)
        d = self.nu+self.aux(K,0.)
        B = d+a*self.H
        J = np.log(B/d)+self.nu*(1/B-1/d)
        value = self.H*np.exp(1-J)
        po.require(np.all((value > 0)&(value < .4*self.H)), 'historical fixture cutoff left the lower layer')
        return value

    def k_diffusivity(self, k, epsilon, auxiliary, outer):
        return super().k_diffusivity(k,epsilon,auxiliary,outer) if outer else self.nu+auxiliary


def native_json(path):
    return json.loads(Path(path).read_text(), parse_int=lambda s: -0.0 if s == '-0' else int(s))


def state(path):
    doc = native_json(path)
    data = doc['state']
    po.require(doc['schema']=='WJH27SteadyCandidateV1' and doc['normal_main'] is False and
               doc['diagnostic_only'] and doc['status']=='CONVERGED' and not doc['qualified_reference'] and
               not doc['TASK3_READY'] and doc['physical_time_advanced']==data['time']==data['step']==0 and
               doc['diagnostic_source']==data['source_identity'], 'invalid diagnostic steady packet')
    return data, doc


def atomic_out(path, value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,allow_nan=False)
        stream.write('\n')


def coefficients(report, ref):
    symbols = 'K r c c_K inverse aux aux_K aux_r E E_K E_r D_inner D_outer'.split()
    errors = {}
    for row in report['reference_actions']:
        po.require(set(row)==set(symbols), 'incomplete native coefficient inventory')
        K,r = row['K'],row['r']
        c, dc = cutoff(ref,K)
        h = historical(ref,K,r)
        A = K**1.5/ref.cl
        E = A/max(c,r)
        eK, er = (E*(1.5/K-dc/c),0.) if r < c else (1.5*E/K,-E/r)
        Do = ref.nu+ref.cd*K*K/(ref.sk*E)
        expected = dict(c=c,c_K=dc,inverse=K,aux=h['D']-ref.nu,aux_K=h['D_K'],aux_r=h['a'],
                        E=E,E_K=eK,E_r=er,D_inner=h['D'],D_outer=Do)
        scale = dict(c=ref.H,c_K=ref.H/ref.Kin,inverse=ref.Kin,aux=h['D'],aux_K=h['D']/ref.Kin,aux_r=h['D']/ref.H,
                     E=E,E_K=E/ref.Kin,E_r=E/ref.H,D_inner=h['D'],D_outer=Do)
        for key,value in expected.items():
            errors[key] = max(errors.get(key,0.),abs(row[key]-value)/scale[key])
    po.require(len(report['reference_actions'])==25 and max(errors.values())<=TOL, 'reference action gate failed')
    return errors


def physical_scales(ref):
    return dict(k=25*ref.A, epsilon=25*(ref.nu*ref.A/ref.H**2 + ref.cd*ref.Kin**2/(ref.se*ref.H)))


def compare_actions(numerical, expected, ref):
    error = {}
    units = physical_scales(ref)
    for name,value in expected.items():
        actual = np.asarray(numerical[name])
        po.require(actual.shape==value.shape and np.all(np.isfinite(actual)), 'incomplete action '+name)
        error[name] = float(np.max(abs(actual-value)))/units[name.split('_')[0]]
    return error


def stale_diffusivity_action(oracle):
    delta = np.zeros_like(oracle.s['k'])
    r0 = oracle.ref
    for cell in oracle.elements:
        if cell['epsilon']:
            continue
        x,r,w = oracle.rule(cell)
        field = oracle.sample(cell,x,r)
        D = r0.nu+r0.cd*field['k']**2/(r0.sk*field['E'])
        _,d1,d2 = po.basis(cell,'k',x,r)
        weight = w*(D-field['Dk'])
        delta[cell['k']] += d1.T@(weight*field['gradk'][:,0])+d2.T@(weight*field['gradk'][:,1])
    return float(np.max(abs(delta)))/(25*r0.A)


def audited_actions(record, reference, order, deadline):
    data, native = record['state'], record['actions']
    oracle = po.PhysicalOracle(data,reference,order,deadline)
    core = oracle.assemble()
    boundary = oracle.boundary_checks()
    expected = {f'k_{label}':oracle.parts[key]['k'] for label,key in
                [('diffusion','diffusion'),('source','source'),('transport','advection'),('natural','natural')]}
    expected.update(k_current=oracle.total['k'], epsilon_current=oracle.total['epsilon'])
    errors = compare_actions(native,expected,reference)
    info = dict(action_errors=errors, accounting=core['accounting_identity_max_error'],
                boundary=boundary, native_dependent_trace=native['dependent_trace_error'],
                independent=current_subset(core),
                native_minima={k:data[k] for k in ('minimum_k','minimum_epsilon','minimum_covariance_over_k','minimum_principal_over_mu')})
    po.require(max(errors.values())<=TOL, f'signed native rows exceed 1e-10: {errors}')
    verify_guards(info['native_minima'],core,boundary)
    po.require(native['dependent_trace_error']<=TOL, 'native current trace differs')
    return oracle, expected, info


def current_subset(core):
    return {k:core[k] for k in ('residuals','steady_residual_status','accounting_identity_status',
                               'accounting_identity_max_error','sampled_minima','maximum_independent_cubic_relative_residual')}


def verify_guards(native_minima, core, boundary):
    po.require(native_minima['minimum_k']>0 and native_minima['minimum_epsilon']>0 and
               native_minima['minimum_covariance_over_k']>1e-10 and native_minima['minimum_principal_over_mu']>1e-8,
               'native admitted positivity, covariance or principal guard failed')
    po.require(core['sampled_minima']['covariance_over_k']>1e-10 and core['maximum_independent_cubic_relative_residual']<=TOL,
               'independent stress or WJ root gate failed')
    po.require(core['accounting_identity_status']=='MET' and boundary['status']=='MET',
               'independent conservation or essential-lift check failed')


def common_flux(report, oracle):
    ref = oracle.ref
    cells = {el['geometry']:el for el in oracle.elements}
    fluxes = {}
    po.require(len(report['H_flux'])==2 and {r['side'] for r in report['H_flux']}=={'lower','upper'},
               'both H physical owners are required')
    for item in report['H_flux']:
        item = dict(item)
        x, r = np.array([item['x']]), np.array([ref.H])
        owner = cells[item['owner']]
        a = oracle.sample(owner,x,r)
        expect = dict(K=a['K'][0],E=a['E'][0],D=a['Dk'][0],gradient=a['gradk'][0,1],Q=(a['Dk']*a['gradk'][:,1])[0])
        scales = dict(K=ref.Kin,E=ref.A/ref.H,D=expect['D'],gradient=ref.Kin/ref.H,Q=ref.A)
        errors = {key:float(abs(item[key]-expect[key])/scales[key]) for key in expect}
        target = .02*(1-abs(item['x'])/2.5)*ref.A
        item.update(action_errors=errors,expected_Q=target,target_error_over_A=abs(item['Q']-target)/ref.A)
        po.require(max([*errors.values(),item['target_error_over_A']])<=TOL, f'H owner flux action differs: {item}')
        fluxes[item['side']] = item
    a,b = fluxes['lower'],fluxes['upper']
    po.require(abs(a['Q']-b['Q'])/ref.A<=TOL and abs(a['Q']/ref.A)>.01 and
               abs(a['D']/b['D']-1)>.005 and abs(a['gradient']-b['gradient'])>1e-6,
               'required nonzero common H flux with unequal D and gradients failed')
    return fluxes


def moving_panels(report, oracle, constant_oracle):
    ref = oracle.ref
    owners = {el['geometry']:el for el in oracle.elements}
    witness = []
    for row in report['moving']['actions']['cutoff_cells']:
        e = owners[row['owner']]
        if not any(abs(value-.30)<1e-7 for value in row['x_cuts']):
            continue
        points = oracle.cutoff_x_splits(e)
        po.require(any(abs(value-.30)<1e-7 for value in points), 'independent same-owner crossing absent')
        for point in row['x_cuts']:
            if abs(point-.30)<1e-7:
                po.require(min(abs(other-point) for other in points)<=TOL*25, 'moving cutoff x root differs')
        po.require(row['points']>12*16,'cutoff witness failed to split native rule')
        witness.append(dict(native=row,independent=points,cell=e['bounds']))
    po.require(len(witness)==2,'moving cutoff must cross shared radial face of two original owners')
    e = next(e for e in oracle.elements if e['bounds'][0]<=.45<=e['bounds'][1] and e['epsilon'])
    x = np.array([.45])
    K0,K1 = [float(o.trace(e,'k',x)[0]) for o in (constant_oracle,oracle)]
    gamma_change = (K1**1.5-K0**1.5)/(ref.cl*ref.A)
    po.require(abs(K1-K0)/ref.Kin>1e-3 and abs(gamma_change)>1e-3,'moving fixture did not change K and Gamma')
    return dict(two_owner_crossing=witness,K=[K0,K1],cutoff=[float(ref.cutoff(K)) for K in (K0,K1)],
                Gamma_change_over_reference=gamma_change)


def f2(args, table, deadline):
    ref = FixtureReference(table)
    new, packet = state(args.parity)
    old, _ = po.load_native(args.v2)
    expected_keys = ('state_sha256','previous_sha256','advection_sha256','previous_advection_sha256')
    same = {k:new[k]==old[k] for k in expected_keys}
    po.require(all(same.values()) and new['reference_identity']==new['boundary_identity']==po.TABLE_SHA and
               new['inner_k_transport']=='ScalarCdV2' and max(new['residuals'])<=new['nonlinear_tolerance'] and
               new['elements']==old['elements'], 'default v2 parity is not byte-identical')
    old_audit = read(args.v2_audit)
    po.require(old_audit['native']['sha256']==po.identity(args.v2)['sha256'] and old_audit['state_sha256']==old['state_sha256'] and
               old_audit['core_steady_equations_status']=='MET','retained independent v2 assessment mismatches parity')
    native = native_json(args.native)
    po.require(native['schema']=='WJInnerKFixtureNativeActionsV1' and not native['TASK3_READY'] and
               native['cross_policy_histories_rejected'] and native['diagnostic_source']==packet['diagnostic_source'],
               'native changed-operator receipt differs')
    try:
        po.PhysicalOracle(native['moving']['state'],po.Reference(table),20,deadline)
    except ValueError as error:
        po.require('identity mismatch' in str(error),'wrong-policy oracle failed for a different reason')
    else:
        raise ValueError('the independent oracle accepted wrong reference identity')
    action, oracles = {}, {}
    for label in ('flux','moving'):
        po.require(native[label]['state']['source_identity']==native['diagnostic_source'],'mixed source in diagnostic')
        o20, e20, r20 = audited_actions(native[label], ref, 20, deadline)
        o32, e32, r32 = audited_actions(native[label], ref, 32, deadline)
        control = compare_actions(e20,e32,ref)
        po.require(max(control.values())<=TOL,'independent split action quadrature unresolved')
        stale = stale_diffusivity_action(o32)
        po.require(stale>1e-8,'stale scalar-Cd action mutation was not detectable')
        oracles[label] = o32
        action[label] = dict(order20=r20,order32=r32,independent_refinement=control,stale_diffusivity_change=stale)
    return dict(schema='WJInnerKFixtureOperatorAssessmentV1', F2='PASS', TASK3_READY=False,
                native=po.identity(args.native),parity=po.identity(args.parity), v2=po.identity(args.v2),v2_audit=po.identity(args.v2_audit),
                v2_same_history_roles=same, v2_current_residuals=new['residuals'],
                diagnostic_source=native['diagnostic_source'],kernel_source=native['kernel_source'], reference_identity=ref.identity,
                boundary_identity=ref.boundary_identity, coefficient_errors=coefficients(native,ref),
                actions=action,common_H_flux=common_flux(native,oracles['flux']),
                moving_cutoff=moving_panels(native,oracles['moving'],oracles['flux']),physical_launch_authorized=False)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('table','native','parity','v2','v2-audit','out'):
        p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--seconds', type=float,required=True)
    args = p.parse_args()
    po.require(sys.platform.startswith('linux') and 0<args.seconds<=135 and not args.out.exists(), 'bounded independent fixture audit required')
    deadline = time.monotonic()+args.seconds
    po.require(po.identity(args.table)['sha256']==po.TABLE_SHA,'fixture boundary table changed')
    out = f2(args,read(args.table),deadline)
    out['oracle'] = po.identity(__file__)
    out['base_oracle'] = po.identity(po.__file__)
    atomic_out(args.out,out)
    print('CHANGED_INNER_K_OPERATOR_GATE_PASS',flush=True)


if __name__=='__main__':
    main()
