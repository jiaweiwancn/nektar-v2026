#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent accepted physical BE/BDF2 + EX1/EX2 audit.

The frozen steady oracle supplies current constitutive, natural/essential, and
owner-aware spatial evaluations. This module replaces its steady advection by
independently rebuilt accepted-history advection and adds independently
integrated time mass. No native operator, cached advection, or reported residual
enters the audited equations. Stored caches are comparison evidence only.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import sys
import time

import physical_oracle as po
import numpy as np
from scipy.sparse.linalg import splu

SCHEMA = 'WJH27IndependentPhysicalTransientOracleV1'
KEYS = ('flow', 'k', 'epsilon')
CACHE_LIMIT = 1e-10


def state_hash(s):
    digest = hashlib.sha256()
    for key in KEYS:
        digest.update(np.asarray(s[key], dtype='<f8').tobytes())
    return digest.hexdigest()


def state(raw, sizes):
    n, np_, nk, ne = sizes
    result = {}
    for key, count in zip(KEYS, (2*n+np_, nk, ne)):
        value = np.asarray(raw[key], dtype=float)
        po.require(value.shape == (count,) and np.all(np.isfinite(value)),
                   'invalid accepted-history '+key+' coefficients')
        result[key] = value
    return result


def combine(a, fa, b, fb):
    return {key: fa*a[key]+fb*b[key] for key in KEYS}


def scheme(data):
    h = data['input_history']
    po.require(type(data['step']) is int and type(h['step']) is int and
               h['step'] >= 0 and data['step'] == h['step']+1,
               'input history is not the immediately preceding accepted step')
    dt = data['dt']
    po.require(math.isfinite(dt) and dt > 0 and h['dt'] == dt,
               'accepted timestep must be positive and unchanged within history')
    for record in (data, h):
        expected = record['step']*dt
        po.require(math.isfinite(record['time']) and
                   abs(record['time']-expected) <= 2e-13*max(1., expected),
                   'accepted/input scientific time differs from integer step times dt')
    startup = h['step'] == 0
    return dict(startup=startup, dt=dt, alpha=(1. if startup else 1.5)/dt,
                mass_coefficients=[1., -1., 0.] if startup else [1.5, -2., .5],
                explicit_coefficients=[1., 0.] if startup else [2., -1.],
                name='BE_EX1_STARTUP' if startup else 'BDF2_EX2')


class TransientOracle:
    def __init__(self, data, reference, order, deadline):
        po.require(data.get('basis') == 'Modified_A', 'explicit Modified_A basis identity required')
        self.data, self.scheme = data, scheme(data)
        self.current = po.PhysicalOracle(data, reference, order, deadline)
        self.ref = reference
        h = data['input_history']
        self.old = state(h['current'], data['sizes'])
        self.previous = state(h['previous'], data['sizes'])
        self.cached = state(h['advection'], data['sizes'])
        self.cached_previous = state(h['previous_advection'], data['sizes'])
        self.cached_new = state(data['current_advection'], data['sizes'])
        accepted_previous = state(data['previous_advection'], data['sizes'])
        po.require(state_hash(self.old) == data['previous_sha256'] and
                   state_hash(self.cached) == data['previous_advection_sha256'] and
                   state_hash(accepted_previous) == data['previous_advection_sha256'] and
                   state_hash(self.cached_new) == data['advection_sha256'],
                   'accepted state/advection byte identity or history shift differs')
        for cache in (self.cached, self.cached_previous, self.cached_new):
            po.require(np.all(cache['flow'][2*self.current.n:] == 0.),
                       'explicit advection cache contains a pressure/continuity load')
        if self.scheme['startup']:
            po.require(state_hash(self.old) == state_hash(self.previous) and
                       state_hash(self.cached) == state_hash(self.cached_previous),
                       'BE startup did not preserve the fresh duplicated state/history')
        self.old_oracle = self.at(self.old)
        self.previous_oracle = self.at(self.previous)
        self.full_factors = None

    def at(self, coefficients):
        data = dict(self.data, coefficients=coefficients, state_sha256=state_hash(coefficients))
        return po.PhysicalOracle(data, self.ref, self.current.order, self.current.deadline)

    @staticmethod
    def fields(o, s, e, x, r):
        b, dx, dr = po.basis(e, 'u', x, r)
        ids = np.asarray(e['u'])
        U = np.column_stack((b@s['flow'][ids], b@s['flow'][o.n+ids]))
        G = np.empty((len(x), 2, 2))
        for j in range(2):
            local = s['flow'][j*o.n+ids]
            G[:, j, 0], G[:, j, 1] = dx@local, dr@local
        values = {'u': U, 'G': G}
        for key in ('k', 'epsilon'):
            if e[key]:
                q, qx, qr = po.basis(e, key, x, r)
                local = s[key][e[key]]
                values[key] = q@local
                values['grad_'+key] = np.column_stack((qx@local, qr@local))
        return values

    def transport(self, o):
        """Conservative scalar/advective velocity forms; current velocity test.

        No constitutive evaluation is needed. Quadrature panels are built from
        this history state's own H trace and cutoff, not current-state panels.
        """
        adv, awf, nonconservative = o.zeros(), o.zeros(), o.zeros()
        work = {region: dict(domain=0., sampled=0., cross_integration_bulk=0.)
                for region in ('whole', 'inner', 'outer')}
        boundary = dict(k=0., epsilon=0., cross_kinetic_flux=0., kinetic_flux=0.)
        minima = dict(k=float('inf'), outer_epsilon=float('inf'))
        rows = []
        for e in o.elements:
            o.clock(); x, r, weight = o.rule(e)
            a = self.fields(o, o.s, e, x, r)
            po.require(np.all(a['k'] > 0) and (not e['epsilon'] or np.all(a['epsilon'] > 0)),
                       'sampled historical scalar positivity failed in owner '+str(e['geometry']))
            minima['k'] = min(minima['k'], float(np.min(a['k'])))
            if e['epsilon']: minima['outer_epsilon'] = min(minima['outer_epsilon'], float(np.min(a['epsilon'])))
            test = self.fields(o, self.current.s, e, x, r)
            b = po.basis(e, 'u', x, r)[0]; ids = np.asarray(e['u'])
            conv = np.einsum('qj,qij->qi', a['u'], a['G'])
            sampled = np.where(r < self.ref.h,
                               2*o.trace(e, 'u', x)*o.trace(e, 'u', x, 1), 0.)
            for j in range(2):
                adv['flow'][j*o.n+ids] += b.T@(weight*conv[:, j])
            awf['flow'][ids] -= b.T@(weight*sampled)
            for key in ('k', 'epsilon'):
                if not e[key]: continue
                qb, qx, qr = po.basis(e, key, x, r); qi = np.asarray(e[key])
                adv[key][qi] -= qx.T@(weight*a['u'][:, 0]*a[key])+qr.T@(weight*a['u'][:, 1]*a[key])
                nonconservative[key][qi] += qb.T@(weight*np.sum(a['u']*a['grad_'+key], axis=1))
            domain = np.sum(test['u']*conv, axis=1)
            sampled_work = test['u'][:, 0]*sampled
            cross = (np.einsum('qi,qj,qij->q', a['u'], a['u'], test['G'])+
                     np.sum(test['u']*a['u'], axis=1)*np.trace(a['G'], axis1=1, axis2=2))
            for region in ('whole', 'outer' if e['epsilon'] else 'inner'):
                for key, value in (('domain', domain), ('sampled', sampled_work), ('cross_integration_bulk', cross)):
                    work[region][key] += float(weight@value)
        for key in ('u', 'k', 'epsilon'):
            for e, side in o.faces[key]:
                o.clock(); x, r, weight = o.face(e, side)
                a = self.fields(o, o.s, e, x, r)
                un = a['u']@np.asarray(po.NORMALS[side])
                if key == 'u':
                    test = self.fields(o, self.current.s, e, x, r)
                    boundary['cross_kinetic_flux'] += float(weight@(un*np.sum(test['u']*a['u'], axis=1)))
                    boundary['kinetic_flux'] += float(weight@(un*np.sum(a['u']**2, axis=1)/2))
                else:
                    flux = un*a[key]
                    adv[key][e[key]] += po.basis(e, key, x, r)[0].T@(weight*flux)
                    value = float(weight@flux); boundary[key] += value
                    rows.append(dict(field=key, geometry=e['geometry'], side=po.SIDES[side], outward_advection=value))
        return dict(advection=adv, awf=awf, total=combine(adv, 1., awf, 1.),
                    nonconservative=nonconservative, work=work, boundary=boundary, faces=rows, sampled_minima=minima)

    def storage(self, startup=None):
        o, dt = self.current, self.scheme['dt']
        if startup is None: startup = self.scheme['startup']
        # Coefficient differences are algebraically the exact BDF combination
        # and avoid cancellation of separate O(1/dt) mass actions.
        increment = combine(o.s, 1., self.old, -1.)
        old_increment = combine(self.old, 1., self.previous, -1.)
        derivative = {key: (increment[key] if startup else
                      1.5*increment[key]-.5*old_increment[key])/dt for key in KEYS}
        derivative['flow'][2*o.n:] = 0.  # No pressure time derivative.
        result = o.zeros()
        names = ('k', 'epsilon', 'velocity_work', 'kinetic_BDF_rate', 'G_energy_rate', 'nonnegative_BDF_dissipation')
        bulk = {region: {key: 0. for key in names} for region in ('whole', 'inner', 'outer')}
        for e in o.elements:
            o.clock(); x, r, weight = o.rule(e)
            a = self.fields(o, o.s, e, x, r)
            change = self.fields(o, increment, e, x, r)
            old_change = self.fields(o, old_increment, e, x, r)
            rate = self.fields(o, derivative, e, x, r)
            b = po.basis(e, 'u', x, r)[0]; ui = np.asarray(e['u'])
            for j in range(2): result['flow'][j*o.n+ui] += b.T@(weight*rate['u'][:, j])
            for key in ('k', 'epsilon'):
                if e[key]: result[key][e[key]] += po.basis(e, key, x, r)[0].T@(weight*rate[key])
            u, du, dv = a['u'], change['u'], old_change['u']
            first_energy_difference = np.sum(du*(2*u-du), axis=1)
            if startup:
                energy_rate = first_energy_difference/(2*dt)
                kinetic_bdf = energy_rate
                dissipation = np.sum(du**2, axis=1)/(2*dt)
            else:
                energy_rate = (first_energy_difference+np.sum((2*du-dv)*(2*u+dv), axis=1))/(4*dt)
                kinetic_bdf = (1.5*first_energy_difference-.5*np.sum(dv*(2*(u-du)-dv), axis=1))/(2*dt)
                dissipation = np.sum((du-dv)**2, axis=1)/(4*dt)
            values = dict(k=rate['k'], epsilon=rate.get('epsilon', np.zeros(len(x))),
                velocity_work=np.sum(u*rate['u'], axis=1), kinetic_BDF_rate=kinetic_bdf,
                G_energy_rate=energy_rate, nonnegative_BDF_dissipation=dissipation)
            for region in ('whole', 'outer' if e['epsilon'] else 'inner'):
                for key, value in values.items(): bulk[region][key] += float(weight@value)
        return result, bulk

    def full_dual(self, residual):
        """Include constrained rows when comparing saved cache inventories."""
        o = self.current
        if self.full_factors is None:
            self.full_factors = {key: splu(o.masses[key]) for key in ('u', 'k', 'epsilon')}
        result = {}
        energy = 0.
        for j in range(2):
            v = residual['flow'][j*o.n:(j+1)*o.n]
            energy += float(v@self.full_factors['u'].solve(v))
        result['momentum'] = math.sqrt(max(0., energy)/o.volume)
        for key in ('k', 'epsilon'):
            v = residual[key]
            result[key] = math.sqrt(max(0., float(v@self.full_factors[key].solve(v)))/(o.volume if key == 'k' else o.evolume))
        result['pressure_load_MAX'] = float(np.max(np.abs(residual['flow'][2*o.n:]), initial=0.))
        return result

    def assemble(self):
        o = self.current
        steady = o.assemble()
        # Keep the instantaneous steady defect as an explicitly separate
        # observation; it cannot qualify or disqualify a transient equation.
        # The physical samples and the accepted test velocity are identical
        # for equal coefficient bytes, including the duplicated BE inputs.
        transported = {}
        for owner in (o, self.old_oracle, self.previous_oracle):
            digest = owner.data['state_sha256']
            if digest not in transported:
                transported[digest] = self.transport(owner)
        current, old, previous = (transported[owner.data['state_sha256']]
            for owner in (o, self.old_oracle, self.previous_oracle))
        f, g = self.scheme['explicit_coefficients']
        explicit = combine(old['total'], f, previous['total'], g)
        mass, storage = self.storage()
        total = {key: o.total[key]-o.parts['advection'][key]-o.parts['awf_convection'][key]+
                       explicit[key]+mass[key] for key in KEYS}
        self.total, self.mass, self.explicit = total, mass, explicit
        self.transport_results = dict(current=current, old=old, previous=previous)
        checks = {}
        for name, cached, independent in (('input_current', self.cached, old['total']),
                ('input_previous', self.cached_previous, previous['total']),
                ('accepted_current', self.cached_new, current['total'])):
            delta = combine(cached, 1., independent, -1.)
            checks[name] = dict(free_dual=o.dual(delta), all_coefficient_dual=self.full_dual(delta))
        cache_ok = all(v <= CACHE_LIMIT for check in checks.values() for v in check['all_coefficient_dual'].values())
        cache_checks = dict(status='MET' if cache_ok else 'NOT_MET', limit=CACHE_LIMIT, comparisons=checks,
            scope='Every signed global coefficient, including constrained rows; independent reconstruction enters equations instead of caches')
        residual = o.dual(total)
        budgets, identities = {}, []
        bulk = steady['bulk']['whole']; boundary = steady['boundary']
        for key, source in (('k', 'k_source'), ('epsilon', 'epsilon_source')):
            fixed, free = o.fixed[key], o.free[key]
            reaction = -float(o.one[key][fixed]@total[key][fixed])
            action = float(o.one[key][free]@total[key][free])
            adv = f*old['boundary'][key]+g*previous['boundary'][key]
            rate = storage['whole'][key]
            prescribed = boundary[key+'_diffusion_prescribed']
            essential = boundary[key+'_diffusion_essential']
            physical = boundary[key+'_diffusion_actual']
            scale = steady['budgets'][key]['fixed_scale']
            weak = rate+adv+prescribed+reaction-bulk[source]
            error = abs(weak-action)/scale; identities.append(error)
            budgets[key] = dict(storage=rate, explicit_advective_outflow=adv,
                current_prescribed_outward_diffusion=prescribed, current_source=bulk[source],
                essential_outward_diffusion_reaction=reaction,
                essential_reaction_scope=steady['budgets'][key]['essential_reaction_scope'],
                current_owner_gradient_outward_diffusion_all=physical,
                current_owner_gradient_outward_diffusion_essential=essential,
                reaction_balance=weak, free_test_action=action, assembly_identity_error_over_scale=error,
                mixed_physical_balance=rate+adv+prescribed+essential-bulk[source],
                all_owner_gradient_balance=rate+adv+physical-bulk[source], fixed_scale=scale,
                actual_production_minus_destruction_not_constrained=True,
                status='OBSERVED_REQUIRES_RESOLUTION_CHECK')
        work = {key: f*old['work']['whole'][key]+g*previous['work']['whole'][key]
                for key in old['work']['whole']}
        # Replace only sampled convective AWF work. Its pressure and known-k
        # terms remain evaluated at the accepted new state.
        current_nonconvective_awf = bulk['roughness_work']-current['work']['whole']['sampled']
        total_awf = current_nonconvective_awf+work['sampled']
        rate = storage['whole']['velocity_work']
        energy_volume = rate+work['domain']+bulk['stress_work']-bulk['pressure_divergence']-total_awf
        fixed_work = sum(float(o.s['flow'][j*o.n+o.fixed[key]]@total['flow'][j*o.n+o.fixed[key]])
                         for j, key in enumerate(('u', 'w')))
        free_work = sum(float(o.s['flow'][j*o.n+ids]@total['flow'][j*o.n+ids]) for j, ids in enumerate(o.ufree))
        balance = energy_volume-boundary['traction_prescribed_work']-fixed_work
        scale = steady['budgets']['mechanical']['fixed_scale']
        error = abs(balance-free_work)/scale; identities.append(error)
        energy_identity = abs(rate-storage['whole']['G_energy_rate']-storage['whole']['nonnegative_BDF_dissipation'])/scale
        identities.append(energy_identity)
        cross_flux = f*old['boundary']['cross_kinetic_flux']+g*previous['boundary']['cross_kinetic_flux']
        cross_error = abs(work['domain']+work['cross_integration_bulk']-cross_flux)/scale
        identities.append(cross_error)
        budgets['mechanical'] = dict(storage_velocity_test_work=rate,
            physical_kinetic_energy_BDF_rate=storage['whole']['kinetic_BDF_rate'],
            G_energy_rate=storage['whole']['G_energy_rate'],
            nonnegative_BDF_dissipation=storage['whole']['nonnegative_BDF_dissipation'],
            discrete_time_energy_identity_error_over_scale=energy_identity,
            G_energy_definition='|u|^2/2 for BE; (|u_n|^2+|2u_n-u_(n-1)|^2)/4 for BDF2',
            explicit_domain_advection_current_velocity_test_work=work['domain'],
            explicit_cross_integration_boundary=cross_flux,
            explicit_cross_integration_bulk=work['cross_integration_bulk'],
            cross_advection_integration_identity_error_over_scale=cross_error,
            explicit_history_kinetic_advective_flux=f*old['boundary']['kinetic_flux']+g*previous['boundary']['kinetic_flux'],
            advection_note='Current velocity tests explicit old-state advection; history kinetic-energy outflow is not substituted for this work',
            current_full_stress_work=bulk['stress_work'], current_molecular_work=bulk['molecular_work'],
            current_raw_turbulent_production_work=bulk['P'], current_isotropic_k_divergence_work=bulk['isotropic_k_divergence'],
            current_modified_pressure_divergence_work=-bulk['pressure_divergence'],
            explicit_sampled_convective_roughness_work=work['sampled'],
            current_pressure_and_known_k_roughness_work=current_nonconvective_awf,
            sampled_roughness_work=total_awf, prescribed_boundary_traction_work=boundary['traction_prescribed_work'],
            essential_traction_reaction_work=fixed_work,
            current_essential_owner_gradient_traction_work=boundary['traction_essential_actual_work'],
            current_all_owner_gradient_traction_work=boundary['traction_actual_work'],
            reaction_balance=balance, free_velocity_test_action=free_work,
            assembly_identity_error_over_scale=error, fixed_scale=scale,
            mixed_physical_balance=energy_volume-boundary['traction_prescribed_work']-boundary['traction_essential_actual_work'],
            all_owner_gradient_balance=energy_volume-boundary['traction_actual_work'],
            status='OBSERVED_REQUIRES_RESOLUTION_CHECK')
        budgets['combined_mechanical_plus_k'] = dict(
            reaction_balance=balance+budgets['k']['reaction_balance'],
            storage=rate+storage['whole']['k'], current_molecular_work=bulk['molecular_work'],
            current_destruction=bulk['E'], raw_current_production_cancels_algebraically=True,
            current_pressure_and_isotropic_k_divergence=-bulk['pressure_divergence']+bulk['isotropic_k_divergence'],
            explicit_momentum_work_plus_scalar_outflow=work['domain']+budgets['k']['explicit_advective_outflow'],
            sampled_roughness_work=total_awf,
            scope='Discrete current-velocity BDF/EX identity; no instantaneous or reference P=E constraint')
        budgets['mass'] = steady['budgets']['mass']
        essential_actions = o.outlet_tangent_reactions(total['flow'])
        outlet = dict(policy=o.outlet_policy,
                      natural_traction=steady['outlet']['prescribed_traction_integral'],
                      natural_traction_work=steady['outlet']['traction_work'],
                      essential_transverse_reaction_rows=essential_actions,
                      essential_transverse_reaction_work=sum(
                          o.s['flow'][o.n+row['dof']]*row['value'] for row in essential_actions))
        identities.extend((steady['accounting_identity_max_error'],))
        # Nonzero fixture tests these controls; actual physical states may make
        # a particular term invisible, which is explicitly not a positive pass.
        deltas = dict(missing_time_mass=combine(mass, -1., mass, 0.),
            current_instead_of_explicit_advection=combine(current['total'], 1., explicit, -1.),
            missing_explicit_AWF_convection=combine(old['awf'], -f, previous['awf'], -g),
            missing_explicit_domain_convection=o.zeros(),
            nonconservative_explicit_scalar_advection=o.zeros(),
            stored_instead_of_reconstructed_advection=combine(combine(self.cached, f, self.cached_previous, g), 1., explicit, -1.))
        deltas['missing_explicit_domain_convection']['flow'] = -f*old['advection']['flow']-g*previous['advection']['flow']
        for key in ('k', 'epsilon'):
            deltas['nonconservative_explicit_scalar_advection'][key] = (
                f*(old['nonconservative'][key]-old['advection'][key])+
                g*(previous['nonconservative'][key]-previous['advection'][key]))
        if not self.scheme['startup']:
            deltas['EX1_instead_of_EX2'] = combine(old['total'], 1., explicit, -1.)
            be_mass, _ = self.storage(startup=True)
            deltas['BE_instead_of_BDF2_mass'] = combine(be_mass, 1., mass, -1.)
        mutations = {}
        for name, delta in deltas.items():
            sensitivity = o.dual(delta)
            mutations[name] = dict(delta_dual=sensitivity, mutated_residual=o.dual(combine(total, 1., delta, 1.)),
                status='DETECTABLE_AT_DECLARED_TOLERANCE' if any(sensitivity[key] > po.RESIDUAL_LIMITS[key] for key in sensitivity)
                else 'INACTIVE_IN_THIS_PHYSICAL_STATE; independent nonzero fixture still required')
        return dict(residuals=residual, residual_limits=po.RESIDUAL_LIMITS,
            transient_residual_status='MET' if all(residual[key] <= limit for key, limit in po.RESIDUAL_LIMITS.items()) else 'NOT_MET',
            history_cache_checks=cache_checks, storage_by_region=storage, budgets=budgets,
            accepted_discrete_outlet=outlet,
            sampled_history_minima={key: value['sampled_minima'] for key, value in self.transport_results.items()},
            accounting_identity_max_error=max(identities), accounting_identity_status='MET' if max(identities) <= 1e-10 else 'NOT_MET',
            explicit_boundary_faces=[dict(**row, extrapolation_coefficient=coefficient, history_role=role)
                for role, coefficient, values in (('input_current', f, old), ('input_previous', g, previous)) for row in values['faces']],
            current_instantaneous_spatial_observation=steady,
            instantaneous_note='Its steady residual/budgets are observed instantaneous defects, not the accepted time-discrete equations',
            mutations=mutations)

    def boundary_checks(self):
        return {name: oracle.boundary_checks() for name, oracle in
                (('accepted', self.current), ('input_current', self.old_oracle), ('input_previous', self.previous_oracle))}


def audit(data, reference, order, refine_order, deadline):
    first = TransientOracle(data, reference, order, deadline)
    coarse_residuals = first.assemble()['residuals']
    coarse_total = first.total
    del first
    final = TransientOracle(data, reference, refine_order, deadline)
    fine = final.assemble()
    difference = final.current.dual(combine(final.total, 1., coarse_total, -1.))
    quadrature_ok = all(difference[key] <= .2*limit for key, limit in po.RESIDUAL_LIMITS.items())
    boundaries = final.boundary_checks()
    boundary_ok = all(row['status'] == 'MET' for row in boundaries.values())
    result = dict(scheme=final.scheme, ex2_exercised=not final.scheme['startup'],
        input_history_identity=dict(step=data['input_history']['step'], time=data['input_history']['time'],
            current_sha256=state_hash(final.old), previous_sha256=state_hash(final.previous),
            advection_sha256=state_hash(final.cached), previous_advection_sha256=state_hash(final.cached_previous),
            provenance='Input current and advection bind the accepted output previous hashes. Earlier input state is evaluated and hashed here; cross-packet lineage requires its earlier artifact.'),
        core=fine, coarse_quadrature_residual=coarse_residuals,
        quadrature_difference_dual=difference, independent_orders=[order, refine_order],
        independent_quadrature_status='MET' if quadrature_ok else 'NOT_MET', boundary_checks=boundaries,
        accepted_history_boundary_status='MET' if boundary_ok else 'NOT_MET')
    result['accepted_time_discrete_equations_status'] = 'MET' if (
        quadrature_ok and boundary_ok and fine['transient_residual_status'] == 'MET' and
        fine['history_cache_checks']['status'] == 'MET' and fine['accounting_identity_status'] == 'MET') else 'NOT_MET'
    return final, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('native', 'table', 'out'): parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--seconds', type=float, required=True)
    parser.add_argument('--order', type=int, default=20)
    parser.add_argument('--refine-order', type=int, default=24)
    parser.add_argument('--compare-finer', type=Path)
    parser.add_argument('--address-space-gib', type=int, default=8)
    args = parser.parse_args()
    po.require(sys.platform.startswith('linux'), 'numerical transient oracle execution is Linux-only')
    po.require(math.isfinite(args.seconds) and args.seconds > 0, 'explicit finite wall cap required')
    po.require(20 <= args.order < args.refine_order <= 48, 'independent quadrature orders must satisfy 20 <= order < refined <= 48')
    po.require(not args.out.exists(), 'immutable transient oracle output already exists')
    po.require(0 < args.address_space_gib <= 20, 'invalid case-specific address-space cap')
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    cap = min([args.address_space_gib*1024**3]+[v for v in (soft, hard) if v != resource.RLIM_INFINITY])
    resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    began = time.monotonic()
    report = dict(schema=SCHEMA, status='NOT_MET', full_H27_qualification=False,
        source=po.identity(__file__), steady_oracle_source=po.identity(po.__file__), basis_source=po.identity(po.fixtures.__file__),
        native=po.identity(args.native), table=po.identity(args.table), wall_cap_seconds=args.seconds,
        memory_cap_bytes=cap, workers=1, unresolved_gates=[
            'Whole-state continuum admissibility and finite perturbations remain separate gates',
            'One accepted timestep does not establish temporal order, long-time stability, or a stationary state',
            'Owner-gradient physical budget convergence needs resolution evidence; reaction identities alone do not qualify it'])
    try:
        ref = po.load_reference(args.table)
        data, metadata = po.load_native(args.native)
        po.require(metadata.get('steady_mode') is not True, 'a steady preparation packet is not an accepted physical timestep')
        report['native_run_metadata'] = metadata
        report['source_identity'], report['state_sha256'] = data['source_identity'], data['state_sha256']
        oracle, result = audit(data, ref, args.order, args.refine_order, began+args.seconds)
        report.update(result)
        report['owner_interfaces'] = oracle.current.interfaces()
        report['profiles'] = oracle.current.profiles()
        if args.compare_finer:
            fd, fm = po.load_native(args.compare_finer)
            po.require(abs(fd['time']-data['time']) <= 2e-13*max(1., data['time']), 'refinement states are not at the same physical time')
            finer = TransientOracle(fd, ref, args.refine_order, began+args.seconds)
            fc = finer.assemble()
            comparison = oracle.current.compare(finer.current)
            comparison.update(input=po.identity(args.compare_finer), native_run_metadata=fm,
                finer_transient_residuals=fc['residuals'], finer_transient_residual_status=fc['transient_residual_status'],
                finer_history_cache_checks=fc['history_cache_checks'], finer_scheme=finer.scheme,
                physical_budget_changes={key: dict(
                    coarse=report['core']['budgets'][key]['mixed_physical_balance'],
                    finer=fc['budgets'][key]['mixed_physical_balance'],
                    difference_over_fixed_scale=(fc['budgets'][key]['mixed_physical_balance']-
                        report['core']['budgets'][key]['mixed_physical_balance'])/report['core']['budgets'][key]['fixed_scale'],
                    status='OBSERVED_REQUIRES_RESOLUTION_CHECK') for key in ('k', 'epsilon', 'mechanical')})
            if (fc['transient_residual_status'] != 'MET' or fc['history_cache_checks']['status'] != 'MET' or
                    fc['accounting_identity_status'] != 'MET'):
                comparison['status'] = 'NOT_MET'
            report['refinement'] = comparison
        else:
            report['refinement'] = dict(status='NOT_MET', reason='No same-time independently refined state supplied')
        report['reason'] = 'Accepted-step independent equations, caches, and budgets recorded; full H27 qualification remains NOT_MET.'
    except Exception as error:
        report['exception'], report['reason'] = type(error).__name__, str(error)
    report['elapsed_seconds'] = time.monotonic()-began
    report['maximum_resident_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False); stream.write('\n')
    print(json.dumps({key: report[key] for key in ('status', 'reason', 'accepted_time_discrete_equations_status',
                     'elapsed_seconds', 'maximum_resident_bytes', 'exception') if key in report}))
    if 'exception' in report: raise SystemExit(2)


if __name__ == '__main__':
    main()
