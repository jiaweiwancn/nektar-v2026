#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Linux-only independent continuous coupled-MMS jet/forcing/flux verifier.

Analytic jets are differentiated symbolically from the prospective formulas.
Physical stress/flux divergence uses two fourth-order finite-difference steps,
independent WJ eigenvalue roots and NumPy Chebyshev cutoff evaluation. No native
AD, sampled native basis, matrix actions or time-discrete residuals are reused.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'BLIS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[variable] = '1'
import numpy as np
from numpy.polynomial.chebyshev import chebval
import sympy as sp
import reference_oracle

JET_ORDER = ['value', 'dx', 'dr', 'dt', 'dxx', 'dxr', 'drr']
LABELS = {'below_cutoff', 'near_cutoff_below', 'near_cutoff_above', 'chi_zero',
          'collar_lower', 'collar_middle', 'roughness_below', 'roughness_above',
          'collar_upper', 'H_minus', 'H_trace', 'H_plus', 'outer_middle', 'outer_high'}
FORWARD_RADIAL = {'collar_lower', 'collar_upper', 'H_trace'}
VARIANTS = ('Historical', 'PolynomialScalarsV1')


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def finite(value, shape=None):
    result = np.asarray(value, dtype=float)
    require(np.all(np.isfinite(result)) and (shape is None or result.shape == shape),
            'native value has wrong shape or nonfinite data')
    return result


def check(errors, name, actual, expected, scale, tolerance=1e-10):
    expected = finite(expected); actual = finite(actual, expected.shape)
    require(scale > 0 and math.isfinite(scale), 'invalid declared comparison scale')
    error = float(np.max(abs(actual-expected), initial=0.)/scale)
    errors[name] = max(errors.get(name, 0.), error)
    require(error <= tolerance, f'{name} error {error:.6g} exceeds {tolerance:.6g}')


class IndependentMMS:
    def __init__(self, table, body_amplitude=0., variant='Historical'):
        require(variant in VARIANTS, 'unsupported independent continuous MMS variant')
        self.variant = variant
        self.body_amplitude = float(body_amplitude)
        require(math.isfinite(self.body_amplitude), 'nonfinite prescribed body amplitude')
        p, c = table['parameters'], table['constants']
        self.Kin, self.H, self.h, self.nu = p['K'], p['H'], p['h'], p['nu']
        self.cd, self.sk, self.se = p['cd'], p['sigma_k'], 1.3
        self.cl, self.qe, self.kappa = c['c_length'], c['q_e'], c['kappa']
        self.d0, self.R = table['inlet_binding']['d0'], .5
        self.Ain, self.tau0 = self.Kin**1.5/self.cl, self.Kin/self.qe
        self.coefficients = table['cutoff_table']['coefficients']
        self.lo, self.hi = table['cutoff_table']['domain']
        x, r, t = sp.symbols('x r t', real=True)
        decimal = lambda value: sp.Rational(str(value))
        H, R, kin, cl = (sp.Float(v, 17) for v in (self.H, self.R, self.Kin, self.cl))
        s = sp.sqrt(kin)*(1+decimal(.012)*sp.sin(decimal(.6)*t)+
                         decimal(.008)*x*sp.cos(decimal(.8)*t))
        eta = decimal(.005)*(1+decimal(.2)*sp.sin(decimal(.7)*t))
        psi = (decimal(.02)*(1+decimal(.1)*sp.cos(decimal(.9)*t))*r*r/2+
               decimal(.05)*(1+decimal(.2)*sp.sin(decimal(.7)*t))*x*x*(1-x)*r*r*(R-r))
        expressions = dict(u=sp.diff(psi, r), w=-sp.diff(psi, x),
            Pi=decimal(.0002)*(1+decimal(.3)*sp.sin(decimal(.8)*t))*(x*x+r*r/R**2-sp.Rational(2,3)),
            k=s*s*(1+eta*r*r*(r-H)**2/(H*H*R*R)), K=s*s,
            epsilon=s**3/(cl*r)*(1+decimal(.01)*sp.sin(decimal(.5)*t)*(r-H)**2/R**2))
        if variant == 'PolynomialScalarsV1':
            expressions['k'] = s*s*(1+eta*r*(r-H)/(H*R))
            expressions['epsilon'] = s**3/(cl*H)*(1+decimal(.01)*sp.sin(decimal(.5)*t)*(r-H)**2/R**2)
        self.functions = {}
        for name, expression in expressions.items():
            jet = [expression, sp.diff(expression, x), sp.diff(expression, r),
                   sp.diff(expression, t), sp.diff(expression, x, 2),
                   sp.diff(expression, x, r), sp.diff(expression, r, 2)]
            self.functions[name] = sp.lambdify((x, r, t), jet, modules='numpy', cse=True)

    def jets(self, x, r, t):
        # Outer epsilon is evaluated only where it is transported; no inner PDE.
        return {name: np.asarray(function(x, r, t), dtype=float)
                for name, function in self.functions.items() if name != 'epsilon' or self.outer(r)}

    def outer(self, r):
        return r >= self.H-1e-14*self.H

    def cutoff(self, K):
        require(self.lo <= K <= self.hi, 'independent K left qualified interval')
        return self.H*math.exp(float(chebval((2*K-self.lo-self.hi)/(self.hi-self.lo), self.coefficients)))

    def auxiliary(self, K, r):
        rough = self.h*math.sqrt(K)/(70*self.nu)
        exponent = max(.5-.4*rough**.7, 1-.79*rough**(-.28))
        intercept = -10.7*math.expm1(exponent*math.log(rough))
        zv = self.nu*intercept/math.sqrt(K)
        require(zv < 0, 'independent roughness branch changed')
        return self.kappa*math.sqrt(K/self.qe)*(r-zv)

    def point(self, x, r, t):
        j = self.jets(x, r, t); K, k = j['K'][0], j['k'][0]
        cutoff = self.cutoff(K); A = K**1.5/self.cl
        outer = self.outer(r); E = j['epsilon'][0] if outer else A/max(r, cutoff)
        auxiliary = self.auxiliary(K, r)
        gradient = np.zeros((3, 3))
        gradient[0, 0], gradient[0, 2] = j['u'][1:3]
        gradient[2, 0], gradient[2, 2] = j['w'][1:3]
        strain = gradient+gradient.T
        n = np.array([0., 0., 1.]); projector = np.eye(3)-np.outer(n, n)
        wall_shear = -auxiliary*projector@strain@n
        isotropic = 2*k*np.eye(3)/3+np.outer(wall_shear, n)+np.outer(n, wall_shear)
        coordinate = min(1., max(0., (r/self.H-.4)/.35))
        chi = coordinate**3*(10-15*coordinate+6*coordinate**2)
        reynolds = ((1-chi)*isotropic+chi*reference_oracle.full_wj(k, E, gradient)
                    if chi > 0 else isotropic)
        turbulent_mu = 0.
        if chi > 0:
            st = (k/E)*(strain/2-np.trace(gradient)*np.eye(3)/3)
            ot = (k/E)*(gradient-gradient.T)/2
            iw = float(np.trace(ot@ot))
            root = reference_oracle.physical_root(float(np.trace(st@st)), iw)
            turbulent_mu = 3*root/(5*(root*root-2*iw))*k*k/E
        mu = self.nu+(1-chi)*auxiliary+chi*turbulent_mu
        stress = self.nu*strain-reynolds+2*k*np.eye(3)/3
        production = -float(np.sum(reynolds*gradient))
        dk = self.nu+self.cd*k*k/(self.sk*E)
        de = self.nu+self.cd*k*k/(self.se*E) if outer else 0.
        flux_k = dk*j['k'][1:3]
        flux_e = de*j['epsilon'][1:3] if outer else np.zeros(2)
        require(np.all(np.isfinite(stress)) and k > 0 and E > 0 and dk > 0,
                'nonfinite/nonpositive independent continuous state')
        return dict(jets=j, gradient=gradient, R=reynolds, V=stress, chi=chi, mu=mu, A=A,
                    cutoff=cutoff, auxiliary=auxiliary, production=production, E=E,
                    Dk=dk, De=de, flux_k=flux_k, flux_e=flux_e, outer=outer,
                    bundle=np.r_[stress.ravel(), flux_k, flux_e])

    def body(self, x, r, t):
        """Prescribed analytic physical force; scalar or broadcast x/r arrays."""
        x, r = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(r, dtype=float))
        if self.body_amplitude == 0:
            return np.zeros((3,)+x.shape)
        return np.stack((self.body_amplitude*2e-4*(1+x+.5*r)*np.cos(.3*t),
                         np.zeros_like(x),
                         self.body_amplitude*3e-4*(1+r)*np.sin(.2*t)))

    def awf_H(self, x, t):
        j = self.jets(x, self.H, t)
        value = np.array([2*j['u'][0]*j['u'][1]+j['Pi'][1]-2*j['k'][1]/3, 0., 0.])
        if self.body_amplitude != 0:
            value[0] -= self.body(x, self.H, t)[0]
        return value

    def region(self, x, r, t):
        K = self.functions['K'](x, r, t)[0]
        return (r <= self.cutoff(K), 0 if r <= .4*self.H else 1 if r < .75*self.H else 2,
                self.outer(r))

    def derivative(self, x, r, t, axis, step, forward=False):
        offsets = [0, 1, 2, 3, 4] if forward else [-2, -1, 1, 2]
        coefficients = np.array([-25, 48, -36, 16, -3] if forward else [1, -8, 8, -1], float)/12
        # At a C2 collar endpoint take the smooth right-hand derivative; at H
        # epsilon requires the outer owner. Central stencils stay in one region.
        expected = self.region(x, r+1e-9*self.H if forward else r, t)
        if axis == 0:
            expected = self.region(x, r, t)
        for _ in range(40):
            sites = [(x+(offset*step if axis == 0 else 0.),
                      r+(offset*step if axis == 1 else 0.)) for offset in offsets]
            safe = all(0 <= xx <= 1 and 0 <= rr <= self.R and
                       (offset == 0 or self.region(xx, rr, t) == expected)
                       for offset, (xx, rr) in zip(offsets, sites))
            if safe:
                values = np.array([self.point(xx, rr, t)['bundle'] for xx, rr in sites])
                return coefficients@values/step, step
            step /= 2
        raise ValueError('unable to construct a smooth-region finite-difference stencil')

    def forcing(self, x, r, t, step, radial_forward):
        state = self.point(x, r, t); j = state['jets']
        steps = (step, step) if np.isscalar(step) else step
        dx, hx = self.derivative(x, r, t, 0, steps[0])
        dr, hr = self.derivative(x, r, t, 1, steps[1], radial_forward)
        divergence_stress = dx[:9].reshape(3, 3)[:, 0]+dr[:9].reshape(3, 3)[:, 2]
        awf = self.awf_H(x, t) if r <= self.h else np.zeros(3)
        u, w = j['u'][0], j['w'][0]; divergence = j['u'][1]+j['w'][2]
        momentum = -divergence_stress-awf
        for component, name in ((0, 'u'), (2, 'w')):
            momentum[component] += j[name][3]+u*j[name][1]+w*j[name][2]+j['Pi'][1 if component == 0 else 2]
        if self.body_amplitude != 0:
            momentum -= self.body(x, r, t)
        def transport(q):
            return q[3]+u*q[1]+w*q[2]+q[0]*divergence
        fk = transport(j['k'])-dx[9]-dr[10]-state['production']+state['E']
        fe = 0.
        if state['outer']:
            e = j['epsilon']
            fe = (transport(e)-dx[11]-dr[12]-(1.44*state['production']-1.82*e[0])*e[0]/j['k'][0])
        return np.r_[momentum, fk, fe], (hx, hr)


def verify(table_path, probe_path, report):
    table = json.loads(table_path.read_text()); probe = json.loads(probe_path.read_text())
    require(table['status'] == 'PASS' and probe['schema'] == 'WJH27ContinuousCoupledMMSProbeV1',
            'wrong/unqualified table or native probe schema')
    require(probe['table_sha256'] == hashlib.sha256(table_path.read_bytes()).hexdigest(),
            'native/reference table identity differs')
    require(probe['formula_identity'] == 'WJH27ContinuousCoupledMMSV1' and
            probe['model_identity'] == 'FiniteNuIsotropicWallFullWJCollarV2' and
            probe['jet_order'] == JET_ORDER and probe['full_coupled_gate_pass'] is False and
            probe['admissibility_qualified'] is False, 'wrong formula/model/scope identity')
    old_formula_keys = {'s', 'k', 'epsilon_outer', 'psi', 'Pi', 'momentum', 'physical_body_at_H'}
    require(old_formula_keys <= set(probe['source_formula']) <= old_formula_keys | {'physical_body'},
            'native prospective source formula metadata missing')
    amplitude = float(probe.get('body_amplitude', 0.))
    if 'body_amplitude' in probe:
        require(probe.get('body_formula_identity') == 'WJH27PrescribedPhysicalBodyMMSV1' and
                'physical_body' in probe['source_formula'], 'unbound prescribed body formula')
    variant = probe.get('mms_variant', 'Historical')
    m = IndependentMMS(table, amplitude, variant); errors = report['maximum_reference_scaled_errors'] = {}
    report['mms_variant'] = variant
    report['variant_provenance'] = ('probe.mms_variant' if 'mms_variant' in probe
                                    else 'Historical pre-variant probe ABI; no polynomial-variant claim.')
    expected_k = ('s^2*(1+.005*(1+.2*sin(.7*t))*r^2*(r-H)^2/(H^2*R^2))'
                  if variant == 'Historical' else 's^2*(1+.005*(1+.2*sin(.7*t))*r*(r-H)/(H*R))')
    expected_e = ('s^3/(cl*r)*(1+.01*sin(.5*t)*(r-H)^2/R^2)'
                  if variant == 'Historical' else 's^3/(cl*H)*(1+.01*sin(.5*t)*(r-H)^2/R^2)')
    require(probe['source_formula']['k'] == expected_k and probe['source_formula']['epsilon_outer'] == expected_e,
            'scalar formula metadata does not match the declared MMS variant')
    report['body_amplitude'] = amplitude
    check(errors, 'parameter_body_amplitude', probe['parameters'].get('body_amplitude', 0.), amplitude, 1., 2e-12)
    body_controls = dict(missing_bulk_body=0., missing_sampled_H_body=0.)
    parameter_map = dict(Kin=m.Kin, H=m.H, h=m.h, d0=m.d0, R=m.R, nu=m.nu,
                         Cd=m.cd, sigmaK=m.sk, sigmaE=m.se, cl=m.cl, qe=m.qe,
                         kappa=m.kappa, tau0=m.tau0)
    for name, value in parameter_map.items():
        check(errors, 'parameter_'+name, probe['parameters'][name], value, abs(value), 2e-12)
    require(len(probe['states']) == 84, 'native continuous state coverage differs')
    expected = {(x, t, label) for x in (.137, .61) for t in (0., .13, .41) for label in LABELS}
    require({(v['x'], v['time'], v['label']) for v in probe['states']} == expected,
            'native continuous state labels/locations/times differ')
    convergence = report['finite_difference_forcing'] = []
    for row in probe['states']:
        x, r, t = row['x'], row['r'], row['time']; state = m.point(x, r, t); j = state['jets']
        c = state['cutoff']
        radii = dict(below_cutoff=.5*c, near_cutoff_below=.95*c, near_cutoff_above=1.05*c,
                     chi_zero=.3*m.H, collar_lower=.4*m.H, collar_middle=.5*m.H,
                     roughness_below=m.h-1e-4, roughness_above=m.h+1e-4,
                     collar_upper=.75*m.H, H_minus=m.H-1e-4, H_trace=m.H,
                     H_plus=m.H+1e-4, outer_middle=.3, outer_high=.47)
        check(errors, 'sample_radius', r, radii[row['label']], m.H, 2e-12)
        require(row['outer'] == state['outer'], 'wrong native outer ownership')
        check(errors, 'physical_z_translation', row['z'], m.d0+r, m.R)
        for name in ('u', 'w', 'Pi', 'k', 'K'):
            scale = m.Kin if name in ('k', 'K') else m.tau0 if name == 'Pi' else 1.
            check(errors, name+'_jet', row[name], j[name], scale)
        epsilon_scale = m.Ain/r if state['outer'] else m.Ain/max(r, m.cutoff(m.Kin))
        if state['outer']:
            check(errors, 'epsilon_jet', row['epsilon'], j['epsilon'], epsilon_scale)
        else:
            require(row['epsilon'] is None and row['epsilon_outward_flux_x_r'] is None,
                    'native probe invented a transported inner epsilon')
        check(errors, 'gradient', row['gradient'], state['gradient'], 1.)
        for name, column in (('gradient_x', (4, 5)), ('gradient_r', (5, 6))):
            g = np.zeros((3, 3)); g[0, 0], g[0, 2] = j['u'][list(column)]
            g[2, 0], g[2, 2] = j['w'][list(column)]
            check(errors, name, row[name], g, 1.)
        for key, value, scale in (
                ('A', state['A'], m.Ain), ('cutoff', state['cutoff'], m.H),
                ('auxiliary', state['auxiliary'], m.auxiliary(m.Kin, r)),
                ('mu', state['mu'], m.nu+m.auxiliary(m.Kin, r)),
                ('chi', state['chi'], 1.), ('R', state['R'], m.Kin),
                ('V', state['V'], m.tau0), ('production', state['production'], m.Ain/m.H),
                ('destruction', state['E'], epsilon_scale),
                ('diffusion_k', state['Dk'], m.nu+m.cd*m.Kin*m.Kin*m.H/m.Ain),
                ('diffusion_e', state['De'], m.nu+m.cd*m.Kin*m.Kin*m.H/(m.se*m.Ain))):
            check(errors, key, row[key], value, scale)
        require(float(finite(row['mu'])) > 0, 'native positive split coefficient is nonpositive')
        body, body_H = m.body(x, r, t), m.body(x, m.H, t)
        require(amplitude == 0 or 'physical_body' in row, 'nonzero-body probe omitted bulk physical body')
        check(errors, 'physical_body', row.get('physical_body', np.zeros(3)), body, 1.)
        check(errors, 'physical_body_at_H', row['physical_body_at_H'], body_H, 1.)
        check(errors, 'awf_at_H', row['awf_at_H'], m.awf_H(x, t), 1.)
        check(errors, 'awf', row['awf'], m.awf_H(x, t) if r <= m.h else np.zeros(3), 1.)
        check(errors, 'k_outward_flux', row['k_outward_flux_x_r'], -state['flux_k'], m.Kin)
        if state['outer']:
            check(errors, 'epsilon_outward_flux', row['epsilon_outward_flux_x_r'], -state['flux_e'], epsilon_scale)
        native_force = np.r_[finite(row['momentum_forcing'], (3,)), row['k_forcing'], row['epsilon_forcing']]
        scales = np.array([1., 1., 1., m.Kin, epsilon_scale])
        coarse, steps0 = m.forcing(x, r, t, 1e-3, row['label'] in FORWARD_RADIAL)
        # Halve the accepted physical steps, not the initial requested values:
        # otherwise branch restrictions could make the two radial steps equal.
        fine, steps1 = m.forcing(x, r, t, tuple(v/2 for v in steps0), row['label'] in FORWARD_RADIAL)
        require(all(abs(b-a/2) <= 1e-16 for a, b in zip(steps0, steps1)),
                'FD refinement did not halve both accepted physical steps')
        error0 = float(np.max(abs(coarse-native_force)/scales))
        error1 = float(np.max(abs(fine-native_force)/scales))
        refinement = float(np.max(abs(fine-coarse)/scales))
        require(error1 <= 1e-8 and refinement <= 2e-7,
                f'continuous forcing mismatch at {row["label"]}: fine={error1}, delta={refinement}')
        # Change just one source term while keeping native manufactured fields
        # and all other body terms fixed. These defects cannot cancel each other.
        missing_bulk = native_force.copy(); missing_bulk[:3] += body
        missing_sampled = native_force.copy()
        if r <= m.h:
            missing_sampled[0] -= body_H[0]
        body_controls['missing_bulk_body'] = max(body_controls['missing_bulk_body'],
            float(np.max(abs(missing_bulk-fine)/scales)))
        body_controls['missing_sampled_H_body'] = max(body_controls['missing_sampled_H_body'],
            float(np.max(abs(missing_sampled-fine)/scales)))
        resolved = error1 <= 1e-10
        require(resolved or error1 <= .6*error0,
                'forcing finite differences neither converge nor reach the declared roundoff band')
        convergence.append(dict(label=row['label'], x=x, r=r, time=t,
            coarse_steps=list(steps0), fine_steps=list(steps1),
            radial_scheme='forward-fourth-order' if row['label'] in FORWARD_RADIAL else 'central-fourth-order',
            coarse_reference_scaled_error=error0, fine_reference_scaled_error=error1,
            reference_scaled_refinement_difference=refinement,
            roundoff_resolved=resolved))
    require(len(probe['boundary_fluxes']) == 24, 'physical boundary flux sample coverage differs')
    expected_boundary = {(x, t, label) for x in (.137, .61) for t in (0., .13, .41)
                         for label in ('ground', 'top', 'epsilon_H')}
    expected_boundary |= {(x, t, label) for x, label in ((0., 'left'), (1., 'right'))
                          for t in (0., .13, .41)}
    require({(v['x'], v['time'], v['label']) for v in probe['boundary_fluxes']} == expected_boundary,
            'boundary sample faces/locations/times differ')
    for row in probe['boundary_fluxes']:
        state = m.point(row['x'], row['r'], row['time']); normal = finite(row['normal'], (3,))
        radius = dict(ground=0., top=m.R, epsilon_H=m.H, left=.3, right=.3)[row['label']]
        normal_expected = dict(ground=[0, 0, -1], top=[0, 0, 1], epsilon_H=[0, 0, -1],
                               left=[-1, 0, 0], right=[1, 0, 0])[row['label']]
        check(errors, 'boundary_radius', row['r'], radius, m.R)
        check(errors, 'boundary_normal', normal, normal_expected, 1.)
        check(errors, 'boundary_normal_length', normal@normal, 1., 1.)
        planar = normal[[0, 2]]
        check(errors, 'boundary_k_flux', row['k_flux'], -state['flux_k']@planar, m.Kin)
        if row['label'] == 'ground' and variant == 'Historical':
            check(errors, 'ground_zero_k_flux', row['k_flux'], 0., m.Kin)
        if row['label'] == 'ground' and variant == 'PolynomialScalarsV1':
            eta = .005*(1+.2*math.sin(.7*row['time']))
            expected_gradient = -state['jets']['K'][0]*eta/m.R
            check(errors, 'polynomial_ground_k_gradient', state['jets']['k'][2], expected_gradient, m.Kin)
            check(errors, 'polynomial_ground_k_flux', row['k_flux'], state['Dk']*expected_gradient, m.Kin)
        if state['outer']:
            check(errors, 'boundary_epsilon_flux', row['epsilon_flux'], -state['flux_e']@planar, m.Ain/row['r'])
        else:
            require(row['epsilon_flux'] is None, 'inner boundary has transported epsilon flux')
        if row['label'] == 'epsilon_H':
            check(errors, 'dependent_H_trace', state['jets']['epsilon'][0], state['A']/m.H, m.Ain/m.H)
            check(errors, 'H_k_trace', state['jets']['k'][0], state['jets']['K'][0], m.Kin)
            if variant == 'PolynomialScalarsV1':
                check(errors, 'polynomial_epsilon_H_normal_derivative', state['jets']['epsilon'][2], 0., m.Ain/m.H)
    report['prescribed_body_negative_controls'] = {
        name: dict(maximum_mutant_forcing_error=error,
                   status='DETECTED' if amplitude != 0 and error > 10e-8 else 'NOT_EXERCISED_ZERO_OR_UNRESOLVED_BODY')
        for name, error in body_controls.items()}
    report['nonzero_prescribed_body_coverage_pass'] = amplitude != 0 and all(
        value['status'] == 'DETECTED' for value in report['prescribed_body_negative_controls'].values())
    require(amplitude == 0 or report['nonzero_prescribed_body_coverage_pass'],
            'nonzero body fixture does not resolve both missing-bulk and missing-H-body mutations')
    report.update(status='PASS_CONTINUOUS_FORCING_ONLY',
        state_count=len(probe['states']), boundary_count=len(probe['boundary_fluxes']),
        maximum_fine_forcing_error=max(v['fine_reference_scaled_error'] for v in convergence),
        derivative_method='Independent SymPy analytic jets; fourth-order FD of complete V and D grad(q), with current K/c/nuA at every displaced point.',
        scales='Fixed physical Lref=tref=Uref=1; k/K use Kin, Pi/V use tau0, outer epsilon uses Ain/r; forcing uses Uref/tref, Kin/tref and Ain/(r*tref). No alpha normalization.',
        tolerance=dict(point_values_jets_fluxes=1e-10, fine_forcing=1e-8,
                       roundoff_resolved_band=1e-10, otherwise_required_error_ratio=.6))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, required=True)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--output', '--out', dest='output', type=Path, required=True)
    args = parser.parse_args()
    report = dict(schema='WJH27ContinuousCoupledMMSCheckV1', status='FAIL_CONTINUOUS_FORCING_CHECK',
        full_coupled_gate_pass=False, admissibility_qualified=False, platform=sys.platform,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        independent_WJ_source_sha256=hashlib.sha256(Path(reference_oracle.__file__).read_bytes()).hexdigest(),
        numerical_libraries=dict(numpy=np.__version__, sympy=sp.__version__),
        table=str(args.table), probe=str(args.probe),
        scope='Continuous field/forcing/flux independence only; no coupled space-time order or application qualification.')
    failure = None
    try:
        require(sys.platform.startswith('linux'), 'authoritative continuous MMS verification requires Linux')
        report['table_sha256'] = hashlib.sha256(args.table.read_bytes()).hexdigest()
        report['probe_sha256'] = hashlib.sha256(args.probe.read_bytes()).hexdigest()
        verify(args.table, args.probe, report)
    except Exception as error:
        failure = error; report.update(exception=type(error).__name__, reason=str(error))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False)+'\n')
    if failure is not None:
        raise failure
    print(json.dumps(dict(status=report['status'], maximum_value_error=max(report['maximum_reference_scaled_errors'].values()),
                         maximum_fine_forcing_error=report['maximum_fine_forcing_error']), sort_keys=True))


if __name__ == '__main__':
    main()
