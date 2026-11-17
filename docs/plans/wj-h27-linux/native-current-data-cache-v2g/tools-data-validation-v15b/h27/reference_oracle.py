#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent H27 finite-nu/collar counterfactual column oracle.

Run numerical work on Linux under the H27 resource runner. This program does
not import any solver implementation, old cutoff, or old reference quadrature.
It is a feasibility diagnostic, not a 2D steady/application solver or a proof
of global uniqueness/principal ellipticity. Failures are retained as JSON.
"""
import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

CP = 1.8
CD = 572.0 / 6561.0
QE = 1.0 / math.sqrt(CD)
KAPPA = math.sqrt(1.3 * (1.82 - 1.44) * math.sqrt(CD))
CL = KAPPA * QE**1.5
Q_PEAK = 2.0 * CP / 0.65


@dataclass(frozen=True)
class Parameters:
    K: float = 0.008809015145515192
    H: float = 0.15
    h: float = 0.10854503542535461
    nu: float = 6.423173803526447e-5
    cd: float = CD
    sigma_k: float = 1.0
    endpoint_flux_over_A: float = 0.0
    inner_diffusion: str = "scalar"
    rtol: float = 2e-11


def physical_root(s, w):
    """Polynomial eigenvalue roots, independently of native bracket solver."""
    if not (math.isfinite(s) and math.isfinite(w) and s >= 0 and w <= 0):
        raise ValueError("invalid WJ invariants")
    roots = np.roots([1.0, -CP, -2.0 * w - 2.7 * s, 2.0 * CP * w])
    selected = [float(z.real) for z in roots
                if abs(z.imag) < 1e-10 * max(1.0, abs(z.real))
                and z.real >= CP * (1 - 1e-12)]
    if len(selected) != 1:
        raise ValueError("physical WJ root absent or ambiguous")
    return selected[0]


def shear_law(q):
    n = physical_root(q*q/2.0, -q*q/2.0)
    denominator = n*n + q*q
    c = 0.6*n/denominator
    fn = 3*n*n - 2*CP*n - 0.35*q*q
    dn = 2*q*(0.35*n + CP)/fn
    dc = 0.6*(dn*denominator - n*(2*n*dn + 2*q))/denominator**2
    return c*q, c + q*dc, n


def full_wj(k, e, gradient):
    """Full 3D explicit WJ tensor with DevStrainRawWorkV1."""
    tau = k/e
    identity = np.eye(3)
    s = tau*((gradient + gradient.T)/2 - np.trace(gradient)*identity/3)
    o = tau*(gradient - gradient.T)/2
    s2, o2 = s@s, o@o
    invs, invw = np.trace(s2), np.trace(o2)
    inv4 = np.trace(s@o2)
    n = physical_root(float(invs), float(invw))
    denominator = (5.0/6)*(n*n - 2*invw)*(2*n*n - invw)
    b1 = -n*(2*n*n - 7*invw)/denominator
    b3 = -12*inv4/(n*denominator)
    b4 = -2*(n*n - 2*invw)/denominator
    b6 = -6*n/denominator
    b9 = 6/denominator
    a = (b1*s + b3*(o2-invw*identity/3) + b4*(s@o-o@s)
         + b6*(s@o2+o2@s-2*inv4*identity/3)
         + b9*(o@s@o2-o2@s@o))
    return k*(2*identity/3+a)


def planar_chi_zero_principal(k, e, g, nu, auxiliary, theta):
    """Exact 2D incompressible rank-one quadratic; no tensor FD used here.

    n=(cos(theta),0,sin(theta)), l=(-sin(theta),0,cos(theta));
    H=l tensor n. Scalars and sampled K are frozen in this principal block.
    """
    q = k*g/e
    n = physical_root(q*q/2, -q*q/2)
    ratio = CP/n
    denominator = 0.7+2.65*ratio-2*ratio*ratio
    a0 = (0.7+0.65*ratio)/denominator
    a1 = (0.7-0.05*ratio-2*ratio*ratio)/denominator
    viscosity = (k/g)*shear_law(q)[0]
    x, s = math.cos(2*theta), math.sin(2*theta)
    value = (nu+auxiliary*x*x+viscosity*s*s
             -viscosity*(q/n)*s*(a0+a1*x))
    return dict(value=value, q=q, n=n, cp_over_N=ratio, a0=a0, a1=a1,
                wj_viscosity=viscosity,
                positive_split_value=nu+auxiliary*x*x)


def principal_obstruction_audit(source):
    """Cheap independent verification using one completed reference result."""
    payload = json.loads(source.read_text())
    case = payload['cases'][0]
    p = Parameters(**case['parameters'])
    column = Column(p)
    if not column.zv < 0:
        raise ValueError('ground/inlet bound requires documented negative-intercept branch')
    accepted = [root for root in case['roots'] if 'audit' in root]
    if len(accepted) != 1:
        raise ValueError('principal witness audit requires one resolved reference root')
    audit = accepted[0]['audit']
    ground = next(row for row in audit['states'] if row['y'] == 0)
    auxiliary = column.auxiliary(0)
    d = p.nu+auxiliary
    # Independent reconstruction of the RETAINED inlet; this historical
    # formula is used only to inspect unchanged boundary/initial data.
    def old_integrated_p_over_A(r):
        return math.log1p(column.a*r/d)+p.nu/(d+column.a*r)-p.nu/d
    old_c = p.H*math.exp(1-old_integrated_p_over_A(p.H))
    def old_Q(r):
        sink = r/old_c if r <= old_c else 1+math.log(r/old_c)
        return column.A*(sink-old_integrated_p_over_A(r))
    inlet_k = p.K-sum(quad(lambda r: old_Q(r)/(d+column.a*r), lo, hi,
                          epsabs=1e-14*p.K, epsrel=1e-12)[0]
                      for lo, hi in [(0, old_c), (old_c, p.H)])
    inlet_g = column.t/d
    states = [dict(label='selected_reference_ground', k=ground['k'], e=ground['e'], g=ground['g']),
              dict(label='retained_inlet_with_new_dependent_E', k=inlet_k,
                   e=column.A/audit['c'], g=inlet_g)]
    results = []
    for state in states:
        k, e, g = state['k'], state['e'], state['g']
        gradient = np.zeros((3, 3)); gradient[0, 2] = g
        witness_rows = []
        for degrees in (15., 20., 25., 30., 35., 45.):
            theta = math.radians(degrees)
            analytic = planar_chi_zero_principal(k, e, g, p.nu, auxiliary, theta)
            normal = np.array([math.cos(theta), 0, math.sin(theta)])
            polarization = np.array([-math.sin(theta), 0, math.cos(theta)])
            direction = np.outer(polarization, normal)
            fds = []
            for relative_step in (1e-4, 2e-5, 1e-5):
                delta = relative_step*e/k
                derivative = (column.tensor(0, k, e, gradient+delta*direction)
                              -column.tensor(0, k, e, gradient-delta*direction))/(2*delta)
                tangent = p.nu*(direction+direction.T)-derivative
                value = float(np.sum(direction*tangent))
                fds.append(dict(relative_step=relative_step, value=value,
                                normalized_error=(value-analytic['value'])/d))
            witness_rows.append(dict(degrees=degrees, analytic=analytic,
                                     normalized=analytic['value']/d, full_3d_tensor_FD=fds))
        grid = np.linspace(0, math.pi, 721)
        function = lambda angle: planar_chi_zero_principal(k, e, g, p.nu, auxiliary, angle)['value']
        index = int(np.argmin([function(angle) for angle in grid]))
        minimum = minimize_scalar(function, bounds=(grid[max(0, index-1)], grid[min(720, index+1)]),
                                  method='bounded', options={'xatol': 1e-14})
        # Closed upper bound for the single theta=30-degree witness.
        # For cp/N<=1/10: f<=3/10, a0>=17/21, a1>=5/7.
        bound = p.nu+auxiliary/4+(k/g)*(9/40-7*math.sqrt(3)/30)
        q = k*g/e
        results.append(dict(**state, q=q, witnesses=witness_rows,
                            theta_min_degrees=math.degrees(float(minimum.x)),
                            planar_minimum=float(minimum.fun), planar_minimum_over_mu=float(minimum.fun)/d,
                            rigorous_theta30_upper_bound=bound,
                            bound_requires_q_at_least=18*math.sqrt(2),
                            bound_condition_met=q >= 18*math.sqrt(2),
                            certified_negative_2d_witness=q >= 18*math.sqrt(2) and bound < 0))
    # For ANY positive reference k diffusion with Q(0)=Q(H)=0 on the
    # increasing-shear branch: P<E on the tail collar (see analytical note).
    # Below .4H, P is fixed by auxiliary shear. This bounds the possible
    # production budget and prevents making c arbitrarily small via diffusion.
    p_upper_over_A = old_integrated_p_over_A(0.4*p.H)+math.log(2.5)
    c_lower = p.H*math.exp(1-p_upper_over_A)
    q_lower = inlet_k*inlet_g*c_lower/column.A
    inlet_bound = p.nu+auxiliary/4+(inlet_k/inlet_g)*(9/40-7*math.sqrt(3)/30)
    uniform = dict(production_upper_bound_over_A=p_upper_over_A,
                   cutoff_lower_bound=c_lower, cutoff_lower_bound_over_H=c_lower/p.H,
                   retained_inlet_q_lower_bound=q_lower,
                   retained_inlet_theta30_upper_bound=inlet_bound,
                   normalized_upper_bound=inlet_bound/d,
                   hypotheses=['same K,A,t,nu,h,H and geometry-only chi',
                               'positive inner k diffusivity, reference source P-E',
                               'zero reference endpoint k fluxes',
                               'selected increasing-shear branch where chi>0',
                               'negative viscosity intercept',
                               'retained inlet k and velocity profile'],
                   bound_conditions_met=c_lower <= 0.4*p.H and q_lower >= 18*math.sqrt(2),
                   diffusion_change_cannot_repair_retained_inlet=(
                       c_lower <= 0.4*p.H and q_lower >= 18*math.sqrt(2) and inlet_bound < 0))
    return dict(schema='WJH27PlanarPrincipalObstructionV1', source=str(source), parameters=asdict(p),
                auxiliary_nu_ground=auxiliary, normalization_mu=p.nu+auxiliary,
                retained_inlet_ground_k=inlet_k, retained_inlet_ground_k_over_K=inlet_k/p.K,
                retained_inlet_old_cutoff=old_c, states=results, uniform_diffusion_bound=uniform,
                restrictions='No changed stress, EP, floor, clipping, roughness or actual-flow budget constraint.')


class Column:
    def __init__(self, p):
        self.p = p
        if not (p.K > 0 and p.H > p.h >= 0 and p.nu > 0
                and p.cd > 0 and p.sigma_k > 0 and p.rtol > 0):
            raise ValueError("invalid positive scientific input or roughness geometry")
        self.A = p.K**1.5/CL
        self.t = p.K/QE
        self.us = math.sqrt(self.t)
        self.a = KAPPA*self.us
        x = p.h*math.sqrt(p.K)/(70*p.nu)
        b = 10.7 if x == 0 else -10.7*math.expm1(
            max(0.5-0.4*x**0.7, 1-0.79*x**(-0.28))*math.log(x))
        self.zv = p.nu*b/math.sqrt(p.K)

    def chi(self, y):
        if y <= 0.4:
            return 0.0
        if y >= 0.75:
            return 1.0
        s = (y-0.4)/0.35
        return s**3*(10-15*s+6*s*s)

    def auxiliary(self, y):
        return self.a*max(self.p.H*y-self.zv, 0.0)

    def state(self, y, ratio, logc):
        p = self.p
        if not math.isfinite(ratio) or ratio <= 0:
            raise ValueError("nonpositive reference k")
        k = p.K*ratio
        e = self.A/(p.H*max(y, math.exp(logc)))
        chi, auxiliary = self.chi(y), self.auxiliary(y)
        alpha = p.nu+(1-chi)*auxiliary
        if chi == 0:
            g = self.t/alpha
            q = k*g/e
        else:
            def traction(q):
                return alpha*e*q/k + chi*k*shear_law(q)[0] - self.t
            # This bounded increasing-shear branch must be declared. No jump
            # to another shear root is allowed on an inadmissible trial.
            if traction(Q_PEAK) < 0:
                raise ValueError("no root on declared q<=q_peak shear branch")
            q = brentq(traction, 0, Q_PEAK, xtol=2e-14, rtol=2e-14)
            g = e*q/k
        f, fp, n = shear_law(q)
        shear = (1-chi)*auxiliary*g + chi*k*f
        production = shear*g
        diffusion = (p.nu+p.cd*k*k/(p.sigma_k*e)
                     if p.inner_diffusion == "scalar"
                     else p.nu+auxiliary/p.sigma_k)
        tangent = alpha+chi*k*k/e*fp
        return dict(k=k, e=e, q=q, n=n, g=g, chi=chi, auxiliary=auxiliary,
                    production=production, diffusion=diffusion,
                    shear_tangent=tangent, shear_function=f, shear_derivative=fp,
                    traction_error=p.nu*g+shear-self.t)

    def solve_trial(self, logc):
        if not -12 <= logc <= 0:
            raise ValueError("cutoff outside declared [exp(-12)H,H] range")
        p = self.p
        breaks = sorted({0.0, 0.4, 0.75, 1.0, math.exp(logc),
                         min(1.0, max(0.0, self.zv/p.H))}, reverse=True)
        # Last two states are analytic derivatives of k/K and Q/A with
        # respect to log(c/H), independently checked by centered differences.
        value = np.array([1.0, p.endpoint_flux_over_A, 0.0, 0.0, 0.0])
        panels = []
        def rhs(y, v):
            state = self.state(float(y), float(v[0]), logc)
            k, e, q = state['k'], state['e'], state['q']
            f, fp, tangent = (state['shear_function'], state['shear_derivative'],
                              state['shear_tangent'])
            de = -e if y < math.exp(logc) else 0.0
            dg = (-state['chi']*(f+q*fp)*p.K*v[3]/tangent
                  +state['chi']*k*q*fp*de/(e*tangent))
            dd = ((state['diffusion']-p.nu)*(2*v[3]/v[0]-de/e)
                  if p.inner_diffusion == 'scalar' else 0.0)
            return [p.H*self.A*v[1]/(p.K*state['diffusion']),
                    p.H*(state['e']-state['production'])/self.A,
                    p.H*state['g']/self.us,
                    p.H*self.A*(v[4]*state['diffusion']-v[1]*dd)
                    /(p.K*state['diffusion']**2),
                    p.H*(de-(self.t-2*p.nu*state['g'])*dg)/self.A]
        for upper, lower in zip(breaks[:-1], breaks[1:]):
            solution = solve_ivp(rhs, (upper, lower), value, method="DOP853",
                                 rtol=p.rtol, atol=p.rtol*0.04,
                                 max_step=0.025, dense_output=True)
            if not solution.success:
                raise ValueError("reference IVP: "+solution.message)
            value = solution.y[:, -1]
            panels.append((lower, upper, solution))
        def sample(y):
            for lower, upper, solution in panels:
                if lower-1e-14 <= y <= upper+1e-14:
                    return solution.sol(y)
            raise ValueError("sample outside column")
        return dict(logc=float(logc), c=p.H*math.exp(logc), ground_flux=float(value[1]),
                    anchor=float(-value[2]*self.us), sample=sample,
                    ground_flux_derivative=float(value[4]),
                    breaks=sorted(breaks))

    def tensor(self, y, k, e, gradient):
        rw = full_wj(k, e, gradient)
        rh = rw.copy()
        shear = -self.auxiliary(y)*(gradient+gradient.T)
        rh[0, 2] = rh[2, 0] = shear[0, 2]
        rh[1, 2] = rh[2, 1] = shear[1, 2]
        return (1-self.chi(y))*rh+self.chi(y)*rw

    def audit(self, trial, principal_directions):
        p = self.p
        logc, sample, breaks = trial['logc'], trial['sample'], trial['breaks']
        def integration(fun):
            return sum(quad(fun, a, b, epsabs=self.A*2e-12,
                            epsrel=2e-11, limit=100)[0]
                       for a, b in zip(breaks[:-1], breaks[1:]))
        integrated_p = p.H*integration(lambda y: self.state(y, sample(y)[0], logc)['production'])
        integrated_e = p.H*integration(lambda y: self.state(y, sample(y)[0], logc)['e'])
        rows = []
        for y in sorted(set(np.linspace(0, 1, 65).tolist()+breaks)):
            state = self.state(y, sample(y)[0], logc)
            gradient = np.zeros((3, 3)); gradient[0, 2] = state['g']
            r = self.tensor(y, state['k'], state['e'], gradient)
            state.update(y=float(y), flux_over_A=float(sample(y)[1]),
                         min_R_over_k=float(np.linalg.eigvalsh(r).min()/state['k']),
                         contraction_error_over_E=float((-np.sum(r*gradient)-state['production'])/state['e']))
            rows.append(state)
        principal = []
        # Finite-difference acoustic tensor on both polarizations transverse
        # to each normal; includes oblique/out-of-plane directions. Sampling
        # supplies witnesses, not a continuum lower-bound certificate.
        for y in sorted(set([0., .2, .4, .5, .6, .75, 1., math.exp(logc)])):
            state = self.state(y, sample(y)[0], logc)
            g = np.zeros((3, 3)); g[0, 2] = state['g']
            scale = p.nu+(1-state['chi'])*state['auxiliary'] + state['chi']*state['k']*shear_law(state['q'])[0]/state['g']
            for step in (2e-5, 1e-5):
                best = (float('inf'), None)
                delta = step*state['e']/state['k']
                for i in range(principal_directions):
                    nz = 1-2*(i+.5)/principal_directions
                    angle = i*math.pi*(3-math.sqrt(5))
                    normal = np.array([math.sqrt(1-nz*nz)*math.cos(angle),
                                       math.sqrt(1-nz*nz)*math.sin(angle), nz])
                    axis = np.eye(3)[int(np.argmin(abs(normal)))]
                    l1 = np.cross(normal, axis); l1 /= np.linalg.norm(l1)
                    basis = [l1, np.cross(normal, l1)]
                    acoustic = np.zeros((2, 2))
                    for j, polarization in enumerate(basis):
                        direction = np.outer(polarization, normal)
                        dr = (self.tensor(y, state['k'], state['e'], g+delta*direction)
                              -self.tensor(y, state['k'], state['e'], g-delta*direction))/(2*delta)
                        action = p.nu*(direction+direction.T)-dr
                        for irow, test in enumerate(basis):
                            acoustic[irow, j] = test@action@normal
                    quality = float(np.linalg.eigvalsh((acoustic+acoustic.T)/2).min()/scale)
                    if quality < best[0]:
                        best = (quality, normal.tolist())
                principal.append(dict(y=y, step=step, min_sampled_normalized=best[0], normal=best[1]))
        delta = 2e-4
        derivatives = []
        for hstep in (delta, delta/2):
            derivative = (self.solve_trial(logc+hstep)['ground_flux']
                          -self.solve_trial(logc-hstep)['ground_flux'])/(2*hstep)
            derivatives.append(dict(log_step=hstep, ground_flux_derivative=derivative,
                                    difference_from_sensitivity=derivative-trial['ground_flux_derivative'],
                                    inverse_absolute_slope=1/abs(derivative) if derivative else None))
        moments = []
        for degree in range(6):
            integral = integration(lambda y: (
                -self.A*float(sample(y)[1])*degree*y**max(0, degree-1)
                +p.H*(self.state(y, sample(y)[0], logc)['production']
                      -self.state(y, sample(y)[0], logc)['e'])*y**degree))
            # Natural endpoint term is Q(H)phi(H)-Q(0)phi(0).
            integral += self.A*p.endpoint_flux_over_A
            integral -= self.A*trial['ground_flux']*(1.0 if degree == 0 else 0.0)
            moments.append(integral/self.A)
        return dict(c=trial['c'], c_over_H=math.exp(logc), anchor_U=trial['anchor'],
                    ground_flux_over_A=trial['ground_flux'],
                    production_integral=integrated_p, dissipation_integral=integrated_e,
                    budget_error_over_A=(integrated_e-integrated_p)/self.A-p.endpoint_flux_over_A,
                    analytic_sink_error_over_A=integrated_e/self.A-(1-logc),
                    epsilon_H=self.A/p.H, root_derivatives=derivatives,
                    analytic_ground_flux_derivative=trial['ground_flux_derivative'],
                    weak_moments_over_A=moments, states=rows, principal_samples=principal)


def run(p, args):
    column = Column(p)
    scans, roots = [], []
    for logc in np.linspace(args.log_min, args.log_max, args.scan_points):
        try:
            trial = column.solve_trial(float(logc))
            scans.append(dict(logc=float(logc), ground_flux=trial['ground_flux'],
                              ground_flux_derivative=trial['ground_flux_derivative'],
                              admissible=True,
                              scope='positive-k bounded-shear IVP only; tensor audit follows at roots'))
        except (ValueError, RuntimeError, OverflowError) as exc:
            scans.append(dict(logc=float(logc), admissible=False, reason=str(exc)))
    for left, right in zip(scans[:-1], scans[1:]):
        if (left['admissible'] and right['admissible']
                and left['ground_flux']*right['ground_flux'] <= 0):
            try:
                root = brentq(lambda c: column.solve_trial(c)['ground_flux'],
                              left['logc'], right['logc'], xtol=2e-12, rtol=1e-13)
                if all('logc' not in prior or abs(root-prior['logc']) > 1e-8 for prior in roots):
                    roots.append(dict(logc=float(root), bracket=[left['logc'], right['logc']],
                                      audit=column.audit(column.solve_trial(root), args.principal_directions)))
            except (ValueError, RuntimeError, OverflowError) as exc:
                roots.append(dict(bracket=[left['logc'], right['logc']], failure=str(exc)))
    valid = [r for r in roots if 'audit' in r]
    return dict(parameters=asdict(p), A=column.A, traction=column.t, zv=column.zv,
                status="ONE_SCANNED_ROOT" if len(valid) == 1 else "NO_ROOT_OR_MULTIPLE_ROOTS",
                interpretation="Finite scan is not a global uniqueness or principal certificate.",
                scan=scans, roots=roots)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--principal-audit-input', type=Path,
                        help='Audit planar obstruction from an existing oracle result; no column solve')
    parser.add_argument('--k-factors', default='1')
    parser.add_argument('--endpoint-flux-over-A', type=float, default=0.0)
    parser.add_argument('--inner-diffusion', choices=['scalar', 'auxiliary'], default='scalar')
    parser.add_argument('--log-min', type=float, default=-4.0)
    parser.add_argument('--log-max', type=float, default=-0.2)
    parser.add_argument('--scan-points', type=int, default=41)
    parser.add_argument('--principal-directions', type=int, default=128)
    parser.add_argument('--rtol', type=float, default=2e-11)
    args = parser.parse_args()
    if args.principal_audit_input:
        payload = principal_obstruction_audit(args.principal_audit_input)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('x') as stream:
            json.dump(payload, stream, indent=2, allow_nan=False)
            stream.write('\n')
        print(json.dumps(dict(output=str(args.output),
                              certified_negative_2d_witnesses=[state['certified_negative_2d_witness']
                                                              for state in payload['states']],
                              uniform_diffusion_obstruction=payload['uniform_diffusion_bound']
                              ['diffusion_change_cannot_repair_retained_inlet'])))
        return
    if not (-12 <= args.log_min < args.log_max <= 0 and args.scan_points >= 3
            and args.principal_directions >= 8):
        parser.error('invalid bounded scan or direction inventory')
    cases = []
    for factor in map(float, args.k_factors.split(',')):
        try:
            cases.append(run(Parameters(K=Parameters.K*factor,
                             endpoint_flux_over_A=args.endpoint_flux_over_A,
                             inner_diffusion=args.inner_diffusion, rtol=args.rtol), args))
        except (ValueError, RuntimeError, OverflowError) as exc:
            cases.append(dict(k_factor=factor, status='REJECTED', reason=str(exc)))
    payload = dict(schema='WJH27IndependentReferenceOracleV1',
                   scope='local counterfactual calibration column, not application steady flow',
                   constants=dict(cp=CP, cd=CD, qe=QE, kappa=KAPPA, c_length=CL, q_peak=Q_PEAK),
                   cases=cases)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(output=str(args.output), statuses=[c['status'] for c in cases])))


if __name__ == '__main__':
    main()
