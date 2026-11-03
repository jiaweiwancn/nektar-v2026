#!/usr/bin/env python3
"""One preregistered M8-R1 phase-D mathematical feasibility screen.

This is a research driver, not an installed closure.  It runs no solver, does
not alter archived receipts, and does not infer state-set invariance from a
constitutive certificate.  Numerical evaluation is restricted to Linux and
one explicitly assigned CPU.  The existing M8 interval and indexed oracles
are imported only from the immutable ab5 archive source.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

sys.dont_write_bytecode = True
THREAD_KEYS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
               'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS')
for _key in THREAD_KEYS:
    os.environ[_key] = '1'
import numpy as np

BASELINE = 'ab5c3c6fee19761e8d88692e8e6c331e5889234e'
RHO = .8
BA = .09*2.55


def stamp():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(source, *args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()


def plain(value):
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, (Decimal, F)):
        return str(value)
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    return value


def save(path, value):
    Path(path).write_text(json.dumps(plain(value), indent=2, allow_nan=False)+'\n')


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def relative_error(left, right):
    return float(np.max(np.abs(np.asarray(left)-np.asarray(right)))/
                 max(1., float(np.max(np.abs(left))), float(np.max(np.abs(right)))))


def auxiliary(K, r, nu, roughness, branch=None):
    """Original physical nuA and all open-branch K/r derivatives.

    Real parts select the same analytic branch during complex-step checks.
    At a max or positive-part knot, callers must request the appropriate
    one-sided action rather than interpreting this branch value as a Jacobian.
    """
    if np.real(K) <= 0 or nu <= 0 or roughness < 0:
        raise ValueError('nuA requires positive sampled K/nu and nonnegative h')
    root = np.sqrt(K)
    if roughness == 0:
        b, bk, active, gap = 10.7, 0., 'smooth', None
    else:
        x = roughness*root/(70*nu)
        p1, p2 = .5-.4*x**.7, 1-.79*x**(-.28)
        active = branch or ('first' if p1.real >= p2.real else 'second')
        p = p1 if active == 'first' else p2
        px = -.28*x**(-.3) if active == 'first' else .2212*x**(-1.28)
        pk = px*x/(2*K)
        b = 10.7*(1-x**p)
        bk = -10.7*x**p*(pk*np.log(x)+p/(2*K))
        gap = float(np.real(p1-p2))
    zv = nu*b/root
    zvk = nu*(bk/root-b/(2*K*root))
    positive = np.real(r-zv) > 0
    value = BA*root*(r-zv) if positive else 0.*K
    derivative_k = BA*(r/(2*root)-nu*bk) if positive else 0.*K
    derivative_r = BA*root if positive else 0.*K
    return dict(value=value, derivative_k=derivative_k, derivative_r=derivative_r,
                zv=zv, zv_k=zvk, branch=active, branch_gap=gap, positive=positive)


def auxiliary_action(K, r, nu, roughness, dK, dr, branch=None, at_positive_knot=False):
    base = auxiliary(K, r, nu, roughness, branch)
    if at_positive_knot:
        return BA*math.sqrt(K)*max(dr-float(np.real(base['zv_k']))*dK, 0.)
    return base['derivative_k']*dK+base['derivative_r']*dr


def inner_value(G, k, nua):
    """Direct tensor law; supports complex perturbations on open branches."""
    G = np.asarray(G)
    S = G+G.T
    c = 2*k/3
    s = np.sqrt(np.sum(S*S))
    if np.real(k) < 0 or np.real(nua) < 0:
        raise ValueError('Negative k or nuA is outside the inner contract')
    if k == 0 or s == 0:
        turbulent = np.zeros((3, 3), dtype=np.result_type(G, k, nua))
        return dict(F=turbulent, R=c*np.eye(3)-turbulent, S=S, s=s,
                    c=c, a=nua*s, t=0., ta=1. if k > 0 else 0., tc=0.,
                    branch='zero-strain' if k > 0 else 'zero-k')
    a = nua*s
    if np.real(a) <= RHO*np.real(c):
        t, ta, tc, active = a, 1., 0., 'inactive'
    else:
        z = (a-RHO*c)/((1-RHO)*c)
        q = np.exp(-z)
        t, ta, tc, active = c*(1-(1-RHO)*q), q, 1-q*(1+(1-RHO)*z), 'active'
    turbulent = (t/s)*S
    return dict(F=turbulent, R=c*np.eye(3)-turbulent, S=S, s=s,
                c=c, a=a, t=t, ta=ta, tc=tc, branch=active)


def inner_action(G, k, nua, H, dk=0., dnua=0.):
    value = inner_value(G, k, nua)
    dS = np.asarray(H)+np.asarray(H).T
    if k <= 0:
        raise ValueError('Use the declared admissible directional limit at zero k')
    if value['s'] == 0:
        dF = nua*dS
    else:
        S, s, ta = value['S'], value['s'], value['ta']
        h = value['t']/s
        dF = (h*dS+(nua*ta-h)*np.sum(S*dS)*S/s**2
              +ta*dnua*S+value['tc']*(2*dk/3)*S/s)
    return dF, 2*dk/3*np.eye(3)-dF


def outer_full_action(G, k, epsilon, H, dk, depsilon, branch):
    """Independent analytic full-state action of the one fixed outer law.

    All C' factors and the homogeneous quadratic/cubic factors are included.
    Gamma and Sref are prescribed spatial data, hence have no k/epsilon action.
    """
    eye = np.eye(3)
    ratio, dratio = k/epsilon, dk/epsilon-k*depsilon/epsilon**2
    s, o = ratio*(G+G.T), ratio*(G-G.T)
    ds = ratio*(H+H.T)+dratio*(G+G.T)
    do = ratio*(H-H.T)+dratio*(G-G.T)
    ss, oo = np.sum(s*s), np.sum(o*o)
    dss, doo = 2*np.sum(s*ds), 2*np.sum(o*do)
    q = math.sqrt((ss if branch == 'strain' else oo)/2)
    dq = (dss if branch == 'strain' else doo)/(4*q)
    aa = .35*q**1.5
    uu = .36*math.exp(.75*q)
    ee = math.exp(-uu)
    C = .3/(1+aa)*(1-ee)
    Cp = -.3*.525*math.sqrt(q)/(1+aa)**2*(1-ee)+.225*uu*ee/(1+aa)
    dC = Cp*dq
    quadratic = -.1*(s@s-eye*ss/3)+.1*(o@s-s@o)+.26*(o@o.T-eye*oo/3)
    dquad = (-.1*(ds@s+s@ds-eye*dss/3)+.1*(do@s+o@ds-ds@o-s@do)
             +.26*(do@o.T+o@do.T-eye*doo/3))
    cubic = -10*(s@s@o-o@s@s)+5*s*(oo-ss)
    dsquare = ds@s+s@ds
    dcubic = (-10*(dsquare@o+s@s@do-do@s@s-o@dsquare)
              +5*ds*(oo-ss)+5*s*(doo-dss))
    base = 2/3*eye-C*s+C*quadratic+.4*C**3*cubic
    return (dk*base+k*(-dC*s-C*ds+dC*quadratic+C*dquad
                      +.4*(3*C*C*dC*cubic+C**3*dcubic))-(H+H.T))


def inner_suite(core, oracle, ledger, check, steps):
    results, aux_results, archived = [], [], []
    Qiv = core.Iv(core.Q_LO, core.Q_HI)
    qref = float((Qiv.lo+Qiv.hi)/2)
    reference_ratio = 3*core.Iv(2).sqrt()/(2*Qiv)
    tracefree_factor = core.Iv(1)-(core.Iv(F(2, 3))).sqrt()
    check('inner_open_reference_margin', reference_ratio.hi < Decimal('0.8'),
          dict(ratio_interval=reference_ratio.pair(), rho='4/5'))
    check('inner_tracefree_covariance_factor', tracefree_factor.lo > 0,
          dict(factor_interval=tracefree_factor.pair()))
    H = np.array([[.2, -.3, .1], [.4, -.1, .2], [-.2, .3, -.1]])
    fixed = []
    for kr in (.5, 1., 1.5):
        K = qref*.4**2*kr
        for r in (.001, .5, 2.):
            for h in (0., .1):
                aux = auxiliary(K, r, 1.5e-5, h)
                nua = float(np.real(aux['value']))
                G = np.zeros((3, 3)); G[0, 2] = (K/qref)/(1.5e-5+nua)
                fixed.append((f'reference_K{kr}_r{r}_h{h}', G, K, K, r, h, .03*K, -.07*K, .013, True))
    fixed.extend([
        ('three_dimensional_active', np.array([[.2, -.3, 3.], [.4, -.1, -.2], [.1, .3, -.1]]),
         .5, .51, .5, .1, -.04, .07, .03, False),
        ('three_dimensional_inactive', np.array([[.02, -.03, .1], [.04, -.01, -.02], [.01, .03, -.01]]),
         .5, .51, .5, .1, -.04, .07, .03, False),
        ('weak_divergence', np.array([[.22, -.3, 3.], [.4, -.1, -.2], [.1, .3, -.1]]),
         .5, .51, .5, .1, -.04, .07, .03, False)])
    for name, G, k, K, r, roughness, dk, dK, dr, reference in fixed:
        nu = 1.5e-5
        aux = auxiliary(K, r, nu, roughness)
        nua = float(np.real(aux['value']))
        dnua = float(np.real(auxiliary_action(K, r, nu, roughness, dK, dr)))
        value = inner_value(G, k, nua)
        action, action_R = inner_action(G, k, nua, H, dk, dnua)
        ladder = []
        for step in steps:
            plus_aux = auxiliary(K+step*dK, r+step*dr, nu, roughness)['value']
            minus_aux = auxiliary(K-step*dK, r-step*dr, nu, roughness)['value']
            plus = inner_value(G+step*H, k+step*dk, plus_aux)['F']
            minus = inner_value(G-step*H, k-step*dk, minus_aux)['F']
            ladder.append(dict(step=step, relative_error=relative_error((plus-minus)/(2*step), action)))
        cs = 1e-30
        pert_aux = auxiliary(K+1j*cs*dK, r+1j*cs*dr, nu, roughness)['value']
        cs_action = inner_value(G+1j*cs*H, k+1j*cs*dk, pert_aux)['F'].imag/cs
        cs_error = relative_error(cs_action, action)
        check(name+'_full_k_K_r_strain_derivative', cs_error < 2e-12 and min(
              point['relative_error'] for point in ladder) < 2e-7,
              dict(complex_step_error=cs_error, difference_ladder=ladder))
        eigen = np.linalg.eigvalsh(np.asarray(value['R'], float))
        S, s, c, t = value['S'], float(value['s']), float(value['c']), float(value['t'])
        general_margin = c-t
        check(name+'_covariance', eigen[0] >= general_margin-2e-13 and eigen[0] > 0,
              dict(eigenvalues=eigen.tolist(), general_lower=general_margin))
        if abs(np.trace(S)) < 1e-13 and s:
            lower = float(tracefree_factor.lo)*c
            check(name+'_tracefree_margin', eigen[0] >= lower-2e-13, dict(lower=lower))
        if reference:
            check(name+'_reference_exact_inactivity', value['branch'] in ('inactive', 'zero-strain') and
                  relative_error(value['F'], nua*S) < 1e-15,
                  dict(branch=value['branch'], a_over_c=float(value['a']/value['c'])))
        # Symmetric tensor coordinates are orthonormal; this checks the full
        # fixed-field tangent, while the radial theorem supplies all directions.
        basis = []
        for i in range(3):
            b = np.zeros((3, 3)); b[i, i] = 1.; basis.append(b)
        for i, j in ((0, 1), (0, 2), (1, 2)):
            b = np.zeros((3, 3)); b[i, j] = b[j, i] = 1/math.sqrt(2); basis.append(b)
        tangent = np.array([[np.sum(b*inner_action(G, k, nua, a/2)[0]) for a in basis] for b in basis])
        tangent_eigen = np.linalg.eigvalsh((tangent+tangent.T)/2)
        radial = nua*float(value['ta'])
        check(name+'_full_symmetric_monotone_tangent', tangent_eigen[0] >= radial-2e-13 and
              relative_error(tangent, tangent.T) < 1e-14,
              dict(eigenvalues=tangent_eigen.tolist(), radial_lower=radial,
                   molecular_principal_lower=nu+radial))
        production = -float(np.sum(value['R']*G))
        expected = t*s/2-c*float(np.trace(G))
        # This is precisely the stress work in momentum, including weak div U.
        turbulent_work = float(np.sum(value['R']*G))
        check(name+'_same_stress_energy_exchange', abs(production-expected) < 2e-12 and
              abs(production+turbulent_work) < 1e-15,
              dict(production=production, weak_divergence_term=-c*float(np.trace(G))))
        n = np.array([0., 0., 1.])
        traction = -1.3*n+nu*S@n-value['R']@n
        shifted = -(1.3+c)*n+nu*S@n+value['F']@n
        check(name+'_same_normal_and_tangential_traction', relative_error(traction, shifted) < 2e-15,
              dict(normal_traction=float(n@traction), traction=traction.tolist()))
        results.append(dict(name=name, k=k, sampled_K=K, r=r, roughness=roughness, nu=nu,
            nu_A=nua, branch=value['branch'], a_over_c=float(value['a']/value['c']),
            minimum_covariance=float(eigen[0]), general_covariance_lower=general_margin,
            radial_tangent=radial, tangential_tangent=float(value['t']/value['s']),
            complex_step_error=cs_error, difference_ladder=ladder))
    # The fixed roughness controls exercise both max branches, with all K
    # dependence differentiated.  No fitting/search of a constitutive parameter.
    for x in (.01, .1, 1., 10., 10000.):
        K, nu = .51, 1.5e-5
        roughness = x*70*nu/math.sqrt(K)
        aux = auxiliary(K, .5, nu, roughness)
        cs = auxiliary(K+1j*1e-30, .5, nu, roughness)['value'].imag/1e-30
        error = abs(cs-aux['derivative_k'])/max(1., abs(cs))
        check(f'roughness_x{x}_sampled_K_derivative', error < 2e-13, dict(error=float(error)))
        aux_results.append(dict(x=x, branch=aux['branch'], branch_gap=aux['branch_gap'],
                                nu_A=float(aux['value']), dnuA_dK=float(aux['derivative_k']), error=float(error)))
    check('both_roughness_max_branches_exercised', {v['branch'] for v in aux_results} == {'first', 'second'})
    # Explicit max-knot one-sided controls: bisection locates a known scalar
    # branch knot; it does not tune rho or the model.
    left, right = .01, 10.
    def gap(x): return .5-.4*x**.7-(1-.79*x**(-.28))
    for _ in range(80):
        midpoint = (left+right)/2
        if gap(midpoint) > 0: left = midpoint
        else: right = midpoint
    x = (left+right)/2; K = .51; nu = 1.5e-5; roughness = x*70*nu/math.sqrt(K)
    knot_records = []
    for direction, branch in ((1., 'second'), (-1., 'first')):
        base = auxiliary(K, .5, nu, roughness, branch)
        analytic = float(auxiliary_action(K, .5, nu, roughness, direction, 0., branch))
        errors = []
        for step in steps:
            numerical = (auxiliary(K+step*direction, .5, nu, roughness)['value']-base['value'])/step
            errors.append(dict(step=step, error=float(abs(numerical-analytic))))
        check('roughness_max_knot_'+branch, errors[-1]['error'] < 2e-5, dict(ladder=errors))
        knot_records.append(dict(direction=direction, limiting_branch=branch, derivative=analytic, ladder=errors))
    # Positive part at the smooth-wall cutoff: the negative side is exactly zero.
    K = .51; nu = 1.5e-5; cutoff = float(auxiliary(K, 1., nu, 0.)['zv'])
    positive_part = []
    for dr in (-1., 1.):
        derivative = auxiliary_action(K, cutoff, nu, 0., 0., dr, at_positive_knot=True)
        errors = [dict(step=h, error=float(abs((auxiliary(K, cutoff+h*dr, nu, 0.)['value']
                  -auxiliary(K, cutoff, nu, 0.)['value'])/h-derivative))) for h in steps]
        check('auxiliary_positive_part_'+str(dr), errors[-1]['error'] < 1e-12, dict(ladder=errors))
        positive_part.append(dict(dr=dr, derivative=float(derivative), ladder=errors))
    # C1 saturation join, zero-strain/nuA/k, and the explicitly nonsmooth joint origin.
    shear = np.zeros((3, 3)); shear[0, 2] = 1.
    k, nua = .5, .1
    join = shear*(RHO*(2*k/3)/(nua*math.sqrt(2)))
    join_actions = []
    for sign in (-1., 1.):
        analytic = inner_action(join, k, nua, sign*shear)[0]
        errors = [dict(step=h, error=relative_error((inner_value(join+h*sign*shear, k, nua)['F']-
                  inner_value(join, k, nua)['F'])/h, analytic)) for h in steps]
        check('saturation_C1_join_'+str(sign), errors[-1]['error'] < 3e-6, dict(ladder=errors))
        join_actions.append(dict(side=sign, ladder=errors))
    zeroS = inner_action(np.zeros((3, 3)), k, nua, shear)[0]
    check('zero_strain_derivative', relative_error(zeroS, nua*(shear+shear.T)) < 1e-15)
    zero_nu = inner_action(shear, k, 0., np.zeros((3, 3)), dnua=1.)[0]
    check('zero_nuA_admissible_derivative', relative_error(zero_nu, shear+shear.T) < 1e-15)
    zero_k = inner_value(shear, 0., nua)
    check('zero_k_R_and_production', np.max(np.abs(zero_k['R'])) == 0)
    zero_k_limit = (2/3)*(shear+shear.T)/math.sqrt(2)
    zero_k_errors = [dict(step=h, error=relative_error(inner_value(shear, h, nua)['F']/h,
                                                       zero_k_limit)) for h in steps]
    check('zero_k_one_sided_limit_fixed_strain', zero_k_errors[-1]['error'] < 1e-12,
          dict(ladder=zero_k_errors))
    # Homogeneity makes the joint-origin directional derivative explicit.
    joint = inner_value(shear, .1, nua)['F']
    separate = inner_value(np.zeros((3, 3)), .1, nua)['F']+inner_value(shear, 0., nua)['F']
    defect = float(np.linalg.norm(joint-separate))
    check('joint_origin_non_Frechet_limitation_reproduced', defect > .01,
          dict(additivity_defect=defect, full_Frechet_derivative_claimed=False))
    # Every archived detailed inner minimum is reconsidered, with source pointer.
    def visit(value, pointer=''):
        if isinstance(value, dict):
            if {'k', 'sampled_K', 'nu_A'} <= set(value) and ('G' in value or 'gradient' in value):
                G = np.asarray(value.get('G', value.get('gradient')), float)
                if G.shape == (3, 3):
                    old = 2*float(value['k'])/3*np.eye(3)-float(value['nu_A'])*(G+G.T)
                    new = inner_value(G, float(value['k']), float(value['nu_A']))
                    eigen = np.linalg.eigvalsh(new['R'])
                    check('archived_inner_'+pointer, eigen[0] >= -1e-14,
                          dict(old_min=float(np.linalg.eigvalsh(old)[0]), new_min=float(eigen[0])))
                    archived.append(dict(pointer=pointer, original=value, old_min=float(np.linalg.eigvalsh(old)[0]),
                        new_min=float(eigen[0]), new_eigenvalues=eigen.tolist(), branch=new['branch'],
                        a_over_c=float(new['a']/new['c']) if new['c'] > 0 else None,
                        new_stress=new['R'].tolist()))
            for key, child in value.items(): visit(child, pointer+'/'+str(key))
        elif isinstance(value, list):
            for index, child in enumerate(value): visit(child, pointer+'/'+str(index))
    visit(ledger)
    check('archived_inner_detailed_witnesses_loaded', len(archived) > 0, dict(count=len(archived)))
    return dict(formal_interval_constants=dict(Qstar=Qiv.pair(), open_reference_a_over_c=reference_ratio.pair(),
                open_reference_margin=(core.Iv(F(4, 5))-reference_ratio).pair(),
                tracefree_covariance_factor=tracefree_factor.pair()),
        fixed_cases=results, auxiliary_derivatives=aux_results, auxiliary_max_knot=knot_records,
        auxiliary_positive_part=positive_part, saturation_join=join_actions,
        zero_k_ladder=zero_k_errors, joint_origin_additivity_defect=defect,
        archived_inner_minima=archived,
        continuum_scope='Fixed positive local k and physical nonnegative nuA: monotone symmetric strain map; all rank-one inner physical momentum forms retain molecular nu; covariance positive at k>0. Cross-field blocks are fully differentiated but are not asserted to form a positive coupled Jacobian.',
        limitations=['No full Frechet derivative at joint k/strain or k/nuA origins.',
                    'No uniform turbulent radial tangent lower bound under unbounded strain or k depletion.',
                    'No full discrete positivity, source-map invariance or nonlinear AWF reconstruction implemented.'])


def source_suite(oracle, check):
    records = []
    weights = np.array([.1, .2, .3, .4])
    EA, EC, PA = np.array([2., 3., 4., 5.]), np.array([4., 2., 4., 5.]), np.array([2., 6., 2., 5.])
    controls = dict(theta=.5, delta=.1, smooth=False)
    for chi in (0., .5, 1.):
        ep, dp = (1-chi)*EA+chi*EC, (1-chi)*EA+chi*PA
        for name, local_k in [('plateau', [1., 1., 1., 1.]), ('depleted', [.01, .03, .05, .02]),
                              ('all_zero', [0., 0., 0., 0.])]:
            data = dict(sampled_k=1., k=local_k, pre_e=ep.tolist(), pre_d=dp.tolist(), metric_weights=weights.tolist())
            v, _ = oracle.capped_source(data, controls)
            B = float(np.dot(weights, dp)); w = np.minimum(1., np.asarray(local_k)/.5)
            b = float(np.dot(weights, dp*w)/B); m = max(b, .1); lam = b/m
            E, D = lam*ep, dp*w/m
            check(f'chi{chi}_{name}_same_budget_and_trace', abs(np.dot(weights, E-D)) < 2e-14 and
                  abs(np.dot(weights, E)-lam*B) < 2e-14 and abs(v['lambda']-lam) < 1e-15 and
                  relative_error(v['e'], E) < 1e-15 and relative_error(v['d'], D) < 1e-15,
                  dict(lambda_value=lam, pre_budget=B, redistribution=float(np.dot(weights, E-D))))
            records.append(dict(chi=chi, case=name, b=b, m=m, lambda_value=lam,
                pre_budget=B, dissipation_budget=float(np.dot(weights, E)),
                destruction_budget=float(np.dot(weights, D)), redistribution_budget=float(np.dot(weights, E-D)),
                actual_outer_epsilon_trace=float(E[-1]), trace_formula='lambda*Epre(H)',
                meets_lambda_guard=lam >= .5, meets_full_outer_box='requires actual physical reference normalization; lambda alone is insufficient',
                high_Re_continuation_supplied=False))
            # Compose the existing source action with changes of sampled K,
            # local k and chi-dependent pre-profiles.  The actual EPS trace
            # includes d(lambda), rather than freezing the projected value.
            direction = dict(sampled_k=.03, k=[.002, .001, .003, .002],
                pre_e=(.2*(EC-EA)).tolist(), pre_d=(.2*(PA-EA)).tolist(),
                metric_weights=[0., 0., 0., 0.])
            _, action = oracle.capped_source(data, controls, direction)
            errors = []
            for step in (1e-3, 1e-4, 1e-5, 1e-6):
                shifted = dict(sampled_k=1.+step*direction['sampled_k'])
                for key in ('k', 'pre_e', 'pre_d', 'metric_weights'):
                    shifted[key] = (np.asarray(data[key])+step*np.asarray(direction[key])).tolist()
                plus, _ = oracle.capped_source(shifted, controls)
                errors.append(dict(step=step, e_error=relative_error((plus['e']-v['e'])/step, action['e']),
                    d_error=relative_error((plus['d']-v['d'])/step, action['d']),
                    trace_error=float(abs((plus['e'][-1]-v['e'][-1])/step-action['e'][-1]))))
            check(f'chi{chi}_{name}_coupled_source_and_actual_trace_action',
                min(max(row['e_error'], row['d_error'], row['trace_error']) for row in errors) < 2e-6,
                dict(ladder=errors, one_sided=True))
            records[-1]['coupled_directional_ladder'] = errors
    return dict(cases=records, interpretation='Algebraic cancellation and dependent-trace controls, including chi=0; no claim that capped destruction repairs the original covariance or keeps outer epsilon in the declared box.')


def outer_witnesses(core, oracle, archive_source, runs, inputs, check):
    """Reevaluate fixed inputs, not a search; keep original rejected models."""
    evidence = archive_source/'docs/plans/cls-epawf-revised-evidence/m8'
    Qiv = core.Iv(core.Q_LO, core.Q_HI)
    a0 = np.zeros((3, 3), dtype=int); a0[0, 2] = 1
    ref = core.interval_array(a0+a0.T)
    records, family_counts = [], {}

    def read(path):
        value = json.loads(path.read_text()); inputs[str(path)] = sha(path); return value

    def array(value):
        return core.interval_array(np.asarray([[F(str(v)) for v in row] for row in value], dtype=object))

    def point(value):
        return float((value.lo+value.hi)/2)

    def evaluate(family, pointer, gradient, k, eps, ratio=None, direction=None, vector=None, dk=0., deps=0., original=None):
        G = array(gradient); K, E = core.Iv(k), core.Iv(eps)
        if ratio is None: ratio = K/Qiv
        else: ratio = core.Iv(ratio)
        depart = max((entry-core.Iv(int(a0[index]))).absmax() for index, entry in np.ndenumerate(G))
        inside = depart <= Decimal('.001') and ratio.lo >= Decimal('.99') and ratio.hi <= Decimal('1.01') and E.lo >= Decimal('.98') and E.hi <= Decimal('1.02')
        Gf = np.array([[point(v) for v in row] for row in G]); kf, ef = point(K), point(E)
        R = core.interval_stress(G, K, E, core.Iv(F(2, 5)), core.Iv(1), ref)
        Rf = oracle.indexed_stress(Gf, kf, ef, .4, 1., a0+a0.T).real
        eigenvalues, eigenvectors = np.linalg.eigh(Rf)
        for index in np.ndindex((3, 3)):
            # A midpoint value is only a consistency check, never an interval proof.
            check(f'{family}{pointer}_indexed_stress_{index}', float(R[index].lo)-1e-10 <= Rf[index] <= float(R[index].hi)+1e-10)
        if vector is None:
            v = eigenvectors[:, 0]
            vector_kind = 'floating minimum-eigenvector nomination at this predeclared state; interval Rayleigh evaluation only'
        else:
            v = np.array([float(F(str(x))) for x in vector])
            vector_kind = 'retained specified covariance direction'
        vi = core.interval_array(np.asarray([F(str(x)) for x in v], dtype=object))
        variance = core.inner(vi, R@vi)/core.squared_norm(vi)
        branches = {}
        for branch in ('strain', 'rotation'):
            L, q = core.interval_tangent(G, K, E, core.Iv(F(2, 5)), branch)
            # Prescribed GammaScale1, independent of the evolving state.
            for column in range(9):
                H = np.zeros((3, 3), dtype=int); H.flat[column] = 1
                L[:, column] += core.interval_array(H+H.T).ravel()
            numerical = oracle.stress_tangent(Gf, kf, ef, .4, 1., a0+a0.T, branch)
            maximum_error = 0.
            for index in np.ndindex((9, 9)):
                maximum_error = max(maximum_error, max(float(L[index].lo)-numerical[index], numerical[index]-float(L[index].hi), 0.))
            check(f'{family}{pointer}_{branch}_full_gradient_tangent', maximum_error < 2e-9,
                  dict(outside_interval_to_roundoff=maximum_error))
            output = dict(q_interval=q.pair(), full_gradient_derivative_crosscheck_error=maximum_error)
            if direction is not None:
                h = array(direction).ravel()
                form = core.inner(h, L@h)/core.squared_norm(h)
                hf = np.array([[point(v) for v in row] for row in array(direction)])
                numeric = float(hf.ravel()@numerical@hf.ravel()/np.sum(hf*hf))
                output.update(retained_direction_principal_interval=form.pair(),
                              independent_indexed_value=numeric,
                              certified_negative=form.hi < 0)
            # State-field derivatives include C' through both k/epsilon and all
            # polynomial factors.  Verify them with a central ladder; this is a
            # consistency calculation, not a coupled coercivity assertion.
            dG = np.array([[.02, -.01, .03], [.01, -.01, .02], [-.03, .02, -.01]])
            dkk = float(dk) if dk else .01*kf
            dee = float(deps) if deps else -.02*ef
            cs = 1e-30
            derivative = oracle.indexed_stress(Gf+1j*cs*dG, kf+1j*cs*dkk,
                ef+1j*cs*dee, .4, 1., a0+a0.T, branch).imag/cs
            analytic = outer_full_action(Gf, kf, ef, dG, dkk, dee, branch)
            analytic_error = relative_error(analytic, derivative)
            check(f'{family}{pointer}_{branch}_independent_full_state_action', analytic_error < 2e-10,
                  dict(relative_error=analytic_error))
            errors = []
            for step in (1e-4, 1e-5, 1e-6):
                plus = oracle.indexed_stress(Gf+step*dG, kf+step*dkk, ef+step*dee, .4, 1., a0+a0.T, branch).real
                minus = oracle.indexed_stress(Gf-step*dG, kf-step*dkk, ef-step*dee, .4, 1., a0+a0.T, branch).real
                errors.append(relative_error((plus-minus)/(2*step), derivative))
            check(f'{family}{pointer}_{branch}_coupled_k_epsilon_derivatives', min(errors) < 2e-6,
                  dict(relative_errors=errors))
            output['full_k_epsilon_directional_ladder_errors'] = errors
            output['independent_full_state_action_error'] = analytic_error
            branches[branch] = output
        record = dict(family=family, source_pointer=pointer, state=dict(gradient=gradient,
            k_interval=K.pair(), k_over_K0=ratio.pair(), epsilon_over_epsilon0=E.pair()),
            maximum_gradient_entry_departure=str(depart), inside_restricted_region=bool(inside),
            beta='2/5', GammaScale='1', GammaUniform='0', branches=branches,
            covariance=dict(direction=v.tolist(), direction_kind=vector_kind,
                normalized_variance_interval=variance.pair(), certified_negative=variance.hi < 0,
                finite_eigenvalues=eigenvalues.tolist(), continuum_certificate_at_point=core.interval_ldl((R+R.T)*F(1, 2)) is not None),
            original_model_record=original,
            interpretation='An out-of-region failure is retained as such; it does not refute or broaden the separately fixed local contract.')
        records.append(record); family_counts[family] = family_counts.get(family, 0)+1

    m1path = runs/'m1-attempt-001/verification/constitutive.json'
    m1 = read(m1path)
    cases = m1.get('cases', m1.get('constitutive_cases', []))
    if isinstance(cases, dict): cases = list(cases.values())
    for index, item in enumerate(cases):
        G = item.get('gradient', item.get('base', {}).get('gradient'))
        if G is not None:
            evaluate('M1_constitutive', f'/cases/{index}/{item.get("name", "")}', G,
                str(item['k']), str(item['epsilon']), direction=item.get('direction'), original=dict(name=item.get('name'), branch=item.get('branch')))
    native = read(runs/'m8/kernel-independent-001/kernel.json')
    for index, item in enumerate(native['closure_cases']):
        evaluate('M8_native_kernel', f'/closure_cases/{index}/{item["name"]}', item['gradient'], str(item['k']), str(item['epsilon']),
            direction=item.get('direction'), dk=item.get('direction_k', 0.), deps=item.get('direction_epsilon', 0.),
            original=dict(beta=item['beta'], gamma=item['gamma'], branch=item['branch']))
    for filename, family in [('proof-fixed-gamma-feasibility.json', 'M8_fixed_gamma_exact'),
                             ('proof-invariant-screen.json', 'M8_invariant_exact')]:
        document = read(evidence/filename)
        for ri, result in enumerate(document['results']):
            for wi, witness in enumerate(result['witnesses']):
                if not witness.get('usable', True): continue
                state = witness['state']; kr = F(state['k_over_k0'])
                evaluate(family, f'/results/{ri}/witnesses/{wi}', state['gradient'], Qiv*core.Iv(kr),
                    F(state['epsilon_over_epsilon0']), ratio=kr, direction=witness.get('rank_one'), vector=witness.get('v'),
                    original=dict(candidate=result.get('candidate'), beta=result.get('beta', result.get('beta0')), kind=witness['kind'], original_certified=witness.get('counterexample_certified', witness.get('adverse_principal_certified'))))
    controls = read(evidence/'proof-invariant-controls.json')
    for index, control in enumerate(controls['controls']):
        for kind in ('transverse', 'covariance'):
            w = control[kind]
            direction = None
            if kind == 'transverse':
                aa, xx = [F(v) for v in w['polarization']], [F(v) for v in w['wave_vector']]
                direction = [[str(a*x) for x in xx] for a in aa]
            evaluate('M8_invariant_analytic', f'/controls/{index}/{kind}', w['gradient'], F(w['k']), F(w['epsilon']),
                direction=direction, vector=w.get('direction'), original=dict(candidate=control['candidate'], beta0=control['beta0'], kind=kind))
    local = read(evidence/'proof-local-joint-11aa.json')
    for index, item in enumerate(local['outside_envelope_controls']):
        evaluate('M8_outside_envelope', f'/outside_envelope_controls/{index}',
            [[v[0] for v in row] for row in item['normalized_gradient']], F(item['k_star'][0]), F(item['epsilon_star'][0]),
            vector=['1', '0', '-1'], original=dict(beta=item['beta'], gamma=item['gamma_over_nut0'], certified_negative=item['certified_negative_variance']))
    for family in ('M1_constitutive', 'M8_native_kernel', 'M8_fixed_gamma_exact', 'M8_invariant_exact', 'M8_invariant_analytic', 'M8_outside_envelope'):
        check('retained_'+family+'_loaded', family_counts.get(family, 0) > 0, dict(count=family_counts.get(family, 0)))
    return dict(family_counts=family_counts, cases=records,
                finite_negative_covariance_count=sum(r['covariance']['certified_negative'] for r in records),
                outside_restricted_region_count=sum(not r['inside_restricted_region'] for r in records),
                original_broad_box='FAILED; necessary fixed-Gamma bounds remain contradictory; this script changes neither the original domain nor those receipts',
                scope='The local all-state interval proof, not this finite witness suite, establishes any restricted continuum margin. Historical candidate failures remain archived and are not erased by reevaluation under a different restricted law.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive-source', required=True, type=Path)
    parser.add_argument('--runs', required=True, type=Path)
    parser.add_argument('--ledger', required=True, type=Path)
    parser.add_argument('--plan', required=True, type=Path)
    parser.add_argument('--contract', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--cpu', default=50, type=int)
    args = parser.parse_args()
    if not sys.platform.startswith('linux'):
        raise SystemExit('Numerical evaluation is authorized only on Linux')
    if set(os.sched_getaffinity(0)) != {args.cpu} or args.cpu != 50:
        raise SystemExit('Run this one-worker suite with taskset -c 50')
    args.output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    plan = json.loads(args.plan.read_text())
    if plan['baseline'] != BASELINE or plan['outer']['ratio_partitions'] != 64:
        raise SystemExit('The preregistered baseline/partition count changed')
    source = args.archive_source.resolve()
    if git(source, 'rev-parse', 'HEAD') != BASELINE or git(source, 'status', '--porcelain', '--untracked-files=no'):
        raise SystemExit('Archive imports must remain clean immutable ab5')
    module_dir = source/'solvers/IncNavierStokesSolver/Tests/RANS'
    core_path, oracle_path = module_dir/'cls_epawf_m8_certificate.py', module_dir/'cls_epawf_m8_oracle.py'
    core = load_module('m8_r1_archived_certificate', core_path)
    oracle = load_module('m8_r1_archived_indexed_oracle', oracle_path)
    inputs = {str(p.resolve()): sha(p) for p in (args.ledger, args.plan, args.contract, core_path, oracle_path)}
    ledger = json.loads(args.ledger.read_text())
    checks, failures = [], []

    def check(name, condition, detail=None):
        item = dict(name=name, passed=bool(condition))
        if detail is not None: item['detail'] = detail
        checks.append(item)
        if not condition: failures.append(item)

    report = dict(schema=1, phase='M8-R1 D', started_utc=stamp(),
        baseline=BASELINE, archive_source=str(source), archive_source_commit=git(source, 'rev-parse', 'HEAD'),
        archive_source_clean_start=True, command=sys.argv, driver_sha256=sha(__file__),
        production_adopted=False, solver_or_build_execution=False, allocated_workers=1,
        cpu_affinity=sorted(os.sched_getaffinity(0)), thread_environment={key: os.environ[key] for key in THREAD_KEYS},
        versions=dict(python=platform.python_version(), numpy=np.__version__, platform=platform.platform()),
        plan=plan, input_sha256=inputs, checks=checks, failures=failures)
    print('M8-R1 kernel: inner differential/constitutive suite', flush=True)
    report['inner'] = inner_suite(core, oracle, ledger, check, plan['inner']['difference_steps'])
    report['source_coupling'] = source_suite(oracle, check)
    save(args.output/'partial-inner.json', report)
    print('M8-R1 kernel: exactly one fixed 64-partition restricted outer certificate', flush=True)
    outer = core.finite_box_certificate(F(2, 5), '.001', '.01', '.02', target='.1',
                                      ratio_splits=64, covariance_gammas=('1',))
    save(args.output/'restricted-outer-interval-proof.json', outer)
    report['outer_proof'] = dict(file='restricted-outer-interval-proof.json',
        sha256=sha(args.output/'restricted-outer-interval-proof.json'),
        certified_principal_lower_over_nut0=outer['fixed_gamma_candidates'][0]['certified_principal_lower_over_nut0'],
        covariance_certified_positive=outer['fixed_gamma_candidates'][0]['certified_positive'],
        joint_constitutive_gate=outer['fixed_gamma_candidates'][0]['joint_constitutive_gate'],
        target='.1', proof_failure_is_inconclusive=True,
        domain='entry radius .001; k/K0 [.99,1.01]; epsilon/epsilon0 [.98,1.02]; both max branches; all incompressible rank-one directions',
        original_broad_domain_status='FAILED, unchanged')
    print('M8-R1 kernel: one fixed retained-witness and derivative suite', flush=True)
    witnesses = outer_witnesses(core, oracle, source, args.runs, inputs, check)
    save(args.output/'retained-witnesses.json', witnesses)
    report['retained_witnesses'] = dict(file='retained-witnesses.json', sha256=sha(args.output/'retained-witnesses.json'),
        family_counts=witnesses['family_counts'], finite_negative_covariance_count=witnesses['finite_negative_covariance_count'],
        outside_restricted_region_count=witnesses['outside_restricted_region_count'])
    membership = ledger['states'][0]['constitutive_range_membership']
    initial = membership.get('0.001', next(iter(membership.values())))
    report['archived_initial_coverage'] = dict(source_pointer='/states/0/constitutive_range_membership/0.001',
        actual=initial, proposed_gradient_radius='.001',
        gradient_covered=initial['maximum_gradient_departure'] <= .001,
        interpretation='Original t0 gradient is outside the proposed region. Improved representation coverage must be checked independently; no state-set invariance is claimed.')
    report['archive_source_commit_end'] = git(source, 'rev-parse', 'HEAD')
    report['archive_source_clean_end'] = not git(source, 'status', '--porcelain', '--untracked-files=no')
    check('archived_source_unchanged', report['archive_source_commit_end'] == BASELINE and report['archive_source_clean_end'])
    report['verification_success'] = not failures
    report['restricted_constitutive_feasible'] = not failures and report['outer_proof']['joint_constitutive_gate']
    report['decision'] = ('One restricted research contract has constitutive feasibility; actual representation coverage, nonlinear inner reconstruction and all full discrete/forward-invariance/terrain gates remain open.'
                          if report['restricted_constitutive_feasible'] else
                          'The one restricted research proposal is unqualified by this bounded screen; preserve failed/inconclusive proof, with no parameter/box adaptation or production adoption.')
    report['lambda_guard'] = dict(lower='1/2', actual_trace_and_outer_box_also_required=True,
        lambda_alone_sufficient=False, lambda_to_zero_high_Re_continuation='unresolved outside applicability; no floor or clipping')
    report['finished_utc'], report['seconds'] = stamp(), time.monotonic()-started
    save(args.output/'validation.json', report)
    print(json.dumps(dict(verification_success=report['verification_success'],
        restricted_constitutive_feasible=report['restricted_constitutive_feasible'],
        checks=len(checks), failures=len(failures), seconds=report['seconds'],
        output=str(args.output)), sort_keys=True), flush=True)
    return 0 if report['verification_success'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
