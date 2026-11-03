#!/usr/bin/env python3
"""Independent dimensional indexed CLS and capped-source M8 oracles.

This module never imports the C++ production kernel and does not run on import.
The actual same stress must be contracted for production and interface traction.
"""
import os
import sys

sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np


def indexed_closure(gradient, k, epsilon, beta=1., gamma=0.,
                    reference_strain=None, branch='physical', distance=1., nu=0.):
    """Dimensional indexed stress; coefficients and production differentiated too."""
    gradient = np.asarray(gradient, dtype=complex)
    strain, rotation, eye = gradient+gradient.T, gradient-gradient.T, np.eye(3)
    ss, oo = np.sum(strain*strain), np.sum(rotation*rotation)
    qs, qo = k/epsilon*np.sqrt(ss/2), k/epsilon*np.sqrt(oo/2)
    q = qs if branch == 'strain' else qo if branch == 'rotation' else (
        qs if qs.real >= qo.real else qo)
    with np.errstate(over='ignore', under='ignore'):
        c = .3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(.75*q)))
    nt = c*k*k/epsilon
    t1 = np.einsum('ik,jk->ij', strain, strain)-ss*eye/3
    t2 = np.einsum('ik,jk->ij', rotation, strain)+np.einsum('jk,ik->ij', rotation, strain)
    t3 = np.einsum('ik,jk->ij', rotation, rotation)-oo*eye/3
    # Independent indexed cubic contraction, rather than the production matrix products.
    t4 = (np.einsum('ki,lj,kl->ij', strain, rotation, strain)
          +np.einsum('kj,li,kl->ij', strain, rotation, strain))
    quadratic = k*nt/epsilon*(-.1*t1+.1*t2+.26*t3)
    original_cubic = nt*k*k/epsilon**2*(-10*c*c*t4+5*c*c*(oo-ss)*strain)
    imposed = np.zeros((3, 3)) if reference_strain is None else np.asarray(reference_strain)
    stress = 2*k*eye/3-nt*strain+quadratic+beta*original_cubic-gamma*(strain-imposed)
    production = -np.sum(stress*gradient)
    length_ratio = k**1.5/(2.5*epsilon*distance)
    yap = epsilon*epsilon/k*.83*(length_ratio-1)*length_ratio**2 if length_ratio.real > 1 else 0.*length_ratio
    return dict(stress=stress, quadratic=quadratic, cubic=beta*original_cubic,
                production=production, k_source=production-epsilon,
                epsilon_production=1.44*epsilon/k*production,
                epsilon_destruction=1.92*epsilon*epsilon/k, yap=yap,
                epsilon_source=1.44*epsilon/k*production-1.92*epsilon*epsilon/k+yap,
                cMu=c, nut=nt, strain=qs, rotation=qo, divergence=np.trace(gradient),
                momentum_diffusion=nu+nt, k_diffusion=nu+nt,
                epsilon_diffusion=nu+nt/1.3)


def indexed_stress(gradient, k, epsilon, beta=1., gamma=0., reference_strain=None,
                   branch='physical'):
    return indexed_closure(gradient, k, epsilon, beta, gamma, reference_strain,
                           branch)['stress']


def stress_tangent(gradient, k, epsilon, beta=1., gamma=0., reference_strain=None,
                    branch='physical'):
    """Complex-step -R tangent, with fixed beta, gamma and prescribed Sref."""
    gradient = np.asarray(gradient, dtype=float)
    output = np.empty((9, 9))
    step = 1e-30
    for column in range(9):
        perturbed = gradient.astype(complex)
        perturbed.flat[column] += 1j*step
        output[:, column] = (-indexed_stress(perturbed, k, epsilon, beta, gamma,
                                            reference_strain, branch).imag/step).ravel()
    return output


def capped_source(data, controls, direction=None):
    """Independent scalar/ray arithmetic and exact one-sided directional action.

    A supplied direction selects the active one-sided branch at exact knots.
    Endpoint extensions are derivative tools, never acceptance of negative k.
    """
    direction = direction or {}
    sampled = float(data['sampled_k'])
    k, ep, destruction, weights = [np.asarray(data[key], dtype=float)
                                  for key in ('k', 'pre_e', 'pre_d', 'metric_weights')]
    dsampled = float(direction.get('sampled_k', 0.))
    derivatives = []
    for name, value in zip(('k', 'pre_e', 'pre_d', 'metric_weights'), (k, ep, destruction, weights)):
        raw = direction.get(name, [])
        derivatives.append(np.asarray(raw, dtype=float) if len(raw) else np.zeros_like(value))
    dk, dep, ddestruction, dweights = derivatives
    theta, delta = controls['theta'], controls['delta']
    ratio = k/sampled
    dratio = dk/sampled-k*dsampled/sampled**2
    limiter, dlimiter = np.zeros_like(k), np.zeros_like(k)
    for index, (value, slope) in enumerate(zip(ratio, dratio)):
        if controls['smooth']:
            if value >= theta:
                limiter[index] = 1
            elif value > 0:
                t = value/theta
                limiter[index] = 10*t**3-15*t**4+6*t**5
                dlimiter[index] = 30*t*t*(1-t)**2/theta*slope
        elif value > theta:
            limiter[index] = 1
        elif value == theta:
            limiter[index] = 1
            dlimiter[index] = min(slope, 0.)/theta
        else:
            limiter[index] = value/theta
            dlimiter[index] = slope/theta
    budget_e, budget_d = np.dot(weights, ep), np.dot(weights, destruction)
    dbudget_e, dbudget_d = np.dot(dweights, ep)+np.dot(weights, dep), np.dot(dweights, destruction)+np.dot(weights, ddestruction)
    numerator = np.dot(weights, destruction*limiter)
    dnumerator = np.dot(dweights, destruction*limiter)+np.dot(weights, ddestruction*limiter+destruction*dlimiter)
    b = numerator/budget_d
    db = dnumerator/budget_d-numerator*dbudget_d/budget_d**2
    # Inputs assembled from decimal constants may represent an intended exact
    # knot to roundoff. Use its exact limiting branch when within 8 ulps.
    at_cap = abs(b-delta) <= 8*np.finfo(float).eps*max(1., abs(b), abs(delta))
    m = max(b, delta)
    dm = max(db, 0.) if at_cap else db if b > delta else 0.
    lam, dlam = b/m, db/m-b*dm/m**2
    e, d = lam*ep, destruction*limiter/m
    de = dlam*ep+lam*dep
    dd = (ddestruction*limiter+destruction*dlimiter)/m-destruction*limiter*dm/m**2
    rate = np.zeros_like(k)
    drate = np.zeros_like(k)
    lipschitz = 15/8 if controls['smooth'] else 1.
    rate_bound = lipschitz*destruction/(delta*theta*sampled)
    drate_bound = lipschitz/(delta*theta)*(ddestruction/sampled-destruction*dsampled/sampled**2)
    positive = k > 0
    rate[positive] = d[positive]/k[positive]
    drate[positive] = dd[positive]/k[positive]-d[positive]*dk[positive]/k[positive]**2
    if not controls['smooth']:
        rate[~positive] = destruction[~positive]/(theta*sampled*m)
        drate[~positive] = (ddestruction[~positive]/(theta*sampled*m)
            -destruction[~positive]*dsampled/(theta*sampled**2*m)
            -destruction[~positive]*dm/(theta*sampled*m**2))
    values = dict(sampled_k=sampled, budget_e=budget_e, budget_d=budget_d,
                  budget_mismatch=budget_e-budget_d, b=b, m=m, **{'lambda': lam},
                  capped=b < delta, dissipation_budget=np.dot(weights, e),
                  destruction_budget=np.dot(weights, d), redistribution_budget=np.dot(weights, e-d),
                  limiter=limiter, e=e, d=d, c=e-d, rate=rate, rate_bound=rate_bound)
    action = dict(sampled_k=dsampled, budget_e=dbudget_e, budget_d=dbudget_d,
                  budget_mismatch=dbudget_e-dbudget_d, b=db, m=dm, **{'lambda': dlam},
                  capped=b < delta, dissipation_budget=np.dot(dweights, e)+np.dot(weights, de),
                  destruction_budget=np.dot(dweights, d)+np.dot(weights, dd),
                  redistribution_budget=np.dot(dweights, e-d)+np.dot(weights, de-dd),
                  limiter=dlimiter, e=de, d=dd, c=de-dd, rate=drate,
                  rate_bound=drate_bound)
    return values, action
