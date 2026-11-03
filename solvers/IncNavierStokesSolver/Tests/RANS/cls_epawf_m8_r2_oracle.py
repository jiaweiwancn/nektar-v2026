"""Independent indexed R2 inner tensor, for Linux numerical verification only.

No production wrapper is imported.  Complex inputs differentiate the open
branches without freezing the local-k/nuA factors.  Exact boundary and joint
origin actions require admissible one-sided differences, not a complex-step
extension outside the physical scalar cone.
"""
import os
import sys

sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np


def indexed_inner(gradient, k, auxiliary_nu, selector='SaturatedR2'):
    if selector not in ('Linear', 'SaturatedR2'):
        raise ValueError('Unknown inner-stress selector')
    if np.real(k) < 0 or np.real(auxiliary_nu) < 0:
        raise ValueError('Negative k/nuA is outside the inner stress domain')
    G = np.asarray(gradient, dtype=complex)
    S = np.empty((3, 3), complex)
    for i in range(3):
        for j in range(3):
            S[i, j] = G[i, j]+G[j, i]
    s = np.sqrt(sum(S[i, j]*S[i, j] for i in range(3) for j in range(3)))
    c, a = 2*k/3, auxiliary_nu*s
    active, ta, tc = False, 1., 0.
    if selector == 'SaturatedR2' and k == 0:
        t, ta, tc = 0., 0., 0.
    elif selector == 'Linear' or np.real(a) <= .8*np.real(c):
        t = a
    else:
        z = (a-.8*c)/(.2*c)
        q = np.exp(-z)
        t, ta, tc, active = c*(1-.2*q), q, 1-q*(1+.2*z), True
    F, R = np.zeros((3, 3), complex), np.zeros((3, 3), complex)
    for i in range(3):
        for j in range(3):
            F[i, j] = (t*S[i, j]/s) if s != 0 else 0.
            R[i, j] = (c if i == j else 0.)-F[i, j]
    production = -sum(R[i, j]*G[i, j] for i in range(3) for j in range(3))
    radial = auxiliary_nu*ta
    tangential = t/s if s != 0 else auxiliary_nu if k != 0 or selector == 'Linear' else 0.
    return dict(stress=R, turbulent=F, strain=S, production=production,
                s=s, a=a, c=c, t=t, ta=ta, tc=tc, radial=radial,
                tangential=tangential, active=active)


def inner_mean(gradient, k, auxiliary_nu, molecular_nu, reduced_pressure,
               selector='SaturatedR2'):
    value = indexed_inner(gradient, k, auxiliary_nu, selector)
    return molecular_nu*value['strain']+value['turbulent']-reduced_pressure*np.eye(3)


def selected_inner(raw, gradient, k, auxiliary_nu):
    return indexed_inner(gradient, k, auxiliary_nu,
                         raw.get('inner_stress_model', 'Linear'))
