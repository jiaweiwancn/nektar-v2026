"""Independent indexed R2 inner tensor, for Linux numerical verification only.

No production wrapper is imported.  Complex inputs differentiate the open
branches without freezing the local-k/nuA factors.  Exact boundary and joint
origin actions require admissible one-sided differences, not a complex-step
extension outside the physical scalar cone.
"""
import os
import sys
from types import FunctionType

sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np

_OUTER_EVALUATIONS = dict(calls=0, rounded_tail_calls=0, maximum_real_q=0.,
                         maximum_omitted_tensor_Frobenius_bound=0.)


def outer_evaluation_statistics():
    return dict(_OUTER_EVALUATIONS, coefficient_tail=coefficient_tail_contract(),
        scope='Finite oracle calls in this Python process; the tensor bound excludes ordinary contraction roundoff and does not certify a state envelope')


def coefficient_tail_contract():
    """Real-axis omission bounds, not a field or invariant truncation.

    For q>8 set F=1 instead of F=1-exp(-x), x=.36 exp(.75q).
    Since exp(6)>403.4, exp(-x)<1e-63. Also x exp(-x)
    decreases for x>1, so |F'|=.75 x exp(-x)<1.1e-61.
    B=.3/(1+.35q^1.5) <=.3 and |B'|<.0072 for q>=8:
    |C-B|<3e-64 and |C'-B'|<4e-62. These bounds cover
    the real-axis derivative of the discarded nested-exponential factor;
    complex steps use the analytic continuation of B in that open branch.
    The actual q and all polynomial tensor terms remain unmodified.
    """
    return dict(threshold=8., factor_value_error_upper=1e-63,
        factor_q_derivative_error_upper=1.1e-61,
        coefficient_value_error_upper=3e-64,
        coefficient_q_derivative_error_upper=4e-62,
        exact_factor='1-exp(-.36*exp(.75*q))',
        stable_tail_factor='1; correctly rounded at binary64 precision',
        proof='exp(6)>403.4; .36*403.4>63*log(10); x*exp(-x) decreases for x>1; B<=.3 and |Bprime|<.0072 for q>=8',
        scope='Real q>=8; no q/field clipping, closure retuning, or finite-state applicability inference')


def indexed_outer(gradient, k, epsilon, branch='physical'):
    """R2-only dimensional CLS index oracle with a stable nested exponential.

    Historical M1/R1/M8 source files remain unchanged. Only the exponentially
    invisible factor is rounded on q>8; C(q), every dimensional prefactor,
    both max branches and the complete quadratic/cubic tensors are retained.
    """
    if branch not in ('physical', 'strain', 'rotation'):
        raise ValueError('Unknown CLS invariant branch')
    if np.real(k) <= 0 or np.real(epsilon) <= 0:
        raise ValueError('Outer indexed stress requires positive k and epsilon')
    g = np.asarray(gradient, complex)
    s, o, eye = g+g.T, g-g.T, np.eye(3)
    ss, oo = np.sum(s*s), np.sum(o*o)
    qs, qo = k/epsilon*np.sqrt(ss/2), k/epsilon*np.sqrt(oo/2)
    q = qs if branch == 'strain' else qo if branch == 'rotation' else (qs if qs.real >= qo.real else qo)
    if not np.isfinite(q):
        raise ValueError('Nonfinite invariant in R2 independent outer oracle')
    tail = bool(np.real(q) > 8.)
    factor = 1. if tail else -np.expm1(-.36*np.exp(.75*q))
    c = .3/(1+.35*q**1.5)*factor
    nt = c*k*k/epsilon
    t1 = np.einsum('ik,jk->ij', s, s)-ss*eye/3
    t2 = np.einsum('ik,jk->ij', o, s)+np.einsum('jk,ik->ij', o, s)
    t3 = np.einsum('ik,jk->ij', o, o)-oo*eye/3
    t4 = np.einsum('ki,lj,kl->ij', s, o, s)+np.einsum('kj,li,kl->ij', s, o, s)
    quadratic_shape = -.1*t1+.1*t2+.26*t3
    cubic_shape = -10*t4+5*(oo-ss)*s
    quadratic = k*nt/epsilon*quadratic_shape
    cubic = nt*k*k/epsilon**2*c*c*cubic_shape
    stress = 2*k*eye/3-nt*s+quadratic+cubic
    if not np.all(np.isfinite(stress)):
        raise ValueError('Nonfinite tensor in R2 independent outer oracle')
    # A conservative real-state Frobenius bound on the omitted tensor term.
    # |C^3-B^3| <=3 B^2 |C-B| and B<=.3. This is distinct from
    # floating-point roundoff in the remaining tensor contractions.
    stress_error = 0.
    if tail and np.isreal(k) and np.isreal(epsilon) and not np.any(g.imag):
        stress_error = 3e-64*(abs(k*k/epsilon)*np.linalg.norm(s)
            +abs(k**3/epsilon**2)*np.linalg.norm(quadratic_shape)
            +.27*abs(k**4/epsilon**3)*np.linalg.norm(cubic_shape))
    _OUTER_EVALUATIONS['calls'] += 1
    _OUTER_EVALUATIONS['rounded_tail_calls'] += int(tail)
    _OUTER_EVALUATIONS['maximum_real_q'] = max(_OUTER_EVALUATIONS['maximum_real_q'], float(np.real(q)))
    _OUTER_EVALUATIONS['maximum_omitted_tensor_Frobenius_bound'] = max(
        _OUTER_EVALUATIONS['maximum_omitted_tensor_Frobenius_bound'], float(stress_error))
    return dict(stress=stress, quadratic=quadratic, cubic=cubic, cMu=c,
        nut=nt, strain=qs, rotation=qo, production=-np.sum(stress*g),
        nested_tail_rounded=tail, nested_tail_stress_Frobenius_error_upper=float(stress_error))


def selected_outer(raw, gradient, k, epsilon, chart):
    """Independent original CLS plus the prescribed M8 beta/Gamma tensors."""
    from validate_cls_epawf_m5_geometry import geometry
    value = indexed_outer(gradient, k, epsilon)
    stress = value['stress'].real.copy()
    if raw['model_variant'] != 'Original':
        beta, uniform, _, _, _, scale = raw['robust_parameters'][:6]
        distance = float(geometry(raw['map_parameters'], chart)['distance'])
        flow, normal = np.asarray(raw['robust_reference_flow']), np.asarray(raw['robust_reference_normal'])
        us, kappa = raw['ustar'], raw['parameters'][7]
        reference = us/(kappa*distance)*(np.outer(flow, normal)+np.outer(normal, flow))
        stress += (beta-1)*value['cubic'].real
        stress -= (uniform+scale*kappa*us*distance)*(gradient+gradient.T-reference)
    return stress, float(value['nut'].real)


def bind_outer_oracle(function):
    """Copy one reused operator's globals; never mutate historical modules.

    The reused Jacobi/quadrature assembly stays identical. Only its explicitly
    named indexed/selected stress dependencies use the R2 stable evaluator.
    """
    bindings = dict(function.__globals__, indexed=indexed_outer, selected_stress=selected_outer)
    return FunctionType(function.__code__, bindings, function.__name__,
                        function.__defaults__, function.__closure__)


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
