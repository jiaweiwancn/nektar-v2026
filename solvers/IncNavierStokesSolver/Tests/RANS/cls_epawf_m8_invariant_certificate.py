#!/usr/bin/env python3
"""Directed interval screen of two optional M8 invariant-response variants.

A: C1 viscosity tail above qc=3.3 with a fixed common cubic factor beta0.
B: the same tail and beta(q)=beta0*[qc*C(qc)/(q*Ctail(q))]^2 above qc.
Neither option is adopted in production. Complete C' and beta' contributions
are included in both formal invariant branches, with physical branch hulls.
Positive one-mode factors are not promoted to all-direction/PDE qualification.

Only this new research module is mutable; the prior fixed-gamma certificate
module supplies its previously verified interval and polynomial primitives.
Run numerical work on the authoritative Linux validation environment.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import platform
import sys

sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np
import cls_epawf_m8_certificate as core

Iv, Q_LO, Q_HI = core.Iv, core.Q_LO, core.Q_HI
interval_array = core.interval_array


def interval_log(value):
    value = Iv(value)
    if value.lo <= 0:
        raise ValueError('Logarithm requires a strictly positive interval')
    with localcontext() as context:
        context.prec = core.PREC
        return Iv(value.lo.ln().next_minus(), value.hi.ln().next_plus())


def interval_power(value, exponent):
    return (Iv(exponent)*interval_log(value)).exp()


def hull(values):
    values = list(values)
    return Iv(min(value.lo for value in values), max(value.hi for value in values))


def array_hull(values):
    values = list(values)
    result = np.empty(values[0].shape, dtype=object)
    for index in np.ndindex(result.shape):
        result[index] = hull(value[index] for value in values)
    return result


@lru_cache(maxsize=1)
def tail_parameters():
    qc = Iv(F(33, 10))
    cc = core.cmu_interval(qc)
    cpc = core.cmu_derivative_interval(qc)
    alpha = -qc*cpc/cc
    p = 1-alpha
    if not (alpha.lo > 0 and alpha.hi < 1 and p.lo > 0):
        raise RuntimeError('The selected C1 tail has no certified positive shear factor')
    return dict(qc=qc, Cc=cc, Cprimec=cpc, alpha=alpha, p=p)


def coefficient_intervals(q, beta0, candidate='A', tail_branch='physical'):
    """Return C,Cprime,beta,betaprime,W=beta*C^3 and full Wprime.

    Explicit original/tail choices are formal branch evaluations. The physical
    choice splits any interval crossing qc and encloses both limiting actions.
    The tail and alpha are fixed declared functions, not iteration-frozen data.
    """
    if candidate not in ('A', 'B'):
        raise ValueError('Candidate must be A or B')
    if tail_branch not in ('physical', 'original', 'tail'):
        raise ValueError('Unknown viscosity-tail branch')
    q, beta0 = Iv(q), Iv(F(str(beta0)))
    if q.lo < 0 or beta0.lo < 0:
        raise ValueError('Nonnegative invariant and cubic factor required')
    parameters = tail_parameters()
    qc = parameters['qc']
    if tail_branch == 'physical':
        parts = []
        if q.lo <= qc.hi:
            parts.append(coefficient_intervals(Iv(q.lo, min(q.hi, qc.hi)), beta0.lo,
                                               candidate, 'original'))
        if q.hi >= qc.lo:
            parts.append(coefficient_intervals(Iv(max(q.lo, qc.lo), q.hi), beta0.lo,
                                               candidate, 'tail'))
        return {key: hull(part[key] for part in parts) for key in parts[0]}
    if tail_branch == 'original':
        c, cp = core.coefficient_bounds_interval(q)
        beta, bp = beta0, Iv(0)
    else:
        if q.lo <= 0:
            raise ValueError('Positive invariant required for formal tail')
        c = parameters['Cc']*interval_power(q/qc, -parameters['alpha'])
        cp = -parameters['alpha']*c/q
        beta = beta0 if candidate == 'A' else beta0*interval_power(q/qc, -2*parameters['p'])
        bp = Iv(0) if candidate == 'A' else -2*parameters['p']*beta/q
    return dict(C=c, Cprime=cp, beta=beta, betaprime=bp,
                W=beta*c*c*c, Wprime=bp*c*c*c+3*beta*c*c*cp)


def invariants_interval(gradient, k, epsilon):
    gradient = interval_array(np.asarray(gradient, dtype=object))
    k, epsilon = Iv(k), Iv(epsilon)
    if k.lo <= 0 or epsilon.lo <= 0:
        raise ValueError('Positive k and epsilon required')
    t = k/epsilon
    s, o = (gradient+gradient.T)*t, (gradient-gradient.T)*t
    return (core.squared_norm(s)/2).sqrt(), (core.squared_norm(o)/2).sqrt()


def _polynomial_values(gradient, k, epsilon):
    gradient = interval_array(np.asarray(gradient, dtype=object))
    k, epsilon = Iv(k), Iv(epsilon)
    if k.lo <= 0 or epsilon.lo <= 0:
        raise ValueError('Positive k and epsilon required')
    dimensionless = gradient*(k/epsilon)
    center = (dimensionless[0, 2].lo+dimensionless[0, 2].hi)/2
    flat = list(dimensionless.flat)
    flat[2] -= Iv(center)
    def evaluate(poly):
        return poly.shifted(2, center).evaluate(flat)
    def evaluate_array(array):
        return np.array([evaluate(poly) for poly in array.flat], dtype=object).reshape(array.shape)
    polynomials = core.polynomial_tensors()
    s, o, _, _, quadratic, cubic = polynomials
    return gradient, k, epsilon, polynomials, evaluate, evaluate_array, evaluate_array(s), evaluate_array(o), evaluate_array(quadratic), evaluate_array(cubic)


def interval_stress(gradient, k, epsilon, beta0, candidate='A',
                    branch='physical', tail_branch='physical'):
    """Complete physical Reynolds tensor, including the isotropic 2k/3 term."""
    _, k, _, _, _, _, s, o, quadratic, cubic = _polynomial_values(gradient, k, epsilon)
    qs, qo = (core.squared_norm(s)/2).sqrt(), (core.squared_norm(o)/2).sqrt()
    if branch not in ('physical', 'strain', 'rotation'):
        raise ValueError('Unknown max-invariant branch')
    q = qs if branch == 'strain' else qo if branch == 'rotation' else Iv(max(qs.lo, qo.lo), max(qs.hi, qo.hi))
    coefficients = coefficient_intervals(q, beta0, candidate, tail_branch)
    eye = interval_array(np.eye(3, dtype=int))
    return (eye*F(2, 3)-s*coefficients['C']+quadratic*coefficients['C']+cubic*coefficients['W'])*k


def interval_tangent(gradient, k, epsilon, beta0, candidate='A',
                     branch='physical', tail_branch='physical'):
    """Return (d(-R)/dA as a 9x9 Iv array, invariant Iv).

    k and epsilon are frozen for this principal-symbol derivative. Every
    derivative of C and beta through the selected invariant remains active.
    A physical ambiguous max branch is enclosed by the formal-branch hull.
    """
    gradient, k, epsilon, polynomials, evaluate, evaluate_array, s, o, quadratic, cubic = _polynomial_values(gradient, k, epsilon)
    qs, qo = (core.squared_norm(s)/2).sqrt(), (core.squared_norm(o)/2).sqrt()
    if branch == 'physical':
        if qs.lo > qo.hi:
            return interval_tangent(gradient, k, epsilon, beta0, candidate, 'strain', tail_branch)
        if qo.lo > qs.hi:
            return interval_tangent(gradient, k, epsilon, beta0, candidate, 'rotation', tail_branch)
        alternatives = [interval_tangent(gradient, k, epsilon, beta0, candidate, side, tail_branch)
                        for side in ('strain', 'rotation')]
        return array_hull(value[0] for value in alternatives), hull(value[1] for value in alternatives)
    if branch not in ('strain', 'rotation'):
        raise ValueError('Unknown max-invariant branch')
    q = qs if branch == 'strain' else qo
    if q.lo <= 0:
        if not all(value.lo == 0 == value.hi for value in gradient.flat):
            raise ValueError('A nontrivial box containing zero invariant needs a separate Lipschitz bound')
        c = coefficient_intervals(Iv(0), beta0, candidate, 'original')['C']
        linear = np.empty((9, 9), dtype=object)
        for column in range(9):
            h = np.zeros((3, 3), dtype=int)
            h.flat[column] = 1
            linear[:, column] = (interval_array(h+h.T)*(k*k/epsilon)*c).ravel()
        return linear, Iv(0)
    coefficients = coefficient_intervals(q, beta0, candidate, tail_branch)
    s_poly, o_poly, ss_poly, oo_poly, q_poly, b_poly = polynomials
    norm_poly = ss_poly if branch == 'strain' else oo_poly
    linear = np.empty((9, 9), dtype=object)
    for column in range(9):
        h = np.zeros((3, 3), dtype=int)
        h.flat[column] = 1
        ds = interval_array(h+h.T)
        dq = evaluate(norm_poly.derivative(column))/(4*q)
        dquad = evaluate_array(np.array([poly.derivative(column) for poly in q_poly.flat],
                                       dtype=object).reshape(3, 3))
        dcubic = evaluate_array(np.array([poly.derivative(column) for poly in b_poly.flat],
                                        dtype=object).reshape(3, 3))
        result = (ds*coefficients['C']+s*(coefficients['Cprime']*dq)
            -dquad*coefficients['C']-quadratic*(coefficients['Cprime']*dq)
            -dcubic*coefficients['W']-cubic*(coefficients['Wprime']*dq))*(k*k/epsilon)
        linear[:, column] = result.ravel()
    return linear, q


def principal_interval(linear, a, xi):
    """An exact incompressible rational rank-one test; norm need not equal one."""
    a, xi = [F(str(value)) for value in a], [F(str(value)) for value in xi]
    if sum(x*y for x, y in zip(a, xi)) != 0:
        raise ValueError('The declared rational polarization must be orthogonal to the wave vector')
    norm = sum(x*x for x in a)*sum(x*x for x in xi)
    if norm <= 0:
        raise ValueError('Nonzero rank-one direction required')
    h = np.array([Iv(x*y) for x in a for y in xi], dtype=object)
    return core.inner(h, linear@h)/Iv(norm)


def covariance_interval(stress, direction):
    values = [F(str(value)) for value in direction]
    norm = sum(value*value for value in values)
    if norm <= 0:
        raise ValueError('Nonzero covariance test vector required')
    vector = np.array([Iv(value) for value in values], dtype=object)
    return core.inner(vector, stress@vector)/Iv(norm)


def broad_box_membership(gradient, k, epsilon):
    """Exact enclosure in the preregistered normalized radius .25, +/-50% box."""
    a0 = np.zeros((3, 3), dtype=int)
    a0[0, 2] = 1
    gradient = interval_array(np.asarray(gradient, dtype=object))
    normalized_k = Iv(k)/Iv(Q_LO, Q_HI)
    return dict(gradient_inside=all(value.lo >= Decimal(str(a0[index]))-Decimal('.25') and
                                   value.hi <= Decimal(str(a0[index]))+Decimal('.25')
                                   for index, value in np.ndenumerate(gradient)),
                k_inside=normalized_k.lo >= Decimal('.5') and normalized_k.hi <= Decimal('1.5'),
                epsilon_inside=Iv(epsilon).lo >= Decimal('.5') and Iv(epsilon).hi <= Decimal('1.5'),
                k_over_k0=normalized_k.pair(), epsilon_over_epsilon0=Iv(epsilon).pair())


def analytic_controls(beta0, candidate):
    beta0 = F(str(beta0))
    gradient = np.zeros((3, 3), dtype=object)
    gradient[0, 2] = F(1)
    transverse = dict(gradient=[['0', '0', '1'], ['0', '0', '0'], ['0', '0', '0']],
        k='4', epsilon='1/2', polarization=['0', '1', '0'], wave_vector=['1', '0', '0'],
        branch_actions={}, invariant_tie='The transverse direction has zero first invariant derivatives and preserves the exact strain/rotation tie; both limiting forms agree.')
    for branch in ('strain', 'rotation'):
        linear, q = interval_tangent(gradient, Iv(4), Iv(F(1, 2)), beta0, candidate, branch, 'tail')
        form = principal_interval(linear, (0, 1, 0), (1, 0, 0))
        transverse['branch_actions'][branch] = dict(q=q.pair(), normalized_principal_form=form.pair(),
                                                   certified_negative=form.hi < 0)
    transverse['box_membership'] = broad_box_membership(gradient, Iv(4), Iv(F(1, 2)))
    transverse['certified_physical_negative'] = all(row['certified_negative'] for row in transverse['branch_actions'].values())
    # At parallel shear the common cubic terms cancel exactly for either
    # beta(q), so this covariance control tests the new C(q) tail itself.
    gradient[0, 2] = F(5, 4)
    stress = interval_stress(gradient, Iv(F(19, 4)), Iv(F(1, 2)), beta0, candidate)
    variance = covariance_interval(stress, (1, 0, 5))
    covariance = dict(gradient=[['0', '0', '5/4'], ['0', '0', '0'], ['0', '0', '0']],
        k='19/4', epsilon='1/2', q='95/8', direction=['1', '0', '5'],
        direction_normalization='divide quadratic form by26',
        variance_interval=variance.pair(), certified_negative=variance.hi < 0,
        stress_intervals=[[value.pair() for value in row] for row in stress],
        box_membership=broad_box_membership(gradient, Iv(F(19, 4)), Iv(F(1, 2))))
    return dict(candidate=candidate, beta0=str(beta0), transverse=transverse, covariance=covariance)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--candidate', action='append', choices=('A', 'B'))
    parser.add_argument('--beta', action='append')
    args = parser.parse_args()
    if platform.system() != 'Linux':
        raise RuntimeError('Numerical research qualification runs only on authoritative Linux')
    report = dict(schema=1, kind='Optional M8 invariant-response mathematical screen',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        fixed_certificate_source_sha256=hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),
        started_utc=datetime.now(timezone.utc).isoformat(), production_adopted=False,
        tail_parameters={key: value.pair() for key, value in tail_parameters().items()},
        original_note_example_only=True, reference_certificates=[], controls=[])
    try:
        for beta in args.beta or ('2/5', '9/20'):
            reference, _ = core.reference_certificate(F(beta))
            if not Iv(Q_LO, Q_HI).hi < tail_parameters()['qc'].lo:
                raise RuntimeError('The new tail changes the neutral reference')
            report['reference_certificates'].append(reference)
            for candidate in args.candidate or ('A', 'B'):
                report['controls'].append(analytic_controls(F(beta), candidate))
        report.update(status='screen_completed',
            interpretation='Only explicit negative principal/variance intervals establish failures. Positive selected directions or inconclusive intervals never establish an admissible constitutive model.')
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(status=report['status'], output=str(args.output), production_adopted=False), indent=2))


if __name__ == '__main__':
    main()
