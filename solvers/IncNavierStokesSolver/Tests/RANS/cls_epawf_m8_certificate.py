#!/usr/bin/env python3
"""Directed-rounding certificates for the explicitly selected M8 CLS variants.

The reference null-form certificate follows the earlier note's embedded listing.
The finite-box calculation differentiates the full polynomial at interval
states, rather than treating a sampled direction/state scan as a certificate.
The analytic fallback also covers zero strain/rotation. These are frozen-k,
epsilon incompressible momentum principal-part statements, never positivity,
covariance, nonlinear CFD, or terrain-validation statements.

Run numerical work on the project's authoritative Linux validation host.
"""
import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
for _key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
             'BLIS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_key] = '1'
import numpy as np

PREC = 60
F = Fraction
Q_LO, Q_HI = '3.19067508162428', '3.19067508162430'
REFERENCE_WEIGHTS = {
    F(1, 5): (F(3, 20), {'strain': [145, 0, -3, 0, 128, 0, -290, 0, 299],
                        'rotation': [106, 0, 9, 0, 273, 0, -158, 0, 183]}),
    F(2, 5): (F(3, 25), {'strain': [137, 0, 4, 0, 212, 0, -270, 1, 366],
                        'rotation': [90, 0, 51, 0, 321, 0, -197, 0, 206]}),
    F(9, 20): (F(1, 10), {'strain': [112, 0, 44, -7, 183, -1, -230, 0, 335],
                         'rotation': [112, 0, 40, 0, 300, 1, -230, 0, 235]})}


def rd(operation, rounding):
    with localcontext() as context:
        context.prec, context.rounding = PREC, rounding
        return operation()


class Iv:
    """Closed Decimal interval; arithmetic endpoints are directed outward."""
    def __init__(self, lo, hi=None):
        if isinstance(lo, Iv):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            def endpoint(value, rounding):
                if isinstance(value, Fraction):
                    return rd(lambda: Decimal(value.numerator)/Decimal(value.denominator), rounding)
                return Decimal(str(value))
            self.lo = endpoint(lo, ROUND_FLOOR)
            self.hi = endpoint(lo if hi is None else hi, ROUND_CEILING)
        if self.lo > self.hi:
            raise ValueError('Reversed interval')

    def __add__(self, value):
        if isinstance(value, np.ndarray):
            return NotImplemented
        value = Iv(value)
        return Iv(rd(lambda: self.lo+value.lo, ROUND_FLOOR),
                  rd(lambda: self.hi+value.hi, ROUND_CEILING))
    __radd__ = __add__

    def __neg__(self):
        return Iv(self.hi.copy_negate(), self.lo.copy_negate())

    def __sub__(self, value):
        if isinstance(value, np.ndarray):
            return NotImplemented
        return self+-Iv(value)

    def __rsub__(self, value):
        return Iv(value)+-self

    def __mul__(self, value):
        if isinstance(value, np.ndarray):
            return NotImplemented
        value = Iv(value)
        lows = [rd(lambda: a*b, ROUND_FLOOR)
                for a in (self.lo, self.hi) for b in (value.lo, value.hi)]
        highs = [rd(lambda: a*b, ROUND_CEILING)
                 for a in (self.lo, self.hi) for b in (value.lo, value.hi)]
        return Iv(min(lows), max(highs))
    __rmul__ = __mul__

    def __truediv__(self, value):
        if isinstance(value, np.ndarray):
            return NotImplemented
        value = Iv(value)
        if value.lo <= 0 <= value.hi:
            raise ValueError('Interval division through zero')
        return self*Iv(rd(lambda: Decimal(1)/value.hi, ROUND_FLOOR),
                       rd(lambda: Decimal(1)/value.lo, ROUND_CEILING))

    def __rtruediv__(self, value):
        return Iv(value)/self

    def square(self):
        upper = rd(lambda: max(self.lo*self.lo, self.hi*self.hi), ROUND_CEILING)
        lower = Decimal(0) if self.lo <= 0 <= self.hi else rd(
            lambda: min(self.lo*self.lo, self.hi*self.hi), ROUND_FLOOR)
        return Iv(lower, upper)

    def exp(self):
        # Decimal exp/sqrt are correctly rounded to nearest, not directed.
        # Expanding each endpoint one representable value encloses the result.
        with localcontext() as context:
            context.prec = PREC
            return Iv(max(Decimal(0), self.lo.exp().next_minus()),
                      self.hi.exp().next_plus())

    def sqrt(self):
        if self.lo < 0:
            raise ValueError('Square root of negative interval')
        with localcontext() as context:
            context.prec = PREC
            return Iv(max(Decimal(0), self.lo.sqrt().next_minus()),
                      self.hi.sqrt().next_plus())

    def absmax(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())

    def pair(self):
        return [str(self.lo), str(self.hi)]


def cmu_interval(q):
    q = Iv(q)
    a = Iv(F(7, 20))*q*q.sqrt()
    u = Iv(F(9, 25))*(Iv(F(3, 4))*q).exp()
    return Iv(F(3, 10))/(1+a)*(1-(-u).exp())


def cmu_derivative_interval(q):
    """Direct C' formula, finite also at q=0 (no log-slope division)."""
    q = Iv(q)
    a, da = Iv(F(7, 20))*q*q.sqrt(), Iv(F(21, 40))*q.sqrt()
    u = Iv(F(9, 25))*(Iv(F(3, 4))*q).exp()
    e = (-u).exp()
    return -Iv(F(3, 10))*da/(1+a).square()*(1-e)+Iv(F(9, 40))*u*e/(1+a)


def interval_array(values):
    result = np.empty(np.asarray(values, dtype=object).shape, dtype=object)
    for index in np.ndindex(result.shape):
        result[index] = Iv(values[index])
    return result


def squared_norm(array):
    return sum((entry.square() for entry in array.flat), Iv(0))


def inner(left, right):
    return sum((a*b for a, b in zip(left.flat, right.flat)), Iv(0))


class Polynomial:
    """Exact rational multivariate polynomial, combining cancellations first."""
    def __init__(self, value=0, variable=None):
        if isinstance(value, Polynomial):
            self.terms = value.terms.copy()
        elif variable is not None:
            power = [0]*9
            power[variable] = 1
            self.terms = {tuple(power): F(1)}
        else:
            self.terms = {} if value == 0 else {(0,)*9: F(value)}

    def __add__(self, value):
        if isinstance(value, np.ndarray):
            return NotImplemented
        value, result = Polynomial(value), Polynomial(self)
        for power, coefficient in value.terms.items():
            result.terms[power] = result.terms.get(power, F(0))+coefficient
            if result.terms[power] == 0:
                del result.terms[power]
        return result
    __radd__ = __add__

    def __neg__(self):
        result = Polynomial()
        result.terms = {power: -coefficient for power, coefficient in self.terms.items()}
        return result

    def __sub__(self, value):
        return self+-Polynomial(value)

    def __rsub__(self, value):
        return Polynomial(value)+-self

    def __mul__(self, value):
        if isinstance(value, np.ndarray):
            return NotImplemented
        value, result = Polynomial(value), Polynomial()
        for left, a in self.terms.items():
            for right, b in value.terms.items():
                power = tuple(x+y for x, y in zip(left, right))
                result.terms[power] = result.terms.get(power, F(0))+a*b
        result.terms = {power: coefficient for power, coefficient in result.terms.items() if coefficient}
        return result
    __rmul__ = __mul__

    def derivative(self, variable):
        result = Polynomial()
        for power, coefficient in self.terms.items():
            if power[variable]:
                child = list(power)
                child[variable] -= 1
                result.terms[tuple(child)] = coefficient*power[variable]
        return result

    def shifted(self, variable, center):
        """Exact univariate translation; used on the nonzero shear coordinate."""
        result = Polynomial()
        center = F(center)
        for power, coefficient in self.terms.items():
            order = power[variable]
            for child_order in range(order+1):
                child = list(power)
                child[variable] = child_order
                child = tuple(child)
                result.terms[child] = result.terms.get(child, F(0))+coefficient*math.comb(order, child_order)*center**(order-child_order)
        result.terms = {power: coefficient for power, coefficient in result.terms.items() if coefficient}
        return result

    def evaluate(self, values):
        result = Iv(0)
        for power, coefficient in self.terms.items():
            term = Iv(coefficient)
            for variable, order in enumerate(power):
                if order == 1:
                    term *= values[variable]
                elif order == 2:
                    term *= values[variable].square()
                elif order == 3:
                    term *= values[variable].square()*values[variable]
                elif order:
                    raise ValueError('Unexpected polynomial degree')
            result += term
        return result


_POLYNOMIAL_TENSORS = None


def polynomial_tensors():
    global _POLYNOMIAL_TENSORS
    if _POLYNOMIAL_TENSORS is not None:
        return _POLYNOMIAL_TENSORS
    a = np.array([Polynomial(variable=i) for i in range(9)], dtype=object).reshape(3, 3)
    s, o = a+a.T, a-a.T
    ss, oo = sum((v*v for v in s.flat), Polynomial()), sum((v*v for v in o.flat), Polynomial())
    identity = np.eye(3, dtype=int)
    q = (s@s-identity*ss*F(1, 3))*F(-1, 10)+(o@s-s@o)*F(1, 10)+(o@o.T-identity*oo*F(1, 3))*F(13, 50)
    b = (s@s@o-o@s@s)*(-10)+s*((oo-ss)*5)
    _POLYNOMIAL_TENSORS = s, o, ss, oo, q, b
    return _POLYNOMIAL_TENSORS


def expanded_interval_tangent(gradient, ratio, T, beta, branch):
    """Same exact tangent, with polynomial cancellation before interval arithmetic.

    A dimensionless B=(k/epsilon)A is formed; tangent = (k^2/epsilon)*L(B).
    This grouping keeps scalar ratio subdivision inexpensive and avoids
    independently bounding canceling matrix products in the cubic polynomial.
    """
    s_poly, o_poly, ss_poly, oo_poly, q_poly, b_poly = polynomial_tensors()
    b_gradient = gradient*ratio
    center = (b_gradient[0, 2].lo+b_gradient[0, 2].hi)/2
    flat = list(b_gradient.flat)
    flat[2] = flat[2]-Iv(center)
    def evaluate(poly):
        return poly.shifted(2, center).evaluate(flat)
    def evaluate_array(array):
        return np.array([evaluate(poly) for poly in array.flat], dtype=object).reshape(array.shape)
    s, o = evaluate_array(s_poly), evaluate_array(o_poly)
    # A norm is evaluated as a sum of nonnegative squares instead of a
    # monomial expansion whose independent intervals can become negative.
    invariant = (squared_norm(s if branch == 'strain' else o)/2).sqrt()
    if invariant.lo <= 0:
        raise ValueError('Finite-box branch includes zero; use analytic fallback')
    c, cp = cmu_interval(invariant), cmu_derivative_interval(invariant)
    q_value, b_value = evaluate_array(q_poly), evaluate_array(b_poly)
    norm_poly = ss_poly if branch == 'strain' else oo_poly
    linear = np.empty((9, 9), dtype=object)
    for column in range(9):
        ds = np.array([Iv(int((i == column//3 and j == column%3)+(i == column%3 and j == column//3)))
                       for i in range(3) for j in range(3)], dtype=object).reshape(3, 3)
        dq = evaluate(norm_poly.derivative(column))/(4*invariant)
        dquad = np.array([evaluate(poly.derivative(column)) for poly in q_poly.flat], dtype=object).reshape(3, 3)
        dcubic = np.array([evaluate(poly.derivative(column)) for poly in b_poly.flat], dtype=object).reshape(3, 3)
        result = (ds*c+s*(cp*dq)-dquad*c-q_value*(cp*dq)
                  -(dcubic*(c*c*c)+b_value*(3*c*c*cp*dq))*beta)*T
        linear[:, column] = result.ravel()
    return linear, invariant


def null_forms():
    forms, names = [], []
    for i in range(3):
        for j in range(i+1, 3):
            for r in range(3):
                for s in range(r+1, 3):
                    matrix = np.full((9, 9), F(0), dtype=object)
                    matrix[3*i+r, 3*j+s] += F(1, 2)
                    matrix[3*j+s, 3*i+r] += F(1, 2)
                    matrix[3*i+s, 3*j+r] -= F(1, 2)
                    matrix[3*j+r, 3*i+s] -= F(1, 2)
                    forms.append(matrix)
                    names.append([i, j, r, s])
    basis = np.zeros((9, 8), dtype=int)
    for column, row in enumerate([1, 2, 3, 5, 6, 7]):
        basis[row, column] = 1
    basis[0, 6], basis[8, 6], basis[4, 7], basis[8, 7] = 1, -1, 1, -1
    return forms, names, basis


def interval_ldl(matrix):
    """Return directed pivots iff this complete interval family is SPD."""
    size = len(matrix)
    lower = np.full((size, size), Iv(0), dtype=object)
    diagonal = []
    for i in range(size):
        lower[i, i] = Iv(1)
        value = matrix[i, i]-sum((lower[i, k].square()*diagonal[k]
                                  for k in range(i)), Iv(0))
        if value.lo <= 0:
            return None
        diagonal.append(value)
        for j in range(i+1, size):
            lower[j, i] = (matrix[j, i]-sum((lower[j, k]*lower[i, k]*diagonal[k]
                                           for k in range(i)), Iv(0)))/value
    return diagonal


def reference_certificate(beta=F(2, 5)):
    """Reproduce both note branches using rational null forms and interval LDL."""
    beta = F(beta)
    margin, weights = REFERENCE_WEIGHTS[beta]
    q = Iv(Q_LO, Q_HI)
    flo = cmu_interval(Iv(q.lo))*Iv(q.lo).square()-1
    fhi = cmu_interval(Iv(q.hi))*Iv(q.hi).square()-1
    if not flo.hi < 0 < fhi.lo:
        raise RuntimeError('CLS root is not bracketed')
    a = Iv(F(7, 20))*q*q.sqrt()
    u = Iv(F(9, 25))*(Iv(F(3, 4))*q).exp()
    e = (-u).exp()
    slope = -Iv(F(3, 2))*a/(1+a)+Iv(F(3, 4))*q*u*e/(1-e)
    if (2+slope).lo <= 0:
        raise RuntimeError('Reference root derivative is not positive')
    s, o, eye = [np.zeros((3, 3), dtype=int) for _ in range(3)]
    s[0, 2] = s[2, 0] = o[0, 2] = 1
    o[2, 0] = -1
    eye[:] = np.eye(3, dtype=int)
    quadratic = (F(-1, 10)*(s@s-F(2, 3)*eye)+F(1, 10)*(o@s-s@o)
                 +F(13, 50)*(o@o.T-F(2, 3)*eye))
    forms, names, basis = null_forms()
    report = dict(beta=str(beta), q_bracket=q.pair(), Cq2_minus_1_at_lower=flo.pair(),
                  Cq2_minus_1_at_upper=fhi.pair(), log_slope=slope.pair(),
                  coercivity_margin_over_nut=str(margin), minor_order=names,
                  scope='neutral-reference frozen-k/epsilon incompressible rank-one momentum principal form',
                  branches={})
    tangent = {}
    for branch, coefficients in weights.items():
        linear = np.empty((9, 9), dtype=object)
        for column in range(9):
            h = np.zeros((3, 3), dtype=int)
            h.flat[column] = 1
            ds, do = h+h.T, h-h.T
            dq = F(-1, 10)*(s@ds+ds@s-F(2, 3)*np.sum(s*ds)*eye)
            dq += F(1, 10)*(do@s+o@ds-ds@o-s@do)
            dq += F(13, 50)*(do@o.T+o@do.T-F(2, 3)*np.sum(o*do)*eye)
            square_variation = s@ds+ds@s
            db = square_variation@o+s@s@do-do@s@s-o@square_variation
            ddiff = 2*np.sum(o*do)-2*np.sum(s*ds)
            tb = int(np.sum((s if branch == 'strain' else o)*h))
            for i in range(3):
                for j in range(3):
                    linear[3*i+j, column] = (Iv(int(ds[i, j]))+slope*tb*(Iv(int(s[i, j]))-q*quadratic[i, j])
                        -q*dq[i, j]+Iv(beta)/q.square()*int(10*db[i, j]-5*ddiff*s[i, j]))
        tangent[branch] = linear
        symmetric = (linear+linear.T)*F(1, 2)
        for coefficient, form in zip(coefficients, forms):
            symmetric += form*F(coefficient, 100)
        for i in range(9):
            symmetric[i, i] -= margin
        pivots = interval_ldl(basis.T@symmetric@basis)
        if pivots is None:
            raise RuntimeError('Reference interval matrix is not certified positive: '+branch)
        report['branches'][branch] = dict(minor_weights=[str(F(c, 100)) for c in coefficients],
            LDL_pivot_intervals=[value.pair() for value in pivots])
    return report, tangent


def interval_tangent(gradient, k, epsilon, beta, branch):
    """Full d(-R_beta)/dA enclosure, differentiating every C and cubic factor."""
    eye = interval_array(np.eye(3, dtype=int))
    t = k/epsilon
    s, o = (gradient+gradient.T)*t, (gradient-gradient.T)*t
    ss, oo = squared_norm(s), squared_norm(o)
    q = ((ss if branch == 'strain' else oo)/2).sqrt()
    if q.lo <= 0:
        raise ValueError('Finite-box branch contains zero invariant; use analytic fallback')
    c, cp = cmu_interval(q), cmu_derivative_interval(q)
    quadratic = (s@s-eye*(ss/3))*F(-1, 10)+(o@s-s@o)*F(1, 10)+(o@o.T-eye*(oo/3))*F(13, 50)
    cubic = (s@s@o-o@s@s)*(-10)+s*(5*(oo-ss))
    linear = np.empty((9, 9), dtype=object)
    for column in range(9):
        h = np.zeros((3, 3), dtype=int)
        h.flat[column] = 1
        ds, do = interval_array(h+h.T)*t, interval_array(h-h.T)*t
        dss, doo = 2*inner(s, ds), 2*inner(o, do)
        dq = (dss if branch == 'strain' else doo)/(4*q)
        dquad = (s@ds+ds@s-eye*(dss/3))*F(-1, 10)+(do@s+o@ds-ds@o-s@do)*F(1, 10)
        dquad += (do@o.T+o@do.T-eye*(doo/3))*F(13, 50)
        square_variation = s@ds+ds@s
        dcubic = (square_variation@o+s@s@do-do@s@s-o@square_variation)*(-10)
        dcubic += ds*(5*(oo-ss))+s*(5*(doo-dss))
        result = (ds*c+s*(cp*dq)-dquad*c-quadratic*(cp*dq)
                  -(dcubic*(c*c*c)+cubic*(3*c*c*cp*dq))*beta)*k
        linear[:, column] = result.ravel()
    return linear, q


def gershgorin_lower(matrix):
    return min(rd(lambda i=i: matrix[i, i].lo-sum(
        (matrix[i, j].absmax() for j in range(len(matrix)) if i != j), Decimal(0)),
        ROUND_FLOOR) for i in range(len(matrix)))


def best_null_lower(symmetric, beta, branch, iterations=28, weight_scale=1):
    """A certified lower bound; bisection is only a search for an LDL witness."""
    forms, _, basis = null_forms()
    margin, weights = REFERENCE_WEIGHTS[F(beta)]
    augmented = symmetric.copy()
    for coefficient, form in zip(weights[branch], forms):
        augmented += form*(F(coefficient, 100)*F(weight_scale))
    projected, metric = basis.T@augmented@basis, basis.T@basis
    def proof(value):
        return interval_ldl(projected-metric*Iv(value))
    low, high = Decimal(-1), Iv(margin).lo
    for _ in range(20):
        if proof(low) is not None:
            break
        low *= 10
    else:
        raise RuntimeError('Could not enclose a finite negative lower bound')
    if proof(high) is not None:
        low = high
    else:
        for _ in range(iterations):
            mid = (low+high)/2
            if proof(mid) is not None:
                low = mid
            else:
                high = mid
    pivots = proof(low)
    if pivots is None:
        raise RuntimeError('Lost selected interval proof')
    return low, dict(bound_over_nut0=str(low), minor_weights=[str(F(c, 100)*F(weight_scale)) for c in weights[branch]],
                     LDL_pivot_intervals=[value.pair() for value in pivots],
                     note='Null forms vanish on every rank-one gradient; trace-free restriction enforces incompressibility.')


def interval_stress(gradient, k, epsilon, beta, gamma=0, reference_strain=None, branch='physical'):
    eye = interval_array(np.eye(3, dtype=int))
    t = k/epsilon
    s, o = (gradient+gradient.T)*t, (gradient-gradient.T)*t
    ss, oo = squared_norm(s), squared_norm(o)
    qs, qo = (ss/2).sqrt(), (oo/2).sqrt()
    q = qs if branch == 'strain' else qo if branch == 'rotation' else Iv(max(qs.lo, qo.lo), max(qs.hi, qo.hi))
    c = cmu_interval(q)
    quadratic = (s@s-eye*(ss/3))*F(-1, 10)+(o@s-s@o)*F(1, 10)+(o@o.T-eye*(oo/3))*F(13, 50)
    cubic = (s@s@o-o@s@s)*(-10)+s*(5*(oo-ss))
    result = (eye*F(2, 3)-s*c+quadratic*c+cubic*(beta*c*c*c))*k
    if reference_strain is not None:
        result -= (gradient+gradient.T-reference_strain)*Iv(gamma)
    return result


def finite_box_certificate(beta, radius, k_radius, epsilon_radius, target='0.1', ratio_splits=8):
    beta, radius = F(beta), F(radius)
    kr, er = F(k_radius), F(epsilon_radius)
    if not (0 <= radius < 1 and 0 <= kr < 1 and 0 <= er < 1):
        raise ValueError('Box radii must be nonnegative and below one')
    if int(ratio_splits) < 1:
        raise ValueError('At least one scalar-ratio subdivision is required')
    a0 = np.zeros((3, 3), dtype=int)
    a0[0, 2] = 1
    gradient = np.empty((3, 3), dtype=object)
    for index in np.ndindex(gradient.shape):
        gradient[index] = Iv(int(a0[index]))+Iv(-radius, radius)
    k, epsilon = Iv(Q_LO, Q_HI)*Iv(1-kr, 1+kr), Iv(1-er, 1+er)
    result = dict(beta=str(beta), gradient_center=a0.tolist(), gradient_entry_radius=str(radius),
                  k_over_k0=[str(1-kr), str(1+kr)], epsilon_over_epsilon0=[str(1-er), str(1+er)],
                  normalization='A*=A/G0; k*=Qstar*k/K0; epsilon*=epsilon/epsilon0; tangent_phys=nut0*tangent_star',
                  G0='u_star/(kappa*r)', K0='Qstar*u_star^2',
                  epsilon0='u_star^3/(kappa*r)', nut0='u_star^2/G0',
                  branch_scope='both complete formal max-invariant branches, hence every active branch and tie convex combination',
                  branches={}, bound_interpretation='The selected gamma is sufficient for this conservative bound; it is not necessary, minimal, or a proof that a smaller gamma fails.',
                  ratio_splits=int(ratio_splits), gradient_entry_radius_float=float(radius))
    lower_bounds = []
    for branch in ('strain', 'rotation'):
        ratio = k/epsilon
        leaves = []
        for index in range(int(ratio_splits)):
            width = (Iv(ratio.hi)-Iv(ratio.lo))/ratio_splits
            left = Iv(ratio.lo)+width*index
            right = Iv(ratio.lo)+width*(index+1)
            local_ratio = Iv(left.lo, right.hi)
            # Both simple products enclose k^2/epsilon on this ratio slice;
            # their intersection preserves correlations at no extra dimension.
            first, second = k*local_ratio, epsilon*local_ratio.square()
            local_T = Iv(max(first.lo, second.lo), min(first.hi, second.hi))
            linear, q = expanded_interval_tangent(gradient, local_ratio, local_T, Iv(beta), branch)
            symmetric = (linear+linear.T)*F(1, 2)
            full_bound = gershgorin_lower(symmetric)
            center_q = Iv((q.lo+q.hi)/2)
            nt = cmu_interval(center_q)*Iv((local_T.lo+local_T.hi)/2)
            weight_scale = (nt.lo+nt.hi)/2
            null_bound, null_proof = best_null_lower(symmetric, beta, branch, weight_scale=weight_scale)
            selected = max(full_bound, null_bound)
            leaves.append(dict(ratio=local_ratio.pair(), T=local_T.pair(), q=q.pair(),
                full_symmetric_gershgorin_lower=str(full_bound), null_form_proof=null_proof,
                selected_rank_one_lower=str(selected),
                tangent_intervals=[[entry.pair() for entry in row] for row in linear]))
        selected = min(Decimal(leaf['selected_rank_one_lower']) for leaf in leaves)
        lower_bounds.append(selected)
        result['branches'][branch] = dict(selected_rank_one_lower=str(selected), ratio_leaves=leaves,
            tangent_intervals=[[[str(min(Decimal(leaf['tangent_intervals'][i][j][0]) for leaf in leaves)),
                                str(max(Decimal(leaf['tangent_intervals'][i][j][1]) for leaf in leaves))]
                               for j in range(9)] for i in range(9)])
    bound = min(lower_bounds)
    required = max(Decimal(0), rd(lambda: Decimal(str(target))-bound, ROUND_CEILING))
    result.update(certified_lower_over_nut0=str(bound), target_over_nut0=str(target),
                  required_gamma_over_nut0=str(required),
                  qualified_property='frozen-k/epsilon incompressible momentum principal-form lower bound only')
    # A covariance witness in the same box, using a fixed rank-one test vector.
    # Negative upper endpoint is a rigorous counterexample, not an eigenvalue scan.
    shear = interval_array(a0)
    shear[0, 2] = Iv(1)+Iv(radius)
    reference = interval_array(a0+a0.T)
    stress = interval_stress(shear, Iv(Q_LO, Q_HI), Iv(1), Iv(beta), required, reference)
    variance = (stress[0, 0]+2*stress[0, 2]+stress[2, 2])/2
    result['regularized_covariance_witness'] = dict(
        normalized_gradient=[[entry.pair() for entry in row] for row in shear],
        k_over_k0=1, epsilon_over_epsilon0=1, unit_direction='(e_x+e_z)/sqrt(2)',
        variance_interval=variance.pair(), certified_negative=variance.hi < 0,
        implication='A negative upper endpoint rejects this sufficient-gamma choice on this state set; it proves neither impossibility of a smaller successful gamma nor an all-state covariance result.')
    return result


def analytic_envelope_certificate(q_max='8', t_max='20', t_min='0', beta=F(2, 5), target='0.1', partitions=64):
    """Norm-bound fallback covering all A, including q=0, at fixed k/epsilon.

    T=k^2/epsilon has viscosity units. q=max(||S~||,||Omega~||)/sqrt(2).
    For incompressible rank-one |H|=1: |dq|<=k/epsilon and
    |delta S~|=|delta Omega~|=sqrt(2)*k/epsilon. Orthogonal deviatoric
    projection and Frobenius submultiplicativity yield Q<=.56 M^2,
    delta Q<=1.12 sqrt(2) M t, B<=25 M^3,
    delta B<=85 sqrt(2) M^2 t, M=sqrt(2)*q_max. The nonnegative
    C*t*deltaS term is retained through C_min*T_min. All remaining
    terms are bounded in absolute value. This is deliberately conservative.
    """
    qm, tm, tn, beta = Iv(q_max), Iv(t_max), Iv(t_min), Iv(F(beta))
    if qm.lo <= 0 or tm.lo <= 0 or tn.lo < 0 or tn.hi > tm.lo:
        raise ValueError('Invalid finite analytic envelope')
    c_upper, cp_upper, c_lower = Decimal(0), Decimal(0), Decimal(1)
    partitions = int(partitions)
    if partitions < 1:
        raise ValueError('Positive partition count required')
    for i in range(partitions):
        lower = Iv(qm.hi)*F(i, partitions)
        upper = Iv(qm.hi)*F(i+1, partitions)
        interval = Iv(lower.lo, upper.hi)
        c, cp = cmu_interval(interval), cmu_derivative_interval(interval)
        c_upper, cp_upper, c_lower = max(c_upper, c.hi), max(cp_upper, cp.absmax()), min(c_lower, c.lo)
    c, cp, cl = Iv(c_upper), Iv(cp_upper), Iv(c_lower)
    root2, m = Iv(2).sqrt(), Iv(2).sqrt()*qm
    negative = cp*m+cp*Iv(F(14, 25))*m.square()+c*Iv(F(28, 25))*root2*m
    negative += beta*(3*c.square()*cp*25*m*m*m+c*c*c*85*root2*m.square())
    lower = (cl*tn-tm*negative).lo
    gamma = max(Decimal(0), rd(lambda: Decimal(str(target))-lower, ROUND_CEILING))
    return dict(beta=str(F(beta.lo)), q_max=str(qm.hi), T_min=str(tn.lo), T_max=str(tm.hi),
        T_definition='k^2/epsilon; dimensional viscosity units, or normalized by physical nut0',
        derivative_scope='all velocity gradients with q<=q_max, positive k,epsilon and T_min<=k^2/epsilon<=T_max; both max branches and q=0 limits',
        coefficient_bounds=dict(C_min=str(c_lower), C_max=str(c_upper), abs_Cprime_max=str(cp_upper),
                                interval_q_partitions=partitions),
        inequalities=dict(M='sqrt(2)*q_max', Q_norm='.56*M^2',
                          delta_Q_norm_over_t='1.12*sqrt(2)*M', B_norm='25*M^3',
                          delta_B_norm_over_t='85*sqrt(2)*M^2',
                          delta_q_abs_over_t='1 for unit incompressible rank-one H'),
        negative_norm_constant=negative.pair(), certified_lower=str(lower),
        target=str(target), required_gamma=str(gamma),
        bound_interpretation='This gamma is sufficient for a conservative analytic bound, not necessary or minimal. Failure of this overdiffusive choice does not rule out another qualified gamma.',
        scope='Sufficient momentum principal-part condition only; large gamma can destroy covariance and physical accuracy.')


def cmu_float(q):
    return .3/(1+.35*q**1.5)*(-np.expm1(-.36*np.exp(.75*q)))


def log_slope_float(q):
    a, u = .35*q**1.5, .36*np.exp(.75*q)
    e = np.exp(-u)
    return -1.5*a/(1+a)+.75*q*u*e/(1-e)


def floating_checks():
    lower, upper = 3.5, 4.
    if not 1+log_slope_float(lower) > 0 > 1+log_slope_float(upper):
        raise RuntimeError('Missing shear threshold bracket')
    for _ in range(64):
        mid = (lower+upper)/2
        if 1+log_slope_float(mid) > 0:
            lower = mid
        else:
            upper = mid
    return dict(q_critical=(lower+upper)/2, factor_at_q5=1+log_slope_float(5),
                interpretation='Independent floating checks; never substituted for interval envelopes.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--beta', action='append', default=None)
    parser.add_argument('--boxes', action='store_true', help='Include finite interval envelopes and gamma requirements')
    parser.add_argument('--radii', default='.01,.05,.1,.25')
    parser.add_argument('--k-radius', default='.5')
    parser.add_argument('--epsilon-radius', default='.5')
    parser.add_argument('--target', default='.1')
    parser.add_argument('--q-max', default='8')
    parser.add_argument('--T-max', default='20')
    parser.add_argument('--ratio-splits', type=int, default=8)
    args = parser.parse_args()
    report = dict(schema=1, kind='M8 constitutive mathematical certificates',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  started_utc=datetime.now(timezone.utc).isoformat(),
                  verified_property='frozen-k/epsilon incompressible rank-one momentum principal-part bounds',
                  excluded=['all-state covariance realizability', 'outer turbulence positivity',
                            'time-discrete positivity', 'coupled-mode stability', 'terrain accuracy'],
                  reference=[], finite_boxes=[], analytic_envelopes=[])
    betas = args.beta or ['2/5', '9/20']
    try:
        for value in betas:
            beta = F(value)
            certificate, _ = reference_certificate(beta)
            report['reference'].append(certificate)
            if args.boxes:
                for radius in args.radii.split(','):
                    report['finite_boxes'].append(finite_box_certificate(
                        beta, radius, args.k_radius, args.epsilon_radius, args.target, args.ratio_splits))
            report['analytic_envelopes'].append(analytic_envelope_certificate(
                args.q_max, args.T_max, beta=beta, target=args.target))
        report['floating_shear_checks'] = floating_checks()
        report['status'] = 'certified_momentum_bounds_only'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['finished_utc'] = datetime.now(timezone.utc).isoformat()
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({key: value for key, value in report.items()
                      if key not in ('reference', 'finite_boxes', 'analytic_envelopes')}, indent=2))


if __name__ == '__main__':
    main()
