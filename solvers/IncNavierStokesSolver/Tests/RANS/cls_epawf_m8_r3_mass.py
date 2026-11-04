"""Exact polynomial tensor mass and certified fixed-vector dual norms.

This module is evaluated numerically only on Linux.  Geometry endpoints and
binary64 load coefficients are interpreted as their exact rationals.  A
Decimal solve supplies an approximation; an exact rational inverse-residual
and a posteriori residual bound enclose its dual quadratic form.  This is a
bound for the stated mass/vector problem, not for residual input errors or
the physical equations.  Native stored mass entries need not be these exact
polynomial integrals.
"""
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
import math


def rational(value):
    if isinstance(value, F):
        return value
    if isinstance(value, str):
        return F(Decimal(value))
    return F(float(value))


def decimal(value, precision=90, rounding=None):
    with localcontext() as context:
        context.prec = precision
        if rounding is not None:
            context.rounding = rounding
        return Decimal(value.numerator)/Decimal(value.denominator)


def interval_sqrt(lower, upper):
    if upper < 0:
        raise ValueError('Negative certified mass-dual square')
    lower = max(F(0), lower)
    with localcontext() as context:
        context.prec = 65
        lo = decimal(lower, 65, ROUND_FLOOR).sqrt().next_minus()
        hi = decimal(upper, 65, ROUND_CEILING).sqrt().next_plus()
    return [str(max(Decimal(0), lo)), str(hi)]


POLYNOMIALS = ([F(1, 2), F(-1, 2)], [F(1, 2), F(1, 2)],
               [F(1, 4), F(0), F(-1, 4)],
               [F(0), F(1, 2), F(0), F(-1, 2)])


def integral_product(a, b):
    return sum((x*y*F(2, i+j+1) for i, x in enumerate(a)
                for j, y in enumerate(b) if (i+j) % 2 == 0), F(0))


def basis_keys(lower, upper, degree):
    return [('v', float(lower).hex()), ('v', float(upper).hex())]+[
        ('e', float(lower).hex(), float(upper).hex(), index)
        for index in range(2, degree+1)]


class RationalFactor:
    """Small Decimal inverse plus an exact rational Banach bound."""
    def __init__(self, matrix):
        self.matrix, self.n = matrix, len(matrix)
        if not self.n:
            raise ValueError('Empty exact mass factor')
        n = self.n
        with localcontext() as context:
            context.prec = 90
            augmented = [[decimal(matrix[i][j]) for j in range(n)]+[
                Decimal(int(i == j)) for j in range(n)] for i in range(n)]
            for column in range(n):
                pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
                if augmented[pivot][column] == 0:
                    raise ValueError('Singular exact polynomial free mass')
                augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
                scale = augmented[column][column]
                augmented[column] = [value/scale for value in augmented[column]]
                for row in range(n):
                    if row == column:
                        continue
                    scale = augmented[row][column]
                    if scale:
                        augmented[row] = [a-scale*b for a, b in zip(augmented[row], augmented[column])]
            self.inverse = [[F(value) for value in row[n:]] for row in augmented]
        residual = [[F(int(i == j))-sum((matrix[i][k]*self.inverse[k][j]
            for k in range(n)), F(0)) for j in range(n)] for i in range(n)]
        self.inverse_residual_inf = max(sum(map(abs, row), F(0)) for row in residual)
        if self.inverse_residual_inf >= 1:
            raise ValueError('Approximate inverse cannot be certified')
        approximate_norm = max(sum(map(abs, row), F(0)) for row in self.inverse)
        self.inverse_inf_bound = approximate_norm/(1-self.inverse_residual_inf)

    def multiply_inverse(self, rhs):
        return [sum((a*b for a, b in zip(row, rhs)), F(0)) for row in self.inverse]


class TensorMass:
    """Canonical physical Modified_A mass, including native maps/signs."""
    def __init__(self, record, fixed):
        self.record, self.n = record, int(record['size'])
        elements = record['mesh']
        self.axes, self.keys, self.one_dimensional = (0, 2), [], []
        for axis in self.axes:
            cells = {(float(e['lower'][axis]), float(e['upper'][axis]), int(e['degree']))
                     for e in elements}
            if any(degree not in (1, 2, 3) for _, _, degree in cells):
                raise ValueError('Exact audit supports only the declared Modified_A degrees1..3')
            keys = sorted({key for lo, hi, degree in cells for key in basis_keys(lo, hi, degree)})
            index = {key: j for j, key in enumerate(keys)}
            mass = [[F(0) for _ in keys] for _ in keys]
            for lo, hi, degree in sorted(cells):
                local = basis_keys(lo, hi, degree)
                length = (rational(hi)-rational(lo))/2
                for i, key_i in enumerate(local):
                    for j, key_j in enumerate(local):
                        mass[index[key_i]][index[key_j]] += length*integral_product(POLYNOMIALS[i], POLYNOMIALS[j])
            self.keys.append(keys)
            self.one_dimensional.append(mass)
        self.native = {}
        indices = [{key: j for j, key in enumerate(keys)} for keys in self.keys]
        for element in elements:
            degree = int(element['degree'])
            local = [basis_keys(element['lower'][axis], element['upper'][axis], degree) for axis in self.axes]
            for j, (native, sign) in enumerate(zip(element['ids'], element['sign'])):
                x, r = j % (degree+1), j // (degree+1)
                entry = (indices[0][local[0][x]], indices[1][local[1][r]], int(sign))
                native = int(native)
                if abs(sign) != 1 or native in self.native and self.native[native] != entry:
                    raise ValueError('Inconsistent native/canonical tensor map')
                self.native[native] = entry
        if set(self.native) != set(range(self.n)):
            raise ValueError('Incomplete native tensor mass map')
        self.free = sorted(set(range(self.n))-set(fixed))
        tuples = {(self.native[i][0], self.native[i][1]) for i in self.free}
        self.fx = sorted({x for x, _ in tuples})
        self.fr = sorted({r for _, r in tuples})
        if tuples != {(x, r) for x in self.fx for r in self.fr}:
            raise ValueError('Free mass mask is not the declared Cartesian product')
        self.x = RationalFactor([[self.one_dimensional[0][i][j] for j in self.fx] for i in self.fx])
        self.r = RationalFactor([[self.one_dimensional[1][i][j] for j in self.fr] for i in self.fr])
        self.volume = sum(((rational(e['upper'][0])-rational(e['lower'][0]))*
            (rational(e['upper'][2])-rational(e['lower'][2])) for e in elements), F(0))

    def canonical(self, value):
        result = [[F(0) for _ in self.fx] for _ in self.fr]
        ix, ir = {x: j for j, x in enumerate(self.fx)}, {r: j for j, r in enumerate(self.fr)}
        for native in self.free:
            x, r, sign = self.native[native]
            result[ir[r]][ix[x]] = sign*rational(value[native])
        return result

    def solution_and_bound(self, value):
        rhs, nx, nr = self.canonical(value), len(self.fx), len(self.fr)
        middle = [[F(0) for _ in range(nx)] for _ in range(nr)]
        for x in range(nx):
            column = self.r.multiply_inverse([rhs[r][x] for r in range(nr)])
            for r in range(nr):
                middle[r][x] = column[r]
        solution = [self.x.multiply_inverse(row) for row in middle]
        # Exact rational a posteriori residual for Mr*solution*Mx^T.
        mr, mx = self.r.matrix, self.x.matrix
        applied_r = [[sum((mr[i][j]*solution[j][x] for j in range(nr)), F(0))
            for x in range(nx)] for i in range(nr)]
        residual = [[rhs[r][x]-sum((applied_r[r][j]*mx[x][j] for j in range(nx)), F(0))
            for x in range(nx)] for r in range(nr)]
        norm_rhs = sum((abs(x) for row in rhs for x in row), F(0))
        norm_residual = max(abs(x) for row in residual for x in row)
        solution_error = self.r.inverse_inf_bound*self.x.inverse_inf_bound*norm_residual
        error = norm_rhs*solution_error
        square = sum((rhs[r][x]*solution[r][x] for r in range(nr) for x in range(nx)), F(0))
        proof = dict(
            approximate_inverse_residual_inf=[str(decimal(self.x.inverse_residual_inf)), str(decimal(self.r.inverse_residual_inf))],
            full_solve_residual_inf=str(decimal(norm_residual)),
            quadratic_absolute_error_upper=str(decimal(error, 70, ROUND_CEILING)),
            coefficient_solution_error_inf_upper=str(decimal(solution_error, 70, ROUND_CEILING)),
            free_tensor_dimensions=[nx, nr], free_count=len(self.free))
        return solution, solution_error, (square-error, square+error), proof

    def norm(self, values, rate):
        _, _, (lo, hi), proof = self.solution_and_bound(values)
        denominator = rational(rate)**2*self.volume
        interval = interval_sqrt(lo/denominator, hi/denominator)
        return dict(interval_decimal=interval, midpoint=float((Decimal(interval[0])+Decimal(interval[1]))/2),
            proof=proof, scope='Exact polynomial physical mass; exact declared rounded residual coefficients; excludes residual input and native stored-mass errors')

    def functional(self, values, native_weights):
        solution, error, _, proof = self.solution_and_bound(values)
        weights = self.canonical(native_weights)
        value = sum((weights[r][x]*solution[r][x] for r in range(len(self.fr))
                     for x in range(len(self.fx))), F(0))
        bound = error*sum((abs(x) for row in weights for x in row), F(0))
        return dict(interval_decimal=[str(decimal(value-bound, 65, ROUND_FLOOR)),
            str(decimal(value+bound, 65, ROUND_CEILING))],
            midpoint=float(value), proof=proof,
            scope='Exact polynomial consistent free mass, exact declared rounded load and exact polynomial evaluation row; no bound on load inputs or native stored-mass entries')

    def action(self, value):
        """Exact rational full native mass action, for a saved double delta."""
        nx, nr = len(self.keys[0]), len(self.keys[1])
        coefficients = [[F(0) for _ in range(nx)] for _ in range(nr)]
        for native, (x, r, sign) in self.native.items():
            coefficients[r][x] = sign*rational(value[native])
        mr, mx = self.one_dimensional[1], self.one_dimensional[0]
        left = [[sum((mr[r][j]*coefficients[j][x] for j in range(nr)), F(0))
            for x in range(nx)] for r in range(nr)]
        applied = [[sum((left[r][j]*mx[x][j] for j in range(nx)), F(0))
            for x in range(nx)] for r in range(nr)]
        return [sign*applied[r][x] for _, (x, r, sign) in sorted(self.native.items())]
