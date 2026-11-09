"""Interval certificates for A93, A95, and E94/E95.
Only the Python standard library is required.
All stress and tangent formulas are given in the main text.
"""
from decimal import Decimal, Context, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction

LOW = Context(prec=60, rounding=ROUND_FLOOR)
HIGH = Context(prec=60, rounding=ROUND_CEILING)

class I:
    """Closed interval with outward-rounded elementary arithmetic."""
    def __init__(self, lo, hi=None):
        if isinstance(lo, I):
            self.lo, self.hi = lo.lo, lo.hi
        else:
            self.lo = Decimal(str(lo))
            self.hi = Decimal(str(lo if hi is None else hi))
        assert self.lo <= self.hi

    def __neg__(self):
        return I(self.hi.copy_negate(), self.lo.copy_negate())

    def __add__(self, other):
        other = I(other)
        return I(LOW.add(self.lo, other.lo),
                 HIGH.add(self.hi, other.hi))

    __radd__ = __add__

    def __sub__(self, other):
        return self + (-I(other))

    def __rsub__(self, other):
        return I(other) + (-self)

    def __mul__(self, other):
        other = I(other)
        pairs = [(self.lo, other.lo), (self.lo, other.hi),
                 (self.hi, other.lo), (self.hi, other.hi)]
        return I(min(LOW.multiply(a, b) for a, b in pairs),
                 max(HIGH.multiply(a, b) for a, b in pairs))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I(other)
        assert not (other.lo <= 0 <= other.hi)
        pairs = [(self.lo, other.lo), (self.lo, other.hi),
                 (self.hi, other.lo), (self.hi, other.hi)]
        return I(min(LOW.divide(a, b) for a, b in pairs),
                 max(HIGH.divide(a, b) for a, b in pairs))

    def __rtruediv__(self, other):
        return I(other) / self

    def sqrt(self):
        assert self.lo >= 0
        # Decimal sqrt is correctly nearest-rounded.
        # Expanding by one representable unit encloses the root.
        lower = max(Decimal(0), LOW.next_minus(LOW.sqrt(self.lo)))
        upper = HIGH.next_plus(HIGH.sqrt(self.hi))
        return I(lower, upper)

    def out(self):
        return "[{}, {}]".format(self.lo, self.hi)


# A95 exact neutral coefficients.
root2 = I(2).sqrt()
q95 = (3 + I(61).sqrt()) / (2 * root2)
alpha_c = 3 / (root2 * q95)
f = 1 - I("4.5") / (q95 * q95)
assert f.lo > 0
CQ = f.sqrt() / (1 + 3 * q95 * q95)
beta95 = CQ * q95 * q95 * q95
c = (1 - f) / f
e = 6 * q95 * q95 / (1 + 3 * q95 * q95)
lp = c * (1 - alpha_c) + e

def form(H):
    """The polynomial P(H), extended to trace-free matrices."""
    norm = sum(x*x for row in H for x in row)
    result = I(norm)
    result -= alpha_c * (H[0][2] - H[1][1]) * (H[0][2] + H[2][0])
    result += beta95 * (
        sum(H[0][j]*H[2][j] for j in range(3))
        - sum(H[i][0]*H[i][2] for i in range(3)))
    result += beta95 * (
        lp*H[0][2] + c*H[2][0] + c*alpha_c*H[1][1]
        ) * (H[0][0] - H[2][2])
    # These two minors vanish on every rank-one matrix.
    result += I("0.2") * (H[0][1]*H[1][2] - H[0][2]*H[1][1])
    result -= I("1.6") * (H[0][0]*H[2][2] - H[0][2]*H[2][0])
    # Desired dimensionless coercivity bound: 1/50.
    result -= I("0.02") * norm
    return result

# Coordinates: xx, xy, xz, yx, yy, yz, zx, zy; zz=-xx-yy.
basis = []
for i in range(3):
    for j in range(3):
        if (i, j) == (2, 2):
            continue
        H = [[0]*3 for _ in range(3)]
        H[i][j] = 1
        if i == j:
            H[2][2] = -1
        basis.append(H)

matrix = [[I(0) for _ in range(8)] for _ in range(8)]
for i in range(8):
    matrix[i][i] = form(basis[i])
for i in range(8):
    for j in range(i):
        H = [[basis[i][a][b] + basis[j][a][b]
              for b in range(3)] for a in range(3)]
        entry = (form(H) - matrix[i][i] - matrix[j][j]) / 2
        matrix[i][j] = matrix[j][i] = entry

L = [[I(int(i == j)) for j in range(8)] for i in range(8)]
pivots = []
for i in range(8):
    pivot = matrix[i][i] - sum(
        (L[i][j]*L[i][j]*pivots[j] for j in range(i)), I(0))
    assert pivot.lo > 0, (i + 1, pivot.out())
    pivots.append(pivot)
    for row in range(i + 1, 8):
        numerator = matrix[row][i] - sum(
            (L[row][j]*L[i][j]*pivots[j] for j in range(i)), I(0))
        L[row][i] = numerator / pivot

print("A95 neutral q:", q95.out())
for i, pivot in enumerate(pivots):
    print("A95 LDL pivot", i + 1, pivot.out())
print("A95: D >= (1/50) nu_t |a|^2 |n|^2 for a.n=0.")


# A93 all-direction sufficient bound on the exact neutral root.
q93 = (57 + I(6249).sqrt()) / 40
v1 = I("2.85") / q93
v2 = I("0.15") / q93
q93cube = q93*q93*q93
beta93 = q93cube / (1000 + q93cube)
linear_loss = 15*beta93 + 45*beta93*beta93
lower93 = 1 - v1 - linear_loss*linear_loss / (4*(v1 - v2))
assert lower93.lo > Decimal("0.036")
print("A93 reference lower bound:", lower93.out())

# E94/E95 sharp simple-shear reference bound.
qE = (3 + I(41).sqrt()) / (2*root2)
lowerE = 4 / (qE*qE)
assert lowerE.lo > Decimal("0.3619")
print("E reference lower bound:", lowerE.out())


# A93 q=8: exact negative eigenvalue and positive covariance.
bad93_exact = (
    Fraction(1)
    + Fraction(3, 25)*(Fraction(304, 329) - 4*Fraction(16, 329))
    - 6*Fraction(47, 360)
    - Fraction(18, 5)*Fraction(47, 360)*Fraction(64, 63))
assert bad93_exact == Fraction(-51301, 296100)
assert bad93_exact < 0
rxx = Fraction(86, 81)
ryy = Fraction(290, 567)
rzz = Fraction(242, 567)
rxz = Fraction(-320, 987)
detxz = rxx*rzz - rxz*rxz
assert min(rxx, ryy, rzz, detxz) > 0
assert detxz == Fraction(247166356, 710169201)
print("A93 q=8 acoustic value / nu_t:", bad93_exact)
print("A93 q=8 xz covariance determinant:", detxz)

# A95 q=1: strict radicand, negative transverse value, positive R.
a0 = I("6.5")
Cbad95 = 1 / (a0 + 3/root2)
fbad95 = 1 - I("4.5")*Cbad95*Cbad95
assert fbad95.lo > 0
CQbad95 = fbad95.sqrt() / 4
beta_bad95 = (a0*a0 + 2*a0*(3/root2)).sqrt() / 4
bad95 = 1 - beta_bad95/2
assert bad95.hi < 0
rmin95 = I(2)/3 - (CQbad95*CQbad95 + Cbad95*Cbad95).sqrt()
assert rmin95.lo > 0
assert Cbad95.lo > 0  # Here P_k=C_mu.
print("A95 q=1 acoustic value / nu_t:", bad95.out())
print("A95 q=1 minimum covariance eigenvalue:", rmin95.out())

# E pure strain: exact eigenvalue, covariance and production.
badE = Fraction(-19, 68)
assert badE < 0
assert min(Fraction(4, 51), Fraction(7, 51), Fraction(91, 51)) > 0
assert Fraction(1084, 34) > 0
print("E pure-strain acoustic value / nu_t:", badE)

# A93 pure rotation violates covariance at tau*omega=8.
rotation_ryy = I(2)/3 - I(38)*64/3000
assert rotation_ryy.hi < 0
print("A93 pure rotation R_yy/k:", rotation_ryy.out())

# A93 can have negative production while R remains positive.
root3 = I(3).sqrt()
Cprod = (I(2)/3) / (I(41)/4 + 2*root3)
Dprod = 1000 + 24*root3
r1 = I(2)/3 - 2*Cprod + I(466)/(3*Dprod)
r3 = I(2)/3 + 4*Cprod - I(932)/(3*Dprod)
production93 = 12*Cprod - I(932)/Dprod
assert r1.lo > 0 and r3.lo > 0
assert production93.hi < 0
print("A93 positive-covariance negative production:", production93.out())
print("PASS: all stated interval and exact-sign checks succeeded.")
