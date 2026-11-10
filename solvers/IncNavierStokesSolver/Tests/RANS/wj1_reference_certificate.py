# SPDX-License-Identifier: MIT
# Extracted from the reviewed WJ note; no production evaluator is imported.
"""Self-contained exact reference certificate and off-reference WJ checks.
Requires Python 3, NumPy, SymPy, and SciPy. No data files are read.
The model uses the full 3D tensor basis and the explicit 2D cubic N root.
The certificate is exact rational arithmetic; optimized values are diagnostic.
"""
import numpy as np
import sympy as sp
from itertools import combinations

F=sp.Rational; r=sp.sqrt(143); Cstar=F(572,6561); qstar=1/sp.sqrt(Cstar)
cp=1.8

def constitutive(B, Nguess=None):
    """Dimensionless anisotropy a, effective coefficient, cubic N."""
    S=(B+B.T)/2; O=(B-B.T)/2; I=np.eye(3)
    s=np.trace(S@S); w=np.trace(O@O); iv=np.trace(S@O@O)
    if Nguess is None:
        rr=np.roots([1,-cp,-2.7*s-2*w,2*cp*w])
        N=max(z.real for z in rr if abs(z.imag)<1e-9)
    else:
        N=Nguess
    for _ in range(12):
        f=N**3-cp*N*N-(2.7*s+2*w)*N+2*cp*w
        df=3*N*N-2*cp*N-2.7*s-2*w
        d=f/df; N-=d
        if abs(d)<1e-14: break
    Q=(5/6)*(N*N-2*w)*(2*N*N-w)
    b1=-N*(2*N*N-7*w)/Q; b3=-12*iv/(N*Q)
    b4=-2*(N*N-2*w)/Q; b6=-6*N/Q; b9=6/Q
    a=(b1*S+b3*(O@O-w*I/3)+b4*(S@O-O@S)
       +b6*(S@O@O+O@O@S-(2/3)*iv*I)
       +b9*(O@S@O@O-O@O@S@O))
    return a,-.5*(b1+w*b6),N

def numeric_tangent(B):
    _,C,N=constitutive(B)
    cols=[]
    for h in np.eye(9):
        aa,_,_=constitutive(B.astype(complex)+1e-30j*h.reshape(3,3),N)
        cols.append(-aa.imag.ravel()/(1e-30*C))
    L=np.column_stack(cols)
    return L,(L+L.T)/2

# Exact derivative of -R divided by nu_eff at the reference; row-major order.
L=sp.Matrix([
 [F(386,143),0,-40*r/1201,0,0,0,10300*r/171743,0,F(100,143)],
 [0,1,0,F(53,28),0,-10*r/143,0,-125*r/4004,0],
 [10*r/143,0,F(86,1201),0,0,0,F(1886,1201),0,-10*r/143],
 [0,1,0,F(53,28),0,-10*r/143,0,-125*r/4004,0],
 [F(775,2002),0,0,0,F(53,14),0,0,0,F(775,2002)],
 [0,10*r/143,0,125*r/4004,0,1,0,F(53,28),0],
 [10*r/143,0,F(86,1201),0,0,0,F(1886,1201),0,-10*r/143],
 [0,10*r/143,0,125*r/4004,0,1,0,F(53,28),0],
 [F(100,143),0,40*r/1201,0,0,0,-10300*r/171743,0,F(386,143)]])
M=(L+L.T)/2
minor_matrices=[]
for i,j in combinations(range(3),2):
    for t,u in combinations(range(3),2):
        V=sp.zeros(9)
        V[3*i+t,3*j+u]+=F(1,2); V[3*j+u,3*i+t]+=F(1,2)
        V[3*i+u,3*j+t]-=F(1,2); V[3*j+t,3*i+u]-=F(1,2)
        minor_matrices.append(V)
weights=[F(289,100),0,0,0,F(57,50),0,0,0,F(289,100)]
Z=[]
for i,j in [(0,1),(0,2),(1,0),(1,2),(2,0),(2,1)]:
    H=sp.zeros(3); H[i,j]=r if i==0 or j==0 else 1
    Z.append(sp.Matrix(list(H)))
Z.extend([sp.Matrix(list(sp.diag(1,0,-1))),sp.Matrix(list(sp.diag(0,1,-1)))])
Z=sp.Matrix.hstack(*Z)
K=sp.simplify(Z.T*(M+sum((w*V for w,V in zip(weights,minor_matrices)),sp.zeros(9))
                         -F(1,250)*sp.eye(9))*Z)
assert all(v.is_Rational for v in K)
ll,dd=K.LDLdecomposition(hermitian=False)
pivots=list(dd.diagonal())
assert all(d>0 for d in pivots)
assert K==ll*dd*ll.T
print('EXACT CERTIFICATE: D >= nu_eff/250 on all incompressible rank-one increments.')
print('Exact positive LDL pivots:',[str(d) for d in pivots])


# Exact adverse-direction sign certificate, independent of the C++ root.
Nsymbol=sp.symbols('N')
polynomial=80*Nsymbol**3-144*Nsymbol**2-343*Nsymbol-1764
lo=F(413,100); hi=F(414,100)
assert polynomial.subs(Nsymbol,lo)<0 and polynomial.subs(Nsymbol,hi)>0
numerator=37338352*Nsymbol**2-156626001*Nsymbol+59032897
assert numerator.subs(Nsymbol,lo)>0 and sp.diff(numerator,Nsymbol).subs(Nsymbol,lo)>0
assert (240*Nsymbol**2-288*Nsymbol-343).subs(Nsymbol,lo)>0
