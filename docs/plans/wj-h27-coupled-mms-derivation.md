# H27 continuous coupled MMS V1 — verification only

This fixes the manufactured fields and amplitudes before native coupled tests.
It supplies verification forcing only; it changes no application data or PDE.
Native admissibility, independent forcing checks, spatial/quadrature refinement
and three time-step levels remain required. No numerical evaluation was run on
the Mac when these formulas and helper were authored.

Use physical `r=z-d0`, `x∈[0,1]`, `r∈[0,R]`, `R=.5`, the model's fixed H and
roughness h, and the qualified immutable reference table. Public helper
coordinates are physical x,z; jet `dr` is the physical z derivative. All tensor
indices use x=0 and z=2. The physical body force is zero; manufactured momentum
forcing is a distinct verification load and must never be fed into sampled AWF.

## Fields and exact dependent trace

Define

```
a(t)=sqrt(Kin)[1+.012 sin(.6t)]
b(t)=sqrt(Kin) .008 cos(.8t)
s(x,t)=a(t)+b(t)x;  K(x,t)=s²
η(t)=.005[1+.2 sin(.7t)]
F(r)=r²(r−H)²/(H²R²)
k(x,r,t)=K(x,t)[1+η(t)F(r)]
ζ(t)=.01 sin(.5t)
ε(x,r,t)=A(K)/r [1+ζ(t)(r−H)²/R²], r≥H
```

`s>0` for every real t and x∈[0,1], and conservatively
`.9604 Kin≤K≤1.0404 Kin`. The positive k shape is at most
`1+.006 R²(R−H)²/(H²R²)` on this fixed domain; with H=.15 it remains
below `1.0744 Kin`. Epsilon's multiplier is at least `.9951`.
These scalar bounds do not certify Reynolds-stress covariance or ellipticity.

`k(H)=K`, `k_r(0)=k_r(H)=0`, and
`ε(H)=A(K)/H=s³/(cl H)` is exactly cubic in x. The outer epsilon x, r and
time derivatives at H equal those of the inner tail `Ew=A(K)/r`; there is no
inlet/interface corner mismatch or independently overwritten epsilon trace.
The scalar left boundary uses these analytic manufactured fields, distinct
from the retained physical application inlet. Epsilon has no inner PDE.

The helper returns analytic value, first derivatives, time derivative and
Hessian for k and outer epsilon. In particular, writing
`B=A(K)`, `g=(r−H)²/R²`, `f=(1+ζg)/r`,

```
B_x=A_K K_x; B_xx=A_K K_xx+A_KK K_x²; B_t=A_K K_t
f_r=ζg_r/r−(1+ζg)/r²
f_rr=ζg_rr/r−2ζg_r/r²+2(1+ζg)/r³
ε_x=B_x f; ε_r=B f_r
ε_t=B_t f+B ζ_t g/r
ε_xx=B_xx f; ε_xr=B_x f_r; ε_rr=B f_rr.
```

## Divergence-free velocity and pressure

```
S(t)=.02[1+.1 cos(.9t)]
V(t)=.05[1+.2 sin(.7t)]
X(x)=x²(1−x); Y(r)=r²(R−r)
ψ=S r²/2+V X Y
u=ψ_r=S r+V X Y_r; w=−ψ_x=−V X_x Y
Pi=.0002[1+.3 sin(.8t)] [x²+r²/R²−2/3].
```

Hence `u_x+w_r=0`, ground velocity is zero, net exterior normal flux is zero,
and Pi has zero physical volume mean. All velocity faces carry analytic
essential values in this MMS fixture; this does not replace the independently
qualified physical mixed boundary cases. Velocity is in tensor Q3 and Pi in
Q2. The quartic r dependence of k and rational outer epsilon still require a
spatial refinement study; exact velocity representation does not remove it.

## Continuous constitutive forcing

At each point use actual manufactured k and `E=Ew(K,r)` inside, transported
analytic epsilon outside. Let `aA=Auxiliary(K,r)` and use the amended v2 stress
`Vstress=nu M−R+2kI/3`, production `P=−R:grad(U)`. This Vstress is distinct
from the scalar velocity amplitude V(t) above. The helper obtains its physical
derivatives from exact `ApplyAmendedDirection` calls with

```
x direction: (k_x,E_x,∂x grad(U),aA_K K_x,dr=0)
r direction: (k_r,E_r,∂r grad(U),aA_r,dr=1).
```

The r call includes the derivative of the geometry collar. The x call includes
the current reference K dependence of both the cutoff and auxiliary viscosity.
The selected table derivative is the derivative of the accepted parameterised
reference law, not a re-fitted manufactured reference or old inlet cutoff.

With exact outer-H ownership, `BodyAtH=0`, and

```
AWF=1_(r≤h)[2u_H (u_x)_H+(Pi_x)_H−(2/3)(k_x)_H] e_x
F_U=U_t+U·grad(U)−div(Vstress)+grad(Pi)−AWF.
```

There is no sampled local time derivative and no recursive sampling of F_U as
the physical body term. For scalars use `Dq=nu+Cd k²/(sigma_q E)`, with
`sigma_k=1`, `sigma_e=1.3`, and the same physical nu. Its gradient is

```
∂i Dq=(Cd/sigma_q)[2k (∂i k)/E−k² (∂i E)/E²].
f_k=k_t+div(Uk)−Dk Δk−grad(Dk)·grad(k)−P+E
f_e=ε_t+div(Uε)−De Δε−grad(De)·grad(ε)
    −(1.44P−1.82ε)ε/k, only outside H.
```

`div(Uq)` is retained as `U·grad(q)+q div(U)` in the helper. No reference
`k_ref` enters actual diffusivity. Molecular viscosity is retained in both
scalar and complete momentum operators, independently of Cd and the momentum
split coefficient. The forcing contains neither a discrete matrix action nor
an inverse time-step term.

## Boundaries, interfaces and regularity

Keep k essential only on the left; its exact ground flux is zero and its
right/top natural loads are `j_k=−Dk grad(k)·n`. Keep epsilon essential on the
left and H, with exact right/top `j_e=−De grad(ε)·n`. The helper returns outward
fluxes; the scalar weak RHS contributes `−∮phi j`. Numerical conservative
advection still uses the current numerical U and q on every exterior face.

The fields are smooth and the value/first-derivative matches at H make complete
traction and k diffusive flux continuous there. No extra interface source is
manufactured. `Ew=A/max(r,c(K))` is continuous and piecewise smooth across the
moving cutoff, with a derivative kink. Dk and its normal flux remain continuous;
the strong scalar forcing can have a bounded jump, with no distributional delta.
At the cutoff use the declared plateau-side derivative only as a point convention;
derivative checks must avoid that measure-zero surface. The cutoff lies below
the chi=0.4H transition, where the amended momentum stress is independent of Ew.
The geometry collar is C2; AWF and its cancelling verification load have the
declared bounded jump at h. Integrate separately at c(x,t), .4H, .75H, h and H.

The analytic solution is prospective verification data, not a successful coupled
run, developed application solution, model admissibility certificate or Task3
qualification.
