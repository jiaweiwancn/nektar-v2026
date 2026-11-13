# SamePhysicalWallDisplacedAcceptedCandidateV1

This is an independently assessed **proposed specification extension**, not the
original Accepted appendix and not an enabled native physical policy. The
existing original Accepted/NeutralCompositeV1 rejects nonzero aerodynamic
displacement. The numerical wall stays at physical r=0 in this candidate.

The experimental fit fixes u_tau and d before the intercept. For natural model
coefficients, K=q u_tau² and A=u_tau³/kappa. The outer reference is a logarithm
in zeta=z-d, constant K and total epsilon=A/zeta, with no roof cap. Compatible
wall coefficients remain q^-2 and kappa*q^(3/2). Neither K nor coefficients
are fitted to RMS data.

The original finite-nu resistance must first match that logarithm at H using
0<=h_s<H. A failed root is not repaired by clipping h_s or moving the wall.
Only then consider the prescribed sink

```
E(r) = A / (max(r,c)-d),        d<c<=H.
```

It is finite throughout the inner region, including r=d. Its flat integral is

```
I(c)/A = c/(c-d) + log((H-d)/(c-d)).
d(I/A)/dc = -c/(c-d)^2 < 0.
```

Its range is [H/(H-d), infinity), so P0>=A*H/(H-d) is necessary and sufficient
for the unique cutoff; equality gives c=H. No cutoff is clipped. Counterfactual
p0(K) and physical nu remain in P0. With Q(r)=integral_0^r(E-p0), define

```
k0(r)=K-integral_r^H Q(s)/D(s) ds,   D=nu+nu_t,A.
```

Then D*k0'=Q, Q(0)=Q(H)=0, k0(H)=K and epsilon(H)=A/(H-d).
Source/diffusion moments cancel together. Budget admissibility does not prove
positive k0 or full covariance; those have separate stationary-point checks.
Pi_inner=0 and Pi_outer=-Qzz give complete normal-traction matching in flat
parallel flow. Shear is u_tau², and both k endpoint fluxes vanish. Local inner
k remains nonconstant and would be transported in a future native policy.
The prescribed sink/cutoff would retain the original frozen one-sweep treatment
and zero local reaction rate; this document does not authorize trace repair.

For a curved frozen column use the actual positive J(r), not the flat cutoff:

```
I_J(c)=integral_0^H J(r)*A/(max(r,c)-d) dr.
dI_J/dc=-A/(c-d)^2 * integral_0^c J(r) dr <0.
minimum I_J=A/(H-d)*integral_0^H J(r) dr.
Q_J=integral_0^r J(E-p0),
k_col(r)=K-integral_r^H Q_J(s)/(J(s)*D(s)) ds.
```

This is a frozen normal-column identity, not a full curved solution or a
curved covariance certificate. Actual Cartesian gradients and full stresses
would still be required. d=0 recovers the original flat/metric formulas.

Placing the modeled wall at the displacement plane is a different option:
the domain starts at physical z=d, and its wall height is H-d. It discards the
physical-wall-to-displacement region from the modeled domain, changes roughness
and observation-distance interpretation, and requires a new resistance fit.
The study computes its limiting height separately; it is not presumed to cure
the original h_s<H restriction. Neither option is silently used in CT7 native
evolution. An interface above the retained 7mm rough RMS traverse would also
change the observational admission contract and is not permitted here.
