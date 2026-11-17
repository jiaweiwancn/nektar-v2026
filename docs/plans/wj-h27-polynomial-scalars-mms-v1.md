# H27 PolynomialScalarsV1 — prospective continuous verification fixture

This optional verification fixture removes the historical scalar target's
quartic/rational spatial representation errors from a BDF2 time-order study.
It changes no application model, inlet, boundary law, numerical threshold or
physical pilot. No numerical execution or qualification is reported here.
The default variant remains `Historical`, including its original formulas and
zero prescribed-body default. The base formula identity remains
`WJH27ContinuousCoupledMMSV1`; the separate explicit variant is part of the
fixture/run identity.

## Exact fields and native spaces

Use the same flat rectangle `x in [0,1]`, `r=z-d0 in [0,R]`, `R=.5`,
`H=.15`, and the same qualified reference table. Define

```
s(x,t) = sqrt(Kin) [1 + .012 sin(.6t) + .008 x cos(.8t)]
K(x,t) = s(x,t)^2
A(K) = K^(3/2)/cl = s(x,t)^3/cl
eta(t) = .005 [1 + .2 sin(.7t)]
zeta(t) = .01 sin(.5t)
F(r) = r(r-H)/(H R)
g(r) = (r-H)^2/R^2
k(x,r,t) = K(x,t) [1 + eta(t) F(r)]
epsilon(x,r,t) = A(K)/(H) [1 + zeta(t) g(r)], r >= H.
```

Retain the historical streamfunction, pressure and optional prescribed body
exactly:

```
psi = .02 [1+.1 cos(.9t)] r^2/2
    + .05 [1+.2 sin(.7t)] x^2(1-x) r^2(R-r)
u = psi_r, w = -psi_x
Pi = .0002 [1+.3 sin(.8t)] [x^2+r^2/R^2-2/3]
b = bodyAmplitude * (2e-4(1+x+.5r)cos(.3t), 0,
                    3e-4(1+r)sin(.2t)).
```

The exact spaces on affine native rectangles are `u in Q3xQ2`,
`w in Q2xQ3`, `Pi in Q2xQ2`, `k in Q2xQ2`, and outer
`epsilon in Q3xQ2`. They are contained in the actual Q3 velocity/scalar and
Q2 pressure spaces. The interface values are exactly `k(H)=K` and
`epsilon(H)=A(K)/H`; since s is positive, the dependent trace is cubic in x.
No coefficient fit or discrete operator is used to define these fields.

For this domain, `F >= -H/(4R)=-.075`, so `1+eta F >= .99955`.
The outer epsilon shape is at least `.9951`. The interface K remains in the
same `[.9604,1.0404] Kin` range as the historical fixture. These analytic
positivity bounds do not replace the native covariance/principal screen.

## Analytic scalar derivatives and boundary data

Use `F_r=(2r-H)/(HR)`, `F_rr=2/(HR)`,
`g_r=2(r-H)/R^2`, and `g_rr=2/R^2`. The implemented k jet is

```
k_x  = K_x (1+eta F)            k_r  = K eta F_r
k_t  = K_t (1+eta F)+K eta_t F
k_xx = K_xx(1+eta F)            k_xr = K_x eta F_r
k_rr = K eta F_rr.
```

Writing `q=(1+zeta g)/H`, epsilon derivatives are
`epsilon_i=A_i q` in x, `epsilon_r=A zeta g_r/H`,
`epsilon_t=A_t q+A zeta_t g/H`,
`epsilon_xx=A_xx q`, `epsilon_xr=A_x zeta g_r/H`, and
`epsilon_rr=A zeta g_rr/H`.

All velocity faces retain analytic essential data. k remains essential only
on the left, with exact outward flux `j_k=-Dk grad(k).n` on natural faces.
At the ground, `k_r=-K eta/R` and therefore
`j_k(ground)=-Dk K eta/R`: it is finite and nonzero. Reusing the historical
zero-ground-flux assumption would be an incorrect fixture. At H,
`k_r=K eta/R`; the actual k diffusivity is continuous because both destruction
traces equal A/H, so the k normal flux also remains continuous.

Outer epsilon remains essential on the left and H. Its exact H normal
derivative is zero, while the dependent inner Ew derivative at H is generally
`-A/H^2`. The two values match, but their normal derivatives need not.
Use exact outer natural fluxes at right/top. The epsilon H flux is a reaction,
not an independently imposed condition. Scalar weak natural loads retain
the sign `-integral(phi*j_out)`.

## Continuous forcing and regularity

The unchanged inner destruction law is `Ew=A(K)/max(r,c(K))`; it is not the
polynomial outer epsilon continued to the wall. Use actual k in
`Dq=nu+Cd*k^2/(sigma_q*E)`, the current qualified reference cutoff and auxiliary
viscosity, and the full amended WJ stress/production. Every scalar derivative
in the strong source follows the selected continuous variant.

```
AWF = I_(r<=h) [2uH uxH + Pi_xH - 2kxH/3 - bxH] e_x
f_M = U_t + U.grad(U) - div(V) + grad(Pi) - b - AWF
f_k = k_t + div(Uk) - div(Dk grad(k)) - P + E
f_e = epsilon_t + div(U epsilon) - div(De grad(epsilon))
      - (1.44P-1.82epsilon)epsilon/k, only outside H.
```

The native helper still differentiates full stress using its exact directional
constitutive derivative. The independent checker instead differentiates the
physical stress and diffusive flux with fourth-order finite differences and
independent WJ/reference evaluation. Its analytic jets come from SymPy
differentiation of the formulas above, not the native jet implementation.

Stress and k flux values match across H; their derivatives and strong sources
may have bounded jumps. The moving cutoff also retains its derivative kink.
There is no manufactured distributional interface source. Split volume and
vertical-face integration at the exact cutoff, collar endpoints, h, H and
element bounds. Polynomial target representability does not make nonlinear
WJ/source coefficients or their quadrature exact.

## Prospective gates and evidence

Run the continuous native probe and independent checker for both `Historical`
and `PolynomialScalarsV1`, recording the explicit variant, body amplitude,
table/source identities and immutable outputs. Keep the existing limits:
value/jet/flux agreement `1e-10`; fine continuous forcing error `1e-8`, with
two halved fourth-order FD steps and the existing convergence/roundoff tests.
The new ground-flux and H-derivative identities are checked explicitly.
Nonzero-body coverage still requires separate bulk and sampled-body controls.

Before a temporal claim, require actual native target projection/essential
trace reproduction at representation accuracy, current coupled residuals
at the unchanged limits (flow/continuity `1e-9`, k/epsilon `1e-8`), independent
quadrature refinement and the existing direct-reference/negative controls.
The physical reference scales remain unchanged, including `Ain/r` for the
outer epsilon norm; the scale is not replaced by the manufactured epsilon.

The inherited contract requires at least three dt levels with U/k/epsilon
order at least 1.8. Freeze a common mesh, terminal time, halved time steps,
startup method and independent sampling/norms in the run ledger before
execution. Preserve the ordinary current BDF2/EX2 path and its declared
startup. Refine quadrature and algebraic tolerance independently until their
observed error contributions are resolved below the temporal differences;
do not label a representation, quadrature or solver floor as temporal order.
Report pressure/gauge and budgets as well. This fixture alone closes none of
the physical development, perturbation, cost or full H27 gates.

## Explicit API and compatibility

Native constructor:
`CoupledMMS(reference, geometry, topRadius=.5, bodyAmplitude=0,
variant=MMSVariant::Historical)`, with enum choices `Historical` and
`PolynomialScalarsV1`; invalid enum values are rejected. `Variant()` and
`VariantName()` expose the selection.

The probe accepts `--variant PolynomialScalarsV1` and optional
`--body-amplitude VALUE` in either order. The independent constructor is
`IndependentMMS(table, body_amplitude=0, variant='Historical')`.
New coupled-run reports must provide `state.mms_variant`; the current coupled
oracle rejects an absent or inconsistent variant. Historical pre-variant
continuous probe artifacts retain an explicitly labelled Historical fallback.
Frozen old source capsules and their recorded results remain unchanged.
