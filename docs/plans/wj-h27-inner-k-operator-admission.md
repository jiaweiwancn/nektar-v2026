# F2 diagnostic inner-k transport and operator checks

The independent F1 row D is supported at Kin. Use the separate
`WJH27InnerKFixture` executable, reusing native Coupling, Scalar,
RuleFactory, Reference, full stress, wall guards and steady solver.
It accepts physical steady sessions at step and time zero. The production
factory constructs ScalarCdV2 by default. The shared implementation has
no production session option for the alternative.

The fixture selects AuxiliaryHistoricalInnerKFixtureV1. The old hashed
v2 table is still validated with the original generator, scalar-law and
zero-reference-flux semantics. It supplies retained geometry, inlet and
all physical boundary data. A separate diagnostic reference identity
binds this boundary table to the analytic $c_A(K)$ and auxiliary inner $D_k$
policy. The compiled shared kernels and diagnostic main have distinct
source hashes. Checkpoint signatures include all of these dependencies.
Only the explicit steady seed interface may import a declared v2 field
across the reference identity. Rebuild its essential lifts, Gamma,
current cutoff panels and advection. Use this seed at zero physical time.
The diagnostic is restricted to steady solves and rejects cross-policy
accepted histories. [Column result](wj-h27-inner-k-fixture-f1.md)

Use the bounded [.9,1.1] Kin interval. The analytic log-cutoff slope is
negative by bounding the endpoint factors of $s=Kd'/d$ and using
$KJ'=(1-\nu/B)(aH/B)(1/2-s)+(\nu/d-\nu/B)s$, with $B=d+aH$.
Certify the single selected roughness branch over the full interval.
Use the same moving-cutoff inverse and native cubic-owner partition with
this enclosure. Check full derivative actions $c_K,E_{w,K},E_{w,r}$,
the $D_{k,w}$ value and dependence, $\nu_A$ derivatives, and inverse radii
at endpoints, Kin and two interior held-out values. Use strictly separated
radii on both sides of the cutoff. Compare to independent analytic and
F1 multiprecision evidence in the existing 1e-10 scales.

F2 uses the retained base v2 28×100 steady coefficients to test default
parity at its same q8×12 recipe. Compare exact source19 exports with the
isolated default-v2 evaluation, including full signed fields/histories and
unchanged boundary tag. Require at most 1e-8 field error and unchanged
acceptance of current equations. Reuse the independent base audit when the
state is byte-identical. The original source19 executable remains pinned.

The changed-action test uses the retained 700-cell geometry with one radial
face added at the historical cutoff at Kin. Use native q12×16 and independent
Gauss20 and 32. Its positive test fields are deliberately nonstationary in
the equation sense. In one state set native constancy of K away from the
inlet and choose two cubic owner profiles with unequal normal slopes and
common $Q=.02A f(x)$ at H, where $f(x)=\max(0,1-|x|/2.5)$.
In the other state change K by $0.04 K_{in}f(x)(x-.30)/2.8$, and verify the
resulting H epsilon lift and moving cutoff crossing near x=.30 on both
physical owners. These diagnostic fields receive no added scalar source.
Assemble the full signed native diffusion, source, advection, natural and
total actions, including essential rows. Compare to independent owner-aware
rows with fixed scales 25A for k and the existing whole-domain epsilon
budget scale. Exercise mismatched and stale-policy controls. The expected
common flux is nonzero and the two slopes/diffusivities are distinct.
Require 1e-10 for small actions and interface trace. Check current covariance,
positive k/epsilon and supported principal response, plus independent weak
account identities for both test states. The steady residual thresholds
apply to converged fields, with unchanged momentum/divergence 1e-9 and
scalar 1e-8. The imposed 2D interface loads contain no counterfactual Q(H).

After verified F1, admit 35 seconds for materialization, 350 seconds for
the four-worker isolated build, 35 seconds for inputs, 110 seconds for
the v2 parity solve, 100 seconds for native changed actions, 140 seconds
for independent checks, and 40 seconds for evidence. The total reserved
cap is 810 seconds under the same 6,000-second fixture account. The serial
2D admission follows actual F2 timings and the prior 5,808-cell peak.
All executable and science artifacts are create-only.
