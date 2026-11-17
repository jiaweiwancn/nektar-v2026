# WJ-H27 amended candidate v2 — authorized contract

The user explicitly approved the exact amended candidate in
`wj-h27-model-amendment-proposal.md` for implementation and full qualification.
The original specified-tensor feasibility result and closed account remain
immutable historical evidence. This approval changes the wall tensor; it does
not waive any H27 closure gate or qualify the proposed model.

## Scientific definition

For fixed physical wall normal n, P=I-nn^T and M=grad(U)+grad(U)^T, use

```
bA = -nuA P M n
Riso = 2k I/3 + bA tensor n + n tensor bA
Rwall = (1-chi) Riso + chi RW(k,Ew,grad(U))
mu = nu + (1-chi) nuA + chi nuW
V = nu M - Rwall + 2k I/3
Trem = V - mu M
Pactual = -Rwall : grad(U)
```

Use the original geometry-only quintic chi on [.4H,.75H], full WJ in the
open collar, physical nu in both domains, the same selected pure WJ root and
DevStrainRawWorkV1, and independent Cd. Old `Adapt` and H26 paths remain
unchanged; the new function/model has a distinct v2 name and identity.

Keep the entire imposed historical inlet, H=.15, hs=.10854503542535461,
d0=.0875, physical calibration and units. Transport k over the whole domain,
epsilon only outside H. Ew is dependent, A(K)=K^(3/2)/cl,
Ew=A/max(r,c), and outer epsilon(H)=Ew(H). Inner scalar diffusion remains
nu+Cd*k^2/(sigma_k*Ew), with the unamended source Pactual-Ew; retain the
stated outer scalar equations, conservative advection and common fluxes.
There is no EP, clipping, hidden floor or actual-flow balance forcing.

The selected counterfactual calibration column retains t=K/qe, k(H)=K,
Q(H)=Q(0)=0. Its parallel shear, production and scalar ODEs are identical
to candidate1, but the amended tensor's values, normal tractions and full
principal properties must be independently rechecked. Qualify the bounded
K/Kin in [.9,1.1] branch/table; reject out-of-range states instead of clamping.
Keep actual-flow storage/transport/interface flux distinct from reference
calibration. The independently solved outer-BVP flux is a diagnostic control,
not silently substituted for the selected zero-endpoint reference definition.

Use prescribed-reference outlet traction
`(Vxx_ref-Vzz_ref+Kin/Ne,0,tau0)`, regenerated with v2. In parallel shear,
it is `(Kin/Ne-2chi*k_ref*d(q_ref),0,tau0)`. Pin its numerical profile before
developing the application. Ground no slip/zero k diffusion; top normal
impermeability/tangential total traction tau0/zero k diffusion/fixed outward
epsilon reference flux including nu; outlet zero scalar diffusion and explicit
backflow rejection. The complete sampled AWF force retains conservative
tangential convection, mechanical pressure, actual Q3 k and body terms;
there is no sampled local time derivative feedback.

## Numerical implementation and qualification

The full original H27 objective remains active. Retain every numerical,
physical-development, finite-perturbation, restart/rotation, monitoring,
storage and complete-cost gate in `wj-h27-contract.md`. Positive sampled
proposal coverage is not a certificate of the coupled solver. At chi=0 the
exact frozen-scalar principal lower bound is nu, but covariance and the
transition/full-WJ collar require current-state checks. The prospective
physical envelopes, perturbation amplitudes/RMS/support and numerical budgets
are unchanged and cannot be chosen from subsequent successful results.

First verify native component-dependent essential spaces, actual restricted
constant-pressure action, and the actual Q3 scalar-to-force action. These
foundation checks are not a substitute for the registered fully coupled
application solver. Then verify evolving coupled continuous MMS, current
residuals, conservation, trace/root/cache/history, the same-inlet developed
2D reference and all physical/accuracy pilots. H27-B/C/D are not closed until
their complete requirements pass on a pinned recipe.

## Resources and initial admission

Use a fresh healthy-ext4 namespace
`/home/jay/build/wj-h27-coupling-r1-20261117`. Preserve the first closed H27
account (`f82e7b89d4d2e68ab177a5e9821703fad888be917cfecb5df3b06ef3565ae0c6`)
and carry its 703.6786409680499-second charge into the original aggregate
43200-second ceiling; there is no new allowance. Protect 3600 seconds for
evidence, cap actual shared-host workers at72, build with4 workers, and use
one science rank/process with all six numerical-library threads1. Retain
8/16GiB science/build RSS, 8GiB source/build and science storage, and 15GiB
free-space reserve.

Admit at most1200 seconds of initial foundation build/verification (two build
caps300s each, fixture/operator/direct/gauge/Q3 checks, focused H26 regression
and assessed native v2 tensor directions). Every failed/overlapping operation
is charged. Expand work from measured evidence; forecast/admit the whole
physical pilot/refinement matrix before expensive runs. Keep T=505.02466521278404,
terminal window [447.8483832973273,505.02466521278404], and the per-trajectory
hard gate 1.25*Ctotal<=86400. Do not execute the full Task3 horizon here.
