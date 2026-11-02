# M6 — Prescribed upstream transition and adapted baseline

Status: implementation and Linux qualification in progress.
Branch: `codex/cls-epawf-revised`; M5 base `71b5e1e`.

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: run; source implementation and independent numerical verification.
- Origin Date: 2026-11-03 (Asia/Shanghai).
- Verification Status: UNVERIFIED until the retained Linux receipts pass.
- Version Label: cls-epawf-revised-m6-v1.

The user authorized execution and closure of the approved M6 scope. The
revised note and contract v1 remain frozen. Tailscale SSH to
`jay@100.105.247.110` is the only transport used. The shared Linux HIT source,
build and installation are preserved; the isolated revised-project source,
build and run directories remain authoritative for numerical conclusions.

The implementation and acceptance work is:

1. Implement analytic C2 spatial, state and admissibility gates, all seven
   independent vector/scalar sensors, their first/second directional actions,
   and a prescribed physical wall coordinate. Keep physical positions fixed
   under refinement; reject correction support without an exactly flat fetch.
   A cosine hill uses signed wall arc length and explicit foot patches.
2. Evaluate the current and stored projected flat reference with identical
   native physical maps at H and fixed independent outer distances r_a,r_b.
   Retain pressure/body and convection separately. Qualify the reference on
   both unit admissibility plateaus. Inactive branches evaluate the original
   profile directly without solving a corrected cutoff.
3. Couple the weight to current K, cutoff, source, dependent epsilon and full
   trace constraint. The outer sensor polynomial contains trace coefficients,
   so the trace remains an implicit coupled constraint in stage Newton.
   Differentiate the complete residual, including those columns and the
   physical trace gradient term (E_c-E_A) grad(chi).
4. Project the baseline trace and correction separately. Use native
   endpoint/edge preserving face transforms and constrain the correction to
   zero on every designated baseline face/shared mode. Require a physical
   buffer before the hill foot and a native face boundary at plateau start.
   Evaluate analytic chi in physical space; do not project the gate itself.
5. Verify equal-input chi=eta=0 identity using two instances of the same full
   curved CG operator. Compare unmasked residuals, constraints, traces and
   fluxes; retain storage, pressure, normal momentum, full-strain stress and
   retained k transport. Different upstream histories are separate comparisons.
6. Qualify kernel arithmetic with an independent high-precision oracle,
   native 2D/full-3D sensors and weak operators with independent Jacobi bases,
   all-mode correction/full-budget identities, trace projections, directional
   refinement, and deliberate omitted-term controls. Run strong manufactured
   transition cases, spatial refinement, full/condensed stages, MPI ownership
   and native full-history restart. Run the appropriate earlier regressions.

Initial acceptance tolerances remain those in the approved plan: isolated
scaled identities 1e-10, coupled constraints/stages 1e-9, MPI/restart full-field
comparisons 1e-8. Directional checks must demonstrate a stable difference-step
range and separately cover gate endpoints and inactive inadmissible branches.
Every terrain/transition run requires eta=0. No constitutive retuning, floors,
clipped cutoffs, independent inner epsilon PDE or spatial residual anchor is
part of this milestone. The original stability and positivity gates remain
failed, as specified; M7 and the conditional M8 campaign remain separate.
