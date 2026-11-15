# H2i: fixed-k joint velocity–epsilon feasibility

15 November 2026, Asia/Shanghai. **Both predeclared fixed-k anchors have a
checked finite sampled obstruction at the actual startup limits.** The joint
experiment used 12 LP rounds and 25 separation scans in 19.41 s. Independent
saved-cut replay preserves positive lower bounds after arithmetic and native
rounding allowances. Neither diagnostic field passes native/seen-grid
qualification. No fresh GL28 audit, certificate, integrator history or physical
hill attempt was created.

This result concerns the selected k anchors, 41 retained velocity directions,
eight epsilon directions, sampled sets and physical bounds. It is not global
physical infeasibility, an ellipticity result or evidence that no hill solution
exists. The calibrated-standard versus regularized-CLS comparison remains open.

## Preserved experiment and provenance

Fetch matched reviewed parent `409f7d98b21ddee999dfe67b1600810b7091b769`.
Source/Git work used `.worktrees/pimple-one-sweep` on
`codex/cls-epawf-pimple-one-sweep`. Other worktrees, the shared Linux checkout's
local changes and historical runtimes were preserved. Tested source blobs and
implementation commits are bound in
[source records](cls-cao-tamura-budget-h2i-evidence/source-revisions.json) and
[commit binding](cls-cao-tamura-budget-h2i-evidence/commit-binding.json).

The [plan](cls-cao-tamura-budget-h2i-plan.md) and
[predeclared policy](cls-cao-tamura-budget-h2i-evidence/policy.json) retain H2g's
calibration, physical boundary roles, absolute pressure convention, source
quadrature and exact 2804-slot mesh. The 3000-slot production guard remains.

| Quantity | Deployed CLS | Calibrated standard |
|---|---:|---:|
| H_A=z_K=z_m, m | .0065 | .0065 |
| u_tau, m/s; nu, m²/s | .07426793792419734; 3.8095238095238094e-6 | same |
| kappa; fitted z0_eff, m | .4405539051866587; .00012222264395414748 | same |
| K0, m²/s² | .017598891430879275 | .01838575534504143 |
| sigmaEpsilon | 1.3 | 1.3478315512167742 |
| Physical nonlinear / regularization / Yap | retained; beta=.4, Gamma=1, gamma0=0 | all absent; Cmu=.09 |
| Deployed h_s, m | .0022472751006336437 | .0022035111852050995 |
| Wall cMu; cLength | .09822789532821126; 2.5108624752197968 | .09; 2.6811256851969465 |
| Independently solved flat cutoff, m | .0008312858196541334 | .0008312858196541342 |

The seven-point fit, weights, viscosity, friction velocity and data are unchanged.
AnalyticalNeutralInflowV1 uses physical incoming height; nearest-floor distance,
inner normal distance and outer chart coordinate remain separate. Both models
retain Linear inner stress, inner molecular diffusion and their own calibrated
accepted budget AWF. The transported source remains P_actual−E_C; p_A0(K)
defines the local cutoff. No hill observation or new stress law was used. The
valid H2g frozen-law table is retained.

Verified standard identities:

- Cold: `038fc076cebbea41bcb7cc349147eacc2b938a9051827bbb386baf7c3cb3716c`.
- Original B: `c83da50b65d2e9980325d365779a9d53394fbf10e6431d1b87f79728366ae13d`.
- A diagnostic: `fc0c75fad0a407a0e8337e5dc02379b51252d21a02847fe3eee184dd2cc1ab9b`.
- B-anchor diagnostic: `f857327e54618ccea62f3a9cd69c1877db4afa68226c5e00b8df9af0644ca2e8`.

Candidates start from exact archived-B copies. B velocity is never reconstructed
from cold search coefficients. Only permitted velocity and epsilon rows receive
search increments; k is fixed to its anchor. Pressure, inner k, endpoint K,
inlet/interface scalars, dependent epsilon, reference, defect and prescribed data
retain their native contracts, including signed-zero checks. These initializer
restrictions introduce no evolution-time clamp.

## Qualified context and two anchors

The new probe-only `FixedKJointContextV1` uses the complete existing conditioned
41/8/8 pool, without the old 32-direction enrichment exporter. Binary64 arrays
retain **450871 native memberships**, plus GL20/24 points, as **263523 evaluation
records**. Native geometry, element, region, radial/tangential sides, field block,
actual/dependent epsilon and family memberships remain explicit.

[Context evidence](cls-cao-tamura-budget-h2i-evidence/context.json) records:

- All 105 regional continuity rows and rank-105 independent Q are checked.
  Direction errors are at most 1.28e-13 for D and 1.18e-12 for Q.
- Three legal nonzero trials reproduce direct native gradients to 1.51e-13
  relative and scalars to 1.52e-15 relative. The independent standard tensor
  minimum agrees to 1.13e-17.
- Physical consistent-mass and GL16/20/24 divergence forms include the actual
  B constant and cross terms. Protected profiles and rows are checked.

The maximal seen-set scalar retention is theta=**.9972828525348357**; the
predeclared inward guard gives **.9972828525338356**. The limiting B/cold k ratio
is .4996404493151861; the retracted value is .5010000000005005. The same theta
contracts permitted B scalar coefficients toward cold in their saved space.
Changing the origin does not reset scalar or velocity allowances.

Anchor A fixes retracted k. Anchor B fixes exact cold outer k, with the same
protected inner/interface values and retracted epsilon starting point. Epsilon
then varies in its original eight-dimensional space. Anchor B ran only after
A's actual-limit obstruction. A's proof includes outer constraints changed by
k, so an unchanged-constraint replay could not replace the second solve. There
was no adaptive anchor or basis search.

GL24 k/cold-k spans [.501147862,1.997602325] for A; integrated TKE is **+1.95193%**
relative to cold. B-anchor k and its integrated TKE are exactly cold. These are
initial-field changes, not positivity-correction supplies. Permitted rows,
anchor hashes and scalar/TKE data are retained in the machine-readable records.

## Convex problem and independently checked obstruction

Variables are 41 signed velocity and eight signed epsilon coefficients. With
S=G+Gᵀ, no trace subtraction or sided averaging, outer constraints are

`[epsilon*(2*k/3−delta)*I − .09*k²*S] / (K0*epsilon_cold) + t*I >= 0`.

Inner constraints are

`[(2*k/3−delta)*I − nu_A*S] / K0 + t*I >= 0`.

Both matrices are affine. The y eigenvalue is retained. Phase I minimizes t
alone: **transformed-tensor slack**, not a raw covariance eigenvalue. Rayleigh
cuts and quadratic tangents form a necessary LP relaxation; full true PSD,
scalar, response and quadratic scans qualify fields. An LP primal cannot do so.

Limits remain total velocity L2<=.12, component increment<=.5, divergence
L2<=.06030809942199109 and scalar ratios [.5,2], all relative to cold. The
existing 1e-8 allowance applies only to independently integrated GL20/24
divergence checks. Finite coefficient boxes follow from the mass ellipsoid and
eight independent scalar rows; no trust boxes were introduced.

Buffered construction uses 2*delta and [.501,1.999], with delta=1e-4*K0. Each
actual-limit check rebuilds delta and [.5,2], discarding stronger covariance
cuts and bounds. No secondary objective ran. A convex segment safeguard retains
true physical caps during construction; it performs no final projection.

| Actual-limit LP relaxation | Anchor A | Anchor B |
|---|---:|---:|
| Primal transformed slack | .0122934929333 | .0795105169555 |
| Maximum scaled primal violation | 1.57e-15 | 4.08e-16 |
| Stationarity infinity norm | 1.98e-13 | 1.48e-13 |
| Native/data allowance | 4.885e-7 | 1.0864e-6 |
| Independent lower bound after all allowances | **.0122930044338** | **.0795094305902** |
| Nonzero dual rows | 48 | 50 |

Independent long-double replay reconstructs supported necessary cuts from the
context with zero recorded discrepancy. A containing-box dual bound handles
imperfect stationarity and omitted bound multipliers. Primal–dual gaps before
allowances are below 1.2e-13. Rows, scales, multipliers, physical coordinates and
sided memberships are in
[active constraints](cls-cao-tamura-budget-h2i-evidence/active-constraints.json).

Each proof combines eight epsilon-bound rows, 38 PSD cuts, response caps and
GL16 divergence tangents. Mass distance is not an active dual constraint.
Different cut banks produce different relaxation bounds; these are not exact
optima of the full PSD/quadratic problem or measures of physical model quality.
The obstruction concerns the required positive margin, not an untested global
zero-margin joint problem.

## Native and sided findings

Two repeatable candidate checks were used. GL20 and GL24 are seen sets. Neither
diagnostic qualifies; unchanged gates prevent fresh sampling and launch.

| State / set | Inner minimum, m²/s² | Outer minimum, m²/s² |
|---|---:|---:|
| A, native | +1.8903444e-6 | −.001568161535 |
| A, GL20 | +3.8954776e-5 | −.001540121924 |
| A, GL24 | **+9.9854222e-7** | −.001559757203 |
| B-anchor, native | +2.1168094e-6 | −.01546239906 |
| B-anchor, GL20 | +3.9134980e-5 | −.01489670514 |
| B-anchor, GL24 | **+1.2096614e-6** | −.01506147933 |

The positive target is **1.8385755345e-6**, so both GL24 inner minima miss it.
A's native outer minimum/K0 is −.0852922; B-anchor's is −.840999. A's native
witness is chart (s,zHat)=(.1,.5400339062), geometry76; its GL24 outer witness
is upstream near physical (x,z)=(−.1404440566,.1997593610). B-anchor's native
witness is chart (0,.2), geometry62; its neighboring GL24 witness is physical
(.0001357570,.2300171163), geometry63. Physical and chart coordinates remain
separate in the receipts.

At the old protected matching witness (x,z)=(.23549580558569316,.0065):

| Side | A lambda_min | B-anchor lambda_min |
|---|---:|---:|
| H-minus, geometry35 | +1.8903444e-6 | +2.1168094e-6 |
| H-plus, geometry41 | **−7.2754537e-5** | **−7.2541313e-5** |

Actual and dependent epsilon agree there but retain separate contexts. Unlike
H2h, H2i permits velocity responses on both sides. No gradient averaging is
used. GL24's inner witness lies farther downstream at
(.2415516231,.0064897663). ExtraSeedAudit now keeps independent inner/outer
covariance witnesses and four scalar-ratio extrema, with full gradients,
tensors, eigenvalues/eigenvectors and provenance. Independent isotropic-plus-
linear tensor replay agrees within roundoff.

Both diagnostic fields pass scalar envelopes, cutoff budgets, prescribed and
endpoint data, weak continuity and response caps. GL24 velocity L2 changes are
.0506007/.0505970; component maxima .467711/.467694. Divergence L2 is essentially
its active .060308099422 bound. Native weak-constraint norms are below 3.6e-14,
while sampled local |div U| reaches **12.63/s**. These are distinct measures.
The defining cutoff-budget relative defect is 3.17e-16, which does not make the
actual transported k budget stationary.

A's GL24 epsilon/cold-epsilon spans [.802940101,1.998431486], time-scale ratios
[.25243,2.38321], and length-scale ratios [.17881,3.36280]. Dimensional data,
alpha=0 residuals and actual source/transport decompositions are in
[native evidence](cls-cao-tamura-budget-h2i-evidence/native.json). No physical
storage history or bound supply was created. H2h's terminal phase-I state was
not reused as an improved seed.

The read-only alpha=0 budgets also distinguish actual production from the cutoff
identity. For A, inner integrated P_actual is .00453604689 and E_C is
.00232620332, while p_A0 integrates to .00232620332. Outer k production/destruction
are .00643400749/.00426858812; the explicit outer outward transport is
.01849618763. These k-rate area integrals have units m⁴/s³ per unit span.
The constant-test outer-k residual is .01633076827, with constrained-row reactions
reported separately. Regional weak-balance identities close within 5.3e-17;
that algebraic reconciliation is not a steady balance or an unconstrained
evolution equation on prescribed rows. Zero diagnostic storage follows from
alpha=0 and does not describe a measured physical trajectory.

## Native accessibility

At six limiting/dual-supported witnesses per field, constrained native
consistent-mass Riesz responses measure retained-span accessibility while
preserving essential and independent continuity rows.
[Response evidence](cls-cao-tamura-budget-h2i-evidence/accessibility.json) gives:

| A witness | eta² |
|---|---:|
| GL24 inner limit, geometry35 | 8.10e-6 |
| Old H-minus match, geometry35 | 7.83e-6 |
| Old H-plus match, geometry41 | 3.72e-5 |
| Native outer minimum, geometry76 | .89168 |
| Roof dual point, geometry75 | .02548 |

The span has substantial response for A's bulk minimum and very little for
limiting interface/inner gradients. Along a single mass-normalized witness-
ascent ray at A's H-plus point, divergence permits about 2.79e-8, versus .00548
from the sampled velocity cap. Several inner/roof ascent rays immediately leave
the active divergence bound. These are local directional diagnostics, not an
additional search or global feasibility statement. No direction was added.

## Paired flat representation diagnostic

Archived H2d already showed analytical free loads near roundoff, improvement
under localized h/p resolution and an unresolved 1e-3 U/k target. Its best
2885-slot reduced diagnostic had residuals about
[.004798,.000691,.001275,.000131]. The six-panel transfer needs 4399–5339 slots.
The earlier single-knot transfer improved traction while worsening gradients
and residuals; it remains rejected. No h/p sweep or evolving knot case repeated.

One paired zero-step **same-space** diagnostic minimizes L2 error squared plus
H_A² times gradient error squared for free U/k. Native compatible prescribed
traces and k's interface trace are retained. W, complete transported epsilon
and both pressures remain unchanged. Existing assembly supplies GL16 matrices.
No balancing load, fit change or new hill reference is introduced.

| Model / projection | U residual | W residual | k residual | epsilon residual |
|---|---:|---:|---:|---:|
| CLS, retained L2 | .1342977 | .00472120 | .0476149 | .00178978 |
| CLS, H1 diagnostic | .0161768 | .00112508 | .00503692 | .00292909 |
| Standard, retained L2 | .1341183 | ~5.0e-16 | .0455783 | .00289101 |
| Standard, H1 diagnostic | .0161437 | .000085832 | .00481859 | .00282737 |

Common U-gradient L2 error falls .2599164→.1746809 and k-gradient error
.0148401→.00947156. U L2 improves; k L2 increases 1.00505e-6→1.17949e-6.
Epsilon coefficients/profile error are unchanged. Sampled traction error moves
+.12496%→−.001632%, yet U/k residuals still miss the target. The dominant remaining
U-gradient panel is .2–.8313 mm.

This scalar H1 projection is not a coupled continuity projection: regional
weak-continuity norms rise to 7.81e-5/2.23e-5 from roundoff baselines. It is
**not adopted**. Shared rows combine before mass-dual contraction; signed
panel/component projections preserve cancellation and reproduce assembled loads
within 8.2e-17. Component norms are not added as percentages. Analytical free
loads stay near 1e-13 or smaller. The actual analytical CLS Yap term remains a
model baseline; standard Yap is zero. Stress/production errors, sided fluxes and
signed panel data are in [flat evidence](cls-cao-tamura-budget-h2i-evidence/flat.json).

## Verification, cost and stopping gates

The unchanged H2g production solver SHA256 is
`844871bb5c70224658c88cca3f8088d7738e6f50a08ae521c74c7d02f4aec0a0`.
The isolated new probe SHA256 is
`f1179bca0b2d71f2c7f1be14209bbf8c20d7432c1c6d854aad4ffdad73891457`.
Flat projection used its original build path; executable and library-content
hashes match the separately frozen probe closure. Original paths are recorded.
No historical executable or library was replaced.

Focused Linux tests pass for exact tensor/LMI equivalence, signed coefficients,
Rayleigh/quadratic necessary cuts, known feasible/obstructed cases, buffered
versus actual limits and budget-exhaustion classification. Driver fault injection
covers every ordinary stage: process status, missing/malformed report, rejection,
wrong steps/checkpoints, missing equation counts and overwritten solved epsilon.
Each result is validated before any later invocation. Deliberate rollback is
distinct. H2i's fresh gate requires actual GL28 order, exact state/base/model/
runtime/basis binding, complete scalars, margin and physical caps. Existing H2g
native calibration and H2h mathematical regressions pass. Completed H2g native
restart/import/rollback evidence is reused, not presented as new hill execution.

| Linux work | Wall time | Peak RSS |
|---|---:|---:|
| Context export, driver-inclusive | 80.50 s | 1.00 GiB |
| Joint search, end to end | 20.72 s | 2.30 GiB |
| Optimizer inside that process | 19.41 s | included above |
| Repeatable native A / B checks | 34.06 / 34.14 s | 1.29 / 1.30 GiB |
| Native accessibility A / B | 7.68 / 7.73 s | not separately measured |
| Paired flat zero-step readback | 57.85 s | .716 GiB |

Native check times include setup, audits and output; they are not exclusive audit
timings. Proof replay and evidence I/O are recorded separately in
[cost evidence](cls-cao-tamura-budget-h2i-evidence/cost.json). Conservative work
accounting includes diagnostic replay beyond the 25 optimizer scans and stays
below 120 scans/600 s. Two of three permitted repeatable checks were used; fresh
audits and physical hill attempts/accepted steps are all zero.

Storage UUID `5AB8A669B8A642FD` is verified at `/media/jay/新加卷1`. Heavy context,
fields, raw samples, logs and runtimes remain under
`/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2i-20261115`.
Four build workers, one native rank and one library thread stayed within the
72-worker limit. See [artifacts](cls-cao-tamura-budget-h2i-evidence/artifact-index.json),
[commands](cls-cao-tamura-budget-h2i-evidence/commands.md),
[verification](cls-cao-tamura-budget-h2i-evidence/verification.json) and
[comparison manifest](cls-cao-tamura-budget-h2i-evidence/comparison-manifest.json).

| Gate | Observed result |
|---|---|
| Calibrations and full native response/nullspace contract | passed |
| Actual-limit sampled obstruction / independent dual replay | verified for both anchors |
| Full covariance and seen-grid positive margin | failed |
| Scalars, profiles, budget, boundary, continuity and caps | passed on diagnostics |
| Fresh GL28 / exact-runtime V2 certificate | unexecuted: prerequisites failed |
| Coarse8/fine16, 4+4 hill restart and curved rollback | unexecuted: no qualified seed |
| Common flat spatial target | open; H1 projection not adopted |
| Temporal accuracy, developed flow and held-out validation | unqualified; no long run |

**Next package:** construct a small native continuity-preserving velocity patch
response around downstream interface geometries35/41, with explicit
strong-divergence compensation and overlap removal against the retained pool.
Measured interface accessibility supports this focused change before another
scalar or height study. Keep physical bounds and calibrations. Common flat
spatial qualification must precede quantitative development, temporal/spatial
sensitivity and held-out experimental claims. EP-AWF and unrelated models remain
deferred.
