# CT-Accepted-5: localized flat-reference traction and observation coordinates

CT5 tests one radial split with unchanged CT4 physical inputs. The primary
A95/6 mm/P4-P2 contrast resolves the ground-traction accuracy failure through
T=.02. All three models pass the transferred recipe through T=.02. A93 uses
the one conditional fresh .0005 history because its original .002/.001
interface-traction comparison is ambiguous at T=.008. Complete final model
decisions and primitive bodies are in the adjacent evidence directory.

The independently corrected observation audit passes. It changes 40 archived
crest distances, adds 36 Fig.7 observations and leaves the rough 7 mm Fig.12
traverse as the limiting known turbulence height. These are corrections to
CT4's height audit, not changes to its frozen observations, fits, numerical
results or failed accuracy decisions. See the explicit
[correction notice](cao-tamura-accepted-5-evidence/CT4-observation-correction.md).

## Contract, source and unchanged physics

The starting revision was `b9882937a5096aa9c73219dd3b6288cd8a346983`, clean in
the Mac implementation worktree and identical to the fetched remote. There
were no intervening commits. The [CT5 contract](cao-tamura-accepted-5-evidence/contract.json)
was committed before new evolution. CT1–CT4 histories, failed preparations,
installations and runtimes remain untouched.

The governing executable and libraries are exactly the frozen CT4 build-v3
runtime, built from `123f31c974e0b75e551e43a5203749a89378fade`. No native C++
equation, kernel, coefficient, boundary policy, projection rule or acceptance
location changes in CT5. Source changes implement read-only load/traction
analysis, exact interface-side replay, observation geometry and bounded
orchestration. The default replay request remains byte-identical in the
constitutive and physical-budget controls; an invalid side inventory rejects
explicitly.

Natural A93/A95/E94E95, model-specific q/kappa and wall calibration, finite
inner viscosity, zero outer molecular diffusion, transported nonconstant
inner k, dependent total epsilon and the Accepted cutoff remain unchanged.
NeutralCompositeV1 and ConstrainedRadialV1 supply the same continuous input.
The one-sweep method retains its stage order, signed production, full 3D
covariance, absolute pressure, free total-traction outlet and current endpoint
K dependence. There is no reference cancellation/restoration, clipping,
source correction, post-k epsilon repair or imposed evolving wall gradient.

For the primary case, u_tau=.07426793792441279, nu=1/10500,
K=.021081109428869052, B=.94768443074701 and h_s=.11096757811954916 Hh
(4.438703124781967 mm). H_A=.15 Hh=6 mm; Hh=40 mm. The streamwise partition
is [0,.125,.25] and the compatible top is z=2 Hh=80 mm. The fit memberships,
weights, normalization and all continuous inlet/top/outlet functions remain
those of CT4.

## Load reconciliation and the one-split decision

Saved CT4 A95/P4 coefficients at 0,.001,.002,.004,.008 are evaluated read-only.
The native integration geometry, basis rows, generating owners, quadrature
coordinates and weights are retained. Independently reconstructed external
advection and prescribed total-traction loads are included with the volume
terms before comparison with the unmasked native residual. Shared interface
rows enter once. Constraint replacement and prescribed-row reactions remain
separate from physical traction.

At all five states the physical-term reconstruction agrees within 5.56e-16
in the maximum row entry, against the predeclared 2e-10 tolerance. Individual
stress, pressure, advection, diffusion, sources, roughness and open loads are
retained in the evidence. Independent continuous-reference weak loads are
evaluated on the same native rows; this diagnostic never enters evolution.
The continuous weak identity passes independently of the nonzero residual of
the represented reference.

Localization uses the signed contribution
`delta_r_slab^T M_free^-1 delta_r_total`, with the same free velocity mass and
normalization for every slab. Negative shares and cancellations are retained;
isolated element stress norms are not ranked as residual errors. The first
slab [0,.0046875] Hh = [0,.1875] mm contributes 2.3664087189e-7, or
62.6502% of the initialized global dual-norm squared. It is the largest
positive signed share. Early-stage traction changes support a gradient
mechanism, so the predeclared first-slab hypothesis is selected before any
new trajectory is constructed.

Exactly one point is added across both columns: H_A/64=.00234375 Hh=.09375 mm.
All original nodes and physical breaks, the outer grading, streamwise mesh,
P4/P2 spaces and quadrature policy remain fixed. Primary native construction confirms
**2957 stored slots**, from 2839 (+118), below the unchanged 3000 guard.
Every space has fresh projections and native zero-step histories. No evolved
CT4 vector is transferred. Same-space timestep controls have identical actual
initialized coefficient hashes. Their represented boundary traces differ from
the old space, while the continuous physical input hashes agree.

## Stage-level physical traction attribution

At each sampled ground point CT5 evaluates
`tau=D(K_H)*(U_z+w_x)` from that stage's actual coefficients and current
polynomial endpoint K trace. A separately labelled frozen-input-K evaluation
is counterfactual only. Full one-sided interface traction includes pressure
and the selected model's complete Reynolds stress.

The exact pointwise decomposition is integrated afterward:
`Delta tau = D0*Delta S + S0*Delta D + Delta D*Delta S`.
For the first CT4 fine step (T=.001), mean traction changes from
.00551572660354443 to .005511260836032797. Its contributions are:

| Contribution | Mean traction change |
|---|---:|
| Gradient, D0 Delta S | -4.46586715129e-6 |
| Endpoint-dependent viscosity, S0 Delta D | +9.97220855157e-11 |
| Cross term | -8.24271966346e-14 |
| Exact summed/actual change | -4.46576751163e-6 |

The decomposition defect is below 1e-18. The momentum predictor gives most of
the first change; pressure correction is smaller and initially opposes it.
Post-epsilon leaves this traction unchanged, and the post-k viscosity update
is tiny. The selected-stage gradient attribution magnitude is 6.80013e-6,
versus 6.03803e-10 for viscosity/cross contributions. This supports resolving
the near-ground gradient; it does not justify freezing K, pressure or flow.
Means, spatial variation, extrema, actual state/history hashes, stage pressure
actions and full interface profiles are retained. Reconstructed outlet
mismatch is a separate observation from the prescribed natural datum.

## Accuracy and uncertainty rules

The inherited ground mean and full interface RMS limits remain 1e-3 on
u_tau². Field reference RMS limits are 1e-3 on the CT4 scales, with vertical
velocity limited to 1e-4. Mass, actual constitutive/covariance/cutoff checks,
source identities and complete scalar budgets remain mandatory. Budgets
require all four controls (native, q48f16, q64f16, q64f32), all four equations
(inner_k, outer_k, outer_epsilon, whole_k), and both reconstructions (k,
outer_epsilon), bound to actual state, geometry and consumed boundaries.

Traction uncertainty uses signed mean differences and the actual ground and
two-sided full interface profiles. Let s be their largest normalized
timestep spread. Admission requires both errors <=1e-3 and
`max(error)+2*s<=1e-3`; a remaining overlap is ambiguous. A spatial improvement
must exceed `4*max(baseline spread, refined spread,2e-7)`. These rules were
fixed before the new runs. Differences of absolute error norms are not used
as a substitute for signed/profile comparisons.

At T=.008 A95's fine-step ground error falls from .00146758984 to
.000163097248. The spatial signal is .00130449259, compared with an uncertainty
limit .0000538356294. The new paired spread is 3.28575e-6 and the robust
ground-error upper bound is .000169668751. Interface traction also passes.
At T=.02 the ground error is .000163311838 and the paired spread is
6.95719e-7. This resolves the particular ground-traction failure with fixed
physical input. It does not establish a formal convergence order or full
spatial convergence of every equation.

The original A93 .002/.001 T=.008 pair has ground errors
.000396914496/.000397469139 and interface errors .000902978380/.000940037404.
Although both endpoint interface errors are below the limit, s=3.51038e-5
makes its upper bound .00101024509: **ambiguous, not admitted**. That receipt
is sealed and retained. It triggers the one permitted fresh .0005 history,
with the same initialized coefficients and continuous input. This is the only
additional timestep control. Continuation uses unchanged-method histories
only after the recomputed .001/.0005 gate passes.

The primary A95 comparison isolates the single split at P4/P2. A93/E94E95
transfer the supported P4/P2 plus split recipe from their frozen physical
cases; CT4's existing six-mm primary spaces for those models were P3/P1.
Their bounded outcomes therefore do not independently isolate the split
from the polynomial-degree change. Initial representation differences,
subsequent drift and continuous-reference errors are reported separately.
The historical CT4 pressure-interpretation uncertainty remains unresolved.

## Completed native pairs

All members below complete the stated physical horizon with passing domain,
field, mass, complete physical-budget and constitutive checks. The listed
errors use the finer timestep; the gate also checks the other member and the
signed/profile spread. All cases have 2957 slots.

| Model | T | Admitted dt pair | Ground error | Full interface RMS | Traction spread s |
|---|---:|---|---:|---:|---:|
| A95 | .008 | .002/.001 | 1.63097e-4 | 8.18157e-5 | 3.28575e-6 |
| A95 | .02 | .002/.001 | 1.63312e-4 | 8.29227e-5 | 6.95719e-7 |
| A93 | .008 | .001/.0005 | 3.97603e-4 | 9.56062e-4 | 9.25325e-6 |
| A93 | .02 | .001/.0005 | 3.98444e-4 | 9.80526e-4 | 4.06108e-6 |
| E94E95 | .008 | .002/.001 | 2.20232e-4 | 1.76361e-4 | 4.01013e-6 |
| E94E95 | .02 | .002/.001 | 2.20610e-4 | 1.76665e-4 | 1.00716e-6 |

A93 is limiting: its T=.02 robust interface upper bound is .000992850582,
only 7.15e-6 below the unchanged .001 limit. This is a bounded pass, not a
large safety margin or permission to extend duration. The original ambiguous
A93 pair remains in its own comparison body. Its .002 member is not extended.
The admitted .001 member continues only after the independently initialized
.0005 companion resolves the short gate.

All 13 distinct new evolution segments retain every accepted-state traction
profile and report no first accuracy crossing. Restart-boundary samples are
retained as segment observations; they do not count as extra evolution steps.
There are no new physical rejections, bound calls, restoration/correction
solves or extra governing solves. Final checkpoint, state, initial/reference,
consumed explicit/implicit boundary and actual mapped-library identities are
checked against the retained native reports.

At T=.02 the fine-step maximum scaled field reference error is 6.164e-6
(across the five fields and three models); vertical RMS is at most 3.770e-8.
Fine-step mass defects are at most 1.43e-16. The largest fine-step scaled
four-control scalar defect is 3.85181e-6, against .02. These physical budgets
retain backward storage, frozen-stage sources, actual advective/diffusive
face fluxes, endpoint-dependent epsilon storage, coefficient splitting and
separate interface contributions. Algebraic row reactions do not replace
physical fluxes. Reconstruction agreement and physical accuracy are separate
statuses.

Native and independent physical sampling both pass. At A95 T=.02 the actual
minimum radicand is about .69194; sampled full covariance remains positive.
The bodies retain actual raw gradients, signed divergence, production and
witness sides/geometry. Weak continuity is not pointwise incompressibility:
at T=.008 on the common physical inventory, the old/new divergence RMS is
1.08689e-5/9.00958e-6 and maximum is 6.38705e-4/4.18693e-4. Full production
continues to include isotropic divergence work.

The common-grid A95 initial U difference is 3.66682e-7; at .008 it is
2.22775e-7. Initial k difference is 2.28987e-7 on K and evolved is
2.20454e-7. The corresponding absolute-Pi difference starts at zero and is
3.64533e-8 at .008. The archived .001, .002 and .004 contrasts separate this
initial representation from subsequent drift. The CT4 baseline is compared
only at common accepted times through .008; its failed trajectory is not
extended to provide a .02 comparator.

## Corrected observation audit

For h(t)=cos²(pi*t/5) on [-2.5,2.5], half the squared-distance derivative is
`F=t-x-h'(t)*(z-h(t))`. With a=pi/5 and c=cos(2at),
`F'=1+a²+a²*(2z-1)c-2a²*c²`. Analytic critical points partition F into
monotone intervals. Exact/sample-node, endpoint and tangent roots are included,
deduplicated with explicit roundoff tolerances, and compared with foot and
flat candidates. No physical coordinate or distance is snapped/clamped.
Near-coalescent binary64 multiplicities are explicitly delimited.

Analytic crest/flat/foot controls, reflection, crest continuity, tangential
roots, competing candidates and the vertical-distance upper bound pass,
alongside independent distance minimization. The old crest witness is
reproduced from the literal reviewed function before verifying its correction.
Audit coverage increases from 760 to 796. Eighty lower/upper extraction-envelope
distances also change and are listed separately. Forty prior nominal distances
change; four smooth observations per height switch from outer to inner:
point0 of f4a_SS_x0, f5a_SS_x0, f8_SS_x0 and f10_SS_x0.

All 36 visible Fig.7 observations retain their RMS/shear quantities, roles,
values, weights and extraction units. Coordinates and extraction envelopes
are converted with delta/Hh=6.25. The rough lower bounds are 9.375017697 mm
(RMS ratio) and 9.334900626 mm (shear ratio). The known Fig.12 rough 7 mm
traverse remains limiting. Fig.11 heights remain unavailable where delta_7.5
is unspecified. Extraction envelopes are not invented experimental or probe
uncertainty. Nothing is refitted and no hill calculation is performed.

CT4's `nu_over_nut_H` is the inner auxiliary ratio. At A95/6 mm it is
.0207671722647; the separately reported outer `nu/(kappa*u_tau*H_A)` is
.0211578311468. Neither changes the stipulated high-Re outer equations.

## Scope retained

CT5 is a bounded prescribed-reference flat-patch accuracy experiment, not a
developed precursor or an experimental validation. Mathematical reference
admissibility, native accuracy, observation geometry and future terrain
admission are distinct. The .1 thin-layer flag (6 mm curvature ratio .118435),
roughness-sublayer adequacy and the interface only 1 mm above 5 mm cubes remain
terrain qualifications. The smooth reference/high-Re/low-observation issues
remain unresolved. EP, 5.5/6.5 mm native repetition, T=.2, hill predictions,
native curved production, caching and backend work are not attempted here.

## Cost, resources and provenance

Mac source/Git work uses the existing isolated implementation worktree.
Authoritative runs use
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-5-20261112` on
`jay-Super-Server`; the actual mounted device is `/dev/sda2`, ntfs3. The
pre-existing relocation bridge and other campaigns are preserved. The source
clone, temporary files, native outputs and evidence are isolated under that
root. No shared installation or old runtime is rebuilt or replaced.

The frozen executable SHA256 is
`9f78fa07ed1eab1f1463bc3a7b546a039fb27d8195c757bd6577590c41bd6496`;
RANSPreservationProbe is
`9291343c45ea984eaa04132e5c9df52bc8b9662e772127a2768bf4bace6c4fc0`.
The reused Release build uses `/usr/bin/c++`, `-O3 -DNDEBUG`, GNU++17 and
MPI-enabled libraries; every CT5 native calculation is serial, with numerical
libraries explicitly one-threaded. This is not MPI qualification. Actual
mapped library paths and bytes, build flags, original source hashes and input
identities are in the sealed provenance, not inferred from a Git commit.

New evolution uses the CT5 driver committed at `b6a5c3d`; the recovered
controller runs at `3d11328`, and the single A93 half-step at `22deba8`.
Final read-only collection is at
`acafc765f24392446a2fea0c4a66fd458152d1bb`. A later orchestration-only guard
(`f8473f0`) enforces the single conditional half-step across all model
transfers; its Linux control rejects a second initialization and verifies
that the physical/classification functions are AST-identical to the qualified
source. It does not resume or regenerate a trajectory. All use the same
native runtime.
The report/evidence commit is separate from these tested source identities.

Thirteen native evolution segments perform 124 accepted steps. Their measured
process CPU totals 731.590 s, elapsed totals 739.191 s and peak native RSS is
443404 KiB (433 MiB). Instrumented exclusive accepted-step CPU is about
5.14–5.34 s/step. For the A95 fine-step .008–.02 segment, 62.784 exclusive CPU
seconds comprise 39.659 domain/constitutive, 2.091 assembly, 1.062 factorization,
1.709 mass/pressure, .906 solves/linear validation, 3.150 snapshot/checkpoint
I/O, 8.993 physical audits and 5.215 remaining orchestration/constraints.
Nested inclusive totals are not added. Initialization and independent
read-only probes are outside these accepted-step figures. No optimization is
implemented.

Every new numerical/bookkeeping family uses the owned 16-worker launcher
with the 72 aggregate guard and one-thread libraries. Completed native,
diagnostic and archival monitoring passes: project peak12, aggregate peak36.
The aggregate peak is at 2026-11-12 15:11:00.794621 UTC. Peak PIDs, parentage,
thread counts, ownership and timestamps are preserved in
`resource-highwater-v1.json`. Allocation wait and recorded suspension events
are zero. Process CPU, elapsed time and their scopes are separate in each
receipt. All archival monitor receipts are closed; the collector's earlier
pending-finalizer snapshot is superseded by completed resource receipts and
the transfer verification. CT1–CT4 resource incidents remain unchanged.

Four technical preparation failures are retained separately from physical
results: the first diagnosis completed five valid packets before a list-abs
reporting error; the first controller completed both A95 short trajectories
before duplicate replay-side arguments stopped analysis; an intended side
negative control raised ValueError outside the original narrower catch; and
one bundle fetch requested a nonexistent branch label instead of its HEAD.
Each recovery uses a new destination and validates any reused sealed packets,
states and histories. No failed physical history is rewritten or restarted.
The original A93 coarse-pair ambiguity also remains explicit. Missing epsilon
budget coverage, missing timestep member, unequal endpoints, forged traction
summaries and wrong side lengths are discriminating negative controls.

## Reproduction and evidence map

The [final decision](cao-tamura-accepted-5-evidence/qualified-v1/final-decision.json)
binds the complete comparison bodies, pressure/stage inventory and completed
resource receipts. See the adjacent
[evidence README](cao-tamura-accepted-5-evidence/README.md) for file identities,
seal checks and exact reproduction commands. Large basis packets, native
coefficients, PIMPLE snapshots, fields and checkpoints remain on Linux.

Reused evidence is limited to CT4's unchanged physical fits/mathematical
references, frozen runtime, existing primary A95 states/stages and frozen
observations. All seven CT5 zero-step histories and all new evolution segments
are newly generated; the one .0005 initialization matches its same-space
companions exactly. Five native diagnostic packets, seven initial probes,
endpoint constitutive/budget controls and the common-grid baseline export are
new read-only evaluations. CT4 trajectories are not rerun or extended.

## Decision and next gate

Implemented and passed: closest-graph/observation correction; complete
volume/open-load reconstruction; signed localization; exact stage traction
decomposition; native one-split count/initialization; A95 representative
accuracy and traction uncertainty; model-specific six-mm transferred recipes
through .02; complete physical budgets, domain checks and resource compliance.

Retained failed/ambiguous evidence: all CT4 accuracy failures and prior
qualifications, the A93 .002/.001 short uncertainty gate, and technical
preparations listed above. Still unresolved: CT4 pressure interpretation,
general spatial convergence, physical roughness-sublayer/interface adequacy,
the .1 thin-layer issue and smooth preparation. Unattempted: native 5.5/6.5 mm
repeat, longer duration, developed precursor, experimental hill prediction,
native curved production, EP and backend/caching changes.

The next implementation package should be a **conservative curved Accepted
manufactured patch**, with independent metric weak balances, one-sided full
traction, transported-k/interface conservation, dependent epsilon and an
explicit zero-curvature regression to this flat case. Keep native curved
production guarded until that verification passes. It must retain the
roughness-sublayer and thin-layer questions as physical terrain gates; this
flat accuracy result does not answer them or validate the hill experiment.
