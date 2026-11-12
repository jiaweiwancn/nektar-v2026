# CT-Accepted-3: bounded evolved space and degree qualification

A95 and E94E95 complete the bounded flat-pilot contract through T=.2.
A93 passes T=.02 and then exits on the original full-tensor covariance gate.
The A95 evolved spatial/degree contrasts remain **unresolved**, despite small
global field differences. This does not admit a developed precursor or hill
prediction. The complete [decision](cao-tamura-accepted-3-evidence/decision.json)
keeps physical completion, spatial uncertainty, implementation diagnostics and
resource compliance separate.

This package starts clean at `db37c6b11882e7a06413a29eb6623ab212c4d4a5`,
identical to the fetched remote. The unrelated main Mac worktree and all CT1/CT2
evidence are preserved. The pre-run
[contract](cao-tamura-accepted-3-evidence/contract.json) fixes four new cases,
two timesteps, three conditional horizons and separate physical, spatial,
diagnostic and resource decisions. The
[spatial interpretation contract](cao-tamura-accepted-3-evidence/spatial-interpretation-contract.json)
was recorded before computing cross-space comparisons.

Authoritative computation uses the actual mounted disk at
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-3-20261112`.
The user-selected CT2 prefix-only relocation bridge remains read-only and
hash checked. Source, temporary files, outputs and the optional timing build
are isolated. The shared installation and existing CLS/WJ projects are untouched.

## Retained scientific contract

The original Natural A93, A95 and E94E95 high-Re outer equations, zero outer
molecular diffusion, finite inner viscosity, Accepted budget cutoff,
nonconstant transported inner k and dependent total epsilon are unchanged.
There is no inner epsilon equation, forcing, reference restoration, clipping,
regularization, bound projection, coefficient/input retuning or post-k trace
repair. The native layer remains H_A=.4 (16 mm).

AcceptedOpenFlatV1 retains its prescribed inlet/compatible top, full-total-
traction free outlet, scalar advection and diffusive data, backflow restriction,
and absolute traction-determined pressure level. No mean-pressure gauge is
reintroduced. Consistent-free-mass one-sweep PIMPLE, stage ordering,
GeneratedVolumeOwnerV1, unsupported-route guards and the 3000-slot limit remain.
Each model's physical FUNCTION strings are copied from its sealed CT1 source;
epsilon's Hermite function is never regenerated from a new mesh.

The CT2 Zminus P3/P1 A95 histories are reused within their scope through exact
executable/library/input identities. Every new model/space/dt has a fresh
zero-step native history. Within a space the two dt initial vectors must be
identical. Different spaces retain their own projections.

## Admission and independent degree diagnostics

The physical-budget classifier now requires exactly four controls
(`native`, `q48f16`, `q64f16`, `q64f32`), all four equations
(`inner_k`, `outer_k`, `outer_epsilon`, `whole_k`) and both reconstruction
records (`k`, `outer_epsilon`). It checks finite primitives, positive scales,
recomputed storage-source-flux defects, reconstruction and quadrature
uncertainty. A paired admission also requires both dt members, the actual
endpoint, matching state/configuration identities and identical initial vectors.
Labels and a collector exit code do not authorize continuation.

All six actual CT2 budget records (both timesteps at .008/.02/.2) still pass.
Negative controls remove outer epsilon, a reconstruction, a quadrature/face
control, a timestep member or a required time; supply a nonfinite primitive,
nonpositive scale or foreign state/configuration; and forge a passing label
with a failed recomputed defect. They are rejected for their intended reasons.
An intentionally incomplete command returns exit 1. The
[coverage controls](cao-tamura-accepted-3-evidence/complete-budget-negative-controls.json)
also exercise reconstruction normalization and stale state bindings. CT2 is
not relabelled as a historical budget failure.

The pressure oracle independently constructs continuous nodal P1/P2 radial
mass projections with the actual default eight-point Functions projection.
It reconstructs tensor polynomials for pressure action. At every rough-volume
point it selects the containing outer geometry at that point's actual x and H,
with declared outer side; it does not reuse one face quadrature row.
Manufactured p=x^2 makes its derivative 2x discriminating. Full pressure action,
mixed rank, conservative weak/matrix rows and absolute-traction offsets are
checked on a small P4/P2 fixture before physical evolution.

The 319-slot P4/P2 fixture has 318 physical unknowns; its unused zero gauge
slot is outside the solve. Pressure action, divergence and consistent-mass
PIMPLE pressure operators each have rank 30. The implemented mixed system
has rank 142 with the intended constraints and no extra pressure null mode.
This is a finite-fixture check, not an inf-sup theorem. The independently
reconstructed raw weak defect is 4.90e-16 and the largest frozen-matrix defect
is 1.12e-14. With nonzero roughness, p=x^2 gives pressure-action defect
9.31e-16, derivative defect 3.11e-14 and a wrong single-row derivative defect
1.96. Absolute traction shift c=.125 gives action defect 3.61e-16 and
two-step coefficient defect 6.21e-14. Full primitives and generating geometry
are in [degree controls](cao-tamura-accepted-3-evidence/degree-pressure-and-coverage-controls.json).

The primary solver and observer remain the pinned CT2 governing runtime.
Its old Functions-projection metadata label is not used as mathematical evidence:
the independent native projection checks establish the actual rule. The
separate timing runtime contains optional telemetry only and is qualified by
bitwise accepted-state/native-history equivalence before its profile is used.

## Numerical results

Every row uses independent dt=.002/.001 histories. Each pair has identical
actual zero-step coefficients, effective initial fields and boundary maps.
The stored-slot counts below are native confirmations, including the unused
zero gauge slot. Exact block counts, mesh knots, basis/point configurations,
projection identities and initial/reference/boundary hashes are in the four
`*-initial.json` receipts. The mesh is fixed throughout each trajectory.

| Model / mesh / space | Slots | T=.008 | T=.02 | T=.2 |
|---|---:|---|---|---|
| A95 / CT2 Zminus / P3/P1 (reused) | 1172 | pass | pass | pass |
| A95 / XZminus / P3/P1 | 1669 | pass | pass | pass |
| A95 / XZminus / P4/P2 | 2981 | pass | pass | pass |
| A93 / Zminus / P3/P1 | 1172 | pass | pass | rejected before endpoint |
| E94E95 / Zminus / P3/P1 | 1172 | pass | pass | pass |

XZminus adds x=.1875 and z=5.125 to CT1. Zminus adds only z=5.125.
P4 uses compatible pressure degree two. Native basis points are 8x8;
wall-normal/tangential/mass/outer-volume/interface orders are 32/8/10/10/10.
Frozen controls retain 48/12/12/12/12 with face16, 64/16/16/16/16 with face16,
and that latter volume rule with face32. They evaluate the same actual
coefficient state read-only. The original cutoff panels and physical functions
are retained.

Every passing milestone has all four budget controls and equations, independent
constitutive comparisons, native and dense physical domain inventories,
physical-face mass flux, source/storage and dependent-trace checks. All 22 new
successful endpoint members pass the physical-budget primitives. At the
last paired qualifying horizon, selected fine-dt results are:

| Case / T / dt=.001 | Scaled mass defect | inner-k budget | outer-k budget | outer-epsilon budget | whole-k budget |
|---|---:|---:|---:|---:|---:|
| A95 XZ P3 / .2 | 5.94e-16 | 5.36e-6 | 5.54e-6 | 3.65e-5 | 4.11e-6 |
| A95 XZ P4 / .2 | 2.14e-16 | 1.51e-6 | 4.63e-6 | 1.48e-5 | 1.87e-6 |
| A93 Zminus / .02 | 8.60e-17 | 1.80e-7 | 7.21e-6 | 3.00e-5 | 1.44e-6 |
| E94E95 Zminus / .2 | 3.30e-16 | 1.32e-7 | 8.10e-6 | 3.94e-5 | 2.95e-5 |

Budget entries are the maximum absolute physically scaled backward-storage /
frozen-source / independently reconstructed physical-flux defect across the
four controls. The limits remain mass 1e-8 and budget .02, with reconstruction
2e-10. These finite-defect passes do not establish convergence. Algebraic row
reactions remain distinct from independently evaluated face fluxes. The
conservative scalar budget has no added f*div(U) correction.

See [physical summaries](cao-tamura-accepted-3-evidence/physical-summary.json)
and the complete `*-T*.json` and `*-derived.json` bodies for both timesteps,
quadrature/face spreads, outlet traction, interface fluxes, trace lag and
stage diagnostics. [Accepted-state extrema](cao-tamura-accepted-3-evidence/accepted-state-extrema.json)
retain all native accepted-state times, component extrema with time/geometry/
side witnesses and final audits separately. Native, fixed physical and common
comparison inventories are separately visible; the inventory with the largest
positive margin is not selected as a substitute for the others.
These are segment-wise inventories; a repeated restart-boundary observation is
not another accepted step. Predictor invariants remain diagnostic: an undefined
intermediate tensor does not introduce a rejection point absent from the
production algorithm.

The first A95 P3 coarse .008 trajectory succeeded, then collection failed on
an empty constitutive-sample request. The adapter-only fix and hash-bound
recovery reuse that exact successful native trajectory and run the missing
fine member. No physical trajectory was replaced or reset. A later compact
collector initially assumed every failed A93 high-quadrature probe supplied
budget primitives; the corrected export retains those failed packets by hash
and supplies no fabricated pass. Both failed collection receipts remain on Linux.

The initial builder also copied transient CT1 dt/step/target/XML-hash template
fields into construction metadata. Actual consumed XML, native dt/time/history,
immutable requests and FUNCTION maps were correct and independently verified.
The writer now records actual execution metadata; existing artifacts are not
edited. Each derived receipt exposes the original metadata and actual XML/report
values explicitly. Its writer-only regression changes no numerical input.

## Evolved spatial and degree sensitivity

The comparison uses one physical grid with six Gauss points per cell of the
union partition, dense delta/inlet probes, the original outlet witness and
explicit inner/outer H sides. Coordinates, weights, region/side identity,
endpoint times and the absolute traction pressure level match. Six points
integrate squared represented P4 fields on each union cell adequately; four
points are retained only in the immutable original P3 packets.
All five compared initial/accepted-state inventories pass the additional common
physical domain checks. This is finite-inventory evidence, not global
realizability or a convergence theorem.

At T=.2 and dt=.001 the scaled physical weighted-RMS differences are:

| Contrast | U | w | k | outer epsilon | absolute Pi |
|---|---:|---:|---:|---:|---:|
| CT2 Zminus P3 -> XZminus P3 | 1.18999e-5 | 2.65796e-6 | 8.22063e-5 | 1.98152e-6 | 6.17751e-6 |
| XZminus P3 -> XZminus P4 | 3.07840e-5 | 3.11548e-5 | 8.21194e-5 | 9.77632e-5 | 3.94992e-5 |

Scales are U(Hh), K0, epsilon(H_A) and U(Hh)^2, with U(Hh)=1.
Both contrasts satisfy the declared 1e-3 global difference limit.
The first contrast is streamwise refinement at fixed degree and vertical
partition. The second changes both the velocity/scalar and compatible pressure
degree on the same partition. Neither is a formal convergence-order estimate.

Initial U/k/epsilon differences are at roundoff scale; initial w is zero.
Streamwise-only initial Pi agrees at roundoff. The P1 -> P2 pressure
projection instead changes initial Pi by scaled RMS 3.75936e-5.
This initial representation difference is retained separately from evolution.
At the original witness the independently verified initial pressure gradients
are +5.00974753e-6 (P1) and -7.38458043e-7 (P2).
The fixed normal-test projection load changes from -2.63745196e-8 to
+4.79086103e-10. Independent 48/96-point integration differs by about 3.1e-16.
This supports improved initial pressure representation; it does not prove a
stationary represented reference or bounded evolved local sensitivity.

The local A95 interpretation fails to establish that spatial changes are small
relative to the actual margins:

| T=.2 contrast | Maximum radicand difference | Smallest margin across both spaces | Difference/margin | Maximum dt radicand difference |
|---|---:|---:|---:|---:|
| Streamwise P3 | .309315 | .0389589 | 7.9395 | .0484756 |
| Degree on XZminus | .453403 | .0664066 | 6.8277 | .0484756 |

The spatial limit is one quarter of the margin. Both contrasts exceed it
substantially. Their largest radicand differences are near the outlet just
above delta, at (x,z)=(.25,6.25000001). The predeclared local dt-uncertainty
limit is max(1e-10,.1*spatial radicand difference); it is also exceeded at .2.
The two limitations are reported separately in the derived classifier rather
than allowing one status to hide the other.

For global field interpretation the limit is
max(1e-7,.1*spatial scaled RMS difference). At .2 the streamwise contrast's
epsilon difference is smaller than its timestep uncertainty, so that observable
remains unresolved; its other fields pass this interpretation check.
All five fields of the degree contrast pass at .008/.02/.2, while local
radicand sensitivity remains unresolved. Earlier common times and their
additional uncertainty limitations remain in the complete receipt.

The independent equal-time dt comparisons pass the original 1e-3 physical-pilot
limit. At .2 their largest scaled RMS field differences are 2.46223e-6 for
XZminus P3, 2.49135e-6 for XZminus P4 and 2.72685e-6 for E94E95 Zminus P3
(each outer epsilon). Passing that absolute pilot limit does not settle the
smaller spatial-difference or local-radicand uncertainty questions.
No temporal order is inferred from two timesteps, and no unequal terminal
times are compared as a spatial error.

## A93's retained covariance failure

Both A93 members pass .008/.02, then stop at the original covariance gate.
This is a forward-defined A93 tensor with a negative covariance eigenvalue,
not an A95 square-root failure or a nonpositive scalar. The failing direction
is the spanwise unit vector: omitting Ryy in a two-dimensional calculation
would miss the failure.

| dt | Last accepted T | Attempted T | Native witness (x,z) | Ryy = minimum eigenvalue |
|---|---:|---:|---|---:|
| .002 | .062 | .064 | (.125, 5.125) | -1.16190220e-5 |
| .001 | .053 | .054 | (.122518116, 5.130961974) | -4.52263963e-5 |

Both witnesses belong to recorded native geometry 24, outer region, radial
side 0. The coarse point lies on element boundaries; explicit geometry is
therefore used for its one-sided derivatives. The fine point is interior.
Independent original-A93 stress evaluation using the native raw gradient and
native audit k/epsilon agrees with all nine native tensor entries within
6.6e-19. Native-audit epsilon minus the raw rejected-field epsilon is zero at
both points, so this rejection is not an ambiguity in the refreshed H trace.

At the .002 witness, k=.00184529615 and epsilon=5.54090888e-6 remain positive.
Signed divergence is .0235771195 and full production is -3.63592392e-5;
the isotropic divergence-work contribution is -2.90045119e-5. At .001 the
corresponding values are k=.00182692975, epsilon=5.44702602e-6,
divergence=.0221973417, full production=-3.40405007e-5 and isotropic work
-2.70353226e-5. Signed production is retained.

The same physical witness is followed through the actual saved stages:

| dt | Last/input covariance minimum | Predictor | Corrected flow | After epsilon | After k/rejected |
|---|---:|---:|---:|---:|---:|
| .002 | 5.99643e-5 | -1.28222e-5 | -1.17189e-5 | -1.15442e-5 | -1.16190e-5 |
| .001 | 1.58212e-6 | -4.54754e-5 | -4.52906e-5 | -4.51787e-5 | -4.52264e-5 |

Loss is already visible in the momentum predictor at this retained point;
pressure correction partly improves it. Scalar transport changes the margin
slightly and does not restore eligibility. This is a diagnostic stage
comparison, not an added predictor rejection rule or proof of the earliest
global continuum failure. One unchanged-method attempt from each retained
last accepted checkpoint reproduces the failure and rejected coefficients
bitwise, with complete accepted-state/native-history rollback.

No source/action inconsistency or ownership defect is inferred from this
forward-defined negative tensor. Approximation, one-sweep coupling and the
inferred input/short-domain choice remain separate possible contributors.
There is no additional A93 mesh/degree search or constitutive modification.

The fine member's retained last accepted .053 state also fails the independent
q64f16/q64f32 audit inventory. Those probes are recorded as failed physical
audits, with no inferred budget-accuracy pass. This additional finite-inventory
observation reinforces the stop; it is not a replacement for the reproduced
native .054 rejection or a new production rejection rule.

## Accepted evolution costs and compatibility

Timing is measured on admitted .008 -> .02 segments with ordinary physical
checks and snapshots enabled. No startup/two-attempt cost fraction is reused.
The ordinary frozen runtime's broad scopes separate accepted evolution from
initialization, independent observation export and rejected attempts. The
[cost receipt](cao-tamura-accepted-3-evidence/accepted-costs.json) contains the
complete exclusive scope totals, per-case zero-step initialization and separate
observation-export CPU/RSS. Snapshot exports can use more memory than the solver;
the solver RSS in this table is not a bound for postprocessing.

Zero-step native preparation and independent observation export are separate
processes. Representative dt=.002 initialization costs are:

| Space | Native zero-step CPU s | Native elapsed s | Observation export CPU s | Observation export RSS KiB |
|---|---:|---:|---:|---:|
| A93 Zminus P3 | 4.781 | 5.358 | 8.186 | 1088648 |
| A95 XZminus P3 | 6.644 | 7.195 | 10.591 | 1454084 |
| A95 XZminus P4 | 7.337 | 7.912 | 13.997 | 2171992 |
| E94E95 Zminus P3 | 4.737 | 5.277 | 8.202 | 1085280 |

| A95 space / dt | Accepted steps | Process CPU s | Elapsed s | Peak solver RSS KiB | Accepted exclusive CPU s |
|---|---:|---:|---:|---:|---:|
| XZ P3 / .002 | 6 | 39.300 | 39.843 | 316536 | 31.936 |
| XZ P3 / .001 | 12 | 70.185 | 70.803 | 347092 | 62.935 |
| XZ P4 / .002 | 6 | 48.295 | 48.875 | 447460 | 40.705 |
| XZ P4 / .001 | 12 | 89.376 | 89.961 | 485208 | 81.856 |

Process CPU includes this segment program's startup/output/helper work. The last
column sums only nonoverlapping accepted-step intervals, excluding startup.
Owned-launcher elapsed time and resource suspensions are separately recorded.
For the six-step coarse-dt segments, accepted exclusive CPU is:

| Exclusive category | P3 CPU s | P4 CPU s |
|---|---:|---:|
| State-dependent domain and constitutive work | 22.366 | 25.819 |
| Volume/interface/external-face assembly | .819 | 1.440 |
| Interface constraints | .206 | .219 |
| Mass and pressure operations | .129 | .958 |
| Factorization | .128 | .632 |
| Solves / linear validation | .089 | .505 |
| Physical audits | 4.901 | 5.804 |
| Snapshot/checkpoint I/O | .547 | 1.691 |
| Other exclusive residual/orchestration work | 2.750 | 3.638 |

The primary binary does not separate geometry-row lookup inside those broad
categories. An isolated Release build adds opt-in buffered geometry, constitutive,
mass and face scopes. The six-step P3 and P4 profiles reproduce **every accepted
state and native history bitwise**. Buffered/nested timing tests pass. Open
uninterrupted/restart, rejected-step rollback, changed-traction history rejection,
absolute-pressure offset, internal reversed flow and ordinary closed CLS/cache
controls pass in the rebuilt runtime. Open ExactStepV1 remains disabled.

Detailed telemetry is intrusive: process CPU is 142.891 s for P3 and 155.301 s
for P4, respectively 3.636x and 3.216x their matched ordinary executions.
Its fine-grained percentages must not be interpreted as ordinary-runtime
fractions. Exclusive detailed measurements and call counts remain useful:
P3/P4 perform 30.76/31.32 million basis lookups, 15.38/15.66 million element-row
evaluations and .536/.568 million closure evaluations during the accepted segment.
All clocks/maps/output overhead remains visible. No estimated overhead is
subtracted to manufacture baseline category costs. No nested inclusive totals
are added to the exclusive tables.

## Provenance and fresh resource status

The primary solver SHA256 is
`858c9cdb8f343d0f68f109d4992123021f4d87d6e1a32d9cf4a7cb902f11984b`,
built from clean `a897896524b8761d3a954d491f96d1b33941d94f` in CT2. The primary
observer is
`43ecbaed68ecc2e6424e05df1b8090646f58cf9914ca82abd2758f6862d17475`.
Every new physical run uses these exact bytes and frozen libraries. CT2's
Zminus trajectories are reused with their sealed native state/configuration/
history identities and checked relocation bridge, not regenerated with new code.

The separate timing solver SHA256 is
`0cc403047b39ea882eb941db83420de2ca0e44fde198acd7a34231d2bc791f4c`,
from clean `421f79e8` and the recorded Release GNU 11.4/C++17 flags. It links
copied frozen CT2 libraries; the observation adapters are explicitly unchanged
CT2 bytes rather than advertised as rebuilt. All C++ changes in CT3 are optional
timing scopes/buffering. No governing equation, field row, physical signature,
input function, ownership rule or admissibility criterion changes.

[Source scope](cao-tamura-accepted-3-evidence/source-scope.json),
[final provenance](cao-tamura-accepted-3-evidence/final-provenance.json), requests,
seals and per-execution actual loaded-image hashes bind source, compiler flags,
runtime, mesh, native coefficients, live explicit/implicit boundaries and native
histories. The final read-only exporter runs from clean `651fa317`. Both source
specifications, all frozen dataset/role/weight files and the constitutive/wall
kernels hash-identically to CT2. Later commits contain reporting/evidence only.
Historical source/runtime receipts remain attached to their own binaries.

Fresh monitored resources pass: observed aggregate peak **22/72**, owned project
peak **13/16**, including both control helpers. Numerical libraries are one
thread before import/launch; only one numerical family runs at a time. The
isolated build uses two compiler jobs and bounded Git settings. Peak timestamps,
PIDs, parent ownership, thread counts, completed command exits, process CPU,
elapsed and pause events are retained. No suspension event is recorded in this
package. The final collector's own completed receipt is supplied separately
because it cannot include its not-yet-ended monitor in its output.
CT1's 79-thread incident and CT2's private-helper/foreign-debugger failures are
unchanged. Fresh compliance does not retroactively repair them.

## Scope carried forward

Two degrees and one bounded mesh contrast do not establish a convergence order.
The supplied finite-depth initial functions are inferred physical inputs, not
an exact steady solution or measurements of k/epsilon. Initial projection,
subsequent numerical drift, temporal differences and local constitutive
sensitivity are reported separately.

CT2's wall study is not repeated. Five millimetres remains a mathematical
frozen-reference candidate; roughness-sublayer adequacy is unestablished.
The native pilot still uses 16 mm, and smooth-reference preparation remains
unresolved. No developed precursor, experimental calibration/prediction,
native curved flow, EP/transition, MPI/3D, ExactStepV1 open caching or backend
expansion is attempted here.
The eight historical unresolved epsilon temporal-order observations remain
unchanged. These two-timestep pilot comparisons are uncertainty controls,
not replacements for those observations or measurements of temporal order.

## Next implementation gate

Implement **immutable geometry/inventory-row reuse for AcceptedOpenFlatV1** as
one bounded package, using these P3/P4 cases as equivalence references. Key rows
by exact mesh/field block, generating owner, explicit side, basis/point
configuration and unchanged physical coordinate; moving profile panels must
retain their current generating geometry. Reuse only immutable rows, never a
state-dependent constitutive or admissibility result. Require cache-off/on
bitwise states/histories, all native/dense domain decisions, pressure-action and
physical-budget equivalence, and a low-intrusion accepted-run cost comparison.
The measured broad domain cost and very high row-call counts support this step;
dense factorization is not the dominant cost in this prototype.

That optimization does not close the scientific spatial gate. Before larger
flat-flow development, resolve the outlet/delta local-radicand sensitivity and
the independent A93 covariance exit on unchanged physics with a separately
predeclared diagnostic. Physical S/R approach-flow preparation still needs the
inferred turbulence/finite-depth choices, roughness-sublayer/interface adequacy
and smooth-reference admissibility resolved. Native curved coupling requires
its own conservative geometry/metric/interface and actual-strain validation.
No precursor, hill prediction or longer response is admitted by this package.

## Reproduction structure

Use a **new** CT3 destination with the recorded source under `source`; never
reuse an existing label as a destination. The runtime paths below refer to the
preserved frozen CT2 binaries, not the shared installation. The original
commands, native-report hashes and immutable requests are exported with the
evidence. The recovery command for the first collector failure is also retained;
a fresh reproduction with the corrected script needs no recovery.

```bash
R=/media/jay/新加卷1/nektar-v2026-runs/NEW-ct-accepted-3
CT2=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-2-20261112
CT1=/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-1-20261111
S=$R/source/solvers/IncNavierStokesSolver/Tests/RANS
RT=$CT2/runtime
C=$R/source/docs/plans/cao-tamura-accepted-3-evidence/contract.json
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 TMPDIR=$R/tmp LD_LIBRARY_PATH=$RT/lib
ct3() {
  tag=$1; shift
  python3 -B "$S/ct_accepted3_command.py" --control "$R/launch-$tag" -- "$@"
}
ct3 controls python3 -B "$S/ct_accepted3_controls.py" \
  --ct2 "$CT2" --runtime "$RT" --output "$R/controls-v1"
ct3 p3-initial python3 -B "$S/ct_accepted3_run.py" --mode initial \
  --case A95-XZminus-P3 --runtime "$RT" --ct1 "$CT1" --contract "$C" \
  --output "$R/A95-XZminus-P3/initial-v1"
ct3 p3-short python3 -B "$S/ct_accepted3_run.py" --mode advance \
  --case A95-XZminus-P3 --runtime "$RT" --ct1 "$CT1" --contract "$C" \
  --initial "$R/A95-XZminus-P3/initial-v1" --target .008 \
  --output "$R/A95-XZminus-P3/T0008-v1"
ct3 p3-follow python3 -B "$S/ct_accepted3_jobs.py" \
  --root "$R" --ct1 "$CT1" --ct2 "$CT2" --case A95-XZminus-P3 \
  --initial "$R/A95-XZminus-P3/initial-v1" \
  --parent "$R/A95-XZminus-P3/T0008-v1" --first-target .02 \
  --output "$R/A95-XZminus-P3/continuation-v1"
ct3 remaining python3 -B "$S/ct_accepted3_suite.py" \
  --root "$R" --ct1 "$CT1" --ct2 "$CT2" --output "$R/remaining-cases-v1"
```

Run commands sequentially and inspect structured gates before any continuation.
`ct_accepted3_suite.py` continues independent cases after another stops; a
nonzero suite exit is not permission to resume the stopped case.
All changed timesteps have independent zero-step histories.

```bash
ct3 gate-tests python3 -B "$S/ct_accepted3_gate_tests.py" \
  --ct2 "$CT2" --output "$R/gate-tests-v1"
ct3 failure-repro python3 -B "$S/ct_accepted3_failures.py" \
  --root "$R" --ct2 "$CT2" --reproduce --output "$R/failures-v1"
ct3 profile-build python3 -B "$S/ct_accepted3_profile.py" --mode build \
  --root "$R" --ct2 "$CT2" --output "$R/profile-build-v1"
# Only admitted .008 -> .02 segments are eligible; use unchanged dt.
ct3 profile-p3 python3 -B "$S/ct_accepted3_profile.py" --mode profile \
  --root "$R" --ct2 "$CT2" --runtime "$R/profile-build-v1/runtime" \
  --parent "$R/A95-XZminus-P3/T0008-v1/0.002" \
  --reference "$R/A95-XZminus-P3/continuation-v1/T0p02/0.002" \
  --output "$R/profile-P3-v1"
ct3 profile-p4 python3 -B "$S/ct_accepted3_profile.py" --mode profile \
  --root "$R" --ct2 "$CT2" --runtime "$R/profile-build-v1/runtime" \
  --parent "$R/A95-XZminus-P4/campaign-v1/T0p008/0.002" \
  --reference "$R/A95-XZminus-P4/campaign-v1/T0p02/0.002" \
  --output "$R/profile-P4-v1"
ct3 profile-controls python3 -B "$S/ct_accepted_open_controls.py" \
  --runtime "$R/profile-build-v1/runtime" --old-runtime "$RT" \
  --output "$R/profile-compatibility-v1"
ct3 compare python3 -B "$S/ct_accepted3_compare.py" \
  --root "$R" --ct2 "$CT2" --output "$R/comparisons-v1"
ct3 compare-classify python3 -B "$S/ct_accepted3_compare_classify.py" \
  --input "$R/comparisons-v1/results.json" \
  --output "$R/comparison-classification-v1"
ct3 costs python3 -B "$S/ct_accepted3_costs.py" \
  --root "$R" --output "$R/costs-v1"
ct3 collect python3 -B "$S/collect_ct_accepted_3.py" \
  --root "$R" --ct1 "$CT1" --ct2 "$CT2" \
  --comparisons "$R/comparison-classification-v1" --output "$R/compact-v2"
ct3 final-supplement python3 -B "$S/ct_accepted3_finish.py" \
  --root "$R" --ct2 "$CT2" --compact "$R/compact-v2" \
  --output "$R/final-supplement-v1"
```

Do not edit a checkpoint's dt/signature, reconstruct a perturbed source, or
copy an evolved coefficient vector to another space. Frozen quadrature probes
use the explicit read-only replay identity, not altered restart metadata.
