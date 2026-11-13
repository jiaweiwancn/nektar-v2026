# Cao–Tamura budget AWF H2f

13 November 2026, Asia/Shanghai. **The inner-covariance defect in the search is
closed, useful native velocity freedom was found, and the bounded enriched
searches were executed. Outer covariance still rejects the best seed. No hill
physical step was attempted.** The changed-candidate V2 affine handoff and the
isolated direct sparse consistent-mass prototype pass.

## Case, provenance and corrected baseline

Work remains on `codex/cls-epawf-pimple-one-sweep`, in the Mac worktree
`.worktrees/pimple-one-sweep`, starting at reviewed head `e93c0fe`. The remote
matched that head. The validated implementation is `e2c7248f9d297d13c7d3342e0f4b85b0b9bf8837`;
[Git binding](cls-cao-tamura-budget-h2f-evidence/commit-binding.json) connects it
to the unchanged native source/runtime manifest. The shared Linux checkout and
unrelated project jobs were preserved. Compilation, all numerical conclusions and timings use the isolated
Linux root:

```
/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2f-20261113
```

UUID `5AB8A669B8A642FD` was verified at `/media/jay/新加卷1`. Four compiler workers,
one native rank and one computational-library thread were used, within the
aggregate 72-worker limit. USB SSH timed out; the established Tailscale path was
used. Frozen historical executables/libraries and candidate files were untouched.

The primary H_A=z_m=z_K=.0065 m, h_s=.0022472751006336437 m, seven-point mean fit,
friction velocity/viscosity, AnalyticalNeutralInflowV1, Regularized CLS
(beta=.4, Gamma=1, gamma0=0), Linear inner stress, budget cutoff, boundary roles,
one-sweep IMEX1 and 2804-slot mesh remain unchanged. There are six streamwise
and thirteen radial panels, cubic velocity/scalars and linear regional pressure.
No H/roughness/closure sweep, EP route, increased floor, stress clipping or
long continuation was introduced. Inlet fitting/held-out experimental data are
unchanged. The roughness-sublayer and flat spatial-accuracy qualifications remain
open.

The immutable PhysicalHeightInterfaceV2 base hash is
`83b1a4fd9a4fcda694f6e6f10bfb789700d1fc6fc83443b91f67bfcb045e79b5`.
Its exact retained H2d native receipt gives inner/outer minima
**+.00160236202103 / −.05289055348758**. The H2e baseline table has been corrected
from that receipt; H2e was not rerun. The Stage-A starting state is the retained
H2e A hash `5fbbd20a…`, with positive inner covariance. H2e B was never used as a
velocity start.

[Policy](cls-cao-tamura-budget-h2f-evidence/policy.json),
[native reduction](cls-cao-tamura-budget-h2f-evidence/native.json),
[verification](cls-cao-tamura-budget-h2f-evidence/verification.json),
[source/runtime identities](cls-cao-tamura-budget-h2f-evidence/source-revisions.json)
and [artifact index](cls-cao-tamura-budget-h2f-evidence/artifact-index.json)
are the compact record. Large typed arrays, states, native packets, checkpoints
and logs remain in Linux storage.

## Versioned context and original-span replay

`NativeOuterScalarContextV2` declares velocity/k/epsilon counts, named observation
blocks with offsets/shapes, global coefficient sizes/offsets and permitted rows,
basis hashes, geometry/reference/inventory identities, physical Gram and rank
policy. The reader verifies dimensions and identity consistency. V1 has one
explicit compatibility adapter; new search code has no `18:210` convention.

A zero-enrichment 32+8+8 export reproduces **every original binary array bitwise**:
coefficients, gradients, scalar/reference data, physical mass/metric forms,
observation records and all membership rows. Named-block interpretations and
conditioning match. The enriched context has 41+8+8 directions, 181975 unique
numerical records and 433207 original memberships. Both radial and tangential
sides, native geometry/element, and actual/dependent epsilon contexts remain
separate. Sharing computations does not remove observations or change their
ordering. Corrupted counts and geometry identities are rejected.

All original 32 physical velocity directions are retained. Their mass
conditioning and coordinate convention are preserved; the nine new directions
are already mass-orthonormal and orthogonal to the old span. Scalar bases remain
the original eight k and eight epsilon directions. No native discretization
DOF was added.

## Native accessibility and the one local enrichment

Ten clustered Rayleigh queries were constructed on the actual H2e-A state,
including slope, roof, matching-plane and protected inner witnesses. One
additional RHS reevaluates the retained H2e-B slope location on that same A
state. Eleven RHS in total are below the declared sixteen-RHS bound. Native
transported epsilon is used for these selected outer queries; their recorded
flags are actual, not dependent epsilon. Full inventory/trace checks retain
both contexts.

Each derivative uses native coefficient gradient rows and the selected full
stress law, at fixed scalars. Global constrained Riesz responses reuse the
native physical mass and continuity Schur factors. An independent reduction
finds projected KKT relative error about 2.44e-14. The accessibility number is
**eta squared**, a mass-energy fraction, not a feasibility certificate.

| Native chart (s,zHat), m | Old 32 eta² | Enriched eta² |
|---|---:|---:|
| (.06796,.10277), slope | .001606 | .999990 |
| (−.05,.4), outer panel face | .000506 | .999999 |
| (.1,.50279), foot-minus outer | .032136 | .999994 |
| (.02247,.14493), slope | .010766 | .999979 |
| (−.05625,.8), roof | .001304 | .999996 |
| (.07262,.0065), outer matching limit | .00000366 | .414153 |

The initial patch is chart [-.1,.1] x [.0065,.4]. It has 374 free coefficient
rows, continuity rank 34, combined continuity/old-span rank 66 and nullspace
dimension 308. It misses upper witnesses and fails the predeclared capture
rule. The **one** permitted enlargement is [-.2,.6] x [.0065,.8]: 738 rows,
continuity rank 56, combined rank 88, nullspace dimension 650.

Only coefficients whose entire native element support lies in the patch are
allowed. All affected global Q rows participate. Old-span mass orthogonality
is imposed inside the patch constraint system; subtracting a global response
afterward would destroy support. Twice-reorthogonalized covariance-gradient
responses produce **nine** independent directions, not sixteen quota fillers.
Their old-span overlap is about 1.4e-15 and homogeneous Q defect about 1.1e-13
at unit mass norm. Their inner velocity gradients are exactly zero throughout
the exported inventory. The matching-plane accessibility remains limited;
this experiment does not exhaust native velocity freedom.

Constraint-aware one-direction tests explain why high accessibility alone is
insufficient. Starting from H2e A, each retained new direction is limited first
by the physical divergence bound. Uncompensated positive steps are only
2.85e-6–3.93e-5 in mass coordinates, versus a maximum trust radius .04. Their
actual local Rayleigh gains are .020–.342 K0; those gains do not remove the
whole outer deficit or guarantee simultaneous improvement at other witnesses.
The full optimizer permits combinations and compensating old directions.

## Hard inner protection and executed searches

With frozen inner k and auxiliary viscosity, R_A=(2/3)k I−nu_A(G+G^T) is affine
in velocity. For each fixed critical eigenvector, its Rayleigh quotient gives
an exact linear cut. **These cuts have no outer phase-I slack.** All nearly
repeated critical eigendirections within 1e-7 K0 are retained. Complete inner
scans independently protect the predictor, proposed accepted trials, line
search, immediate best-state saving and terminal selection. A repeated-eigenvalue
fixture demonstrates why canonical-axis cuts alone are insufficient.

The inner target is lambda_min/K0=1e-4, with the existing normalized arithmetic
allowance 1e-10; native hard acceptance tolerances are unchanged. The lowest
accepted inner margins are about 9.99999e-5, safely positive. Initializer-only
strong-divergence and strain/roof controls do not alter the production momentum
equation or roof boundary policy.

All bounds concern total departure from the immutable base: physical velocity
L2 <=.12, sampled component increment <=.5 m/s, GL16 divergence L2
<=.0603080994, independent GL20 allowance 1e-8, outer scalar ratios [.5,2].
The pressure, inner k, inlet/interface scalars, endpoint K and dependent epsilon
remain fixed exactly during initialization. None becomes an evolution clamp.
The native validator confirms fixed rows, moving-profile locations, cutoffs and
trace data independently.

| Metric | Enriched A, frozen scalars | Enriched B, outer scalars free |
|---|---:|---:|
| Major iterations / stop | 20 / budget | 40 / budget |
| Search wall seconds | 42.06 | 143.05 |
| Native inner minimum | +1.759889e-6 | +1.759889e-6 |
| GL20 inner minimum | +3.530622e-5 | +3.298460e-5 |
| Native outer minimum | −.034267843 | −.029804047 |
| GL20 outer minimum | −.034385617 | −.030013454 |
| Velocity L2 increment | .04699357 | .04846096 |
| Maximum sampled component increment | .43143691 | .43115801 |
| Divergence L2 | .06030809942 | .06030809099 |
| Conservative full-scan work | 103 / 180 | 270 / 320 |
| Active derivative batches | 160 / 3000 | 480 / 6000 |
| Maximum retained inner cuts | 51 | 30 |

A preliminary Python name collision stopped before the optimizer; its failed
scan is charged to A's aggregate work and the log is retained. No major
iteration was repeated. At A's twentieth iteration a relative-progress branch
initially labelled the stop as stagnation even though a small improvement was
accepted. Budget exhaustion takes precedence; the original receipt is retained
beside the corrected label. Coefficients/hashes were not edited. B likewise
exhausts forty majors, with improvement at iteration 38 and rejection at 39.
Local subproblem limits are recorded; only fully checked improving states are
accepted. This is not a global incompatibility certificate.

The single permitted hybrid uses exact H2e-A velocity and H2e-B values only on
allowed outer scalar rows. Its native outer minimum −.03728150 is worse than
the best enriched A result −.03426784, so it is not the B start. Its inner and
GL20 margins are positive. This was one combination and one native audit.

At B, divergence and the inner margin bind. k reaches its upper ratio 2;
epsilon is approximately [.55,2]. The .12 mass-distance and .5 component caps
are **not** active at the accepted final field, although response/scalar rows
enter the local trial working sets (up to 139/1092). There is no extra response
allowance on top of the starting state.

B's integrated initial TKE increases from .0112015608 to .0124298530 (+10.97%).
This is an initializer departure, not correction supply or a time-integrated
source. Exact scalar envelopes and k/epsilon time- and length-scale departures
are in the native reduction. Inner fields and endpoint K remain unchanged.

The B outer witness is at physical (x,z)≈(−.0250901,.1649363), native
(s,zHat)≈(−.0226247,.1359198), geometry/element 62. Its k=.03489343 and
epsilon=.00842951. Along the independent critical eigenvector, isotropic,
linear, quadratic, beta-weighted cubic and regularization contributions are
+.02326229, −.00605107, −.00673405, −.00060696 and −.03967425. Their signed sum
is −.02980405. Independent tensor replay errors are below 1e-17. These are local
stress contributions, not additive percentages of a flow budget.

Finite-inventory gradient maxima are retained by region. B has gradient
Frobenius maxima approximately 2653.8/s inner and 55.46/s outer; the large inner
wall gradient must be distinguished from its increment (about 386.6/s). Sampled
local divergence reaches 16.46/s inner and 2.13/s outer despite small native weak
continuity defects. GL20 divergence and energy measures are separate. No
all-element covariance theorem or loss-of-ellipticity conclusion follows.

## V2 startup handoff, with a genuinely changed candidate

The new V2 runner uses the **already selected executable and library directory**;
it does not refreeze/move them after certificate issuance. It sets the actual
startup-version information and validates receipt/state/runtime identities.
A changed runtime requires revalidation and a new certificate. The conditional
hill runner forms independent coarse/fine histories, stops each leg at its
first rejection, supports coarse 4+4 restart and post-assignment rollback.

A positive reference-compatible affine fixture changes permitted native k/epsilon
rows 1960/2522 by 1.759889e-6 / 1.430513e-5, respectively. It is independent of
the open mapped-hill qualification and is not an experimental precursor.
Independent complete native/GL20 checks pass. Its candidate hash differs from
the immutable base/reference. Exact imported seed, fixed reference, defect,
prescribed data, endpoint K and protected epsilon retain their intended distinct
identities; import performs zero Lift/projection.

Two direct steps agree bitwise with 1+1 split evolution. Every accepted step
has one U/p/epsilon/k solve and one logical bound call; solved epsilon is retained.
A deliberate post-assignment rejection restores state, histories, time labels,
reference/defect and samples, preserving parent checkpoint bytes. Ordinary
changed-dt restart and a changed state presented with its old receipt are
rejected. V1 h0/history still reproduces the retained regression bitwise.

Three complete hill-candidate certifications were used: A, hybrid and B, below
six. All fail outer covariance. **No startup certificate was issued for a hill,
no hill attempt or step occurred, and curved 8/16 stepping and 4+4 restart remain
unexecuted.** An affine passing fixture cannot waive this gate.

## Direct sparse physical-mass prototype

The opt-in diagnostic assembles global CSR/block mass directly from 78 native
physical element contributions, using signed shared CG IDs and physical
quadrature. The sparse factor receives those local contributions, never the
exported dense reference. The current native fixture still constructs its dense
baseline for comparison; this is not a whole-solver storage replacement.

| Current free velocity mass | Measured result |
|---|---:|
| Rows / nonzeros | 1386 / 31164 |
| Dense matrix bytes | 15,367,968 |
| CSR bytes | 379,516 |
| Explicit L/U bytes / fill ratio | 1,764,908 / 4.69 |
| Direct local assembly | .0234 s |
| Sparse direct factor, three trials | .00496–.00500 s |
| Dense reference factor | .04645 s |
| Sparse multi-RHS solve | .00102–.00113 s |
| Matrix relative max / row-scaled error | 1.34e-15 / 2.43e-15 |
| Worst per-RHS inverse / native constrained response | 2.25e-13 / 2.47e-13 |
| Mass-dual relative error | 3.71e-15 |

Essential elimination error is about 6.94e-18. Unit-mass normalized Q response
error is about 1.57e-12; unscaled response defects are also reported. Explicit
L/U storage excludes permutation/workspace overhead. Peak process RSS before
loading the dense reference is 49152 KiB and the full comparison peaks at
112268 KiB. These timings do not imply a campaign speedup: production already
reuses its unchanged dense mass factor, and inventories/audits remain costly.

Bounded symbolic support graphs reproduce the current nonzero count. The retained
4399/5339-slot allocations imply free-velocity masses of 2250/2682 rows,
51324/61404 structural entries and approximately .625/.748 MB CSR storage,
versus 40.5/57.5 MB dense free masses. Fill-only surrogate estimates are labelled
explicitly; no physical target factorization or >3000-slot evolution was run.
Pressure Schur handling stays dense. Future integration must preserve
B=(alpha M_free)^-1, complete pressure/roughness/outlet action and the exact
consistent-mass k-bound objective. Flow/scalar dense assembly and the unclosed
flat spatial gate still require a separate implementation package.

## Costs, gates and next decision

Zero/enriched context construction takes 37.42/43.80 s, or 41.15/47.52 s inclusive,
with roughly .765-GiB native peaks. A/B search peaks are approximately
1.13/1.14 GiB; exact runtime/receipt hashes are recorded separately. The changed affine V2 zero-step import costs
19.87 s and its two-step block 33.96 s inclusive, with peak about 1.42 GiB.
These include initialization, full audits and output; there is no physical
hill-step timing. At the final 16:11 Asia/Shanghai process check, no H2f native/build job or
continuation remained running.

| Gate | Outcome |
|---|---|
| Corrected receipt-based baseline and V2 zero-enrichment replay | passed |
| Native accessibility and one patch enlargement | useful additional freedom found; nine directions |
| Inner protection throughout accepted states and native/GL20 certification | passed |
| Outer covariance qualification | failed |
| Scalar/budget/frozen-row/boundary/weak-continuity checks | passed |
| Changed-candidate V2 import, direct/split, rollback and V1/negative checks | passed |
| Direct sparse physical mass and symbolic scaling | passed diagnostic; production unchanged |
| Conditional hill pair / curved restart / long-time or experimental validation | unexecuted / unqualified |

**Next recommendation:** a bounded divergence-aware initializer response or
compensating-direction construction, retaining the present inner/scalar/velocity
bounds and physical model. All nine isolated new responses are divergence-limited;
adding more mass-optimal covariance responses alone is not justified. The current
41-dimensional search does not exhaust the 650-dimensional outer patch nullspace
or prove the bounds incompatible. No scalar-envelope widening, endpoint-K
release, closure retuning or automatic optimization extension follows.

## Reproduction

Use a fresh run root; never overwrite the indexed artifacts. The build was:

```sh
cmake -S "$H2F_ROOT/source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$H2F_ROOT/build" -DCMAKE_BUILD_TYPE=Release \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS' \
  '-DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake' \
  -DRANS_BUILD_SOLVER=ON \
  '-DRANS_PREBUILT_SPATIAL_LIBRARY=/media/jay/新加卷1/nektar-v2026-runs/cls-cao-tamura-budget-h2b-20261112/curved-lazy-final/runtime/lib/libSpatialDomains.so.5.9.0'
cmake --build "$H2F_ROOT/build" --parallel 4 \
  --target RANSWallResolutionProbe IncNavierStokesSolver RANSKBoundQPChecks
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
```

H2F_TESTS denotes `source/solvers/IncNavierStokesSolver/Tests/RANS` inside that
root. H2F_PREFLIGHT is the retained H2c `preflight-final.json`, H2F_OLD is the
H2e run root, and H2F_H1 is H1 `compatibility-resolved`. The main sequence is:

```sh
python3 -B "$H2F_TESTS/run_cao_tamura_budget_h2f.py" select \
  --context "$H2F_OLD/context-v1/context.json" --preflight "$H2F_PREFLIGHT" \
  --start "$H2F_OLD/stage-A-final/best-candidate.json" --output "$H2F_ROOT/selection"
python3 -B "$H2F_TESTS/run_cao_tamura_budget_h2f.py" native \
  --probe "$H2F_ROOT/build/RANSWallResolutionProbe" \
  --context "$H2F_OLD/context-v1/context.json" \
  --selection "$H2F_ROOT/selection/selection.json" --output "$H2F_ROOT/native-data"
```

The `retained-slope` operation adds the one recorded H2e-B location, using
`--selection`, `--context`, `--start`, `--preflight`, `--output`; a second
`native` call evaluates that one RHS. `combine --native ... --extra ... --output
...` joins only those already evaluated arrays. Regenerating the extra query
and combination reproduces their recorded witnesses/arrays bitwise, without
another native solve.

```sh
python3 -B "$H2F_TESTS/enrich_cao_tamura_budget_h2f.py" \
  --native "$H2F_ROOT/native-combined/native-data.json" \
  --context "$H2F_ROOT/zero-context/context.json" --output "$H2F_ROOT/enrichment"
python3 -B "$H2F_TESTS/run_cao_tamura_budget_h2f.py" context \
  --probe "$H2F_ROOT/build/RANSWallResolutionProbe" \
  --context "$H2F_OLD/context-v1/context.json" --preflight "$H2F_PREFLIGHT" \
  --basis "$H2F_ROOT/enrichment/velocity.f64" --output "$H2F_ROOT/enriched-context"
python3 -B "$H2F_TESTS/search_cao_tamura_budget_h2f.py" \
  --context "$H2F_ROOT/enriched-context/context.json" --preflight "$H2F_PREFLIGHT" \
  --start "$H2F_OLD/stage-A-final/best-candidate.json" --stage A \
  --prior-full-scans 1 --output "$H2F_ROOT/stage-A-final"
```

The zero context uses the same `context` command with
`selection/velocity.f64`, before enrichment. A native/GL20 certification uses
`certify_cao_tamura_budget_h2f.py --context ... --candidate ... --preflight ...
--probe ... --output ...`. The recorded hybrid audit and B-start selection
precede B:

```sh
python3 -B "$H2F_TESTS/search_cao_tamura_budget_h2f.py" \
  --context "$H2F_ROOT/enriched-context/context.json" --preflight "$H2F_PREFLIGHT" \
  --start "$H2F_ROOT/stage-A-final/best-candidate.json" --stage B \
  --output "$H2F_ROOT/stage-B"
```

Before final certificate issuance the intended solver/probe closures are frozen
in `pilot-solver`/`pilot-probe`. B certification supplies that exact solver via
`--solver`. Its failed gates prevent certificate creation. The conditional
`run_cao_tamura_budget_h2f_startup.py` requires a real V2 certificate and the
already frozen solver; it was **not** run on a hill.

`verify_cao_tamura_budget_h2f_startup.py` executes the changed affine fixture;
`check_cao_tamura_budget_h2f_startup.py` adds the narrow V1/negative checks.
`prototype_cao_tamura_sparse_mass.py --native ... --context ... --output ...`
runs the direct local mass comparison. `symbolic_cao_tamura_mass.py --output
...` constructs only the labelled support graphs. `check_cao_tamura_budget_h2f.py`
and the retained mass-QP/branch checks pass. `analyze_cao_tamura_budget_h2f.py`
and `collect_cao_tamura_budget_h2f.py` produce the compact reduction and manifest.
Native input packets, receipts, source overlays and exact commands are indexed.
