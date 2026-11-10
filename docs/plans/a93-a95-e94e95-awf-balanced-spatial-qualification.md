# Balanced spatial qualification

Actual starting revision: `3b0961e12142d7eca56d75962feebfffa34630cc`.
The existing implementation branch and isolated Mac worktree are clean, with
no subsequent local commits. Preserve all historical receipts, including the
eight unresolved epsilon-order observations. Specifications retain SHA256
`056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64`
and `17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.

The authoritative new Linux root is
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-balanced-20261110`,
with separate source/build/dependency/frozen runtime/tmp/run directories.
The primary Linux checkout is dirty on an unrelated branch and remains untouched.
The previous pilot root and active CLS campaigns are read-only. Initial relevant
process census is 15 threads; allocate eight build workers and at most four
serial native controllers after recounting. Numerical library thread counts are
one, including Python, and the aggregate hard limit remains 72.

## Frozen contract, before new experiments

No governing closure, coefficient, wall source, pressure convention, one-sweep
ordering or native history policy changes. High-Re outer molecular diffusion
stays zero, physical inner molecular viscosity stays finite, all three normal
stresses and full production are retained. Sources are uncapped, corrections
remain disabled, eta is zero. No reference cancellation, repair or resetting.
The existing 3000-unknown evolution guard stays in place. Unsupported routes
remain guarded.

Supplement the previous export into a new directory with all twelve complete
required comparison bodies, endpoint/request/runtime/history/consumed-input
identities and fieldwise limits. Missing comparison coverage and a scaled
outer-epsilon constant-test defect above 1e-10 must return nonzero. Its scale is
outer volume * epsilon(H) * u_star/H, independent of timestep and mesh. Historical
printed values and thresholds are not overwritten.

NeutralV1 is eligible only for canonical positive-x neutral flow. Its existing
D array entries are explicitly labelled change from the fixed neutral polynomial
representation; Python's independently sampled change_from_initial continues to
mean the actual initialized field, including a declared perturbation.

Residual diagnosis is read-only, model-aware alpha=0. Export the existing native
physical terms, free equations, constrained-row reactions, explicit constraint
replacement, native mass matrices, element support and one-sided fluxes. Require
term sum versus unmasked native residual and constraint replacement versus the
constrained residual within 2e-11 times maximum(1, sum of absolute term loads).
This is an accumulation comparison, not an equilibrium tolerance. Independent
analytical weak balances, represented-state residuals, evolved spatial residuals
and stage/splitting receipts remain separate.

Keep full native mass-dual norms with the existing physical normalization.
Additional common-test results use an unchanged base p3 native test space and
its mass. Verify test prolongation by value/gradient evaluation, mass congruence,
constraints and basis orientation; restrict fine loads by P^T. Any mesh that
omits an inactive auxiliary breakpoint uses the explicitly documented common
physical p3 test subspace without that breakpoint, for every comparison in that
control. Newly introduced modes stay visible in full native norms.

Screen h1/p4 and localized p5 hierarchies before advancement. After localization,
record exact selected partitions and their rationale before qualifying evolved
pairs. For dominant momentum-x and inner-k defects above a 1e-8 scaled floor,
require at least 15% reduction in both full native and common-test norms for the
selected genuinely nested pair, beyond five times quadrature/evaluation
uncertainty. Selected finest initial defects must also improve over base p5;
other equations and new modes must be reported, not dropped. A nonnested h1/p4
comparison is a candidate screen, not an h convergence measurement.

Retain field weighted RMS limits U,k,epsilon,Pi <=1e-3, w<=1e-4,
wall/interface traction <=1e-3, on the existing fixed physical scales. Timestep
uncertainty for each member must satisfy delta_t <= max(1e-7,.1*delta_h).
Use identical actual inputs within same-space temporal controls and valid
independently initialized native histories. Compare on the union of physical
partitions and active breaks with one global pressure gauge.

Start E94E95 natural both walls at T=.0016/.016; continue other model/wall and
selected prescribed controls independently where supported. T=.5 (theta=.1)
is conditional on the balanced spatial and timestep gates, and requires a
finer-space and dt/2 companion at that same horizon. Preserve the first failure,
last accepted history and all diagnostic evidence. No broad longer campaign or
new backend is authorized by this package.

## Results

Pending native term diagnosis and candidate screening; this is not a completion
receipt. Results, exact reproduction commands and remaining gates follow in
logical commits.

## Initial screen and localized candidates (before advancement)

Both E94E95 walls pass the native term/reaction reconstruction and common-test
inclusion checks on base p5, Interface8, Nearwall8, global h1/p5, h1/p4 and the
inactive-breakpoint mesh-only control. Actual global h1/p4 count is 2471;
global h1/p5 is 3801 and remains ineligible for evolution. h1/p4 increases full
momentum and accepted inner-k norms relative to base p5, despite smaller loads
on shared p3 tests, so it is not selected merely for its DOF count.

Signed native dual-energy localization assigns about 62--63% of base momentum
energy to the tiny original-epsilon panel, and about 91% of accepted k energy
to that panel. EP k energy is mainly .02--.04 (63%), .04--.1 (20%) and .1--.25
(12%). Nearwall8 removes the accepted near-wall defect but EP's new .07--.1
modes dominate its increased full norm. Both walls improve on fixed shared
tests under Nearwall8. Omitting the inactive original-epsilon mesh breakpoint
reduces full norms but leaves common momentum nearly unchanged and slightly
increases accepted common k. Thus removal alone is not physical consistency
improvement and does not pass the balanced gate.

The next bounded initialization screen is a nested p5 hierarchy on the base
mesh with only that inactive epsilon mesh breakpoint omitted. Active viscosity,
cutoff, roughness/forcing and interface breaks remain; profile quadrature is
unchanged. `Balanced5V1` splits inner [0,.0025], [zc,.02], [.02,.04], [1,2] and
outer [2,2.25]. `Balanced11V1` also splits [.0025,.005], [.04,.1], [.1,.25],
[.25,.5], [.5,1] and outer [2.5,3]. Expected counts are 2355 and 2997 and must be
verified natively. Both use velocity/k p5 and pressure p3, quadrature
48/12/12/12/12. The exact per-model/policy cutoff is computed from its descriptor;
no mesh motion is allowed. This allocation follows assembled-load support
energy rather than raw element stress loads that cancel between neighbors.
The predeclared 15% full/common dominant-residual reduction, field and temporal
uncertainty gates remain unchanged. No trajectory is admitted until this
initial screen establishes a meaningful pair.

The E94E95 pair passes at both .0016 and .016, including independent dt/2
controls on **both** members. At .016 fine/coarse reference-error ratios are
about .040 for U, .034 for w, .027--.045 for k, .325 for epsilon and .058 for Pi.
These are spatial/uncertainty receipts, not measured temporal orders.
Before A93/A95 screening, require their non-floor vertical-momentum defects to
reduce on both full and common tests as well as the already declared U/k gate.
Use the identical selected mesh strategy, including a base-p5 benchmark.

First read-only evolved replays exposed two diagnostic mistakes, retained in
`replay-previous*.log` and their complete native packets. An exact time equality
rejected an archived time .016-1.04e-17; use the existing 1e-13 endpoint tolerance
and retain both actual and nominal time. More substantially, near moving
cutoffs the native field lookup allows roundoff-outside points and its basis
evaluator clamps local coordinates. Side-zero lookup and raw coordinates do
not always describe the fine test row's actual polynomial branch. The new
**read-only** test-basis export selects the exact parent element; common tests
also use the fine basis's effective (clamped) coordinate. Export the maximum
raw/effective displacement. Do not change native lookup, clamping, equations,
quadrature weights, physical source coordinates or identity tolerances.
Separate frozen solver and revised diagnostic-probe runtimes preserve all
ongoing trajectories. The first failed common-test receipts remain failures.

## Duration-cost contract (before duration-cohort runs)

E94E95 .016 costs are 5.9--7.5 seconds/step for Balanced5 and 7.0--8.5 for
Balanced11. Extending the .0002 fine companion would take several hours.
To keep the conditional experiment affordable, use a **fresh**, independently
initialized duration cohort at dt=.004/.002, on both members of the already
selected pair. Its milestones are .016,.1,.5. This changes neither a stored
history nor its timestep/signature; within each new trajectory continuation
through .1 to .5 uses genuine unchanged-method native restart. The original
.0004/.0002 short histories and receipts remain immutable.

Require the same field/traction/admissibility/conservation gates and
`delta_t <= max(1e-7,.1*delta_h)` on both spaces at **each** new horizon. No
new temporal order is inferred. Admit this cohort only for an evidence-selected
limiting combination after its full/common spatial and short uncertainty gates
pass. If the coarser duration pair fails uncertainty while admissible, one
separately initialized .002/.001 cohort is the predeclared targeted refinement;
do not relax the criterion or alter an existing checkpoint. Any hard accuracy
or physical failure retains its rejected state and last accepted history.

## Verification corrections and supplemented evidence

`collect_rans_neutral_pilot.py` now exports all twelve required detailed bodies
from `natural-schedule`, including both A93 .016 bodies. Each carries three
case/request/space/dt/runtime identities, actual consumed raw/effective/reference
and both live-boundary hashes, endpoint/report/sample/state/history identities,
fieldwise errors, differences, limits and classifications. The supplement is a
new export; historical receipts are untouched. An isolated omission of
A93/Accepted/.016 returns exit 1 for the intended missing-file reason. Printed
outer-epsilon load/reaction defects now have an explicit scaled classifier;
a deliberately doubled tolerance witness returns exit 1, while the recognized
negative-control suite returns 0.

`NeutralV1` rejects a negative-x reference specifically as diagnostic
ineligibility. Forward closures retain their existing supported scope. The
native report explicitly labels legacy D entries as change from the fixed
neutral polynomial, not drift from a perturbed initial state. Both positive-x
E94E95 wall controls preserve every accepted state, native history and numeric
state-audit receipt bitwise against 3b0961e's runtime. A zero-step A95 velocity
perturbation separately has zero drift from its actual initialized vector and
nonzero change from the fixed neutral reference. The first fixture failed
because it pre-created the generator's exclusive output directory; its log is
preserved and the corrected fixture uses a fresh directory.

Explicit-parent native test exports have valid, unknown-geometry, wrong-count
and conflicting-interface-side controls. All intended results are observed.
The targeted CTest selection passes 3/3, including 724 spectrum/domain/ground
assertions in the native helper. Unchanged constitutive/wall kernels,
evolution/history/cache/rollback sources and independent analytical-oracle
contents are checked by content against the reviewed baseline. Reuse is scoped
to those unchanged components; new observation paths receive fresh checks.

## Why the wall refinements differ

Element stress loads individually contain large internal-boundary contributions
that cancel after assembly. They must not be interpreted as additive local
residual norms. The export retains those physical terms and reactions, and
additionally partitions the assembled dual energy r^T M^-1 r by exact native
basis-support mass fractions. Signed contributions sum to the original squared
norm; coordinates, support geometry and field/equation are retained.

On the base mesh the inactive original-epsilon panel contributes approximately
62--63% of momentum dual energy and 91% of accepted inner-k energy. EP k is
instead dominated by .02--.04, .04--.1 and .1--.25. Nearwall8 removes most of the
accepted wall's small-panel defect, but the newly available EP modes on .07--.1
then carry about 79% of its k dual energy. Its full native norm increases even
though the same fixed common tests improve. Interface8 resolves the EP bulk
panels and interface while leaving almost all accepted small-panel energy.

Removing only the inactive epsilon mesh breakpoint lowers full momentum from
.108/.111 to .0659/.0658 (EP/accepted), and accepted inner k from .1008 to .0275.
Common momentum is essentially unchanged; accepted common k slightly increases
(.007652 to .007885). Thus that mesh-only control demonstrates sensitivity to
new test modes, not a cure for the physical spatial defect. Every active
viscosity, cutoff, roughness/forcing and interface break remains represented
and native profile integration retains its physical panels.

The selected 2355/2997 hierarchy resolves both the small-panel/budget region and
the dominant bulk/interface contributions. The fixed common test space is p3
with pressure p1 on the active-break base mesh, unchanged across all candidates
within one model/wall/policy. P is constructed from native basis cross mass;
value/gradient, mass congruence, exact free-ID constraints and orientation checks
precede P^T r restriction. Full native norms remain alongside those restrictions.
Raw-to-effective coordinate shifts in evolved replays are at most 1e-11 and
are explicit. Native physical/source coordinates and weights are never changed.

## E94E95 residual results

All entries below use the existing dimensionless physical mass-dual scales.
Common k includes the assembled shared equation once; native inner-k excludes
the separately reported shared rows. Small differences between them are not
silently identified as the same norm.

| Wall | Space/time | Native U | Native inner k | Common U | Common k |
|---|---|---:|---:|---:|---:|
| EP | B5 initial | .131300 | .00183904 | .0199993 | .000241191 |
| EP | B11 initial | .0134396 | .000184628 | .00153787 | .0000126857 |
| EP | B5 .016 | .000560082 | .000155467 | .000392345 | .000125627 |
| EP | B11 .016 | .0000126612 | .00000216634 | .00000357527 | .000000740115 |
| Accepted | B5 initial | .130617 | .0563665 | .0200098 | .00826489 |
| Accepted | B11 initial | .0127152 | .00803083 | .00156286 | .00107652 |
| Accepted | B5 .016 | .000572601 | .000159732 | .000399322 | .000127750 |
| Accepted | B11 .016 | .0000129984 | .00000234189 | .00000366510 | .000000755087 |

The term sum reproduces the native **unmasked** operator within the predeclared
accumulation tolerance. Constraint replacement is checked independently against
the constrained result; prescribed-row loads remain reactions. Inner/outer
stress, conservative advection, k diffusion/production/prescribed destruction,
outer k, model-specific epsilon diffusion/production/destruction, continuity,
gauge and shared rows are separate. The diagnostic never feeds evolution.

The source/profile and outer/mass/interface quadrature checks increase
48/12/12/12/12 to 64/16/16/16/16 with identical coefficients and live boundaries.
Both initial spaces and both evolved .016 states pass the 1%-of-spatial-change
criterion with absolute floor 1e-10. GL10/8 versus GL14/10 physical evaluation
also passes on identical coefficients. The actual rules, projection quadrature,
profile panels, frozen state hashes and target boundary identities are recorded.

Independent analytical values/gradients evaluated against native tests verify
the continuum weak balances. EP plateau source moments vanish; accepted source
moments are nonzero and cancel with diffusion and boundary terms. These checks
remain distinct from the represented and evolved spatial residuals above and
from time-discrete stage residuals.

The outer-epsilon constant test is a limited receipt on this prescribed-side,
single-column patch: all its nonzero constant coefficients are constrained.
Consequently load minus prescribed-row reaction is bookkeeping, not an
independent free-equation conservation theorem. Its finite scaled tolerance
still classifies failures. Nonconstant free tests, source/storage identities,
full native/common residuals and interface balances provide separate evidence.
No zero dependent-epsilon boundary rate is assumed; no semidiscrete scalar-rate
claim is made.

## Duration-admission outcome so far

The scoped admission selects E94E95 Natural Accepted, whose nonconstant wall has
the larger qualified E94E95 scalar residual and reference error. The .004/.002
fresh-history cohort passes execution, physical and field-accuracy checks at
.016 but **fails pressure timestep uncertainty**: 4.1630e-7 exceeds 1.9186e-7.
Its nonzero process result and all four trajectories remain archived. The
predeclared .002/.001 independent-history refinement is the next experiment;
no .5 result is inferred from the successful short .0004/.0002 study.

## Source mapping and reproduction

| Responsibility | Implementation / evidence |
|---|---|
| Existing physical term assembly, with optional per-element export | `RANSWallEvolutionResolution.cpp`, `ResidualResolution` and `ResidualResolutionComponents`; production default does not collect elements |
| Model-aware immutable replay and native test selectors | `library/UnitTests/SolverUtils/RANSPreservationProbe.cpp`; legacy CLS replay guard remains |
| Diagnostic eligibility and D semantics | `RANSEquilibriumWallEvolution.cpp`, `NeutralAccuracy` in `RANSWallEvolutionDiagnostics.cpp` |
| Full/native and fixed common-test analysis | `rans_spatial_resolution.py`, `explain_rans_spatial_residual.py`, `audit_rans_balanced_evolved.py` |
| Complete previous bodies and epsilon budget classifier | `rans_pilot_evidence.py`, `collect_rans_neutral_pilot.py`, `check_rans_balanced_evidence.py` |
| Selected fixed mesh / stage gates | `make_rans_balanced_meshes.py`, `classify_rans_balanced_spatial.py`, `compare_rans_balanced_evolution.py` |
| Same-space inputs, histories and validated reuse | `run_rans_neutral_pilot.py`, immutable request manifests, `run_rans_balanced_duration.py` |
| Scoped duration admission and failed pressure control | `admit_rans_balanced_duration.py`, `duration-admission.json`, the separate .004 and .002 cohort receipts |
| Source reuse and complete export | `record_rans_balanced_source.py`, `collect_rans_balanced_spatial.py` |

The solver runtime is frozen at source snapshot 577b5f1 with executable SHA256
`389dee61618f51c9d451be2ddbdc0e83ee049351683da2eae3b59ab0758afcb4`.
The revised explicit-parent diagnostic probe is separately frozen at 0a5496f
under `runtime-resolution-v2`; its solver copy is byte-identical and is not
used to replace any trajectory executable. Manifests include relevant source
file contents/hashes, compiler/options, CMake cache, executable, direct linked
libraries, actual loaded libraries and thread environment. Later changes are
Python observation/orchestration and documentation. The source proof checks
34 governing/history/geometry files and three independent-oracle files against
3b0961e. Ordinary positive-x histories and numeric diagnostics are tested for
bitwise equivalence; unsupported routes retain their guards.

All heavy native states, histories, field data, basis/residual packets, logs,
failed inputs and timing streams remain outside the checkout under the Linux
root. `runtime` and `runtime-resolution-v2` are immutable. Source/build/tmp/
run storage is isolated. Neither a shared installation nor any active CLS
controller/library/checkpoint is replaced.

These commands use a **new** Linux destination and the existing MPI-enabled SDK
without installing into it. Recount other jobs before allocating workers.
Python >=3.10, NumPy and SciPy are checked at CMake configuration.

```bash
export B=/media/jay/新加卷/nektar-v2026-runs/balanced-reproduction
export ORIGINAL=/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-balanced-20261110
export OLD=/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-pilot-20261110
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$B/dependency" "$B/tmp"
git -C /media/jay/新加卷/nektar-v2026 worktree add --detach "$B/source" TESTED_HEAD
cp "$ORIGINAL/runtime/lib/libSpatialDomains.so.5.9.0" "$B/dependency/"
export T="$B/source/solvers/IncNavierStokesSolver/Tests/RANS"
export TMPDIR="$B/tmp"
cmake -S "$T" -B "$B/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DRANS_PREBUILT_SPATIAL_LIBRARY="$B/dependency/libSpatialDomains.so.5.9.0" \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$B/build" --parallel 8 --target IncNavierStokesSolver \
  RANSPreservationProbe RANSNewModelsProbe RANSNativeDiagnosticChecks
python3 -B "$T/freeze_rans_pilot_runtime.py" --source "$B/source" \
  --build "$B/build" --output "$B/runtime"
python3 -B "$T/freeze_rans_pilot_runtime.py" --source "$B/source" \
  --build "$B/build" --output "$B/runtime-resolution-v2"
python3 -B "$T/check_rans_balanced_evidence.py" --root "$B" --previous-root "$OLD"
python3 -B "$T/check_rans_balanced_diagnostics.py" --root "$B" --previous-root "$OLD"
for wall in EP Accepted; do
  python3 -B "$T/diagnose_rans_balanced_spatial.py" --root "$B" \
    --wall "$wall" --label "initial-$wall"
  python3 -B "$T/make_rans_balanced_meshes.py" --wall "$wall" \
    --output "$B/meshes/E-N-$wall"
  python3 -B "$T/diagnose_rans_balanced_spatial.py" --root "$B" \
    --wall "$wall" --label "balanced-$wall" \
    --specifications "$B/meshes/E-N-$wall/spaces.json"
  python3 -B "$T/classify_rans_balanced_spatial.py" --root "$B" \
    --wall "$wall" --output "$B/initial-gate-$wall.json"
done
python3 -B "$T/schedule_rans_balanced_evolution.py" --root "$B" --time .0016
python3 -B "$T/schedule_rans_balanced_evolution.py" --root "$B" --time .016
python3 -B "$T/qualify_rans_balanced_models.py" --root "$B"
```

Use `audit_rans_balanced_evolved.py --root "$B" --campaign-root "$B" --label
E-N-Accepted-B11 --output-label balanced-E-N-Accepted-B11 --analysis-runtime
"$B/runtime-resolution-v2"` for a representative complete .0016/.016 term
replay; repeat for both members and all six natural combinations with their
matching common directory. E defaults to `diagnosis/initial-WALL/common`;
A93/A95 use `--common-directory diagnosis/MODEL-Natural-WALL/common`.
The complete immutable command arrays, XMLs, artifact hashes and exact history
paths are in each campaign/request/receipt directory.

Frozen quadrature controls:

```bash
python3 -B "$T/check_rans_pilot_quadrature.py" --root "$B" \
  --coarse-suffix B5 --fine-suffix B11 --output-label balanced-quadrature
python3 -B "$T/check_rans_pilot_quadrature.py" --root "$B" \
  --coarse-suffix B11 --fine-suffix B5 --output-label balanced-fine-quadrature
python3 -B "$T/check_rans_balanced_state_quadrature.py" --root "$B" --time .016
python3 -B "$T/check_rans_balanced_evaluation.py" --root "$B"
python3 -B "$T/check_rans_common_test_sides.py" --root "$B"
python3 -B "$T/check_rans_observer_baseline.py" --root "$B"
```

After the scoped prerequisites and complete term replays pass, create the
immutable admission with `admit_rans_balanced_duration.py --root "$B"`. The
original experiment verifies its existing admission instead of replacing it.
The failed .004 cohort is a separate request and can be reproduced only into
its own fresh destination/label. Passing duration continuations use:

```bash
python3 -B "$T/run_rans_balanced_duration.py" --root "$B" \
  --model E94E95 --wall Accepted --dt .002 --max-time .016
python3 -B "$T/run_rans_balanced_duration.py" --root "$B" \
  --model E94E95 --wall Accepted --dt .002 --max-time .1 --resume
# .5 remains conditional on the .1 four-member gate.
python3 -B "$T/run_rans_balanced_duration.py" --root "$B" \
  --model E94E95 --wall Accepted --dt .002 --max-time .5 --resume
```

The controller validates an immutable current request before using completed
work, calls the native runner's exact-request validator even on completed
members, and never edits checkpoint signatures, dt or histories. Larger worker
allocations are not implied by these commands: at most four trajectory
controllers plus two serial diagnostic slots share the 72-worker total with
existing campaigns and numerical-library threads one.
