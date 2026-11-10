# Neutral persistence pilot

Starting revision: `1297301156c703fd1e4fe3b7776dc3676a5f3267`, clean existing
implementation worktree. The original Mac checkout and Linux primary checkout
retain their unrelated work. Two active CLS controllers and two native serial
campaigns are read-only. New Linux root:
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-pilot-20261110`.
Separate source, build, dependency, runtime, tmp and case directories are used.
Build allocation starts at eight, later four; each native run has one MPI rank
and numerical library thread counts one. Recount other work before concurrency;
never exceed 72 aggregate workers.

The existing specification hashes and mathematical contract are unchanged.
The earlier preservation receipt, including eight unresolved epsilon orders,
is immutable. No retrospective threshold change or reclassification is allowed.

## Predeclared experiment gates

FrozenReferenceV1 is opt-in, fingerprinted, and reconstructs live external data
from validated stored reference coefficients in both data objects before
reference restoration and initial residual receipts. Legacy signatures remain
byte-compatible. Strict frozen comparisons require matching raw/effective
initial coefficient, reference and live boundary hashes. No evolved signature
or history is edited. Interface requests at H require an explicit inner/outer
side; observations carry side and native element provenance.

Campaigns default to a new destination. Before any mutable receipt, an immutable
request records physics, runtime/dependencies, spaces, integration, initial and
boundary identities, perturbation and evaluation grid/sides/weights. Resume is
explicit and exact-request only; completed artifacts are hash-verified.

Select E94E95 natural EP and accepted spaces using the existing h2/p3 and h0/p5
results; first test h0/p5 and a combined h1/p5 comparator. If needed, one targeted
additional refinement is permitted with its reason recorded. Pressure degree
is velocity degree minus two. Meshes stay fixed. Boundary data prescribe the
neutral reference on the finite patch; this is not an open-outflow ABL test.

For every pilot use fixed physical scales u*=.4, H=2, K0=q*u*^2,
epsilon(H)=u*^3/(kappa*H), pressure/traction=u*^2, k flux=K0*u*.
Physical errors E compare to the analytical reference; D compares to the initial
numerical field. Both use weighted RMS; true L2 and one global Pi gauge are
retained. Targets: E_u,E_k,E_epsilon,E_Pi <=1e-3; RMS vertical velocity/u* <=1e-4;
mean wall-traction error and resolved interface traction RMS/u*^2 <=1e-3.
Never divide by a vanishing local reference.

A material spatial reduction means fine/coarse <=0.85 for non-floor errors
above 1e-10, or classification as unresolved. Same-time spatial comparisons
must be consistent with the reference-error triangle bound. A targeted dt/2
control must be <=10% of the interpreted spatial difference, with absolute
scaled floor 1e-7; otherwise the evolution gate is temporally unresolved and
requires a smaller dt. These new pilot uncertainty thresholds do not change
the prior epsilon temporal-order receipt.

Advance satisfactory matched pairs first to T=.0016, then .016. Extend passing
selected spaces to .05 and .1 by genuine same-method restart. theta=T*u*/H is
.0032,.01,.02 at the latter milestones. Stop only the affected case on the
first sampled hard-domain or declared accuracy failure, preserving state and
history; continue independent cases. Retain accepted/stage extrema, actual
radicands, cutoff applicability, frozen trace and refreshed lag, source-storage
identities, conservative budgets, correction calls/integrals and scheduled
alpha=0 residuals. Corrections must remain zero. Storage identities <=1e-9 and
continuity/gauge <=1e-10 retain existing scaled tolerances. No claim that every
transient residual decreases is made.

Profile a refined case before optimizing. Any diagnostic partition/mass-factor
cache is scoped to immutable mass/space and exact free-ID lists and must give
identical diagnostics and state/history with caching disabled.

After the neutral evolved gate passes, test a small divergence-free analytic
streamfunction perturbation on at least A95 or A93. It vanishes at external
boundaries and preserves interface traces. Compare against the unperturbed
numerical trajectory at identical times; measure actual discrete divergence,
pressure and normal-stress response. No monotone-decay or global-stability
claim. Include a two-column streamwise check before extending the single-column
patch claim. All longer/curved/3D/MPI/nonlinear/physical qualification is outside
this package.

## Results and provenance

In progress. Gates, coverage, first failures and the next justified experiment
will be completed after the bounded evolution work.

### First p5 advance and diagnostic extension

The first E94E95 EP h0/p5 advance failed at time zero with the inherited
`WallTraceR4` cubic-only basis guard, before any governing solve. The frozen
p studies did not exercise this stage check. Its raw report, initial field
errors and native history are retained in `campaigns/E-N-EP-h0p5` under the new
Linux root. This is an architectural qualification gap, not a continuum failure.

The new-model p4/p5 continuous ground check uses the exact Modified_A polynomial
and all stationary points, isolated recursively using Rolle partitions, including
repeated roots. Endpoints are retained. Arithmetic uncertainty scales with the
coefficient norm. The cubic CLS implementation is unchanged; other unsupported
degrees remain rejected. Independent known-root and narrow-negative-well tests
check this extension before a new runtime is frozen. The failed campaign is
never resumed with changed source/runtime identities.

### Evolved space selection, before extended pilots

E94E95 natural EP and accepted h0/p4 versus h0/p5 pass the evolved-space gate
at T=.0016 with dt=.0004. Field-error fine/coarse ratios are .126--.402;
both independent dt=.0002 controls satisfy the predeclared uncertainty gate.
The p4 continuations first exceed reconstructed interface traction tolerance
at T=.0028 (last committed T=.0024). Field errors remain within targets; this
is a pointwise numerical-accuracy failure, not loss of weak conservation or
constitutive admissibility. Both full candidate states/histories are retained.

Full h1/p5 would exceed the existing 3000-unknown dense-backend guard and remains
rejected. The next targeted comparator is `Interface8V1`: split the six inner
intervals [.02,.04],[.04,.1],[.1,.25],[.25,.5],[.5,1],[1,2] and outer intervals
[2,2.25],[2.5,3]. Thus both neighbors of H are refined, with fixed mesh and the
same continuum boundary specification. Base p5 remains the candidate, and this
partial h refinement is the genuinely finer comparator. No backend or accuracy
limit is increased. Compare them and a base-p5 dt/2 control at .0016 and .016
before advancing .05/.1.

The profiled base p5 four-step jobs cost 6.8/8.5 seconds per requested step
including initialization, with exclusive domain computation 8.0/11.5 seconds
across 16 stage checks. The entire physical audit costs only about 1.4 seconds
per job; restricted-mass caching cannot remove the dominant mandatory domain
work. No diagnostic optimization is introduced without a measured benefit.

### Boundary specification and runtime identities

The sand-grain input is h_r=.1, not a separately fixed aerodynamic z0. The
original roughness-intercept correlation is retained. Each model/wall/policy
has U_H=u_star^2 J0(H) and outer U=U_H+u_star/kappa log(z/H), equivalently
z0_m=H exp(-kappa*U_H/u_star), as allowed by the appendix's roughness compatibility
condition. All compared spaces within a configuration use that same continuum
reference. Mesh-dependent prescribed traces are its documented constrained
polynomial representations. Both streamwise sides and the top prescribe the
reference velocity/k/outer-epsilon traces; the ground has no slip and natural
k flux. Pressure has one global gauge. This bounded patch cannot demonstrate
open-outflow downstream ABL preservation.

The first new runtime (source 8b91714, C++ content compiled at 4b76bc0) has solver
SHA256 `411a6eb88e3b5aba89a589e687ce565fce8f7550b5c38cf27312ebd833543e25`.
It records the cubic-only p5 failure. The qualified p4/p5 runtime is separately
frozen under `p5/runtime`, source snapshot 99c6bf0 with identical C++ content to
compiled ce11c23, solver SHA256
`6f16a3366f67d2964f5a2ebb395d4f2f06f025058ab6679c24ac5f304dc08950`.
Later changes are Python verification/orchestration and documentation only.
Full source-file hashes, compiler flags, linked libraries and 122 actually
loaded shared libraries/plugins are recorded. All numerical library thread
counts are one. The eight unresolved historical epsilon observations and the
historical overall unresolved receipt remain unchanged.

### Verification receipts already completed

* FrozenReferenceV1 accepted Q0/Q2 controls: all three models pass with identical
  actual raw/effective/reference and explicit/implicit boundary hashes, zero
  ExternalLift difference, native fresh histories, and matched final replay maps.
  Descriptor, identity, nonzero time, unequal state/reference and boundary-map
  negative controls reject for their intended reasons.
* One-sided A95 interface values, different gradients, regional pressures and
  independent full tractions pass. An ambiguous H observation rejects.
* The immutable-request unit suite passes. CLI changes to runtime, quadrature,
  perturbation, horizon or boundary identity reject without modifying the
  aggregate receipt. Identical resume completes only the two missing native
  cases and preserves the first case's sealed evidence. A two-point/T=.0001
  reuse fixture fails its Pi temporal-order check and correctly exits 1; that
  numerical failure is retained separately from verified reuse behavior.
* The targeted CTest selection passes 4/4; extended spectrum/domain/ground-root
  checks pass 724 assertions. The 22 affected restart/cache/rollback,
  calibration/model rejection and unsupported-route controls pass, including
  bitwise continuation of existing ordinary new-model histories.

Base-p5 and Interface8V1 pass the E94E95 natural evolved-space and dt/2 uncertainty
gates at .0016 and .016. At .016, base/fine velocity errors are approximately
2.1e-5/7.4e-7 on u_star; both walls pass every field gate. Both base spaces and
the finer accepted-wall case reach .05 with zero bound/correction calls. These
are bounded accuracy observations, not measured temporal orders or long-time
stability results. The subsequent .1 and remaining-model results follow below.

## Reproduction commands

Use a fresh Linux destination, with a clean source snapshot and available
worker budget. The original `p5/runtime` and every campaign/checkpoint are
immutable. `TESTED_HEAD` below means the final branch revision reported with
this document; its governing C++ hashes match compiled ce11c23.

```bash
export REPRO=/media/jay/新加卷/nektar-v2026-runs/a93-pilot-reproduction
export ORIGINAL=/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-pilot-20261110
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
mkdir -p "$REPRO/dependency" "$REPRO/tmp"
git -C /media/jay/新加卷/nektar-v2026 worktree add --detach "$REPRO/source" TESTED_HEAD
cp "$ORIGINAL/p5/runtime/lib/libSpatialDomains.so.5.9.0" "$REPRO/dependency/"
export TMPDIR="$REPRO/tmp"
export T="$REPRO/source/solvers/IncNavierStokesSolver/Tests/RANS"
cmake -S "$T" -B "$REPRO/build" \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DRANS_PREBUILT_SPATIAL_LIBRARY="$REPRO/dependency/libSpatialDomains.so.5.9.0" \
  '-DCMAKE_CXX_FLAGS=-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$REPRO/build" --parallel 8 --target IncNavierStokesSolver \
  RANSPreservationProbe RANSNewModelsProbe RANSNativeDiagnosticChecks
python3 -B "$T/freeze_rans_pilot_runtime.py" --source "$REPRO/source" \
  --build "$REPRO/build" --output "$REPRO/runtime"
python3 -B "$T/check_rans_frozen_initial_history.py" --root "$REPRO"
python3 -B "$T/check_rans_interface_sides.py" --root "$REPRO"
python3 -B "$T/schedule_rans_neutral_pilot.py" --root "$REPRO" --concurrency 4
```

The scheduler validates immutable requests and staged comparisons, restricts
numerical libraries to one thread, accounts for existing Nektar/Python/compiler
jobs, and reserves ten worker threads before starting another controller. A
separate native control/probe slot must share the same 72-worker total. It
never modifies other controllers. A failed/unresolved gate is retained and
prevents only that continuation.

A single admissible reference and genuine restart continuation can be reproduced
without the full scheduler:

```bash
python3 -B "$T/run_rans_neutral_pilot.py" --root "$REPRO" \
  --label E-natural-accepted-example --model E94E95 --wall Accepted --dt .0004
python3 -B "$T/run_rans_neutral_pilot.py" --root "$REPRO" \
  --label E-natural-accepted-example --model E94E95 --wall Accepted --dt .0004 \
  --resume --max-time .016
```

Use `--refinement Interface8V1` for the finer 2793-unknown comparator, and a new
label with `--dt .0002` for its independent-history temporal control. Keep the
same dt when continuing a checkpoint. The `--max-time` operational limit does
not change the immutable requested milestone set (.0016,.016,.05,.1).

The velocity control uses `--model A95 --wall EP --velocity-amplitude .00002`;
its actual zero-step native initialization and history hashes are bound before
advancement. `analyze_rans_velocity_pressure.py` compares it to the unperturbed
numerical trajectory at identical times. Two-column p4 sensitivity uses
`--degree 4 --columns 2`; comparisons are resampled on x panels [-1,0],[0,1].
The selected coefficient-policy/cache/CLS controls are reproduced by
`run_rans_pilot_controls.py --root "$REPRO" --cls-baseline
/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-awf-20261109/cls-reviewed-whole4`.
Exact inputs, commands, request hashes and source/native-history locations are
retained in the evidence index.
