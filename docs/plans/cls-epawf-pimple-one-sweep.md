# CLS–EP-AWF PIMPLE single sweep

`RANSWallTemporalCoupling=PimpleSingleSweep` selects an opt-in candidate
first-order sequential, linearly implicit temporal coupling method. It uses
native IMEX1 and one implicit callback per physical step. It does not solve
the fully implicit backward-Euler stage to nonlinear convergence. The
existing default `ConvergedStage` and all nonlinear policies remain available.

This report records the uncorrected `RANSWallPimpleKBound=Off` baseline.
The later opt-in mass projection, its modified step map and continued
Linux experiments are documented in
[the focused k-bounding report](cls-epawf-pimple-k-bound.md).

The branch is `codex/cls-epawf-pimple-one-sweep`. The fetched remote base was
`0f2574b7f558963aafbc814b479c840e60880893`; the actual base was the newer local
`6156937022124575e67376799487accb9b271069`, containing the repaired pressure
verification and five subsequent development commits. Source editing moved
to the isolated Mac worktree `.worktrees/pimple-one-sweep` when concurrent
reduced-search edits appeared in the original checkout. Those edits were
preserved in the original checkout; the net PIMPLE branch diff leaves the
reduced-search implementation unchanged from the actual base.

## Scheduled operation

All four field updates have fraction one, so a fixed relaxation does not
rescale the leading physical-time increment. Equation under-relaxation, extra
pseudo-time mass, equilibrium-defect anchoring, nonlinear guess changes,
reduced/continuation/Newton/ground-projection searches and non-k equilibration
are absent. The selected native history is copied directly into the stage,
including any native explicit-history contribution. Prescribed body forces
are on the native explicit side, with zero prescribed body and scalar
forcing in the inherited outer6 experiment.

| Operation, once per step | Frozen coefficients and sources | Latest fields |
| --- | --- | --- |
| Momentum predictor | Oseen transport, symmetric-strain secant, full deferred constitutive stress, wall profiles/normalization and mechanical k/convection forcing from the committed state, externally lifted at stage time | Incoming native differential history; current separate regional pressures |
| Pressure correction | Same frozen full native G, both D blocks and gauge/pressure C as momentum, including roughness pressure forcing; consistent free-velocity mass | Momentum predictor |
| Outer epsilon | Transport, diffusion, signed production/Yap, implicit `Cepsilon2 * epsilon/k` reaction and projected trace target from the pre-epsilon state | Corrected flow; incoming k and epsilon |
| Full-domain k | Transport/diffusion/production, actual inner wall destruction rate with its zero-k extension; outer positive/negative production split and implicit destruction | Corrected flow and solved outer epsilon; incoming k |
| Refresh and diagnostics | No further transport/correction/restoration solve or trace overwrite | Retained flow, epsilon and k; recomputed wall profiles, normalization, trace target and original residual |
| Acceptance and commit | Strict R2LocalV1 or explicitly labelled R2DiagnosticV1 | Commit coupled state, native history and physical time once after all gates |

Velocity is a single vector equation: its two components are coupled in one
980-unknown predictor solve, and each component is updated once. No inner
epsilon transport equation exists. The outer epsilon/trace system has 254
unknowns and k has 504. Initialization may use the existing nonlinear `Lift`;
its calls/iterations are reported separately and must not increase during
physical advancement. Constitutive evaluations, boundary projections and
native diagnostic Riesz norms are evaluations, not transport solves.

## Approximate pressure correction

The production path factors A once and solves `A du* = -r_u` once. It
separately factors the actual consistent mass on the free velocity space,
uses `B = alpha^-1 M_free^-1`, and forms `S = D B G - C`. One direct,
generally nonsymmetric 85-unknown pressure/gauge solve gives
`S dz = D du* + r_a`; velocity receives the algebraic correction
`du = du* - B G dz`. The same stored BG columns enter S and the velocity
correction. There is no repeated full-momentum inverse in S construction.

For `A=alpha M+L`, this B is the leading term in the small-dt inverse
expansion; its absolute action error is O(dt²) for fixed L and pressure
action. Essential-boundary lifting, two pressure regions, relative regional
pressure offset and common gauge are retained. The identity checked is
`r_u + A du + G dz = (I-A B) G dz`. Equality with the exact frozen saddle
solution is deliberately not an acceptance condition.

`RANSWallPimpleVerify=On` runs a separately counted exact-flow reference and
four frozen-L pressure-action comparisons. Production counts exclude that
verification-only work: 86 full-momentum inverse applications and one
pressure solve for the retained exact correction, four additional momentum
action solves and one additional mass action. The production pressure
construction uses 85 mass-inverse applications, including the zero gauge
column; its final correction combines the stored columns algebraically.

## Matching, admissibility and rollback

The epsilon solve imposes the target evaluated at its assembly state, with
trace columns retained in every free transport row, including physical
mass/history. The following k update can change the nonlinear matching
target. Both the imposed-target and refreshed-target coefficient defects
are recorded. The solved epsilon coefficients survive unchanged.
`RoughnessBudget` and `EnergyBudget` use this retained state in the new mode;
no diagnostic calls `Lift` to restore the trace. Native original residual
and block norms are recorded after every completed scheduled sweep.

Numerical acceptance requires successful accurate scheduled linear solves,
finite fields, complete declared scalar/sampling inventory, continuous
cubic ground-k admissibility, boundaries, both continuity blocks, gauge,
covariance and normalization validity. There is no clipping, field floor,
residual-descent requirement, retry or nonlinear residual tolerance gate.
The original tolerance remains enforced for converged-stage Picard both in
its driver and in `RANSEquilibriumWallEvolution::Implicit`.

`R2LocalV1` stops before commit at its first qualification exit.
`R2DiagnosticV1`, available only for this mode, records the unchanged box's
first exit and may continue numerically while retaining all hard validity
checks. Such states are not physically qualified M8 evolution. It does not
change any qualification bound. A failed callback restores every native
history slot bitwise; candidate coupled state stays local until acceptance.
The deliberate rejection control acts after scheduled solves, before commit.

## Linux evidence

All source edits and Git commits were made on the Mac. The shared Linux
checkout at `f4b70baec10016dd098005bc5a5507b87e4a57bf`, its pre-existing
working-tree changes, installation and five ongoing solver experiments were
left untouched. `linux-usb` timed out; Tailscale SSH to
`jay@100.105.247.110` succeeded. Heavy files, isolated source/build, all
attempts, checkpoints, native reports, binary/library hashes and logs are in:

`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-pimple-one-sweep`

The build uses Release, MPI-enabled installed dependencies, eight explicit
compiler jobs, `-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`.
Numerical runs use one rank and one thread in each OpenMP/BLAS runtime;
CTest concurrency is one. Build and test allocations stayed below 72
workers while the existing five numerical experiments ran independently.
The production executable is never replaced during an experiment.

The unchanged outer6 problem has 2,165 coupled unknowns, full-domain
velocity/k, outer epsilon, p3 velocity/scalars, p1 regional pressures,
source rules 8/16 and native GLL8. Regularized CLS, SaturatedR2 inner stress,
capped destruction and its inherited smoothing, molecular viscosity,
transition/normalization laws, mesh, boundary data and all zero prescribed
loads are retained. The actual represented perturbed h0 has SHA256
`b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71`
and initial native k residual norm `62.60959506799984`. The fixed reference
is a different state; it is not substituted for h0.

The first strict qualification exit occurs on the attempted first step at
`t=.0025`, gradient departure `.002357770290111416 > .001`, witness
`(x,z)=(1.6046496089512394,2)`, component `(0,2)`, outer geometry 34.
Both covariance minima are positive (`.03293658`, `.17689637`). The strict
run commits zero physical steps and keeps time/history at zero. The legacy
error text says R1, but the selected guard and complete audit are R2LocalV1.

The diagnostic dt=.0025 trajectory completes two steps to .005. This does
not establish a converged trajectory at that horizon: dt=.00125 rejects
its third attempted step at .00375. Its solved k polynomial has continuous
ground minimum `-.01770974049885598` at `x=1.6022028132`, geometry 4. A
prescribed-source evaluation at `z=4.8478017123e-7` is also negative
(`-.017417516873183554`). Frozen k and epsilon linear residuals are
`1.11e-16` and `1.95e-18`, respectively. Thus the first observed hard failure
is scalar inadmissibility of an accurately solved update, not insufficient
linear convergence. The dt=.000625 run similarly rejects its fourth
attempt at .0025, retaining three accepted steps to .001875. No field
increment is reduced and no mesh/source modification is made. Refinement
therefore uses a shorter common physical horizon.

The final build passes seven focused CTests (40.40 s total), including the
six existing tests and the new `RANSHistoryRejectionChecks`. The latter uses nonzero explicit
forcing, checks IMEX1 and IMEX2 rollback bitwise in every history slot and
time label, then compares resumed and uninterrupted evolution. The retained
Picard verification passes: verified/unverified two-sweep states and norms
are identical, while both runs still reject the physical stage at the old
nonlinear convergence limit. The constitutive kernel passes 28,434 checks,
including the defined zero-k rate extension and signed source handling.

The first common-time study at `.00125`, using dt `.000625`, `.0003125`
and `.00015625` (2/4/8 steps), did **not** show convergence. Successive
velocity L2 differences were `2.16160e-4` and `2.15405e-4`; epsilon and
pressure differences increased. This failed comparison is retained.
The refreshed trace defect nevertheless decreased
`3.07826e-6 -> 1.53666e-6 -> 7.49085e-7`. There is no observed amplification
of that trace defect in this comparison. On the finest trajectory, ground
k decreased to `.0510055` and minimum B changed rapidly to `.950578`,
where the coarsest endpoint still had ground k `.165387` and B=1. This
associates the poor comparison with rapid near-ground/source-normalization
evolution; it does not establish a unique causal explanation. No linearized
trace variant or nonlinear restoration was introduced on that evidence.
An earlier common-time study at `.000625` tests the unchanged baseline
before that rapid normalization change.

## Final comparison and checks (November 6, 2026)

The final tested implementation/test revision is
`7baa9f9f274e0ede7bf953a3e7a7adfe835d88f3`. The main implementation is
`76844da`; the independent worktree cleanup/compile fix is `63418a5`,
focused experiment `fcbee66`, history regression `5879e63`, harness fix
`dc2f6c5`, complete coupled rollback `69f323d`, and final restart/horizon
harness `7baa9f9`. Tracked producing sources are clean. Some intermediate
runner manifests record an untracked Python bytecode cache; it was removed
after testing. Every changed source and executable/library identity is bound
in [compact native evidence](cls-epawf-pimple-one-sweep-evidence/native.json).
The earlier `.00125` comparison records its own producing revision `dc2f6c5`
and binary hash. An accepted-state prefix from that epoch and the final
rollback-hardened executable is bitwise identical.

The final solver SHA256 is `cf592ab290f82486df4a3e582ddcd5c8d811690a8cb4c0d1a0445376bf0f1193`.

These are independent fixed-dt, zero-forcing runs from the same h0 to
`t=.000625`. Norms are physical quadrature L2 norms on the retained fields.
Pressures are aligned by one common weighted mean across the two regions;
the regional pressure offset remains in the difference. Ground minima are
continuous polynomial minima, not coefficient signs. Weak continuity/gauge
and strong divergence are reported separately.

| dt | Accepted steps | Ground k min | Outer k min | Outer epsilon min | Max weak continuity/gauge | Refreshed trace defect | Solver wall seconds | R2 qualified at final time |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 0.00015625 | 4 | 0.236459269 | 0.510507393 | 0.0363179166 | 1.042e-16 | 7.560736e-07 | 100.82 | no |
| 7.8125e-05 | 8 | 0.214026843 | 0.510507391 | 0.0363179166 | 1.311e-16 | 3.771311e-07 | 182.97 | no |
| 3.90625e-05 | 16 | 0.200292708 | 0.51050739 | 0.0363179166 | 9.197e-17 | 1.883226e-07 | 348.39 | no |

All three final states have B=lambda=1, and their imposed trace defects are
at most `2.12e-20`. Their first R2 qualification exit is at `.000625`.
The inner/outer strong-divergence L2 norms are respectively
`(2.33461e-3,1.08488e-5)`, `(2.39303e-3,1.07711e-5)` and
`(2.42203e-3,1.07646e-5)`. The lower-degree pressure tests enforce weak
continuity; pointwise divergence is a spatial diagnostic. The fully
refreshed original stage maxima are `8.63035`, `5.13719`, `2.86797`.
They are recorded, not used as nonlinear convergence gates. Wall seconds
include initialization, the complete acceptance inventories and output on
the NTFS run volume; these runs are not a performance benchmark.

| Physical L2 difference | dt vs dt/2 | dt/2 vs dt/4 | Observed difference order |
| --- | ---: | ---: | ---: |
| Velocity | 3.88799230e-05 | 2.08650194e-05 | 0.898 |
| k | 1.06278139e-03 | 6.45200394e-04 | 0.720 |
| Outer epsilon | 7.61263348e-08 | 3.86139156e-08 | 0.979 |
| Pressure, common gauge | 9.04593911e-03 | 6.94062769e-03 | 0.382 |
| Ground shear stress | 1.58265682e-02 | 1.01433492e-02 | 0.642 |
| Wall sampled K | 1.23830715e-08 | 9.21989369e-09 | 0.426 |

Wall B and lambda differences are zero in this early comparison. All
measured evolving-field differences decrease here. The rates, especially
pressure, do **not** establish uniform first-order accuracy. This is
short-interval empirical convergence of a candidate first-order method,
not a convergence claim at `.00125` or `.005`.

A separate strict run at dt `.00015625` accepts three R2-qualified steps to
`.00046875`, then rejects its fourth attempt at `.000625` for the unchanged
gradient qualification bound. Every retained state in that strict run is
qualified. The diagnostic endpoints in the tables are outside that narrow
reference box and are not physically qualified M8 endpoints.

Each accepted step executes `[outer, momentum, pressure, epsilon, k,
implicit callback] = [1,1,1,1,1,1]`. Additional correction/restoration solves
are zero. Pressure construction applies the consistent free-velocity mass
inverse 85 times; no momentum inverse enters its pressure matrix. All
production runs retain initialization work `[Lift calls, iterations]=[2,3]`
with no later increase. The exact-flow comparisons occur only in the
separate two-step verification trajectory; verified and unverified states
are identical.

The first final-epoch pressure verification gives:

| Check | Result |
| --- | ---: |
| Essential boundary defect | 0.00000000e+00 |
| Pressure linear residual | 8.41805166e-20 |
| Frozen momentum splitting-identity error | 4.01154804e-17 |
| Native full G action error | 6.07153217e-18 |
| Difference from exact frozen flow, max coefficient | 8.11295247e-05 |

The nonzero exact-flow difference is expected. For a fixed representative
native pressure action and fixed L, the verification varies alpha without
advancing the trajectory:

| 1/alpha | Absolute mass-norm inverse-action error | Relative error |
| ---: | ---: | ---: |
| 0.00015625 | 1.23537477e-07 | 3.18896415e-02 |
| 3.90625e-05 | 8.73067243e-09 | 8.95175769e-03 |
| 9.765625e-06 | 5.65514108e-10 | 2.31467231e-03 |
| 2.44140625e-06 | 3.56758426e-11 | 5.83785077e-04 |

The absolute error approaches quadratic decay with dt, as expected for
this leading mass inverse. The same-state scalar source/storage comparison
has maximum coefficient defect below `9e-16` on the verification steps.
Their frozen scalar linear residuals pass the scheduled solve tolerances.
Matching-trace mass-column magnitude reaches `355.5555556`; on the second
verification step one retained trace-column physical mass/history term
reaches `1.2961241e-4` in a free epsilon transport row. Removing those columns would
change the equation actually solved.

Checkpoint/restart runs two steps, saves, and resumes two more. Final
coupled state, physical time and every native history slot/time label are
bitwise identical to the uninterrupted four-step run. A separate run
accepts two steps, deliberately rejects its third after all scheduled
solves, and is bitwise identical to the two-step checkpoint state/history
at time `.0003125`. Its three callbacks and three attempted sweeps are
reported separately from its two accepted steps.

| Outcome | Result |
| --- | --- |
| Implementation and Linux build | Passed |
| Frozen-operator, scalar and pressure checks | Passed |
| Single-sweep numerical evolution | Completed 4/8/16-step common-time runs |
| Time-step convergence | Decreasing physical differences at `.000625`; failed earlier comparison at `.00125` retained |
| Uniform first-order accuracy | Not demonstrated by these rates |
| Original constitutive qualification | Three-step strict prefix passes; diagnostic comparison endpoints fail the unchanged box |
| Original fully converged M8 stage | Not solved; not this temporal method's acceptance contract |

The retained initial compiler namespace error, test link error and Python
provenance/ASLR/restart assertion fixes are development/infrastructure
failures, separately recorded from the numerical k rejections and
qualification exits. Some serial MPI executable logs contain PMIx dstore
permission messages at teardown on the NTFS temporary volume. Completed
runs have zero process exit status, complete native receipts and passed
readback/checks. No distributed-MPI claim follows from these serial tests.

## Reproduction

The following uses the retained clean producing source at revision
`7baa9f9f274e0ede7bf953a3e7a7adfe835d88f3`; use fresh output directories.
The compact evidence contains exact commands, all case XML hashes,
executable/library hashes and per-attempt solve counts. Raw fields,
checkpoints, complete state inventories and readbacks stay under the Linux
run root. Nothing is installed over the shared dependency prefix.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
run='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-pimple-one-sweep'
src="$run/source"
export TMPDIR="$run/tmp"
cmake -S "$src/solvers/IncNavierStokesSolver/Tests/RANS" -B "$run/build" \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/home/jay/.local/opt/nektar-v2026-mpi \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
  -DCMAKE_CXX_FLAGS='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$run/build" --parallel 8 --target IncNavierStokesSolver \
  RANSWallResolutionProbe RANSEquilibriumTimeProbe RANSEquilibriumEvolutionProbe \
  RANSEquilibriumPatchProbe RANSReducedModelChecks RANSHistoryRejectionChecks \
  RANSRobustVariantProbe
ctest --test-dir "$run/build" --parallel 1 --output-on-failure \
  -R '^(RANSReducedModelChecks|RANSHistoryRejectionChecks|RANSEquilibriumTimeNative|RANSEquilibriumEvolutionNative|RANSEquilibriumPatchNative|RANSEquilibriumWallEvolutionSmoke|RANSEquilibriumWallStageSmoke)$'
python3 -B "$src/solvers/IncNavierStokesSolver/Tests/RANS/test_cls_epawf_pimple.py" \
  --solver "$run/build/IncNavierStokesSolver" \
  --resolution-probe "$run/build/RANSWallResolutionProbe" \
  --comparison-dt .00015625 --comparison-steps 4 --output "$run/replay-short"
python3 -B "$src/solvers/IncNavierStokesSolver/Tests/RANS/run_cls_epawf_pimple.py" \
  --solver "$run/build/IncNavierStokesSolver" --dt .00015625 --steps 4 \
  --output "$run/replay-strict-early"
python3 -B "$src/solvers/IncNavierStokesSolver/Tests/RANS/test_cls_epawf_picard_verification.py" \
  "$run/build/IncNavierStokesSolver" "$run/replay-old-picard"
"$run/build/RANSRobustVariantProbe" --output "$run/replay-kernel.json"
```

The original requested dt comparison and its retained admissibility failures
can be reproduced with the dedicated experiment runner's default
`--comparison-dt .0025 --comparison-steps 2`, using a fresh output. To repeat
the nonconvergent `.00125` common-time comparison directly, use
`--comparison-dt .000625 --comparison-steps 2`. These are independent runs
from h0, not adaptive retries inside the integrator.
