# Segregated Picard for the CLS/EP-AWF physical stage

The segregated backend completed Linux operator and history checks and
ran the unchanged M8 startup stage. The original convergence problem
remains open. The default run stopped at a full residual of 13.69368.
A k update projected against continuous ground constraints reduced the
residual further, with the k block still dominant. No Picard physical
step was accepted.

The branch is `codex/cls-epawf-segregated-picard`. Its base
`0ed7790450bb7e81f17ff1d0b577f12fc95f9cce` matches the earlier reviewed tip.
The implementation and harness are in `0486db8c39ce325bc0c2d004382bfc8881f6ea07`.
The ground projection is in `d3bc3f5f516b7ecb957debac52ed1ac9001f47bd`.
The Mac checkout had only unrelated untracked files at startup.

## Discrete stage

The state uses full-domain continuous velocity and k, outer continuous
epsilon, two separate continuous regional reduced pressures, and one
common gauge multiplier. For the archived outer6 case the respective
coefficient counts are 592, 592, 592, 304, 42, 42 and 1. Velocity and k
use p3 and pressure uses p1. The source rules are 8/16 and the native
expansions use GLL8. Epsilon has independent outer coefficients on the
matching trace. Inner dissipation is reconstructed from the profiles and
complete-ray normalization. Shared velocity and k modes supply the weak
interface flux and traction coupling, with separate regional pressure
constants and one combined pressure mean constraint.

The implicit callback supplies alpha, time and the actual IMEX history.
The stage uses physical mass applied to current lifted differential
coefficients minus the incoming differential history. Pressure and gauge
receive algebraic rows. State-dependent epsilon trace coefficients
participate in the transport rows. Matching trace rows enforce their
difference from the reconstructed and projected target. For the target,
eta is zero and all prescribed loads are zero. The original `Residual`
and `Norms` compute the final mass-dual block norms and scaled constraint
defect, with the historical threshold of `2e-12` on their maximum.

## Sweep

1. Refresh the lagged profile, quadrature, cap normalization, closure and
   forcing. The velocity predictor uses conservative Oseen transport and
   the symmetric-strain operator. Its inner viscosity is the molecular
   viscosity plus the selected inner stress secant. Its outer viscosity is
   the CLS eddy viscosity plus the prescribed regularization coefficient.
   The complete stress at the lagged state supplies the deferred load,
   including the nonlinear tensor remainder and the reference-strain
   contribution. The pressure gradient part of the roughness forcing
   enters the pressure matrix. Convection and `-2 grad(k)/3` are lagged
   within this subsolve.
2. With the velocity factorization, form `D A^-1 G - C`, where `C` is the
   augmented pressure/gauge coupling. Solve for the pressure increment and
   apply the velocity correction with the same factorization. Both
   continuity blocks and the common gauge retain their original rows.
   All solves use increments on free coefficients, so prescribed values
   and their contributions are present in the initial defect. Damping
   applies the same fraction to velocity and pressure together and
   preserves the linear constraints between feasible endpoints.
3. Assemble k using the corrected flow. Freeze advection, diffusivity,
   signed production and complete-ray inputs. Treat inner destruction with
   the model's point `rate` and its zero-k limit. For outer k use
   `(epsilon - min(P,0))/k` as the lagged reaction and `max(P,0)` as
   production on the source side. Refresh the wall profiles, normalization
   and mechanical pressure forcing after the k update.
4. Assemble outer epsilon with the newest k. Freeze diffusivity, advection,
   production, Yap and the reconstructed matching target for this
   subsolve. Implicit destruction uses `Cepsilon2 * epsilon_lag/k_latest`.
   The free transport rows include the matching trace columns. The trace
   rows move toward the current projected target, which is recomputed from
   the resulting fields. A transient trace mismatch is permitted after k
   and epsilon updates.
5. Recompute the complete native residual. The squared sum of scaled block
   norms selects a sweep or triggers backtracking. The maximum original
   block norm determines convergence. The existing native qualification
   runs before publication of the accepted state or advancement of IMEX
   history.

The frozen defect plus block matrix times increment is the affine equation
solved by each block. At zero increment the volume, lifting, mass/history,
gauge and matching loads recover the original residual. In particular
`rate * k_lag` is the point destruction and the outer split recovers the
signed native sources. `RANSWallPicardVerify=On` compares the separately
assembled same-state blocks with the native residual at initialization and
on an admissible nonuniform perturbation. On the first sweep it compares
the Schur correction with the existing direct dense solve of the frozen
saddle system, checks the affine pressure action independently, and records
native predictor/corrector continuity, gauge and interface observations.

## Configuration

Set `RANSWallNonlinearPolicy=SegregatedPicard` in the native evolution
session. The initial scope is flat, 2D, cubic modified-basis ground traces,
one rank, uncondensed upstream CLS with capped destruction. The backend
uses the existing 3000-unknown guard. It accepts the actual native stage
arguments across the IMEX startup and subsequent stages. The `Legacy`
default and historical stage-test policies remain selectable.

| Session parameter | Default | Role |
| --- | ---: | --- |
| `RANSWallPicardFlowRelaxation` | 1 | Joint velocity and pressure fraction |
| `RANSWallPicardKRelaxation` | 0.5 | k fraction |
| `RANSWallPicardEpsilonRelaxation` | 0.5 | Epsilon transport and trace fraction |
| `RANSWallPicardPseudoAlpha` | 0 | Optional coefficient multiplying mass times differential increment |
| `RANSWallPicardGuessPerturbation` | 0 | Signed deterministic guess-only perturbation for comparisons |
| `RANSWallMaxIterations` | 20 | Completed-sweep budget, with zero available for a forced rejection |
| `RANSWallNonlinearTolerance` | 2e-12 | Original full residual threshold |

Relaxation lies in (0,1]. The pseudo-time term vanishes at the fixed point
and has no pressure term. The implementation backtracks scalar updates
against the declared finite scalar inventory, including the sampled K and
moving normalization rays and the projected trace, and checks the
continuous cubic ground k trace. These checks cover the existing scalar
evaluation domain plus the ground trace. The fixture applies the original
`R2LocalV1` constitutive qualification at native stage acceptance.
Trial sweeps can leave that constitutive envelope. No clipping or field
floors are applied.

`RANSWallPicardGroundProjection=On` adds an optional search policy for k.
If a relaxed k solve crosses the continuous ground constraint, the solver
projects its increment against cuts at the trial polynomial extrema.
It minimizes the change in frozen mass-dual k residual and uses the
free-coefficient metric inverse `A^-1 M A^-T`. The projected increment
retains all original k equations for the acceptance test. Other scalar
inventories and the full residual still gate the sweep. The default
remains backtracking, with the projection selector off.

## Reproduction on Linux

Inspect the load and the shared source revision and status via `ssh
linux-usb` or the documented Tailscale address. Transfer the committed
branch to a fresh Linux source checkout. Verify the commit and clean
working tree there before building. The example uses the existing
MPI-enabled dependency installation and places both fresh build products
and all test artifacts outside the shared checkout.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
run='/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard'
src="$run/validation-source"
prefix=/home/jay/.local/opt/nektar-v2026-mpi
mkdir -p "$run/tmp"
export TMPDIR="$run/tmp"
revision=$(git -C "$src" rev-parse HEAD)
python3 "$src/docs/plans/cls-epawf-revised-evidence/m8-r3/build-epoch.py" \
  --source "$src" --build "$run/build" --prefix "$prefix" \
  --output "$run/build-receipt" --expected-source "$revision" --jobs 8 \
  --target IncNavierStokesSolver --target RANSWallResolutionProbe \
  --target RANSRobustVariantProbe --target RANSEquilibriumTimeProbe \
  --target RANSEquilibriumEvolutionProbe --target RANSEquilibriumPatchProbe
ctest --test-dir "$run/build" --parallel 1 --output-on-failure \
  -R '^(RANSEquilibriumTimeNative|RANSEquilibriumEvolutionNative|RANSEquilibriumPatchNative|RANSEquilibriumWallEvolutionSmoke|RANSEquilibriumWallStageSmoke)$'
"$run/build/RANSRobustVariantProbe" --output "$run/robust-kernel.json"
python3 "$src/solvers/IncNavierStokesSolver/Tests/RANS/run_cls_epawf_picard.py" \
  --solver "$run/build/IncNavierStokesSolver" \
  --resolution-probe "$run/build/RANSWallResolutionProbe" \
  --output "$run/unchanged-first-stage" --verify --sweeps 120
```

Run the build and tests sequentially, using eight compiler jobs only when
the measured server workload leaves room below the 72-worker cap. The
native harness uses one rank and one thread per configured runtime. It
binds binary and linked libraries by path and SHA256. The XML is generated
from the tracked unforced session via the original R2 outer6 generator.
The harness checks the startup coefficient hash, incoming guess, actual
stage history values, alpha 400 and time .0025. It records signed-zero
history differences separately via the raw and normalized-zero hashes.
The deliberate zero-sweep rejection compares all IMEX history slots and
time labels with the zero-step startup. The readback probe rebuilds the
original native stage residual from saved coefficient and history data
in a separate process. On a reported first-stage failure it can also
read the last feasible trial for the existing k-term decomposition.

For the projection comparisons, build the solver from the separate
`projection-source` into `projection-build` with the same Release flags.
Use that solver and the first build's `RANSWallResolutionProbe` in the
harness. Set `--ground-projection --sweeps 120` and use a fresh output
directory for each relaxation setting. The readback probe uses the original
residual, whose assembly is shared by both producing revisions.

After a qualified first stage, run with `--steps 3` in a fresh output
directory for the focused multi-step fixture through .0075. Use
`--perturbation 1e-7`, `--perturbation -1e-7` and
`--perturbation 1e-6` in separate first-stage runs. A second relaxation
setting can use `--k .25 --epsilon .25` while keeping the same physical
target. Choose later continuation budgets and pseudo-time settings from
measured sweep, block and readback evidence. The serial Picard path has
no MPI comparison at present.

## Linux results on November 4, 2026

Tailscale access recovered after the initial connection failures. The clean
shared Linux checkout was on `codex/hit-inflow-generator` at
`f4b70baec10016dd098005bc5a5507b87e4a57bf`. The baseline producing source
at `f6e8ff839d28247bd5c5f9b7195493ffbdd67543` was checked out in
`validation-source` using the retained R7 epoch-002 objects and a thin Git
bundle from the Mac. The projection source at `d3bc3f5f516b7ecb957debac52ed1ac9001f47bd`
was checked out separately. Both Release builds used eight jobs against
the existing MPI-enabled prefix. The numerical runs used one rank and one
thread per configured runtime. The build receipts record the exact compile
flags, commands and compiler inputs. The runtime manifests record paths,
hashes and resolved libraries. Both solvers resolved Nektar++ 5.9.0
libraries from the MPI-enabled prefix, Open MPI and the system OpenBLAS
pthread LAPACK and BLAS libraries.

| Executable under the run root | Producing revision | SHA256 |
| --- | --- | --- |
| `build/IncNavierStokesSolver` | `f6e8ff83` | `6d9521f1923e7fc79a808192a89726fde6b2851ec48b4abe507c7278cf76d252` |
| `projection-build/IncNavierStokesSolver` | `d3bc3f5f` | `d69ff5ebf8202141c978494467940108f96c0858956007f3e8943b9c20f9f93d` |

The first-stage startup and incoming guess matched the archived h0 hash
`b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71`.
The native callback history matched h0 coefficient values. Its raw hash
`f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`
differs by signed zero bits from the stored guess. Both the zero-sweep
deliberate rejection and the failed physical attempts retained their
initialized coupled state, native time labels and all four history slots.
The accepted step count is zero and the physical clock remains at 0.

### Operator and regression checks

| Native check | Result |
| --- | --- |
| Same-state flow, k and epsilon block agreement | Maximum scaled difference `3.33e-13` across the reference, nonuniform perturbation and post-update checks |
| Frozen flow versus direct dense saddle solve | Maximum coefficient difference `3.61e-16`, including both regional pressure spaces |
| Native affine pressure action | Maximum momentum row difference `2.01e-18` |
| Predictor/corrector continuity and common gauge | From `(8.54e-4, 2.53e-4, 1.17e-5)` to `(6.87e-16, 8.68e-16, 6.22e-17)` |
| Velocity correction and frozen flow defect | Largest velocity coefficient change `5.60e-4` and maximum frozen row defect `7.24e-17` |
| Interface and relative pressure | Integrated reduced-pressure jump `-0.4557867404` in both paths and maximum interface observation difference below `3e-17` |
| Existing native regressions | 5/5 passed, including equilibrium evolution and the older backend smoke cases |
| Constitutive kernel | 28,434 checks passed on Linux, including active/capped rays, the zero-k rate, signed production and the epsilon reaction |

The frozen systems had 980 velocity unknowns, 85 pressure and gauge
unknowns, 504 free k modes, and 254 epsilon modes including the matching
coefficients. Interface observations compare two frozen-flow solutions
and include the physical traction and flux terms. The discontinuous
regional pressure result is retained in the native representation.

### Unchanged first-stage target

All comparisons used the original h0, alpha 400, time .0025, zero
prescribed forcing and tolerance `2e-12`. The reported full residual is
the maximum original scaled block norm from a fresh readback. The
projection comparisons used `--ground-projection`, with
`--k .75 --epsilon 1` for the second relaxation setting. The
pseudo-time comparison used `--pseudo-alpha 400`.

| Run | Completed sweeps | Full residual | Solver seconds | Stop |
| --- | ---: | ---: | ---: | --- |
| Default | 21 | 13.69367965 | 239.90 | Admissible residual-decreasing sweep unavailable |
| Pseudo alpha 400 | 23 | 13.94347378 | 188.91 | Admissible residual-decreasing sweep unavailable |
| Ground projection | 120 | 3.17771676 | 205.18 | Sweep budget |
| Ground projection, higher scalar relaxation | 74 | 3.17771675 | 145.79 | Admissible residual-decreasing sweep unavailable |

The default terminal relaxation for k was `2.98e-8`, and the
continuous ground minimum was `2.21e-10` near x=1.596. The finite
evaluation inventory still had minimum k `1.07e-4`. This difference
localizes the step restriction to the continuous ground condition.
Pseudo-time produced similar k relaxation collapse. Projection completed
all 117 triggered projections in its first run and continued through the
ground contact. Its k residual changed by about `1.01e-6` between
sweeps 90 and 120, with the ground minimum at the arithmetic margin.
The higher-relaxation run reached the same plateau with a final ground
minimum of `2.21e-16` and a trace coefficient defect of `7.22e-15`.
Its fresh readback gave the following full norm vector in native order.

| $u$ | $v$ | $k$ | $\varepsilon$ | Inner continuity | Outer continuity | Gauge | Constraints |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1.36e-5 | 6.10e-8 | 3.17771675 | 1.47e-9 | 7.40e-16 | 9.65e-16 | 1.05e-17 | 9.94e-14 |

Each fresh readback reproduced the complete terminal norm vector
exactly as serialized by the producing run.

The default terminal k residual decomposes into signed projections
along the original mass-dual residual direction. Inner destruction
contributes `+57.2282`, storage `-23.0076`, diffusion `-17.6386`,
inner production `-2.88685`, and advection `-0.001499`. Outer terms
nearly cancel, and prescribed scalar forcing is zero. The decomposition
uses an independently integrated cubic mass and agrees with the native
dual norm. The full term reassembly differs from the native k normed
rows by less than `8e-13` in the default case. This is a ground-limited
k balance failure at order-one residual. The independently rebuilt
trace constraint is much smaller and continues to converge while k
stalls.

The saved endpoints pass the scalar evaluation-domain check and have
minimum normalization lambda one. Normalized outer-gradient departures
are `.00288665` for the default run, `.00288513` for pseudo-time, and
`.00288904` for both projected runs, against the original `.001`
bound. They fail the original constitutive qualification. The retained R7
audit executable rebuilt the finite-set qualification from each saved
candidate in an audit-only mode on the original selected stress and source.

The focused multi-step and perturbed-guess trajectory checks are
conditional on a qualified first stage. The relaxation and pseudo-time
comparisons are first-stage diagnostics. The evidence supports a
ground-constrained k residual plateau for the tested split updates.
The next numerical experiment is to seek an admissible residual descent
direction that accounts for the feedback of the flow and epsilon
corrections on the active ground k modes, while retaining the block
solves and the original qualification. The current experiments leave
the existence of a qualified root unresolved. The earlier R7 coupled method
reached 2.9261 on an in-box nonroot from a different starting state.
Its result remains a separate historical experiment.

### Artifacts

The Linux root is
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-segregated-picard`.
The run directories contain `target/solver.log`, `target/sweeps.json`,
`target/case.evolution.json` and `target/readback.json`. Build receipts
are in `build-receipt` and `projection-build-receipt`. The five existing
test results are in `regression.log`. The two smoke tests also passed on
the second build and are recorded in `projection-regression.log`.
Candidate qualification and its
input packets are under `terminal-audit`.

The repository holds the [validation receipt](cls-epawf-segregated-picard-evidence/validation.json),
[compact Linux histories](cls-epawf-segregated-picard-evidence/linux-sweeps.json)
and [collection script](cls-epawf-segregated-picard-evidence/summarize-linux.py).
