# Segregated Picard for the CLS/EP-AWF physical stage

The backend and its native test harness are implemented. Linux validation is
blocked by access to `jay-Super-Server`. The original M8 convergence problem
remains open. No Picard physical stage has been run or accepted.

The branch is `codex/cls-epawf-segregated-picard`. Its base
`0ed7790450bb7e81f17ff1d0b577f12fc95f9cce` matches the earlier reviewed tip.
The implementation and harness are in `0486db8c39ce325bc0c2d004382bfc8881f6ea07`.
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
src="$run/source"
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

After a qualified first stage, run with `--steps 3` in a fresh output
directory for the focused multi-step fixture through .0075. Use
`--perturbation 1e-7`, `--perturbation -1e-7` and
`--perturbation 1e-6` in separate first-stage runs. A second relaxation
setting can use `--k .25 --epsilon .25` while keeping the same physical
target. Choose later continuation budgets and pseudo-time settings from
measured sweep, block and readback evidence. The serial Picard path has
no MPI comparison at present.

## Evidence as of November 4, 2026

| Check | Result |
| --- | --- |
| Mac syntax check of changed solver translation units | Passed with Clang C++17 and warnings as errors |
| Mac standalone robust constitutive probe | 28,434 checks passed, including active/capped rays, the zero-k rate, signed production and the epsilon reaction |
| Archived input preparation | Startup, deliberate rejection and first-stage fixtures generated |
| Linux source sync and build | Blocked. Direct and Tailscale SSH timed out, and Tailscale reported the target offline |
| Native consistency, pressure, rejection and regression tests | Pending Linux |
| Original M8 residual, constraints, qualification and readback | No Picard measurement |
| Accepted steps, time and perturbation runs | No native run |

See the [validation receipt](cls-epawf-segregated-picard-evidence/validation.json)
for local fixture identities and access results. Once Linux is reachable,
run the unchanged first stage with verification enabled. If it stalls,
the saved readback provides the k-term breakdown and fresh trace defect
needed to choose the next numerical change.
