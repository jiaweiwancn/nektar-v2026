# CLS–EP-AWF optional post-k bounding

This experiment continues `codex/cls-epawf-pimple-one-sweep` from the actual
local and fetched remote HEAD `c166f4441b5f89beed2874f267e2d5df1c75c335`.
The separate Picard checkout and its subsequent work are preserved.

`RANSWallPimpleKBound=Off|MassProjection` is an opt-in, algebraic post-k
operation. Off remains the default. The bound changes the numerical
step map; the repaired state is not an exact solution of the original
frozen k transport system. The U, pressure/gauge, outer epsilon and k
transport equations are each solved once. There is no inner-epsilon
transport equation, repeated governing solve, outer iteration, epsilon
trace restoration or field under-relaxation.

## Pinned OpenCFD reference

The official OpenCFD OpenFOAM-v2606 tag resolves to
`481094fdf34f11ed6d0d603ee59a858a0124236d`. All eight requested files were
read from Git objects at that exact revision; their SHA256 hashes are in
[reference provenance](cls-epawf-pimple-k-bound-evidence/openfoam-reference.json).
No OpenFOAM code or dependency is compiled into this implementation.

In [pimpleFoam.C](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/applications/solvers/incompressible/pimpleFoam/pimpleFoam.C),
turbulence correction follows the momentum predictor and pressure loop,
inside the outer loop and behind `turbCorr()`. The
[pimpleControl defaults](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/finiteVolume/cfdTools/general/solutionControl/pimpleControl/pimpleControl.C)
are one outer corrector, one pressure corrector and turbulence on the final
outer iteration. Its [inline control](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/finiteVolume/cfdTools/general/solutionControl/pimpleControl/pimpleControlI.H)
permits other configured schedules; OpenFOAM is not inherently restricted
to one sweep.

[kEpsilon.C](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/TurbulenceModels/turbulenceModels/RAS/kEpsilon/kEpsilon.C)
solves and bounds epsilon, then solves and bounds k, then updates turbulent
viscosity. Its frozen epsilon/k destruction coefficients enter implicit
`Sp` terms. That implicit reaction treatment and the subsequent field
bound are distinct operations. The equations also call relaxation and
finite-volume option hooks. Our temporal unit updates and physical source
laws remain the previous baseline.

[LienCubicKE.C](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/TurbulenceModels/incompressible/turbulentTransportModels/RAS/LienCubicKE/LienCubicKE.C)
also solves/bounds epsilon before solving/bounding k, including implicit
sink coefficients. It is a comparator for nonlinear turbulence-model
practice, not the same closure as CLS–EP-AWF.

[bound.C](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/finiteVolume/cfdTools/general/bound/bound.C)
activates when the field minimum is below its lower bound. For nonpositive
internal k, it takes a positive-neighbour-type target computed from the
field already clamped to the lower bound. It retains the maximum of raw
k and that target, then clamps every internal value to the lower bound.
A positive value below kMin therefore receives the lower clamp, whereas
nonpositive values can receive a substantially larger average. Boundary
values are also lower-bounded, followed by evaluation of coupled patches.

The [average implementation](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/finiteVolume/finiteVolume/fvc/fvcAverage.C)
linearly interpolates the clipped volume field to faces and then takes a
face-area-weighted cell average. It includes contributions involving the
cell itself; it is not simply an unweighted average of positive neighbours.
For fixed cell volumes this internal operation is upward-only, so its
change in integrated k is the nonnegative sum of cell volume times the
value increase. The printed global average is not that integral budget.

[RASModel.C](https://gitlab.com/openfoam/core/openfoam/-/blob/481094fdf34f11ed6d0d603ee59a858a0124236d/src/TurbulenceModels/turbulenceModels/RAS/RASModel/RASModel.C)
reads/adds kMin in the RAS dictionary, with velocity-squared dimensions
`[0 2 -2 0 0 0 0]` and a default SMALL value. Runtime dictionary rereading
can update it. This numerical default is not copied into our physical
model scale. A post-solve repair is not a positivity theorem for the
uncorrected transport discretization.

## Native map

The pilot minimizes `delta_f^T M_ff delta_f / 2` relative to the original
raw post-transport k candidate, subject to native physical evaluation
rows `k_raw(x)+E(x) delta_f >= floor + arithmetic reserve`.
`delta_c=0` at every actual essential coefficient. M_ff is the actual
physical consistent k mass restricted to the shared free CG coefficients.
A Cholesky whitening and dual active set with two-pass QR solve this small
convex problem. No k transport inverse enters the corrector. Algebraic
active-set, triangular mass-factor and KKT work is counted separately.
The projection can redistribute k downward elsewhere; unlike OpenFOAM's
upward FV repair, neither the sign nor conservation of its integral is
assumed.

The inner law admits nonnegative k, including its existing zero-k
extension. This pilot uses one small positive target at all protected k
locations, a stronger condition compatible with this case's essential
data. It does not redefine that physical zero-k extension.

The constraints include static k evaluation sets, strict sampled K and
outer-k sets, available pre-k normalization rays, and exact extrema of
each continuous cubic ground trace. Corrected states are scanned again;
new ground extrema and moving-ray k witnesses add linear cuts to the
same problem relative to the original raw state. The existing ScalarDomain
ordering checks all newly moving k values before normalization. An epsilon,
covariance or normalization failure is not converted into another field
repair. Finite volume samples do not certify positivity throughout every
interior volume polynomial. The claim is the declared physical inventory
plus the continuous ground trace.

The single positive floor is `RANSWallPimpleKFloorRelative * K0`, where
`K0=Q*u_star^2` has k dimensions (velocity squared). The default relative
floor is `1e-8`. A stated reserve `4096*double_epsilon*K0` is added to the
QP targets, and separation requires at least half that reserve above the
physical floor. Existing hard domain checks are unchanged. Essential data
have precedence; incompatible floor constraints fail. An admissible raw
state above the floor/reserve is returned bitwise unchanged. On startup
and restart the represented state must pass the corrector as a no-op.
The reference and actual perturbed h0 are checked separately.

Only k free coefficients are written. Solved epsilon coefficients are
retained; their target discrepancy is recomputed after correction. The
raw transport residual, added `A_k delta`, corrected frozen residual,
`alpha M_ff delta`, mass norm, signed integral budgets and active-set
certificates are all distinct diagnostics. The native time/history
transaction commits once after the final hard checks.

| Operation | Frozen/latest state |
|---|---|
| Momentum predictor | Incoming coupled coefficients and actual native callback history |
| Pressure/gauge | Same complete frozen G, including roughness forcing; `B=(alpha M_u,ff)^-1` |
| Outer epsilon | Corrected flow, incoming k/epsilon coefficients; matching target evaluated once before this solve |
| k | Corrected flow and solved epsilon; frozen transport/source coefficients, actual inner destruction rate and its zero-k extension |
| Optional bound | Original raw k candidate, physical mass metric; pre-k rays plus newly discovered k witnesses |
| Final audit/commit | Corrected k and retained solved epsilon; refresh dependent fields, report matching lag, commit once |

The baseline remains native IMEX1, a candidate first-order sequential,
linearly implicit method. Eta and additional Picard pseudo-time are zero;
all field relaxation factors are one. The selected CLS variant,
SaturatedR2 inner stress, capped destruction/smoothing, molecular
viscosity, normalization/transition laws, original outer6 mesh/spaces,
boundary data and zero prescribed forcing are unchanged. The bound does
not turn the method into an exact backward-Euler or converged M8 solve.

On checkpoint signatures explicitly bind mode, algorithm version,
consistent metric, constraint policy, effective default/explicit floor,
K0 and arithmetic reserve. The corrector is stateless per step. Off omits
this new signature suffix and ignores inactive bound parameters in the
signature, preserving old Off checkpoints. Correction integrals are
reported per attempt; experiment summaries accumulate accepted attempts
only and combine restart segments explicitly.

## Linux validation

The [compact native evidence](cls-epawf-pimple-k-bound-evidence/native.json)
contains every attempt's raw/corrected hashes, minima/locations, correction
budgets, frozen residuals, counters, qualification exits and artifact pins.
Heavy artifacts are under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-pimple-k-bound/` on
`jay-Super-Server`. `artifact-index.json` pins 554 files, including raw and
corrected coefficient snapshots, incoming histories, frozen matrices,
fields, checkpoints and logs. Its SHA256 is
`3d98a3ea46c7d62ad1bde5132d7da465aa2b9989118732db710534bf47b0c814`.

The authoritative solver was compiled from
`c7ce4cfb1a624771a64d9b1dc3e73845791de339`; its SHA256 is
`7c0d132c40a68c7c01fec5552dcb2029dfbe4016312d12ee189e3bbdd8758db9`.
The readback probe was compiled from
`33bba9680ff492afb68c65c42959ae17749e4dbf`; its SHA256 is
`e49bbdf42507b9b424081a43dced5de8d77b72c30c2486811282e2bc8727a902`.
The later commits add experiment/readback tests; final analysis ran from
`443e99e05dd77c121d159cdff85db7cbb897c3d7` with clean Linux source.
The evidence pins changed source files, linked libraries, XML files,
commands and log hashes. Runtime binary/library identities were checked
unchanged by the runner.

`linux-usb` was unavailable, so validation used the documented Tailscale
fallback. Source/build/run directories were isolated; the shared Linux
checkout, installation and concurrent Picard experiments were preserved.
The Release build used eight compiler workers and
`-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS`. Runtime used one
rank, all six recorded thread settings equal to one, and at most three
trajectory jobs concurrently. The batch audit found 36 existing OS
worker/helper threads; compilation and additional serial checks stayed
within the 72-worker limit.

Eight focused CTests passed in 40.72 seconds: QP checks, reduced-model
checks, native history rejection, native time integration, equilibrium
evolution, regional patch, evolution smoke and stage smoke. The QP
known-answer checks also passed ASan/UBSan. Existing Picard verification
passed: verified/unverified retained states were identical, exact saddle
comparison differed by only `3.61e-16`, and the unconverged stage remained
rejected under its original acceptance rule.
The existing fully refreshed `2e-12` nonlinear acceptance rule remains at
both enforcement locations for converged Picard stages. The PIMPLE runner
also checks frozen scalar systems, signed source/storage agreement and
retained epsilon trace mass columns on every scheduled completion.

Native corrector checks passed for exact no-ops on a constant, smooth
polynomial, the reference, and actual represented h0 separately. A cubic
ground fixture with positive native interpolation-node values had an
actual between-node minimum `-0.00170453648`; repair made its continuous
minimum positive. Essential coefficients and non-k coefficients stayed
bitwise unchanged, shared CG jumps remained below `1e-12`, and independent
quadrature reproduced the mass norm/integral. An incompatible prescribed
floor failed cleanly. Over the production trajectories, maximum primal
floor defect was zero, mass-dual stationarity below `1.9e-19`,
complementarity below `8.7e-23`, and active multipliers nonnegative.

For this case `K0=0.5105080130598861`, so the default physical floor is
`5.105080130598862e-9`, with arithmetic reserve `4.643043330905913e-13`.
These are k units, not dimensionless floors. No volume-wide polynomial
positivity certificate is claimed.

## Fixed-time evolution and qualification

Every case starts from the actual h0, SHA256
`b8ba7917424057a9f6f16a026d82d53ede43985871537facea170ef6e1841b71`.
The initial state is not replaced by the reference and no equilibrium
cancellation source is added. Off runs request T=0.0025 independently.
Bounded diagnostic runs request T=0.005 independently; T=0.0025 values
below are their actual committed 4/8/16-step prefixes, read back without
interpolation or another PDE solve. An independent coarse T=0.0025
checkpoint matches that prefix bit for bit. No failed horizon was halved.

| Mode | dt | Accepted steps | Final accepted time | First rejection attempted time | Raw ground minimum at rejection | Total runtime s |
|---|---:|---:|---:|---:|---:|---:|
| Off | 0.000625 | 3 | 0.001875 | 0.0025 | -0.00827757827 | 85.75 |
| Off | 0.0003125 | 5 | 0.0015625 | 0.001875 | -0.00530833752 | 127.50 |
| Off | 0.00015625 | 9 | 0.00140625 | 0.0015625 | -0.00730911267 | 212.50 |
| MassProjection | 0.000625 | 8 | 0.005 | none | — | 195.96 |
| MassProjection | 0.0003125 | 16 | 0.005 | none | — | 376.14 |
| MassProjection | 0.00015625 | 32 | 0.005000000000000002 | none | — | 703.36 |

The Off failures are near x=1.604–1.605 on geometry 4 at the ground.
The raw solves are accurate to roundoff. On and Off accepted prefixes,
and their raw candidates/incoming history through the first rejection,
are bitwise identical. Inactive On and Off runs to T=0.000625 also have
identical coupled state, time and all native history.

All long runs use **R2DiagnosticV1**. They retain hard scalar, covariance,
normalization and sampling checks, but are not physically qualified M8
evolution after the reference-box exit. First qualification exits are
0.00125, 0.0009375 and 0.000625 respectively, all initially from the outer
gradient box. A separate strict coarse run accepts one step to 0.000625,
then stops at its first qualification exit at 0.00125, before bounding
activates. No qualification bounds were widened.

| T | dt | Bound activations / steps | Continuous ground min | Sampled K min | B min | Refreshed epsilon trace defect | Cumulative integral delta k |
|---|---:|---:|---:|---:|---:|---:|---:|
| .0025 | .000625 | 1/4 | 5.10548e-9 | .5105022921 | .829717 | 3.23414e-6 | 2.86748e-6 |
| .0025 | .0003125 | 3/8 | 5.10552e-9 | .5105025793 | .819184 | 2.21511e-6 | 2.90838e-5 |
| .0025 | .00015625 | 7/16 | 5.10545e-9 | .5105035182 | .814465 | 2.79866e-6 | 7.88205e-5 |
| .005 | .000625 | 5/8 | 5.10537e-9 | .5104906512 | .742501 | 1.58771e-5 | 1.09701e-4 |
| .005 | .0003125 | 11/16 | 5.10546e-9 | .5104547174 | .718938 | 3.37479e-6 | 2.61375e-4 |
| .005 | .00015625 | 23/32 | 5.10550e-9 | .5104321973 | .704261 | 1.76533e-6 | 4.97802e-4 |

At these six states, both corrected regional continuity blocks and the
common gauge are below `1.5e-16`; the boundary defect is below `1e-13`.
Strong physical divergence is not zero: the final inner/outer L2 values
at T=.005 are `.01427/.0001103`, `.01611/.0001364`, and
`.01738/.0001530`. Outer k minima are .510484–.510507, outer epsilon
minimum is .03631791663, and lambda is one. The imposed frozen epsilon
trace defect is at most `1.03e-18` in this table. The refreshed defect
is the retained trace-coefficient infinity norm against the newly
projected target; it is not forced back to zero. It is nonmonotone at
T=.0025 under refinement.

The fully refreshed original residual maxima at T=.005 are 17.2092,
16.6706 and 18.5702. These are diagnostics of a single sweep with bounding,
not nonlinear acceptance tolerances. Per-attempt block norms are retained
in the Linux reports and readback summary.

## Refinement and correction size

These are physical quadrature L2 differences between adjacent dt runs at
the same time, with one common pressure-gauge alignment over both regions.
No separate regional mean is removed. Gauge shifts are below `8e-18`.

| T | dt pair | Velocity | k | Outer epsilon | Pressure | Wall stress | Wall B | Sampled wall K |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| .0025 | .000625 / .0003125 | .00109044 | .00471120 | 4.68929e-5 | .0686865 | .0506418 | .178418 | 1.68829e-5 |
| .0025 | .0003125 / .00015625 | .000789223 | .00314989 | 4.21316e-5 | .0340411 | .0235918 | .0667662 | 2.22879e-5 |
| .005 | .000625 / .0003125 | .00205388 | .00722887 | 2.85556e-5 | .0479688 | .0650910 | .103134 | 7.09805e-5 |
| .005 | .0003125 / .00015625 | .00129308 | .00391180 | 7.52505e-6 | .0230634 | .0163626 | .0557083 | 4.43308e-5 |

Wall lambda differences are zero. At T=.005 the observed successive
difference orders for velocity/k/epsilon/pressure/wall stress/B/K are
`.67/.89/1.92/1.06/1.99/.89/.68`. These three resolutions show decreasing
differences for the bounded map at this horizon. They do not establish
uniform first-order accuracy: T=.0025 has weak epsilon reduction and an
increasing wall-K difference. No exact or fully converged M8 root is used
as a reference.

The first coarse repair has mass norm `1.4259115883e-4`, 1.9536% of the raw
k-step norm, and integral `2.8674772076e-6`. The raw continuous minimum
`-0.00827757827` becomes `5.10547802e-9`. Raw/repaired state hashes begin
`4085794f237a7115` / `bd390aa1df48b54a`. Its raw frozen linear residual is
`2.22e-16`; the corrected frozen residual is `0.00600624410`, with added
`A_k delta` maximum `0.00600624410` and storage contribution maximum
`0.00156937910`. These are deliberately different residuals.

Independent saved-matrix readback checks reconstruct raw and corrected
frozen residuals for every saved attempt without another solve. Their
raw-vector discrepancy is at most `2.22e-16`. Across the three long runs,
the maximum corrected frozen residual is .1164/.1209/.1232, while raw
linear residuals stay below `5.56e-16`. The maximum correction/raw-step
mass-norm ratios rise to 13.4%, 23.3% and **41.8%**. Thus correction does
not become a rare, vanishing safeguard over this refinement range.

| dt at T=.005 | Net added integral | Positive contribution | Negative contribution | Near-ground signed integral (z<=.01) |
|---:|---:|---:|---:|---:|
| .000625 | 1.09700652e-4 | 2.36707195e-4 | -1.27006543e-4 | 1.34481562e-4 |
| .0003125 | 2.61375060e-4 | 5.53769664e-4 | -2.92394604e-4 | 3.20418572e-4 |
| .00015625 | 4.97801618e-4 | 1.04785498e-3 | -5.50053358e-4 | 6.10252882e-4 |

These are physical 2D integrals (per unit span). Almost all net addition
is in the inner region; outer totals are `1.13e-9`, `2.70e-9`, `5.13e-9`.
The negative contributions demonstrate redistribution by the consistent
mass metric. Cumulative correction budgets have not converged or vanished
under this refinement, even though the final field differences decrease.

Each accepted step has one outer sweep, one vector momentum solve, one
pressure/gauge solve, one epsilon solve, one k solve, one implicit callback,
and zero additional governing-equation solves. Pressure correction uses
85 mass-inverse applications per step. At T=.005, bound calls/activations
are 8/5, 16/11 and 32/23; active-set iterations are 161/367/742,
KKT solves 160/367/742, and triangular mass-factor applications
192/452/936. Each activated corrector factors M_ff once. Two mass-matrix
applications per bound call measure norms; a separate application measures
the storage residual. Initialization reports two Lift calls and three
iterations, separately from the physical schedule. There is no physical
step Lift or nonlinear epsilon trace restoration.

## Mechanism, floor sensitivity, and transactions

At the original coarse failure, the last accepted ground minimum is
0.05989979 at x=1.6095023 on geometry 4. Consistent-mass instantaneous
k rates at that point are:

| State | Advection | Diffusion | Inner production | Net inner destruction | Total rate |
|---|---:|---:|---:|---:|---:|
| Last accepted | -.038796 | 227.450438 | 18.490400 | -658.076083 | -412.174040 |
| After flow | -.052632 | 227.450438 | 18.931691 | -658.062434 | -411.732936 |
| After epsilon, before k | -.052632 | 227.450438 | 18.931691 | -658.062434 | -411.732936 |

At the fixed x=1.60 witness, corresponding totals are -412.037199,
-411.589665 and -411.589665. Outer production/destruction contributions
are approximately +/-8.7e-7, other k terms zero. Net destruction is counted
once. Component sums reproduce the independent native rate within
`1.8e-12`; mass-solve residuals are below `2.3e-16`. Each selected-state
diagnostic uses 19 mass-rate inverse applications, separately from all
production counts. This identifies strong represented near-ground
depletion, not an inaccurate raw transport solve. It is an instantaneous
spatial budget, not a clean full-sweep h-to-zero consistency experiment.
It does not by itself prove which spatial/implicit coupling destroys
positivity, or that the continuous physical law should cross zero.
The epsilon trace defects at these three phases are `3.14e-6`, `.001708`,
and `2.45e-9`, so their different trace/storage contexts remain explicit.

Saved raw failing-candidate replay gives:

| Relative floor | Physical floor | Correction mass norm | Integral delta k |
|---:|---:|---:|---:|
| 1e-10 | 5.10508e-11 | .0001425910716 | 2.867475448e-6 |
| 1e-8 | 5.10508e-9 | .0001425911588 | 2.867477208e-6 |
| 1e-6 | 5.10508e-7 | .0001425998797 | 2.867653122e-6 |

All three repair the continuous ground and complete the declared inventory.
Independent active trajectories at dt=.000625 reach T=.005 for all three
floors, with five activations each. Their cumulative integrals are
`1.097006120e-4`, `1.097006521e-4`, `1.097046574e-4`. Compared with the
default, the high floor changes velocity/k/epsilon/pressure physical norms
by `5.92e-9 / 3.39e-8 / 2.55e-10 / 2.75e-7`; wall-stress/B/K changes are
`5.14e-7 / 2.72e-6 / 7.59e-10`. Low-floor differences are about 100 times
smaller. These are much smaller than the time-step differences.

The near-floor field is stiff: maximum evaluated actual inner destruction
rates at T=.005 are about 17623, 18200, 18579 in the three dt runs; the
coarse rate changes only from 17623.147 to 17623.112 across the floor
range. Lambda remains one and no additional hard model-domain exit occurs
through T=.005. A positive-neighbour-target variant was therefore not
needed to resolve an observed floor/stiffness failure in this experiment.

Checkpoint/restart across the first active repair (T=.0025 to .005) is
bitwise identical to uninterrupted evolution in coupled state, time and
every native history vector/time. The accepted correction integral also
agrees when restart segments are combined. A deliberate hard-audit test
rejects step four after actual bound activation: the last accepted state,
time=.001875 and history match Off's accepted prefix bit for bit. The
rejection record contains the actual repaired candidate, complete incoming
history, measured audit and explicit test reason. It does not manufacture
a covariance failure. A changed mode or floor is rejected by checkpoint
provenance; explicit default 1e-8 is accepted. A historical Off checkpoint
from the earlier branch also restarts with bitwise-identical state/history.

An early development diagnostic omitted the independent `spatial_unmasked`
load and threw `map::at`; this adapter error was fixed in 62c21e1. Its
`initial/` runs are retained and are not labelled physical k failures.
An initial readback budget omitted the final finest step because its
binary64 time exceeded .005 by 2e-18; f4e52bc corrected prefix selection.
`readback-verified/` includes that final step and independently verifies
the saved matrices. These fixes did not modify the production step map.

## Reproduction

Use an isolated Linux source at this branch, build and run root as below;
do not overwrite the shared checkout/install. Existing commands and runtime
XML are also pinned per case in the evidence and Linux manifests.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
run=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-pimple-k-bound
src="$run/source"
tests="$src/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export TMPDIR="$run/tmp"
cmake -S "$tests" -B "$run/build" -DRANS_BUILD_SOLVER=ON \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH=/home/jay/.local/opt/nektar-v2026-mpi \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake \
  -DCMAKE_CXX_FLAGS='-fno-fast-math -ffp-contract=off -D_GLIBCXX_ASSERTIONS'
cmake --build "$run/build" --parallel 8 --target IncNavierStokesSolver \
  RANSWallResolutionProbe RANSKBoundQPChecks RANSHistoryRejectionChecks \
  RANSEquilibriumTimeProbe RANSEquilibriumEvolutionProbe \
  RANSEquilibriumPatchProbe RANSReducedModelChecks
ctest --test-dir "$run/build" --parallel 1 --output-on-failure \
  -R '^(RANSKBoundQPChecks|RANSReducedModelChecks|RANSHistoryRejectionChecks|RANSEquilibriumTimeNative|RANSEquilibriumEvolutionNative|RANSEquilibriumPatchNative|RANSEquilibriumWallEvolutionSmoke|RANSEquilibriumWallStageSmoke)$'
solver="$run/build/IncNavierStokesSolver"
probe="$run/build/RANSWallResolutionProbe"
python3 -B "$tests/check_cls_epawf_k_bound.py" --solver "$solver" \
  --probe "$probe" --output "$run/reproduce-corrector"
python3 -B "$tests/run_cls_epawf_pimple.py" --solver "$solver" \
  --output "$run/reproduce-checkpoint" --dt .000625 --final-time .0025 \
  --diagnostic --k-bound MassProjection --checkpoint
python3 -B "$tests/experiment_cls_epawf_k_bound.py" --solver "$solver" \
  --output "$run/reproduce-fixed" --jobs 3 \
  --checkpoint "$run/reproduce-checkpoint/checkpoint"
python3 -B "$tests/analyse_cls_epawf_k_bound.py" \
  --experiments "$run/reproduce-fixed" --probe "$probe" \
  --output "$run/reproduce-readback"
python3 -B "$tests/check_cls_epawf_k_bound_provenance.py" --solver "$solver" \
  --checkpoint "$run/reproduce-checkpoint/checkpoint" \
  --checkpoint-report "$run/reproduce-checkpoint/case.evolution.json" \
  --output "$run/reproduce-provenance"
python3 -B "$tests/test_cls_epawf_picard_verification.py" \
  "$solver" "$run/reproduce-old-picard"
```

Output directories must be unused. Inspect current server load before
building/running. For a single explicit horizon, use the dedicated runner's
`--final-time`, `--k-bound`, `--floor-relative`, `--snapshots` and optional
`--budget-time`; rejection remains a recorded result. `--hard-reject-step 4`
exercises actual-candidate hard-audit recording after repair.

## Interpretation

| Outcome | Result |
|---|---|
| Implementation/build | Passed; Off remains default and compatible |
| Frozen-operator checks | Passed; repaired k residual is explicitly nonzero |
| Single-sweep bounded trajectory | Completed at three dt values through T=.005 |
| Time-step convergence/accuracy | Decreasing differences at T=.005; uneven at T=.0025; uniform first-order accuracy not established |
| Original constitutive qualification | Exits early; strict run stops, diagnostic continuations are unqualified |
| Fully converged original M8 stage | Not solved or claimed; not this task's objective |

Bounding permits useful continuation and materially modifies the
trajectory. It is floor-insensitive over the tested small range, but its
increasing activation and cumulative supply prevent treating survival as
accuracy for the uncorrected equations. Keep this option experimental.
The next numerical step should target the near-ground semi-discrete
destruction/diffusion balance and the nonvanishing correction budget,
alongside a qualified temporal fixture; further refinement must assess
that budget as well as field differences. Increasing the floor or widening
qualification bounds is not supported by these results.
