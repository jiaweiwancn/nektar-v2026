# Neutral ABL closure and analytical wall-function kernels

Milestone 2 is **CLOSED**. Reusable kernels now evaluate standard k–epsilon,
the full high Reynolds number CLS cubic stress and transport sources, the
original smooth/rough analytical wall function (AWF), and the supplied note's
equilibrium dissipation correction. Linux tests verify their local mathematical
contracts independently. Coupling these kernels into an unsteady RANS solver
remains Milestones 3–5.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | implementation and numerical verification |
| Origin date | 2026-10-29 |
| Version label | neutral_abl_rans_m2_v1 |
| Milestone status | CLOSED on 2026-10-29 |
| Verification status | Closure/AWF kernels and native column regression VERIFIED on Linux |
| Tested source parent | `fd82c378cd4aba02f9d77946f79beb548ceece45` plus the exact source hashes in the validation manifest |
| Project branch | `codex/neutral-abl-rans-awf` |
| Linux CPU limit | 16 concurrent project workers, including compilation, MPI, and threads |

## Deliverables and API contracts

The implementation lives in
[`RANSModels.h`](../../library/SolverUtils/Turbulence/RANSModels.h),
[`RANSModels.cpp`](../../library/SolverUtils/Turbulence/RANSModels.cpp),
[`AnalyticalWallFunction.h`](../../library/SolverUtils/Turbulence/AnalyticalWallFunction.h),
and [`AnalyticalWallFunction.cpp`](../../library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp).
The main SolverUtils CMake lists register both implementation files and their
installed headers. The isolated validation build compiles these same sources
into a static library, without replacing libraries used by ongoing HIT runs.

All quantities are per unit density, using physical dimensional inputs. The
gradient convention is `gradient[i][j] = dU_i/dx_j` in a fixed, nonrotating
Cartesian frame. The three-dimensional tensors use `S = G + G^T` and
`Omega = G - G^T`, with no factor of one half. A two-dimensional caller must
embed its gradient in the corresponding three-dimensional tensor, retaining
the third normal stress.

| Interface | Contract |
| --- | --- |
| `EvaluateRANSClosure` | Returns Reynolds covariance, separate quadratic/cubic stresses, production `-R:G`, k and dissipation source components, invariants, divergence, eddy viscosity, and momentum/k/dissipation diffusion coefficients. |
| `RANSParameters` | Configurable standard-model coefficients, transport constants, and algebraic Yap constants; CLS retains its published variable Cmu relation and cubic coefficients. |
| `RANSTotalDissipation` | Explicit conversion `epsilon = epsilonTilde + nu*|grad(k)|^2/(2k)` when a caller chooses that convention. |
| `RANSNeutralReference` | Immutable reference constants solved from the selected coefficients; anchored log profiles, compatible derivatives/gradient flux, standard shifted-origin option, and standard log-wall traction. |
| `AnalyticalWallFunction` | Original AWF integrals, viscosity and dissipation distributions, tangential traction, velocity reconstruction/gradient, and volume-mean production/dissipation for current K. |
| `CLSEquilibriumWallFunction` | Current-K counterfactual reference, cutoff/admissibility diagnostics, fixed-weight or independently sensed correction, fallback/strict policies, and the common blended dissipation distribution. |

Inputs and coefficients must be finite; k and transported dissipation must be
positive. CLS requires positive **geometric** wall distance. AWF requires
positive physical molecular viscosity and `0 <= h < H`, `0 < P < H`; its
reconstruction accepts `0 <= z <= H`. Invalid arguments throw exceptions;
unrepresentable arithmetic results raise range errors. This makes invalid
states visible to the future time integrator. No solver stabilization policy,
k/epsilon floor, production clipping, or realizability projection is hidden
inside these kernels. A test explicitly retains negative CLS production in
a strong strain case. The original closure is not globally realizable.

The raw velocity gradient is used and its divergence is returned. Incompressible
states have `trace(R) = 2k`; a small numerical divergence is not silently removed
by a pointwise projection. Future assembly must monitor divergence and handle
positivity through its chosen transport algorithm.

### Resolved closure and reference constants

Standard k–epsilon uses the usual linear stress with constant Cmu. CLS uses
the complete quadratic and cubic relation, including the `S^2 Omega - Omega S^2`
term and the invariant-weighted cubic terms; the published c5 is exactly zero.
Production contracts the complete covariance with the supplied gradient. The
algebraic Yap source uses physical wall distance, independent of AWF roughness
or its dissipation cutoff.

CLS uses the note's high Reynolds number outer convention: damping is one,
Cepsilon2 is 1.92 by default, the low Reynolds number E term is zero, and the
transported modified dissipation equals the k sink in this limit. The optional
`resolvedNu` adds molecular diffusion only. It does not enable low Reynolds
number damping, the full original dissipation split, or a preservation guarantee
at finite molecular diffusion. AWF always retains its separate physical nu.

The default neutral constants are computed, rather than inserted as rounded
values:

| Model | q in `k = q u*^2` | kappa | Reference Cmu |
| --- | --- | --- | --- |
| Standard k–epsilon | 3.3333333333333335 | 0.43266615305567874 | 0.09 |
| CLS with algebraic Yap | 3.1906750816242875 | 0.4405539051866587 | 0.09822789532821127 |

CLS solves `Cmu(q) q^2 = 1` and
`q kappa^2/sigmaEpsilon + HY(kappa q^(3/2)/yapLength) = Cepsilon2-Cepsilon1`.
Changing transport/Yap constants recomputes kappa. The standard reference uses
its corresponding coefficient relation. Both references expose anchored
profiles for positive height and u*. Only the standard model accepts the
Richards–Norris shifted-origin option. Its log-wall traction accepts two
tangential components and reverses sign with velocity.

### Original AWF and equilibrium correction

The original wall object is constructed from **current** K, geometry, and wall
coefficients. Its viscosity slope uses `cMuWall*cLength*sqrt(K)`; its smooth or
rough intercept follows the AWF paper. The resolved CLS viscosity never replaces
that auxiliary viscosity. Reconstruct the wall object whenever K changes.

The integrals are
`J0 = integral(1/D)`, `J1 = integral(max(z-h,0)/D)`, and `J2 = integral(1/D^2)`,
where `D = nu + nu_t,A`. Traction is `U_H/J0 - (J1/J0) F`, in the note's positive
shear convention. The momentum force on the fluid has the opposite sign.
Reconstruction satisfies no slip at zero and the supplied matching velocity
at H. Tangential reversal and two-component forcing use the same operator.

Integrals use exact piecewise formulas with `log1p` and small-argument series.
Production is evaluated through positive moments of `nu_t,A/D^2` and a completed
square. This avoids subtracting `J0 - nu J2` near an entirely molecular layer
and avoids cancellation between large forcing terms. The only numerical clamp
removes roundoff-sized negative variance in a positive moment matrix; a larger
loss of positivity raises a range error. It does not clip modeled production.

For the correction, the reference uses current `u0 = sqrt(K/q)` and **zero**
forcing, regardless of the actual forcing. It computes `U_H,0 = (K/q) J0`,
`P0 = (K/q)^2 integral(nu_t,A/D^2)/H`, `A = u0^3/kappa`,
`lambda = H P0/A`, and `zc = H exp(1-lambda)`. The reference cutoff must be
representable and satisfy `zc <= P`. The cutoff is never clipped to P.

One cell-uniform effective weight blends the original and corrected
distributions. Its mean, P value, and H trace all come from that same
distribution; at weight one the mean is P0. All weights retain the original
traction and production **bit for bit** for identical inputs, and weight zero
recovers all original outputs exactly. The API reports requested/effective
weights and a reason for activation or fallback. An inadmissible cutoff or
caller-reported departure from the intended high Reynolds number wall regime
either selects the original wall law or throws under `RequirePreservation`.
The caller must supply the wall-regime decision; no unvalidated Reynolds-number
threshold is invented here. A baseline fallback has no preservation guarantee.

Automatic sensing uses independent matching velocity, the P-to-N logarithmic
velocity difference, kN, epsilonN, pressure forcing, convection forcing, and
normal velocity at P. Pressure and convection are sensed separately even when
their sum cancels. The incoming direction is prescribed, then normalized;
it is never inferred from local velocity. N must lie above H. The switch has
the note's exact plateaus and compact C2 transition. Thresholds and seven
nonnegative sensor weights are explicit experimental inputs. Any departure
that must trigger recovery independently needs a positive weight.

## Independent checks and completion evidence

[`RANSKernelProbe.cpp`](../../library/UnitTests/SolverUtils/RANSKernelProbe.cpp)
and [`validate_m2.py`](../../solvers/IncNavierStokesSolver/Tests/RANS/validate_m2.py)
provide complementary checks. The C++ executable performs **11,250 assertions**,
including 24 rejected inputs/policies. It evaluates 78 tensor cases and 60
wall cases. The Python oracle independently evaluates the dimensional tensor
relation in index form with NumPy and integrates the original wall differential
relations using SciPy quadrature. It does not call the C++ moment primitives.

The tensor cases cover zero gradient, simple shear, rigid rotation, pure strain,
mixed three-dimensional flow, 32 deterministic general gradients for each
closure, raw nonzero divergence, custom coefficients, symmetry, trace,
coordinate-frame covariance, production/source components, and diffusion.
Separate checks cover roots for default/custom coefficients and a deliberately
zero-Yap variant, the reference stress ratios, compatible epsilon gradient flux,
standard shifted heights/log traction, extreme invariants, and the explicit
dissipation conversion.

Wall geometries cover smooth, roughness below/above the viscous intercept,
fully rough, the `h*=70` zero-intercept transition, entirely molecular layers
on both sides of the dissipation cutoff, a weak turbulent segment, a very thin
forcing layer, custom constants, and changed current K. Each uses multiple
forcing values/components and reversed velocity. Tests check reconstruction,
no slip, pressure-response identity, positive production, dissipation integrals,
the molecular/turbulent power balance, weights 0/0.25/0.7/1, uniform-profile
consistency, cutoff equality and immediate failure below it, strict/fallback
behavior, independent departures, cancelling forcings, and switch endpoints.

| Measured error | Maximum |
| --- | --- |
| Full stress versus independent index form, scaled | 1.75e−16 |
| Quadratic/cubic stress versus independent index form, scaled | 1.90e−16 |
| Stress frame covariance, scaled | 4.71e−16 |
| AWF production versus independent SciPy quadrature, relative | 8.62e−16 |
| AWF production versus native 16-point Gauss quadrature, relative | 1.20e−12 |
| Reference power balance, relative | 1.43e−16 |
| Corrected mean versus independently integrated distribution, relative | 4.86e−16 |

The validator defines each scale and tolerance; the manifest records maximum
errors and hashes that validator. Zero-crossing
intercepts use the smooth intercept/height as their scale; very small eddy
viscosities use molecular viscosity. The independent quadrature uses local
segment coordinates with extended-precision interpolation so a thin `z-h`
is not lost by subtracting nearly equal double-precision heights. Initial
failures in these oracle comparisons were resolved at those numerical scales;
the final source passed all recorded checks.

For the note's `nu=1.5e-5`, `H=2`, `h=0.1`, `P=1`, `u*=0.4` example,
the kernels give `K=0.510508013059886`, `U_H,0=5.780254367996081`,
`P0=0.46109260899813753`, and `zc=0.009515269139914385`.
The original mean dissipation is `0.7749439664570105`; the corrected mean
at weight one equals P0. Original traction remains `0.16` in both cases.

The Milestone 1 native column now calls the **general production CLS kernel**
at its projected fields and native gradients. Its geometry, residual norms,
source-ownership tests, h/p/quadrature studies, incompatible-boundary controls,
and frozen-linear timestep retain their original acceptance criteria. The
degree-14 normalized residuals are momentum `1.0110e-10`, k `2.1353e-11`,
dissipation `3.3207e-9`, and normal momentum `2.9986e-10`. Serial/four-rank
results differ by at most `2.835e-17`. N is still owned by rank 1, the wall
strip by rank 0, with other ranks owning no strip.

Retained evidence:

- [Validation manifest](neutral-abl-rans-awf-milestone-2/validation.json): exact source/binary hashes, commands, constants, and comparisons.
- [Kernel report](neutral-abl-rans-awf-milestone-2/kernels.json): inputs and outputs for all tensor, wall, and sensing cases.
- [Serial column report](neutral-abl-rans-awf-milestone-2/column-serial.json) and [four-rank report](neutral-abl-rans-awf-milestone-2/column-mpi.json).
- [Environment report](neutral-abl-rans-awf-milestone-2/environment.json): compiler flags, installed native library hashes, allocation, and build/test log hashes.

GNU 11.4.0 Release compilation passed with warnings treated as errors, followed
by both registered CTests and the independent audit. The tested Linux worktree
is under `/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m2-kernels`.
The manifest identifies parent `fd82c37` plus exact source bytes, which are
checked against the milestone commit. The primary Linux checkout remained
clean at `812a2c84946d1fdad2b042c265e3e4ed32914138`. Its installed native
dependencies are unchanged from Milestone 1 and have no source changes between
that primary revision and the tested parent. The shared build/installation
was not rebuilt while HIT simulations were active.

Builds used four compiler jobs; tests ran sequentially with at most four MPI
ranks and one OpenMP/BLAS thread per rank, within CPU affinity 112–127. This
stays below the project's 16-worker limit. Reproduction commands are in the
[test README](../../solvers/IncNavierStokesSolver/Tests/RANS/README.md).

## Remaining implementation

These checks close the independent-kernel completion condition in the
[proposal](neutral-abl-rans-awf-proposal.md). No new session solver option is
exposed yet. Milestone 3 must add standard-model URANS transport, changing
diffusion and stiff sources, pressure consistency, controlled temporal
convergence, and the Richards–Norris control. Milestones 4–5 must integrate the
full CLS stress and both AWF policies with the conservative strip. The full
multidimensional inlet-preservation, perturbation, and restart study remains
Milestone 6. Local reference compatibility and convergent column residuals do
not establish exact finite-order analytical preservation or nonlinear stability.

The mathematical sources remain the supplied
[CLS/AWF note](../literature/cls_neutral_abl_equilibrium_awf_note.tex) and the
three primary papers linked in the [proposal's references](neutral-abl-rans-awf-proposal.md#source-references).
