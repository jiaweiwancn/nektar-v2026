# Neutral ABL wall discretization and residual prototype

Select a conservative control-volume wall strip coupled to the CG outer domain.
The strip retains one representative tangential velocity and k state per wall
patch and replaces its dissipation transport equation with the note's algebraic
condition. This avoids introducing unproved higher-order equations inside the
auxiliary wall layer. The outer spectral elements continue to solve the
specified high Reynolds number equations.

## Material Passport

| Field | Value |
| --- | --- |
| Origin skill | academic-research-suite, experiment-agent |
| Origin mode | implementation and numerical verification |
| Origin date | 2026-10-29 |
| Version label | neutral_abl_rans_m1_v1 |
| Milestone status | CLOSED on 2026-10-29 |
| Verification status | Strip/column prototype VERIFIED on Linux; complete URANS integration remains later work |
| Source parent | `d8411938909169c0b30651fb6727a43c47cf0a19` |
| Project branch | `codex/neutral-abl-rans-awf` |
| Project Linux CPU limit | 16 concurrent workers, including compilation, MPI, and threads |

## Decision and scope

Milestone 1 in the [implementation proposal](neutral-abl-rans-awf-proposal.md)
allows either a verified constrained spectral wall element or the conservative
strip alternative with its additional scope recorded. Select the strip for the
initial implementation. A native prototype demonstrates its column assembly
and a limitation of a direct pointwise CG interpretation of the wall sources.

For the note's fixed example, integrate the auxiliary production and corrected
dissipation over the wall layer. Their difference has zero mean, but its moment
against `2z/H - 1` is nonzero. Thus the integrated k budget does not close all
the equations introduced by a pointwise higher-order wall reconstruction.
This rejects that direct mapping; it is not an impossibility proof for every
possible constrained spectral element. A uniform mean source would eliminate
this particular source moment, but would still need separate justification of
its inner velocity/dissipation representation and equation constraints.

This milestone does not expose a new IncNavierStokesSolver solver type. The
general CLS cubic tensor, general AWF kernels and switch, nonlinear coupled
URANS integration, terrain, and downstream profile tests retain their planned
later milestones. The native residual prototype uses the simple-shear CLS
reduction, without molecular terms in the outer PDE, and fixed constants from
the note. It evaluates constitutive quantities from projected fields and
native derivatives; it does not insert analytical velocity gradients or a
reference residual cancellation.

## Geometry and unknowns

Use physical geometric wall distance z, with the stationary wall at z=0.
For each flat wall patch, `0 < zP < H < zN < Z`. The version 1 column uses
zP=1 m, H=2 m, zN=4 m and Z=30 m. The strip occupies `0 < z < H` and the
CG mesh occupies `H <= z <= Z`. Their volumes do not overlap. In multiple
dimensions multiply strip balances by the physical patch area.

P is the representative strip state, not a CG node or a quadrature point.
N is an explicitly located outer sample. The initial prototype places the N
plane on a CG vertex and assigns ownership to the element to its right.
All three locations remain fixed under p refinement and under subdivision of
the outer elements. General sample location and curved layers need separate
validation before widening this initial geometry contract.

Keep these unknowns distinct:

- Tangential velocity `UP` and turbulence energy `KP` in the strip.
- Dissipation `EP` at P, constrained algebraically, with no independent strip
  dissipation PDE or dissipation time history.
- Outer CG velocity, k, dissipation and reduced pressure.
- Interface momentum and k diffusion fluxes, eliminated through the strip
  balances or retained as equivalent coupling unknowns.

The sampled outer logarithm at P is a reference for the strip's representative
state. It is not the integral mean of that logarithm down to the wall.
The inventory `H UP` is the selected finite-volume approximation, as is `H KP`.

## Matching and the inlet anchor

Let uH and uN be outer values and choose

\[
w_U=\frac{\log(H/z_P)}{\log(z_N/z_P)},\qquad
u_H=(1-w_U)U_P+w_Uu_N.
\]

This constant-preserving interpolation reproduces a logarithm exactly, so its
input is the actual CG trace at H. With the note's zero-forcing family,
set `uH = t0 J0`, and use

\[
U_0(z)=t_0J_0+\frac{u_*}{\kappa_c}\log(z/H).
\]

The selected physical-height weights for k are

\[
w_k=\frac{H-z_P}{z_N-z_P},\qquad
k_H=(1-w_k)K_P+w_kk_N.
\]

Eliminate the representative states with these exact linear constraints:

\[
U_P=\frac{u_H-w_Uu_N}{1-w_U},\qquad
K_P=\frac{k_H-w_kk_N}{1-w_k}.
\]

Their time derivatives must be eliminated by the same operators. Dropping the
resulting strip inertia would change the transient equations even if a steady
reference still passed. The initial P/H/N geometry gives `wU=1/2` and `wk=1/3`.
No additive offset is selected using a different interpolation.

## Conservative strip rows

Use gradient fluxes positive in the increasing-z direction. Let tH be the
outer tangential momentum diffusion flux at H, tA the original AWF wall
traction in the note's positive-shear convention, qH the outer k gradient
flux, and qw the modeled wall k flux. The per-unit-area, per-unit-density rows
are

\[
R_U=H(\dot U_P+F)-t_H+t_A,
\]
\[
R_k=H(\dot K_P+A_k)-q_H+q_w
-H(\overline P_A-\overline\epsilon_{EP}).
\]

Here F is the mean resolved tangential convection/pressure/body-force
contribution in the strip momentum row, with the solver's consistent signs,
and Ak is mean k advection. The AWF auxiliary momentum forcing has its
published cutoff below roughness height h. Its reconstructed upper derivative
is not prescribed as an additional independent resolved interface flux away
from equilibrium. The resolved strip balance determines tH.

The CG domain's bottom outward flux is `-tH` or `-qH`. After moving diffusion
to the residual left side, the CG constant-test row contains `+tH` or `+qH`.
Adding that row to the strip row cancels the interface flux exactly. Only the
physical wall/top fluxes, inventories, transport, and volume sources remain.
The strip mean production and dissipation are supplied once. There is no
ordinary CG production/dissipation integration over `0 < z < H`.

On the reference, `F=Ak=0`, `tH=tA=t0`, `qH=qw=0`, and
`PbarA=EbarEP=P0`, so both rows vanish. Original AWF supplies its original
dissipation mean instead, giving a measurable nonzero reference k row.

The reusable header
[`RANSWallStrip.h`](../../library/SolverUtils/Turbulence/RANSWallStrip.h)
implements these scalar rows, matching constraints, elimination and integrated
normal mass balance. Two tangential components use the same scalar geometric
operators componentwise; the multidimensional assembly remains later work.

## Dissipation constraint and inner representation

Use the same auxiliary dissipation distribution to produce all three values:
its strip mean, its value at P, and its value at H. In the corrected plateau,

\[
\epsilon_{c,w}(z)=A_c/\max(z,z_c),\quad
\overline\epsilon=P_0,\quad E_P=A_c/z_P,\quad E_H=A_c/H.
\]

The row `EP - prescribedEP = 0` replaces the strip dissipation PDE. The outer
boundary row at H is `epsilonH - prescribedEH = 0`, replacing that boundary's
transport row. Its diffusive reaction is an output, not a second imposed
epsilon flux. Deriving EH from the common auxiliary distribution is an explicit
choice for the CV-to-CG interface, additional to the note's original single
cell-centred row. Preserve this choice in the original-AWF baseline as well
as in the corrected variant.

Both values can be eliminated before the coupled solve or imposed with
equivalent algebraic constraints. A post-solve overwrite of one quadrature
sample is not equivalent. The switch must use independent outer dissipation
information, such as at N, so prescribing EP does not define its own switch.

There is no polynomial resolved field on `0 < z < H`. The finite auxiliary
AWF reconstruction enforces tangential no slip at z=0 and supplies traction;
its capped dissipation supplies the mean and point values. The resolved and
auxiliary velocities at P are intentionally different. Neither log(z), 1/z,
nor the resolved geometric-distance Yap formula is evaluated at z=0.

## Normal velocity and pressure

Physical-wall impermeability is `Wwall=0`. Integrated strip continuity is

\[
W_H-W_{wall}+H\nabla_t\cdot\boldsymbol U_P=0.
\]

Consequently WH=0 on the homogeneous column. For a later multidimensional
unsteady flow, WH must follow the integrated strip continuity equation; imposing
impermeability independently at H would discard the strip's tangential mass
transport. Pressure is the outer incompressibility multiplier. Its normal
boundary treatment must use the normal momentum equation and this prescribed
interface normal velocity, including its time derivative. Do not prescribe
an independent pressure value there as well.

Store `pi = p/rho + 2k/3` after removing hydrostatic pressure. In the column the
remaining CLS nonlinear normal stress is

\[
R^{NL}_{zz}=-\frac{11}{75}k\frac{\nu_tG^2}{\epsilon},
\qquad \partial_z(\pi+R^{NL}_{zz})=0.
\]

It is constant on the reference and pi can be constant. The native prototype
assembles this weak normal-momentum balance, weak continuity, and departures
with nonzero pressure gradient or divergence. It does not run Nektar++'s
velocity-correction pressure projection. Coupling strip continuity to that
projection, tangential transport across neighboring patches, and the normal
pressure reaction are required integration work in milestones 3–5.

## Prototype and residual definitions

[`RANSWallLayerProbe.cpp`](../../library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp)
uses native `LocalRegions::SegExp` geometry, endpoint-constrained forward
transforms, backward transforms, physical derivatives, and weak inner products.
It assembles CG rows using shared endpoint maps and interior modes. The outer
domain starts at H; four base intervals are `[2,4]`, `[4,8]`, `[8,16]`, and
`[16,30]` metres. Tests use degrees 2 through 14, repeated subdivision at degree
4, and additional quadrature at degree 10.

Evaluate C(q), viscosity, production and Yap from the projected U and epsilon
fields and the native derivative of U. Use constant reference k. Momentum has
the original-AWF traction at H and reference velocity at Z. Dissipation has the
common-profile value at H and reference gradient flux at Z. k has zero modeled
wall/top gradient flux. The normal stress uses the same simple-shear reduction.

Residuals are mass-dual L2 norms `sqrt(R^T M^-1 R)` on the unconstrained test
space, divided by the relevant scale times `sqrt(Z-H)`. Scales are
`u*^2/H` for tangential/normal momentum, `u*^3/H` for k, `u*^4/H^2` for
dissipation, and `u*/H` for continuity. The replaced top velocity and bottom
dissipation rows are excluded from their differential norms; algebraic matching
is checked separately. This avoids comparing basis-dependent raw row maxima.

The stationary and timestep experiment freezes outer viscosity to
`kappa_c u* z` and K/J0 to the example inputs. Let eH select the H vertex and
EN evaluate the physical N sample. After eliminating UP and tH, the scalar
momentum matrices are

\[
M_c=M+\frac{H}{1-w_U}e_H(e_H-w_UE_N)^T,\qquad
L_c=L+J_0^{-1}e_He_H^T.
\]

The new mass contribution is generally nonsymmetric. Solve its frozen linear
steady problem and one backward-Euler step with dt=0.125 s. Recover strip
inventory changes and the top Dirichlet reaction to verify the total momentum
budget. This demonstrates exact algebraic elimination and timestep consistency
for that linear experiment; it does not establish nonlinear URANS stability
or the order of the future IMEX solver.

## Completion evidence

The final Release build passed with prototype warnings treated as errors.
The registered `RANSWallLayerPrototype` CTest passed, followed by the serial
and four-rank audit. The [validation manifest](neutral-abl-rans-awf-milestone-1/validation.json)
records the tested parent plus exact source hashes and commands. The
[environment report](neutral-abl-rans-awf-milestone-1/environment.json) records
compiler/build settings and installed library hashes. Full
[serial](neutral-abl-rans-awf-milestone-1/serial.json) and
[MPI](neutral-abl-rans-awf-milestone-1/mpi.json) numerical summaries are retained.

The tested source is parent `d8411938909169c0b30651fb6727a43c47cf0a19` plus
the six source files identified by SHA256 in the manifest. Those source bytes
are checked against the milestone commit before delivery. Native dependency
sources in LibUtilities, StdRegions, SpatialDomains, and LocalRegions are
unchanged between the primary Linux revision and this parent. The primary
Linux checkout remained clean at
`812a2c84946d1fdad2b042c265e3e4ed32914138`; no shared build or installation
was modified.

Measured wall-source values for the note's example are:

| Quantity | Measured value |
| --- | --- |
| Corrected source mean | −2.6421e−16 m2/s3 |
| Corrected source moment against `2z/H−1` | 0.008677571954 m2/s3 |
| Original AWF strip k residual at the same reference | 0.6277027149 m3/s3 |
| Difference between auxiliary and representative UP | About 0.04439 m/s |

The outer high Reynolds number residuals converge without a balancing volume
source. Representative normalized mass-dual results are:

| Polynomial degree | Momentum | k | Dissipation | Normal momentum |
| --- | --- | --- | --- | --- |
| 2 | 2.0508e−3 | 2.5740e−2 | 1.9455e−2 | 5.6936e−3 |
| 4 | 1.5116e−4 | 7.5980e−4 | 2.3370e−3 | 4.4152e−4 |
| 8 | 9.3345e−7 | 7.2450e−7 | 1.9062e−5 | 2.7787e−6 |
| 12 | 2.3305e−9 | 6.9248e−10 | 6.6379e−8 | 6.9186e−9 |
| 14 | 1.0102e−10 | 2.1349e−11 | 3.3212e−9 | 2.9949e−10 |

Degree-4 subdivision by eight reduced the respective residuals to 4.4063e−7,
3.6076e−7, 4.8834e−6, and 1.3044e−6. Increasing degree-10 quadrature from 28
to 40 points changes each residual by less than 1e−8. Continuity and matching
constraints are zero to the reported precision.

Across all cases, the largest normalized frozen-linear backward-Euler momentum
inventory error is 1.40e−12 and the largest change of its numerical stationary
state is 1.83e−13. At degree 14, one step from the projected analytical velocity
changes its coefficient vector by 1.00e−11 relative to u*. This timestep check
uses frozen viscosity and does not stand in for a nonlinear URANS step.

At degree 14, omitting Yap produces a normalized dissipation residual of
9.0020e−4; setting top epsilon flux to zero gives 7.7641e−2; and a one-percent
wall traction error gives a momentum residual of 4.0030e−2. The controlled
pressure gradient and divergent crossflow produce normalized residuals 0.01
and 0.005. These controls distinguish a sensitive residual diagnostic from a
hard-coded zero result.

Serial and four-rank results differ by at most 4.45e−17 across the compared
numerical fields. In MPI the N sample is on rank 1, the strip on rank 0, and
other ranks have no strip. Tests ran sequentially with one OpenMP/BLAS thread
per process, at most four MPI ranks, and a 16-CPU inherited affinity. The build
allowed at most four compiler workers. A warning in OpenMPI's legacy C++
headers encountered during the strict-build check was resolved by marking
external dependency include directories as system includes; warnings from
the prototype remain errors.

These results close the milestone's permitted conservative-strip alternative:
geometry, matching, volume/source ownership, dissipation row replacement,
normal mass/pressure interpretation, and algebraic elimination are fixed and
the selected column assembly contracts have passed native tests.

The acceptance checks require convergence of all four outer dynamic residuals
under p and h refinement, independent-quadrature agreement, balanced strip
rows/constraints, conservative frozen-linear stepping, and preservation of
that numerical linear equilibrium. Controls omit Yap, change the top epsilon
flux, perturb traction, apply a pressure gradient, and introduce divergent
crossflow. They must create detectable residuals.

MPI partitions native elements round-robin. All ranks participate in assembly
and sample communication, including ranks with no strip. The N-plane owner
differs from the strip owner in the four-rank case. Compare the serial and
four-rank reports with one OpenMP/BLAS thread per rank.

## Reproduction and remaining work

The [test README](../../solvers/IncNavierStokesSolver/Tests/RANS/README.md)
contains build/run commands. The external CMake project links the existing
Release MPI installation and leaves its shared build and installed libraries
unchanged. Validation uses a detached Linux worktree of the same repository
because two existing HIT runs are using the primary checkout and installation.
Record both revisions and preserve that primary checkout.

The project is limited to 16 concurrent CPU workers. Build with at most four
compiler workers for this prototype, run serial and four-rank tests sequentially,
and limit their Linux affinity to 16 CPUs. Store builds, TMPDIR, logs and full
results under `/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf`.

Later milestones must implement the full tensor and wall kernels, switch and
admissibility policies, nonlinear boundary derivatives, patch advection and
mass/pressure coupling, multidimensional geometry, and full URANS timestep
and restart. Exact finite-order retention of analytical samples remains an
open numerical question. This column prototype demonstrates convergent
residuals and the selected wall assembly, without claiming that later goal.
