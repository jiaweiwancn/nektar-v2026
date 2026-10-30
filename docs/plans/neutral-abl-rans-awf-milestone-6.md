# Neutral-ABL preservation and robustness study

The compatible standard k–epsilon model and the fixed and switched corrected
CLS/AWF policies retain their analytical neutral-ABL profiles within the
specified 0.1% engineering target for five flow-throughs in 2D and full 3D.
The original AWF leaves that reference and reaches the nonlinear iteration
limit. These results establish finite-resolution accuracy and perturbation
recovery for the tested flat neutral problem; they do not establish exact
discrete preservation.

| Material passport | Value |
| --- | --- |
| Origin | academic-research-suite / experiment-agent; implementation and numerical study |
| Version | neutral_abl_rans_m6_v2 |
| Status | CLOSED on 2026-10-30; accepted-resolution gates pass, higher-degree counterexample retained |
| Tested source parent | `88e2dafe2664bc6a7b185b645d7b382f7138f637` plus the source hashes in the evidence manifest |
| Branch | `codex/neutral-abl-rans-awf` |
| Authoritative environment | Linux, isolated worktree and Release build against the installed MPI libraries |
| Run directory | `/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m6-preservation` |
| Project allocation | At most 16 workers; numerical threads explicitly one |

## Equations, boundaries and measurements

The study uses the evolving URANS equations and the wall-strip contracts from
Milestones 1–5. It introduces no reference residual subtraction, field reset,
additional physical source, or frozen turbulence coefficients. The
`ExactSolution` functions supply passive diagnostics and the prescribed
reference inlet/boundary data.

The empty domain has streamwise length L=300 m, top Z=30 m and, in full 3D,
periodic span 20 m. The reference mesh has 12 geometric vertical layers,
velocity/turbulence degree 4, pressure degree 2, degree+4 GLL points per
direction, one streamwise element and one spanwise element. The copy test
uses two streamwise elements. Reference and small-perturbation 3D data are
span independent, while the solver and mesh are fully three dimensional.
Separate manufactured regressions exercise all three velocity components.

For the standard model, U(10 m)=10 m/s and z0=0.02 m. Its configured constants
give q=10/3, kappa=sqrt(0.48*1.3/q) and the compatible friction velocity.
The CLS reference uses u*=0.4 m/s, q=3.190675081624288,
kappa=0.440553905186659, H=2 m, P=1 m and N=4 m. The independently integrated
original-AWF matching relation fixes U(H)=5.780254367996082 m/s, hence
U(10)=7.241540166157099 m/s. Above H,
U=U(H)+(u*/kappa) log(z/H), k=q u*² and epsilon=u*³/(kappa z).
The high-Reynolds-number transported and total dissipation coincide.

The top carries compatible shear and dissipation fluxes and zero k flux.
The inlet prescribes the reference primary fields. The flat x+ outlet uses
total traction, zero scalar normal gradients and the coupled pressure marker.
For CLS, Sx=-(0.92-2/3)k and Sz=u*²; the standard normal traction is zero.
Streamwise-periodic columns provide a separate vertical-balance control.

The timestep is 0.5 s. Standard cases reach 150 s, or five L/U(10) times;
CLS cases reach 207.5 s, or 5.00873 times. The independent manufactured tests
use IMEX2 with dt=0.025, 0.0125 and 0.00625 s at the same final time.

`RANSReferenceMonitor=1` records initial projection error and every completed
step. The normalized L2 error is ||f-fref||/||fref||; the pointwise maximum
over quadrature points
uses |f-fref|/|fref|. Zero-reference velocities use u*, and zero-reference
stress/pressure uses u*², including the corresponding volume factor in L2.
The first-layer and upper-domain norms use the GLL masks H<=z<=N and z>N.
These masks are quadrature diagnostics, rather than cut-cell integrations.

The 0.1% gate applies to the maximum over time in both norms for U, k,
dissipation, viscosity and the relevant modeled stresses. Perturbations
intentionally exceed this gate initially and have a separate decay test.
The field audit also samples x/L=0.25, 0.5, 0.75 and 0.95 at 65 logarithmic
heights plus common physical heights 4, 8, 10, 16, 20, 25 and 30 m where
inside the domain. Derived stresses and viscosity are projected output fields;
that representation error is included in the reported errors.

## Implementation completed for the study

`RANSLinearSolver=MatrixFree` applies the current weak variable-diffusion and
mass operator with native backward transforms, derivatives and inner products.
It includes the same Robin strip operator used by the assembled residual
audit. Native static-condensation matrices provide frozen preconditioners,
with at most the two inertia values of constant-step IMEX2. Every physical
coefficient, nonlinear stress, source and wall datum remains current.
Normalizing each preconditioner right-hand side prevents native successive-RHS
magnitude smoothing from discarding small Krylov residuals. `Native` remains
the default solver. Some native Helmholtz collections omit variable mass in
their matrix action, and the full iterative 3D path lacks the required Robin
trace action; the new implementation therefore applies the weak action using
the elementary native operations.

`RANSOutflow` marks all velocity Neumann entries at the flat x+ outlet.
Its XML value is total traction (T-pi I)n, where
T=nu_eff(grad U+grad U^T)-R_nonlinear and pi=p/rho+2k/3.
The current pi*n completes the scalar Helmholtz load each nonlinear iteration.
The prescribed traction fixes the pressure constant. The approximate pressure
inverse uses a positive outlet Robin impedance
beta=sqrt(alpha/mean(nu_eff)) and a separate full-space mass inverse. This
impedance is only a preconditioner: the final multiplier and continuity keep
their full weak test space and have no added pressure boundary equation.
Domains with prescribed normal velocity everywhere retain the zero-mean gauge.

Both closures now expose the same derived viscosity, production, total
dissipation and complete modeled stress fields. The passive reference monitor,
common-height field sampling, small outlet/checkpoint fixtures and the
reproducible study/plot/evidence scripts accompany the solver changes.

Full IMEX1/2 restart stores the native solution and explicit derivative
histories at quadrature points, time labels, all primary coefficients including
pressure, element IDs, basis/points signatures and element coordinates.
Checksummed rank payloads and their XML manifest are bound to the imported
field state. Repartitioning maps by global element ID. A restart does not
repeat the Euler startup. Each publication reserves a fresh directory and
atomically publishes its payloads, preserving histories referenced by backup
fields. The strip has no independent history unknowns; its elimination uses
the complete restored predictor.

`RANSRestartMode=Auto` restores referenced history when present and requires
its files. Legacy fields without history use cold startup with a log message.
`History` requires complete metadata; `Cold` explicitly requests startup.
`RANSHistoryPath` relocates a sidecar directory. Keep field and history files
together, with the same constant timestep, IMEX order, geometry, quadrature,
model, wall parameters, boundaries and forcing. Boundary/forcing expressions
are not fingerprinted and must be kept consistent by the caller. Checkpoints
use step intervals; adaptive timesteps and truncated final steps are outside
the validated continuation protocol.

## Preservation, convergence and perturbations

| Peak relative error over five flows, percent | Standard, 2D/3D | CLS fixed/switched, 2D/3D |
| --- | ---: | ---: |
| U, L2 | 0.00000364 | 0.00000287 |
| k, L2 | 0.00000436 | 0.00000405 |
| dissipation, L2 | 0.00012084 | 0.00012084 |
| viscosity, L2 | 0.00009550 | 0.00113358 |
| viscosity, pointwise maximum | 0.0003032 | 0.0048784 |
| largest relevant stress pointwise error | 0.0061766 | 0.0020797 |

All six unperturbed standard/fixed/switched cases meet the target. The corrected
CLS switch stays at chi=1 with unchanged on/off thresholds 0.02/0.1.
The completed cases activate no positivity guards or cutoff/regime/range
fallbacks. Maximum weak-divergence residuals remain below 1e-8. Wall traction
has a maximum fractional error about 4.3e-8. In four-rank 3D, all 64 N samples
are owned away from the wall-owning rank, including ranks without wall faces.

Final dissipation L2 errors are approximately 0.0001114% in H<=z<=N and
0.0001108% above N for both models. The corresponding viscosity errors are
0.00009055%/0.00008680% for standard and 0.0010264%/0.0010152% for CLS.
The full U/k/stress band data are retained in the metrics. Maximum normal
velocity/u* is below 3.64e-8 in the six reference cases. Their maximum strong
divergence L2 norms are 2.51e-7 in 2D and 1.01e-6 in 3D; strong divergence
need not vanish pointwise in the mixed approximation. Final corrected strip
production/dissipation means are about 0.461092611/0.461092609 m²/s³.

Fixed-H refinement from 6 to 12 to 24 layers gives final L2 rates
4.97/4.99 for dissipation, 4.04/4.01 for viscosity and 4.98/5.00 for U.
Degree 3/4/5 reduces the 2D errors. Full-3D coarse-layer and degree-three
controls also finish five flows; degree-five results are collected separately.
Halving dt to 0.25 s in 2D and 3D changes the nearly stationary result by
much less than its spatial representation error. The unsteady manufactured
sequence independently gives standard primary-field rates 1.97–2.13 and
CLS rates 1.91–2.22, approaching second order. Pressure rates are about 2.00.

Changing H changes the physical matching problem. Recomputing its anchor gives
U(10)=7.1971536, 7.2415402 and 7.2872273 m/s for H=1, 2 and 4 m;
effective z0 is 0.00360927, 0.00343707 and 0.00326840 m. Each anchored case
meets the engineering target. Their common-height samples support comparison
without conflating changing physical H with mesh refinement. The additional
molecular-diffusion case also finishes five flows with small finite errors;
the strict analytical compatibility statement remains the high-Re limit.

The smooth perturbation uses an analytically divergence-free streamfunction
and a separate k amplitude 0.01. Finite-degree projection need not preserve
its pointwise divergence exactly. Both models recover after five flows:

| Final normalized L2 errors, 2D/3D approximately | U | k | dissipation | viscosity |
| --- | ---: | ---: | ---: | ---: |
| Standard | 2.98e-6 | 2.58e-6 | 5.07e-6 | 1.40e-6 |
| CLS switched | 1.70e-6 | 2.32e-6 | 4.78e-6 | 1.74e-5 |

The initial domain k error is 0.0035502; its decay exceeds a factor of 1000.
This small departure remains on the correction plateau; Milestone 5 supplies
the separate switch-off and recovery tests.
This tests the selected smooth departures in a flat domain. It does not
establish stability for arbitrary disturbances or terrain flows.

### Higher-degree full-3D robustness limit

The additional degree-five, four-rank 3D experiment does not supply another
preservation pass. At dt=0.25 s and linear tolerance 1e-13, it develops growing
strong divergence and reaches the coupled-stage iteration limit after its last
completed t=172.75 s. Its strong divergence L2 norm reaches 0.0306 despite a
weak-divergence norm below 1e-8. Tightening the linear tolerance to 1e-15 and
raising the iteration cap delays the growth; the nonlinear and continuity
targets remain 1e-8. The tighter repeat exceeds the 0.1% pointwise viscosity
target at t=165.25 s (3.98888 flows), with measured error 0.102706% and
strong divergence L2 norm 0.00165. It is stopped at that measured failure:
subsequent recovery could not restore the maximum-over-time criterion.
Its exact command, fixture/binary hashes, threshold, time and termination
receipt are retained. It does not complete or pass the five-flow test.

The accepted degree-four cases and 2D h/p convergence are separate findings.
The full-3D coarse-to-degree-four controls improve, but this additional result
prevents a claim of uniform long-time convergence with increasing degree in
3D. The cause of the growing mode has not been isolated between the numerical
method and closure. Full-3D degree five at these settings needs further
stability work before application; tighter linear precision alone is not a
demonstrated general remedy.

## Controls, budgets and the numerical column

The zero top dissipation-flux control completes five flows but reaches peak
L2 errors 2.8144% in dissipation and 11.8579% in viscosity. The inconsistent
kappa control develops large drift and stops at the iteration limit. The
original-AWF inlet/outlet cases in both dimensions stop after the last completed
t=9 s, far short of five flows. Its initial mean production 0.461092609 and
mean dissipation 0.774943967 already expose the specified source mismatch.
A stopped nonlinear solve is a measured numerical failure; it does not by
itself establish physical blow-up or long-time behavior. At dt=0.25 s, the
original periodic controls reach the coupled-stage iteration limit before
their first accepted step in both dimensions; they supply no completed-step
preservation or budget claim.

The compatible periodic CLS case has final integrated momentum and k budget
defects about -7.76e-8 and 6.62e-8; peak absolute defects are 4.99e-6 and
1.81e-6. Its `budget_available` flag is one. The existing natural-flux budget
does not recover inlet Dirichlet reactions, so inlet/outlet cases have
`budget_available=0`; their raw momentum column cannot be read as a
conservation defect. Subtracting the independently known reference inlet
normal-stress reaction, 3.6212035 in 2D and 72.4240701 in 3D, leaves final
raw-column differences about -1.37e-6 and -2.91e-5. General perturbed-inlet
reactions have not been recovered.

The periodic corrected column is continued for another 10.0054 flows, reaching
t=622 s (about 15 flows in total). Its final stage change is 6.87e-9 and weak
divergence norm is 6.61e-10. Piecewise degree-four reconstruction from interior
samples reproduces its saved fields within 1.34e-14 absolute error; derived
DG traces are fitted from their element interiors. The reconstructed primary
fields are copied into two streamwise elements in 2D and full 3D, with inlet
values and total outlet traction derived from that numerical column.

Both copies complete 41.5 s, or 1.00175 flows, at dt=0.1 s. Peak normalized
L2 changes relative to the copied column are at most 2.59e-8 for U,
1.41e-7 for k, 1.82e-7 for dissipation and 6.07e-7 for viscosity.
This separate finite-time consistency test retains small settling, projection
and boundary-trace errors. It is neither a machine-exact stationary solution
claim nor the five-flow analytical-profile test reported above.

## Verification evidence and reproduction

The final isolated solver passes eight CTests and the 27-case standard,
27-case original-CLS/AWF and 34-case corrected-AWF regression suites. Native
and matrix-free outlet comparisons cover both models in 2D and four-rank 3D;
the largest complete-field absolute differences are 2.41e-10, 4.68e-9,
8.98e-10 and 6.71e-9 for standard-2D, CLS-2D, standard-3D and CLS-3D.
Separate pilots check current assembled stage residuals and complete fields.

The 38-case continuation/manufactured suite verifies standard and CLS serial
to four-rank continuation, full-3D four-to-two-rank continuation, first-step
checkpoint recovery, IMEX1 and IMEX2, repeated publication without overwriting
backup histories, and rejection of changed model/timestep/order, bad checksums
and missing payloads. The temporal sequence above uses this same executable.
The largest complete-field restart difference is 1.37e-10 absolute.

The optional profile sampler initially transformed local coordinates twice.
Only the audit tool was rebuilt after correcting that interpolation. Saved
solver fields were re-audited; time and complete-field error/comparison JSON
must be identical, while profile CSVs are regenerated. Original reports and
an explicit refresh receipt remain on Linux. The solver hash is unchanged.
A fresh reproduction uses the corrected sampler directly.

Final compact plots, profile metrics, source/binary hashes, environment,
case counts and report hashes are in
[the M6 evidence directory](neutral-abl-rans-awf-milestone-6-figures/).
The parent plus recorded source hashes identify the tested patch before the
completion commit. Large fields, restart payloads and logs remain on Linux.
The primary Linux checkout, shared build and installation are unchanged.
All 153 recorded source hashes match the Mac delivery checkout. The tested
solver SHA-256 is `aa59534b39daf85415b73bc0bbb4b8e821c529995a10df2b68559b18e0c8f298`;
the corrected audit SHA-256 is
`c1f91be0c170dab3cbe03b2f493e5dc6abb4f50826ab925817650161615d28b2`.
The environment is GCC 11.4, CMake 3.22.1, Open MPI 4.1.2 and Python 3.10.12;
package versions and every cohort's report hash are in the manifest.

The closure manifest checks six degree-four reference targets, both models'
2D/3D perturbation decay, 2D spatial refinement, temporal order, numerical-column
copies, executable consistency and the regression reports. It records the
higher-degree experiment separately. Its `passed` status means these closure
gates passed, rather than that every compared wall policy or resolution passed.

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/m6-reproduction
mkdir -p "$RANS_RUN/tmp"
export TMPDIR="$RANS_RUN/tmp"
cmake -S solvers/IncNavierStokesSolver/Tests/RANS -B "$RANS_RUN/build" \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$RANS_RUN/build" --parallel 4
taskset -c 112-127 bash solvers/IncNavierStokesSolver/Tests/RANS/run_m6.sh \
  "$RANS_RUN/build" "$RANS_RUN" 2
```

The reproduction wrapper runs phases sequentially and permits at most four
study jobs of four MPI ranks each; its default uses at most eight workers.
The recorded run overlaps independent phases within 16 workers total, uses
at most four ranks per case, and sets OMP/OpenBLAS/MKL/BLIS threads to one.
Build parallelism is explicit and does not overlap the 16-worker allocation.

Supported scope remains full-dimensional CG on flat aligned layers with the
specified high-Re closures and wall policies, constant-step IMEX1/2 and a flat
streamwise traction outlet. Backflow/curved outlets, Homogeneous1D, terrain,
thermal or stratified ABL, low-Re CLS and empirical terrain calibration need
separate implementation or validation. Milestone 7 remains a separate delivery
and documentation task.
