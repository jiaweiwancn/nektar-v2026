# M7 — Admissibility and perturbation diagnosis

## Material Passport

- Mode: experiment-agent run and reproducibility validation.
- Scope: unchanged revised-note CLS and equilibrium-preserving AWF.
- Branch: `codex/cls-epawf-revised`; baseline M6 `f41c7e46986928d380866f16df284c791b4d5dac`.
- Authority: approved implementation plan, frozen M0 contract, revised note.
- Execution: isolated Linux source/build and run directories, Tailscale SSH
  `jay@100.105.247.110`; numerical library threads one; compiler jobs 16.
- Closure date: 2026-11-03.
- Status: **M7 CLOSED — faithful-model implementation and diagnosis VERIFIED.**
  Original-model robust stability remains **failed**; no M8 variant is adopted.

## Experiment contract fixed before execution

The original CLS and scalar AWF source and the frozen note/contract stay
unchanged. Diagnostics are opt-in (`RANSWallDiagnostics=1`). The accepted
state must have nonnegative sampled inner k and strictly positive sampled K,
outer k and outer epsilon. A rejected candidate is recorded as a rejection,
never a valid trajectory. Covariance and tangent signs are measured independently.
No clipping, filtering, SVV, constitutive retuning or positivity intervention
is introduced. There is no claim about polynomial positivity everywhere.

Every initial/accepted native stage is audited at the actual split ray source
quadrature, outer source quadrature, fixed prescribed-source quadrature,
physical mass quadrature, both matching traces, physical K and transition
samplers, native GLL reconstruction and a denser equispaced reconstruction.
Each element uses `2*p+3` points per coordinate for the dense set.
The independent Jacobi oracle reconstructs every accepted k polynomial on that
set. A deliberately negative between-node polynomial must be rejected even
though native GLL and source GL subsets are positive.

The native nonlinear campaign separately perturbs velocity, k and epsilon,
then their combination, a near-interface disturbance and a disturbance crossing
the prescribed upstream transition on a cosine-foot map. Amplitudes, wavelengths,
dt, p and mesh are changed independently. Isolated flat verification uses
eta=1 with the original frozen defect; unanchored paired controls and all
transition/terrain runs use eta=0. An unperturbed native benchmark spans five
flow-through times; bounded nonlinear perturbations are short diagnostic runs,
not arbitrary-time stability tests. A separate matching-height study fixes
physical roughness, nu, rK and rm and reports the resulting wall-model dependence.

These bounded diagnostic families use separately recorded small domains
(flat L=8 m/Z=8 m and transition L=6 m/Z=4 m). The transition samplers are
ra=2.5 m/rb=3.5 m. They do not replace the frozen M4 reference fixture.
Within each dt/p/mesh refinement, the physical domain and matching/sampling
inputs remain fixed. The H study is explicitly a separate model study.

The three-component solenoidal transverse test uses a clearly labelled reduced
linearized Fourier/native-CG system supported on `4<z<12`, strictly above H=2.
The compiled unchanged constitutive law supplies its tangent. The complete
transverse stress block and base advection are implicit in native IMEX1/2;
the explicit prescribed body term is zero. Independent indexed stress,
quadrature, matrix exponential and time recurrences check this experiment.
Wavenumbers 1/2/4, amplitudes 1e-5/1e-4/1e-3 and independent dt/p/h changes
test the adverse growth described in the note. This is not a nonlinear 3D
terrain trajectory.

Acceptance remains the frozen contract: scaled coupled residual <=1e-9,
isolated identities <=1e-10, MPI physical differences <=1e-8. Strong divergence
and reconstructed traction/flux jumps are measured separately from weak
constraints. Native stage energy uses the actual physical mass operator,
the unmasked velocity/constant-k test, prescribed and roughness work,
essential-boundary reactions and the separate anchor work. Trapezoidal
integrals of sampled accounting defects are labelled as such. They do not
replace the native IMEX storage identity. Negative principal diffusion,
growth, or loss of sampled admissibility leaves robust-model qualification
failed even when implementation checks pass.

## Verified native and independent results

The final complete native campaign at
`dea75dab62fb59f5264868f3d2033d13e0b4ef57` passed 5,849 checks: 25 successful
reports and three required initial-state rejections. It covers 135 audited
accepted initial/stage records. The compiled solver/kernel physics remains
`bf202be017bab7aceeedf9a48d86bcfeeca0580f`; no C/C++ changes occur between these
revisions. The original forward campaign and its 2,515 checks are retained
separately. The audit-only
independent final replay passed 9,852 checks at
`dea75dab62fb59f5264868f3d2033d13e0b4ef57`, after proving that all 99 runtime
images and every original numerical source other than the named M7 audit
drivers remained unchanged. The stricter transverse qualification passed
101,736 checks at `dea75dab62fb59f5264868f3d2033d13e0b4ef57`; its compiled
constitutive and time-integration code is the unchanged `bf202be` build.

| Native result | Measurement |
|---|---|
| Largest scaled accepted-stage coupled residual | 1.755e-12 |
| Largest scaled weak/direct energy discrepancy | 1.566e-13 |
| Largest independent final energy/traction audit discrepancy | 1.620e-13 |
| Dimensionless dt differences at fixed p/mesh | 4.18885e-8 and 1.04769e-8; ratio 3.99816 |
| MPI parity | Complete physical fields and diagnostics pass the 1e-8 gate |
| Perturbed restart | Complete coupled state and native histories pass |
| Five-flow stationary benchmark | 25 actual native IMEX2 steps, five declared flow-through times; all states audited |

Both the source/volume sets and the actual boundary/trace transfer samplers
are checked. The latter use their own native GLL tangential points at
physical rK/rA/rB and H. They are not assumed to be subsets of the source
Gauss or dense reconstruction grids. The independent Jacobi reconstruction
checks every accepted k polynomial on the declared dense finite grid.

The final temporal comparison uses the frozen per-field physical mass norms,
including each field scale and domain volume. Pressure-block differences are
reported separately. The first native receipt's unscaled combined differences
are retained as historical output; the dimensionless final audit supplies the
comparison used here.

The focused dense-only counterexample keeps every other audited
admissibility set positive, including native GLL reconstruction, physical
samplers, both interface traces and source/mass quadrature. One additional
dense point gives k=-5.10508e-4, so the candidate is rejected and no initial
state is accepted. Its initial-condition support ends at the first radial
element, leaving the K sampler and outer epsilon unchanged. The `max` in
that initial support expression is neither a transported-field limiter nor
a closure modification. A separate invalid-K-sampler control preserves the
original failure before native history exists. These controls establish the
rejection mechanism; they do not establish positivity invariance or a bound
over all polynomial points.

`energy` contracts the actual unmasked weak stage loads with velocity and
constant k. It separately reports essential-boundary reactions and the fixed
anchor work. Independent dimensional AWF primitives, indexed CLS stress,
Jacobi fields and physical mass matrices reconstruct the direct Fixed-mode
components. The unforced Upstream preflight additionally reconstructs all
seven physical sensor groups, chi, initial mixed-chi source budgets and the
direct energy terms. The zero-velocity manufactured transition receives a
clearly labelled mapped mass/storage/work audit; its prescribed manufactured
forcing is not guessed by the flat neutral-data oracle. The current M6
regression independently checks its complete coupled weak operator.

Strong divergence and reconstructed total-traction/k-flux jumps are measured
separately from converged weak pressure/flux constraints. They are not forced
to zero by a reporting convention. A small coupled residual does not imply
pointwise continuity, pointwise matching of all reconstructed gradients, or
perturbation stability.

## Adverse transverse dynamics

The original tensor gives a0/nu_t=-0.964557906564225, with zero first-order
turbulence-variable forcing in the transverse subsystem. The finest native
stress difference agrees with the independently indexed derivative within
1.74e-12; complete native weak matrices agree within 3.07e-12. Covariance
remains positive at the actual finite-perturbation samples.

| Wavenumber xi_x [1/m] | Native energy ratio at T=0.05 s, p=7/eight radial cells |
|---|---:|
| 1 | 1.09258702038 |
| 2 | 1.64509274569 |
| 4 | 8.62184899621 |

The duration is 0.01 turnover times using H/u_star=5 s. It spans
0.0560/0.1120/0.2241 wavelength flow-through times for the three modes.
Independent p and mesh refinement converge to this growing solution. Against
the independently assembled matrix exponential, step halving reduces IMEX1
field errors by 2.025/2.014 and IMEX2 errors by 3.916/3.960. Amplitudes
1e-5, 1e-4 and 1e-3 retain the linear energy ratios. The restricted chi/eta
controls leave the outer operator unchanged. These are explicitly reduced
linearized Fourier/native-CG results, not nonlinear terrain trajectories.

The observed growth is a failed physical stability gate. A successful
nonlinear stage, a positive sampled covariance and a stationary anchored
benchmark do not repair it. No remedy is adopted in this milestone.

## Matching-height dependence

This is a wall-model study separate from mesh/p refinement. Physical nu,
roughness h=0.1 m, rK=rm=1 m, u_star and the closure coefficients remain fixed.
The native matching trace and the independent fully turbulent formula for
z0,c agree to better than 1e-10, without prescribing or retuning z0,c.

| H [m] | Matching velocity [m/s] | Output z0,c [m] |
|---|---:|---:|
| 1.5 | 5.50041976422 | 0.00350833606268 |
| 2.0 | 5.78025436800 | 0.00343706659146 |
| 2.5 | 5.99746331957 | 0.00338221660291 |

This verifies `M7.matching_height_dependence` and retains the revised note's
height dependence. Exact preservation at one matching height does not impose
the same aerodynamic roughness at another height.

## Unforced transition feasibility

The extra physical-flow experiment removes all manufactured/body/scalar
functions, retains the actual neutral inflow, prescribed mask and original
state thresholds, and uses eta=0. Its paired one-step preflights pass. A
localized analytical streamfunction pulse, k and epsilon perturbations are
used; the ordinary boundary/trace lifting remains coupled.

The first long control stops at the 35-iteration limit after 0.03 s. Increasing
only that budget to 80, with the same 2e-12 nonlinear tolerance, extends the
accepted prefix to 0.055 s, before every admissible Newton trial is rejected
at the next step. All recorded accepted states remain positive. Neither stop
is a demonstrated positivity failure or a continuum-instability proof, and
neither supplies a qualified transition-passage trajectory. The independent
transverse experiment supplies the demonstrated constitutive failure.

Additional paired/dt observations and the earlier-regression closure are
recorded separately from trajectory qualification. Both controls and pulses
with dt=0.005/0.0025 s reach 0.055 s and then stop, leaving 12/23 actual
accepted records per trajectory. At common accepted epochs, the perturbation
norm centroid drifts downstream by about 0.264/0.360 m, but it does not pass
the transition. The final common-time paired field difference is 0.346;
no temporal-order or full-duration convergence claim is made.

The stopped prefixes also show negative outer modeled covariance eigenvalues:
about -0.138 for the coarse control and -0.152 for the finer control at their
last accepted state. Independently indexed stress and eigenvalue evaluation
confirm the recorded values. These are finite-resolution observations,
separate from turbulence positivity and from the exact-reference transverse
failure witness. The coarse projected initial inner stress can already be
indefinite; it must not be conflated with the positive continuum reference
covariance. No failed or extrapolated state is used for physical validation.

The final closure receipt distinguishes verified implementation identities,
verified observations of stopped trajectories and the failed robust-model
gate. All 31 earlier registered tests pass at the compiled-physics revision,
including full M3/M4/M5/M6 qualification, MPI and restart checks. Together with
the two M7 driver qualifications, the evidence covers all 33 registered gates;
this is combined source-bound coverage, not a single 33-test CTest invocation.

## Closure and reproducibility

The integration audit passes 241 checks, covers all 33 registered gates and
closes both frozen M7 verification owners: reference-only matching-height
dependence and full energy accounting. The latter is an accounting identity,
not a stability or nonnegative-energy bound. M7 is closed by the faithful
diagnosis required in the approved plan. The original-law growth witness and
the unqualified unforced-flow observations remain explicit failed physical
gates. M8 requires a separately identified approved model decision; M9 is open.

The primary isolated Linux source is synchronized to the qualified audit
revision. Incremental Release compilation succeeds with 16 compiler workers,
and all 99 runtime images remain byte-identical. The full fingerprints also
match loader bindings, numerical environment, installed NumPy/SciPy/SymPy
code and tools against the original campaign. The user Linux checkout remains
clean at `f4b70baec10016dd098005bc5a5507b87e4a57bf` on its original branch.
The final documentation-only commit is checked against this qualified
fingerprint; its receipt is retained in the Linux `m7-final-sync` directory.

The tests use one numerical thread per rank, a 12-rank CTest reservation and
at most three overlapping M7 ranks. Compilation is separate from numerical
tests. The combined allocation stays below the hard 72-worker limit.

The [evidence index](cls-epawf-revised-evidence/m7/README.md) links the final
native, independent, covariance/energy, initial-rejection and regression
receipts, equation-owner mapping, full raw-artifact hashes and exact runners.
Large fields, checkpoints, native histories and failed attempts remain under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised`, outside both source
checkouts. The [scientific figure](cls-epawf-revised-evidence/m7/m7-transverse.png)
and [PDF](cls-epawf-revised-evidence/m7/m7-transverse.pdf) use the final strict
transverse data and retain their plotting/data hashes.

Reproduce on Linux over Tailscale, in a fresh run directory, with the installed
MPI environment and the isolated `build/geometry-lib` overlay:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1
export LD_LIBRARY_PATH=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build/geometry-lib:$LD_LIBRARY_PATH
python3 /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m7.py --build /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build --output FRESH_NATIVE_DIRECTORY
python3 /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m7_transverse.py --build /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build --output FRESH_TRANSVERSE_DIRECTORY
```

The earlier-regression command is retained verbatim in
`regression-receipt.json`: CTest with `--parallel 12 -E ^RANSCLSEPAWFM7` and
one numerical thread. Solver-only iteration-budget follow-ups preserve both
the original 35-iteration stop and the 80-iteration results; no failed receipt
is replaced. The first build/pilot corrections and all subsequent raw evidence
are included in the artifact index. The note, frozen contract, original CLS
and scalar AWF hashes remain unchanged.
