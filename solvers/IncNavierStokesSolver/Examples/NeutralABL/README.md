# Neutral atmospheric boundary layer with URANS

This example prepares complete, self-contained XML sessions for standard
k–epsilon and the high Reynolds number Craft–Launder–Suga (CLS) cubic model.
Both closures can use the original rough-wall analytical wall function (AWF), its fixed
equilibrium correction, or the correction with active departure sensors.
The preparation scripts use the same fixture generator as the verified study.
They need Python 3.9 or newer and the tracked source tree; no literature files,
NumPy or SciPy are required to prepare an example.

The accepted benchmark uses degree 4 in full 2D or 3D. The standard log-law
control and corrected CLS profiles remained below the 0.1% engineering drift target for five flows
in [Milestone 6](../../../../docs/plans/neutral-abl-rans-awf-milestone-6.md).
The [standard-AWF extension](../../../../docs/plans/standard-kepsilon-equilibrium-awf.md)
adds fixed/switched standard-model AWF preservation and 2D refinement evidence.
The original CLS/AWF drifts. A degree 5, four-rank CLS 3D trial develops a growing mode
and fails the target even with tighter linear solves. Higher degree is not a
demonstrated remedy. The validated flat-domain scope does not include terrain,
thermal/stratified flow, curved/backflow outlets or Homogeneous1D.
The separate [terrain-wall extension](../TerrainHill/README.md) adds a plane,
compact hill and spanwise extrusion with physical normal-coordinate sampling,
integrated AWF strip continuity and a prescribed upstream correction ramp.
Its short pilots and C1-foot limitations are documented separately from the
five-flow flat-domain preservation benchmark.

## Build and run a short case

Run these commands on Linux from a checkout of this branch. Use a fresh run
directory on the large volume. The external build links an existing Nektar++
MPI installation and builds the production incompressible-solver source list
without replacing that installation.

```bash
RANS_SOURCE=$(git rev-parse --show-toplevel)
RANS_EXAMPLE="$RANS_SOURCE/solvers/IncNavierStokesSolver/Examples/NeutralABL"
RANS_TESTS="$RANS_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS"
RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/neutral-abl-rans-awf/tutorial
RANS_BUILD="$RANS_RUN/build"
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
mkdir -p "$RANS_RUN/tmp"
export TMPDIR="$RANS_RUN/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$RANS_TESTS" -B "$RANS_BUILD" \
  -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$RANS_BUILD" --parallel 4
ctest --test-dir "$RANS_BUILD" --parallel 1 --output-on-failure

python3 "$RANS_EXAMPLE/generate.py" "$RANS_RUN/cls/neutral.xml" \
  --model cls_switched --steps 4 --checkpoint 2 --linear-tolerance 1e-15
cd "$RANS_RUN/cls"
"$RANS_BUILD/IncNavierStokesSolver" neutral.xml > solver.log 2>&1
"$RANS_BUILD/RANSFieldAudit" neutral.xml neutral.fld neutral.fld profiles.csv > audit.log
```

The four steps are a startup example, not a five-flow preservation test.
The generator prints the physical reference, duration and XML SHA-256.
It embeds the mesh, initial/reference functions and boundary conditions.
Both preparation commands require an `.xml` output path; the continuation
command also refuses to overwrite its original session or field.
Output files are written in the solver's working directory. Run from the case
directory so diagnostics, checkpoints and history sidecars stay together.
`audit.log` contains the complete-field JSON audit; `profiles.csv` samples four
downstream stations at common physical heights.

This project's limit is 16 CPU workers total, including simultaneous builds
and MPI jobs. The commands above run sequentially, using at most four workers.
On the recorded server, optional `taskset -c 112-127` restricts the allocation;
use allowed CPU IDs on another machine. Account for other project work before
running multiple cases. Keep numerical-library threads at one per rank.

## Select the closure and wall policy

Pass `--model` to `generate.py`; all seven choices are available in 2D and 3D.

| Example choice | `TurbulenceModel` | `RANSWallModel` | Lower boundary |
| --- | --- | --- | --- |
| `standard` | `StandardKEpsilon` | `None` | Tangential `RANSLogLaw`, normal Dirichlet, zero k flux, prescribed epsilon |
| `standard_original` | `StandardKEpsilon` | `AWFOriginal` | Conservative original AWF strip |
| `standard_fixed` | `StandardKEpsilon` | `AWFEquilibrium` | Standard-model corrected strip with fixed chi=1 |
| `standard_switched` | `StandardKEpsilon` | `AWFEquilibrium` | Standard-model corrected strip with seven active sensors |
| `cls_original` | `CLSHighRe` | `AWFOriginal` | Conservative original AWF strip |
| `cls_fixed` | `CLSHighRe` | `AWFEquilibrium` | Corrected strip with fixed chi=1 |
| `cls_switched` | `CLSHighRe` | `AWFEquilibrium` | Corrected strip with seven active sensors |

The standard AWF uses the configured standard closure's
q=1/sqrt(Cmu) and kappa=sqrt((Cepsilon2-Cepsilon1)*sigmaEpsilon/q).
There is no Yap contribution or cubic stress. The original wall momentum and
production are retained; one capped dissipation distribution supplies the
mean k sink and both P/H epsilon values. Its admissibility check is the same
budget condition used by CLS, evaluated with the standard constants.
See the [standard-AWF derivation and validation](../../../../docs/plans/standard-kepsilon-equilibrium-awf.md).
For example:

```bash
python3 "$RANS_EXAMPLE/generate.py" "$RANS_RUN/standard-awf/neutral.xml" \
  --model standard_switched --flow-times 5
```

This generates the AWF-matched analytical reference. Its default U(10) is
approximately 7.141675 m/s, with effective z0 approximately 0.004417 m.
The matching anchor determines z0 from physical sand-grain h, H and current K;
the `standard` log-law control's independently chosen z0=0.02 m is a different
reference. Specifying an independent `RANSZ0` with `AWFEquilibrium` is rejected.

The solver selection is explicit:

```xml
<I PROPERTY="EQTYPE" VALUE="UnsteadyNavierStokes" />
<I PROPERTY="SolverType" VALUE="RANSVelocityCorrectionScheme" />
<I PROPERTY="TurbulenceModel" VALUE="CLSHighRe" />
<I PROPERTY="RANSWallModel" VALUE="AWFEquilibrium" />
<I PROPERTY="RANSResolvedDiffusion" VALUE="HighRe" />
<I PROPERTY="RANSLinearSolver" VALUE="MatrixFree" />
<I PROPERTY="Projection" VALUE="Galerkin" />
```

Use `u,v,k,eps,p` for standard in 2D and `u,v,k,eps_tilde,p` for CLS.
Insert `w` after `v` in 3D. In 2D, mesh (x,y) means physical (x,z), and solver
v is the normal velocity; in 3D, the axes coincide. Stress output names always
refer to physical x,y,z. The saved p is reduced pressure
Pi=p_mechanical/rho+2k/3, not mechanical pressure. All-normal-velocity domains
use a zero volume mean; an open traction outlet fixes the pressure constant.

The examples use degree 4 velocity/scalars, degree 2 continuous pressure,
eight GLL points per direction, 12 geometric vertical layers and constant-step
IMEX2. All fields share the same physical quadrature grid. The domain is
300 m long, reaches 30 m height, and has a 20 m periodic span in full 3D.
The generator's `--degree`, `--layers`, `--height`, `--dt`, `--perturbation`
and `--boundary outflow|periodic` options support controlled experiments.
`--backend native` selects the default native linear solver instead of the
matrix-free backend used in the preservation study.

The `standard` log-law reference has U(10 m)=10 m/s and aerodynamic z0=0.02 m,
with u* derived from the configured model kappa. All AWF choices use u*=0.4 m/s
and the independent original-AWF integral to fix U(H); they do not add a separately
chosen log-law offset. Its equivalent sand-grain roughness h=0.1 m is not z0.
The default strip locations are P=1 m, H=2 m and N=4 m. H is the start of the
resolved CG domain; no slip is enforced by the auxiliary wall solution below
H. Setting tangential velocity to zero at H would change the model.

## Configuration reference

These are solver settings, distinct from the analytical `ExactSolution`
functions used only for diagnostics. The example generator uses the original
model constants. If constants are changed by editing XML, recompute q, kappa,
the matching anchor and every inlet/top reference consistently; the inlet
expressions do not update themselves from the closure parameters.

| Setting | Core default or rule | Delivered example |
| --- | --- | --- |
| `RANSResolvedDiffusion` | `Molecular`; `HighRe` omits resolved molecular diffusion | `HighRe` |
| `Kinvis` | Positive physical viscosity, also needed in HighRe AWF | 1.5e-5 m²/s |
| `RANSLinearSolver` | `Native` or `MatrixFree` | `MatrixFree` |
| `RANSCmu`, `RANSCepsilon1/2` | 0.09, 1.44, 1.92 | Defaults |
| `RANSSigmaK/Epsilon` | 1, 1.3 | Defaults |
| `RANSYapCoefficient/Length` | CLS algebraic Yap: 0.83, 2.5 | Defaults |
| `RANSNonlinearTolerance` | 1e-9; state change and weak continuity | 1e-8 |
| `RANSPressureTolerance` | Nonlinear tolerance unless specified | 1e-8 |
| `RANSMaxIterations` | 200 standard; 800 CLS | 2000 |
| `RANSPicardRelaxation` | 0.5 standard; 0.1 CLS | Defaults |
| `RANSPressureRelaxation` | 0.8 standard; 0.25 CLS | Defaults |
| `RANSMinimumK/Epsilon` | Zero rejection thresholds; no floor is applied | Defaults |
| `RANSZ0`, `RANSHeightShift` | Standard log wall: 0.02, 0 | No independent z0 with corrected AWF; every AWF requires zero height shift |
| `RANSWallOrigin` | Physical wall coordinate, default 0 | 0 |
| `RANSStripPointP/Height/PointN` | Require 0<P<H<N; default 1,2,4 | Defaults |
| `RANSRoughnessHeight` | AWF equivalent sand-grain h, default 0.1 | `--roughness` |
| `RANSAWFCmu/CLength/SmoothIntercept` | Original AWF kernel: 0.09,2.55,10.7 | Defaults |
| `RANSAWFWeight` | `Switched` or `Fixed` | Chosen by `--model` |
| `RANSAWFFixedWeight` | Fixed chi in [0,1], default 1; zero recovers original | 1 |
| `RANSAWFOn/Off` | Required for switched correction; 0<on<off | 0.02/0.1, validated here only |
| `RANSAWFWeightU/DeltaU/K/Epsilon/Pressure/Convection/Normal` | Nonnegative sensor weights, default one; at least one positive | Defaults |
| `RANSAWFDirectionX/Y` | Incoming horizontal direction, default (1,0); Y=0 in 2D | Defaults |
| `RANSAWFMinimumRe` | sqrt(KP) H/nu gate, default 100 | 100 |
| `RANSAWFFailurePolicy` | `UseOriginal` default; `RequirePreservation` rejects an inadmissible reference | `RequirePreservation` |
| `RANSReferenceMonitor`, `RANSStageAudit` | Optional 0/1 diagnostics | Monitor 1, stage audit 0 |
| `RANSAWFPointDiagnostics` | Optional 0/1 per-rank strip receipts | 0 |

`RANSCmu` is the standard-model eddy-viscosity constant. CLS evaluates its
variable coefficient; `RANSAWFCmu` is a separate wall-kernel parameter.
For small-step AWF transients, use `--linear-tolerance 1e-15` and verify
coupled-stage convergence. At dt=0.001 s an original standard-AWF control
reached the iteration limit with 1e-13 but completed with 1e-15. The five-flow
degree 4 reference cases use 1e-13. A smaller timestep alone does not remove
linear-solver error from the pressure iteration.
The stronger 50% epsilon departure also reaches the iteration limit at
dt=0.001 s with 1e-15; switch-off recovery is validated at dt=0.025 s.

The correction requires its cap zc<=P and an admissible high-Re wall regime;
it never clips an inadmissible cutoff. Fixed chi=0 recovers the original AWF.
The switched policy senses velocity matching, the P/N increment, kN, epsilonN,
pressure minus body acceleration, convection and normal velocity separately.
Thresholds 0.02/0.1 describe this benchmark, not terrain calibration.
`Molecular` adds resolved molecular diffusion to either high-Re closure;
it does not implement a low-Re CLS model or its dissipation split.

## Boundary data and flux signs

For either closure with AWF at the lower interface, tangential velocity and k need Robin tags
matching the selected wall policy. Normal velocity and transported dissipation
need Dirichlet entries; the coupled strip determines their current traces.
Use a positive initial dissipation value. Every nonperiodic pressure entry
uses the `RANSPressure` Neumann marker. The generated 2D corrected wall is:

```xml
<R VAR="u" USERDEFINEDTYPE="AWFEquilibrium" VALUE="0" PRIMCOEFF="0" />
<D VAR="v" VALUE="0" />
<R VAR="k" USERDEFINEDTYPE="AWFEquilibrium" VALUE="0" PRIMCOEFF="0" />
<D VAR="eps_tilde" VALUE="Ustar^3/(Kappa*y)" />
<N VAR="p" USERDEFINEDTYPE="RANSPressure" VALUE="0" />
```

The top uses tangential traction +u*², zero normal velocity and zero k flux.
The reference dissipation flux is -u*⁴/(sigma_epsilon Z), not zero. At a flat
x+ outlet, tag every velocity Neumann entry `RANSOutflow`; its XML value is
the total traction (T-Pi I)n. With
T=nu_eff(grad U+grad U^T)-R_nonlinear, the current Pi n is supplied internally.
For CLS, normal traction is -(0.92-2/3)k and vertical traction is +u*²;
for standard the normal traction is zero. Scalars have zero normal gradients.
An ordinary scalar Neumann value is its diffusion-gradient flux; do not
multiply it by the diffusivity a second time.

## Full 3D and longer comparisons

```bash
python3 "$RANS_EXAMPLE/generate.py" "$RANS_RUN/cls3/neutral.xml" \
  --model cls_switched --dimension 3 --steps 2 --dt 0.05 \
  --perturbation 0.01 --linear-tolerance 1e-15
cd "$RANS_RUN/cls3"
mpirun --bind-to none --map-by slot -np 4 \
  "$RANS_BUILD/IncNavierStokesSolver" neutral.xml > solver.log 2>&1
```

Use `--flow-times 5 --dt 0.5` and omit `--steps` in a new directory for the long
comparison. Omit `--perturbation` for the unperturbed analytical case. This is
300 steps for `standard`, 421 for the standard AWF choices, and 415 for CLS.
With `--perturbation 0.01`, a smooth streamfunction velocity perturbation
and a k perturbation are added to the initial state; the reference stays
unperturbed. These full-3D reference/perturbation data are span independent.
Independent manufactured regressions exercise all three velocity components;
the study does not establish general three-dimensional disturbance stability.

## Continue from a field or checkpoint

Keep the entire native field file/directory and its referenced
`.rans-<completed-step>[.<generation>]` history directory. MPI field output
requires all partitions and Info.xml. Copying only one rank's file is not a
restart. Histories contain solution and explicit derivative quadrature arrays,
primary coefficients including pressure, time labels and geometry/basis
signatures. They use a checked native-endian 64-bit-double schema.

The short 2D example above finishes at t=2 s. Prepare four additional steps:

```bash
python3 "$RANS_EXAMPLE/restart.py" "$RANS_RUN/cls/neutral.xml" \
  "$RANS_RUN/cls/neutral.fld" "$RANS_RUN/resume/resume.xml" --steps 4
cd "$RANS_RUN/resume"
mpirun --bind-to none --map-by slot -np 4 \
  "$RANS_BUILD/IncNavierStokesSolver" resume.xml > solver.log 2>&1
```

To continue from the first step-interval checkpoint at t=1 s, use
`neutral_1.chk` instead of `neutral.fld` and `--steps 6` to reach the same
t=4 s. The helper preserves the model, mesh, quadrature, timestep and forcing;
it changes InitialConditions, the additional NumSteps and restart mode,
removes FinalTime, and follows the new field's history metadata. These
generated cases are self-contained. For other XML sessions, check any external
mesh/function paths before moving a session into a new directory.

`--mode History` is the helper default and requires a complete continuation.
Core `RANSRestartMode=Auto` restores referenced history when present and uses
cold startup for legacy fields without history metadata. Missing/corrupt
referenced history fails; Auto does not silently downgrade that case.
`--mode Cold` explicitly requests native startup without multistep history,
so it is not an uninterrupted-run equivalence claim. Stored physical
configuration checks still apply. A moved complete sidecar can be selected
with `--history-path /absolute/new/history-directory`.

Keep the mesh, basis, quadrature, constant dt, IMEX order, closure/wall
parameters, boundary definitions and forcing consistent. MPI repartitioning
is supported. Boundary/forcing expressions are not fingerprinted by the core;
the caller must preserve them. Use step-based `IO_CheckSteps` scheduling as
in the supplied examples. Time-based/filter scheduling, adaptive timesteps
and shortened final steps are outside the validated continuation protocol.

## Read the results

| Output | Interpretation |
| --- | --- |
| `.fld`, `.chk` | Primary fields and projected modeled viscosity, production, dissipation, six stresses and quadratic/cubic norms |
| `.rans.csv` | Stage iterations, state/pressure changes, strong/weak divergence, extrema, guards, cache and Krylov counts |
| `.awf.csv` | Standard/CLS strip inventories, traction, production/dissipation, switch, admissibility and MPI ownership |
| `.rans-errors.csv` | Initial and per-step independent reference errors; diagnostics do not reset the solution |
| `.awf-points-r*.csv` | Optional per-rank inputs/sensors with `RANSAWFPointDiagnostics=1` |
| `.rans-residual.csv` | Optional assembled stage-coefficient residuals with `RANSStageAudit=1` |

Reference L2 errors use the reference norm. Zero-reference velocities use u*,
and zero-reference stress/pressure uses u*². Pointwise maxima are sampled at
quadrature points, not exact continuous L-infinity norms. The first/upper-layer
bands split at N and include polynomial representation error. Weak continuity
does not imply pointwise divergence is zero: inspect both diagnostics.
The degree 5 counterexample demonstrates why this distinction matters.

Check `budget_available` before interpreting wall budget columns. Inlet
Dirichlet reactions are absent from the natural-flux inventory, so inlet/outlet
cases have flag zero; their raw momentum column is not a conservation defect.
The periodic inventory has flag one. Modeled `rans_R_*` stresses are distinct
from covariance of resolved unsteady fluctuations.

For the field-level continuation comparison with the four-step CLS prefix:

```bash
python3 "$RANS_EXAMPLE/generate.py" "$RANS_RUN/full/neutral.xml" \
  --model cls_switched --steps 8 --linear-tolerance 1e-15
cd "$RANS_RUN/full"
"$RANS_BUILD/IncNavierStokesSolver" neutral.xml > solver.log 2>&1
"$RANS_BUILD/RANSFieldAudit" "$RANS_RUN/full/neutral.xml" \
  "$RANS_RUN/full/neutral.fld" "$RANS_RUN/resume/resume.fld" \
  > "$RANS_RUN/restart-comparison.log"
```

Do not confuse a cold restart with this full-history comparison.

## Reproduce and extend the evidence

The [test index](../../Tests/RANS/README.md) lists kernels, manufactured
solutions, wall assembly, MPI, history and the complete study. After building,
the delivery checks exercise these public scripts and restart commands:

```bash
python3 "$RANS_TESTS/validate_m7.py" --build "$RANS_BUILD" \
  --output "$RANS_RUN/delivery-check"
```

For the full numerical study in a separate run directory:

```bash
taskset -c 112-127 bash "$RANS_TESTS/run_m6.sh" \
  "$RANS_BUILD" "$RANS_RUN/full-study" 2
```

Study status `recorded` means the comparisons were collected, including
failed controls. The closure collector's `passed` status applies to its
degree 4 acceptance gates and retains the higher-degree failure separately.
The [M6 report and evidence](../../../../docs/plans/neutral-abl-rans-awf-milestone-6.md)
identify the exact settings, source/binary hashes and measured limits.
The [M1–M7 plan](../../../../docs/plans/neutral-abl-rans-awf-proposal.md) links
the equations, conservative strip contract and development evidence.

The original model sources are Craft, Launder and Suga (1996),
[cubic eddy-viscosity model](https://doi.org/10.1016/0142-727X(95)00079-6);
Suga, Craft and Iacovides (2006),
[rough-wall AWF](https://doi.org/10.1016/j.ijheatfluidflow.2006.03.011);
and Richards and Norris (2011),
[neutral-ABL boundary conditions](https://doi.org/10.1016/j.jweia.2010.12.008).
The supplied equilibrium correction is implemented through the tracked
milestone specifications; the example does not require untracked paper files.
