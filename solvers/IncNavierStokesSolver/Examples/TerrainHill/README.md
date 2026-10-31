# Terrain walls and an upstream AWF transition

This example extends the [neutral-ABL URANS example](../NeutralABL/README.md)
to a plane, a two-dimensional cosine-squared hill, and its spanwise extrusion.
Standard k–epsilon and high-Re CLS support original and equilibrium-preserving
AWF. The compatible log law supports standard k–epsilon. The spatial AWF
correction finishes upstream of the hill; it is zero over the hill and wake.

The geometry and transition enable a later Cao–Tamura experiment campaign.
The supplied short runs are numerical functionality tests. Their approach
profile, equivalent sand-grain roughness, domain and span are not calibrated to
the wind-tunnel measurements. Two timesteps cannot determine an appropriate
transition width, a statistically stationary separation bubble, or experimental
accuracy. The implementation/validation record is
[the terrain plan](../../../../docs/plans/hill-wall-geometry-transition.md).

The [WG4 evidence](../../../../docs/plans/hill-wall-geometry-evidence/wg4.json)
records 23 hill pilots, five neutral controls, legacy regressions and ten passed
CTests. The [complete-field rank check](../../../../docs/plans/hill-wall-geometry-evidence/wg4-rank-check.json)
compares four/eight-rank 3D CLS checkpoints. The plan also records the C1-foot
maximum-error limitation and the native-inverse runtime/precision diagnostics.

## Prepare, build and run on Linux

Use the Linux environment and large-volume run directory specified in AGENTS.md.
Keep all OMP/BLAS threads at one; this project's limit is **16 workers total**.
The isolated external build also builds the source SpatialDomains library, so
XML partitioning retains full precision for curved edges/faces. It leaves the
shared installation unchanged. Its geometry library must precede installed
libraries in `LD_LIBRARY_PATH` (the validation runner and CTest do this).

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
RANS_SOURCE=$(git rev-parse --show-toplevel)
RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/my-terrain-study
RANS_TEST="$RANS_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS"
RANS_EXAMPLE="$RANS_SOURCE/solvers/IncNavierStokesSolver/Examples/TerrainHill"
mkdir -p "$RANS_RUN/tmp" "$RANS_RUN/cls-width2"
export TMPDIR="$RANS_RUN/tmp"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
cmake -S "$RANS_TEST" -B "$RANS_RUN/build" -DCMAKE_BUILD_TYPE=Release \
  -DRANS_BUILD_SOLVER=ON \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$RANS_RUN/build" --parallel 4
export LD_LIBRARY_PATH="$RANS_RUN/build/geometry-lib:${LD_LIBRARY_PATH:-}"
python3 "$RANS_EXAMPLE/generate.py" "$RANS_RUN/cls-width2/neutral.xml" \
  --kind hill --model cls_fixed --width-ratio 2 --layers 6 --dt 0.00002 --steps 2
python3 "$RANS_TEST/terrain_initial.py" "$RANS_RUN/cls-width2/neutral.xml" \
  --build "$RANS_RUN/build"
cd "$RANS_RUN/cls-width2"
mpirun --bind-to none -np 4 "$RANS_RUN/build/IncNavierStokesSolver" neutral.xml
```

The initial-field preparation is required for hill pilots. Its streamfunction
uses the normal-distance log profile through N and blends into a roof-compatible
profile by 2N. Near H/N the raw analytical tangential traces are independent of
wall foot, so zero raw initial transpiration follows from integrated continuity. The
subsequent wall condition solves for normal velocity. The 2D native boundary
projection and the 3D conforming reference tensor projection preserve boundary
polynomial traces. The 3D preparation checks shared face/edge coefficients and
span invariance; an independent oracle checks the analytical formula, geometric
reconstruction and turbulence positivity. Projection does not establish discrete
equilibrium, and the pilot includes the initial pressure/flow adjustment. This
is a synthetic start, not a reset during time advancement or a measured profile.

Generate an identical empty-floor control with `--empty`. Width ratios
`1,2,4,8` mean transition length divided by hill height. All cases use the same
streamwise mesh; the downstream endpoint stays at x=−0.20 m, before the upstream
hill foot x=−0.10 m. The starts are −0.24, −0.28, −0.36 and −0.52 m.
`--diffusion Molecular` retains resolved molecular transport; `HighRe` omits
it. Molecular viscosity remains positive in the AWF auxiliary model in both
modes. Finite resolved viscosity does not give exact preservation of the
inviscid log-law reference.

`--model` accepts `standard` (compatible log law), `standard_original`,
`standard_fixed`, `standard_switched`, `cls_original`, `cls_fixed`, and
`cls_switched`. `--dimension 3` produces full hex elements with periodic span;
it is not Homogeneous1D. `--kind plane --angle 0.35` rotates the entire plane
test, including profiles and boundary data. `--refinement 2` or `4` refines
streamwise columns and every resolved normal layer of the hill fixture. The matching height
and timestep are separate controls (`--height`, `--dt`).

## Geometry and boundary configuration

```xml
<I PROPERTY="RANSWallGeometry" VALUE="CosineHill" />
<I PROPERTY="RANSLinearSolver" VALUE="MatrixFree" />
<I PROPERTY="RANSWallModel" VALUE="AWFEquilibrium" />
<I PROPERTY="RANSAWFWeight" VALUE="Fixed" />
<I PROPERTY="RANSAWFSpatialWeight" VALUE="Ramp" />
```

The physical terrain is z=origin+height/2·[1+cos(pi·(x−centre)/halfLength)]
inside its support and z=origin outside it. Parameters are `RANSWallOrigin`,
`RANSTerrainHeight`, `RANSTerrainCentre`, `RANSTerrainHalfLength`; a plane uses
`RANSTerrainSlope`. In 2D, mesh y and solver v represent physical z and w.
In 3D, terrain is extruded along physical y. General terrain files and
non-extruded three-dimensional surfaces are not supported.

The resolved lower boundary must be the actual physical normal offset H of
the terrain, with P<H<N on the same normal ray. Use `RANSStripHeight`,
`RANSStripPointP`, `RANSStripPointN`. Heights are physical normal distances,
not global vertical offsets. CLS uses geometric closest-wall distance for its
Yap correction. Aerodynamic z0 and sand-grain roughness do not alter that distance.
The adapter checks interface coordinates/normals, uniqueness of normal rays,
positive metric J=1−curvature·eta, and |curvature·H|≤0.25 by default
(`RANSMaximumCurvatureRatio`). A folded or overly thick strip is rejected.

Use Robin conditions on **all Cartesian velocity components** at the terrain
wall, tagged `AWFOriginal`, `AWFEquilibrium` or `RANSLogLaw`. AWF k is Robin;
log-wall k is zero Neumann. Epsilon is Dirichlet; pressure uses `RANSPressure`.
The generator supplies complete working examples. The current vector solver
requires matching Cartesian velocity boundary types on each outer region and
identical scalar assembly maps. Fully prescribed inlet/top velocities, vector
traction outlets and periodic spans satisfy this restriction. Arbitrary mixed
Cartesian conditions on an inclined outer boundary do not.

The solver applies the tangential Robin operator and a wall-normal multiplier
in Cartesian coordinates. Let P be the L2 projector onto the assembled
continuous free scalar velocity trace, u=P1, a=1−u and c=(1,u). The conservative
constraint C f=P f+u(a,f)/c preserves the integral also at prescribed inlet/wall
corners. Enforce C(U_normal−W_H)=0; use its adjoint C* for the reaction and
beta C*C for augmentation. This gives an SPD velocity block and conserves
integrated strip mass. Independent pointwise and side-flux audits remain
necessary because the constraint is enforced in a finite-dimensional space.
`RANSNormalAugmentation` (default 2) controls the augmented iteration;
`RANSNormalRelaxation` (0.5), `RANSPressureRelaxation` (terrain default 0.1)
and `RANSPicardRelaxation` (terrain default 0.1) damp the coupled iteration.
The terrain default `RANSNonlinearSolver=NewtonKrylov` solves the coupled
velocity, turbulence, pressure and normal-reaction stage using native GMRES,
with positivity-constrained line search. Its CLS Jacobian probes hold the
active max-invariant branch fixed; accepted physical residuals always use the
original model. `RANSNewtonMaxIterations` defaults to 60. The legacy default
remains `Picard`. `RANSNewtonKrylovDimension` defaults to 240 (20–512);
`RANSNewtonDifferenceStep` defaults to 1e−6 for CLS and 1e−7 for standard
k–epsilon. Pressure and reaction scales balance the implicit inertia blocks;
physical residual tolerances remain separate. `RANSNormalDirectLimit=0`
selects the iterative trace projection instead of cached Cholesky factors.
The converged problem is independent of augmentation. The velocity block is
solved together; separately lagging its Cartesian cross terms breaks rotation
covariance. The legacy flat adapter remains the default (`RANSWallGeometry=Flat`) for
constant-wall regressions. Spatial ramps/properties require `Plane` (also for
horizontal ground) or `CosineHill` and its conservative normal coupling.

## Spatial transition and properties

Set `RANSAWFTransitionStart` < `RANSAWFTransitionEnd`. The latter must be
strictly upstream of the cosine hill foot. With r clipped to [0,1], the C2
mask is 1−10r³+15r⁴−6r⁵. Fixed mode uses mask·fixedWeight; switched mode uses
mask·sensorWeight. All seven existing sensors remain active inside this
envelope. The same effective weight selects the strip dissipation sink and
both P/H epsilon values. The original traction and production are retained.
Zero mask bypasses corrected-reference cutoff/high-Re admissibility checks
and selects original AWF exactly, including under `RequirePreservation`.

Optional stationary expressions are evaluated at the physical terrain foot:

```xml
<FUNCTION NAME="RANSWallProperties">
  <E VAR="roughness" VALUE="0.002+0.001*(1+tanh(x/0.1))/2" />
</FUNCTION>
```

AWF uses `roughness` (equivalent sand-grain height); the compatible log law uses
`z0`. Time-dependent and file-based properties are rejected. The corrected AWF
derives its effective aerodynamic roughness from its matching anchor, so an
independent z0 is not permitted. Align discontinuous property changes with wall
element boundaries; one-sided element-interior traces are used at shared endpoints.

AWF and CLS reject an independent `RANSHeightShift`; their delivered reference
uses unshifted physical height. Fitting a displaced measured inlet exactly is a
further modelling extension. The standard compatible log wall supports its
documented height shift.

Referenced expression parameters, geometry, transition, wall policy and normal
constraint controls enter the history-restart signature. Native history checks
also compare basis/quadrature and physical element coordinates. Restarting with
different MPI rank counts is supported; changed model/geometry/property/mask
data are rejected. The multiplier is a converged-iteration warm start, not a
physical history variable.

## Auxiliary-model limits and diagnostics

The original auxiliary AWF remains locally parallel. Curved coupling uses the
exact normal metric: strip volume V=H−curvature·H²/2, interface factor J(H),
metric-weighted production/dissipation, and the rotating-shear-basis contribution.
AWF retains integrated strip continuity:
J(H) W_H = -H d_s(U_P) - V(H) d_span(U_span,P).
The corresponding relation at P supplies the normal-velocity sensor. The
pulled-back scalar strip flux is differentiated on the native curved surface;
interface normal velocity is not set to zero. Independent GLL differentiation
and strip-side/interface-flux budgets audit this relation. Only the compatible
log wall uses zero normal velocity and has no AWF inventory. Tangential
diffusion and normal-stress gradients inside the locally parallel auxiliary
strip remain approximations. The equilibrium reference is supported on the
flat approach region; its mask must be zero on curved terrain. Matching-height
sensitivity is required when comparing these policies.

`neutral.terrain-wall.csv` records projected/pointwise normal errors, physical
mass defects, actual/physical-wall area, resolved/unresolved volume, N-sampling
inverse-map error, projection audits, positivity, chi and strong interface flux
residuals. `neutral.terrain-points-r*.csv` records local frames, roughness, mask,
requested/effective chi, all seven sensors and source/traction data.
`RANSStageAudit=1` adds independent assembled free-row residuals. Strong flux and
pointwise normal residuals are discretization errors; a small weak residual
alone does not establish pointwise enforcement. The paper-shaped hill is C1
with a curvature jump at each foot. Keep its global maximum errors visible,
and inspect global strong L2 and regular-region errors under refinement as well.
At a shared C1 foot the differentiated flux has two one-sided targets; a
continuous velocity trace cannot match both when their jump is nonzero. The
compact-hill suite reports that jump, its half-jump lower bound for the maximum
error, and the measured global maximum. It requires strong L2 convergence;
maximum-error convergence at these feet is a documented limitation. The smooth
periodic-curved suite separately requires decreasing global maximum and strong
L2 normal errors. This preserves the failed maximum evidence and does not
equate a small weak residual with pointwise enforcement.

AWF hill verification/pilots use dt=2e−5 s; the log-wall spatial verification
uses dt=2e−4 s. The larger AWF timestep produced a difficult CLS startup and
is retained as a failed development case. These choices do not establish a
general timestep stability limit or time-converged experimental predictions.

## Reproduce verification and the width pilots

Run suites sequentially below, with a fresh output directory for each. Each
suite permits two cases concurrently. The full-hex CLS hill uses eight MPI
ranks; all other cases use at most four. Supply at least eight available CPUs
to the pilot (sixteen for two concurrent cases). Every
subprocess has a timeout; retained receipts include commands, source/binary
hashes, mesh fingerprints, actual CPUs and thread/rank configuration.
Full-hex hill cases allow up to 7,200 seconds; their centered CLS Jacobian
probes can be costly. This runtime allowance does not relax the solve gates.

```bash
python3 "$RANS_TEST/validate_terrain.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/planes" --suite planes --jobs 2
python3 "$RANS_TEST/validate_terrain.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/smooth-curved" --suite periodic-curved --jobs 2
python3 "$RANS_TEST/validate_terrain.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/curved" --suite curved --jobs 2
python3 "$RANS_TEST/validate_terrain.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/transitions" --suite transitions --jobs 2
python3 "$RANS_TEST/validate_terrain.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/preservation" --suite preservation --jobs 1
python3 "$RANS_TEST/validate_terrain.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/pilot" --suite pilot --jobs 1
python3 "$RANS_TEST/validate_standard_awf.py" --build "$RANS_RUN/build" \
  --output "$RANS_RUN/legacy" --suite checks --jobs 1
python3 "$RANS_TEST/record_terrain.py" --build "$RANS_RUN/build" \
  --suite "$RANS_RUN/pilot/evidence.json" \
  --suite "$RANS_RUN/preservation/evidence.json" \
  --legacy "$RANS_RUN/legacy/validation-checks.json" \
  --ctest "$RANS_RUN/legacy/ctest.log" --output "$RANS_RUN/wg4.json"
python3 "$RANS_TEST/plot_terrain.py" \
  --wg2 "$RANS_SOURCE/docs/plans/hill-wall-geometry-evidence/wg2.json" \
  --wg4 "$RANS_RUN/wg4.json" --output "$RANS_RUN/figures"
```

The `preservation` suite advances five 2D periodic plane controls for five
flow-through clocks: standard/log wall and both fixed/switched corrected AWF
with standard and CLS closures. On four single-threaded ranks, 40 steps,
degree 4 and 12 normal layers,
their final pointwise relative errors are at most 4.53e−8 in U, 3.16e−8 in k
and 6.74e−7 in epsilon; |normal velocity|/u_tau is below 5.12e−8. This tests
stationarity of the new terrain adapter with resolved `HighRe` diffusion.
All four AWF controls retain chi=1 throughout the accepted steps; complete
saved fields are also compared between fixed and switched modes.
The inlet/outlet preservation and refinement studies of the original adapter
remain separately documented in the neutral-ABL example. Finite polynomial
degree does not represent a logarithm or reciprocal-height profile exactly.

`validate_standard_awf.py --suite checks` also runs the ten CTests, original/zero
weight and switch-off recovery, full-history restart, flat strip budgets and
independent kernel/short CLS regressions. `record_terrain.py` checks executable
hashes and the completed timestep/saved-field times before collecting evidence.

The pilot writes paired floor/hill receipts for four transition lengths, complete
field changes from width 1, and whole-height incident/crest/wake `profiles.csv`.
An optional `--suite pilot-2d --jobs 1 --publish-to "$RANS_RUN/pilot"` run in a
separate output directory can prepare the 2D cases while the full-hex case runs.
Allocate four additional available CPUs and keep the combined total within 16.
Only completed solver cases are published, using atomic cache symlinks; existing
destinations are untouched and original input/output paths remain valid. The
full pilot checks fixture/executable hashes before reuse and applies its
independent wall/continuity oracles to every case.

Before claiming an experiment reproduction, calibrate equivalent roughness and
displacement from the measured approach layer, match its finite thickness and
facility boundaries, converge mesh/timestep/H and placement/width independently,
run to stationarity, and average long enough. The supplied pilot uses a 40 mm
hill with 100 mm half-length, a 0.30 m roof, 2 mm sand-grain height, and a
0.032 m periodic span in 3D, with reference u_tau=0.3 m/s. These are not the
paper's 5 mm cube array,
0.25 m approach-layer thickness, or complete wind-tunnel geometry. URANS
comparison of turbulence intensities must combine modelled and resolved
covariances; it does not resolve all measured spectra or intermittency.
