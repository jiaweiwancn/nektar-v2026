# Terrain wall geometry and upstream AWF transition

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: implementation plan and deterministic validation
- Origin Date: 2026-10-31
- Verification Status: WG1 VERIFIED; WG2--WG4 in progress
- Version Label: hill_wall_geometry_v1

## Objective and scope

Extend the delivered standard k-epsilon / high-Re CLS URANS implementation
from a flat wall to a two-dimensional terrain graph and its spanwise extrusion.
Provide a prescribed upstream transition from equilibrium-preserving AWF to
original AWF, independently of the existing departure sensors. This is enabling
work for the Cao--Tamura (2006), doi:10.1016/j.jweia.2005.10.001, hill experiment;
it is not yet an experiment-calibrated hill reproduction or a claim that URANS
resolves measured broadband spectra/intermittency.

Baseline: `cba0df5d985ad863d438c48dd2e094769874e1b6`.
Implementation branch: `codex/hill-wall-geometry-transition`.
Keep the existing flat adapter as the default and retain its regression cases.
The compatible log wall remains a standard-k-epsilon policy; both AWF policies
support both closures. General non-extruded three-dimensional terrain, low-Re
CLS and a resolved cube array are outside this extension.

## Execution environment

Source changes and Git work use the Mac checkout. Builds and solver/numerical
validation use an isolated checkout/build under
`/media/jay/新加卷/nektar-v2026-runs/hill-wall-geometry-transition` on Linux.
The primary Linux checkout is on `codex/hit-inflow-generator` at
`6307bcec190adf0af0416096ca2b261f34fd90d3`; preserve it and the shared install.
Connect with `ssh jay@100.105.247.110` while `linux-usb` is unavailable.
Never exceed this project's **16 total workers**. Set OMP/BLAS threads to one,
bound compiler/test/MPI concurrency explicitly, and record commands, source
hashes, binary hashes, timeouts, exit status and failed development attempts.

## Geometric and numerical contract

Use the physical terrain curve r(s), unit tangent t and inward unit normal n.
With signed curvature c defined by dt/ds=c n, normal coordinates have metric
J(eta)=1-c eta. The wall strip has physical normal heights 0<P<H<N,
volume per unit physical-wall area V=H-c H^2/2 and interface area J(H).
Sample N along the physical normal, not the global vertical. Compute geometric
closest-wall distance for CLS separately from aerodynamic roughness/displacement.
Reject folded/ambiguous offset geometry and strips outside a stated thin-layer
curvature range. Support a plane and the paper's compact cosine-squared hill.

The original AWF remains a locally parallel auxiliary model. Extend its strip
inventory, source integration and interface exchange with the actual metric;
keep the curvature contribution of the auxiliary normal shear. Tangential
diffusion/normal-stress gradients inside the unresolved layer remain modelling
approximations and must be documented, not described as an exact integration of
the full curved RANS equations. The equilibrium correction is intended for the
flat approach region; the spatial gate must be zero over the hill and wake.

Terrain velocity boundaries require a vector tangential condition and a normal
constraint. Use an augmented boundary multiplier with a lower-order trace space
for the normal constraint, converged with the coupled implicit stage; the
augmentation is a preconditioning parameter, not a permitted finite-penetration
wall law. Retain Cartesian velocity unknowns and verify normal constraint,
boundary traction and refinement independently. The terrain path initially
requires the current-coefficient matrix-free linear backend.

For a spatial ramp with endpoints xa<xb, define r=clip((x-xa)/(xb-xa),0,1)
and w=1-10r^3+15r^4-6r^5. Fixed mode uses chi=w*chi_fixed; switched mode uses
chi=w*chi_flow. The same effective chi selects mean dissipation and both P/H
epsilon values. Traction and production retain their original-AWF definitions.
Zero w must bypass corrected-reference admissibility checks and exactly select
original AWF. Ramp coordinates refer to the physical wall foot. Optional
stationary wall-property expressions provide spatial roughness and log-wall z0.
All new physics/configuration must enter restart compatibility checks.

## Milestones and acceptance gates

| Milestone | Deliverable | Required evidence | Status |
|---|---|---|---|
| WG1 | Geometry, metric/source integration and spatial-ramp primitives | Linux kernel tests: closest points/normals, rotations, curvature/Jacobian, invalid geometry, endpoint/derivative limits and independent quadrature; existing kernel tests pass | CLOSED |
| WG2 | Terrain solver adapter, CLS distance and compatible standard log wall | Actual Linux solves: flat reduction, inclined-plane covariance, curved-wall refinement/normal constraint, source/interface budgets and both closure/wall families; positivity and coupled residual checks | PLANNED |
| WG3 | Spatial transition/property integration and reproducible workflows | Prescribed-mask receipts; inactive-correction/original recovery; spatially varying properties; 1/2/4-rank agreement including remote samples; history restart and mismatch rejection; empty-floor transition controls | PLANNED |
| WG4 | Hill-enabling pilot, regression and delivery | Paper-shaped terrain pilot for standard log wall and original/corrected AWF with both closures; four upstream ramp lengths with paired flat controls; relevant legacy regression; documented errors, limits and reproducible commands | PLANNED |

Commit and push after each completed milestone. A gate is closed only when its
evidence exists; failed attempts and any scientific limitations stay visible.
Initial tolerances: kernel relative errors <=5e-10; converged stage and projected
normal constraint <=1e-7 in the verification cases; serial/MPI/restart field
differences <=1e-6 relative to max(1, field magnitude). Curved-wall pointwise
normal error must decrease under refinement; weak enforcement alone is not
evidence of pointwise impermeability. Record measured orders rather than assume
them. Pilot runs are numerical functionality checks, not accuracy acceptance
against experimental data. Threshold changes require an evidence-based reason
in this record.

## Transition-size study design

Use x=0 at the hill crest, H_h=0.04 m and upstream foot x=-0.10 m. Keep the
downstream ramp endpoint xb fixed strictly upstream of that foot; vary width
Delta/H_h in {1,2,4,8} without moving xb. Use a C2 ramp and sample it densely
enough to distinguish a physical-width change from mesh resolution. For every
width run a corresponding empty-floor control. Record incident U/k/epsilon,
skin friction, effective chi and fallback reasons at the hill foot, in addition
to hill/wake profiles and separation indicators. A separate placement test and
wall-strip-height / mesh / timestep studies are required for a later calibrated
experiment campaign.

## Monitoring and outputs

Each command has a hard timeout and retained stdout/stderr/exit receipt. Inspect
progress while processes are alive; retain a failed attempt before fixing its
cause and rerunning a distinct case. Store all generated fields, meshes, history
and detailed logs on the Linux large volume. Keep only small reproducible case
generators, regression fixtures and compact evidence in Git. Final delivery
includes public preparation commands, configuration/restart documentation and a
report separating numerical verification from experimental validation.

## WG1 closure (2026-10-31)

Implemented physical plane/compact-hill geometry, whole-hill closest-point
search, signed wall distance, orthonormal frames, normal-ray/curvature guards,
exact metric volume/interface measures, metric-weighted AWF source integration
and a bounded C2 spatial mask. Zero mask recovers all original wall outputs,
including when corrected-reference high-Re checks would otherwise reject.

Linux: three CTests pass, including the two existing kernel/strip tests.
The geometry probe performs 3,759 checks and five invalid-input rejections.
Independent SciPy quadrature checks 54 smooth/rough, standard/CLS,
negative/zero/positive-curvature and partially/fully corrected combinations;
maximum relative source/shear integration error is 5.76e-16.
Compact receipt: `hill-wall-geometry-evidence/wg1.json`.

Tested source is plan commit `0fb850f` plus the WG1 changes. Geometry .cpp SHA256:
`c081dd17008e475394ab9bd010f6b94b26c513b523e211f739861706095de67c`;
probe .cpp SHA256:
`9576db509a3f16e264cf89309d7f1f11fa5e3c7fa0bc8bc4529f0bf756414dbd`.
Full evidence is in the Linux run root's `wg1` directory. Configure/build/test:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
RANS_RUN=/media/jay/新加卷/nektar-v2026-runs/hill-wall-geometry-transition
export TMPDIR="$RANS_RUN/tmp"
cmake -S "$RANS_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$RANS_RUN/build" -DCMAKE_BUILD_TYPE=Release -DRANS_BUILD_SOLVER=OFF \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$RANS_RUN/build" --parallel 8
ctest --test-dir "$RANS_RUN/build" --output-on-failure -j 1
python3 "$RANS_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_wall_geometry.py" \
  --probe "$RANS_RUN/build/RANSWallGeometryProbe" --output "$RANS_RUN/wg1/quadrature"
```

Retained development attempts: two CMake discovery failures before using the
documented package directory, and one C++ mixed-auto-declaration compile error,
corrected before the successful run. These supplied no numerical results.
