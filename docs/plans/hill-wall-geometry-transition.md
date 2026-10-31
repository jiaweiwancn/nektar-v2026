# Terrain wall geometry and upstream AWF transition

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: implementation plan and deterministic validation
- Origin Date: 2026-10-31
- Verification Status: WG1--WG2 VERIFIED; WG3--WG4 pending
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
constraint. Use an augmented boundary multiplier with the assembled continuous
free velocity trace space for the normal constraint, converged with the coupled implicit stage; the
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
| WG2 | Terrain solver adapter, CLS distance and compatible standard log wall | Actual Linux solves: flat reduction, inclined-plane covariance, curved-wall refinement/normal constraint, source/interface budgets and both closure/wall families; positivity and coupled residual checks | CLOSED |
| WG3 | Spatial transition/property integration and reproducible workflows | Prescribed-mask receipts; inactive-correction/original recovery; spatially varying properties; 1/2/4-rank agreement including remote samples; history restart and mismatch rejection; empty-floor transition controls | PLANNED |
| WG4 | Hill-enabling pilot, regression and delivery | Paper-shaped terrain pilot for standard log wall and original/corrected AWF with both closures; four upstream ramp lengths with paired flat controls; relevant legacy regression; documented errors, limits and reproducible commands | PLANNED |

Commit and push after each completed milestone. A gate is closed only when its
evidence exists; failed attempts and any scientific limitations stay visible.
Initial tolerances: kernel relative errors <=5e-10; converged stage and projected
normal constraint <=1e-7 in the verification cases; serial/MPI/restart field
differences <=1e-6 relative to max(1, field magnitude). Smooth curved-wall
pointwise normal error must decrease under refinement; weak enforcement alone
is not evidence of pointwise satisfaction. For the compact C1 hill, require
decreasing global strong L2 normal error and record its global maximum and both
one-sided foot targets, as justified below. Record measured orders rather than
assume them. Pilot runs are numerical functionality checks, not accuracy acceptance
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

## WG2 development and verification decisions

The terrain adapter now uses physical normal rays, metric strip measures and
weighted sources, a Cartesian vector Robin operator, and an assembled continuous
velocity-trace multiplier. Its conservative constraint and adjoint use the
orthogonal trace projection as defined below, including constant-test
conservation at prescribed-velocity corners. Measure the unprojected normal
error and the actual integrated normal/mass defect independently. The original
auxiliary AWF remains locally parallel, as specified above.

Several failed development attempts are retained under the Linux run root's
`wg2` directory:

- Separately lagging Cartesian Robin cross terms failed rotated AWF cases.
  A simultaneous SPD velocity-block solve fixes this coordinate dependence.
- A pressure-order normal trace left high velocity-normal modes uncontrolled.
  Two-mesh normal errors did not decrease, and a 3D inclined case lost positivity.
  A facewise full-order trace instead imposed conflicting endpoint constraints.
  The assembled continuous free velocity trace removes both problems.
- Native curved-element inverse coordinates needed Newton polishing for exact
  physical N samples. XML partitioning also rounded curved edges/faces to eight
  decimal digits, changing MPI normals by up to 3.5e-6. Curved geometry now uses
  the same full-precision serialization as vertices. An isolated source
  SpatialDomains library verifies the fix without changing the shared install.
- The new surface-integral diagnostic initially reduced an already collective
  ExpList integral again. Local element integration followed by one collective
  fixes the MPI volume accounting, including ranks without a wall.
- A zero initial normal reaction perturbed the CLS normal stress and could send
  the nonlinear iteration to a different branch. Initialize the reaction with
  minus the current resolved normal stress. This is an iteration warm start,
  not a reset to an analytical profile.
- Derivative-based transpiration made the coupled iteration difficult.
  Multisecant acceleration and a surface filter did not resolve its refinement
  failure. Neither experimental iteration option is delivered. The retained
  terrain stage uses a positivity-constrained Newton--Krylov iteration with
  native GMRES. A fixed active CLS invariant branch defines each generalized
  Jacobian; physical residuals retain max(strain,rotation). Kernel checks verify
  active-branch agreement and directional-derivative additivity. Redundant
  pressure/reaction values are projected to their native CG spaces to avoid
  spurious off-space null modes.
- Cache a Cholesky factor of the exact trace mass for up to 1,024 trace modes;
  retain the iterative projection for larger systems or `RANSNormalDirectLimit=0`.
  Runtime idempotence/self-adjoint projection audits verify the same operator.

Three refinement levels are required. Earlier development fixtures refined
the streamwise mesh and H--N layers but kept the outer pressure mesh fixed;
these are not a whole-domain h-convergence test. They exhibited nonmonotonic
strong flux/normal errors and are retained as failed development evidence.
The accepted study must refine every resolved normal layer as well as every
streamwise column. Keep the C1 hill-foot curvature jumps and global maxima
visible, together with global L2/regular-region errors; do not replace a failed
maximum by a weak norm.

The original-CLS plane test at dt=0.002 s met the 1e-6 primary-field comparison
but its derivative-sensitive outputs differed by up to 1.39e-5. Tightening the
nonlinear target to 1e-10, with linear targets 1e-14 and 1e-15, did not converge
in 4,000 iterations (last weak divergence about 7.8e-9). Resolve its rapid
unbalanced original-AWF dissipation transient using the hill pilot timestep
dt=0.0002 s: the observed complete-field rotation difference is then at most
2.34e-9. **The 1e-6 complete-field gate is retained.** Other plane tests retain
dt=0.002 s. The two-step tests exercise IMEX2 startup and continuation.

Curved verification targets nonlinear/pressure iteration at 1e-8 and independently
requires assembled free-row residuals, weak divergence, projected normal error
and physical mass defects below 1e-7. This uses the originally declared gates;
the earlier 1e-9 iteration target was unnecessarily tighter and caused prolonged
stagnation. All failed/overly tight attempts remain available. The final accepted
case receipts and measured refinement orders will be recorded at WG2 closure.

The volume-gradient reconstruction of curved strip continuity failed an early
maximum-normal refinement gate (finest error 0.0320 m/s). An impermeable variant
was investigated but is **not the requested model and is not delivered**.
Integrated strip continuity remains required. Differentiate the pulled-back
scalar strip flux in the native surface space, rather than applying a continuous
chain rule to separately interpolated volume gradients and rotating frames.
For streamwise/normal/span coordinates, retain
J(H) W_H + H d_s(U_P) + V(H) d_span(U_span,P) = 0,
and the corresponding expression at P. Constant scalar flux has exactly zero
surface derivative even on a curved mesh; face contributions telescope under
integration. The physical/weak normal errors and strip side/interface flux
budget remain independent acceptance gates.

A compatible synthetic initial streamfunction is G(normal distance) through N,
blended to the roof-compatible G(e) by 2N. Thus its H/N tangential log traces
are initially independent of wall foot; zero initial transpiration follows
from continuity, without constraining subsequent normal velocity. Numerical
verification is rerun with this initial condition and conservative reconstruction.

Prescribed inlet/wall corners require a conservative dual trace. Let P denote
the orthogonal L2 projector onto the native free scalar velocity trace V_f,
u=P1, a=1-u and c=(1,u)=(u,u)>0. Define

```
C f = P f + u (a,f)/c,          C* v = v + a (u,v)/c, v in V_f.
```

C is an idempotent map into V_f and preserves the integral of every f; its
adjoint lifts the test space to include the constant. Enforce
C(U_normal-W_H)=0. Use C*lambda for the normal reaction and beta C*C for
augmentation, with the same beta C*C W_H on the right-hand side. This retains
an SPD velocity block and makes augmentation vanish at the converged constraint.
In particular, integrated interface normal flux balances the strip side flux
also at prescribed-velocity corners. Runtime audits check P's orthogonality,
C's idempotence/integral preservation, the adjoint identity and symmetry of
C*C. Independent GLL/endpoint quadrature checks the H/P reconstruction and
actual interface-plus-side-flux defect from saved accepted wall data.

The periodic Newton state also needs an explicit zero-mean pressure gauge;
projection into P-2 alone leaves its constant null mode. Outflow cases retain
the pressure level determined by total traction. Difference probes are bounded
in physical field scales, not solely by a volume L2 norm that can hide a large
fine-grid vertex perturbation. Full native GMRES storage is bounded by
RANSNewtonKrylovDimension (default 240, allowed 20--512); difference-step and
Krylov-dimension controls enter restart compatibility. Failed pre-gauge and
pre-conservative-corner solves, including 3D CLS line-search failures, are
retained in the Linux development directories. They are not accepted results.

The original unconditional maximum-normal refinement gate needs a regularity
qualification. At a compact C1 hill foot, the two one-sided differentiated strip
fluxes give different targets, while the resolved CG velocity has one continuous
trace. If the target jump is delta, every continuous trace has maximum error at
least |delta|/2 there. The saved corner audit reports both targets, this lower
bound, the actual trace jump and the measured maximum at every refinement;
the failed global maximum is retained, not called a convergence pass. Add a
smooth periodic cosine hill, with matching one-sided curvature across the
periodic seam, to require decreasing global maximum **and** strong L2 normal
errors. Preliminary AWF maximum errors on this smooth problem decrease
2.29e-3, 1.11e-3, 2.38e-4 m/s, whereas the compact-foot maximum grows. The C1
case must still pass decreasing global strong L2 normal and interface-flux
errors, coupled/free-row tolerances and independent physical flux budgets.
This is an explicit geometry-regularity qualification, not a relaxed nonlinear
or conservation tolerance. Maximum-norm convergence at the compact feet is a
recorded limitation of this continuous-trace discretization.

The AWF startup verification uses dt=2e-5 s after the dt=2e-4 CLS case failed
to converge; all three AWF meshes use the same timestep. The compatible log
wall retains dt=2e-4 s, resolving its rapid wall-flux transient more effectively.
Different wall policies are not compared as an equal-time convergence pair.
The source guard for terrain-wall presence is global: partitions with no wall
faces are valid participants in pressure, trace projection and N sampling.

## WG2 closure (2026-10-31)

Closed with integrated AWF strip continuity retained. All 47 Linux cases pass:
21 flat/rotated/legacy/full-3D cases, 13 compact C1-hill cases and 13 smooth
periodic curved-wall cases. All ten CTests pass. The 32-file tested source
manifest matches the Mac checkout; source baseline is `cfb870a` plus those
manifest changes. Compact evidence: `hill-wall-geometry-evidence/wg2.json`.

Maximum weak divergence is 6.53e-9, conservative normal constraint 7.46e-9,
assembled free-row residual 2.68e-9, interface normal-integral defect 1.02e-11
and resolved mass defect 2.19e-10. Complete-field rotation/flat-reduction
agreement is at most 2.71e-7; augmentation, trace-backend and CLS difference-step
comparisons are below 7.0e-9. Independent GLL/strip-side quadrature verifies
nonzero AWF transpiration and actual interface-plus-side-flux conservation.

| Normal-condition error (m/s) | h | h/2 | h/4 | Last measured order |
|---|---:|---:|---:|---:|
| Smooth periodic standard AWF, global maximum | 2.295e-3 | 1.114e-3 | 2.375e-4 | 2.23 |
| Compact standard log wall, global maximum | 2.594e-4 | 1.126e-4 | 3.903e-5 | 1.53 |
| Compact standard AWF, global strong L2 | 1.116e-4 | 9.739e-5 | 7.503e-5 | 0.38 |
| Compact standard AWF, global maximum (limitation) | 2.303e-3 | 3.232e-3 | 4.804e-3 | -0.57 |

The C1 AWF maximum is **not converged**. Corner receipts show the one-sided
target jump and its half-jump lower bound for every continuous velocity trace.
The regularity-qualified gates and their rationale above are part of this
closure. Both tangential and k strong flux errors decrease on both geometries.
This establishes the stated numerical adapter contracts, not a mesh/time-
converged prediction of the experiment.

The final commands use `validate_terrain.py --build BUILD --suite SUITE
--output RUN --jobs JOBS`: `planes` at `wg2/wallless-planes` (one case at a
time), `curved` at `wg2/final-c1` (two), and `periodic-curved` at
`wg2/final-smooth` (one). The evidence records every exact command, CPU lease,
MPI rank/thread setting and binary/fixture hash. Four compiler workers and at
most sixteen total project workers were used; all OMP/BLAS threads were one.
The isolated geometry library was selected through `LD_LIBRARY_PATH`.
The last smooth cases were moved off CPU64, an SMT sibling used by the separate
HIT study; unfinished logs are retained in `wg2/final-smooth-reallocated` and
accepted earlier cases retain their actual allocations. The primary Linux
checkout and shared installation were left untouched.
