**M3 closure for the revised CLS/EP-AWF native flat wall stage**

M3 is closed on 2 November 2026 on codex/cls-epawf-revised. The tested
Linux code revision is 84442541ad3c4ff014b251f314c22b08ae84cbfe. The complete external
production-solver build and all 22 registered tests pass. The independent
M3 qualification completes 8,085 checks, including actual full and
condensed nonlinear stages. M4 and subsequent milestones remain open.

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: implementation and native computational verification.
- Origin Date: 2026-11-02.
- Verification Status: VERIFIED for the M3 native serial flat stage, its
  spatial residual, coupling and Newton condensation. Original constitutive
  stability and general inner-k positivity remain FAILED AS SPECIFIED.
- Version Label: cls-epawf-revised-m3-v1.

The [revised note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex)
and [frozen contract](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-contract-v1.json)
retain SHA-256 d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f
and 4e7bd75022409b8dc2ef31cc73425779b646a97535760b31eeea4ec49087de8a.
The original CLS and scalar AWF kernels, coefficients and prior solver
selections are unchanged. The contract's historical availability entries
remain frozen; this report records delivery of its nine M3 verification IDs.

**Implemented stage and field layout**

[RANSEquilibriumWallSystem](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallSystem.cpp)
is registered in the complete IncNavierStokesSolver executable and its
SolverType enumeration. Its [runnable guide](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Examples/RevisedEquilibriumWall/README.md)
and [small session](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS/EquilibriumWallFlat.xml)
document the opt-in selection, functions, physical parameters and output.
The native compatibility gate is qualified against actual Nektar++ 5.9.0;
the note's 5.10.0 interface discussion is not assumed to describe this build.

[RANSWallFieldLayout](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp)
creates real subdomain ExpList/ContField/AssemblyMapCG objects. A full native
CG injection represents the common inner/outer velocity traces; k has one
full-domain CG space. Epsilon has only outer elements and PDE rows. Inner
and outer pressure spaces remain separate, including their relative
constant mode. One volume-mean constraint fixes pi=p_mechanical/rho+2k/3
across both regions; no independent pressure Dirichlet values are imposed.
The two-element p=4/p=2 fixture has 179 unknowns: 45 each for U,W,k, 25 for
outer epsilon, 9 each for the two pressures and one gauge multiplier.

Native physical basis values and derivatives use LocalRegions APIs and
native assembly maps/signs. Physical element bounds come from vertices,
since the framework's padded search bounding box is unsuitable for domain
classification. Independent physical Gauss quadrature is split at native
inner cell boundaries as well as wall-profile breakpoints. This adds only
optional physicalBreakpoints to the M2 profile quadrature, with no source-law
change; all prior M2 native profile results remain identical.

[RANSWallPatch](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp)
assembles conservative momentum, full-domain k, outer epsilon and both
continuity blocks. Every differential field retains alpha*(state-history),
including free inner U,W,k. Inner stress is -pi I+D_A S; its actual production
is the signed contraction -R_A:G, including the isotropic contribution when
discrete divergence is nonzero. The outer operator uses the complete
original three-by-three CLS stress, production and Yap source, with physical
2D mean axes (x,z) and no added outer molecular diffusion.

The inner force is b+(1-a_h)*F_A. Neither mass nor the residual is multiplied
by a_h. Each current K is evaluated at the physical r_K, and the ray-constant
F_A is the one-sided outer value of d_x(U^2)+d_x(pi-2k/3)-b_x at H. Profiles,
cutoffs, traces and sampler dependencies are recomputed for changed states.
The source uses P_A-(1-chi)E_A-chi*p_A0 once, with current sampled K in p_A0;
forced production remains distinct from that reference production.

Ground velocity is strongly zero, and ground modeled k flux is naturally
zero. Shared native velocity/k traces and one common flux give opposite
internal dual loads. Epsilon at H is projected from the dependent E_EP(H);
its flux is determined by the outer equation. External side/top data must
agree with the dependent epsilon trace at corners. Output consists of eight
native field files with their own element sets and a complete stage/history
JSON; no full-domain epsilon mask is used.

**Independent operator and pressure verification**

[RANSEquilibriumPatchProbe](/Users/jay/Documents/nektar-v2026/library/UnitTests/SolverUtils/RANSEquilibriumPatchProbe.cpp)
exports native maps, basis/mesh metadata, residuals and matrices.
[validate_cls_epawf_m3.py](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m3.py)
reconstructs Jacobi bases independently and uses dimensional indexed CLS,
separate physical quadrature and manufactured strong derivatives. It does
not call the C++ stress/source operator to construct expected residuals.

All seven native field roles pass physical polynomial value/derivative,
CG dual-assembly and nonpolynomial L2 projection moment checks. The largest
independent six-block residual error across chi=0,0.37,1 is 4.170e-16,
using the componentwise scale max(1,max(abs(expected load))).
Those arbitrary states include full strain, nonzero discrete divergence,
normal transport, body forcing, separate pressure traces and all storage
terms. The pressure/body-potential sampler cancellation error is
6.661e-16.

The mixed divergence matrix has rank 17 of
18 pressure rows, giving exactly one common
constant null mode. The separate interface constraint matrix has full
rank 5 of 5 rows.
The globally gauged 61-row pressure/interface
system is nonsingular. Its prescribed-reference-stress probe produces
zero velocity to 1.728e-14, global mean
1.683e-19, and the measured matching pressure jump
-0.074874508582116653=-11*K0/75. The normal traction mismatch is
6.939e-17. The jump is determined by the retained modes
and full traction, rather than prescribed as a pressure boundary value.
This isolates the reference stress/pressure example; the complete neutral
ABL reference and timestep are M4 work.

The compatible top-reference check verifies nu_t U'=u_star^2, k'=0,
epsilon'=-u_star^3/(kappa_c Z^2), and the gradient flux
-u_star^4/(sigma_epsilon Z), including native Dirichlet trace values.
The physical outward diffusive flux has the opposite sign. A nonzero
traction placed at the strongly fixed ground gives exactly zero load on
free velocity tests, reproducing the note's strong-no-slip comparison.

| Frozen verification ID | Delivered evidence |
| --- | --- |
| M3.outer_normal_pressure | Independent normal-momentum residual, full stress/pressure convention, mixed pressure rank and normal traction |
| M3.compatible_outer_boundary | Native compatible top values and analytic stress/epsilon gradient-flux factors |
| M3.physical_domains_and_sampler | Seven native field roles, actual subdomain maps, physical r_K, one-sided H+ forcing, invalid-H rejection |
| M3.strong_no_slip_comparison | Ground strong constraints and exactly zero ground traction load on free tests |
| M3.flat_analytical_momentum | Native restricted CG weak solve, independent FE assembly, zero/nonzero forcing, h/p refinement, reaction tractions and energy |
| M3.developing_flat_mass_pressure_forcing | Nonzero shared normal flow, regional mass budgets, full normal transport and physical pressure/body cancellation |
| M3.retained_inner_k | Full-domain k modes, separate measured inner storage, full-strain production, conservative transport and zero ground k flux |
| M3.dependent_dissipation | Outer-only epsilon rows, current-K dependent H trace, compatible corners and independent interface checks |
| M3.reference_normal_traction_example | Native pressure/interface solve reproducing -11*K0/75 with one gauge |

**Steady analytical AWF reduction**

For frozen K and supplied matching velocity, the native inner 2D CG space
is restricted to constant tangential dependence and solves (D_A U')'=a_h F.
An independent one-dimensional Jacobi FE assembly verifies every retained
coefficient, weak reaction, physical derivative and the energy identity.
This is the specified steady scalar reduction, not a coupled neutral-ABL
transport solution. K, geometry and boundary data are fixed during each
refinement study. All three forcing values are retained.

| p | Inner elements | RMS velocity error, F=-0.02 | F=0 | F=0.03 | Maximum weak traction error |
| --- | --- | --- | --- | --- | --- |
| 4 | 19 | 1.496e-06 | 1.430e-06 | 1.331e-06 | 5.595e-10 |
| 6 | 19 | 7.247e-09 | 6.928e-09 | 6.449e-09 | 2.920e-14 |
| 8 | 19 | 3.899e-11 | 3.727e-11 | 3.470e-11 | 3.830e-15 |

RMS velocity errors are dimensional sqrt(integral(error^2)/H), in m/s.
The fixed graded partition includes h and resolves the molecular/virtual
length. Weak traction means the Galerkin endpoint reaction; pointwise
D_A U'_h is also retained in the JSON and converges separately. The
coarse two-element p=4 reduction is deliberately unresolved and has a
large profile error; it is retained as a refinement control.

The p=4 h study uses a fixed graded coarse grid and bisects each interval:

| Inner elements | RMS velocity error, F=0 |
| --- | --- |
| 8 | 1.710e-04 |
| 16 | 1.088e-05 |
| 32 | 4.826e-07 |

The measured two refinement rates are
3.9747, 4.4947.
The weak energy identity agrees to 1.905e-14
after the declared scaling. All three forcing cases and pointwise traction
errors are available in the compact verification data.

**Evolving transport, storage and condensation**

The independent [manufactured-case generator](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS/cls_epawf_m3_cases.py)
provides evolving, divergence-free velocity, nonconstant transported k,
outer epsilon, normal flow across H and separate pressure functions. Its
physical body/scalar sources come from symbolic strong derivatives.
Inner k has degree six; velocity and pressure have lower degrees. The
three-element mesh has a first outer cell of thickness 0.02 and samples
both active and inactive Yap branches. Every resolved audit checks that
both branches are actually sampled.

History is the supplied exact state minus its time derivative divided by
alpha. This verifies the complete stage storage operator at several times
and alpha values without asserting that alpha is a native IMEX time map.
The inner storage norms at t=0.3 are nonzero for U,W,k. A separate regional
mass integration verifies the shared nonzero W_H and the cancelling inner
side/normal fluxes. Its maximum net regional boundary mass flux is
1.752e-16. No residual anchor, solution reset, source floor or model
coefficient change is used; all manufactured forcing is explicit session data.

Residuals use sqrt(r_free^T M_free^-1 r_free) divided by
(field_scale/time_scale)*sqrt(volume). Continuity uses u_star/H; state
errors use field_scale*sqrt(volume). The scales are u_star=0.4, K0,
epsilon0=u_star^3/(kappa_c H), pressure u_star^2 and time H/u_star.
Native JSON retains raw block norms, the global gauge and removed-row
defects. The oracle also checks scaled state/storage errors and complete
interface Gauss-rule maxima, which bound the corresponding positive-weight, unit-area
physical surface quadrature norms. The frozen acceptance threshold is 1e-9.

| p | Evaluation time | alpha | Maximum scaled residual | Maximum scaled field L2 error |
| --- | --- | --- | --- | --- |
| 4 | 0.3 | 5 | 4.074e-01 | 1.308e-05 |
| 6 | 0.3 | 5 | 9.179e-12 | 2.193e-16 |
| 8 | 0.3 | 5 | 4.190e-11 | 2.314e-16 |
| 6 | 0 | 5 | 7.312e-12 | 2.125e-16 |
| 6 | 0.9 | 5 | 1.760e-11 | 2.348e-16 |
| 6 | 0.3 | 2 | 9.179e-12 | 2.193e-16 |

The unresolved p=4 k defect is visible. Degrees six and eight reproduce
the manufactured fields and satisfy the residual threshold; their remaining
roundoff floor is reported rather than treated as a convergence order.
The maximum resolved scaled storage-rate error is 3.011e-15;
the largest scaled resolved interface/ground check is 5.522e-13.

Actual nonlinear stage solves start from perturbed velocity fields:

| Stage | Newton iterations | Maximum scaled physical residual |
| --- | --- | --- |
| stationary-full | 3 | 2.055e-14 |
| stationary-condensed | 3 | 2.068e-14 |
| evolving-full | 3 | 6.538e-14 |
| evolving-condensed | 3 | 6.506e-14 |

The full and condensed evolving maximum coefficient difference is
1.939e-15 with comparison scale max(1,max(abs(coefficients))), in the same pressure
gauge. Each solution also passes the separate physical field-mass checks. The stationary comparison is
4.362e-16.

Condensation is a Newton Schur reduction of free inner velocity only.
The isolated 179-row stage eliminates
18 velocity rows and retains
161 rows, including every pressure mode, full-domain k
and outer epsilon. Inner state/history and storage remain in the recovered
stage. NumPy independently reconstructs the Schur matrix, residual, full
step and recovered step. The recovered-step discrepancy is
6.749e-15. True physical residual
differences along the implicitly recovered direction give errors
6.194e-09, 5.523e-10, 5.781e-11
for steps 0.001, 0.0003, 0.0001.
This separates finite-difference truncation from the tighter algebraic
Schur comparisons; it is not an analytic-Jacobian claim.

**Retained development findings**

[attempts.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m3/attempts.json)
retains all failed and successful attempts and exploratory diagnostics.
Early failures exposed native array/type aliases, padded physical bounds,
compiler warnings and missing SolverType registration. They are preserved
with their source commits and command/log hashes.

The 5.9 expression interpreter mishandles max(expression, literal-zero)
in the manufactured epsilon source. A small residual alone then accepted
the wrong manufactured target. Independent exact-field/source readback
caught it. The fixture now uses the equivalent unary positive-part form
(a+fabs(a))/2; a separate native control reproduces the original expression
failure. Shared interpreter and turbulence kernels are unchanged.

The first coarse evolving fixture also prescribed incompatible epsilon
corner data after projection of degree-six k onto p=4. Its k side traces
now remain exactly K0 at every p; incompatible data are still rejected.
Increasing quadrature did not resolve the subsequent thin-cell roundoff
floor. Native constant-centred value/gradient evaluation and physical L2
defect correction remove background cancellation while retaining the same
projection equations, verified independently for every basis moment. The
original thin mesh remains in the final qualification. Endpoint evaluation
also exposed and corrected an inner-output sampler applying the H-minus
side at the ground. These are numerical/API corrections, not model remedies.

**Linux execution and provenance**

The closing run is
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m3-attempt-017.
The preceding direct qualification and complete 22-test run at
95f203efd33ae5ab3d74fad402b79dcea2d0530f are retained as attempt 016.
The closing revision strengthens only acceptance checks/documentation;
compiled source, fixtures, solver binary and native probe results are
verified identical to that successful direct qualification.

The exact closing commands are in the
[execution receipt](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m3/execution-receipt.json):

```bash
cmake -S "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS" \
  -B "$CLS_EP_RUN/build" -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
  -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
cmake --build "$CLS_EP_RUN/build" --parallel 16
taskset -c 0 ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure
```

Tests use one rank and one numerical thread, with OMP_NUM_THREADS,
OPENBLAS_NUM_THREADS, MKL_NUM_THREADS and BLIS_NUM_THREADS all set to 1.
The registered M3 oracle executes the complete qualification once in the
closing suite; a duplicate direct solve is not needed. Compile allocation
is 16 workers and tests run sequentially, below the 72-worker cap.
The complete suite takes 708.93 seconds in this run.

[integration.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m3/integration.json)
binds 25 source/input pins,
12 executables, the static kernel archive and
37 resolved shared libraries. The resolved
shared-library set is unchanged from M2. Compiler flags include strict
warnings with Werror for the new native components. The isolated geometry
library is selected explicitly through LD_LIBRARY_PATH and checked with
ldd, resolving the external build's retained RPATH-cycle warning.
The M1 kernel/constitutive and M2 native-profile outputs remain identical
to their earlier accepted reports. The separate native CTest and oracle
reports agree across 70,188 numeric values to
1.765e-14 scaled error
(tolerance 1e-12), with exact matching structure. Their small roundoff
difference is retained with both file hashes; bytewise equality is not
claimed for that comparison. The earlier direct-oracle native report is
byte-identical to the closing oracle's native report.

The Linux HIT checkout remains clean at
f4b70baec10016dd098005bc5a5507b87e4a57bf. The shared build/installation is
preserved. Both source checkouts and final source hashes are checked;
user-owned untracked literature and unrelated Mac files are preserved.
The initial Tailscale route was used for development; a later interruption
was diagnosed read-only and the restored USB route completed qualification
and collection, without changing network settings.

Raw fields, matrices, sessions, logs, rejected controls and the complete
artifact hash index remain on Linux, outside the source checkout. Compact
[verification.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m3/verification.json),
receipts and source-bound manifests are committed with this report.

**Remaining qualification**

This native backend is serial, affine, axis-aligned and two-dimensional,
with a dense mixed-stage limit of 3,000 unknowns. The larger steady
restricted helper skips mixed-stage preparation and has a separate limit
of 2,000 radial unknowns. Passing these fixtures establishes M3's native
flat spatial stage and condensation, not scalability or time integration.
M4 owns the constrained neutral reference, native IMEX stages/startup,
repeated steps, full-history restart, three-dimensional fields and MPI.
M5 owns curved geometry. Original adverse momentum diffusion and failure
of general inner-k positivity remain unchanged failed physical gates;
robust long-time terrain qualification remains conditional on M8.
