**M0 contract and baseline for the revised CLS and EP-AWF project**

Status: contract v1 frozen; closure audit pending. The implementation branch
is codex/cls-epawf-revised, based on
6415398d856323c3116e646eaf7c882c66617602. Branch setup and the exact-base
Linux regression run were completed during M1. M0 adds the complete
machine-readable contract, equation mapping and baseline provenance audit.

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: implementation contract and reproducibility audit.
- Origin Date: 2026-11-02.
- Verification Status: contract frozen; native baseline already passed;
  M0 closure audit pending.
- Version Label: cls-epawf-revised-m0-v1.

The normative artifact is
[cls-epawf-revised-contract-v1.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-contract-v1.json).
Its executable checker is
[validate_cls_epawf_m0.py](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m0.py).
This pair specifies and checks the implementation contract. It does not
provide the coupled wall-patch solver reserved for M3.

**Source identity and authority**

[source-intake.json](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m0/source-intake.json)
records the Mac workspace, exact reference tips and both note hashes:

| Source | Frozen identity |
| --- | --- |
| Governing revised note | d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f |
| Comparison original note | 8b7f67035e553318a2b0896c9e31dd74d94bf0cc4864013ea209a48eddf98819 |
| codex/neutral-abl-rans-awf | cba0df5d985ad863d438c48dd2e094769874e1b6 |
| codex/hill-wall-geometry-transition | c48d9fb14afeca6caa80cf5b3f8ce0a34a02d2d7 |
| codex/cao-tamura-calibrated-validation | 6415398d856323c3116e646eaf7c882c66617602 |

The revised note is tracked and unchanged. The original comparison note
remains an unchanged local reference; its appended remedies are outside
the selected model. The three reference commits form an ancestor chain.
The existing scalar kernels and old CV/CG strips are references, while
the revised multiregion wall element remains new implementation work.

**Frozen model and field contract**

Use original CLSHighRe with the complete cubic tensor, original C(q),
max strain/rotation invariant, transport constants and algebraic Yap term.
Tensor components use physical x/y/z and G_ij=dU_i/dx_j; S and Omega have
no half factor. A 2D mesh uses its x/y as physical x/z, retaining all
three modeled stress components. The fixed AWF constants remain
0.09, 2.55, 0.2295, 10.7 and 5.1 with their distinct roles.

The outer approximation omits molecular diffusion, damping, low-Re E and
the finite-Re dissipation split. Transported outer epsilon is identified
with eps_tilde only within that declared approximation. Physical positive
molecular viscosity remains in the AWF. Negative virtual intercepts are
retained. The model includes no coefficient retuning, production clipping,
turbulence floor, filter or SVV.

Outer velocity/pressure/epsilon and full-domain k are independent fields.
The full patch solves internal velocity/pressure stage unknowns and
retains their storage/history through condensation. Its inner epsilon
is dependent, with no independent PDE row. Velocity and k have shared
interface traces; pressure may differ across the model interface.
The initial implementation uses matching interface meshes and CG
subdomain spaces.

The pressure variable is pi=p_mechanical/rho+2k/3, in kinematic units.
Mechanical pressure and the body force enter the forcing sampler
consistently. Full normal traction requires the flat-reference jump
(p_w-p_o)/rho=-11K0/75. Use one physical-volume mean gauge across the
combined closed domain, or the pressure level supplied by physical
traction in an open domain. Do not impose a second independent gauge
or mechanical-pressure continuity.

Inner viscosity uses sampled K, while isotropic stress uses local
transported k. Production uses the same stress as momentum and traction,
including the actual isotropic/divergence contribution where discrete
incompressibility is weak. The inner source is assembled once as

    P_A - (1-chi)*E_A - chi*p_A0

p_A0 is the counterfactual ZERO-forcing production at current K. The
dependent profile E_EP supplies the single outer epsilon value at H.
The outer equation determines its diffusive flux. Source loads are
dual vectors B_Q^T W_J s_Q; coefficient conversion requires a
constrained mass solve.

**Initial numerical and admission choices**

The initial reference fixture has u*=0.4 m/s, nu=1.5e-5 m2/s,
H=2 m, sand-grain h=0.1 m, r_K=r_m=1 m, independent outer
sensor distances r_a=3 m and r_b=4 m, Z=30 m, L=300 m and
span=20 m. Its incoming direction is physical x, with zero body
acceleration. K is evaluated at a fixed physical point. The matching
anchor selects the outer velocity offset and aerodynamic z0; no
independent z0 is imposed. These are illustrative verification
inputs, rather than calibrated terrain parameters.

Begin with degree-four velocity/k/outer-epsilon and degree-two
pressure, Modified_A bases and eight GLL points per direction.
For this basis the corresponding NUMMODES values are five and three.
The pressure pairing must pass the M3 nullspace, constraint-rank and
manufactured-convergence gates.

Split positive wall quadrature at h and positive viscosity/dissipation
breakpoints, including the corrected cutoff. Begin with 16-point
Gauss-Legendre panels, resolving small molecular/reciprocal scales,
and compare independently with refined quadrature. The metric
cutoff and redistribution use the same discrete budget quadrature.
Under discretization refinement hold H, h, physical sampling
distances and transition positions fixed. Changes of those inputs
are separate wall-model sensitivity studies.

The reference target uses a constrained physical L2 projection
preserving constant k, matching velocity and dependent epsilon
trace, followed by explicit algebraic-constraint checks. Sensors
compare the stored projected inlet reference through the same
evaluation and forcing operators. This stored target is distinct
from the current-K counterfactual used in p_A0.

Use the native empty-variant IMEX family, first-order startup and
constant-step order two. The starter dt is 0.05 s. State-dependent
conservative advection, diffusion, reactions, roughness forcing and
interface/pressure coupling belong to the converged implicit block;
prescribed state-independent body forcing is explicit. This
partition is a future wall-element implementation choice; it
does not change the recorded legacy baseline run. Preserve full
native value/explicit-derivative histories and extend the history
layout for the inner/outer fields. M4 must verify the complete
stage equations and stationary time map.

The correction is chi=G_sp G_adm G_st, constant along a normal ray.
Use seven independent sensor weights, all initially one. Diagnostic
thresholds are d_on/off=0.02/0.1, margin_off/on=0.05/0.1 and
Re_off/on=100/200. They are initial declared choices without
universal calibration. Terrain cases must separately supply a
physical mask coordinate, upstream-flat endpoints, trace-support
buffer and geometry certificate. The mask stays zero through
hill/wake/outlet.

Admission requires positive outer k/epsilon and sampled K, valid
physical distances, positive normal stretch/Jacobian and global
normal-map injectivity. The initial minimum normal stretch is 0.1;
coordinate validity does not establish physical thin-layer adequacy.
The corrected cutoff must satisfy r_c<=r_m and the metric budget
condition with smooth interior margins. An inactive inadmissible
correction evaluates the adapted baseline directly. Nonpositive
turbulence or invalid geometry is rejected, not repaired by chi=0.

Eta defaults to zero and is zero for terrain. The optional flat
verification anchor uses a fixed global eta=1, one stored target
and free differential loads only. Algebraic constraints are
unchanged; spatial eta and solution resets are excluded. Stage
consistency and the anchor's budget/work remain explicit M4 gates.

Residuals use the physical dual mass norm of free differential
loads; state changes use physical mass norms. Scales are
u*, K0, u*^3/(kappa_c H) and u*^2 for velocity, k, epsilon and
pressure, with time scale H/u*. Normalize by the corresponding
field/rate scale and square root of physical volume, using
surface norms for interface data. Audit removed/constrained
rows separately. Initial scaled residual/time-map and coupled
constraint tolerances are 1e-9, with 1e-10 for isolated interface
identities.

**Equation to implementation and verification mapping**

The contract maps all 178 equation labels in the revised note
exactly once, in 56 groups, to 11 components and 35 verification
identifiers. Existing kernels, known failures, comparisons and
planned work have distinct statuses. The complete matrix is
stored in the JSON artifact; its checks reject missing, duplicated
or unknown labels and unresolved component/test references.

| Main groups | Component and verification ownership |
| --- | --- |
| Original tensor, high-Re sources, Yap, neutral roots/profiles, scalar AWF and integrals | Existing RANSModels/AnalyticalWallFunction; verified M1 kernel oracle |
| Finite-Re terms, MOST/FV and scalar Robin comparisons | Scope/reference entries; no new low-Re or FV implementation is claimed |
| Physical sampling, cutoff, redistribution, weak/modal loads | Planned equilibrium profiles; M2 |
| Multiregion spaces, inner k/epsilon, flat storage and pressure coupling | Planned field layout, patch and interface; M3 |
| Composite target, constrained residual, anchor and stationary time map | Planned reference and multiregion history components; M4 |
| Normal geometry, variable H, conservative momentum/mass/k flux and pressure | Planned geometry, coupled patch and interface; M5 |
| Spatial mask, independent sensors, derivatives and exact adapted-baseline recovery | Planned transition and trace support; M6 |
| Full energy and matching-height sensitivity | System/reference diagnostics; M7 |
| Constitutive adverse directions | Known failure verified in M1 |
| Sampled-K nonnegative-cone obstruction | Mathematical failure; native counterexample remains M2 work |

Planned file names and verification identifiers are implementation
contracts, not available XML settings or registered native tests.
Changing the model, pressure convention, sampling, time partition,
geometry assumptions or discretization requires an identified
contract revision and rechecking the affected milestones.

**Baseline evidence and closure audit**

The exact-base Release build and all 16 baseline tests passed
during M1, before the new kernel checks were built. Their
[original receipt](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m1/baseline-receipt.json)
retains commands, exits, times and raw-log hashes. M0 reuses this
result because it changes only the contract, documentation and
provenance checker.

The M0 checker verifies those raw logs and the retained exact-base
worktree, inventories the 16 registered tests without executing
them again, and records the source tree, nine native binary hashes,
resolved library hashes, compile flags, reference XML spaces/time
scheme and actual toolchain. Fingerprints captured at M0 audit
are labelled separately from the earlier numerical execution.
It also verifies that the tested M1 sources and common runtime
dependencies still match their committed evidence.

The actual framework is 5.9.0, Release with MPI enabled. The note
discusses 5.10.0 interfaces; native availability for the new
multiregion solver remains an M3 gate. The retained baseline used
16 compiler workers and sequential single-process tests pinned
to CPU 0, with OMP/OpenBLAS/MKL/BLIS threads set to one.
The M0 read-only audit uses one CPU and launches no numerical
solver. The project hard limit remains 72 total Linux workers.

All generated data is under
/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised. The
unrelated primary checkout and installation are preserved.
The audit requires a fresh output directory:

    source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
    CLS_EP_RUN=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised
    export TMPDIR=$CLS_EP_RUN/tmp
    taskset -c 0 python3 "$CLS_EP_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m0.py" \
      --output "$CLS_EP_RUN/m0-audit-001"

The M1 constitutive failure remains explicit. M0 defines admission
to the formula domain, not a stable state envelope or a general
positivity guarantee. Robust terrain qualification still requires
the separately identified model decision in conditional M8.
