**CLS and EP-AWF revised formulation implementation and validation plan**

Status: Approved by the user on 2 November 2026; M1 closed.
Implementation branch: codex/cls-epawf-revised, based on
6415398d856323c3116e646eaf7c882c66617602.

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: plan.
- Origin Date: 2026-11-02.
- Verification Status: M1 kernel identities and constitutive counterexamples
  VERIFIED on Linux; later implementation and numerical milestones UNVERIFIED.
- Version Label: cls-epawf-revised-plan-v1.

The objective is a native implementation of the high-Re CLS model and the
spectral/hp EP-AWF specified by the revised note, with independent verification
of its equations, assembly, interface coupling and time integration. The
principal work is a wall element with transported inner k, dependent inner
dissipation and a coupled inner velocity/pressure patch. The existing branches
supply reference kernels, solver components and experimental infrastructure.

The [M1 closure report](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-m1.md)
records the exact tested revision, independent errors, documented negative
constitutive directions and passing 17-test Linux suite. M2 and later
implementation milestones remain open.

The revised note retains two demonstrated limitations: adverse incremental
momentum diffusion and failure of general inner-k positivity. A faithful
implementation must reproduce these results. Stable long-time terrain
qualification needs an explicitly selected further remedy; success of
equilibrium or nonlinear-solver tests cannot establish that qualification.

**1. Governing evidence and inspection**

The source specification is
[the revised note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex).
Its SHA-256 is
d35e437c2f2f8b25a46f7a3f285c6f1ad8b48d7b554aa1d161315f59eeb6de0f.
The comparison document is
[the original note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note.tex),
SHA-256
8b7f67035e553318a2b0896c9e31dd74d94bf0cc4864013ea209a48eddf98819.

The revised file shares its main derivation with the original file but removes
the appended common-cubic-scaling and capped-destruction remedies and their
interval certificate. Those remedies are outside this implementation contract.
The revised note's stability section discusses possible changes without
selecting a complete stabilized model.

[README.md](/Users/jay/Documents/nektar-v2026/README.md) and
[AGENTS.md](/Users/jay/Documents/nektar-v2026/AGENTS.md) govern the project and
execution environments. The inspected reference tips form an ancestor chain:

| Branch and inspected tip | Existing capability | Consequence for this project |
| --- | --- | --- |
| codex/neutral-abl-rans-awf, cba0df5d985ad863d438c48dd2e094769874e1b6 | Original high-Re CLS tensor/Yap/reference kernels; original and corrected AWF; integrated CV strip, outer CG solver, histories and flat-domain studies | Reuse verified kernel formulas and regression controls; implement the revised wall-layer representation separately |
| codex/hill-wall-geometry-transition, c48d9fb14afeca6caa80cf5b3f8ce0a34a02d2d7 | Normal-coordinate terrain sampling, metric strip integrals, integrated strip continuity, vector boundary coupling and upstream mask | Reuse geometric and transfer ideas after reconciling curvature conventions; extend to the full coupled patch and variable H |
| codex/cao-tamura-calibrated-validation, 6415398d856323c3116e646eaf7c882c66617602 | Native physical observations, residual-map/Newton diagnostics, provenance and calibration workflows | Use as branch base and diagnostic infrastructure; its CLS development gate remains unmet |

The physical CLS kernel and scalar AWF formulas remain the original laws across
these tips. The terrain branch adds branch-selected constitutive linearization
for nonlinear solves; it does not retune the physical CLS stress.

The current implementation stores eliminated strip states through P/H/N
interpolation and couples mean production minus mean dissipation to the outer
boundary. It has no transported polynomial inner-k field or the revised
all-mode redistribution. The terrain adapter remains a thin-layer CV/CG
formulation, rather than the note's full-strain momentum/pressure patch.
These differences are visible in
[RANSWallStrip.h](/Users/jay/Documents/nektar-v2026/library/SolverUtils/Turbulence/RANSWallStrip.h:17),
[the flat adapter](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSAnalyticalWallStrip.cpp:391)
and [the terrain adapter](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSTerrainWall.cpp:826).

The existing solver requires matching velocity/k/epsilon meshes and quadrature;
the history schema also assumes this layout. The revised multiregion spaces
therefore require solver and checkpoint work, beyond adding a boundary tag.
See [field initialization](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSVelocityCorrectionScheme.cpp:98)
and [history serialization](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems/RANSHistory.cpp:90).

Prior studies document degree-four flat preservation and selected perturbation
controls, a degree-five 3D growing mode, and short terrain pilots. The
[CV3-R1 record](/Users/jay/Documents/nektar-v2026/docs/plans/cao-tamura-cv3-r1-reference-feasibility.md:14)
records developed standard controls and rejected or retired CLS cases. Those
are historical results for their producing binaries; they do not validate the
new formulation or establish the cause of CLS nonlinear rejection.

**2. Model and numerical choices proposed for approval**

| Component | Selected contract |
| --- | --- |
| Outer closure | Original CLSHighRe: full cubic stress, original C(q), max strain/rotation invariant, algebraic Yap correction and original transport constants |
| Dissipation convention | Transport outer total dissipation identified with modified dissipation under the stated high-Re approximation; retain physical molecular viscosity in the AWF, omit outer molecular terms in the reference-preserving model |
| Geometry and samples | Physical matching height H, roughness h, independent k-sampling distance r_K and dissipation-tail start r_m; outer sensor distances r_a,r_b greater than H |
| k sampler | Start with fixed physical point evaluation k(r_K); a normalized-average sampler would be a separately identified wall-model option |
| Independent state | Outer velocity/pressure/dissipation; k throughout inner and outer regions; full inner velocity/pressure stage unknowns retained through patch condensation |
| Inner dissipation | Dependent profile, with one matching value supplied to outer dissipation; no independently solved inner dissipation equation or additional prescribed epsilon flux |
| Inner momentum | Full-strain auxiliary-viscosity stress, conservative transport, roughness force, normal momentum and storage; analytical columns provide the flat quasi-steady reduction |
| Discretization | Begin with matching interface meshes and subdomain CG spaces, a shared velocity/k trace, and separate inner/outer pressure spaces permitting the required pressure offset |
| Pressure pairing | Begin with the existing p / p-2 continuous velocity-pressure pairing within each subdomain; qualify its nullspace, constraint rank and manufactured convergence before accepting it |
| Time integration | First-order stages for isolation, then constant-step IMEX2 with complete inner/outer histories and a converged coupled stage residual |
| Correction | Ray-constant chi, physical upstream mask, independent state sensors and smooth interior admissibility margins |
| Discrete anchor | Default eta=0; optional fixed global eta=1 only for isolated flat-reference verification; eta=0 throughout terrain calculations |
| Invalid states | Report and reject states outside the declared positive-k/epsilon and geometric regime; record failure location and inputs |

No cubic scaling, C(q) continuation, production clipping, hidden turbulence
floors or positivity repair is included in this selection. A trial-state
line search may protect residual evaluation, but it does not establish
positivity of an accepted time-dependent solution.

Use p=4 as the initial comparable configuration, then p=2,4,6,8 where the
chosen mixed spaces are valid. Record actual basis degrees and quadrature;
NUMMODES is not assumed to equal degree. Refine quadrature independently of
degree. Keep H, h, r_K, r_m and the physical sensor and transition locations
fixed during discretization studies.

**3. Implementation architecture**

Introduce an opt-in equation system, provisionally RANSEquilibriumWallSystem,
with an explicit inner/outer block layout. Keep existing solver selections
and their regressions available. The proposed name is a development choice,
not an existing XML capability.

Reuse
[RANSModels](/Users/jay/Documents/nektar-v2026/library/SolverUtils/Turbulence/RANSModels.cpp)
for the original CLS/reference law and
[AnalyticalWallFunction](/Users/jay/Documents/nektar-v2026/library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp)
for stable original-AWF integrals and reconstruction. Add separate components
for equilibrium profiles/modal loads, normal-tube geometry, coupled patch
residual/condensation, common interface transfers and multiregion histories.
Factor reusable native matrix-free and pressure operations where needed;
keep the model-dependent wall residual distinct from its preconditioner.

The mathematical core of the inner k source is

    p_A0 = nu_t,A * (t0 / D_A)^2
    E_EP = (1 - chi) * E_A + chi * E_c
    C_k  = chi * (E_c - p_A0)
    S_k  = P_A - (1 - chi) * E_A - chi * p_A0

Assemble the simplified net source once. The reference p_A0 uses zero forcing
and the same current sampled K. Current forced production must remain P_A.
On curved rays solve the metric-weighted cutoff equation; the old flat
exponential cutoff is insufficient there. Use the same physical integration
measure for the cutoff, zero-integral redistribution and source loads.

Within the patch, use the full-strain stress and its consistent production.
If discrete divergence is nonzero, retain the isotropic term in the actual
stress contraction and account for the chosen pressure convention.
Represent source loads as dual vectors B_Q^T W_J s_Q; use a mass solve only
when converting a load to a coefficient field. Do not integrate a load twice
through a physical-source API.

At the matching surface use one mass flux, one full momentum traction and
one k flux with equal and opposite weighted contributions. Distinguish wall
traction t_A from matching traction T_H=t_A+F(H-h). Mechanical pressure is
reconstructed from the solver pressure. At the flat CLS reference, enforce

    (p_w - p_o) / rho = (0.52 - 2/3) * K0 = -11/75 * K0

through full traction matching and retained relative pressure modes. Apply
one compatible global gauge. The dependent epsilon trace is E_EP(H);
the outer equation determines its flux.

**4. Milestones and completion evidence**

1. **M0 — Freeze the contract and establish the working branch.**
   After approval, create codex/cls-epawf-revised from the inspected latest
   reference tip, preserving all existing local files. Record note hashes,
   branch tips, compiler/library configuration, pressure convention, field
   spaces, time stages, model parameters and admissibility rules. Create an
   equation-to-component-to-test matrix using the revised note's labels.
   Run the existing small Linux regression suite once against the exact base
   revision in an isolated build. Resolve any baseline failures before using
   that build as a comparison. Deliver the executable specification and
   baseline manifest.

2. **M1 — Verify the original CLS and scalar AWF independently (closed).**
   Check the complete 3-by-3 tensor, rotation/sign conventions, production,
   Yap source and geometric distance, including 2D mean flow with all three
   modeled stress components. Reproduce q*=3.190675081624288,
   Cmu*=0.0982278953282113 and kappa_c=0.440553905186659. Verify the covariance
   entries and eigenvalues, AWF smooth/rough branches, negative virtual
   intercept, closed integrals, no slip, energy identity, signed/reversed
   traction and finite forcing increments against independent quadrature.
   Test both limiting invariant branches and the indexed cubic contraction.
   Reproduce the adverse transverse coefficient -0.9645579066*nu_t and the
   nearby smooth 2D adverse direction. These tests pass by reproducing the
   counterexamples while recording the physical stability gate as failed.

3. **M2 — Implement the spectral/hp wall profiles and source loads.**
   Separate r_K and r_m from the historical P/N interpolation. Implement
   E_c, the monotone cutoff/admissibility test, chi and the simplified S_k.
   Add a bracketed metric cutoff solve and split physical quadrature at
   roughness, viscosity and dissipation breakpoints. Verify every active
   native basis moment, source/flux integration-by-parts equivalence and
   zero redistribution endpoint fluxes. Reproduce the note's nonzero
   mean-only higher moments as a control. Test chi=0, intermediate chi and
   chi=1, forcing variations and inadmissible inactive branches. Include
   the smooth inner-k local-zero counterexample away from r_K; record its
   loss of positivity without changing the source law.

4. **M3 — Build the final field layout and an uncondensed flat patch.**
   First demonstrate native creation, evaluation and assembly of the
   subdomain spaces with full-domain k and outer-only dissipation. Establish
   the pressure nullspace and interface constraint rank. Implement ground
   no slip, zero ground k flux, shared traces/fluxes, the pressure offset,
   physical forcing sampler and inner storage. Use the one-sided outer
   conservative tangential forcing specified by the note, with consistent
   pressure/body-force cancellation. Compare a steady, parallel patch with
   the analytical AWF at zero and nonzero forcing. Then use manufactured
   evolving fields to verify full transport and storage. Condense only
   after the uncondensed residual is verified; compare residuals and
   directional actions at identical states and pressure gauges. Deliver a
   native flat wall-element solver and small reproducible fixtures.

5. **M4 — Verify the constrained flat reference and actual timestep.**
   Construct the declared finite-dimensional composite target, with exact
   matching trace and constant k, and recompute dependent fields from it.
   Check divergence, gauge, ground values, interface constraints and all
   free momentum/k/outer-epsilon residual blocks. Report eta=0 defects
   first and demonstrate their refinement behavior. If exact target
   preservation needs anchoring, subtract only fixed free differential
   residual loads and report their integrals and energy work. Derive the
   corresponding reference balance for every time stage and pressure
   boundary path. Verify one step and repeated steps with consistent
   histories, several dt values, and IMEX startup. Check 2D/full 3D,
   rotated flow, serial/MPI and uninterrupted/full-history restart
   equivalence. Exact anchoring is a verification mode with separately
   reported changes to the discrete operator.

6. **M5 — Extend and verify the conservative curved patch.**
   Progress from a rotated plane to a smooth circular patch, constant-H
   curved patches, variable H and a genuinely 3D smooth surface. Reconcile
   the old graph-curvature sign with B=grad_Gamma(n) before reusing geometry.
   Verify positive normal stretches/Jacobians, tube injectivity, between-node
   mapping validity and matching geometry. Use the actual tilted matching
   normal and area. In particular verify
   m_H=J_H*(w_H-v_H dot grad_Gamma H), not J_H*w_H alone.
   Test discrete metric identities, constant Cartesian flux divergence,
   full normal momentum, normal traction, conservative momentum/k/mass
   budgets and roughness-force work. Check lateral fluxes across patches,
   patch partition independence and MPI ownership. Repeat M4 reference and
   stage checks after curved assembly is enabled. Treat C1 hill feet by
   explicit patches and interface conditions; use a separately labelled
   smooth hill for initial verification.

7. **M6 — Implement the prescribed upstream transition and baseline.**
   Use a prescribed physical wall coordinate and C2 spatial mask ending
   on exactly flat terrain before the hill, with no downstream
   reactivation. Use independent outer velocity/increment, normal-motion,
   k, epsilon, pressure and convection indicators with positive weights.
   Compare sensors with the same evaluated discrete reference. Add smooth
   cutoff-margin and high-Re admissibility gates. Include K, cutoff,
   sensor, chi, trace-lifting and condensation dependencies in directional
   actions, including the trace term (E_c-E_A)*grad(chi). Use a correction
   projection constrained to zero on designated baseline faces, allowing
   a physical transition buffer. Test equal-input identity with the same
   new patch operator at chi=0 and eta=0. Compare complete residuals,
   traces and fluxes; a solution from a different upstream history is a
   separate comparison. Keep eta=0 for every terrain/transition run.

8. **M7 — Characterize admissibility and perturbation behavior.**
   Check accepted states at source/volume quadrature, interfaces, physical
   samplers and reconstruction locations. Distinguish sampled positivity
   from any claim about the minimum of the polynomial everywhere.
   Perturb velocity, k and epsilon independently and jointly; include the
   three-component transverse mode supported outside the layer,
   near-interface disturbances and transition-crossing disturbances.
   Compare amplitudes and wavelengths with independent dt/p/mesh
   refinement. Record true coupled residuals, weak/strong continuity,
   pressure/normal traction, covariance eigenvalues, tangent diagnostics
   and mass/momentum/energy budgets. Use durations in both flow-through
   and turnover units; a five-flow benchmark is a comparable control,
   not an arbitrary-time stability criterion. Reproduction of adverse
   high-frequency growth or positivity failure closes the faithful-model
   diagnosis and leaves robust-model qualification failed.

9. **M8 — Qualify a robust variant and validate terrain conditionally.**
   Before a production-scale or physical-accuracy campaign, resolve the
   demonstrated constitutive and positivity failures with a separately
   identified model decision. Prepare the concrete remedy and its
   changed equations/reference/energy/baseline contract for approval
   if that extension is requested. Re-run all affected earlier gates
   after such a change. For a qualified candidate, compare paired
   empty-floor/hill cases, surface pressure, wall stress/reversal,
   separation/reattachment, velocity, k and modeled Reynolds stresses.
   Compare spatial EP/transition with an all-baseline calculation using
   the same qualified patch and physical inputs. Reuse the frozen
   Cao–Tamura data/observation conventions, but version the old
   original-first selector/freezer before any new calibration. Fit flat
   precursor data only and retain reserved flat/hill scores for
   validation. No result from an unresolved or rejected trajectory
   enters a developed-flow average or an accuracy claim.

10. **M9 — Deliver the implementation and evidence.**
    Provide documented model selections, runnable sessions, checkpoint
    and observation workflows, compact source/binary-bound evidence,
    regression results and refinement figures. Update the repository
    guides with the new wall representation and distinguish baseline
    verification, known failed model properties, bounded numerical
    behavior and physical validation. M9 can deliver a faithful
    research implementation with explicitly failed stability/positivity
    gates; a robust terrain-validation deliverable additionally needs
    the conditional M8 work.

The practical dependency order is M0, M1/M2, M3, M4, M5, M6, M7 and M9.
M8 follows a resolved model decision and repeats the affected gates.
This ordering keeps geometry, source assembly, pressure coupling and time
integration failures individually diagnosable.

**5. Validation matrix and proposed acceptance rules**

These are proposed project tolerances, rather than claims made by the note.
M0 will bind each test to units, reference scales, norm domains and the exact
free/constrained rows it measures. Use a physical dual norm for residual
loads and a physical mass norm for state changes; do not normalize by a
reference residual that should vanish.

| Test family | Initial acceptance rule |
| --- | --- |
| Roots, stress identities and original AWF oracle | Relative error at most 1e-11 with stated absolute treatment for zero quantities; reproduce the appendix example's tighter 2e-12 to 3e-12 absolute modal identities where applicable |
| Wall loads, cutoff and redistribution | Independent split quadrature agreement at most 1e-10 after dimensional scaling; all active modes checked, including a demonstrably failing mean-only control |
| Shared interface and constraint identities | Scaled algebraic mismatch at most 1e-10 in isolated probes; coupled constraints at most 1e-9 |
| Native reference residual and stationary step | Every reported scaled field block at most 1e-9, with tighter solve tolerances; pressure compared in one gauge; separate eta=0 and anchored results |
| Condensation and residual directional actions | Agreement with uncondensed residuals to solve tolerance; stable difference-step convergence and separate limiting-branch checks |
| Manufactured spatial/temporal behavior | Correct discretization consistency and expected temporal order on at least three dt levels; independent quadrature and p/mesh studies without simultaneous model-input changes |
| Serial/MPI and full-history restart | Initially at most 1e-8 scaled full-field/flux difference for small fixtures; any tolerance adjustment justified by reduction/solve accuracy before campaign acceptance |
| Engineering flat-profile drift | Retain the historical 0.1 percent over five flow-throughs as a comparison metric, separate from exact-reference and perturbation gates |
| Stability and positivity | Report the original-law counterexamples as failed physical gates; no numerical tolerance can turn their signs into a stability or positivity proof |
| Physical accuracy | Freeze observables, data uncertainty, coverage, calibration split and acceptance thresholds before fitting; numerical errors must be assessed separately |

The independent oracles should differ from production arithmetic: indexed
versus matrix stress, quadrature versus closed integrals, source versus flux
weak loads, analytic manufactured derivatives, and uncondensed versus
condensed assembly. Checking two calls to one kernel is insufficient.

The proposed small-test matrix includes smooth and rough AWF branches,
original/EP wall policies, fixed/switched weights, 2D/full 3D, fixed/variable H,
and 1/2/4 MPI ranks, with an 8-rank ownership/restart case once justified.
Use standard k-epsilon as a control to isolate wall/pressure errors from CLS
constitutive behavior. Larger matrices and durations follow measured small
pilots; elapsed solver runtime is recorded, not predicted as validation.

**6. Linux execution, synchronization and resource plan**

At inspection the Mac is on the proposed base at 6415398 with untracked
literature and unrelated local files. The authoritative Linux checkout is
clean on codex/hit-inflow-generator at
f4b70baec10016dd098005bc5a5507b87e4a57bf. The USB SSH connection timed out;
the Tailscale connection succeeded. No source synchronization, build or
solver run was performed during this planning inspection.

After approval, compare revisions and working trees again. Transfer committed
project changes by an explicit Git object transfer/fetch and check out the
exact commit in a dedicated Linux worktree, without resetting either
existing checkout. Verify source manifests before compilation and record
any staged/uncommitted patch if validation ever intentionally includes one.
Use this run root:

    /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised

Keep its source worktree, build, optional isolated installation, cases,
checkpoints, logs, figures and tmp directories there. Begin with the existing
external native build pattern against the installed Release MPI libraries.
If implementation requires changed core libraries, use a complete isolated
build/install and matching library paths. Preserve the shared HIT checkout,
build and installation.

The initial build/smoke command pattern is:

    source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
    CLS_EP_RUN=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised
    CLS_EP_SOURCE=$CLS_EP_RUN/source
    export TMPDIR=$CLS_EP_RUN/tmp
    mkdir -p "$TMPDIR"
    cmake -S "$CLS_EP_SOURCE/solvers/IncNavierStokesSolver/Tests/RANS" \
      -B "$CLS_EP_RUN/build" -DRANS_BUILD_SOLVER=ON -DCMAKE_BUILD_TYPE=Release \
      -DNektar++_DIR=/home/jay/.local/opt/nektar-v2026-mpi/lib64/nektar++/cmake
    cmake --build "$CLS_EP_RUN/build" --parallel 16
    ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure

Use 16 compiler workers initially and 1/2/4/8-rank solver tests with one
numerical thread per rank, scheduled separately by default. Check concurrent
server work before each allocation; all builds, MPI ranks, test runners and
their threads together must remain at or below the AGENTS.md hard limit of
72 workers. Historical README allocations of 16 describe earlier campaigns.

Each native case must retain the exact command, commit and source hashes,
binary/library hashes, mesh/basis/quadrature, boundary and forcing data,
pressure/gauge convention, dt/history configuration, MPI/CPU/thread allocation,
exit status, wall time and acceptance results. Monitor owned process status,
per-stage progress and declared metrics. Set and record pilot-informed wall
timeouts for longer cases; preserve the last accepted state and all rejected
stage evidence. Run-generated data stays outside the repository.

**7. Source changes and deliverables**

The expected change areas are:

- New kernel/profile and geometry helpers under
  /Users/jay/Documents/nektar-v2026/library/SolverUtils/Turbulence.
- New block equation-system, patch/interface and history components under
  /Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/EquationSystems.
- Native probes under
  /Users/jay/Documents/nektar-v2026/library/UnitTests/SolverUtils, and small
  reproducible definitions/oracles under
  /Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Tests/RANS.
- A separate revised-formulation example and updated user/developer guidance,
  with compact milestone receipts under
  /Users/jay/Documents/nektar-v2026/docs/plans.

Register production sources and focused tests in the corresponding CMake
files. Extend the complete build rather than validating only an ad hoc kernel
binary. Keep large meshes/fields/logs and historical source copies on Linux;
commit only compact fixtures, manifests, summaries and necessary figures.

Every milestone ends with reviewable changes, its source revision, commands,
results and remaining gate status. Implementation does not automatically
unlock a failed model qualification or the existing Cao–Tamura score gates.

**Approved scope**

The user approved the branch/base, the faithful revised-note model and the
M0–M7/M9 implementation/verification sequence, and requested execution and
closure of M1. This includes recording the
original CLS and sampled-K failures as validation findings. A robust-model
extension and the conditional physical-validation campaign in M8 require
a concrete, separately identified resolution of those failures; the removed
remedies supplement is not adopted by this approval.
