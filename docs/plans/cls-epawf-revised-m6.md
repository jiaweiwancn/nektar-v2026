# M6 — Prescribed upstream transition and adapted baseline

Status: M6 closed; implementation, Linux qualification and evidence audit passed.
Branch: `codex/cls-epawf-revised`; M5 base `71b5e1e`.
Final numerical revision: `8a80594a47eb34edbb76a2fe827a200227b3d34b`.
Complete main qualification/regressions: `2d5d175f98dffa543892e453569b57109d6d55c7`.

## Material Passport

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: run; source implementation and independent numerical verification.
- Origin Date: 2026-11-03 (Asia/Shanghai).
- Verification Status: VERIFIED for M6 within the declared native backend and test scope; original physical failure gates unchanged.
- Version Label: cls-epawf-revised-m6-v1.

The user authorized execution and closure of the approved M6 scope. The
revised note and contract v1 remain frozen. Tailscale SSH to
`jay@100.105.247.110` is the only transport used. The shared Linux HIT source,
build and installation are preserved; the isolated revised-project source,
build and run directories remain authoritative for numerical conclusions.

The implementation and acceptance work is:

1. Implement analytic C2 spatial, state and admissibility gates, all seven
   independent vector/scalar sensors, their first/second directional actions,
   and a prescribed physical wall coordinate. Keep physical positions fixed
   under refinement; reject correction support without an exactly flat fetch.
   A cosine hill uses signed wall arc length and explicit foot patches.
2. Evaluate the current and stored projected flat reference with identical
   native physical maps at H and fixed independent outer distances r_a,r_b.
   Retain pressure/body and convection separately. Qualify the reference on
   both unit admissibility plateaus. Inactive branches evaluate the original
   profile directly without solving a corrected cutoff.
3. Couple the weight to current K, cutoff, source, dependent epsilon and full
   trace constraint. The outer sensor polynomial contains trace coefficients,
   so the trace remains an implicit coupled constraint in stage Newton.
   Differentiate the complete residual, including those columns and the
   physical trace gradient term (E_c-E_A) grad(chi).
4. Project the baseline trace and correction separately. Use native
   endpoint/edge preserving face transforms and constrain the correction to
   zero on every designated baseline face/shared mode. Require a physical
   buffer before the hill foot and a native face boundary at plateau start.
   Evaluate analytic chi in physical space; do not project the gate itself.
5. Verify equal-input chi=eta=0 identity using two instances of the same full
   curved CG operator. Compare unmasked residuals, constraints, traces and
   fluxes; retain storage, pressure, normal momentum, full-strain stress and
   retained k transport. Different upstream histories are separate comparisons.
6. Qualify kernel arithmetic with an independent high-precision oracle,
   native 2D/full-3D sensors and weak operators with independent Jacobi bases,
   all-mode correction/full-budget identities, trace projections, directional
   refinement, and deliberate omitted-term controls. Run strong manufactured
   transition cases, spatial refinement, full/condensed stages, MPI ownership
   and native full-history restart. Run the appropriate earlier regressions.

Initial acceptance tolerances remain those in the approved plan: isolated
scaled identities 1e-10, coupled constraints/stages 1e-9, MPI/restart full-field
comparisons 1e-8. Directional checks must demonstrate a stable difference-step
range and separately cover gate endpoints and inactive inadmissible branches.
Every terrain/transition run requires eta=0. No constitutive retuning, floors,
clipped cutoffs, independent inner epsilon PDE or spatial residual anchor is
part of this milestone. The original stability and positivity gates remain
failed, as specified; M7 and the conditional M8 campaign remain separate.

## Implemented operator

`EquilibriumAWFTransition` implements the prescribed C2 decreasing polynomial,
its first/second directional actions, all seven independent physical sensors,
their normalization, strictly positive weights, and the margin/Reynolds gates.
The wall coordinate is independent of velocity and physical normal distance.
It is signed wall arc length for the extruded cosine hill, with a prescribed
positive incoming orientation. An affine-coordinate helper also covers a
general 3D plane. An active upstream selector rejects geometries without an
exactly flat fetch, a variable upstream H, a missing plateau/baseline face,
unsafe Reynolds thresholds and buffers that reach the hill foot.

`RANSWallEvolutionTransition` evaluates the actual stored flat native target
and current state with the same samplers and forcing operators. The fixed
reference body force is zero. Pressure/body and convection remain independent,
including in an IMEX stage that moves prescribed body loads to the explicit
side. Prescribed body quadrature remains fixed; it is not differentiated as
a spurious state-dependent partition. Profiles reuse a qualified cutoff only
when K and geometry agree, and recompute the state-dependent weight on every
residual evaluation. Inactive branches use the baseline directly and report
zero corrected-cutoff iterations, including the positive low-K/low-Re control.

The matching epsilon coefficients remain unknowns because the independent
outer sensor sample contains those same coefficients. The physical residual
lifts ordinary essential data but retains these dependent trace columns;
their rows impose the nonlinear projected trace. Newton and Schur actions
differentiate this complete coupled residual. Initialization solves the trace
lifting fixed point. The current dependent inner epsilon API requires time
and forcing data in Upstream mode; its earlier two-argument form is restricted
to Fixed mode. No extra inner epsilon transport or flux equation is introduced.

Baseline and correction traces use separate native constrained transforms.
Segment transforms preserve endpoints; 3D uses commuting segment transforms
to preserve shared edges/corners. The correction is explicitly zero on every
designated baseline face and shared coefficient. A flat physical buffer allows
polynomial support beyond the analytic mask endpoint, while plateau start and
baseline start lie on native face boundaries. The matching trace takes
priority over side-function epsilon values at shared H corners; this policy
is recorded in output and the user guide. The gate itself is never projected.

The energy diagnostic now contracts storage and stored energy with the actual
physical stage mass matrix. Constitutive quadrature can use a different order;
using it for storage caused a diagnostic defect on the unsplit curved panel,
despite independently matching residual rows. This change corrects accounting
and leaves the stage equations intact.

## Verification design and observed results

The complete main implementation qualification is at
`2d5d175f98dffa543892e453569b57109d6d55c7`. It was compiled in Release with
16 compiler workers in the isolated Linux build. All numerical threads are
one; native cases use 1/2/4 ranks. A supplementary quadrature run uses one
worker on CPU 4 alongside at most four main qualification ranks. No compiler
job overlaps that supplement. Exact commands, wall times and source/runtime
hashes are retained with the evidence. Earlier regressions use 12 CTest worker
reservations on CPUs 0-15, respecting each test's MPI PROCESSORS allocation.
Compilation finishes before this test allocation begins. All settings stay
below the 72-worker server limit.

The independent kernel oracle uses 80-digit arithmetic and independent
numerical differentiation, not the production gate derivative. It covers
all 18 input and 15 control directions, all 40 physical-normalization input
components, C2 endpoint approaches, independent sensor activation, pressure/
convection cancellation controls, signed physical coordinates, and actual
profile K/cutoff/source/trace differentiation on three step levels. Its
1,654 checks pass; the native probe also rejects 34 invalid inputs.

The native operator oracle reconstructs geometry, Jacobi bases, field values,
the seven sensors, dimensional AWF primitives/cutoff, indexed CLS stresses,
weak loads and the constrained GLL trace projection independently. Native
maps/signs and coefficient vectors specify the FE space; native stresses,
source values or weights do not construct the oracle. Both 2D geometries and
full 3D pass. The maximum independent complete weak-row absolute discrepancy
is 8.44e-15. Current/reference sensors agree within 2.4e-16 after scaling;
the actual projected reference departure is independently zero. Designated
baseline correction rows and the equal-input native baseline identity are
exactly zero on common compatible boundary data satisfying the matching-trace
constraint. The independently reconstructed baseline agrees within 7e-15.
Upstream at chi=0 defines the adapted baseline with the same implicit trace
rows and corner priority. Its additional comparison with Fixed mode is on that
common admissible manifold, or after trace elimination. Fixed strongly lifts
the dependent trace; Upstream retains its coefficients in physical loads.
Their residual extensions away from that constraint manifold are different
and are not claimed identical.

Complete 2D K, outer epsilon and pressure actions are independently rebuilt
on two difference steps, retaining dependent H coefficients. Agreement with
the finest native action is at most 1.26e-9. Frozen-chi controls detect the
missing implicit trace derivative. The 3D case independently checks complete
weak rows and tensor trace projections; its time evolution is not part of
this new M6 matrix. Actual new transition time/Schur/restart tests are 2D.

Strong manufactured cases use U=pi=0, spatially constant k(t), and
epsilon=B(ell,t)/r. Their prescribed differential sources are not native weak
loads. Deliberately broad state thresholds isolate the spatial transition;
these verification thresholds are not physical calibration. The source
expressions are checked in the actual installed interpreter against separate
formula evaluations before solver qualification.

| Geometry/refinement | epsilon physical RMS | matching-trace tangential-gradient RMS |
| --- | ---: | ---: |
| Flat p=3 | 1.57514e-4 | 2.57027e-4 |
| Flat p=4 | 2.05868e-5 | 1.25487e-4 |
| Flat p=5 | 7.08771e-6 | 1.00359e-4 |
| Flat p=3, h refined by two | 1.22098e-5 | 5.81042e-5 |
| Cosine feet p=3 | 1.57521e-4 | 2.57032e-4 |
| Cosine feet p=4 | 2.05891e-5 | 1.25490e-4 |

Physical mask positions, H, roughness and sampling distances remain fixed.
The trace-gradient control omitting (E_c-E_A) grad(chi) has RMS 6.338e-4;
the resolved gradient approaches the correct nonzero value under refinement.
These data establish consistency, not a universal asymptotic p-rate for a C2
mask whose endpoint lies inside its projection-buffer face.

The full/condensed physical state difference is below 2.5e-17. MPI 1/2/4
physical-field differences are below 5.9e-17 and full native history differences
below 1.7e-16. Four-rank restart restores every native history slot, pressure
and gauge. IMEX1 and the standard k-epsilon control pass. The full M6 entry
point passes 2,753 checks in 33 commands, with 154 additional independent weak
checks. Its registered kernel test passes all 1,654 checks at the same revision.

At fixed p=3, mesh, physical inputs and final time 0.34 s, actual native IMEX2
scalar mass-norm differences to dt=0.0025 s are 2.70489e-6, 6.83422e-7 and
1.31722e-7 at dt=0.02, 0.01 and 0.005 s. Their successive reduction factors
are 3.96 and 5.19, consistent with second order while separating spatial
error from the temporal comparison. No fitted/extrapolated CFD result is used.
All 11 invalid-session controls are rejected as required: nonzero anchor,
sensor at H, absent buffer, crossed plateau/baseline face, zero weight,
unsafe Reynolds threshold, inadmissible reference, reversed direction,
absence of an exactly flat fetch and correction support overlapping the hill.

Doubling tangential quadrature from 12 to 24 on the same p=3 mesh changes
epsilon by 4.22992e-8 physical RMS, 0.0269 percent of its measured spatial
error. The supplement preserves the exact source and solver hash and uses
one numerical worker; its stage constraints pass. This bounds that fixture's
quadrature sensitivity and does not certify every arbitrary switched state.

The source/runtime fingerprint covers 3,891 numerical/configuration inputs
and 98 executables/libraries, plus numerical package/tool code and the fixed
thread environment. All 28 earlier registered regressions pass, including the
complete M3 (557.97 s), M4 (2050.14 s) and M5 (1416.87 s) qualifications.
The suites run concurrently within the declared 12-worker reservations; their
elapsed times are measurements under that allocation, not performance claims.

The final orientation correction is verified at
`8a80594a47eb34edbb76a2fe827a200227b3d34b`. Four constructors are accepted and
two rejected according to independently computed physical wall-coordinate
progression. The three complete 2D/3D native operator reports are replayed,
including fields, weak rows, constraints, histories, directional actions and
budgets. All 16,212 orientation/replay checks pass.

The final audit mechanically verifies that the only production change after
the main qualification is the incoming-orientation validation. The optional
constructor-only probe and its test registration/oracle are the only other
numerical input changes. No residual, profile, time integrator, condensation
or Fixed-mode equation path changes. The previous 28 regression cases do not
call the changed Upstream guard. Their passing evidence is retained at its
actual revision; operator replay qualifies the admitted M6 cases after the
guard correction. This is an explicit source transition, not a claim that all
31 tests ran in one CTest invocation.

Final registered coverage is 31: 28 earlier CTests, the registered M6 kernel
CTest, direct execution of the registered complete M6 qualification entry
point, and the registered orientation entry point with the additional three
operator replays. The integration manifest records each execution separately.
Seventeen actual stage reports have maximum scaled residual 1.98e-12 and
maximum energy closure error 9.60e-16. Seven physical output sets independently
agree in Cartesian coordinates within 1.78e-15 and field values within 3.34e-16.
Both frozen M6 equation owners, `scalar_forcing_response` and
`transition_and_adapted_baseline`, are delivered by these receipts.

The final guard fingerprint covers 3,892 numerical/configuration inputs and
98 runtime artifacts. The final documentation-only revision must match that
fingerprint without numerical-source or runtime exceptions before sync closes.

## Evidence and reproduction

The [curated evidence](cls-epawf-revised-evidence/m6/README.md) includes the
qualification, regression/kernel XML, quadrature and orientation receipts,
source/runtime fingerprints, integrated equation/coverage/stage/output audit,
and a hash index of 445 raw artifacts. Large raw data remain on Linux.
The [PNG](cls-epawf-revised-evidence/m6/figures/transition-verification.png) and
[PDF](cls-epawf-revised-evidence/m6/figures/transition-verification.pdf) show
measured p/h/time behavior and the one-way mask/baseline buffer. They were
rendered with Matplotlib on Linux and visually checked.

Activate the installed MPI environment, use the isolated build and set all
numerical threads to one. Exact executable paths/arguments, MPI bindings,
timeouts, source revisions and wall times are in the receipts. Entry points:

```bash
python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m6.py --build /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build --output /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/new-m6-qualification
ctest --test-dir /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build --parallel 1 -R '^RANSCLSEPAWFM6KernelOracle$' --output-on-failure
python3 solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m6_orientation.py --build /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/build --output /media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/new-m6-orientation
```

Use the saved runner scripts for source synchronization, complete build and
the reserved regression allocation. The shared HIT source remains clean at
`f4b70baec10016dd098005bc5a5507b87e4a57bf`; its build/installation are preserved.
No push or physical calibration was part of M6. M7 and subsequent milestones
remain open.

## Retained diagnostics and limitations

The first strict build stopped on two misleading-indentation warnings in new
diagnostics; both were corrected. The first manufactured trajectory was
rejected. Direct expression sampling isolated an installed interpreter defect:
a binary function with one constant operand does not load that operand before
evaluation. Expressions such as max(0,variable) therefore returned incorrect
values. Mathematically equivalent dynamic zero*t operands avoid that path;
the actual evaluated sources now independently agree within 3.5e-18 after
scaling. The framework interpreter and original model laws were not changed.

A subsequent curved audit exposed an 8.94e-8 storage-accounting defect from
using constitutive quadrature instead of the stage mass matrix. Independent
weak rows already agreed within 7e-15. The diagnostic was corrected and the
original acceptance tolerances retained. Failed logs and receipts remain on
Linux. A diagnostic-generated Python cache was preserved outside the isolated
source when the synchronization guard detected it; neither checkout's source
work was overwritten. A package-provenance helper was adjusted for Ubuntu
packages that omit their installed-file manifest.

Final review identified an oblique, doubly tilted 3D plane edge case: positive
direction coefficient X need not imply positive physical progression along
q1 when the second orthonormal tangent also has an x component. The guard is
corrected to check the actual Cartesian incoming flow. Targeted plane
accept/reject tests and replay of all three qualified operators verify this
initialization-only change; model equations and time integration are unchanged.

This is a bounded native CG reference backend with replicated global
coefficients and a dense coupled solve limited to 3,000 unknowns. General
imported curved Cartesian meshes, scalable sparse condensation and arbitrary
3D terrain masks are not qualified by M6. The physical stability and general
inner-k positivity failures remain unchanged. No engineering five-flow or
physical hill-accuracy claim follows from these short manufactured runs.
