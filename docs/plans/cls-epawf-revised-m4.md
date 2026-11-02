**M4 closure for the revised CLS/EP-AWF reference and native timestep**

M4 is closed on 2 November 2026 on codex/cls-epawf-revised. The tested Linux
runtime revision is 722d55c5d3c2055b9c0c4e878931350b92d47350. Its complete M4
qualification and the other 25 registered regressions pass, covering all 26
registered tests without a duplicate full M4 solve. M5 and later milestones
remain open.

This milestone implements the constrained flat composite target and its actual
native IMEX evolution. Its authority is the [revised note](/Users/jay/Documents/nektar-v2026/docs/literature/cls_neutral_abl_equilibrium_awf_note_revised.tex), especially
`eq:hp-discrete-target`, `eq:hp-constraints`, `eq:hp-defect-correction` and
`eq:hp-time-fixedpoint`, together with the unchanged
[approved M0 contract](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-contract-v1.json).

**Material Passport**

- Origin Skill: academic-research-suite / experiment-agent.
- Origin Mode: implementation and native computational verification.
- Origin Date: 2026-11-02.
- Verification Status: VERIFIED for the constrained flat reference, native
  IMEX map, bounded 2D/full-3D MPI backend and complete history restart. Physical
  stability and general inner-k positivity remain FAILED AS SPECIFIED.
- Version Label: cls-epawf-revised-m4-v1.

**Implementation and scope**

`RANSEquilibriumWallEvolution` is a new equation-system selector. The M3
`RANSEquilibriumWallSystem` remains its original serial single-stage path.
The new path creates native Modified_A CG fields on affine, axis-aligned,
canonically ordered quadrilaterals or hexahedra. It stores full-domain velocity
and k, outer-only epsilon, two regional pressure spaces and one common physical
volume-mean pressure gauge. Velocity/scalar degree p is at least three;
pressure degree is p−2 in the qualified fixtures.

The MPI implementation replicates globally numbered native coefficient spaces.
Each physical cell owns volume/source quadrature on `geometry ID % rank_count`.
Each rank assembles its owned loads, reductions form the global residual, and
rank zero solves the dense coupled system and broadcasts the update. Physical
mass and gauge terms enter once. Nonlinear solves are bounded to 3,000 coupled
unknowns; larger reference audits are permitted. This is a bounded verification
backend, with no claim of scalable distributed sparse solution.

The original full XML input is recovered from Nektar's positional command-line
arguments. The loader's current filenames refer to partition files after MPI
mesh creation and cannot be reused as a global mesh. Rank-local full-mesh
sessions preserve effective parameter values and use native DirectFull maps;
the world input uses the parallel-compatible IterativeStaticCond entry setting.
All generated mesh partitions, cloned sessions, fields and checkpoints stay in
the experiment directory, outside the source checkout.

The implementation comprises `RANSWallGlobalField`, `RANSWallEvolutionPatch`,
`RANSWallTimeIntegration` and `RANSEquilibriumWallEvolution`. Native time and
spatial probes accompany independent physical Jacobi/mass/constitutive oracles.
The [user guide](/Users/jay/Documents/nektar-v2026/solvers/IncNavierStokesSolver/Examples/RevisedEquilibriumWall/README.md)
describes session selection, output and restart.

**Declared target and anchor**

The default neutral case fixes u_star=0.4, H=2, Z=30, length 300,
roughness 0.1, r_K=r_m=1 and nu=1.5e-5. Its two-dimensional coarse fixture has
two horizontal cells and six radial cells, p=3/p=1 and 515 coupled unknowns.
The full 3D reference adds span 20 and has 2,541 unknowns.

The velocity is a constrained physical L2 projection of the original AWF below
H and its matched logarithm above H. Constant k and horizontal invariance are
represented exactly. The projection preserves ground velocity zero, U(H), the
one-sided reference velocity gradients at H and compatible top value/gradient.
The outer A/z epsilon projection preserves its dependent H trace and compatible
top value/gradient. The actual stored polynomial, differentiated by native
operators, satisfies these trace/gradient constraints. No evaluated derivative
is replaced by an analytic value. The two pressures retain the required
pi_w−pi_o=−11*K0/75 and one common volume-mean gauge.

Independent one-dimensional Jacobi KKT assembly reconstructs the constrained
L2 minimizer. Independent native-map evaluation also checks exact horizontal
invariance, discrete divergence, ground no-slip, matching velocity/k/epsilon,
full vector traction, normal pressure offset, zero k flux and top epsilon flux.

Chi=1 selects the physical EP-AWF source. Eta=0 is the default unanchored
operator and does not disable EP-AWF. Its polynomial target generally has
nonzero free differential residuals. The qualification evaluates eta=0 first,
records every block, and studies polynomial refinement of all regions and
inner h refinement at fixed outer mesh and physical wall parameters. Separate
unanchored native one/repeated steps demonstrate the resulting target drift.

For isolated neutral verification only, eta=1 stores one reference and subtracts
its fixed free momentum/k/outer-epsilon loads. Divergence, gauge, essential and
dependent-field constraints are checked first and remain unmodified. The
complete unmasked reference loads, removed rows, free defects, constant-test
integrals and energy work are retained. The primary two-horizontal-cell case
has nonzero free constant-test loads and energy work; the audit is not made
vacuous by removing every constant-bearing coefficient as a boundary row.

**Native stage balance**

The unchanged native IMEX implementation supplies first-order startup and
orders 1/2 with an empty variant. Differential coefficients are flattened for
that interface, while each field keeps its distinct native space, physical mass
matrix, domain volume and constraints. Pressure and the gauge remain algebraic
coupled stage unknowns. There is no separate velocity-correction pressure path.

Write B(t) for the fixed physical prescribed momentum load. Its inner
contribution is the projection of b_i−1_(r<=h)*b_o,t(H) in tangential momentum
and b_i in normal momentum; its outer contribution is b_o. The same owned
physical load and the same quadrature are used by the full residual and by
`ExplicitBody=M_ff^−1*B_f`. In particular, prescribed load integration does not
change when current sampled K moves the adaptive wall quadrature. The
state-dependent part of F_A, transport, original full CLS stress, production,
reaction, current-K profiles and dependent epsilon trace stay implicit.
Manufactured scalar source functions are evaluated at the implicit stage time.

At a stage with native history h and coefficient lambda, the free rows solve

    M*(Q_stage−h)/lambda + R_implicit(Q_stage,t_stage) − eta*d = 0,

along with the original continuity, gauge and boundary/trace constraints.
For steady reference data, d=R_implicit(Q0)−B. Consistent native histories give
h=Q0+lambda*M_ff^−1*B on free rows. Substitution cancels both the spatial
reference balance and B, including the same coupled pressure/constraint path.
The zero-body reference has h=Q0 and zero explicit derivative histories.
Nonzero prescribed-body controls separately exercise this stage balance.

Native first-order startup generates valid second-order histories before use.
Unused IMEX2 master slots at zero steps retain their native zero values; no
synthetic backward-time state is inserted. The native recurrence/history probe
checks those slots, the startup conversion, every subsequent history slot and
checkpoint continuation at steps 0/1/2/3.

Every accepted coupled stage is checked using its actual full scaled physical
mass-dual residual. The fixtures explicitly select 2e-12, stricter than the
frozen 1e-9 reference residual/state/constraint gate; isolated interface
identities use the frozen 1e-10 tolerance. A previous state is only
accepted when it already meets that stage equation. No free-field reset,
reference projection, turbulence floor, production clipping or stabilization
is applied after the native step.

**Qualification coverage**

- Independent 2D/full-3D native basis and gradient evaluation; complete free
  residual blocks at chi=0, 0.37 and 1; one pressure null mode; and essential
  Jacobian-column controls.
- Independent physical mass verification of the exact explicit/implicit body
  partition for nonpolynomial prescribed data and distinct sampled-K states.
- Constrained target, pressure/traction/flux and gauge checks; raw eta=0 defects,
  h/p studies, removed-row audit, free-load integrals and anchor energy work.
- Unanchored native one/repeated steps, anchored IMEX1/2 at three dt values,
  100 steps over H/u_star=5 seconds, and stationary nonzero explicit-body controls.
- Full 3D straight and horizontally rotated flow, direction (0.6,0.8), on
  1/2/4 MPI ranks. Tilted and curved walls remain M5 work.
- Independently manufactured 2D unsteady solution with inner storage and
  normal transport, outer k/epsilon transport and active Yap source; observed
  differential-field temporal orders from dt=0.1, 0.05 and 0.025.
- Full/condensed and genuinely evolving MPI 1/2/4 agreement, with actual root
  coupled solves and retained inner k/pressure/history.
- Complete checkpoint continuation through startup and established histories,
  including serial-to-MPI4 and MPI4-to-serial continuation; incompatible source,
  dt and eta controls plus corruption/shape/publication controls.

The native checkpoint contains all differential value and explicit-derivative
history slots and time labels, the full coupled state including pressure/gauge,
fixed reference, fixed defect and eta. A mesh/basis/physics/function signature
and checksum protect the canonical payload. Publishing uses an atomic final
file in a unique directory and never overwrites an earlier generation. Restart
input retains the original Time epoch; payload time and steps determine the
continuation. NumSteps requests additional steps.

**Retained development findings**

The Linux run directory retains failed attempts, logs, raw states and receipts.
Strict compilation caught a diagnostic vector allocation warning, fixed using
explicit fixed-size row storage. The native StdQuadExp constrained transform
reverses its top/left physical edges against the current increasing-order trace
map. The wall adapter uses two commuting native segment projections, preserving
all edges/corners and the interior-bubble L2 equations; the framework routine
and original CLS/AWF kernels are unchanged.

The first 3D absolute-state IMEX2 startup transition exposed an inherited M3
1e-13 tolerance. Its 1.62e-13 physical scaled residual triggered unnecessary
large Jacobians. M4 fixtures now explicitly use the intended 2e-12 physical
stopping rule. Exact lifted boundary columns are identity columns; only free
columns need full coupled differentiation, including the dependent epsilon
trace's variation with current free k.

MPI development also identified the need for a parallel-compatible world
solver selection and the original full input filenames after mesh partitioning.
A separate nonpolynomial body-load audit identified and removed a quadrature
mismatch between the full residual and native explicit derivative. The final
partition uses one fixed physical prescribed load in both paths.

An exploratory unanchored IMEX2 run at dt=0.001 reached a 2.62e-12 residual and
rejected further Newton trials against a 2e-12 stopping rule. This failure is
retained rather than described as a physical instability. The qualified
unanchored baseline uses IMEX1 at dt=0.01/0.005; the anchored and unsteady IMEX2
dt studies use 0.1/0.05/0.025. No claim is made for arbitrary dt or state ranges.

**Execution and results**

The full M4 qualification passes on Linux revision
`722d55c5d3c2055b9c0c4e878931350b92d47350`, with 12,863 recorded checks
and its final source-immutability check (12,864 in the command output).
The independent native time probe passes 6,809 checks and rejects 24 invalid
controls. These counts are reported separately. The isolated-interface recheck
passes 8,272 checks over all 28 neutral cases; its largest traction error is
1.125e-12 against the frozen 1e-10 tolerance. The all-slot history recheck uses
physical field masses and scales, including dt*field_scale/(H/u_star) for native
dt*RHS slots. All 216 comparisons pass, with a maximum error of 2.905e-15.

The coarse eta=0 differential norms (U,W,k,epsilon) are
(5.54367, 0.0163108, 0.930097, 0.00468988). The stored free-load constant-test
integrals are (9.02463, 0.0118606, −33.95895, −0.00299031). These integrals use
the declared zero-constrained-row extension of the dual load; the unmasked
removed-row reactions are audited separately. At eta=1, the coarse zero-body
anchor energy work is +0.0254621. The nonzero prescribed-body verification
control has work −26.3425. The anchor is explicitly changing the discrete
balance and its energy contribution is not hidden.

| p | Coupled coefficients | U defect | W defect | k defect | epsilon defect | Velocity profile RMS | Epsilon profile RMS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 | 625 | 0.440874 | 0.00620732 | 0.0145008 | 0.00236872 | 0.000167797 | 4.90336e-06 |
| 5 | 1641 | 0.0706874 | 0.00019748 | 0.0014219 | 0.000218344 | 4.93858e-06 | 6.40775e-08 |
| 7 | 3153 | 0.00767394 | 5.25196e-06 | 9.75963e-05 | 9.48866e-06 | 1.84227e-07 | 8.58662e-10 |

The physical inputs and mesh are fixed in this p study. The p=7 case is an
audit with 3,153 coefficients and does not invoke the bounded nonlinear solve.
The unanchored polynomial target has measurable defects even at p=7; it does
not meet the exact stationary-reference gate. The native eta=1 map is qualified
separately. Inner h refinement reduces the dominant U/k defects and velocity
projection error. The outer mesh remains fixed in that h study, so its epsilon
defect stays fixed; the W dual norm also reflects the changing full test space.
All-region p refinement reduces both outer momentum and epsilon defects.

![Unanchored reference refinement](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/figures/reference-refinement.png)

All anchored zero-body tests, including IMEX1/2 at dt=0.1/0.05/0.025, the
100-step H/u_star run, full 3D and rotation (0.6,0.8), have exactly zero reported
state drift. Their accepted stage balances remain below the frozen limit.
The nonzero explicit-body controls also have zero state drift. Unanchored
IMEX1 instead changes the stored target: after one dt=0.01 step its scaled
U/k errors are 7.697e-3 and 1.533e-3. All unanchored accepted stages satisfy
the same residual requirement; the largest is 1.130e-12. No stationary claim
is made for that unanchored coarse target.

| Method | U observed orders | W observed orders | k observed orders | Epsilon observed orders |
| --- | --- | --- | --- | --- |
| IMEX1 | 0.9233, 0.9612 | 0.9471, 0.9726 | 0.9475, 0.9728 | 0.9524, 0.9754 |
| IMEX2 | 2.0004, 2.0125 | 1.9885, 2.0031 | 2.3693, 2.1667 | 2.1439, 2.0699 |

These are physical mass-scaled errors against independent strong manufactured
fields from t=0.3 to 0.7, with dt=0.1/0.05/0.025. Pressure errors and their
observed rates are retained in the JSON in one common gauge; the table reports
the differential-field temporal orders.

![Native temporal orders](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/figures/temporal-order.png)

The largest full/condensed final-state discrepancy is 1.698e-16. The evolving
MPI 2/4 discrepancy is at most 2.360e-16; those cases perform three or four
actual nonlinear iterations per stage. Full-history restart discrepancies are
at most 2.383e-15, including MPI4-to-serial continuation. Native slots, stored
reference/defect, pressure and gauge are compared as well as the final
current state. The frozen MPI/restart allowance is 1e-8.

The M4 oracle is one of 26 registered tests. It is executed once directly,
with identical source, binaries and test arguments apart from its output
directory. The remaining 25 registered tests pass in 724.24 seconds at the same source
and binary, including all 8,085 M3 checks. The final integration audit records
both commands and their complete registered coverage. All 37 resolved shared
libraries are byte-identical to M3, and the source manifest, original kernels,
revised note, frozen contract and unrelated main checkout are verified unchanged.

Across all accepted M4 stages, the largest scaled physical residual is
1.206e-12 and the largest constraint defect is 3.802e-15. The raw artifact
index contains 283 hashed input, output, log and checkpoint artifacts.

The numerical build uses the complete production incompressible-solver source
list in the isolated Release build. All new kernels/probes/solver sources use
Wall/Wextra/Wpedantic/Werror checks. Compiler parallelism is 16; tests are
sequential, with MPI 1/2/4 and OMP/OPENBLAS/MKL/BLIS threads set to one. The
72-worker limit is respected. The server checkout, shared build/installation
and unrelated HIT work remain preserved. All connections use Tailscale.

**Evidence and reproduction**

The [qualification JSON](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/verification.json),
[integration audit](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/integration.json),
[regression output](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/regression-output.txt),
[isolated-interface audit](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/isolated-interfaces.json),
[physical history audit](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/history-norms.json)
and [raw artifact index](/Users/jay/Documents/nektar-v2026/docs/plans/cls-epawf-revised-evidence/m4/raw-artifact-index.json)
retain exact commands, inputs, logs, results and hashes. Native meshes, fields,
checkpoints and larger data remain below
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised`. Failed and successful
development attempts remain there with their original receipts.

The principal numerical commands, after the source-pinned isolated Release
build, are:

```bash
source /home/jay/.local/opt/nektar-v2026-mpi/env.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1
CLS_EP_RUN=/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised
export TMPDIR="$CLS_EP_RUN/tmp"
export LD_LIBRARY_PATH="$CLS_EP_RUN/build/geometry-lib:$LD_LIBRARY_PATH"
python3 "$CLS_EP_RUN/source/solvers/IncNavierStokesSolver/Tests/RANS/validate_cls_epawf_m4.py" \
  --build "$CLS_EP_RUN/build" --output "$CLS_EP_RUN/manual-m4-qualification"
ctest --test-dir "$CLS_EP_RUN/build" --parallel 1 --output-on-failure \
  -E '^RANSCLSEPAWFM4Oracle$'
```

Use a fresh output directory. The standalone oracle is the complete registered
M4 test; running all 26 directly through CTest is also supported. Later verifier
changes tighten isolated comparisons to 1e-10 and use physical masses for all
native history slots. The archived postchecks apply precisely those stronger
norms to every corresponding tested case, without rerunning or modifying the
numerical trajectories. The final source-immutability guard also runs before
publishing a successful verifier report. Runtime C++/fixtures remain identical
to the tested revision. These reporting/verifier changes do not alter the
model, numerical solve or stored states.

**Frozen verification IDs**

| Contract ID | Evidence |
| --- | --- |
| M4.projected_composite_target | Independent constrained physical L2 minimizer, exact k/invariance, native traces/gradients and full traction/pressure/flux identities |
| M4.assembled_equilibrium | Eta=0 free defects and refinement first; separate eta=1 physical mass-dual reference residuals with unchanged constraints |
| M4.constrained_discrete_anchor | Fixed free differential loads, removed-row audit, constant-test integrals, energy work and exact prescribed-load time partition |
| M4.stationary_time_map | Actual native IMEX startup/stages/histories, repeated dt studies, nonzero-body balance, full 3D rotation, MPI and complete restart equivalence |

The contract's historical `availability: planned` entries remain frozen; this
report records delivery of these four IDs.


**Physical and implementation limits**

The original CLS adverse incremental diffusion and the inner-k positivity
counterexample remain failed physical gates. M4 preserves and verifies the
specified equations; a stationary anchored target does not establish stability
or constitute a physical repair. No constitutive coefficients, molecular-term
policy, production rule, turbulence floor or frozen contract are retuned.

The qualified implementation is flat, affine, bounded and uses compatible
external Dirichlet data. Its 3D time tests concern the stationary reference and
horizontal rotation; independent 3D residual checks exercise all components and
derivatives. Nonstationary temporal-order tests use the independently
manufactured 2D case. Curved geometry, localized transition, scalable parallel
solvers, nonlinear perturbation stability and terrain/wind-tunnel validation
remain outside M4.
