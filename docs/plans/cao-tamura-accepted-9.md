# CT-Accepted-9: pressure consistency and a serial scalable one-sweep backend

CT9 completes the independent backend implementation/equivalence/capacity work.
The one permitted finer dense history still fails absolute-pressure timestep
uncertainty. **Full-domain physical preservation is not admitted.** Backend
success does not replace that failed physical gate.

The reviewed and actual starting revision was
`0c28610e685aef95f7e6f21b9a4fdcd305dbf319`, clean in the dedicated Mac
implementation worktree. The initial remote check confirmed that head; there
were no intervening user changes. Existing main/CLS/WJ worktrees, installations,
sealed CT4–CT8 records and native histories were preserved.

The [prospective contract](cao-tamura-accepted-9-evidence/contract.json) separates
pressure consistency, backend qualification and physical preservation. The
[primitive decision](cao-tamura-accepted-9-evidence/qualified-v2/decision.json)
and [publication scope](cao-tamura-accepted-9-evidence/publication-v3/source-scope.json)
identify the exact tested data and revisions. Old successful and failed runtimes
remain attached to their own evidence.

## Outcomes

| Gate | Result and scope |
|---|---|
| Central physical reference | CT8 mathematical feasibility retained; the eight failed extraction corners and observation classifications are unchanged. |
| Initial pressure reconstruction | Pass, using all actual free pressure tests, consistent mass, frozen blocks and lifting/history. |
| New dense `.001/.0005`, T=.008 | **Fail: temporal:Pi only.** Domain, field accuracy, traction, mass and complete physical budgets pass. |
| Sparse operator and histories | Pass on d=0/displaced A95 P3/P1 and P4/P2 fixtures, with strict backend/component identity and rollback. |
| Sparse physical pair | Pass dense/scalable equivalence through T=.008; reproduces the old dense `.002/.001` pressure failure. |
| Flat capacity | Pass actual 8,864, 33,267 and 60,907 stored-slot zero-time cases. |
| Resource compliance | Pass closed numerical and publication receipts; no historical resource failure is reclassified. |
| Full-domain preservation | Unattempted: failed pressure prerequisite; no target timestep campaign assigned. |
| Terrain | Unqualified: mapped numbering/operators, joined-map and historical curved temporal/spatial questions remain open. |

## Frozen physics and environment

The physical definition remains Natural A95, the original high-Re closure,
`SamePhysicalWallDisplacedAcceptedV1`, `DisplacedNeutralCompositeV1`, physical
floor0, H=7.5mm, displacement3.5mm and the physical-height three-point fit.
The full-precision CT8 manifest remains authoritative. K0=q*u_tau², surface
h_s is fixed, and the actual evolving endpoint K enters the governing AWF.
Inner viscosity is physical and finite; outer molecular diffusion remains zero.
Nonconstant inner k, dependent epsilon, original signed production, full 3D
covariance, full-height inlet, prescribed Richards–Norris top fluxes, nonlinear
absolute outlet traction and one-sweep stage order are unchanged.

No pressure reset/regauging, reference cancellation, restoration, clipping,
coefficient/input tuning, post-k trace repair or extra governing sweep was added.
The model/profile kernels, specification files and CT8 manifest have unchanged
hashes in the publication source audit.

Linux is authoritative. The isolated root is
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-9-20261114`;
`/dev/sda2` was verified mounted there with fuseblk. Old dependency/history paths
and their existing relocation bridge were retained. Intermittent IPv4 SSH loss
was worked around with the same Tailscale peer's IPv6 address, verifying the
existing trusted host key. A shared authenticated connection avoided repeated
new connections. No numerical family was duplicated after a lost connection.

## A. Bounded dense pressure study

Exactly one fresh finer history used the sealed CT8 RT3 governing executable:
dt=.0005,16 steps,T=.008. Actual initial coefficient bytes, effective lifting,
reference and live boundary/load hashes match the old `.001` member. The first
collection's float-equality test omitted temporal fields. Its incomplete receipt
is retained; `pressure-comparison-v2` reaggregates the unchanged histories with
the declared1e-13 endpoint tolerance.

The old `.002/.001` pressure spread remains
4.804457393280823e-6 >8.411880056796802e-7. The new pair gives:

| Observable | Scaled weighted RMS timestep spread | Unchanged limit |
|---|---:|---:|
| U | 6.55838468e-8 | 2.80274272e-7 |
| w | 1.02514685e-8 | 1e-7 |
| k | 4.57692207e-9 | 1e-7 |
| Outer total epsilon | 2.17918259e-8 | 1.49210684e-7 |
| Absolute Pi | **1.67715161e-6** | **9.14165462e-7** |

The Pi physical RMS spread is9.250709748368297e-9; its scale is
u_tau²=.005515726603544432. Its small absolute magnitude does not replace the
failed uncertainty criterion. Both traction uncertainty measures and signed
mean traction pass. No additional halving or T=.02 continuation was launched.
All accepted-state field/domain/traction inventories and the endpoint four
physical-budget controls are retained in
[pressure-pair.json](cao-tamura-accepted-9-evidence/qualified-v2/pressure-pair.json).

### Initial weak balance and frozen coupling

The independent reconstruction uses the actual first frozen A,G,D,C blocks,
consistent free M, r_m/r_c, all free pressure IDs, input and history/lifting.
C is exactly zero; the unused gauge is excluded. It reproduces

`deltaPi = [D (alpha M)^-1 G - C]^-1 [-D A^-1 r_m + r_c]`.

For the `.002` first stage:

| Primitive | Measured value |
|---|---:|
| Maximum initial weak-continuity load | 2.1752556577563494e-17 |
| Continuity-driven pressure RMS | 5.966832365066653e-16 |
| Momentum-driven pressure RMS | 8.778447476731445e-8 |
| Maximum reconstructed pressure error / u_tau² | 6.5453565206926284e-18 |
| Exact frozen coupled pressure RMS | 9.001063009391185e-8 |
| Split-minus-exact pressure RMS | 1.4339262870818112e-8 |

Thus significant alpha-amplification of initial continuity error is not the
observed startup mechanism. The compatible-history/C=0 alpha-limit diagnostic
has momentum pressure RMS1.1590856036462652e-7 and negligible continuity
contribution. Exact frozen coupled correction also produces nonzero pressure:
spatial weak momentum imbalance and finite-step splitting are distinct issues.
These extra diagnostic solves never enter the governing history.

At identical initial and `.001` endpoint coefficients, normal/tangential/
mass/volume/interface orders32/8/12/12/64,48/12/12/16/64 and64/16/12/20/64
were compared. The largest diagnostic pressure-image change / u_tau² is
3.590308673615297e-12, far below the observed timestep spread. This does not
support governing quadrature as the dominant cause on those states. The
pressure observation rule separately integrates affine P2 pressure squared
exactly, using six-point physical Gauss panels.

Whole/inner/outer signed means, variation and common-time vectors are exported.
At T=.008 successive whole-domain pressure differences have alignment.954892;
||a-2b||_L2=2.32605e-8 versus||b||_L2=2.06852e-8. The inner signed means are
2.41242e-7 and7.57055e-8, much larger than their whole-domain means. The regional
pattern merits an interface/normal-stress load investigation. These are
absolute-pressure diagnostics, not mean-subtracted admission or proof of an
asymptotic order. See
[pressure-diagnosis.json](cao-tamura-accepted-9-evidence/qualified-v2/pressure-diagnosis.json)
and [pressure-vectors.json](cao-tamura-accepted-9-evidence/qualified-v2/pressure-vectors.json).

## B. Serial scalable backend

`AcceptedSerialSparseV1` is distinct from the physical descriptor. The qualified
component is `CSR1:EigenLU-LDLT-IR3-F40M:FGMRES80x1000`. Native checkpoints and
preservation/replay identities bind it. Cross-backend/component histories are
rejected; no signature was edited to resume an old trajectory.

Scope is serial affine flat2D, Original Natural A95, displaced Accepted/RN,
isotropic P3/P1 or P4/P2, IMEX1 single sweep, cacheOff and no bound/action/
nonlinear route. An explicit d=0 verification flag is limited to1000 stored
slots. Unsupported mapped/MPI/3D/research actions remain guarded. The old dense
backend and its3000-slot guard are unchanged.

Direct sorted-row assembly finalizes to CSR. Signed maps, startup mass,
general consistent projection, frozen flow/scalars, free restrictions,
mass-dual norms and required physical inventories use sparse storage. Local
radial KKT systems remain dense with dimension at most2048. Optional dense
exports/per-element global-vector localization are bounded to small fixtures;
larger packets contain CSR/actions/vectors/hashes. No global dense staging,
M-inverse-G columns or dense Schur matrix is constructed.

Eigen3.4.0 supplies SparseLU for predictor/scalars/preconditioning and positive
LDLT for consistent mass. Installed Nektar supplies flexible right-preconditioned
GMRES, restart80,max1000. The actual governing pressure action remains

`p -> D[(alpha M_free)^-1(Gp)] - Cp`.

`D diag(M_free)^-1 G/alpha-C` is solely the preconditioner. The contract's
unscaled shorthand has the same scalar-factor role when C=0; this records the
actual implementation. Every direct solve and final pressure solve recomputes
its true unscaled residual, requiring
`max|r| <= 2e-12 max|rhs| + 1e-25`. RHS normalization avoids the installed GMRES
small-RHS stopping floor. Diagonal mass never replaces governing consistent M.

The [allocation/method audit](cao-tamura-accepted-9-evidence/backend-method.md)
covers startup, projection, assembly, pressure, norms, diagnostics and snapshots.
The route caps each CSR matrix at20M entries, reported factor fill at40M, stored
slots at70000 and process address space at16GiB. IMEX1 requests are bounded to64
steps per segment, so retained states and attempted-stage records remain O(N)
with a fixed segment bound. Longer campaigns require genuine same-method restart.
Sparse fill can grow substantially; this is a measured bounded prototype, not
an asymptotic complexity theorem or a distributed backend.

### Preserved failures and narrowly justified changes

`capacity-v1` retains the intermediate pressure-preconditioner true-residual
failure and target P4's original30M factor cap failure. At identical coefficients,
one factor-reused residual correction reduces the intermediate288-row
preconditioner defect from1.28184e-12 to2.01658e-13, below the unchanged
7.49159e-13 bound. The IR3 component permits at most three such linear residual
corrections. These are not additional governing sweeps.

The P4 recheck measures38,165,601 LU entries in the29,636-row momentum block
with2,081,794 original nonzeros. A separately declared40M resource cap gives
about4.8% measured fill headroom under the unchanged16GiB ceiling. It changes
neither a physical tolerance nor the dense guard. Fresh component-bound controls,
histories and capacity runs follow both amendments. The final capacity ledger
records exactly one residual-refinement event, in the intermediate preconditioner.

Failed builds v1/v2/v6, the initial control-collector key error and incomplete
first pressure export remain archived. The first publication audit checked the
wrong SDK file for PETSc; the installed build cache explicitly records OFF.
That metadata failure is retained separately from numerical qualification.

### Equivalence and regression

All four tiny operator cases pass. Worst complete row/mass/frozen-update
comparison is5.21805e-15, below2e-10; this raw bound also implies the established
max(1,expected-magnitude)-scaled bound. Independent Nektar all-mode value/
derivative/face/projector checks and complete pressure-action/rank tests pass.
The zero-epsilon-flux, shear-sign, mass-replacement, rough-pressure omission and
free-mask controls are discriminating. Absolute traction shift.125 gives
pressure shift.125 with maximum coefficient defect2.44182e-12.

Same-backend uninterrupted/restart histories match exactly; rejected-step rollback
is exact. Default dense tiny histories remain bitwise equal to CT8. Wrong model,
fit/boundary identity, cache, method, backend and old sparse-component histories
are rejected for their intended reasons. Both focused CTests pass, including
intentionally failed required-check process-exit controls.

Fresh sparse physical `.002/.001` histories reach T=.008 and retain complete
four-control/four-equation budgets, k/epsilon reconstruction, domain/constitutive
checks and zero physical corrections/extra governing sweeps:

| dt | Maximum scaled field difference from dense, all common times | Maximum scaled ground/interface traction difference |
|---|---:|---:|
| .002 | 8.79622e-14 | 3.31562e-13 |
| .001 | 1.62990e-12 | 8.59902e-13 |

Both are well below the1e-8 trajectory budget. The sparse pair reproduces the
old dense failed `temporal:Pi` classification. This proves cross-backend
numerical equivalence on the stated cohort, not physical preservation.

## Target capacity and measured cost

All three fixed flat layouts pass fresh initialization/domain/projection binding
and complete frozen mass, predictor, pressure and scalar operations. Counts
include constrained coefficients and the unused fixed-zero gauge slot.

| Layout | Stored slots | Flow rows / nnz | Pressure iterations | Startup wall / CPU s | Frozen wall / CPU s | Peak RSS KiB |
|---|---:|---:|---:|---:|---:|---:|
| Intermediate P3/P1 | 8864 | 4584 /247825 | 42 | 36.18 /35.52 | 36.94 /36.27 | 753288 |
| Target P3/P1 | 33267 | 17671 /993162 | 66 | 133.58 /132.82 | 146.50 /145.72 | 2656656 |
| Target P4/P2 | 60907 | 33542 /2814858 | 150 | 144.87 /144.09 | 179.43 /178.67 | 3662944 |

Each target has930 elements. P4 predictor/mass/preconditioner factor counts are
38,165,601 /1,129,586 /1,778,418. Its true pressure/continuity defects are
2.84157e-22 /2.84206e-22, and independent complete-native G action defect is
7.14684e-22. Matrix, RHS and direct/mass residual primitives are retained in
[capacity.json](cao-tamura-accepted-9-evidence/qualified-v2/capacity.json).
Per-case output is7,989,754 /20,037,120 /29,309,775 bytes. These are static
costs, not measured target accepted-step times.

Native-quadrature representation errors against the same continuous reference:

| Quantity | Target P3 RMS / gradient RMS | Target P4 RMS / gradient RMS |
|---|---:|---:|
| U | 5.95814e-6 /1.27491e-3 | 7.32900e-7 /1.88154e-4 |
| k | 6.60630e-7 /8.94181e-5 | 2.16313e-7 /4.36326e-5 |
| Outer epsilon | 6.14378e-8 /1.76829e-5 | 6.17598e-9 /2.13532e-6 |

w and reference Pi representation errors are zero. RMS is weighted physical
RMS; true L2 is retained separately, with full volume1100 and outer volume
1089.6875. Gradients use the retained H_h coordinates. This is initial
representation characterization on native mass rows, not evolved spatial
convergence. The target retains the actual wall cutoff/roughness integration
breaks even where they are not mesh nodes. Its detailed integration orders and
breakpoint audit are in the case bodies. Additional target quadrature/traction
and temporal admission remain prerequisites to physical advancement.

On the actual accepted strip `.001` steps, native exclusive CPU ranges
6.07–7.02s/step. Across that run, domain evaluation takes26.77s, flow assembly
10.99s, residual work4.98s and complete/physical audits4.47/2.58s. Target P4's
**detailed-instrumentation** frozen preflight spends38.48s in element-row work,
21.57s in basis lookup/evaluation and13.14s in evaluation-row construction.
These are exclusive timers; nested inclusive totals are never added. Timing
instrumentation has overhead, and static costs are not extrapolated into a
claimed target trajectory. Immutable geometry-row reuse is a later measured
engineering opportunity; mandatory state/gradient-dependent domain work remains.

## Full-domain preservation and remaining admission

The matched flat/terrain definition retains x=[-25,30],roof20,H=.1875 and the
explicit CT8 radial grid. Flat counts are measured; the current affine
height classifier does not provide an admitted terrain field count.

Independent inner resistance integration plus the outer logarithm gives bulk
reference velocity1.3539582080070292 and transit time40.621637857609905. The
prospective.01/.1/1-transit endpoints are.40621638 /4.06216379 /40.62163786.
**None was launched.** Pressure uncertainty fails, target evolved-reference and
timestep/quadrature checks remain, and a finite accepted-run resource plan must
be declared from relevant measured costs. The dense3000 guard was never raised.

Physical experimental applicability is unchanged: roughness-sublayer/interface
adequacy, the exceeded historical.1 thin-layer screen, smooth preparation and
complete approach-layer fidelity remain open. CT6/CT7 curved temporal/degree
failures and all eight historical unresolved epsilon-order observations remain
unchanged. Sparse flat equivalence does not close joined-map Cartesian gradients,
cofactors, transmission, region-aware numbering or variable-height evolution.
EP, MPI/3D, physical terrain production and hill validation remain outside scope.

## Reproduction and provenance

The tested C++ implementation revision is
`1e73bf53353490e7117f36799250ef79208e8df6`. Later commits change drivers,
publication and documentation. The exact controller revisions are recorded in
requests/publication scope. Linux uses GCC11.4.0, Release `-O3 -DNDEBUG`, C++17,
Eigen3.4.0-2ubuntu2, SDK MPI ON with one rank, PETSc OFF and one-thread numerical
libraries. At most two compiler jobs and one numerical family ran at once.

The frozen sparse runtime is `build-v8/runtime`. SHA256:

- IncNavierStokesSolver: `5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7`.
- RANSPreservationProbe: `3038ff84241d493c069baeea19037b6a484956ae05498529ebc689e10db4ac33`.
- The sole new dense physical history uses CT8 RT3 solver
  `a95fc6037946ce088185337fbf92d3d143a1e12b0dfa2d62ed591bf2f9fcdfa1`.

Exact compiler flags, header/library hashes, actual mapped-file catalogs,
commands, configuration/state/history hashes, high-water PIDs/threads and
suspension events are in the compact evidence and closed publication records.
Peak listed numerical resources are15 project /18 aggregate workers, within
16/72. Publisher/collector closure is bound after process exit; elapsed time,
process CPU and suspension intervals retain separate meanings. The capacity
family uses698.34 elapsed seconds and979.81 process CPU seconds including its
control/monitor/postprocessing work. Historical resource failures are untouched.

[Reproduction commands and archive conventions](cao-tamura-accepted-9-evidence/README.md)
cover each immutable destination. Complete matrices, fields, histories and first
failure witnesses stay in the isolated Linux archive; compact initial arrays
are represented by counts/hashes. Neither compact receipts nor diagnostic
increments are restart inputs.

The smallest next package is a **single localized strip pressure/representation
contrast**, selected using signed initial weak-momentum loads and their pressure
image. Preserve absolute pressure, the same physical inputs, one-sweep method
and existing uncertainty limits. A
[ready-to-use prompt](cao-tamura-accepted-9-evidence/next-package-prompt.md)
predeclares that bounded study. Full-domain preservation may follow only after
its physical gate and the target-specific preflights pass.
