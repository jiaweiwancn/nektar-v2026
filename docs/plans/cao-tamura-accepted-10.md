# CT-Accepted-10: signed pressure attribution and one radial contrast

The single radial contrast resolves the retained strip pressure-uncertainty
gate. The refined A95 pair passes through **T=.02**. The principal initial
driver is represented **outer nonlinear normal stress**, with important signed
cancellations. The actual target's initial/operator/quadrature preflight passes,
but its required physical-scalar observer fails the prospective 16-GiB
address-space capacity test. Target evolution remains **unattempted, T=0**.

Starting/fetched remote revision: `af87414abc2b583492c3a9f3466b25cb0fe5285a`;
the dedicated worktree was clean with no intervening remote commits. Source/Git
work used the existing Mac branch worktree. Authoritative native work used a
new isolated Linux root on the verified `/dev/sda2` fuseblk mount:
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-10-20261114`.
Existing worktrees, shared installations and sealed CT4–CT9 records are preserved.

## Contract and implementation

[The prospective diagnostic contract](cao-tamura-accepted-10-evidence/contract.json)
fixes the central [CT8 physical manifest](cao-tamura-accepted-8-evidence/physical-reference-manifest.json),
Natural A95, displaced Accepted/composite policies, H=7.5mm, d=3.5mm, the
three-point fit, finite inner viscosity, high-Re outer equations, nonconstant
inner k, dependent epsilon, full signed production/covariance, RN top and
absolute nonlinear outlet traction. Surface calibration, continuous functions,
ConstrainedRadialV1 projection, IMEX1 single sweep and cacheOff are unchanged.
No governing solver was rebuilt or modified for these trajectories.

`MomentumAttributionPacket` adds opt-in read-only component vectors and sparse
generating-element row supports. It separates augmented pressure, molecular,
linear and nonlinear normal/shear stress, conservative advection, rough
pressure/inertia forcing, prescribed body load, natural traction and natural
advection. It exports both analytical and represented evaluations on identical
native owners, coordinates, weights and tests. Free rows, replaced constraints,
prescribed reactions and shared interface IDs remain distinct. One-sided
interface fluxes are diagnostic; they are not added again to the shared rows.
There is no element-count times global-size dense array. The original
1,000-slot dense localization guard remains and its rejection is tested.

A separate `export_faces_only` option exposes the existing physical face rows
without the full velocity/scalar volume catalog. Its face/owner/side/weight
arrays are bitwise equal to the sealed CT9 arrays; state, reference and boundary
data are unchanged. Native audit differences satisfy the tight operator bound.
The old scalar observer still constructs its own full internal volume packet;
this is the later target-capacity blocker.

The final `CT10SignedMomentumV2` observer explicitly records the effective-input
policy used by the governing `Residual()`, including any `ExternalLift` changes.
It does not change that policy or modify the supplied state/history. Baseline
and refined initial replays reproduce the sealed component loads, and a refined
T=.02 replay agrees with the old original residual to 8.89e-15. The three
tested raw/effective inputs have zero lift difference. Raw accepted-state
traction and refreshed epsilon-trace lag remain separate diagnostics evaluated
from stored coefficients. This supplement leaves the sealed initial attribution
and all governing trajectories unchanged.

## Signed loads and pressure images

Independent finite-nu profile integration and spectral constitutive formulas
reconstruct native component loads to 5.30e-15 or better on the baseline.
Their sum reconstructs the original unmasked momentum residual to 1.07e-14;
agreement with the actual first frozen momentum load is 7.74e-15. Modified_A
constant tests use signed vertex coefficients, not ones in every modal slot.
Whole/inner/outer volume, gradients, native sides and element supports are retained.

The diagnostic uses the actual captured free consistent mass, G and D, with
C=0 and the unused gauge excluded:
`Pi_j = -W^-1 D M^-1 r_j`, `W = D M^-1 G`.
Mass solves are sparse direct; W is applied without storing M^-1G or a dense
Schur matrix. True residuals, a tighter repeat, signed component/band sums and
cross-products are checked. Recomputed coefficient-sum uncertainties are
4.72e-10 / 5.98e-10 on the u_tau² scale, below 1e-8. Held-state 48/64 controls
also pass. These extra solves never enter a governing state or history.

For the represented-minus-continuous baseline loads, define
`a_j=<Pi_j,Pi_total>_physical_mass/||Pi_total||²`. Signed alignments can exceed
one because contributions cancel:

| Contribution | Signed alignment |
|---|---:|
| Outer outlet advection | +7.434782 |
| Outer volume advection | -7.434782 |
| Outer nonlinear normal stress | +1.218412 |
| Inner linear shear | -0.099137 |
| Outer linear shear | -0.069332 |
| Inner molecular shear | -0.049942 |

The normal-momentum image contributes +1.083027; the streamwise image opposes
it by -0.083027. The large advection images almost cancel. Ranking unsigned
stress or advection norms would misidentify the driver. The analytical total
image is near roundoff. The represented total image has whole-domain RMS
1.15909e-7 and inner RMS 7.74025e-7, illustrating why a large inner pressure
response does not locate its forcing inside the AWF.

The outer FE gradient and epsilon representation do not maintain the exact
parallel-flow constitutive relation pointwise; the resulting varying nonlinear
normal stress is not balanced by the exactly represented constant reference
pressure. This explains the dominant signed spatial load. The contrast changes
the joint velocity/epsilon representation; it does not isolate either field as
the unique cause. Finite-step splitting remains a separate CT9 finding.

## Prospective contrast and fresh histories

[The sealed selection](cao-tamura-accepted-10-evidence/single-contrast-v1.json)
captures 80.4629% of the absolute signed-alignment signal in the smallest
contiguous band: the first two outer slabs `[.1875,.28125]` and
`[.28125,.421875]`. Factor-two subdivision adds `.234375` and `.3515625`
(9.375 and 14.0625mm). Inner cells, x=[0,.125,.25], roof20, P4/P2, physical
interface, all functions, projection and governing quadrature remain fixed.
Native stored slots change **2,957 -> 3,265**, within the existing sparse guard;
the dense 3,000-slot guard is unchanged.

Static admission precedes trajectories. Both candidate timestep members have
the identical actual initial hash
`8d47e251c6e3324f246073185502bd94548a3aa2201f71d92ea89d369c711aa9`.
Baseline and candidate initial coefficients are intentionally different spaces;
their continuous physical inputs are identical. Live explicit/implicit boundary,
effective state and fixed-reference receipts are retained. Initial scaled U
reference error falls from 1.77668e-5 to 9.60916e-6; outer-epsilon error from
1.07564e-5 to 1.40641e-6. The initial momentum-driven pressure-image RMS falls
to 4.45858e-8, ratio **.384663**. Pressure initialization is unchanged.

Exactly three new primary histories were generated: sparse baseline .0005,
candidate .001 and candidate .0005. The sealed sparse baseline .001 is reused.
The fresh sparse .0005 reproduces the sealed dense .0005 within 3.27e-12 in
scaled fields and 1.45e-12 in traction, below the 1e-8 equivalence budget.

All following pressure numbers are scaled weighted physical RMS, with absolute
pressure and the unchanged `max(1e-7,.1*drift)` rule:

| Space / endpoint | dt pair | Pressure spread | Limit | Decision |
|---|---|---:|---:|---|
| Baseline, .008 | .001/.0005 | 1.677152e-6 | 9.141655e-7 | **FAIL**, retained |
| Single contrast, .008 | .001/.0005 | 2.928135e-7 | 5.026636e-7 | **PASS** |
| Same contrast, .02 | .001/.0005 | 1.513382e-7 | 5.311233e-7 | **PASS** |

At .008 the pressure spread is .02928% of the absolute .001 preservation
tolerance, versus .16772% for the baseline. This additional interpretation
does not replace either frozen relative classification. Common observations
resolve the union partition and analytical breaks, with explicit interface
sides. Value/gradient projection error, drift, absolute error, spatial
differences and regional pressure means/variation are separate in the evidence.
Two spaces establish no spatial order; timestep pairs establish no temporal order.

At T=.02 the fine member has scaled drifts U=2.1572e-6, w=3.7141e-7,
k=3.3746e-7, outer epsilon=5.1150e-7 and Pi=5.3112e-6. Ground error is
9.0551e-5; full interface-traction error is 1.16362e-4; mass defect is
2.01e-15. All field/traction temporal gates, native/independent domain and
covariance/radicand/cutoff inventories, and all four controls/four equations
plus both scalar reconstructions pass. Native scaled budget defects are about
1.4e-7 for inner k and 1.4e-6 for outer epsilon, below .02. Corrections and
extra governing solves are zero. Actual current endpoint K and gradients are
used in traction; prescribed data and row reactions are not substituted.

Continuation uses unchanged-method native histories and the original time-zero
drift baseline. A controller exact-float check stopped the fine continuation
before launch after the coarse member completed .02. A new receipt validates
and reuses that completed coarse member, applies the established 1e-13 endpoint
tolerance, and advances only the fine history. The failed preparation remains
archived. Earlier probe-build/API and array-adapter preparations are likewise
retained; no physical trajectory was rerun to erase them.

## Actual target readiness and the remaining gate

The unchanged matched target is x=[-25,30], roof20, H=.1875, P4/P2,
930 elements and **60,907 stored slots**. Its actual initialized vector, source
reference, complete free-row loads, natural faces and independent samples pass.
Scaled initial U/k/outer-epsilon errors are 9.85156e-6 / 1.03679e-5 /
6.05198e-7; Pi error is 3.02e-16. Initial traction errors are near roundoff.
The independent covariance minimum is .00611043 and A95 radicand minimum
.69165083. Raw momentum reconstruction is 8.18e-15. Governing-quadrature
pressure changes are at most 2.46e-10 on u_tau², below 1e-8. The native frozen
flow has 33,542 rows and predictor fill 38,314,871, below the unchanged 40M cap.
These are initial/frozen results, not evolved spatial convergence.

[The target plan](cao-tamura-accepted-10-evidence/target-preflight-v1.json)
prospectively fixes a .001/.0005 paired stepping check to .002, then conditional
stages toward .407 (at least .01 transit). Actual target minimum dx=.5 and
minimum dz=.0029296875 support the strip-informed initial timestep choice;
target uncertainty must still pass independently. The first-stage allocation
is 24,000 process-CPU seconds, 28,800 active wall seconds, 24GiB output and
2,048 pair steps, with at most 32 steps per native restart segment. No stage is
admitted from collector completion alone.

Inspection identified a separate required diagnostic-capacity prerequisite.
[Its contract](cao-tamura-accepted-10-evidence/target-audit-capacity-v1.json)
uses the actual initialized target field in a read-only frozen transport-load
fixture. It is explicitly not consecutive accepted-state budget evidence.
The sealed CT9 observer fails its **first native control** with `std::bad_alloc`:
159.06 process-CPU seconds, 159.82 elapsed seconds, peak RSS11,387,496KiB,
peak virtual address space16,777,216KiB at the declared16-GiB bound. Later
controls are unattempted under the predeclared first-failure rule.

`PhysicalScalarBudget()` constructs a complete
`patch.PreservationPacket(preK,time,data,true,{},{})`, then parses its volume
rows into a property tree. This global serialized catalog is the diagnosed
resource concern; no new global dense Schur or M^-1G allocation was introduced.
The exact failing allocation site has not been instrumented; this is a failed
required-observer execution-capacity test, not an allocation-site proof.
Sparse governing capacity retains its CT9 scope. Required physical-observer
capacity does not pass. Therefore the target paired trajectory, accepted-step
profile, evolved fill and .01-transit stage are **unattempted**, actual T=0.
Full target physical budgets and temporal uncertainty remain unqualified.

## Costs, resources and identities

The refined .001/.0005 T=.008 solver processes use 58.77/113.25 CPU seconds
and 60.00/114.99 elapsed seconds. The .008->.02 segments use86.96/167.01 CPU
and88.45/169.30 elapsed seconds. Instrumented exclusive accepted-step CPU is
about6.16–6.81s; startup/output are separate. For the 24-step fine continuation,
exclusive domain computation is79.20s and flow assembly35.10s. Nested inclusive
totals are not added. These are strip measurements, not target accepted costs.
The target initial/frozen preflight uses166.85 CPU /170.26 elapsed seconds and
peak RSS10,086,752KiB. Detailed timing, process memory, suspensions, commands and
loaded-library catalogs are in the compact evidence.

Every new numerical/build/postprocessing family uses the owned launcher,
one-thread numerical libraries, serial execution, at most two compiler jobs and
one family at a time. The scientific collection peaks at 14 project /14 aggregate;
the final read-only build supplement records **14 project /18 aggregate**,
including concurrent external activity, within 16/72. Collector and publication
receipts close separately and pass.
The controlled address-space rejection is a capacity failure, not a worker-limit
violation. Historical resource incidents are unchanged. A transient SSH outage
interrupted a source transfer after numerical preflight closure; no trajectory
or sealed artifact was altered.

Governing executable SHA256:
`5830a2cd8a4b0b5b12095555171eb28bc0f63e1daa97652d63da244c3fd808c7`,
sealed CT9 component `CSR1:EigenLU-LDLT-IR3-F40M:FGMRES80x1000`.
Read-only probeV2/V3 were built at `124006f7` / `a3250a81`, SHA256
`40ce227b3cb460858b18e9d74c8a873f4713dc5bc115e0543a0a5e84222cb9ae` /
`21fd68d440539ecb40bd39a401885e71fcb82175d789efd01a313483882fb0d1`.
The final input-policy probeV4 was built at `ec48cb35`, SHA256
`41f5e9d66b61d4b39188816060495f9543307151c384b0982d099cdd2f2ca42f`;
its sealed build/input controls and separately closed resources are in
`observer-input-v2/`.
The numerical collector ran at `7e8d4258`; publication at `baf29937`.
Complete governing/diagnostic/controller roles, compiler settings, MPI/PETSc
flags, physical-file hashes, state/history/configuration/boundary identities and
actual mapped libraries are in `qualified-v1/source-scope.json`, `runtime.json`
and `publication-v1/provenance.json`. Production physics/PIMPLE/sparse sources
and the authoritative physical manifest are unchanged from the reviewed head.

## Decision and next package

The signed-load diagnosis, one representation contrast, matched strip
qualification through .02, backend-transfer control, target initial/frozen
preflight and fresh resource compliance are demonstrated. General spatial
convergence, target evolution and physical experiment applicability are not
established. The CT8/CT9 failures, eight extraction failures, curved CT6/CT7
temporal/degree gates, CT4 pressure interpretation, A93 narrow CT5 margin,
roughness-sublayer/thin-layer issues and smooth preparation remain unchanged.

The next bounded package is **a streamed physical-scalar observer**, with
small-case all-row/face/actual-stage equivalence and target memory qualification,
followed by the already contracted paired target stepping/cost check. This is
the immediate gate to full-domain flat preservation. Mapped field numbering,
joined Cartesian gradients/cofactors, conservative transmission and curved
temporal/spatial qualification still precede a physical rough-hill pilot.
See the [ready-to-use prompt](cao-tamura-accepted-10-evidence/next-package-prompt.md)
and [evidence/reproduction index](cao-tamura-accepted-10-evidence/README.md).
