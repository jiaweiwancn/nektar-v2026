# CT14: short full-domain A95 flat-pilot qualification

Status: **PASS — short numerical qualification only**. All five CT14 closure
criteria pass. Two fresh .00025/.000125 histories reach T=.002 in 8/16 accepted
steps, with every unchanged physical/uncertainty check passing. Absolute-pressure
spread **2.287441372732e-07 <= 3.144090599826e-07**, leaving 8.566492e-08
margin (27.25% of the allowance). The initial .001/.0005 failure is preserved.
The source/evidence checkpoint containing this report closes delivery; the exact
pushed commit and remote containment are reported with the final handoff.

The fetched starting head was `f121dd1e5dc14ed147714876ea182a913e6fcae5` on
`codex/a93-a95-e94e95-awf`; no newer remote progress existed. Development uses
the dedicated Mac worktree. Compilation and all numerical conclusions use
Linux, in the isolated archive
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-14-20261115`.
The verified mount is `/dev/sda2`, `fuseblk`. Unrelated main checkouts, dirty
work, shared installations, historical failures and histories are preserved.

| Closure criterion | Fixed requirement | Measured result | Status |
|---|---|---|---|
| Required frozen/evolved systems | 70k slots;20M matrix nnz;40M entries/factor;16GiB address space;original true infinity residual<=2e-12 maxabs(rhs)+1e-25;IR<=3 | 62533 slots;final full-flow maximum 2911762 nnz;qualified-factor maximum 37266316;all direct and24 iterative pressure checks pass;actual direct IR<=1 | PASS |
| Component/operator/update/history | actions/reconstruction<=2e-10;state/absolute-pressure/traction equivalence<=1e-8;native restart/full rollback/incompatible rejection | maximum action4.313e-15;original field1.355e-14;traction4.071e-14;all history controls pass | PASS |
| Complete selected static admission | Native/alpha1000/2000/full mass/G/D/C/scalars;reference, gradients, signed loads, boundaries, domain, covariance, radicand, cutoff, quadrature | All required selected-space checks and all four initialized-frozen scalar controls pass;changed-alpha4000/8000 controls also pass | PASS |
| Fresh pairedT=.002 qualification | Every unchanged CT11 physical and uncertainty criterion | DeltaPi2.2874413727322336e-7<=3.144090599826331e-7;all 9/17 accepted states, required stages and all physical controls pass | PASS |
| Reproduction and finite resources | Pushed source;runtime/input/actual-history seals;24kCPU s/28.8kactive-wall s/16GiB new output | 24 closed owned families;20350.165 measured CPU s plus explicit source/transport reserve;closed accounting and compact sealed evidence linked below | PASS |

The [prospective contract](cao-tamura-accepted-14-evidence/contract.json)
was pushed before costly execution. The physical manifest remains SHA256
`004bc0fc82054296a68d56de77751fc1f0a41570ee16b0a174ae8087f262d296`.
The frozen experimental dataset, literature files, calibration and original
CT11 classifier/budget source are unchanged; the digitization audit was reused.

The physical floor remains r=0, Hh=40mm, HA=7.5mm, d=3.5mm. This is the accepted
SamePhysicalWallDisplacedAcceptedV1/DisplacedNeutralCompositeV1 extension,
with the same three physical-z<=25mm fit points0/1/2, fixed measured friction
velocity and normalization, natural q/kappa and fixed roughness. Actual evolving
K=k(HA) drives wall coefficients, cutoff and sink;hs remains fixed. Nonconstant
inner k, dependent inner epsilon, finite inner viscosity, high-Re outer equations,
full signed Cartesian stress/production and conservative physical transmission
remain. Inlet, Richards–Norris roof, absolute nonlinear outlet traction using
D(r;KH), IMEX1 single sweep and no-cache policy are unchanged. No clipping,
residual subtraction, reference reset, pressure mean/gauge change, trace repair,
extra governing sweep or source/model retuning occurs.

The fixed solver reference scales are ustar=.07426793792441279,
K0=.021081109428869052, outer epsilonH=.01013811394236224 and
Pi scale ustar²=.005515726603544432. Compatible wall constants remain
c_mu^w=q^-2 and c_l^w=kappa*q^(3/2).

The selected space remains 62533 stored slots/961 elements, P4/P2,
x/Hh=[-25, 30], roof/Hh=20, HA/Hh=.1875 and the first-inner-slab insertion at
HA/128=.00146484375. Its initial shear-pressure-image reduction96.27% is a
static diagnostic; it never supplied the evolved-pressure decision.

## Capacity, component and static evidence

CT14-VectorCOLAMDV1 is an opt-in predictor ordering. It builds the directed
quotient of the actual canonical matrix pattern by common native velocity FE
IDs, applies unchanged Eigen COLAMD defaults and expands each group in canonical
order. Both matrix indices are permuted internally, every coefficient is checked
exactly, and solutions/residuals return to canonical field IDs. Scalar and
preconditioner COLAMD, consistent-mass LDLT, pivots, scales, explicit-zero policy,
complete G/D/C actions and residual/refinement bounds are unchanged. Governing
pressure uses D[(alpha*Mfree)^-1 G]-C;diagonal mass is only a preconditioner.

| Frozen system | Matrix nnz | Factor entries | Margin to40M | Complete ordering-probe CPU s |
|---|---:|---:|---:|---:|
| Original initial, alpha1000 | 2081794 | 35438135 | 4561865 | 2.787 |
| Original evolved, alpha2000 | 2082125 | 35399996 | 4600004 | 2.844 |
| Selected initial, alpha1000 | 2152435 | 37266312 | 2733688 | 2.571 |
| Selected initial, alpha2000 | 2152435 | 37225528 | 2774472 | 2.529 |
| Selected initial, alpha4000 | 2152435 | 37199662 | 2800338 | recorded in timestep admission |
| Selected initial, alpha8000 | 2152435 | 37162966 | 2837034 | recorded in timestep admission |

These are complete-command times, including input, permutation checks and
actual/tiny/cancelling RHS controls, not isolated LU timings. Required original
actual-RHS residuals are <=2.885e-19;tiny residuals<=3.885e-31. New alpha4000/8000
actual residuals are<=3.854e-19, tiny<=3.760e-31 and cancelling<=2.004e-25;
each passes its own unchanged RHS-scaled bound. Inverse permutations and the
independent nonsymmetric48-row manufactured control pass. Selected systems
use exact canonical coefficients/RHS/masks/history and independent controls;
no impossible old selected solve above its fill guard is required.

Three measured failed policies are retained: scalar Scotch ND40862187 entries;
native-vector Scotch ND40119892;keeping its external ND order40119892 again.
The identical third result rules out Eigen postordering as the fill cause.
The directed COLAMD quotient then passes all required systems. Each attempt has
a prospective amendment and immutable witness;no seed/threshold sweep or
known-failed rerun occurs. CT13's unfinished Natural trial remains work-limited,
not a measured fill result.

Four complete small zero/displaced P3/P4 operator controls pass. Maximum scaled
independent action defect is 4.3125225612783424e-15. Fresh original full-domain
.001 two-step evolution agrees with the qualified original update at every
common state:maximum scaled field difference1.3541903104511208e-14, including
absolute pressure, and traction4.070102570243455e-14. Its historical pressure
failure remains 3.824614134304275e-6>5.71117429897371e-7. Genuine same-component
restart and full rejected-step rollback reproduce all native histories bitwise.
Old-component, new-to-old, foreign-fit histories and unknown ordering reject.
No old history is relabelled. The new component has its own identity suffix;
default histories retain their old identity.

The complete selected static state reproduces CT12/13's actual coefficients,
SHA256 `dd0af02695e58934d6469580aec40249700b404b7f079c0abe8db10c99d3c77a`.
Both full-command factor-free exports have zero sparse factorizations and
exact canonical matrix/RHS/history/mask/lift equality to the CT13 inputs.

| Selected initial/diagnostic quantity | Measured | Criterion/role |
|---|---:|---|
| Scaled reference RMS, U/k/outer-epsilon | 9.85149e-6/1.03679e-5/6.05198e-7 | each<=1e-3 |
| Scaled reference RMS, w/absolute Pi | 0/3.01921e-16 | <=1e-4/1e-3 |
| Maximum scaled physical-gradient reference RMS | 4.74436e-4 | U;all five fields retained |
| Ground/full-interface physical traction error | 1.41528e-15/2.16173e-15 | each<=1e-3 |
| Physical outlet reconstruction mismatch | 8.44214e-5 | natural-face diagnostic, using D(r;KH) |
| All-mode signed momentum reconstruction | 9.74901e-15 | <=2e-10 |
| q48/q64 maximum scaled pressure-action change | 1.28385e-10/2.44748e-10 | each<1e-8 |
| Independent initial inner k minimum | .021081109428869045 | >0 |
| Independent initial full3D covariance minimum | .006110429044108287 | >=-1e-10 K0 |
| Independent initial A95 radicand minimum | .6916508330449576 | >=0 |
| Static scalar-control maximum VmPeak, KiB | 15894660 | <=16777216 |

All native hard/domain/cutoff flags pass. Native initial regional Pi-reference
errors are zero;independent Python accumulation gives the roundoff global RMS
above, with no mean removal. All wall test modes, physical gradients, signed
boundary/volume loads, flow/scalar/consistent-mass/full-pressure actions are
retained. Native, q48f16, q64f16, q64f32 initialized-frozen scalar controls pass
capacity and full source/physical-flux reconstruction;their role is static,
not an evolved budget pass. Complete bodies are external and SHA-bound by the
[compact static evidence](cao-tamura-accepted-14-evidence/static-v1/).

The [changed-timestep admission](cao-tamura-accepted-14-evidence/timestep-v1/results.json)
passes both alpha4000/8000 canonical exports and all six native/q48/q64 full
action controls. G/D/C/M and initial/history identities are exact;independent
A(alpha)-A(1000)=(alpha-1000)M defects are 1.56146e-16/1.33840e-16.
Maximum q64 pressure-action change is 1.95452e-9<1e-8;maximum address space is
15344516KiB<16777216. Its complete result SHA256 is
`ed814ff3b1e085ebd588e59f9bac15b269a70209aeb55e7275f776ec2a52548b`.

## Fresh physical histories

The first contract was pushed at `eca8f3ab` before two fresh independent
.001/.0005 histories. Both reachedT=.002, with all physical, budget, domain and
uncertainty checks passing except absolute-pressure timestep uncertainty:

| First-pair quantity | Measured | Criterion | Status |
|---|---:|---:|---|
| Own-initial normalized Pi drift, .001/.0005 | 1.7107716728451436e-6/2.519311058540182e-6 | each<=1e-3 | PASS |
| Absolute Pi pair RMS | 9.147839144604025e-7 | <=2.519311058540182e-7 | FAIL |
| Maximum physical-scalar budget defect | 1.189617257349346e-6 | <=.02 | PASS |
| Clock-off native CPU s, .001/.0005 | 366.035/662.695 | measured complete commands | RECORDED |

The full failure/independent collection SHA256s are
`dc11f77a2b60133d4a501f9ddf18929144c146dd0871ef444c927f2c08edbda3` and
`68cf158ba29eb553a0415e4bc5943930d33b9c88d140406f923b42cf837e3673`.
Absolute levels, regional means/variations, all common times, actual consecutive
states/stages and every factor/solve remain archived. No endpoint is retimed.

Backend/update differences and quadrature uncertainty are far below that
spread;all original-system solves pass with at most one refinement. The
single-sweep momentum remainder and separately retained refreshed post-k epsilon
lag support testing finite-step error. The lag falls from1.14539e-7 to6.05308e-8
when dt halves, while imposed traces remain roundoff-scale. A first-order
planning estimate predicted a .0005/.00025 pair still fails (4.574e-7 spread
against2.947e-7 allowance). The smallest dyadic pair predicted to pass was
.00025/.000125 (2.287e-7 against3.164e-7);these estimates are not acceptance or
formal convergence evidence. The prospective amendment was pushed before new
admission or evolution. No physical space/model/coupling change was made.

The [retry contract](cao-tamura-accepted-14-evidence/short-pair-contract-v2.json)
was pushed at `b4867df7` after admission. Two new zero-step histories were created
before either advanced, with identical actual initialized coefficients and their
own dt/native checkpoint. Both consume only their own histories. CT11 source
is unchanged;the request specialization changes exactly its required-dt list
and member-name set. Every scientific check, tolerance, scale, T=.002 and actual-dt
coverage rule is unchanged. Actual first-pair parity, physical/request negatives,
explicitly synthetic new-request positive, foreign-old-cohort rejection and
coverage negative pass. Synthetic guard/classifier fixtures are never physical
histories or actual resource expenditure.

The final pair passes the complete classifier and independent collection. Both
fresh initial histories have actual state SHA256
`dd0af02695e58934d6469580aec40249700b404b7f079c0abe8db10c99d3c77a`.
All 9 coarse/17 fine accepted-state records pass; every required native stage
passes admissibility. Fine raw time is .002000000000000001, within the unchanged
1e-13 time coverage rule; no state is retimed. All 9 common times and whole/inner/
outer absolute-pressure content remain in the complete evidence.

| Field | Coarse own-initial drift | Fine own-initial drift | Coarse absolute reference error | Fine absolute reference error |
|---|---:|---:|---:|---:|
| U | 1.879743081e-06 | 1.904741645e-06 | 9.224864555e-06 | 9.221536295e-06 |
| w | 2.724403875e-08 | 2.809720643e-08 | 2.724403875e-08 | 2.809720643e-08 |
| k | 1.290425755e-06 | 1.304011932e-06 | 1.050568217e-05 | 1.050800056e-05 |
| outer_epsilon | 1.311885184e-07 | 1.342032373e-07 | 5.925773992e-07 | 5.910252838e-07 |
| Pi | 2.945177709e-06 | 3.144090600e-06 | 2.945177709e-06 | 3.144090600e-06 |

These are endpoint values, normalized on the original scales. Initial
representation error is separately tabulated above. All accepted-state drift
limits remain 1e-3 (U/k/epsilon/Pi) and1e-4 (w).

| Common-endpoint uncertainty | Normalized spread | Original allowance | Status |
|---|---:|---:|---|
| U | 2.743636751e-08 | 1.904741645e-07 | PASS |
| w | 1.974208720e-09 | 1.000000000e-07 | PASS |
| k | 1.438715889e-08 | 1.304011932e-07 | PASS |
| outer_epsilon | 4.314681850e-09 | 1.000000000e-07 | PASS |
| Pi | 2.287441373e-07 | 3.144090600e-07 | PASS |
| Traction ground | 2.199237602e-08 | 9.477054909e-07 | PASS |
| Traction interface | 7.318313695e-06 | 5.259686946e-05 | PASS |
| Traction signed_mean | 2.182916594e-08 | 9.459510498e-07 | PASS |

The absolute Pi norm includes its level. AtT=.002 the whole/inner/outer raw
pressure-difference RMS values are 1.2616901234e-9/2.5923570313e-9/1.2423074144e-9;
whole signed mean9.6806856554e-10 is retained, not subtracted. The normalization
is the unchanged ustar²=.005515726603544432. The pressure spread falls by almost
exactly4 relative to the first pair when both dt values are quartered; this
establishes this declared temporal uncertainty test, not spatial convergence.

| Endpoint physical witness | .00025 | .000125 | Criterion/role |
|---|---:|---:|---|
| Ground traction error | 9.459510498e-06 | 9.437681331e-06 | <=1e-3 |
| Full-interface traction RMS | 8.319604851e-04 | 8.446814744e-04 | <=1e-3 |
| Physical outlet reconstruction RMS | 7.299570258e-05 | 7.289851635e-05 | natural-face diagnostic D(r;KH) |
| Normalized net mass | 2.211851361e-14 | 3.237952319e-15 | <=1e-8 |
| Maximum physical-scalar budget defect | 1.102024936e-06 | 1.088841258e-06 | <=.02 |
| Maximum signed-defect quadrature spread | 5.363106808e-15 | 1.237213973e-14 | <=.02 |
| Imposed epsilon trace residual | 1.845739477e-21 | 2.442792593e-21 | reported separately |
| Refreshed post-k trace lag | 3.119507546e-08 | 1.584719243e-08 | retained single-sweep lag;no repair |

All native/q48f16/q64f16/q64f32 budgets use the actual consecutive accepted
states and actual flow/pre-k/history snapshots. Every equation and physical
source/flux reconstruction passes; prescribed-row reactions are not substituted.
The maximum native momentum splitting-identity defect is 4.060e-19, while the
expected uncancelled single-sweep momentum remainder is recorded separately.

| Minimum across required native accepted states/stages | Value | Criterion |
|---|---:|---|
| inner_k | 0.021081013846 | >0 |
| outer_k | 0.0210810139336 | >0 |
| outer_epsilon | 5.0913315467e-05 | >0 |
| inner_covariance | 0.00862714441413 | >=-1e-10 K0 |
| outer_covariance | 0.00611015427944 | >=-1e-10 K0 |
| A95_radicand | 0.691581016567 | >=0 |
| required_cutoff_ratio | 1.46864739378 | >=1;d<c<=HA |

Independent physical endpoint inventories and constitutive comparisons also
pass. Located extrema retain time, stage, coordinates, native side and geometry
in the complete bodies. Finite output is not used as an admissibility proxy.

| Clock-off native trajectory | Steps | Factor receipts | Direct-solve receipts | Pressure solves | CPU s | Wall s | VmPeak KiB |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0.00025 | 8 | 73 | 2729 | 8 | 1269.812781 | 1279.968576 | 10062772 |
| 0.000125 | 16 | 137 | 5399 | 16 | 2514.342541 | 2531.764635 | 10154516 |

Both zero-step histories additionally have9 factors/35 direct-solve receipts
each, with maxIR0. All evolved direct solves satisfy the original true-residual
bound, with maxIR1 and largest residual/bound fraction.997173909. Independent
recovery from every actual native flow snapshot checks all 24 unscaled pressure
RHS/residual pairs;maximum residual/bound fraction.001730891. No extra solve
is performed by collection. Final full-flow maximum is 2911762 nonzeros;final
predictor maxima37199666/37162970. The complete-command times include startup,
mandatory checks, history and output, with detailed timing off;they are not
exclusive LU timings, steady-step rates or proof of long-run affordability.

Full passing pair SHA256:
`dfeb102cba472160f0f3e0718a85353710735d759130c438fc9275acad438a60`.
Independent collection SHA256:
`f2be1f0835ef654a21a2945d1417cd61208b14e77671bb40f85a884d82df4b03`.
The [sealed compact closure views](cao-tamura-accepted-14-evidence/closure-v1/)
retain every aggregate primitive and decision and bind complete external bodies.

## Runtime, accounting and reproduction

Native tested source is `6477d9d988e969686a418c8e333fadfb85bb3f9d`.
Governing executable SHA256:
`9ae7d3c63492c82a4c9f979616f673d8e4909a5c41ee8069fc1592f2fb23da06`.
Probe SHA256:
`95d6a0258926ab250b69ea0940976f43a993bf4bdccc3001a5819f4f36835824`.
The private runtime is `build-runtime-v1/runtime` under the archive. Exact
compiler/build flags, loaded libraries and hashes are in the
[capacity/component runtime evidence](cao-tamura-accepted-14-evidence/capacity-component-v1/).
Component/static/first-pair/changed-alpha/retry controller scopes are respectively
6477d9d9/087fa720/eca8f3ab/ae18441f/b4867df7;immutable requests bind each source
file, configuration, runtime and predecessor seal. Later publication edits do
not rebuild or relabel these executables.

All compilation and numerics use the authoritative Linux environment: MPI
ranks=1, OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1 and the other four numerical
library thread variables=1. Compiler jobs are at most2, with one numerical
family at a time. Recorded owned-family worker high water is 13 and aggregate
high water 23, including launchers, monitors and MPI helpers, within 16/72.

All 24 controlled families are closed at **20,350.164785 process-CPU seconds**.
The explicit 120 CPU-second source/Git/transport/inspection reserve gives
**20,470.164785<24,000**. Conservative whole-package elapsed wall at final
accounting is 20,046.45s, including source/Git/inspection; adding 1800s for final
Git delivery gives **21,846.45<28,800s**. Linux retained output is 9,896,853,512
bytes; known Mac publication/transfer files add 52,948,738 bytes. An additional
256MiB Git/transfer reserve gives **10,218,237,706<17,179,869,184 bytes**.
All failed attempts, controllers, audits, builds and publishers remain charged.
The live package guard completed without reaching its allocation stop.

Direct owned-native observations retain a largest sampled VmPeak of
**16,122,596KiB (15.38GiB)**,654,620KiB below the16GiB address-space guard. This
observed peak is not claimed as every process's final peak; all native jobs
inherit the hard16GiB limit. Exact closed receipts, publication supplement,
source/transport reserves and the memory witness are in
[closed-publication-accounting.json](cao-tamura-accepted-14-evidence/closed-publication-accounting.json).
Python3.10.12, NumPy1.21.5 and SciPy1.8.0 are recorded with the Linux environment;
long double storage is 16bytes. These metadata and source/transport reserves
are separate from the native command CPU benchmarks.

The [archive/reproduction index](cao-tamura-accepted-14-evidence/README.md)
contains commands and archive roles. Use fresh destinations, matching private
runtime/library bytes, request/completion seals and the owned worker launcher.
Retry collection/publication includes --admission for the exact changed-alpha
seal. Full CSR, states, histories, stage packets and per-point faces remain outside
Git. Compact views bind them by path/SHA/count and retain aggregate primitives
and classifications;they are not classifier/checkpoint/replay inputs. The
initial static publication copied a large initialized vector;the normal
additive correction preserves that checkpoint while current views use external
references. Historical publication/collector seals keep their original scope.
The first CT14 publisher's missing private-library environment was rejected and
charged;the matching environment passed without changing native data/runtime.

This passingT=.002 pair establishes only this short numerical qualification,
not spatial convergence, experimental validation or affordable preservation.
Long flat preservation still requires its own longer-time physical/uncertainty
contract and measured whole-campaign CPU/wall/output affordability, including
startup, mandatory audits, controller costs and qualified retention/reuse if
optimized. Lower factor fill is not a throughput qualification. The first
rough hill still requires mapped region-aware numbering, joined-map Cartesian
cofactor/gradient checks, conservative transmission, the outstanding CT6/CT7 curved
spatial/temporal gates, and roughness-sublayer/thin-layer/approach-flow physical
admission. Positive geometry alone supplies none of these. No meaningful-duration
preservation, terrain, other-model displaced/smooth production, EP, MPI/3D, variable
height or new calibration is launched by CT14.
