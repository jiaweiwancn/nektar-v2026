# Nonlinear duration and response qualification

A95 Natural / Accepted passes the complete four-member neutral qualification
through T=.5 (theta=.1). A93 reaches the same duration with passing physical
accuracy, admissibility and timestep bounds, but fails the fixed-common vertical
momentum reduction; its response is not admitted. A95's T=.016 controlled
response resolves pressure, vertical velocity, nonlinear shear stress and full
production with passing space/timestep uncertainty. Selected A95/EP and
half-amplitude controls complete T=.016. Their field comparisons pass, but
native weak-balance failures prevent full response qualification or extension.

The Accepted response's T=.1 extension stops at a read-only native weak-balance
gate. A frozen control identifies incorrect containing-element selection at
moving-cutoff integration points; the original response remains unqualified.

Starting revision: `dfd94161fd9d97f7f1a8ce2966e9cc9461b6d686`, clean implementation
worktree, no newer local or remote commits. The previous evidence and eight
unresolved epsilon temporal-order observations remain immutable. Specifications
retain SHA256 `056a85c2aea5ce4ca168caf561613e662416bb6b7d7177451ba5038d2cf68a64`
and `17c7d6a241ecf8f929c9232ef34e3a6f48823cae2c0b9357911a377d28768a40`.

New authoritative Linux root:
`/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-nonlinear-20261110`.
Source, build, copied immutable runtime, diagnostic runtime, tmp and trajectories
are isolated. The primary Linux checkout and concurrent CLS/WJ jobs are untouched.
The initial relevant thread census is 39; numerical libraries use one thread,
initial native concurrency is two, and the aggregate limit is 72. Resource
availability is checked again at each launch. The existing 3000-unknown guard
and unsupported-route guards remain unchanged.

## Contract declared before new evolution

No governing turbulence/wall equations, coefficients, sources, molecular policies,
one-sweep stage order, pressure convention or native histories change. Eta and
corrections remain zero. No reference cancellation, resets, trace repair,
radicand clipping or new regularization. B5/B11 retain 2355/2997 unknowns, p5/p3
spaces, fixed model-specific meshes and 48/12/12/12/12 integration.

Admit A93 and A95 Natural Accepted independently. Validate both .0016/.016
short comparisons, both model-specific common-test/spatial receipts, scalar
observation coverage, actual consumed inputs and byte-level runtime equivalence.
E94E95 success cannot satisfy either admission. New four-member cohorts use
.002/.001 and milestones .016,.1,.5. If uncertainty fails while admissible,
allow one fresh .001/.0005 cohort; further refinement needs a recorded limiting
observable and a new bounded contract. No checkpoint timestep is edited.

At every milestone all four field/completion/identity/source gates and
`max(delta_t_B5,delta_t_B11) <= max(1e-7,.1*delta_h)` must pass before continuing.
Neutral field RMS limits remain U,k,epsilon,Pi <=1e-3, w<=1e-4, mean wall/full
interface traction <=1e-3, on the existing fixed physical scales. Full native
and fixed common-test momentum-x/inner-k residuals, and non-floor nonlinear-model
vertical momentum, require fine/coarse <=.85 above the existing 1e-8 residual
floor and five times evaluation/quadrature uncertainty (absolute 1e-10).
Analytical balance, reconstruction, reduction, field/time accuracy and execution
are separate statuses. Other equations and new fine modes remain visible.

Final requirements include both spaces' alpha-zero audits, all four trajectories'
scalar observations through .5, frozen quadrature on both, and independent full
R/Q/production/normal-traction checks on actual native samples. Accepted states
are counted uniquely as well as with repeated restart endpoints. The one-column
epsilon constant test remains a prescribed-row reaction identity.

After a model passes neutral duration, use StreamfunctionPolyV1 amplitude 2e-5,
and a selected half-amplitude 1e-5 control, at .016 then conditionally .1. The
streamfunction and first normal derivative vanish at the essential boundaries
and interface. Bind the actual perturbed initialized vector/history separately
from the neutral reference and boundary map. Same-space dt controls require
identical actual initialization. Compare response to a matching unperturbed
numerical trajectory, never to a regenerated neutral initializer.

Response spatial differences must be <=max(5e-7,.1*response RMS); timestep
differences <=max(1e-7,.05*response RMS) on the existing physical scales. A
component below these absolute floors is unresolved. The half-amplitude test
reports departure from linear scaling with limit max(5e-7,.1*half-response RMS),
without assuming exact linearity or monotone decay. Retain the existing hard
physical and absolute field/traction accuracy gates. A93/EP is the preselected
short wall-source/production comparison, with all four companions at .016;
extend only for a specific resolved wall question. No temporal order, global
stability or open-outflow claim follows from these experiments.

Pressure, vertical velocity, Qxz and full-production responses must be resolved
for extension to .1; individual smaller scalar/stress responses retain their
unresolved labels. Stress response uses scale u_star^2 and production response
uses epsilon(H). Initial coefficients and native zero-step history are sealed
separately from the neutral reference; their stable artifact identities are
checked on every same-space temporal comparison and exact resume.

The declared response observables also include all three Q normal components,
mechanical pressure `pi=Pi-2k/3`, and full vertical-plane normal/shear traction.
Mechanical pressure and normal traction use the same single gauge shift as Pi;
they are not independently regauged. Traction response RMS over volume samples
is distinct from the native interface-traction-jump accuracy gate.

For exactly divergence-free planar flow, A93's Q:S vanishes because the two
nonzero strain eigenvalues are opposite and Omega^2 is a multiple of the
identity on that plane. A95's commutator has zero contraction with S generally.
This special kinematic fact does not permit dropping full production: the
variable viscosity responds, Qxz can be nonzero, and the discrete raw gradient
can have divergence. Report the nonlinear and isotropic production contributions
separately and keep the general three-dimensional constitutive evidence intact.

The optional physical scalar-budget diagnostic reconstructs the actual IMEX1
flow/epsilon/pre-k stages from consecutive native accepted vectors. It integrates
storage and the actual frozen-stage sources separately from physical diffusive
and advective face fluxes. Current-coefficient sources/fluxes are also reported,
so frozen-operator splitting is visible. Epsilon storage includes its moving
dependent boundary coefficients. Prescribed-row reactions do not enter this
calculation. Nonzero physical finite-element budget defects are classified by
quadrature, space and timestep comparisons, not forced to zero. This diagnostic
does not establish a mean-plus-turbulent energy theorem.

Before the response experiments, the physical-budget controls are declared as
volume orders 48/64 and independent physical-face orders 16/32. Both primary
spaces and both timestep companions are observed. Budgets are signed integrated
rates divided by volume times `K0*u_star/H` (k) or `epsilon(H)*u_star/H`
(epsilon), rather than RMS field norms. Quadrature changes must be at most
`max(1e-9, .01*spatial_difference)` to interpret their spatial differences.
Finite-defect spatial reduction uses the same .85 ratio above a 5e-9 floor and
five times quadrature uncertainty. Timestep differences in those defects are
compared with `max(1e-9,.1*spatial_difference)`. These diagnostic classifications
remain separate from the established field-response gates; unresolved physical
flux or temporal contributions are reported explicitly, without forcing a zero
budget or identifying prescribed-row reactions with physical face fluxes.

## Boundary configuration and interpretation

The physical patch is `-1 <= x <= 1`, `0 <= z <= 4`, with interface `H=2`,
roughness `.1`, `u_star=.4` and physical molecular viscosity `1.5e-5`.
Velocity is prescribed on the two streamwise sides and top, and vanishes at
the wall. Scalars have prescribed side/top traces; the wall k condition is
natural. The outer epsilon trace at H is dependent and frozen during its
declared scalar stage. Pressure uses one global mean gauge. The common velocity
and k spaces couple the interface weakly; reconstructed one-sided traction and
flux jumps are additional diagnostics rather than extra interface equations.

`FrozenReferenceV1` binds the external map to the stored neutral reference.
Actual perturbed initial coefficients and their native zero-step history are
bound separately. Across a restart, repeated endpoint observations must have
the same full state hash, including the retained epsilon trace. Native legacy
D entries remain differences from the fixed neutral polynomial reference;
the Python physical comparisons additionally measure drift from the actual
perturbed initial field. Neither quantity replaces analytical-reference error.

The whole-patch physical k budget subtracts the outward sum of the two internal
H-face fluxes from the regional face totals before forming the external-boundary
budget. The internal flux mismatch remains visible separately. Independent
integral fixtures verify this accounting and verify that changing a prescribed
row reaction cannot change the physical budget. In this one-column patch the
epsilon constant-test coefficients that are nonzero are all constrained: its
native load/reaction identity is not an independent conservation theorem.

## Verification changes

Admission and finishing now bind model/policy/wall, the effective descriptor,
both exact meshes, both timesteps, required time sets, the model-specific common
packet, runtime/library identities, and sealed state/history/boundary artifacts.
Every accepted scalar observation must carry the matching state hash. Native
report intervals, step counts and timestep increments are checked, with unique
accepted states distinguished from duplicated restart endpoints.

Required CLI negative controls cover missing times, foreign model/wall receipts,
missing spatial/timestep companions, truncated scalar observers, and failed
dominant residual reduction despite passing analytical/reconstruction checks.
A real partial export also preserves a complete report and exits 1. Passing
mock receipts do not qualify actual numerical cases. A separate real control
rejects unperturbed neutral trajectories presented as perturbed inputs.

The final collector also requires conditional response/audit/physical-budget
coverage, the selected EP and half-amplitude results, and current loaded-library
maps rather than inherited maps. Exact-request resume checks completed hashes;
changed requests reject before replacing results. Completed earlier evidence is
not overwritten by operational recovery or supplemented exports.

The added forward-tensor probe initially used an exact H comparison whereas
the observation path uses a 1e-12 interface tolerance. A first fixture attempted
to sample both sides at H +/- 5e-13. It failed because an inner wall profile is
only defined on [0,H]; that failed fixture and runtime are preserved. The final
probe scope requires exact H with an explicit side for interface observations,
and rejects nearby ambiguous coordinates without extending or clipping the
wall profile. Tensor and test-basis exports use consistent side selection.
Separate before/after regressions exercise all six natural model/wall
initializations and unchanged exact/interior sample values. The primary
campaign grids use exact H or interior points and retain their frozen runtime.

The old A95 B11 initialization sampler hit its legacy 600-second wall deadline
while resource suspensions were in use, before any T=.1 advancement. The failed
comparison/log remain immutable. `init-timeout-recovery-v1` uses the same
trajectory requests and native histories with a separately recorded operational
deadline and new comparison receipts. It does not change timestep or accuracy
criteria. Required Linux Python is 3.10 or later; the Mac's default Python 3.9
cannot execute `zip(strict=True)` and supplies no numerical qualification claim.

## Completed preparation and first milestone

The stronger admission initially found missing scalar observers on the old dt/2
short runs. New read-only receipts cover those sealed histories without changing
the old campaign. The first loader-identity launch checks also rejected an
unpinned environment before any advancement; new immutable admissions bind the
explicit trajectory and diagnostic loader environments. Original failed logs
and admissions remain archived.

Both A93 and A95 .002/.001 four-member cohorts pass T=.016, including complete
scalar/time coverage, independent native R/Q/P/eigenvalue checks, frozen
quadrature and full/common dominant-residual reductions. Pressure timestep
differences are respectively 8.07688e-8 and 7.67047e-8, below unchanged limits
1.69744e-7 and 1.73943e-7. Their T=.1/.5 continuations remain conditional on each
preceding complete milestone.

The historical global common-k norm is retained. An additional sealed check
uses only common free tests supported strictly inside the inner k region;
outer and shared tests remain separately reported. It rechecks the two reused
short times and all three new duration milestones on the same fixed mass and
physical scale. The final export and perturbation admission require this
inner-only check. A regression demonstrates that a passing global-k norm cannot
substitute for a failed inner-k reduction. These supplementary receipts do not
alter the previously sealed global-k results.

## A93 first duration gate failure

All four A93 Accepted trajectories reach T=.5 (theta=.1) with unchanged native
histories. Their physical accuracy, admissibility, source/storage, scalar
coverage, frozen-quadrature, and independent constitutive checks pass. The
complete duration qualification fails its fixed-common vertical-momentum
reduction. A93 perturbation and selected A93/EP response runs are consequently
not admitted. This is a post-step qualification failure in accepted physical
states; there is no invented rejected timestep or constitutive-domain failure.

| T | Native u ratio | Native inner-k ratio | Native w ratio | Common u ratio | Common global-k ratio | Common w ratio |
|---|---:|---:|---:|---:|---:|---:|
| .016 | .023740 | .016009 | .054212 | .010210 | .005602 | .023143 |
| .1 | .034640 | .031849 | .120642 | .017310 | .008058 | .052528 |
| .5 | .631766 | .489268 | .547848 | .832472 | .610110 | **1.013463** |

Ratios are B11/B5 on fixed physical scales. The unchanged required maximum is
.85 above the 1e-8 floor. At T=.5 the common w norms are 1.0800020643e-7 and
1.0945424092e-7, so the failing observation is above the floor. The actual
same-time half-timestep states give 1.1025359676e-7 and 1.1133651181e-7, ratio
1.009822: timestep halving does not remove this failure. The common residual
load difference across spaces is 1.7347749779e-7; timestep load differences are
1.3770209677e-8 and 7.9255652265e-9, below the separately declared diagnostic
uncertainty bound .1 times that spatial load difference. This diagnostic does
not replace the original ratio gate.

At the primary T=.5 endpoints, scaled reference RMS errors on B5/B11 are:

| Field | B5 | B11 |
|---|---:|---:|
| U | 1.83999e-5 | 7.22311e-7 |
| w | 2.83147e-7 | 2.93619e-8 |
| k | 3.61579e-6 | 1.60314e-7 |
| outer epsilon | 6.56483e-9 | 2.12770e-9 |
| Pi | 2.37958e-6 | 4.15255e-7 |

These passing field errors do not override the failed spatial-residual gate.
Pressure timestep differences are at most 7.09847e-9, against a 2.30733e-7
limit. The largest assembled common vertical-momentum rows have support around
z=.04.. .5, particularly across the roughness height .1. Large separate
element stress norms contain cancelling constant-flux boundary contributions;
they cannot by themselves be called localized physical defects. The read-only
failure supplement retains those full terms and additionally separates actual
augmented pressure, the constitutive-stress remainder, and advection, and
localizes the varying-flux buckets with their constant bucket explicitly kept.
It changes no evolution equation or residual used by the solver.

## Conditional A95 response and EP control selection

A95's T=.5 evolved primary pair passes all six native/common dominant-equation
checks: u, inner k and w on native tests; u, global k and w on the fixed common
tests. The w ratios are .375064 native and .463096 common. Remaining final
scalar coverage, quadrature and constitutive checks must pass before response
admission, together with the additional inner-only common-k qualification.

A93's failed neutral gate prevents its preselected EP response. The conditional
replacement is A95/EP at T=.016, using the same B5/B11 spaces, .002/.001 histories,
StreamfunctionPolyV1 amplitude 2e-5, physical observation grid and response
uncertainty criteria. An A95 Accepted half-amplitude 1e-5 control is limited to
.016. Selection will bind the observed Accepted short response before the EP
response advances; EP neutral preparation is independently read-only/short.
The comparison asks how constant versus nonconstant inner-reference k and the
two prescribed wall-source constructions affect the baseline-subtracted
pressure, nonlinear shear stress and production responses under one A95 outer
closure. It adds no EP T=.5 or open-boundary claim. These controls retain the
existing accuracy/domain criteria and cannot admit an A93 perturbation.

A95's complete neutral final audit and the additional inner-only qualification
subsequently pass. At T=.016 the Accepted response's scaled RMS signals on the
primary pair are Pi=3.77347e-4, w=3.32463e-5, Qxz=1.42792e-4 and P=1.91826e-4.
Their spatial differences are respectively 5.42e-8, 2.08e-9, 4.00e-8 and
4.58e-8. Maximum timestep differences are 2.71e-6, 9.18e-8, 7.05e-7 and
1.28e-6, all below the predeclared response-specific limits. The measured
Accepted receipt is hashed into the A95/EP selection before its advancement.

The first EP admission correctly rejects older primary observers without
accepted-state hashes. New read-only scalar supplements are selected with the
explicit `--scalar-observer-root` interface; the original observers and rejected
admission log remain unchanged. All four supplemented short histories and
frozen quadrature pass before the new EP admission is written.

## First response failure and geometric witness

All four A95 Accepted responses reach T=.016 with passing physical accuracy,
domain and response space/timestep comparisons. Their first additional spatial
audit fails analytical vertical-momentum weak balance: maximum row load
1.2308036225e-9, scaled 2.9848536942e-9 against the unchanged 1e-9 criterion
(with the existing 1e-12 absolute floor). Native physical-term reconstruction
agrees to 2.49e-14; the epsilon constant-test reaction identity also passes.
These checks test different properties and retain separate outcomes.

The frozen B5 packet contains 156 integration points outside the bounds of the
element chosen for their test basis, by at most 9.8861908243e-12. For example,
z=.014912396178320015 lies above the mesh node .01491239617830955, but the native
locator selects geometry 4 rather than containing geometry 5. The shared
`RANSWallGlobalField::LocateElement` admits neighboring elements within 1e-11
and chooses the lowest geometry ID. `Refresh` uses this locator for the moving
wall-profile quadrature. Near a fixed cell boundary, the new profile panels
therefore acquire basis rows from the adjacent cell. The counterexample is
inside the wall region, separate from exact-interface side sampling.

The read-only control keeps every coordinate, quadrature weight, coefficient,
time and consumed boundary value fixed. It exports native test bases using
the explicit containing geometry. The Stokes defect in the radial direction
falls from 1.4848854007e-8 to 5.8684761338e-15. The analytical momentum-z weak
defect falls to 5.0401848493e-16 and the original weak-balance criterion passes
on those counterfactual test bases. The original constant-normal-traction
load's scaled native free mass-dual norm is 4.1313245389e-6; its counterfactual
value is 2.5597805389e-13. This identifies a geometric ownership defect without
feeding the counterfactual into production assembly or modifying a trajectory.

The frozen control returns nonzero for the original failed qualification. Its
successful diagnosis does not relabel the response or resume it to T=.1.
The original sealed first-failure receipt is retained. A separately named
supplement completes independent spatial, scalar and frozen-quadrature checks
after that failure. The response auditor now gathers independent checks while
preserving a nonzero overall exit; its status output also reports the actual
failed numerical check instead of printing an unconditional "pass".

The A93 T=.5 common-w failure remains separate. Its measured constant-flux
quadrature effect is much smaller than its dominant evolved common load;
the current control does not establish that ownership alone explains A93's
failed reduction. Both limitations must remain visible in the next package.

Concurrent external starts produced a conservative schedule census of 75; the
immediate active census was 67. An owned process family was suspended promptly.
The record is retained rather than claiming the original census stayed below
72. A process-family guard now reserves additional startup capacity and pauses
only owned families, targeting at most 56 active threads including external
work. Suspensions preserve in-memory numerical state; waiting Python parents
resume after their native children finish so suspension does not consume the
subprocess timeout. Gross wall timings therefore include resource waits and
are distinguished from instrumented process-CPU cost. No existing campaign,
executable, library or checkpoint is modified or signalled.

## Final A95 neutral duration receipt

The A95 Natural Accepted cohort passes all three required milestones, both
final spatial-term audits, all four scalar observers, both frozen-quadrature
controls, independent constitutive checks and the additional common inner-k
checks. Each primary observer contains 251 unique accepted states (253 records,
including two repeated restart boundaries); each half-timestep observer contains
501 unique states (503 records). Time intervals, state hashes and duplicate
identities agree with the sealed native reports. All four histories retain
their original timestep and IMEX1 method. No finer trial cohort is needed for
the declared fieldwise uncertainty gate.

| T=.5 scaled mass-dual residual | B5 | B11 | Ratio |
|---|---:|---:|---:|
| Native momentum-x | 1.126894e-6 | 3.764294e-7 | .334042 |
| Native inner k | 1.350593e-6 | 4.143729e-7 | .306808 |
| Native momentum-z | 7.182219e-7 | 2.693788e-7 | .375064 |
| Fixed-common momentum-x | 6.783273e-7 | 3.242799e-7 | .478058 |
| Fixed-common global k | 8.974764e-7 | 3.750535e-7 | .417898 |
| Fixed-common momentum-z | 2.501352e-7 | 1.158366e-7 | .463096 |

Every resolved ratio is below .85. The separately partitioned common inner-k
receipt also passes; the historical global-k value cannot replace it. Other
equations remain visible: outer epsilon native norms are 9.36418e-10 and
1.11968e-9, below the declared residual floor, so these do not establish a
spatial rate. Shared k rows, pressure, continuity, prescribed-row reactions and
full reconstructed traction remain distinct from these six reduction checks.

| T=.5 scaled reference RMS | A95 B5 | A95 B11 |
|---|---:|---:|
| U | 1.947736e-5 | 7.642459e-7 |
| w | 3.115986e-7 | 2.742108e-8 |
| k | 3.631649e-6 | 1.599723e-7 |
| outer epsilon | 6.501799e-9 | 2.071413e-9 |
| Pi | 2.583396e-6 | 4.104976e-7 |

These are weighted physical RMS divided by the fixed scales declared above,
not pointwise relative errors. True scaled L2 and actual-initial drift are
also retained. For example B11 drift RMS is U=5.23170e-7, w=2.74211e-8,
k=1.11803e-7, epsilon=1.56288e-9 and Pi=4.10498e-7. Source/storage identities
are at most 2.54e-16 for k and 7.65e-18 for epsilon on the two primary final
segments. Refreshed trace lag remains visible (up to 7.82e-12 on those segments)
and is not repaired. Native correction calls and signed/absolute correction
integrals are zero.

The B11 final segment's native all-state/stage A95 radicand minimum is
.6919428274, outer covariance minimum .1772842540, inner covariance minimum
.2476860343, outer epsilon minimum .03959798195, and required cutoff budget
ratio minimum 5.898709634. Coordinates, side, stage and time accompany each
extremum. These values describe their declared native sampling set; the
independent 800-point constitutive grid has its own margins and receipts.
Across both models' final independent grids, R, Q, full signed production,
model-specific epsilon sources, covariance and full normal/shear traction agree
with independent formulas within 3.34e-16 (absolute criterion 5e-11). No stress
component or isotropic raw-divergence work is omitted.

## Completed bounded response controls

All three A95 response cohorts use four native trajectories at T=.016
(theta=.0032), with their own immutable actual perturbed zero-step histories.
Accepted full amplitude, Accepted half amplitude and EP full amplitude each
pass their field/traction/domain and baseline-subtracted space/timestep checks.
Their native analytical weak momentum-z audits fail; these remain failed
response qualifications. The complete supplement for the first Accepted
failure gathers both spatial audits, all four scalar observers and both frozen
quadrature controls without overwriting the first failed receipt. The EP and
half-amplitude audits likewise have complete independent coverage.

| Scaled primary response RMS at .016 | Accepted, 2e-5 | EP, 2e-5 |
|---|---:|---:|
| w | 3.324634e-5 | 3.319142e-5 |
| Pi | 3.773466e-4 | 3.450708e-4 |
| full normal traction Tzz | 3.724946e-4 | 3.417200e-4 |
| Rxx | 1.476798e-4 | 1.475258e-4 |
| Ryy | 7.423326e-6 | 7.262230e-6 |
| Rzz | 1.482754e-4 | 1.480380e-4 |
| Rxz | 1.899707e-4 | 1.892548e-4 |
| Qxz | 1.427917e-4 | 1.427989e-4 |
| full production | 1.918261e-4 | 1.909692e-4 |

All 17 observables retain their four-member comparisons. k and outer epsilon
responses are below the predeclared 5e-7 spatial signal floor; Qyy is identically
zero here. They are unresolved signals, not measured convergence orders.
Other listed core responses are resolved. Independent R/Q/P/traction and
radicand comparisons remain within 4.45e-16; actual A95 radicands on those
sampling grids remain above .6917828. Neutral analytical Qxz=0 is kept separate
from the actual perturbed nonlinear shear response and from Rxz=-u_star^2.

The initialized Accepted full-amplitude B5 divergence RMS is 2.87045e-16,
maximum 1.60821e-15. Its actual/effective state hash is
`f482078dcda8ec5403975a2f9b17e820aeb24f7b16c7e55d7e0880175a9d54f3`,
different from neutral reference
`b1743ac5799f15d30684aedafcef4d1e97662d646597f59f7e2a0f513a37590e`.
Explicit and implicit boundary hashes both equal
`efa95f6e5e2190ebbcebfb94c9d273db0107dab5a8ff92df4ea59f718ece2b8e`;
raw/effective and source/target boundary differences are zero. Fresh/resume
zero-step histories are byte identical. All other member identities are in the
complete comparison bodies rather than inferred from this representative case.

The half-amplitude test passes all 68 observable/member checks. Maximum
departure from half the full response is 7.68644e-9, below the declared
observable-specific limits. This is a bounded linear-scaling control, not a
claim of linear governing turbulence equations. On the common wall-comparison
grid, EP versus Accepted differences are Pi=3.44672e-5, Tzz=3.44384e-5,
w=7.69576e-7, Qxz=3.27634e-6 and P=6.13491e-6 on the primary pair. Same-wall
space/time uncertainty checks pass; different wall laws are not required to
predict equal responses. The grid resolves both fixed partitions and analytical
profile breakpoints, with exact-H inner/outer samples and one pressure gauge.
Moving cutoff panels remain part of native source quadrature and are recorded
separately; the computational mesh is not moved.

The independent EP frozen ownership control also diagnoses the same locator
failure: 127 misowned points, analytical weak momentum-z 9.60712e-10 falling to
5.04018e-16 with explicit containing-element bases. The one-point neutral
initialization fixture in the evidence directory reproduces geometry 4 being
selected for a strictly interior point in geometry 5, with no time stepping or
perturbation. Its forward closure remains valid; its ownership checker exits 1.
This is distinct from the earlier near-H diagnostic fixture, whose initial
failure and subsequent explicitly restricted scope are both archived.

## Independent physical scalar budgets

All 16 required Accepted full-response physical-budget reconstructions complete:
perturbed and neutral baselines, both spaces, both timesteps, volume orders
48/64 and physical-face orders 16/32. Reconstruction, finiteness, exact request
resume and quadrature/face-order controls pass. They include the dependent
epsilon boundary coefficient's storage rate and exclude the internal interface
from whole-patch external flux. Prescribed-row reaction values do not enter.

Their finite-defect accuracy does **not** pass the separately declared criteria.
For perturbed B5/B11, signed scaled frozen-stage defects are approximately
1.22755e-5/1.23242e-5 (inner k), 1.42630e-5/1.42644e-5 (outer k),
1.38970e-5/1.38974e-5 (outer epsilon), and 1.32304e-5/1.32508e-5 (whole k).
Spatial reductions fail or are unresolved and temporal uncertainty is unresolved.
On the baseline, outer k is below the floor and whole-k temporal uncertainty
passes; neither fact establishes full physical budget accuracy. A separate
derived classifier recomputes these numerical conditions rather than trusting
printed pass labels, preserves the earlier reconstruction receipt, and exits 1.
Four regression fixtures distinguish failure, uncertainty, floor and success.
Physical budget accuracy for EP/half amplitude is unattempted after the common
native weak-balance gate; no independent conservation theorem is claimed.

## Provenance, verification and cost

The trajectory runtime is unchanged from the byte-bridged prior numerical
snapshot. No governing turbulence, wall, patch, stage or acceptance source file
differs from the reviewed baseline. Thirty-four governing file hashes and all
loaded dependencies are checked by each admission. New C++ code is confined
to the read-only probe. Final reporting source is recorded independently from
each runtime's actual compilation source and compiler/link identities.

| Artifact | SHA256 |
|---|---|
| Frozen trajectory solver | `389dee61618f51c9d451be2ddbdc0e83ee049351683da2eae3b59ab0758afcb4` |
| Frozen descriptor probe | `a050e8af62360391588667ad168e2daa4dc433bc39679f6d326bf6ff33ae9279` |
| Spatial/response probe v3 | `6a24d6eaaeee1c129b338ed90526f41fddfece01b69c7636bc5efb5bd7ad86df` |
| Physical-budget probe v5 | `71501637551c05e5c948d81c2deb7665a5b5ea209b1881fc283821c431eb9dab` |
| Exact-interface diagnostic regression probe v7 | `4e461b4000dade3e3a033d197451f265b82e765149257c4e37e39f4626fbb67f` |

v3 was compiled from `b48fc2c7f2f3ea04848af2d2c0bb31aee335f4b8`, v5 from
`9034f1efa2353896069fa32036a267b4613bf829`, and v7 from
`c5ca6b597539b219389c48c20962867552129e23`. Actual current loader maps, source
hashes, compiler flags, CMake/compile/link hashes and runtime lineage accompany
the export. GCC 11.4 Release uses `-fno-fast-math -ffp-contract=off` and
`_GLIBCXX_ASSERTIONS`. Native work is rank 1, flat 2D, all numerical thread
variables explicitly 1. No shared installation is replaced.

Six targeted Python suites pass on Linux Python 3.10.12, as do CTest's duration,
evidence, physical-budget accounting, forward-interface and budget-accuracy
controls. The first CTest attempt without a built descriptor target remains
failed and archived; the target was then built in isolation and the affected
test passed. Command-line negative controls reject missing times, foreign
configuration, missing companions, truncated observers, failed residual
reductions, changed sealed evidence, existing output reuse, and neutral vectors
misrepresented as perturbed inputs for the intended reasons. The real final
scientific collector exits 1 and preserves four full accepted-state/history
failure witnesses. Reconstruction success never replaces spatial reduction.

Unchanged constitutive/wall, analytical weak-oracle, CLS, cache/rollback and
ordinary restart evidence is reused only through its byte/source bridges, not
rerun as a new broad campaign. Current neutral continuations and same-space
perturbation initialization controls additionally exercise genuine native
histories. Historical receipts and the eight epsilon temporal-order observations
are byte checked and remain unchanged.

Forty completed native segments contain 3192 accepted steps: A93 Accepted 1500,
A95 Accepted 1596, and A95 EP 96. Sum of instrumented exclusive process CPU is
22066.554799 seconds (6.13 CPU hours); diagnostic probes and orchestration are
additional. Native domain/physical/audit scopes account for about 65.5% of that
instrumented cost. Final-duration process CPU per step is approximately
6.1–6.6 seconds on B5 and 7.4–7.5 seconds on B11. Gross wall times include
resource suspensions and initialization and are not substituted for CPU cost.
No backend or state-dependent validity-cache optimization was introduced.
The owned guard is stopped only after checking there are no live or paused
owned numerical families; the final census and the earlier exception remain
in the evidence. Concurrent existing jobs are untouched.

## Reproduction and retained artifacts

The adjacent `a93-a95-e94e95-awf-nonlinear-duration-response-qualification-evidence`
directory contains `decision.json`, full `complete-results.json.gz`, compressed
extended summary, admissions, all four qualification failure witnesses, runtime
manifests, source/consumed-input/history identities, controls, selected complete
response comparisons, budget classifications, geometry diagnosis and the small
one-point fixture. `artifact-inventory.json` maps every exported artifact to its
Linux source and original hash. `SHA256SUMS` covers committed evidence. Heavy
packets, quadrature arrays, binary/checkpoint files and full logs remain under
the isolated Linux root above. No prior evidence file is replaced.

Reproduction on Linux, with those preserved frozen runtimes and artifacts:

```bash
R='/media/jay/新加卷/nektar-v2026-runs/a93-a95-e94e95-nonlinear-20261110'
S="$R/source/solvers/IncNavierStokesSolver/Tests/RANS"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 TMPDIR="$R/tmp"
export LD_LIBRARY_PATH="$R/runtime/lib"
export NEKTAR_RANS_PROBE_TIMEOUT=14400

# Exact resume: validates immutable requests/seals before reusing completed work.
# Do not increase concurrency without a fresh aggregate worker census.
python3 -B "$S/run_rans_balanced_duration.py" --root "$R" \
  --admission "$R/admissions/A95-N-Accepted-loader-v2.json" \
  --dt .002 --max-time .5 --concurrency 1 --resume \
  --receipt-tag init-timeout-recovery-v1 \
  --recovery-from "$R/comparisons/A95-N-Accepted-duration-0.002-0.1.json"
# This selects the recorded recovery request and preserves the original timeout.

python3 -B "$S/collect_rans_nonlinear_qualification.py" --root "$R" \
  --admissions "$R/admissions/A93-N-Accepted-loader-v2.json" \
    "$R/admissions/A95-N-Accepted-loader-v2.json" \
  --regression-receipt "$R/targeted-regressions-final-v6.json" \
  --forward-interface-receipt "$R/forward-interface-exact-scope-status.json" \
  --output "$R/export-reproduction-new"
# Expected exit 1: A93 common-w and A95/EP response-audit gates still fail.

python3 -B "$S/classify_rans_physical_budget_accuracy.py" \
  --source "$R/physical-budget-analyses/A95-N-Accepted-response-2e-05-dt-0.002-0.016/results.json" \
  --output "$R/physical-budget-accuracy-reproduction-new"
# Expected exit 1, with complete failed/unresolved numerical classifications.
```

The geometry diagnosis commands and immutable requests are copied under
`linux-controllers/run-A95-volume-ownership.py` and
`run-A95-EP-volume-ownership.py`. A new destination is required for a fresh
diagnosis; exact resume validates every completed artifact. The minimal fixture
can be rerun with the frozen v3 probe on its `case.xml request.json new-packet.json`
and checked with `check_minimal_ownership.py --packet new-packet.json --output new-result.json`.
The checker must exit 1 for `native_geometry_not_containing_point`; recognizing
that intended negative outcome is a successful reproduction, not a passing
native quadrature qualification. Full response/EP/half commands and the
separately initialized timestep identities are in `linux-controllers` and the
complete result bodies. Existing destination mismatches fail before overwrite.

## Scope and next capability gate

Demonstrated: A95 Natural Accepted neutral qualification through T=.5/theta=.1
on B5/B11 with both timestep companions; A93 physical evolution through the
same duration with a failed final common-w reduction; short A95 Accepted and
EP nonlinear responses with controlled field uncertainty and half-amplitude
scaling. Earlier E94E95 Accepted T=.5 and the six short natural configurations
remain qualified only within their prior scope. No new prescribed-kappa
duration or response campaign is claimed.

Failed or unresolved: A93's final common-w reduction; native analytical weak
balance for all three A95 short response cohorts; physical scalar-budget
spatial/temporal accuracy. A93 response is not admitted. A95 response beyond
.016, EP neutral duration beyond .016, and EP/half physical-budget accuracy are
unattempted after their stated gates. There is no hard-domain rejected timestep
in these failures; the preserved witnesses are accepted physical states and
their histories at failed post-step qualification checks.

The next bounded package should correct native integration-point ownership:
retain the generating geometry for volume points, distinguish strict interior
containment from tolerance-based interface selection, and verify Stokes,
partition, side and constraint identities on tiny moving-cutoff panels. This
changes the numerical operator and needs an explicit method/signature version,
fresh histories, ordinary CLS controls, and repeated short weak/space/time
qualification before longer response work. It must not rewrite these failures
or continue a changed operator from these old histories. Recheck A93's common-w
load on the corrected operator; if it persists, use the retained pressure,
stress, advection and common-row decomposition before another spatial choice.

These data do not justify another duration increase now. After ownership and
the residual/budget gates close, repeat bounded nonlinear response and then
introduce a multi-column neutral inflow/outflow test with an independent
physical face-flux budget before curved-terrain validation. MPI, true 3D,
nonlinear/action solvers, terrain, open-outflow preservation, global stability
and realizability/positivity theorems remain outside this evidence.
