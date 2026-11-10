# Nonlinear duration and response qualification

A95 Natural / Accepted passes the complete four-member neutral qualification
through T=.5 (theta=.1). A93 reaches the same duration with passing physical
accuracy, admissibility and timestep bounds, but fails the fixed-common vertical
momentum reduction; its response is not admitted. A95's T=.016 controlled
response resolves pressure, vertical velocity, nonlinear shear stress and full
production with passing space/timestep uncertainty. Its conditional extension
and selected EP/half-amplitude controls are reported separately below.

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
