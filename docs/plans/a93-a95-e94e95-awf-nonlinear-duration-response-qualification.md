# Nonlinear duration and response qualification

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
