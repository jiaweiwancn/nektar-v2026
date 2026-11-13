# CT-Accepted-6: conservative curved Accepted manufactured patch

Starting source: `7d32b43dc9d176628cfef1dbae2aee02a86ce330`, clean and equal to fetched remote.
The [versioned contract](cao-tamura-accepted-6-evidence/contract.json) precedes native execution.

This package implements a guarded manufactured curved route. Qualification of
the complete evolving formulation is pending the declared degree, timestep and
physical-budget gates. CT4/CT5 physical fits, observations, sealed artifacts and
historical decisions remain unchanged. No ordinary curved production or
experimental admission is implied by a manufactured test.

## Formulation and implementation

`CT6ManufacturedV1` permits only original Natural A93/A95/E94E95, Accepted AWF,
serial 2D, constant physical height and the declared Gaussian/zero map. The
boundary policy is `AcceptedManufacturedOpenV1`; the continuous provider is
`CT6ContinuousV1`. General curved production, new-model curved EP, MPI, 3D,
nonlinear/action routes, ExactStepV1 caching and full variable-height evolution
remain guarded. The 3,000 stored-slot limit is unchanged.

The manufactured graph is `h(q)=.005 exp(-q²)`, with ground coordinate
`-1.2 <= q <= 1.2`, physical interface height `.5` and mapped top `r=1`.
It has varying curvature of both signs and nonzero tangential metric
derivatives. With `B_g=grad_Gamma n`, `gamma=-h''/(1+h'^2)^(3/2)` and
`J=1+r gamma`; the profile kernel's older `1-kappa*r` notation uses
`kappa=-gamma`. Positive stretch and global normal-map reach cover the outer
region too. This geometry is not the Cao–Tamura hill.

The kernel retains finite inner molecular viscosity and the model-compatible
wall constants. Its cutoff solves the monotone metric budget
`integral(J E_w)=integral(J p0)`; `p0` is counterfactual zero-forcing production.
Actual production remains `-R:grad(U)` with the complete Cartesian tensor and
raw numerical divergence. Positive column budget, positive reference k and
actual-flow covariance are separate checks. `Profile()` now changes Accepted
`sampleHeight`, `tailStart` and `height` together when local height changes.

The evolving fixture uses a time-dependent streamfunction with
`u_t=psi_r`, `w=-psi_q/(a J)`. There is no extra `1/J` in `u_t`.
The inner Hermite coefficient matches tangential total traction; the pressure
jump matches normal traction. K is common, the two diffusion coefficients
match a nonzero k flux, and outer epsilon has the dependent trace `A(K)/H`.
There is no inner epsilon equation or extra epsilon-flux condition. The local
transported inner k is nonconstant and differs from the frozen reference column.

Manufactured sources come from continuous Cartesian derivatives, using
multivariate Taylor arithmetic independently checked by arbitrary-precision
differentiation and a Lambert-W cutoff. They never use the native discrete
residual. Inner body forcing explicitly removes the counterfactual roughness
completion so that adding the actual native completion recovers the intended
continuous equation. All corresponding work is retained. This forcing is
available only through the named manufactured capability.

Mapped Cartesian gradients and oriented cofactor measures are used once in the
original native residual and one-sweep assembly. The roughness pressure action
now includes its projected Cartesian normal component on the named curved
route; the flat arithmetic path is retained. Mechanical pressure sampling uses
`grad(Pi)-2 grad(k)/3` at fixed physical r before the one-sided outer trace.
No mean-pressure constraint or extra outlet pressure is introduced.

`GeneratedVolumeOwnerV1` remains intact. Generating geometry, field and explicit
interface sides survive through evaluation rows and caching. Curved capability,
map and provider identities are fingerprinted separately from closure identity.

## Completed independent checks

The immutable [contract](cao-tamura-accepted-6-evidence/contract.json) and
[pre-evolution addendum](cao-tamura-accepted-6-evidence/pre-evolution-addendum-v1.json)
define the tolerances. The [execution notes](cao-tamura-accepted-6-evidence/execution-notes-v2.json)
record diagnostic repairs without changing those tolerances.

| Check | Completed evidence |
|---|---|
| Independent continuous fields/sources, all models | 18 samples; maximum divergence/transmission/gradient defect `3.76e-37`; positive sampled k, epsilon and covariance, A95 radicand `>=.97556` |
| Metric Accepted kernels | Nine model/signed-curvature combinations and 15 intended boundary/rejection controls; maximum scaled native/oracle difference `6.40e-13` |
| Native conservative all-row operator | P3/P1: `1.09e-15`; P4/P2: `1.08e-15`, including separate prescribed-row replacement |
| Continuous weak equations on native tests | Maximum `2.09e-13`; shared-interface cancellation `1.83e-17` |
| Geometry, cofactors and every-mode cell Stokes checks | Maximum Stokes defect `2.29e-12`, against `2e-11` |
| Actual small consistent-free-mass pressure action | Full rank for 12/30 pressure unknowns, no additional null mode; constant-pressure action is nonzero |
| Zero-curvature profiles, stresses, production, raw loads and frozen matrices | All three models pass the `2e-10` comparison against the ordinary flat route on identical states/data |
| Static variable-height transmission | All models in P3/P1 and P4/P2 pass independent full-traction, nonzero conservative k-flux, trace-chain, pressure recovery and dependent-epsilon checks |

The static test uses `H(q)=.5+.03 sin(pi q/2)`. It checks arbitrary polynomial
geometry probes separately from constructive compatible local Hermite collar
jets. It detects both omitted tilted diffusion and fixed-chart differentiation.
A zero base-normal derivative produces nonzero tilted-interface diffusion.
These are static transmission results, not full variable-height FE evolution.

The small operator fixtures have 188/319 stored slots. The evolving P3/P1 and
P4/P2 fixtures have 1,273/2,283 slots. Native counts are measured, not predictions.
Native chart-basis checks are distinguished from independent mapped physical
gradient and Stokes checks.

## Evolving qualification in progress

Fresh histories use the prescribed degree/timestep matrix through `T=.016`.
A95 P4/P2 at `.002` has completed all eight accepted steps. Its independent
endpoint scaled RMS errors are approximately `1.22e-5` (U), `1.10e-5` (w),
`1.09e-5` (k), `3.00e-5` (outer epsilon) and `.001598` (absolute Pi).
Ground traction error is `.000742`; the two interface-side reference errors
are `.000910/.000906`, and their full jump is `.000310`. All pass the CT6
limits. These are bounded endpoint accuracy results, not a convergence claim.

The companion histories, degree comparison and full frozen-quadrature coverage
are still required. A native-only physical-budget diagnostic explicitly exits
3 and cannot supply full qualification or admit the A93/E94E95 transfers.
Time-discrete stage residuals, unchanged alpha-zero spatial loads and the
read-only consistency load `r_spatial + load(analytical time derivative)` are
reported separately. The last expression never enters evolution.

## Retained failures and diagnostic corrections

The first attempted evolving step stopped before any governing solve on an
inherited flat-only continuous ground-k guard. Extrema of the ground chart
polynomial are unchanged by monotone physical arclength. The CT6 route now
uses that continuous check, including a narrow negative interval between
ordinary sample points. Its arithmetic allowance and admissibility floor
are unchanged.

An evolved replay originally checked the immutable time-zero reference against
time-T live boundary data. The diagnostic now validates the reference at zero
and observes the unchanged trajectory at T. A native JSON-array reporting bug
was also corrected. Failed preparations and original packets remain archived.

The first scalar reconciliation matched its own reassembled frozen system at
roundoff, but differed from the actual captured governing k matrix by
`3.3069e-6`. Old coefficients must include current-time essential boundary data,
while dependent epsilon trace rows retain their existing corner precedence.
Physical storage still uses the actual previous accepted field. Binding the
original post-flow/pre-k snapshots and captured matrix reduces the difference
to `9.42e-16`. This changes diagnostics only. Solver executables after this
correction are byte-identical to the sealed governing runtime. CT1–CT5 used
time-independent inputs and their historical budget decisions are not relabelled.

## Provenance and remaining scope

Linux is authoritative. All heavy data are isolated under
`/media/jay/新加卷1/nektar-v2026-runs/ct-accepted-6-20261113/recovered-v1`.
The original interrupted CT6 preparations remain beside that directory.
The user-coordinated filesystem recovery changed the mount driver to ntfs-3g;
it did not format, repair or delete experiment data. CT4/CT5 archives, frozen
runtimes, shared installations and unrelated campaigns remain untouched.

Builds use two jobs; numerical libraries use one thread. Every numerical family
uses the owned 16-worker project guard under the 72-worker aggregate limit.
Final resource closure, source/runtime/library hashes, commands, accepted-run
costs, failure witnesses and the complete qualification decision will accompany
the completed numerical campaign.

Ordinary curved production, experimental predictions and developed precursors
remain unadmitted. A93's narrow CT5 value `.000992850582/.001`, general spatial
convergence/long-time preservation, the historical CT4 pressure interpretation,
roughness-sublayer/interface adequacy, the `.1` thin-layer criterion, smooth
preparation and approach-layer fidelity remain open. EP, optimization, full
variable-height evolution and backend expansion remain outside CT6.
