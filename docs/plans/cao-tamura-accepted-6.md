# CT-Accepted-6: conservative curved Accepted manufactured patch

Starting source: `7d32b43dc9d176628cfef1dbae2aee02a86ce330`, clean and equal to fetched remote.
The [versioned contract](cao-tamura-accepted-6-evidence/contract.json) precedes native execution.

This package implements a guarded manufactured curved route. The complete
evolving formulation is **not qualified**: the declared degree and timestep
gates fail in the A95 primary experiment. Independent quadrature receipts are
being completed separately. CT4/CT5 physical fits, observations, sealed artifacts and
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

## Bounded A95 evolution and failed consistency gates

All five fresh A95 histories reach `T=.016`, passing the accepted-state domain
and manufactured field/traction limits. Same-space timestep members have
identical actual initial coefficient vectors. Each has complete native
physical-budget equations and actual captured-stage reconciliation. Full
four-control quadrature qualification is a separate receipt and is not supplied
by native-only completion.

The following are scaled weighted physical RMS errors against the continuous
time-dependent solution, on one common physical inventory. Pi is absolute;
there is no pressure mean subtraction. Scales are `.2` for U and w,
`K0=.6115200681722618` for k, `.3167838556291072` for outer epsilon and `.1`
for Pi. The limits are `.005` for the four transported fields and `.01` for Pi.

| Space | Stored slots | dt | Accepted steps | U | w | k | Outer epsilon | Pi |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| P3/P1 | 1273 | .002 | 8 | 1.52414e-5 | 1.40619e-5 | 1.09086e-5 | 1.20792e-4 | .00158203 |
| P3/P1 | 1273 | .001 | 16 | 8.97399e-6 | 6.91036e-6 | 5.49718e-6 | 1.29468e-4 | .000775087 |
| P4/P2 | 2283 | .004 | 4 | 2.56345e-5 | 2.04337e-5 | 2.17874e-5 | 6.09391e-5 | .00263361 |
| P4/P2 | 2283 | .002 | 8 | 1.22157e-5 | 1.10066e-5 | 1.09031e-5 | 2.99731e-5 | .00159824 |
| P4/P2 | 2283 | .001 | 16 | 5.25802e-6 | 4.98739e-6 | 5.45498e-6 | 1.67812e-5 | .000775897 |

The original gates are retained. They apply to resolved differences rather than
successful process exits or finite values.

| Observable | P4/P3 error ratio at dt=.001 | Required <=.85 | P4 triplet order | Required .8--1.2 |
|---|---:|---|---:|---|
| U | .585918 | pass | 1.208038 | **fail** |
| w | .721726 | pass | 1.234440 | **fail** |
| k | .992324 | **fail** | .998575 | pass |
| Outer epsilon | .129616 | pass | .997092 | pass |
| Absolute Pi | 1.001045 | **fail** | .357693 | **fail** |

Every temporal difference is above the declared `1e-6` signal floor. Neither
the small U/w band misses nor the pressure failure is rounded into a pass.
The k and Pi errors exceed the spatial floor `1e-7`; they are not relabelled as
unresolved-floor observations.

For k, the timestep spread is about `5.45e-6` on either space, while the
P3/P4 difference is `6.55e-7`. For Pi, these quantities are about `.00082`
and `.000125`. Thus the raw degree contrast is substantially contaminated by
time/splitting error. Read-only first-order extrapolation and second differences
are exported as diagnostics; they do not replace the failed gates or alter any
state. Two degrees do not establish a spatial order.

P4/P2 initial projection errors are `2.11e-9`, `3.00e-8`, `2.43e-12`,
`1.89e-5` and `1.17e-7` on these same field scales. Initial projection,
subsequent numerical drift and time-dependent reference change are reported
separately. Initial pressure value error alone is much smaller than the evolved
pressure error; this does not exclude an initial weak-continuity or splitting
effect. No violated governing discrete identity has been established from the
temporal failure.

At P4/P2 dt=.002, ground-traction reference error is `.000742`;
the interface inner/outer reference errors are `.000910/.000906` and the full
jump is `.000310`, all below `.005`. These are reconstructed physical
tractions using current FE gradients and endpoint K, not prescribed loads or
row reactions.

The last time-discrete scaled stage residuals at dt=.002/.001 are about
`.0130/.00389`. They are not steady-state residuals. Unmodified alpha-zero
spatial loads and the separate read-only consistency load
`r_spatial + load(analytical time derivative)` remain visible; the latter never
enters the governing equations.

Because the A95 consistency contract fails, A93 and E94E95 manufactured
**evolutions are not admitted**. Their successful kernel, continuous-source,
initialized native weak-balance and static-transmission results retain their
own scope. No additional timestep, mesh or manufactured input has been used to
search for a passing result.

## Compatibility and rejection controls

The rebuilt runtime reproduces fresh CT5 A95 physical states at `.002/.004`
and its native history payload bit for bit. The current curved same-method
restart is also bitwise identical. Deliberate rejection at step two retains the
last accepted state and every native history slot; its rejected candidate equals
the corresponding uninterrupted next state. Changed geometry is rejected by
the checkpoint signature. These tests do not extend a CT5 trajectory.

Both registered CT6 CTests pass, including the admission suite's process-level
negative checks. Missing members/equations/controls, invalid primitives,
incorrect actual-stage binding, failed independent traction and failed shared
prerequisites cannot admit a model transfer merely through enclosing pass labels.

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
