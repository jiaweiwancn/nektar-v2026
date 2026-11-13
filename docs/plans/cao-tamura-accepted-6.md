# CT-Accepted-6: conservative curved Accepted manufactured patch

Starting source: `7d32b43dc9d176628cfef1dbae2aee02a86ce330`, clean and equal to fetched remote.
The [versioned contract](cao-tamura-accepted-6-evidence/contract.json) precedes native execution.

This package implements a guarded manufactured curved route. The complete
evolving formulation is **not qualified**: the declared degree and timestep
gates fail in the A95 primary experiment. All five endpoint physical-budget
audits completed with four integration controls; their accuracy and integration
uncertainty remain distinct from reconstruction. CT4/CT5 physical fits, observations, sealed artifacts and
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
residual. Inner body forcing explicitly removes the analytically evaluated
manufactured roughness completion so that adding the actual native completion recovers the intended
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
by native-only completion. Each completed audit contains `native`, `q48f16`,
`q64f16` and `q64f32`, each with inner k, outer k, outer epsilon and whole k,
plus both k and outer-epsilon reconstruction records. Whole-patch k excludes
the internal-interface flux from the external total. Physical fluxes use actual
one-sided FE gradients, including the ground; prescribed-row reactions never
replace them. Dependent-epsilon boundary storage and current-versus-frozen
coefficient splitting remain explicit.

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
The primary `.002` endpoint's refreshed dependent-epsilon trace lag is
`5.35e-5` absolute (about `1.69e-4` on its epsilon scale). The relevant frozen
stage constraints and this refreshed lag are distinct; there is no post-k
trace repair.

All five four-control endpoint audits pass the physical budget limit `.02`.
The largest scaled defect is outer epsilon: `.00350792` for P3/P1 and
`.000404766` for P4/P2. Actual captured-stage matrix reconstruction is at most
`2.95e-15`, against `2e-10`. Independently integrated net mass defects are at
most `2.02e-16` on the declared `.2 * physical outlet length` scale.
Each member also passes all 45 frozen integration-uncertainty comparisons.
These compare current unanchored mass-dual residual norms, complete frozen
scalar load vectors and physical budget defects, without changing effective
coefficients, immutable reference or live essential boundaries. Momentum norm
agreement alone is not presented as a direction-sensitive load comparison;
the separate all-row operator check supplies that reconstruction test.

Weak mass conservation does not imply pointwise incompressibility. At the
P4/P2 `.001` endpoint the independent physical grid gives raw-divergence RMS
`1.78117e-5`, maximum `1.09094e-4`, and isotropic production-work RMS
`7.43328e-6`. The signed `-2 k div(U)/3` term is retained. Its coordinate,
region, native field owners and side are exported; no trace removal is applied
to production.

For the primary P4/P2 `.002` member, aggregation over every recorded accepted
state and audited stage gives inner/outer k minima `.587161/.605405`, outer
epsilon `.156022`, inner/outer covariance minima `.370790/.270410`, and A95
radicand `.956606`. Final-state extrema are retained separately. These are
finite native inventories with witnesses, not a global positivity or
realizability theorem.

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
The derived integration gate recomputes the original uncertainty bound from
the measured degree difference; a forged large limit or a passing enclosing
label cannot override it.

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

The first outlet-orientation negative control changed an observation inventory
that the all-row reconstruction did not consume. That ineffective test is
retained as a failed harness attempt. The corrected control reverses the actual
natural-boundary cofactor rows, verifies that rows were changed, and observes
the intended reconstruction failure. Its new receipt does not overwrite the
original attempt. All 17 corrected discriminating controls pass.

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
costs and failure witnesses are in the
[sealed compact export](cao-tamura-accepted-6-evidence/qualified-v1/decision.json)
and [final primitive verification](cao-tamura-accepted-6-evidence/final-audit-v2/results.json).
The admission CLI exits **1** for the actual failed convergence primitives;
the completed collector exits **3**, by design, rather than implying admission.
The final verification reproduces the classification exactly and validates
all sealed compact hashes. There are no missing required compact receipts.

The compiled diagnostic source is
`0cce5b53716db58e39ee9ee0eb5c14f29f62e668`. The governing executable is
byte-identical to the ground-guard-corrected runtime built at
`99aac0832ad009715f386040985267e99042444c`. Complete endpoint audits and the
comparison use analysis source `764180e0c4d6d8e67a5fdce9bab6fa874d3a5904`;
the corrected negative controls, 18 admission tests and final primitive export
use `e8e49ae527ca73870fd65c16f04e0917dfa91ed6`. Later report/evidence commits
do not change tested C++ sources. Runtime SHA256 values are:

| Frozen runtime file | SHA256 |
|---|---|
| `IncNavierStokesSolver` | `c4ad1b18f4a57e880700bfb585cdb12a4f2ca4fb12401eb988153c7ca6e31601` |
| `RANSPreservationProbe` | `c4ee55a14a8ee254a97f7ce69364c01ac150cb3a2260cd428017780c9b03fb60` |
| `RANSCurvedAcceptedKernelProbe` | `720b862961479adc5e24327742b981fcbed37e42ba227ccca463132f4eb4fd80` |

Build requests retain compiler settings (GCC 11.4.0, Release `-O3 -DNDEBUG`,
GNU C++17 and the existing warning/error policy); each command binds consumed XML,
provider/boundary functions, state and native history. Actual mapped-library
paths and hashes are in the compact library registry. MPI is enabled in the
build, but every numerical test is serial rank one; this is no MPI qualification.
`RANSModels.{cpp,h}`, the two specifications, CT5 evidence and frozen observation
data have no changes from the reviewed start.

The recovered campaign's closed monitors record a maximum of **14 project**
and **30 aggregate** workers, including control helpers, below 16/72. The final
verification family closes at 3/3; its own completed receipt accompanies the
otherwise closed 25-family inventory. CPU time, elapsed time and suspension
intervals are recorded separately, with PID/thread high-water witnesses.
Historical resource failures are unchanged.

The following are actual solver-process costs for initialization plus accepted
evolution, with ordinary checks enabled. Separate read-only endpoint probes
are excluded from these rows; the full exclusive instrumentation and analysis
costs remain in the archive. Nested inclusive timers are not added.

| Space | dt | Steps | Process CPU seconds | Elapsed seconds |
|---|---:|---:|---:|---:|
| P3/P1 | .002 | 8 | 1013.40 | 1015.56 |
| P3/P1 | .001 | 16 | 1769.72 | 1773.67 |
| P4/P2 | .004 | 4 | 615.43 | 617.88 |
| P4/P2 | .002 | 8 | 1035.86 | 1040.27 |
| P4/P2 | .001 | 16 | 1838.33 | 1846.09 |

These synthetic-source costs are not an experimental hill cost prediction.
For example, the P4/P2 `.002` solver execution peaks at 598,504 KiB RSS.
Per-execution RSS, initialization/projection and exclusive timer categories are
retained with the native member receipts. No backend or performance
optimization was introduced.

## Required-scope decision and next gate

| Scope | Decision |
|---|---|
| Curved Accepted kernel, all three models | **Pass**, signed curvature and intended rejection controls |
| Conservative spatial operator | **Pass** in the stated small P3/P1 and P4/P2 native fixtures |
| Full traction and transported-k transmission | **Pass** for the constructed continuous/native/static tests |
| Constant-height evolving MMS, A95 | Bounded completion and physical checks pass; **degree/time consistency fails** |
| Constant-height evolving MMS, A93/E94E95 | **Not attempted** because complete A95 admission fails |
| Static variable-height transmission | **Pass**; evolving variable height remains rejected |
| Zero-curvature/flat regression | **Pass**, including bitwise fresh CT5 A95 state/history compatibility |
| Resources | **Pass** for the recovered isolated campaign |
| Experimental applicability / ordinary curved production | **Not admitted** |

The smallest next package is a bounded diagnosis of absolute-pressure time
level, initial weak continuity and frozen-boundary splitting on this same MMS,
followed by an independently initialized finer timestep contrast that separates
time error from the degree comparison. Existing state/row/matrix packets are
the first diagnostic inputs. A governing change needs a demonstrated discrete
identity defect; no such defect is inferred merely from a failed error ratio.
Only after the evolving consistency gate passes should the remaining model
transfers precede physical approach-flow/interface admission and a first bounded
rough-hill pilot. That later pilot also needs a reach-safe outer geometry, rather
than extrapolating this shallow Gaussian normal map to the experimental hill.

Ordinary curved production, experimental predictions and developed precursors
remain unadmitted. A93's narrow CT5 value `.000992850582/.001`, general spatial
convergence/long-time preservation, the historical CT4 pressure interpretation,
roughness-sublayer/interface adequacy, the `.1` thin-layer criterion, smooth
preparation and approach-layer fidelity remain open. EP, optimization, full
variable-height evolution and backend expansion remain outside CT6.
The original CT5 A93 coarse-pair ambiguity and all eight historical unresolved
epsilon-order observations are unchanged. CT6's time-dependent manufactured
orders use their own declared horizon and signal floor.
