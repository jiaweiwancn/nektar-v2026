# Approved-v2 inverse-cache and physical-history checkpoint

This is a numerical/infrastructure checkpoint after `6e290f5a82e8b6106885262e3245b0ea8d1563e0`.
It does not close H27-C/D or authorize the proposed v4 physical model.
`TASK3_READY=NO`; the full Task3 horizon has not been launched.

## Equations and numerical policy

The active model remains `FiniteNuIsotropicWallFullWJCollarV2`: physical nu
in both regions, isotropic lower-wall normal stress plus auxiliary shear,
the unchanged geometry-only C2 blend into full WJ below H, actual raw stress
production, and inner `Dk=nu+Cd*k*k/(sigma_k*Ew)`. The dependent amplitude,
v2 reference cutoff, inlet, conservative scalar transport, single outer
epsilon trace, finite-nu boundaries and full current admissibility gates are
unchanged. No historical-cutoff/auxiliary-diffusion v4 implementation is active.

`CoupledJFNKV1` now has an opt-in `OperatorDriftV1` inverse policy. It compares
all three current variable inverse-operator matrices with their prepared
matrices, using the maximum relative Frobenius change. It refreshes for cold
start, the actual BE/BDF2 mass change, drift>.02 or eight prepared Newton-base
uses. Current residuals, finite-difference Jacobian actions, forcing and final
physical guards never use stale equations. Reused epsilon factors retain
their matching full operator for the dependent-Gamma cross action. A typed
linear retry starts cold and refreshes at every Newton base. The default
`RebuildEveryNewtonV1` control is retained. See the prospective
[admission](wj-h27-newton-cache-admission-v1.md) and
[history checker specification](wj-h27-physical-cache-history-foundation-v1.md).

## Pinned execution

Linux namespace: `/home/jay/build/wj-h27-coupling-r1-20261117`.
Source13 scientific identity:
`11b6087b6e090b5cac0290ff3229cf2b83ad94ab4da57eb1eaee25087cfd715c`.
Runtime ELF SHA256:
`808364eb743afb60b9cdaa1c4e08c3b97ca9ba2d5bcde7390fcb7d0a87e3a8f5`.
Reference-table SHA256 remains
`b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.
The source manifest records both preserved primary checkouts. All1717 locally
present C++/header/CMake source files in the executed capsule match the Mac
working tree at this checkpoint. Compilation used the standard solver main,
Release, GCC11.4, `-ffp-contract=off -fno-fast-math`, four compiler workers;
the build took200.845039s. Science used one serial process/rank, all six
numerical-library thread counts1. Exact configure/build/run commands,
runtime libraries, XML/option files, hashes and per-operation resource receipts
are in the evidence bundle.

The physical fixture is the original14x50=700 cells, Q3/Q2, dt=.005,
nonlinear3e-12, linear1e-12, no Anderson, cached max uses8. All final physical
steps have four Newton updates and seven coupled Krylov iterations. This is
only T=.06 startup. It is not the required >=26.13350407707765 physical pilot
and .005 is not an accuracy-qualified production timestep.

## Independent numerical results

PolynomialScalarsV1/body1 (four steps, dt=.02) and Historical MMS (eight
steps, dt=.0025) both pass the fixed<=1e-8 fresh-versus-reused current/previous/
older field, explicit-history, essential-lift, Gamma/Ew/cutoff comparisons.
Their independent higher-quadrature current-equation audits pass. These runs
use max uses3 and exercise age refresh. The threshold-zero Polynomial control
also passes: every nonzero reported operator change refreshes, and its
independent current equations pass. No threshold was relaxed.

The physical fresh-versus-cached comparison passes all fixed<=1e-8 gates.
At step12 the scaled coefficient MAX differences are U1.08715e-14,
Pi1.47635e-9, k8.60936e-15 and epsilon3.27591e-15. Pressure is not shifted.
The cached uninterrupted trajectory records9 triple refreshes and39 reused
Newton bases; cold, mass-alpha, drift and age reasons are all observed.

| Advancing endpoint | Independent momentum | Continuity | k | epsilon |
| --- | ---: | ---: | ---: | ---: |
| Cached uninterrupted | 8.93282e-15 | 6.85934e-16 | 1.23969e-13 | 1.18000e-14 |
| Cold step3→12 | 1.03190e-14 | 6.64463e-16 | 1.28836e-13 | 1.22257e-14 |
| Retained recovery9→12 | 8.97632e-15 | 6.64942e-16 | 1.27089e-13 | 1.12224e-14 |
| Retry at step4→12 | 9.95182e-15 | 6.77462e-16 | 1.27994e-13 | 1.19564e-14 |

The independent physical transient reports keep overall `NOT_MET`: their
accepted time-discrete equation, history-boundary, cache/accounting and
quadrature subgates are `MET`, while physical/reference and refinement gates
remain unresolved. A passed subprocess receipt is not overall physical closure.

## Cold restart and actual retention

The independent physical history checker passes its explicitly limited
foundation gate. A three-step prefix plus cold nine-step process continuation,
retained recovery9 plus three steps, and the injected retry all reach step12
within the fixed<=1e-8 field/history/dependent-state limits. Current Pi scaled
coefficient MAX differences from uninterrupted are respectively1.27593e-9,
7.43504e-10 and1.13540e-9; U/k/epsilon differences are near machine precision.
These different inverse paths are not claimed byte-identical.

Recovery writes occur at3/6/9/12. Files3/6 were actually deleted, files9/12
remain, and scientific JSON at all four steps remains. Both retained files
were actually cold-imported. Recovery12's zero-step reload preserves all four
complete histories and its binary byte-for-byte; its empty step log, exact
`NumSteps=0` input and zero inverse work are verified. That reload does not
supply an advancing residual test. Recovery9's advancing endpoint does.
The retry performs the injected triple rebuild, discards it, and refreshes
each second-attempt base from unchanged accepted histories. Complete callback
timing survives; its explicitly incomplete engine-statistics flag is retained.

The geometry checker derives the actual14x50 mesh from the session's vertices,
edges and quadrilaterals. A static pre-execution correction removed an
incorrect assumed20x35 layout before the checker was first run. This test-only
overlay is hash-bound in `tools-cache-validation-v13/overlay.json`; no compiled
code, field, signed map, owner gradient or clock was changed. Its classification
allows only64-epsilon axis-span coordinate roundoff. The final checker SHA256 is
`4fc3ce8a66e911429b8ae2a139ccccb9b76a1d42317b0d4d9407a1755f493dac`.

## Complete short-run observations

All values below include scheduled step logging, recovery and scientific output.
They are measured physical progress, with the full per-step variation retained.

| Step block | Fresh factors, seconds/step | Reused factors, seconds/step |
| --- | ---: | ---: |
| Startup1–3 | 11.22612378 | 8.73781749 |
| 4–6 | 11.06795812 | 7.12910307 |
| 7–9 | 11.06821319 | 7.48358179 |
| 10–12 | 11.06248966 | 7.11894653 |
| All12 | 11.10619619 | 7.61736222 |

Whole timed native operations are137.798089s fresh and97.312824s cached,
a1.45801x mean advancing-step ratio. Cold9-step continuation costs73.424373s;
recovery9's3-step continuation29.420360s; zero-step recovery12 import5.872341s;
the12-step retry trajectory102.789771s. The cache's last step has zero factor
rebuilds and6.707825s native advancement, versus12 block rebuilds and10.999469s
for the fresh control. Its disjoint seconds are constitutive2.180888,
wall.180312, assembly1.815266, flow.103014, scalar.026577, guards1.851128,
residual.523868, nonlinear.024834 and other.001938. Complete elapsed time
also includes scientific/recovery output outside native advancement.

The three post-startup blocks are useful cache diagnostics. They do not sample
a qualified developed reference or finite perturbation and cannot support the
required complete per-trajectory/campaign forecast. No one-day PASS is claimed.

## Evidence and open gates

Published evidence: [native-newton-cache-v2e2](wj-h27-linux/native-newton-cache-v2e2/inventory.json).
The immutable Linux archive is14,222,544bytes, SHA256
`5f9e84c784f5ba1c688bd33c7006cc73be3d525071fc9c4874733f0e56674fea`.
Full raw physical fields, owners, audit rows and histories remain on Linux;
compact copies retain their exact raw identities. The first publisher failed
after4.336520s on a missing supplemental-account path; its input/log/receipt
and partial output remain. The corrected create-only publisher passed in
7.296528s. Both are charged. The post-publication account is still OPEN at
9722.40879069583 aggregate operation seconds, including earlier carried-forward
work; no budget or failed outcome was reset. The same43200s ceiling and3600s
protected evidence reserve apply.

The v2 same-inlet reference still fails the unchanged inner-k engineering
envelope and its numerical refinement budget. This checkpoint does not change
that evidence. Final h/p/q/dt qualification, required evolved physical duration,
baseline/startup/both-sign perturbations, monitoring escalation, terminal
handling, storage and complete coupled cost remain open. The next independent
terminal-output foundation is separately admitted; v4 science remains pending.
The original Task3 horizon505.02466521278404 and terminal window
[447.8483832973273,505.02466521278404] are unchanged. No qualified Task3 launch
command exists yet.
