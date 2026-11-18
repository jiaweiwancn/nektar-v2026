# H27 source19 inverse guards and numerical-evidence reuse

The changed source19 reference-inverse input path passes its focused native
and independent tests. The audit also binds reusable mixed-boundary/gauge/Q3-k
and constitutive/admissibility evidence to their unchanged numerical sources.
**This is progress, not full H27 closure. TASK3_READY=NO.** The active model
remains the user-approved v2; v4 remains inactive.

## Exact native test

Parent is remotely verified `d205335c80458009387430d2fa54f3c3e39ef560`.
The application and its source19 libraries are unchanged. The separate probe
links those exact existing libraries, uses source19 headers and numerical
identity `56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf`,
and builds with C++17 Release, no-fast-math and fp-contract=off. The selected
application ELF remains SHA256
`37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b`.
The table remains SHA256
`b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d`.

The new diagnostic ELF is 240448 bytes, SHA256
`a0ef00e11367e765a75f962102890dc3bd8431ccc95491bd2018ef62c6ad986d`.
It is a reference-API probe, not a replacement application. Source, existing
static libraries, loader, compiler commands, options and link identities
are in [probe-runtime.json](wj-h27-linux/native-reference-inverse-v2l/probe-runtime.json)
and the frozen project manifest. All numerical execution was on Linux.
The small standalone [CMake project](../../solvers/IncNavierStokesSolver/Tests/RANS/h27/reference_inverse/CMakeLists.txt)
allows replay without modifying the immutable application build.

The [prospective admission](wj-h27-reference-inverse-guards-admission-v1.md)
fixes the input controls, independent checks, 1e-10 limits and <=400-second
matrix before native outcomes. The [probe](../../solvers/IncNavierStokesSolver/Tests/RANS/h27/reference_inverse_probe.cpp)
queries 129 admitted K values and 129 distinct radius values. The
[independent checker](../../solvers/IncNavierStokesSolver/Tests/RANS/h27/reference_inverse_check.py)
uses direct Chebyshev definitions, polynomial extrema and SciPy brent roots,
without native bisection code or native inverse levels.

| Quantity | Maximum scaled error | Limit |
|---|---:|---:|
| c(K) / H | 6.938893903907228e-17 | 1e-10 |
| Inverse K / Kin, including independent roots | 1.3784814909780096e-15 | 1e-10 |
| Inverse residual / H | 6.938893903907228e-17 | 1e-10 |
| Relative cutoff derivative | 4.512506311378652e-16 | 1e-10 |

Exact endpoints and an endpoint query after visiting the opposite endpoint
pass. The independent whole-interval negative slope enclosure agrees with
the native enclosure and contains the independently evaluated polynomial
derivative extrema. This is a bounded table certificate, not a global
physical-reference root theorem or an all-system cache-return test.

Seven radius controls reject: nextafter outside each endpoint, negative,
zero, NaN, positive infinity and negative infinity. Four synthetic tables
exercise the new constructor path:

| Synthetic input | Default application identity | Diagnostic explicit-hash constructor |
|---|---|---|
| Constant cutoff | Reject | Reject: no certified decreasing inverse |
| Increasing linear cutoff | Reject | Reject: no certified decreasing inverse |
| Quadratic with two inverse roots | Reject | Reject: no certified decreasing inverse |
| Stale derivative coefficient | Reject | Reject: analytic derivative coefficients differ |

The quadratic has two independently found distinct roots inside the admitted
K interval at its declared radius level. Synthetic tables are deliberately
unqualified controls. Only the probe passes their matching explicit SHA to
reach a deeper guard; the application retains its default golden SHA. These
controls are not alternative model data or new physical reference solutions.

## Reused evidence and scope

The [Linux reuse audit](wj-h27-linux/native-reference-inverse-v2l/science/reference-inverse-v19/reuse-audit.json)
checks actual immutable source7/source8/source19 bytes against their source
manifests and binds original native results to their independent reports.
It verifies 11 relevant numerical source/header files for the source8 mixed
tests and 14 for the source7 constitutive/admissibility tests. The selected
source19 application runtime remains separately pinned.

* All three closed/slip/physical-outlet tests retain independent small actions
  below 1e-10. Closed/slip each have one unaugmented pressure null mode;
  physical outlet has none. The slip/outlet native v maps have 11/13 permuted
  ids. The full Q3 known-k action, -2/3 sign and outer H owner remain tested;
  Q2 truncation, missing transpose and flipped signs remain negative controls.
* The 67-state planar audit retains four negative-covariance and one negative-
  principal witnesses, spanwise rejection and independent tensor/tangent/minimum
  checks. This supports the declared flat planar application; it is not an
  arbitrary 3D flow qualification.
* The unchanged optimized WJ kernel retains native value/direction parity
  7.346672369403494e-15/1.1368683772161603e-13 and independent joint-FD
  error 1.0800093702989166e-10 under its separately frozen 2e-8 FD tolerance.
  Three value and three direction fallback controls remain. This finite-
  difference tolerance is not substituted for the 1e-10 small-action target.
  Kernel speed measurements are not transferred into coupled forecasts.

Changed integration and trace paths use the newer source19 proofs:
[temporal/spatial/current continuous MMS](wj-h27-source19-mms-checkpoint-v2.md),
[native cutoff panel actions](wj-h27-cutoff-panels-checkpoint-v2.md), and
[actual owner transition with inverse reuse/cold history](wj-h27-cutoff-owner-transition-checkpoint-v2.md).
This avoids rerunning successful unchanged kernels merely because their
surrounding executable grew. The aggregate numerical gate remains explicitly
unclosed: this report does not establish an explicit complete transported-
scalar A-to-B-to-A action/cache control or replace the final full-requirements
review. Existing prescribed-family H26 tests and the pure inverse query above
must not be described as that complete coupled control.

## Complete tranche cost and preservation

| Operation | Whole charged seconds |
|---|---:|
| Fixtures | 1.7924106270074844 |
| Configure thin project | 1.7670248359791003 |
| Compile/link probe | 3.4715040100272745 |
| Native inverse and rejection controls | 1.8016111379838549 |
| Independent inverse assessment | 1.785138117033057 |
| Source/evidence reuse audit | 1.784675469971262 |
| Publication | 3.5916875419788994 |
| Total | 15.994051739980932 |

All seven operations pass; no failed attempt, overlap or allowance reset is
omitted. Native analysis times and kernel timings are not production cost.
The exact [post-publication account](wj-h27-linux/native-reference-inverse-v2l/accounting-post-publication-open.json)
is **OPEN16991.556598995754/43200 seconds**, with no RESERVED operation and
3600 seconds protected for evidence. Remaining non-protected allowance is
22608.443401004246 seconds. The goal remains active and is not marked complete.

The phase uses one science rank, six library thread limits 1, four compiler
workers, at most 8 actual shared-host workers and reported maximum tree RSS
314953728 bytes. Source/build/runtime including the thin project uses
3002450622 bytes; science 2726725508 bytes; free 69430677504 bytes.
The thin project/build adds 675145 bytes. Both unrelated primary working
trees and older capsules/accounts remain preserved.

The compact archive is 150779 bytes, SHA256
`45cd47b58f0f54a1f0008b39ae24dc38919e1522ee49d8ee61984d4e93830daf`.
All 50 inventoried files were hash-verified and retained. The final publication
receipt/log and exact post-publication account are separate additions after
the immutable archive. Large existing libraries, native fields and old
results remain on Linux at their bound identities.

The physical reference still has inner-k RMS/MAX
.13279468209091205/.44168324293315253 against .05/.10. Its selected-source19
physical h/p/q/dt qualification, four required-duration pilots, monitoring
escalation, persisted statistics/storage and complete costs remain open.
The unselected .005 short-work proxy still exceeds the one-day target.
No additional physical refinement was admitted: the existing 5600-cell steady
already used 5.470GB RSS, so a blind doubled grid lacks a safe memory forecast.
This is not a proof that finer solutions cannot be computed within the cap.

The [Task3 launch hold](wj-h27-task3-launch-hold-v2.json) remains unchanged in
all global gate statuses, horizon/window and null production commands.
No full Task3 run was launched. The [v4 scalar/reference proposal](wj-h27-historical-inner-transport-proposal-v4.md)
remains inactive pending a separate model decision; the approved tensor
amendment does not authorize it.
